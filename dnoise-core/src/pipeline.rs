//! The per-frame denoising pipeline: [`Denoiser`] holds one run's settings and
//! prebuilt state, and runs every enabled stage on one frame at a time.
//!
//! # Run state versus frame state
//!
//! A [`Denoiser`] is built once per run from plain metadata ([`FrameMeta`] for
//! every frame, in `Frames.Id` order) and the parameters. Optional per-run
//! state is attached with the `set_*` methods before the first frame:
//!
//! | State | Built from | Needed for |
//! |---|---|---|
//! | [`MsmsKeep`] | every ddaPASEF MS/MS frame, grouped by precursor ([`crate::msms::MsmsKeepBuilder`]) | ddaPASEF MS/MS denoising |
//! | [`PrmWindows`] | `PrmFrameMsMsInfo` | prm-PASEF MS/MS denoising |
//! | [`DiaWindows`] / [`DiaRegions`] | `DiaFrameMsMsInfo` / `DiaFrameMsMsWindows` | diaPASEF window gates and per-window filtering |
//! | [`DiaMs1Gate`] / [`PolygonGate`] | isolation windows or the IMS polygon, plus calibration | the MS1 gates |
//! | [`CropGate`] and the RT mask | crop bounds, plus calibration | cropping |
//! | [`NeighborIndex`] | frame compatibility metadata | temporal neighbor support |
//!
//! Everything else is per frame. After setup the denoiser is read-only
//! (`&self`), so frames can be processed in any order and in parallel.
//!
//! Two stages read **other frames** while denoising frame `i`:
//!
//! - ddaPASEF MS/MS denoising pools each precursor's scans across all the frames
//!   that isolated it. That pooling happens up front, in [`MsmsKeep`]; the
//!   per-frame call only looks up the result.
//! - Temporal neighbor support ([`NeighborIndex`]) sums up to `radius`
//!   compatible frames on each side of frame `i`. [`Denoiser::process`] asks the
//!   caller's frame source for them; [`Denoiser::denoise_ms1`] and
//!   [`Denoiser::denoise_msms`] refuse when a neighbor would be needed.
//!
//! The order of the stages is an internal detail of this module and may change
//! between releases; callers never sequence stages themselves.

use crate::box_centroid::box_centroid;
use crate::crop::CropGate;
use crate::dia_ms1::DiaMs1Gate;
use crate::dia_window::{filter_per_window, in_window_mask};
use crate::error::{Error, Result};
use crate::filter::filter_iterated;
use crate::frame::FlatFrame;
use crate::halo::horizontal_halo_keep_mask;
use crate::msms::MsmsKeep;
use crate::neighbor::{NeighborIndex, NeighborUsage};
use crate::params::{CropParams, FilterParams, HaloParams, MsmsFilterParams, Stages};
use crate::polygon::PolygonGate;
use crate::smooth::box_average;
use crate::watershed::watershed_centroid;
use crate::windows::{DiaRegions, DiaWindows, FrameMeta, PrmWindows};
use std::sync::atomic::AtomicBool;

/// Upper bound on watershed groups formed per frame — a guard against
/// pathological frames. Real MS1 frames centroid to far fewer than this.
const MAX_CENTROIDS: usize = 100_000;

/// Acquisition scheme of a run, detected from frame types and checked PRM events.
/// Drives the `--preset auto` gate selection (see the CLI): ddaPASEF wants the MS1
/// selection-polygon gate, diaPASEF the isolation-window gates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Acquisition {
    /// ddaPASEF (`MsMsType` 8 present): data-dependent PASEF.
    DdaPasef,
    /// diaPASEF (`MsMsType` 9 present): data-independent PASEF.
    DiaPasef,
    /// prm-PASEF (`MsMsType` 10) with checked target/event metadata.
    PrmPasef,
    /// More than one nonzero frame type. Acquisition gates are disabled.
    Mixed,
    /// Only MS1 frames — no MS/MS in the run.
    Ms1Only,
    /// MS/MS frames present but of an unrecognised `MsMsType`.
    Unknown,
}

impl std::fmt::Display for Acquisition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::DdaPasef => "ddaPASEF",
            Self::DiaPasef => "diaPASEF",
            Self::PrmPasef => "prm-PASEF",
            Self::Mixed => "mixed acquisition",
            Self::Ms1Only => "MS1-only",
            Self::Unknown => "unknown acquisition",
        })
    }
}

/// Check parameter values every front end must reject (neighbor gap, internal
/// gap, centroider choice, halo and gate pads). [`Denoiser::new`] calls it.
pub fn validate(params: &FilterParams, stages: &Stages) -> Result<()> {
    let invalid = |message: &str| Error::InvalidInput(message.into());
    if !stages.neighbors.max_rt_gap_seconds.is_finite()
        || stages.neighbors.max_rt_gap_seconds <= 0.0
    {
        return Err(invalid(
            "neighbor_max_rt_gap must be finite and positive (seconds)",
        ));
    }
    if params.max_internal_gap >= u32::MAX as usize || params.num_iterations > u32::MAX as usize {
        return Err(invalid(
            "max_internal_gap and num_iterations must be below 2^32",
        ));
    }
    if stages.watershed.is_some() && stages.box_centroid.is_some() {
        return Err(invalid("choose only one centroider"));
    }
    if let Some(h) = stages.halo {
        if !h.peak_fraction.is_finite()
            || !(0.0..=1.0).contains(&h.peak_fraction)
            || h.scan_half_width > u32::MAX as usize
        {
            return Err(invalid(
                "halo fraction must be finite in [0,1], with a valid scan width",
            ));
        }
    }
    for (mz, im) in stages
        .ms1_polygon
        .map(|p| (p.mz_pad, p.im_pad))
        .into_iter()
        .chain(stages.dia_ms1.map(|p| (p.mz_pad, p.im_pad)))
    {
        if !mz.is_finite() || !im.is_finite() || mz < 0.0 || im < 0.0 {
            return Err(invalid("gate padding must be finite and nonnegative"));
        }
    }
    Ok(())
}

/// The stages that actually run for `acquisition`. Targeted, mixed and unknown
/// runs may use the MS1 filter, but never discovery geometry or the whole-frame
/// MS/MS fallback. Cropping remains an explicit operation. Errors when the
/// requested MS/MS denoising cannot run on this acquisition.
pub fn effective_stages<'a>(
    acquisition: Acquisition,
    stages: &Stages<'a>,
    crop_only: bool,
) -> Result<Stages<'a>> {
    let mut effective = *stages;
    if matches!(
        acquisition,
        Acquisition::PrmPasef | Acquisition::Mixed | Acquisition::Unknown
    ) {
        if !crop_only && acquisition == Acquisition::PrmPasef && stages.filter_all_frames {
            return Err(Error::InvalidInput(
                "--all-frames (filter_all_frames) is unsupported for prm-PASEF; use --denoise-msms for experimental per-event fragment filtering".into()
            ));
        }
        if !crop_only
            && acquisition != Acquisition::PrmPasef
            && (stages.denoise_msms.is_some() || stages.filter_all_frames)
        {
            return Err(Error::InvalidInput(format!(
                "MS/MS denoising is not supported for {acquisition}; use MS1-only processing (disable denoise_msms and filter_all_frames)"
            )));
        }
        if crop_only || acquisition != Acquisition::PrmPasef {
            effective.denoise_msms = None;
        }
        effective.filter_all_frames = false;
        effective.ms1_polygon = None;
        effective.dia_ms1 = None;
        effective.dia_window = None;
        effective.dda_window = None;
        effective.dia_per_window = false;
    }
    if acquisition == Acquisition::DiaPasef && stages.neighbors.dia_radius > 0 {
        effective.dia_per_window = true;
    }
    Ok(effective)
}

/// A frame after every denoising / crop stage has run, but *before* it is
/// re-encoded. The `dnoise` crate encodes it into a `.d` record; an in-process
/// caller (e.g. a feature finder) can use the surviving points directly.
#[derive(Debug, Clone)]
pub struct DecodedFrame {
    /// Bruker frame `Id`.
    pub frame_id: usize,
    /// True for MS1 frames.
    pub is_ms1: bool,
    /// Scan count for this frame (encode size; also maps scan -> mobility).
    pub num_scans: usize,
    /// Retention time in **seconds** (`Frames.Time`).
    pub rt_seconds: f64,
    /// Surviving points as integer `(scan, tof_idx, intensity)`, after every
    /// enabled stage (vertical filter, halo, gates, smoothing, centroiding).
    pub survivors: Vec<(u32, u32, u32)>,
    /// Input point count before filtering, including RT-cropped frames.
    pub raw_points: u64,
    /// Summed decoded input intensity, including RT-cropped frames.
    pub raw_summed: u64,
    /// Actual temporal evidence used for this frame.
    pub neighbor_usage: NeighborUsage,
    /// True when this frame was emptied by the retention-time crop.
    pub cropped: bool,
}

/// One run's denoising settings and prebuilt state; denoises one frame per call.
/// See the [module docs](self) for which state crosses frames.
///
/// ```
/// use dnoise_core::{Acquisition, Denoiser, FilterParams, FlatFrame, FrameMeta, Stages};
///
/// // One MS1 frame: a real streak at TOF 1000 and an isolated noise point.
/// let meta = vec![FrameMeta { id: 1, num_scans: 700, num_peaks: 9, ms_ms_type: 0, rt: 0.5 }];
/// let params = FilterParams::default();
/// let stages = Stages::default();
/// let denoiser = Denoiser::new(&params, &stages, Acquisition::Ms1Only, meta, false)?;
///
/// let mut scan: Vec<u32> = (10..18).collect();
/// scan.push(300);
/// let mut tof = vec![1000; 8];
/// tof.push(50_000);
/// let frame = FlatFrame { frame_id: 1, num_scans: 700, scan, tof, intensity: vec![100; 9] };
///
/// let out = denoiser.denoise_ms1(0, &frame)?;
/// assert_eq!(out.survivors.len(), 8); // the streak survives, the lone point does not
/// # Ok::<(), dnoise_core::Error>(())
/// ```
#[derive(Debug)]
pub struct Denoiser<'a> {
    params: FilterParams,
    stages: Stages<'a>,
    meta: Vec<FrameMeta>,
    rt_keep: Vec<bool>,
    crop_only: bool,
    msms_keep: Option<MsmsKeep>,
    prm_windows: Option<PrmWindows>,
    dia_msms: Option<&'a MsmsFilterParams>,
    dia_windows: Option<DiaWindows>,
    dia_regions: Option<DiaRegions>,
    dda_windows: Option<DiaWindows>,
    dia_ms1: Option<DiaMs1Gate>,
    polygon: Option<PolygonGate>,
    crop: Option<CropGate>,
    neighbors: Option<NeighborIndex>,
}

impl<'a> Denoiser<'a> {
    /// Validate `params` and `stages` ([`validate`]), apply the acquisition
    /// policy ([`effective_stages`]) and hold `meta` (every frame, ordered by
    /// `Frames.Id`; index `i` below is a position in it). With `crop_only`,
    /// no denoising runs and only the crop applies.
    ///
    /// No optional state is attached yet: set it with the `set_*` methods.
    pub fn new(
        params: &FilterParams,
        stages: &Stages<'a>,
        acquisition: Acquisition,
        meta: Vec<FrameMeta>,
        crop_only: bool,
    ) -> Result<Self> {
        validate(params, stages)?;
        let stages = effective_stages(acquisition, stages, crop_only)?;
        Ok(Self {
            params: *params,
            stages,
            rt_keep: vec![true; meta.len()],
            meta,
            crop_only,
            msms_keep: None,
            prm_windows: None,
            dia_msms: None,
            dia_windows: None,
            dia_regions: None,
            dda_windows: None,
            dia_ms1: None,
            polygon: None,
            crop: None,
            neighbors: None,
        })
    }

    /// The stages that run for this acquisition (after [`effective_stages`]).
    /// Use these, not the requested stages, to decide which state to build.
    pub fn stages(&self) -> &Stages<'a> {
        &self.stages
    }

    /// The run's frame metadata.
    pub fn meta(&self) -> &[FrameMeta] {
        &self.meta
    }

    /// Number of frames.
    pub fn len(&self) -> usize {
        self.meta.len()
    }

    /// True when the run has no frames.
    pub fn is_empty(&self) -> bool {
        self.meta.is_empty()
    }

    /// True when frame `i` is an MS1 frame.
    pub fn is_ms1(&self, i: usize) -> bool {
        self.meta[i].is_ms1()
    }

    /// True in crop-only mode.
    pub fn crop_only(&self) -> bool {
        self.crop_only
    }

    /// ddaPASEF per-precursor keep sets for MS/MS denoising. Ignored (stored as
    /// `None`) when MS/MS denoising does not run for this acquisition
    /// ([`effective_stages`]), as the `dnoise` CLI never builds them then.
    pub fn set_msms_keep(&mut self, keep: Option<MsmsKeep>) {
        self.msms_keep = keep.filter(|_| self.stages.denoise_msms.is_some());
    }

    /// Checked prm-PASEF isolation events. With `denoise_msms` set, each event
    /// is filtered on its own, and a nonempty PRM frame without events is an error.
    pub fn set_prm_windows(&mut self, windows: Option<PrmWindows>) {
        self.prm_windows = windows;
    }

    /// Filter each diaPASEF MS/MS frame (or, with `dia_per_window`, each of its
    /// isolation windows) with the `denoise_msms` parameters. Used when MS/MS
    /// denoising is on and the run has no ddaPASEF isolation events.
    pub fn set_whole_frame_msms(&mut self, on: bool) {
        self.dia_msms = if on { self.stages.denoise_msms } else { None };
    }

    /// diaPASEF isolation windows, for the MS/MS out-of-window gate and
    /// per-window filtering.
    pub fn set_dia_windows(&mut self, windows: Option<DiaWindows>) {
        self.dia_windows = windows;
    }

    /// Checked DIA filter regions (static windows and scanning regions); when
    /// set they take precedence over [`Self::set_dia_windows`] intervals.
    pub fn set_dia_regions(&mut self, regions: Option<DiaRegions>) {
        self.dia_regions = regions;
    }

    /// ddaPASEF isolation-event intervals for the MS/MS out-of-window gate
    /// ([`DiaWindows::from_pasef`]).
    pub fn set_dda_windows(&mut self, windows: Option<DiaWindows>) {
        self.dda_windows = windows;
    }

    /// The diaPASEF MS1 out-of-window gate. Ignored in crop-only mode and when
    /// the gate does not run for this acquisition ([`effective_stages`]).
    pub fn set_dia_ms1_gate(&mut self, gate: Option<DiaMs1Gate>) {
        self.dia_ms1 = gate.filter(|_| self.stages.dia_ms1.is_some() && !self.crop_only);
    }

    /// The MS1 selection-polygon gate. Ignored in crop-only mode and when the
    /// gate does not run for this acquisition ([`effective_stages`]).
    pub fn set_polygon_gate(&mut self, gate: Option<PolygonGate>) {
        self.polygon = gate.filter(|_| self.stages.ms1_polygon.is_some() && !self.crop_only);
    }

    /// The point-level region-of-interest crop (every frame).
    pub fn set_crop_gate(&mut self, gate: Option<CropGate>) {
        self.crop = gate;
    }

    /// Retention-time crop: frames outside `[rt_min, rt_max]` (minutes) are
    /// emitted empty and marked `cropped`. `None`, or bounds without RT, keep all.
    pub fn set_rt_crop(&mut self, crop: Option<&CropParams>) {
        self.rt_keep = match crop {
            Some(cp) if cp.has_rt() => {
                let lo = cp.rt_min.map(|m| m * 60.0).unwrap_or(f64::NEG_INFINITY);
                let hi = cp.rt_max.map(|m| m * 60.0).unwrap_or(f64::INFINITY);
                self.meta.iter().map(|m| m.rt >= lo && m.rt <= hi).collect()
            }
            _ => vec![true; self.meta.len()],
        };
    }

    /// Temporal neighbor support. Ignored in crop-only mode.
    pub fn set_neighbors(&mut self, neighbors: Option<NeighborIndex>) {
        self.neighbors = neighbors;
    }

    /// The ddaPASEF keep sets, if set.
    pub fn msms_keep(&self) -> Option<&MsmsKeep> {
        self.msms_keep.as_ref()
    }

    /// The prm-PASEF events, if set.
    pub fn prm_windows(&self) -> Option<&PrmWindows> {
        self.prm_windows.as_ref()
    }

    /// True when whole-frame diaPASEF MS/MS filtering is on.
    pub fn whole_frame_msms(&self) -> bool {
        self.dia_msms.is_some()
    }

    /// The diaPASEF windows, if set.
    pub fn dia_windows(&self) -> Option<&DiaWindows> {
        self.dia_windows.as_ref()
    }

    /// The DIA regions, if set.
    pub fn dia_regions(&self) -> Option<&DiaRegions> {
        self.dia_regions.as_ref()
    }

    /// The ddaPASEF event intervals, if set.
    pub fn dda_windows(&self) -> Option<&DiaWindows> {
        self.dda_windows.as_ref()
    }

    /// The diaPASEF MS1 gate, if set.
    pub fn dia_ms1_gate(&self) -> Option<&DiaMs1Gate> {
        self.dia_ms1.as_ref()
    }

    /// The polygon gate, if set.
    pub fn polygon_gate(&self) -> Option<&PolygonGate> {
        self.polygon.as_ref()
    }

    /// The neighbor index, if set.
    pub fn neighbors(&self) -> Option<&NeighborIndex> {
        self.neighbors.as_ref()
    }

    /// Denoise frame `i`, reading it (and any temporal neighbors) through
    /// `read_frame(j)`, which returns 0-based frame `j`. Empty frames
    /// (`num_peaks == 0`) are never read. `cancel` is checked while reading
    /// neighbors. This is the call the `dnoise` CLI makes for every frame.
    pub fn process(
        &self,
        i: usize,
        read_frame: &dyn Fn(usize) -> Result<FlatFrame>,
        cancel: Option<&AtomicBool>,
    ) -> Result<DecodedFrame> {
        let meta_i = self.frame_meta(i)?;
        // Empty input frames are not read: a reader cannot decode an absent
        // payload. The `dnoise` writer still writes them as an all-zero scan
        // table so the output stays readable.
        if meta_i.num_peaks == 0 {
            return Ok(DecodedFrame {
                frame_id: meta_i.id,
                is_ms1: meta_i.is_ms1(),
                num_scans: meta_i.num_scans,
                rt_seconds: meta_i.rt,
                survivors: Vec::new(),
                raw_points: 0,
                raw_summed: 0,
                neighbor_usage: NeighborUsage::default(),
                cropped: false,
            });
        }
        let flat = read_frame(i)?;
        self.denoise_frame(i, &flat, read_frame, cancel)
    }

    /// Denoise MS1 frame `i` from its decoded points. Errors when frame `i` is
    /// not MS1, or when temporal neighbors are enabled for it (use
    /// [`Self::process`] with a frame source then).
    pub fn denoise_ms1(&self, i: usize, frame: &FlatFrame) -> Result<DecodedFrame> {
        if !self.frame_meta(i)?.is_ms1() {
            return Err(Error::InvalidInput(format!(
                "frame {} is not an MS1 frame",
                self.meta[i].id
            )));
        }
        self.denoise_frame(i, frame, &no_neighbors, None)
    }

    /// Denoise MS/MS frame `i` from its decoded points. Errors when frame `i` is
    /// MS1, or when temporal neighbors are enabled for it (use [`Self::process`]
    /// with a frame source then).
    pub fn denoise_msms(&self, i: usize, frame: &FlatFrame) -> Result<DecodedFrame> {
        if self.frame_meta(i)?.is_ms1() {
            return Err(Error::InvalidInput(format!(
                "frame {} is not an MS/MS frame",
                self.meta[i].id
            )));
        }
        self.denoise_frame(i, frame, &no_neighbors, None)
    }

    /// `meta[i]`, or an error for an index outside the run.
    fn frame_meta(&self, i: usize) -> Result<&FrameMeta> {
        self.meta.get(i).ok_or_else(|| {
            Error::InvalidInput(format!(
                "frame index {i} is outside the run ({} frames)",
                self.meta.len()
            ))
        })
    }

    /// Check a caller's frame against `meta[i]`: same id and scan count, and
    /// every point inside the scan range and the column arrays equal in length.
    fn check_frame(&self, i: usize, flat: &FlatFrame) -> Result<()> {
        let m = self.frame_meta(i)?;
        let n = flat.scan.len();
        if flat.frame_id != m.id || flat.num_scans != m.num_scans {
            return Err(Error::InvalidInput(format!(
                "frame {} ({} scans) passed for index {i}, which is frame {} ({} scans)",
                flat.frame_id, flat.num_scans, m.id, m.num_scans
            )));
        }
        if flat.tof.len() != n || flat.intensity.len() != n {
            return Err(Error::InvalidInput(format!(
                "frame {}: scan, tof and intensity columns differ in length",
                m.id
            )));
        }
        if flat.scan.iter().any(|&s| s as usize >= m.num_scans) {
            return Err(Error::InvalidInput(format!(
                "frame {}: a point's scan is outside 0..{}",
                m.id, m.num_scans
            )));
        }
        Ok(())
    }

    /// The whole per-frame pipeline, in its one fixed order.
    fn denoise_frame(
        &self,
        i: usize,
        flat: &FlatFrame,
        read_frame: &dyn Fn(usize) -> Result<FlatFrame>,
        cancel: Option<&AtomicBool>,
    ) -> Result<DecodedFrame> {
        self.check_frame(i, flat)?;
        // The stages that act per frame; the polygon/dia_ms1/denoise_msms knobs were
        // already consumed into the gates and keep sets held by `self`.
        let &Stages {
            filter_all_frames,
            halo,
            smooth,
            watershed,
            box_centroid: box_centroid_params,
            dia_window,
            dda_window,
            dia_per_window,
            ..
        } = &self.stages;
        let params = &self.params;
        let prm_msms = self.stages.denoise_msms.zip(self.prm_windows.as_ref());
        let msms = self.msms_keep.as_ref();
        let dia_msms = self.dia_msms;
        let dia_regions = self.dia_regions.as_ref();
        let crop_only = self.crop_only;
        let neighbors = self.neighbors.as_ref();
        let meta_i = &self.meta[i];
        let is_ms1 = meta_i.is_ms1();

        let raw_points = flat.len() as u64;
        let raw_summed: u64 = flat.intensity.iter().map(|&it| it as u64).sum();
        let num_scans = flat.num_scans;
        let frame_id = flat.frame_id;

        // Decode RT-excluded frames too: QC denominators must include removed intensity.
        if !self.rt_keep[i] {
            return Ok(DecodedFrame {
                frame_id: meta_i.id,
                is_ms1,
                num_scans: meta_i.num_scans,
                rt_seconds: meta_i.rt,
                survivors: Vec::new(),
                raw_points,
                raw_summed,
                neighbor_usage: NeighborUsage::default(),
                cropped: true,
            });
        }

        // Sum compatible local observations only to decide which native points survive.
        let neighbor_params = if is_ms1 {
            *params
        } else {
            self.stages
                .denoise_msms
                .map(|p| p.as_filter_params())
                .unwrap_or(*params)
        };
        let neighborhood_mask = match neighbors.filter(|_| !crop_only) {
            Some(n) => n.keep_mask(
                read_frame,
                &self.meta,
                i,
                flat,
                &neighbor_params,
                halo,
                &self.rt_keep,
                cancel,
            )?,
            None => None,
        };
        let (neighborhood_mask, neighbor_usage) = match neighborhood_mask {
            Some((mask, usage)) => (Some(mask), usage),
            None => (None, NeighborUsage::default()),
        };
        let to_filter: &FlatFrame = flat;

        // diaPASEF isolation-window scan intervals for this frame (None for MS1, for
        // ddaPASEF, or when neither DIA feature is enabled). Drives both per-window
        // MS/MS filtering and the out-of-window gate below.
        let dia_iv = dia_regions
            .and_then(|r| r.intervals(meta_i.id))
            .or_else(|| {
                self.dia_windows
                    .as_ref()
                    .and_then(|dw| dw.intervals(meta_i.id))
            });

        // MS1 frames: vertical filter then (optional) horizontal-halo on the survivors.
        // MS/MS frames: pruned by the precursor keep sets when MS/MS denoising is on,
        // otherwise re-encoded unchanged (or vertical-filtered if `filter_all_frames`).
        // In `crop_only` mode no denoising runs at all — every point survives to the
        // crop below, which is the sole filter.
        let mut keep = if crop_only {
            vec![true; to_filter.len()]
        } else if meta_i.is_ms1() {
            match neighborhood_mask {
                // Temporal support already ran streak and halo on the summed
                // neighborhood; only the gates remain.
                Some(mut keep) => {
                    self.apply_ms1_gates(to_filter, &mut keep);
                    keep
                }
                None => self.ms1_keep(to_filter),
            }
        } else if let Some(mask) = neighborhood_mask {
            mask
        } else if let Some((mp, windows)) = prm_msms {
            let intervals = windows.get(&meta_i.id).ok_or_else(|| {
                Error::InvalidInput(format!(
                    "nonempty PRM frame {} has no checked isolation events",
                    meta_i.id
                ))
            })?;
            filter_per_window(to_filter, intervals, &mp.as_filter_params(), halo)
        } else if let Some(mk) = msms {
            mk.keep_mask(to_filter, meta_i.id)
        } else if let Some(mp) = dia_msms {
            // diaPASEF MS/MS: run the same MS/MS filter on each whole frame. With
            // `dia_per_window`, filter each isolation window's scan slice on its own
            // instead, so a mobility run cannot be fused across a window boundary
            // (cross-talk between unrelated isolation events). Temporal support, when
            // enabled, has already supplied the mask above; the `msms_*` knobs apply via
            // FilterParams just like the ddaPASEF path.
            let fp = mp.as_filter_params();
            match dia_iv {
                Some(iv) if dia_per_window => filter_per_window(to_filter, iv, &fp, halo),
                _ => {
                    let mut keep = filter_iterated(to_filter, &fp);
                    if let Some(hp) = halo {
                        apply_halo(to_filter, hp, &mut keep);
                    }
                    keep
                }
            }
        } else if filter_all_frames {
            match dia_iv {
                Some(iv) if dia_per_window => filter_per_window(to_filter, iv, params, halo),
                _ => {
                    let mut keep = filter_iterated(to_filter, params);
                    if let Some(hp) = halo {
                        apply_halo(to_filter, hp, &mut keep);
                    }
                    keep
                }
            }
        } else {
            vec![true; to_filter.len()]
        };

        // Out-of-window gate: drop any MS/MS point whose scan falls outside every
        // isolation window for this frame. Independent of the streak filter, so it
        // also trims mobility-edge noise when no MS/MS filtering runs. (Per-window
        // filtering already excludes these points, making this a no-op there.)
        if !crop_only {
            if let (Some(dp), Some(iv)) = (dia_window, dia_iv) {
                let mask = in_window_mask(&to_filter.scan, iv, dp.scan_pad);
                for (slot, keep_pt) in keep.iter_mut().zip(mask) {
                    *slot &= keep_pt;
                }
            }
            // The ddaPASEF twin of the gate above, driven by PasefFrameMsMsInfo
            // isolation events. Standard timsTOF ddaPASEF files record no
            // out-of-event scans, so this is expected to change nothing there; it
            // enforces the invariant rather than trusting the acquisition.
            let dda_iv = self
                .dda_windows
                .as_ref()
                .and_then(|w| w.intervals(meta_i.id));
            if let (Some(dp), Some(iv)) = (dda_window, dda_iv) {
                let mask = in_window_mask(&to_filter.scan, iv, dp.scan_pad);
                for (slot, keep_pt) in keep.iter_mut().zip(mask) {
                    *slot &= keep_pt;
                }
            }
        }

        // Region-of-interest crop: AND the `(m/z, 1/K0, intensity)` box into the keep
        // mask. Applies to every frame regardless of MS level (a subset of the raw
        // acquisition), and is the only active filter under `crop_only`.
        if let Some(cg) = &self.crop {
            cg.apply(
                &to_filter.scan,
                &to_filter.tof,
                &to_filter.intensity,
                &mut keep,
            );
        }

        let survivors = to_filter.survivors(&keep);

        // Optional intensity smoothing then watershed centroiding, both applied only
        // to frames the vertical filter actually processed — MS1 always, MS/MS only
        // under `filter_all_frames` (and never on the separate per-precursor MS/MS-
        // denoise path, which has its own keep logic). Smoothing runs first so the
        // watershed seeds on the stabilised intensities.
        let filtered_here = !crop_only
            && (meta_i.is_ms1() || dia_msms.is_some() || (msms.is_none() && filter_all_frames));
        // PRM postprocessing stays inside each event as well: smoothing and
        // centroiding must not mix adjacent targets after the keep mask is built.
        let finish = |points: Vec<(u32, u32, u32)>| {
            let points = match smooth {
                Some(sp) => box_average(&points, num_scans, sp),
                None => points,
            };
            let points = match watershed {
                Some(wp) => watershed_centroid(&points, wp, MAX_CENTROIDS),
                None => points,
            };
            match box_centroid_params {
                Some(bp) => box_centroid(&points, bp),
                None => points,
            }
        };
        let postprocess_intervals = if !crop_only && !is_ms1 {
            neighbors
                .map(|n| n.intervals(i))
                .filter(|iv| !iv.is_empty())
                .or_else(|| {
                    dia_regions
                        .and_then(|r| r.intervals(meta_i.id))
                        .map(|iv| iv.to_vec())
                })
                .or_else(|| prm_msms.map(|(_, w)| w[&meta_i.id].clone()))
        } else {
            None
        };
        let survivors = if let Some(intervals) = postprocess_intervals
            .filter(|_| smooth.is_some() || watershed.is_some() || box_centroid_params.is_some())
        {
            intervals
                .iter()
                .flat_map(|&(begin, end)| {
                    finish(
                        survivors
                            .iter()
                            .copied()
                            .filter(|p| p.0 >= begin && p.0 < end)
                            .collect(),
                    )
                })
                .collect()
        } else if filtered_here {
            finish(survivors)
        } else {
            survivors
        };

        Ok(DecodedFrame {
            frame_id,
            is_ms1,
            num_scans,
            rt_seconds: meta_i.rt,
            survivors,
            raw_points,
            raw_summed,
            neighbor_usage,
            cropped: false,
        })
    }

    /// MS1 keep mask without temporal support, in stage order: m/z pre-cut,
    /// streak filter, MS1 gates, halo.
    ///
    /// The pre-cut runs only with exactly one MS1 gate. It drops points whose TOF
    /// is more than `(iterations + 1) x mz_half_width` outside the gate's TOF
    /// hull. The streak filter decides a point from points within
    /// `iterations x mz_half_width` TOF of it, so every point within one
    /// `mz_half_width` of the hull is decided exactly as on the whole frame. A
    /// point gate keeps nothing outside the hull, so there the cut is exact. A
    /// feature-level gate can keep a feature that runs out of the hull; any such
    /// feature leaves the hull through a point within `mz_half_width` of it, which
    /// the cut decides exactly, so it shows up as a kept point outside the hull.
    /// The frame is then redone without the cut. Output is therefore identical to
    /// running the stages on the whole frame.
    fn ms1_keep(&self, frame: &FlatFrame) -> Vec<bool> {
        let params = &self.params;
        let span = match (&self.dia_ms1, &self.polygon) {
            (Some(g), None) => g.tof_span(),
            (None, Some(g)) => g.tof_span(),
            _ => None,
        };
        if let Some((lo, hi)) = span {
            let margin = (params.num_iterations as u32)
                .saturating_add(1)
                .saturating_mul(params.mz_half_width);
            let (cut_lo, cut_hi) = (lo.saturating_sub(margin), hi.saturating_add(margin));
            let idx: Vec<usize> = (0..frame.len())
                .filter(|&k| (cut_lo..=cut_hi).contains(&frame.tof[k]))
                .collect();
            if idx.len() < frame.len() {
                let sub = FlatFrame {
                    frame_id: frame.frame_id,
                    num_scans: frame.num_scans,
                    scan: idx.iter().map(|&k| frame.scan[k]).collect(),
                    tof: idx.iter().map(|&k| frame.tof[k]).collect(),
                    intensity: idx.iter().map(|&k| frame.intensity[k]).collect(),
                };
                let mut sub_keep = filter_iterated(&sub, params);
                self.apply_ms1_gates(&sub, &mut sub_keep);
                let escaped = sub_keep
                    .iter()
                    .zip(&sub.tof)
                    .any(|(&kept, t)| kept && !(lo..=hi).contains(t));
                if !escaped {
                    if let Some(hp) = self.stages.halo {
                        apply_halo(&sub, hp, &mut sub_keep);
                    }
                    let mut keep = vec![false; frame.len()];
                    for (&k, kept) in idx.iter().zip(sub_keep) {
                        keep[k] = kept;
                    }
                    return keep;
                }
            }
        }
        let mut keep = filter_iterated(frame, params);
        self.apply_ms1_gates(frame, &mut keep);
        if let Some(hp) = self.stages.halo {
            apply_halo(frame, hp, &mut keep);
        }
        keep
    }

    /// AND the MS1 acquisition gates into `keep`: the diaPASEF window gate, then
    /// the selection-polygon gate. With `overlap`, a gate decides per streak
    /// feature of the current survivors ([`crate::overlap`]).
    fn apply_ms1_gates(&self, frame: &FlatFrame, keep: &mut [bool]) {
        let gates = [
            self.dia_ms1
                .as_ref()
                .map(|g| (g.keep_mask(&frame.scan, &frame.tof), g.overlap)),
            self.polygon
                .as_ref()
                .map(|g| (g.keep_mask(&frame.scan, &frame.tof), g.overlap)),
        ];
        for (mut mask, overlap) in gates.into_iter().flatten() {
            if overlap {
                mask = overlap_mask(frame, keep, &mask, &self.params);
            }
            for (slot, inside) in keep.iter_mut().zip(mask) {
                *slot &= inside;
            }
        }
    }
}

/// The frame source of [`Denoiser::denoise_ms1`] / [`Denoiser::denoise_msms`]:
/// temporal neighbors need other frames, which those calls do not have.
fn no_neighbors(index: usize) -> Result<FlatFrame> {
    Err(Error::InvalidInput(format!(
        "temporal neighbor support needs frame {index}; use Denoiser::process with a frame source"
    )))
}

/// Feature-level gate mask: link surviving points with the streak filter's own
/// adjacency (column half-width, bridged gap + 1) and keep features that touch
/// the gate anywhere, over their whole extent. See [`crate::overlap`].
fn overlap_mask(
    frame: &FlatFrame,
    keep: &[bool],
    inside: &[bool],
    params: &FilterParams,
) -> Vec<bool> {
    crate::overlap::extend_to_features(
        &frame.scan,
        &frame.tof,
        keep,
        inside,
        params.mz_half_width,
        params.max_internal_gap as u32 + 1,
    )
}

/// Run the horizontal-halo filter on the currently-kept points of `frame` and
/// turn off `keep` for any the filter removes. Operates in integer
/// `(scan, TOF index)` space — no calibration needed.
fn apply_halo(frame: &FlatFrame, hp: &HaloParams, keep: &mut [bool]) {
    let idx: Vec<usize> = (0..frame.len()).filter(|&i| keep[i]).collect();
    if idx.is_empty() {
        return;
    }
    let scan: Vec<u32> = idx.iter().map(|&i| frame.scan[i]).collect();
    let tof: Vec<u32> = idx.iter().map(|&i| frame.tof[i]).collect();
    let inten: Vec<u32> = idx.iter().map(|&i| frame.intensity[i]).collect();

    let hmask = horizontal_halo_keep_mask(&scan, &tof, &inten, frame.num_scans, hp);
    for (k, &i) in idx.iter().enumerate() {
        if !hmask[k] {
            keep[i] = false;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dia_ms1::TofScanBox;
    use crate::params::DiaMs1WindowParams;

    fn meta(ms_ms_type: i64) -> Vec<FrameMeta> {
        vec![FrameMeta {
            id: 1,
            num_scans: 700,
            num_peaks: 1,
            ms_ms_type,
            rt: 1.0,
        }]
    }

    fn lone_point() -> FlatFrame {
        FlatFrame {
            frame_id: 1,
            num_scans: 700,
            scan: vec![5],
            tof: vec![5],
            intensity: vec![7],
        }
    }

    #[test]
    fn ms_level_is_checked() {
        let p = FilterParams::default();
        let s = Stages::default();
        let d = Denoiser::new(&p, &s, Acquisition::DdaPasef, meta(8), false).unwrap();
        assert!(d.denoise_ms1(0, &lone_point()).is_err());
        // No MS/MS stage enabled: the frame passes through unchanged.
        assert_eq!(
            d.denoise_msms(0, &lone_point()).unwrap().survivors,
            vec![(5, 5, 7)]
        );
    }

    #[test]
    fn empty_frames_are_not_read() {
        let p = FilterParams::default();
        let s = Stages::default();
        let mut m = meta(0);
        m[0].num_peaks = 0;
        let d = Denoiser::new(&p, &s, Acquisition::Ms1Only, m, false).unwrap();
        let out = d
            .process(0, &|_| panic!("empty frame must not be read"), None)
            .unwrap();
        assert!(out.survivors.is_empty() && !out.cropped);
    }

    #[test]
    fn rt_crop_empties_frames_but_counts_their_input() {
        let p = FilterParams::default();
        let s = Stages::default();
        let mut d = Denoiser::new(&p, &s, Acquisition::Ms1Only, meta(0), false).unwrap();
        d.set_rt_crop(Some(&CropParams {
            rt_min: Some(1.0),
            ..CropParams::default()
        }));
        let out = d.denoise_ms1(0, &lone_point()).unwrap();
        assert!(out.cropped && out.survivors.is_empty());
        assert_eq!((out.raw_points, out.raw_summed), (1, 7));
    }

    /// Streaks and scattered noise across TOF 0..2000, plus a staircase of
    /// linked streaks running from inside a TOF 800..=1200 gate out to TOF 1300.
    fn gated_frame() -> FlatFrame {
        let (mut scan, mut tof, mut intensity) = (Vec::new(), Vec::new(), Vec::new());
        let mut push = |s: u32, t: u32, i: u32| {
            scan.push(s);
            tof.push(t);
            intensity.push(i);
        };
        let mut x: u64 = 12345;
        let mut next = move || {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (x >> 33) as u32
        };
        for _ in 0..40 {
            let (t, s0, len) = (next() % 2000, next() % 90, 3 + next() % 10);
            for s in s0..(s0 + len).min(100) {
                push(s, t, 50 + next() % 500);
                if next() % 3 == 0 {
                    push(s, t + 1 + next() % 60, 5 + next() % 40); // halo-like flank
                }
            }
        }
        for _ in 0..400 {
            push(next() % 100, next() % 2000, 1 + next() % 100);
        }
        for step in 0..40 {
            let t = 1180 + 3 * step;
            for s in 30..40 {
                push(s, t, 300);
            }
        }
        FlatFrame {
            frame_id: 1,
            num_scans: 100,
            scan,
            tof,
            intensity,
        }
    }

    #[test]
    fn bad_frames_are_errors_not_panics() {
        let d = Denoiser::new(
            &FilterParams::default(),
            &Stages::default(),
            Acquisition::Ms1Only,
            meta(0),
            false,
        )
        .unwrap();
        assert!(d.denoise_ms1(0, &lone_point()).is_ok());
        assert!(d.denoise_ms1(1, &lone_point()).is_err());
        let mut f = lone_point();
        f.scan = vec![700];
        assert!(d.denoise_ms1(0, &f).is_err());
        let mut f = lone_point();
        f.frame_id = 2;
        assert!(d.denoise_ms1(0, &f).is_err());
        let mut f = lone_point();
        f.tof.push(9);
        assert!(d.denoise_ms1(0, &f).is_err());
    }

    #[test]
    fn precut_matches_the_whole_frame() {
        let p = FilterParams::default();
        let halo = HaloParams::default();
        let dia_ms1 = DiaMs1WindowParams::default();
        let s = Stages {
            halo: Some(&halo),
            dia_ms1: Some(&dia_ms1),
            ..Stages::default()
        };
        let frame = gated_frame();
        for overlap in [false, true] {
            let mut d = Denoiser::new(&p, &s, Acquisition::DiaPasef, meta(0), false).unwrap();
            let boxes = [TofScanBox {
                scan_lo: 20,
                scan_hi: 80,
                tof_lo: 800,
                tof_hi: 1200,
            }];
            let mut gate = DiaMs1Gate::build(&boxes, 100).unwrap();
            gate.overlap = overlap;
            d.set_dia_ms1_gate(Some(gate));

            let mut want = filter_iterated(&frame, &p);
            d.apply_ms1_gates(&frame, &mut want);
            apply_halo(&frame, &halo, &mut want);
            assert_eq!(d.ms1_keep(&frame), want, "overlap = {overlap}");
            // The staircase escapes past the cut only in overlap mode, which is
            // the case the whole-frame fallback exists for.
            let escaped = (0..frame.len()).any(|k| want[k] && frame.tof[k] > 1200 + 9);
            assert_eq!(escaped, overlap);
        }
    }

    #[test]
    fn invalid_parameters_are_refused() {
        let p = FilterParams {
            max_internal_gap: usize::MAX,
            ..FilterParams::default()
        };
        let s = Stages::default();
        let e = Denoiser::new(&p, &s, Acquisition::Ms1Only, meta(0), false).unwrap_err();
        assert_eq!(
            e.to_string(),
            "max_internal_gap and num_iterations must be below 2^32"
        );
    }
}
