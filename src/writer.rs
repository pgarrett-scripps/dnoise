//! Orchestration: copy the source `.d`, rewrite `analysis.tdf_bin` with filtered
//! frames (re-encoded as type 2), and fix up the `analysis.tdf` SQLite database.

use crate::box_centroid::box_centroid;
use crate::codec::encode_frame_type2;
use crate::crop::CropGate;
use crate::dia_ms1::{DiaMs1Gate, TofScanBox};
use crate::dia_window::{filter_per_window, in_window_mask};
use crate::error::{DnoiseError, Result};
use crate::filter::filter_iterated;
use crate::frame::FlatFrame;
use crate::halo::horizontal_halo_keep_mask;
use crate::msms::{MsmsKeep, build_msms_keep};
use crate::params::{
    CropParams, DiaMs1WindowParams, FilterParams, HaloParams, Ms1PolygonParams, MsmsFilterParams,
    Stages,
};
use crate::polygon::PolygonGate;
use crate::smooth::box_average;
use crate::tdf::{self, DiaWindows, FrameUpdate};
use crate::watershed::watershed_centroid;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use timsrust::converters::ConvertableDomain;
use timsrust::readers::{FrameReader, MetadataReader};
use tracing::{debug, info, warn};

/// Frames are read+filtered+encoded in parallel batches of this size, then the
/// batch is written sequentially (so offsets stay ordered) before the next.
/// Bounds peak memory to roughly this many encoded frames.
const CHUNK: usize = 2048;

/// Upper bound on watershed groups formed per frame — a guard against
/// pathological frames. Real MS1 frames centroid to far fewer than this.
const MAX_CENTROIDS: usize = 100_000;

/// Summary returned by [`denoise`].
#[derive(Debug, Clone, Copy, Default)]
#[non_exhaustive]
pub struct DenoiseStats {
    /// Total frames in the run (MS1 + MS/MS + empty).
    pub frames: usize,
    /// MS1 frames.
    pub ms1_frames: usize,
    /// MS/MS frames.
    pub msms_frames: usize,
    /// Frames emptied by the retention-time crop (subset of `frames`).
    pub cropped_frames: usize,
    /// Frames actually processed. Equals `frames` unless a dry-run `sample` was
    /// requested, in which case it is the sampled subset.
    pub processed_frames: usize,
    /// Total input points across all processed frames.
    pub raw_points: u64,
    /// Total points kept after filtering + crop.
    pub kept_points: u64,
    /// Input points in MS1 frames only.
    pub raw_ms1_points: u64,
    /// Kept points in MS1 frames only.
    pub kept_ms1_points: u64,
    /// Summed intensity of all input points (processed frames).
    pub raw_summed_intensity: u64,
    /// Summed intensity of all kept points.
    pub kept_summed_intensity: u64,
    /// True when this was a dry run (no output written).
    pub dry_run: bool,
}

/// Progress update passed to the callback of [`denoise_with_progress`].
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct Progress {
    /// Frames written so far.
    pub frames_done: usize,
    /// Total frames to process.
    pub frames_total: usize,
}

/// Dry-run frame sampling: process only a pseudo-random subset of frames to
/// estimate the data reduction quickly, without touching output. Selection is
/// deterministic in `seed`, so a run is reproducible and comparable across
/// parameter sweeps.
#[derive(Debug, Clone, Copy)]
pub struct SampleSpec {
    /// Fraction of frames to process, in `(0, 1]`.
    pub fraction: f64,
    /// Seed for the deterministic frame selector.
    pub seed: u64,
}

/// Run-level options orthogonal to the filter itself: overwrite behaviour, the
/// region-of-interest crop, crop-only mode, dry-run / sampling, and an optional
/// cancellation token. Bundled so the `denoise*` entry points stay to a few
/// arguments.
#[derive(Debug, Clone, Copy, Default)]
pub struct RunOptions<'a> {
    /// Overwrite an existing `output` folder.
    pub force: bool,
    /// Compute statistics without writing any output `.d`.
    pub dry_run: bool,
    /// Region-of-interest crop applied to every frame (`None` = no crop).
    pub crop: Option<&'a CropParams>,
    /// Skip all denoising (vertical filter, halo, gates, centroiders) and only
    /// apply the crop — carve a subset `.d` without altering retained signal.
    /// Requires `crop` to be set; ignored otherwise.
    pub crop_only: bool,
    /// Dry-run frame sampling for a fast reduction estimate (`None` = all frames).
    /// Only honoured together with `dry_run`.
    pub sample: Option<SampleSpec>,
    /// Cooperative cancellation token. When set and flipped to `true`, the run
    /// stops at the next frame-chunk boundary and returns [`DnoiseError::Cancelled`];
    /// the caller should discard any partial output. `None` = never cancelled.
    pub cancel: Option<&'a AtomicBool>,
}

/// Deterministic per-frame selector for dry-run sampling: hash `(seed, index)`
/// with SplitMix64 and keep the frame when the hash falls below `fraction` of the
/// u64 range. Order-independent, so any frame subset is reproducible from the seed.
fn frame_sampled(index: usize, seed: u64, fraction: f64) -> bool {
    let mut z = seed
        .wrapping_add(index as u64)
        .wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    let threshold = (fraction.clamp(0.0, 1.0) * u64::MAX as f64) as u64;
    z <= threshold
}

/// Denoise `input` (.d) into a new `output` (.d).
///
/// The core vertical-IM filter ([`FilterParams`]) runs on MS1 frames; `stages`
/// selects every optional stage layered on top (halo, MS/MS denoising, smoothing,
/// centroiding, the diaPASEF window gates, and the ddaPASEF selection-polygon
/// gate) — see [`Stages`] for the per-stage semantics. `force` overwrites an
/// existing `output`.
///
/// This reports no progress; use [`denoise_with_progress`] to receive
/// [`Progress`] updates as frames are written. For the crop, crop-only, and
/// dry-run options use [`denoise_with_options`].
pub fn denoise(
    input: &Path,
    output: &Path,
    params: &FilterParams,
    stages: &Stages,
    force: bool,
) -> Result<DenoiseStats> {
    let opts = RunOptions {
        force,
        ..RunOptions::default()
    };
    denoise_with_options(input, output, params, stages, &opts, |_| {})
}

/// Like [`denoise`] but takes a full [`RunOptions`] (crop, crop-only, dry-run,
/// sampling) and a progress callback. This is the most general entry point; the
/// others are thin wrappers over it.
pub fn denoise_with_options<F: FnMut(Progress)>(
    input: &Path,
    output: &Path,
    params: &FilterParams,
    stages: &Stages,
    options: &RunOptions,
    progress: F,
) -> Result<DenoiseStats> {
    run(input, output, params, stages, options, progress)
}

/// Like [`denoise`], but invokes `progress` once before processing and again
/// after each frame is written, so callers (e.g. a CLI) can drive a progress bar
/// without the library depending on any UI crate.
pub fn denoise_with_progress<F: FnMut(Progress)>(
    input: &Path,
    output: &Path,
    params: &FilterParams,
    stages: &Stages,
    force: bool,
    progress: F,
) -> Result<DenoiseStats> {
    let opts = RunOptions {
        force,
        ..RunOptions::default()
    };
    run(input, output, params, stages, &opts, progress)
}

/// The full pipeline behind every public entry point. Reads the input `.d`, builds
/// the per-run gates (including the crop), filters + crops each frame in parallel
/// chunks, and — unless `options.dry_run` — writes the rewritten `analysis.tdf_bin`
/// and fixes up `analysis.tdf`.
fn run<F: FnMut(Progress)>(
    input: &Path,
    output: &Path,
    params: &FilterParams,
    stages: &Stages,
    options: &RunOptions,
    mut progress: F,
) -> Result<DenoiseStats> {
    let &RunOptions {
        force,
        dry_run,
        crop,
        crop_only,
        sample,
        cancel,
    } = options;
    // Unpack the stages this function builds gates from; the per-frame stages
    // (smoothing, centroiding, etc.) are forwarded to `process_frame` via `stages`.
    let &Stages {
        halo,
        denoise_msms,
        dia_window,
        dia_per_window,
        dia_ms1,
        ms1_polygon,
        ..
    } = stages;

    let in_tdf = input.join("analysis.tdf");
    let in_bin = input.join("analysis.tdf_bin");
    if !in_tdf.is_file() || !in_bin.is_file() {
        return Err(DnoiseError::NotADotD(input.to_path_buf()));
    }
    // A dry run never writes, so the output folder is left completely untouched.
    if !dry_run {
        if output.exists() {
            if force {
                fs::remove_dir_all(output)?;
            } else {
                return Err(DnoiseError::OutputExists(output.to_path_buf()));
            }
        }
        // Copy everything except the binary (we regenerate that), so
        // calibration.sqlite, analysis.tdf, etc. come along.
        copy_dir_except(input, output, "analysis.tdf_bin")?;
    }

    let reader = FrameReader::new(input).map_err(|e| DnoiseError::OpenFrames(e.to_string()))?;
    let n_frames = reader.len();
    // Frame metadata (ordered by Id == timsrust index). Empty frames are handled
    // without timsrust, which cannot decode their absent payload.
    let meta = tdf::read_frame_meta(&in_tdf)?;

    // Frame inventory + acquisition scheme, logged up front so a caller (human or
    // agent) can see what kind of run this is and how much work it entails before
    // any frames are written. Scheme is read straight off `MsMsType` (8 = ddaPASEF,
    // 9 = diaPASEF) — no extra table reads.
    let n_ms1 = meta.iter().filter(|m| m.is_ms1()).count();
    let n_empty = meta.iter().filter(|m| m.num_peaks == 0).count();
    let scheme = if meta.iter().any(|m| m.ms_ms_type == 9) {
        "diaPASEF"
    } else if meta.iter().any(|m| m.ms_ms_type == 8) {
        "ddaPASEF"
    } else {
        "MS1-only/unknown"
    };
    info!(input = %input.display(), output = %output.display(), "denoise: starting");
    info!(
        scheme,
        frames = n_frames,
        ms1 = n_ms1,
        msms = n_frames - n_ms1,
        empty = n_empty,
        "denoise: frame inventory"
    );

    // MS1 subsequence used by the running-average pre-filter: `ms1_indices` maps a
    // position in the MS1-only stream to a global frame index; `ms1_pos` is the
    // reverse (None for MS/MS frames). The window skips interleaved MS/MS frames.
    let ms1_indices: Vec<usize> = (0..n_frames).filter(|&i| meta[i].is_ms1()).collect();
    let mut ms1_pos: Vec<Option<usize>> = vec![None; n_frames];
    for (p, &gi) in ms1_indices.iter().enumerate() {
        ms1_pos[gi] = Some(p);
    }

    // MS/MS denoising splits by acquisition scheme, both driven by the same
    // `denoise_msms` params:
    //   * ddaPASEF — each precursor is re-isolated across several frames, so we
    //     build per-precursor keep sets up front (PasefFrameMsMsInfo) and combine
    //     a precursor's fragment scans across frames before filtering.
    //   * diaPASEF — has no PasefFrameMsMsInfo (each isolation window is sampled
    //     once per cycle, nothing to combine). We detect this by an empty
    //     ddaPASEF window table and instead run the same MS/MS filter on each
    //     whole MS/MS frame as-is (the `dia_msms` path in `process_frame`).
    let (msms_keep, dia_msms) = match denoise_msms {
        Some(mp) => {
            let windows = tdf::read_pasef_msms(&in_tdf)?;
            if windows.is_empty() {
                info!("MS/MS denoise: diaPASEF whole-frame path (no PasefFrameMsMsInfo)");
                (None, Some(mp))
            } else {
                let keep = build_msms_keep(&reader, &meta, &windows, mp, halo)?;
                info!(
                    isolation_events = windows.len(),
                    "MS/MS denoise: ddaPASEF per-precursor path"
                );
                (Some(keep), None)
            }
        }
        None => (None, None),
    };
    let msms_ref = msms_keep.as_ref();

    // diaPASEF isolation windows, read once and shared. Needed by either the
    // out-of-window gate (`dia_window`) or per-window MS/MS filtering
    // (`dia_per_window`). Empty for ddaPASEF, so both features no-op there.
    let dia_windows: Option<DiaWindows> = if dia_window.is_some() || dia_per_window {
        let w = tdf::read_dia_windows(&in_tdf)?;
        if w.is_empty() {
            debug!("diaPASEF window feature requested but no windows found (ddaPASEF?) — skipped");
            None
        } else {
            info!("diaPASEF isolation-window scheme loaded");
            Some(w)
        }
    } else {
        None
    };
    let dia_windows_ref = dia_windows.as_ref();

    // diaPASEF MS1 out-of-window gate: build the padded `(scan, TOF)` lookup once
    // from the isolation windows + calibration. `None` for ddaPASEF (no windows).
    let dia_ms1_gate = match dia_ms1 {
        Some(mp) => build_dia_ms1_gate(&in_tdf, mp, &meta)?,
        None => None,
    };
    if dia_ms1.is_some() {
        match &dia_ms1_gate {
            Some(_) => info!("diaPASEF MS1 out-of-window gate active"),
            None => debug!("diaPASEF MS1 gate requested but no isolation windows — skipped"),
        }
    }
    let dia_ms1_ref = dia_ms1_gate.as_ref();

    // MS1 selection-polygon gate: build the per-scan TOF lookup once from the
    // run's IMS PolygonFilter + calibration. `None` when the run stores no polygon.
    let polygon_gate = match ms1_polygon {
        Some(pp) => build_polygon_gate(&in_tdf, pp, &meta)?,
        None => None,
    };
    if ms1_polygon.is_some() {
        match &polygon_gate {
            Some(_) => info!("MS1 selection-polygon gate active"),
            None => {
                debug!("MS1 polygon gate requested but run stores no usable polygon — skipped")
            }
        }
    }
    let polygon_ref = polygon_gate.as_ref();

    // Region-of-interest crop: convert the physical `(m/z, 1/K0)` bounds to integer
    // `(TOF, scan)` once via the run calibration (RT bounds are applied per frame
    // below). Applies to every frame — this is a subset of the acquisition, not a
    // signal/noise decision. `None` when no crop is requested or it is RT-only.
    let crop_gate = match crop {
        Some(cp) if !cp.is_empty() => {
            let md =
                MetadataReader::new(&in_tdf).map_err(|e| DnoiseError::Metadata(e.to_string()))?;
            let num_scans = meta.iter().map(|m| m.num_scans).max().unwrap_or(0);
            let g = CropGate::build(
                cp,
                num_scans,
                |mz| md.mz_converter.invert(mz),
                |k0| md.im_converter.invert(k0),
            );
            info!(
                point_crop = g.is_active(),
                rt_crop = cp.has_rt(),
                crop_only,
                "crop: region-of-interest gate built"
            );
            g.is_active().then_some(g)
        }
        _ => None,
    };
    let crop_ref = crop_gate.as_ref();

    // Per-frame retention-time keep mask (crop bounds are in minutes; `Frames.Time`
    // is in seconds). Frames outside the window are emitted empty rather than
    // deleted, so the frame axis stays valid. All-true when no RT bound is set.
    let rt_keep: Vec<bool> = match crop {
        Some(cp) if cp.has_rt() => {
            let lo = cp.rt_min.map(|m| m * 60.0).unwrap_or(f64::NEG_INFINITY);
            let hi = cp.rt_max.map(|m| m * 60.0).unwrap_or(f64::INFINITY);
            meta.iter().map(|m| m.rt >= lo && m.rt <= hi).collect()
        }
        _ => vec![true; n_frames],
    };

    // Frames to process. In a dry run with `sample` set, this is a deterministic
    // pseudo-random subset (for a fast reduction estimate); otherwise every frame,
    // in order (so the sequential offsets written below stay consistent).
    let selected: Vec<usize> = match sample {
        Some(s) if dry_run => (0..n_frames)
            .filter(|&i| frame_sampled(i, s.seed, s.fraction))
            .collect(),
        _ => (0..n_frames).collect(),
    };
    match sample {
        Some(_) if !dry_run => {
            warn!("--sample ignored without --dry-run (a real run must process every frame)")
        }
        Some(s) => info!(
            sampled = selected.len(),
            total = n_frames,
            fraction = s.fraction,
            "dry-run: processing a frame sample"
        ),
        None => {}
    }
    let n_process = selected.len();

    // Preserve the leading header that precedes the first frame (Bruker reserves a
    // block at the start of the .tdf_bin), and start writing frames after it so the
    // new TimsId offsets land in the same layout the Bruker reader expects. A dry
    // run opens no output file.
    let header_len = tdf::binary_header_len(&in_tdf)?;
    let mut bin = if dry_run {
        None
    } else {
        let mut b = BufWriter::new(fs::File::create(output.join("analysis.tdf_bin"))?);
        if header_len > 0 {
            let mut header = vec![0u8; header_len as usize];
            fs::File::open(&in_bin).and_then(|mut f| f.read_exact(&mut header))?;
            b.write_all(&header)?;
        }
        Some(b)
    };

    progress(Progress {
        frames_done: 0,
        frames_total: n_process,
    });

    // Per-run context shared by every frame: the prebuilt MS/MS keep sets and
    // gates derived above, the crop, plus the MS1-stream index maps. Bundling these
    // keeps `process_frame` to a handful of arguments.
    let ctx = FrameCtx {
        msms: msms_ref,
        dia_msms,
        dia_windows: dia_windows_ref,
        dia_ms1: dia_ms1_ref,
        polygon: polygon_ref,
        crop: crop_ref,
        crop_only,
        rt_keep: &rt_keep,
        ms1_indices: &ms1_indices,
        ms1_pos: &ms1_pos,
    };

    let mut offset: u64 = header_len;
    let mut updates: Vec<FrameUpdate> = Vec::with_capacity(n_process);
    let mut raw_points: u64 = 0;
    let mut kept_points: u64 = 0;
    let mut raw_ms1: u64 = 0;
    let mut kept_ms1: u64 = 0;
    let mut raw_summed: u64 = 0;
    let mut kept_summed: u64 = 0;
    let mut cropped_frames: usize = 0;
    let mut frames_done: usize = 0;

    for chunk in selected.chunks(CHUNK) {
        // Cooperative cancellation: check once per chunk (a real run's partial
        // output is incomplete, so the caller discards it on Cancelled).
        if let Some(c) = cancel {
            if c.load(Ordering::Relaxed) {
                return Err(DnoiseError::Cancelled);
            }
        }
        let processed: Vec<ProcessedFrame> = chunk
            .par_iter()
            .map(|&i| process_frame(&reader, &meta, i, params, stages, &ctx))
            .collect::<Result<_>>()?;

        for pf in processed {
            raw_points += pf.raw_points;
            kept_points += pf.num_peaks;
            raw_summed += pf.raw_summed;
            kept_summed += pf.summed_intensities;
            if pf.is_ms1 {
                raw_ms1 += pf.raw_points;
                kept_ms1 += pf.num_peaks;
            }
            if pf.cropped {
                cropped_frames += 1;
            }
            if let Some(b) = bin.as_mut() {
                b.write_all(&pf.record)?;
                updates.push(FrameUpdate {
                    frame_id: pf.frame_id,
                    tims_id: offset,
                    num_peaks: pf.num_peaks,
                    max_intensity: pf.max_intensity,
                    summed_intensities: pf.summed_intensities,
                });
                offset += pf.record.len() as u64;
            }
            frames_done += 1;
            progress(Progress {
                frames_done,
                frames_total: n_process,
            });
        }
    }
    if let Some(b) = bin.as_mut() {
        b.flush()?;
    }
    drop(bin);

    // Only a real run rewrites the database (offsets, peak counts, compression type).
    if !dry_run {
        tdf::update_metadata(&output.join("analysis.tdf"), &updates)?;
    }

    let kept_pct = if raw_points > 0 {
        // Round to 2 decimals so the log field is readable (e.g. 36.64, not
        // 36.635067948941554).
        ((10_000.0 * kept_points as f64 / raw_points as f64).round()) / 100.0
    } else {
        0.0
    };
    info!(
        dry_run,
        processed_frames = n_process,
        raw_points,
        kept_points,
        kept_pct,
        "denoise: complete"
    );

    Ok(DenoiseStats {
        frames: n_frames,
        ms1_frames: n_ms1,
        msms_frames: n_frames - n_ms1,
        cropped_frames,
        processed_frames: n_process,
        raw_points,
        kept_points,
        raw_ms1_points: raw_ms1,
        kept_ms1_points: kept_ms1,
        raw_summed_intensity: raw_summed,
        kept_summed_intensity: kept_summed,
        dry_run,
    })
}

struct ProcessedFrame {
    frame_id: usize,
    record: Vec<u8>,
    raw_points: u64,
    num_peaks: u64,
    max_intensity: u32,
    summed_intensities: u64,
    /// Summed intensity of the frame's input points (before filtering/crop).
    raw_summed: u64,
    /// Whether this is an MS1 frame (for the per-level stat split).
    is_ms1: bool,
    /// Whether this frame was emptied by the retention-time crop.
    cropped: bool,
}

/// Per-run context for [`process_frame`]: the MS/MS keep sets and gates built
/// once in [`run`], plus the crop and the MS1-stream index maps. Lets the
/// per-frame worker take the run's derived state as a single value.
struct FrameCtx<'a> {
    /// ddaPASEF per-precursor keep sets (`None` unless MS/MS denoising on ddaPASEF).
    msms: Option<&'a MsmsKeep>,
    /// diaPASEF MS/MS filter knobs (`None` unless MS/MS denoising on diaPASEF).
    dia_msms: Option<&'a MsmsFilterParams>,
    /// diaPASEF isolation windows (`None` for ddaPASEF or when unused).
    dia_windows: Option<&'a DiaWindows>,
    /// Built diaPASEF MS1 out-of-window gate (`None` when disabled / ddaPASEF).
    dia_ms1: Option<&'a DiaMs1Gate>,
    /// Built MS1 selection-polygon gate (`None` when disabled or no polygon).
    polygon: Option<&'a PolygonGate>,
    /// Built region-of-interest crop (`None` when no point-level crop is requested).
    crop: Option<&'a CropGate>,
    /// Skip all denoising and apply only the crop.
    crop_only: bool,
    /// Per-frame retention-time keep mask (`false` = emit this frame empty).
    rt_keep: &'a [bool],
    /// MS1-stream position -> global frame index (running-average pre-filter).
    ms1_indices: &'a [usize],
    /// Global frame index -> MS1-stream position (`None` for MS/MS frames).
    ms1_pos: &'a [Option<usize>],
}

fn process_frame(
    reader: &FrameReader,
    meta: &[tdf::FrameMeta],
    i: usize,
    params: &FilterParams,
    stages: &Stages,
    ctx: &FrameCtx,
) -> Result<ProcessedFrame> {
    // The stages that act per frame; the polygon/dia_ms1/denoise_msms knobs were
    // already consumed into the gates and keep sets held by `ctx`.
    let &Stages {
        filter_all_frames,
        frame_half_width,
        halo,
        smooth,
        watershed,
        box_centroid: box_centroid_params,
        dia_window,
        dia_per_window,
        ..
    } = stages;
    let &FrameCtx {
        msms,
        dia_msms,
        dia_windows,
        dia_ms1,
        polygon,
        crop,
        crop_only,
        rt_keep,
        ms1_indices,
        ms1_pos,
    } = ctx;
    let meta_i = &meta[i];
    let is_ms1 = meta_i.is_ms1();
    // Empty frames: timsrust cannot decode their absent payload, so emit the
    // canonical empty record directly (Bruker stores these too).
    if meta_i.num_peaks == 0 {
        return Ok(ProcessedFrame {
            frame_id: meta_i.id,
            record: crate::codec::encode_empty_frame_type2(meta_i.num_scans),
            raw_points: 0,
            num_peaks: 0,
            max_intensity: 0,
            summed_intensities: 0,
            raw_summed: 0,
            is_ms1,
            cropped: false,
        });
    }

    // Retention-time crop: a frame outside the window is emitted empty without
    // decoding its payload. Its input points still count toward the raw total (from
    // `NumPeaks`) so the reported reduction reflects the crop.
    if !rt_keep[i] {
        return Ok(ProcessedFrame {
            frame_id: meta_i.id,
            record: crate::codec::encode_empty_frame_type2(meta_i.num_scans),
            raw_points: meta_i.num_peaks,
            num_peaks: 0,
            max_intensity: 0,
            summed_intensities: 0,
            raw_summed: 0,
            is_ms1,
            cropped: true,
        });
    }

    let frame = reader.get(i).map_err(|e| DnoiseError::FrameRead {
        index: i,
        message: e.to_string(),
    })?;
    let flat = FlatFrame::from_frame(&frame);
    let raw_points = flat.len() as u64;
    let raw_summed: u64 = flat.intensity.iter().map(|&it| it as u64).sum();
    let num_scans = flat.num_scans;
    let frame_id = flat.frame_id;

    // Optional cross-frame combine for the keep/drop DECISION only (frame_half_width
    // > 0): sum this MS1 frame with its MS1-frame neighborhood into one spectrum,
    // filter the combined spectrum, and keep this frame's NATIVE points whose
    // (scan, tof) survive there. Neighbors raise signal-to-noise so a faint but
    // persistent feature survives, exactly as the ddaPASEF MS/MS path combines a
    // precursor's re-isolated scans before deciding (msms.rs::combine_and_filter).
    // Unlike the old running-average smoother it never merges or averages points,
    // so output stays a subset of the native frame with native intensities.
    let neighborhood_keys: Option<HashSet<u64>> =
        if frame_half_width > 0 && meta_i.is_ms1() && meta_i.num_peaks > 0 {
            let p = ms1_pos[i].expect("MS1 frame must have an MS1-stream position");
            let lo = p.saturating_sub(frame_half_width);
            let hi = (p + frame_half_width).min(ms1_indices.len() - 1);
            let mut neighbors: Vec<FlatFrame> = Vec::with_capacity(hi - lo);
            for &gi in &ms1_indices[lo..=hi] {
                if gi == i || meta[gi].num_peaks == 0 {
                    continue;
                }
                let nf = reader.get(gi).map_err(|e| DnoiseError::FrameRead {
                    index: gi,
                    message: e.to_string(),
                })?;
                neighbors.push(FlatFrame::from_frame(&nf));
            }
            let mut window: Vec<&FlatFrame> = neighbors.iter().collect();
            window.push(&flat);
            Some(neighborhood_keep_keys(num_scans, &window, params, halo))
        } else {
            None
        };
    // Filtering, the DIA gates, survivors and any centroiding all operate on the
    // native frame; the neighborhood only informs the MS1 keep mask above.
    let to_filter: &FlatFrame = &flat;

    // diaPASEF isolation-window scan intervals for this frame (None for MS1, for
    // ddaPASEF, or when neither DIA feature is enabled). Drives both per-window
    // MS/MS filtering and the out-of-window gate below.
    let dia_iv = dia_windows.and_then(|dw| dw.intervals(meta_i.id));

    // MS1 frames: vertical filter then (optional) horizontal-halo on the survivors.
    // MS/MS frames: pruned by the precursor keep sets when MS/MS denoising is on,
    // otherwise re-encoded unchanged (or vertical-filtered if `filter_all_frames`).
    // In `crop_only` mode no denoising runs at all — every point survives to the
    // crop below, which is the sole filter.
    let mut keep = if crop_only {
        vec![true; to_filter.len()]
    } else if meta_i.is_ms1() {
        let mut keep = if let Some(keys) = &neighborhood_keys {
            // Prune native points by the combined-spectrum decision (mirrors MS/MS).
            (0..to_filter.len())
                .map(|j| keys.contains(&frame_key(to_filter.scan[j], to_filter.tof[j])))
                .collect()
        } else {
            let mut keep = filter_iterated(to_filter, params);
            if let Some(hp) = halo {
                apply_halo(to_filter, hp, &mut keep);
            }
            keep
        };
        // diaPASEF MS1 out-of-window gate: drop surviving points whose (scan, TOF)
        // is in no padded isolation window. Composes as an AND on the keep mask.
        if let Some(gate) = dia_ms1 {
            for (slot, in_win) in keep
                .iter_mut()
                .zip(gate.keep_mask(&to_filter.scan, &to_filter.tof))
            {
                *slot &= in_win;
            }
        }
        // MS1 selection-polygon gate: drop surviving points outside the run's IMS
        // PolygonFilter region (never-selected precursor space). Also ANDed in.
        if let Some(gate) = polygon {
            for (slot, inside) in keep
                .iter_mut()
                .zip(gate.keep_mask(&to_filter.scan, &to_filter.tof))
            {
                *slot &= inside;
            }
        }
        keep
    } else if let Some(mk) = msms {
        mk.keep_mask(to_filter, meta_i.id)
    } else if let Some(mp) = dia_msms {
        // diaPASEF MS/MS: run the same MS/MS filter on each whole frame. With
        // `dia_per_window`, filter each isolation window's scan slice on its own
        // instead, so a mobility run cannot be fused across a window boundary
        // (cross-talk between unrelated isolation events). No cross-frame combine
        // (each DIA window is sampled once per cycle); the `msms_*` knobs apply via
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
    }

    // Region-of-interest crop: AND the `(m/z, 1/K0, intensity)` box into the keep
    // mask. Applies to every frame regardless of MS level (a subset of the raw
    // acquisition), and is the only active filter under `crop_only`.
    if let Some(cg) = crop {
        cg.apply(&to_filter.scan, &to_filter.tof, &to_filter.intensity, &mut keep);
    }

    let survivors = to_filter.survivors(&keep);

    // Optional intensity smoothing then watershed centroiding, both applied only
    // to frames the vertical filter actually processed — MS1 always, MS/MS only
    // under `filter_all_frames` (and never on the separate per-precursor MS/MS-
    // denoise path, which has its own keep logic). Smoothing runs first so the
    // watershed seeds on the stabilised intensities.
    let filtered_here = !crop_only
        && (meta_i.is_ms1() || dia_msms.is_some() || (msms.is_none() && filter_all_frames));
    let survivors = match smooth {
        Some(sp) if filtered_here => box_average(&survivors, num_scans, sp),
        _ => survivors,
    };
    let survivors = match watershed {
        Some(wp) if filtered_here => watershed_centroid(&survivors, wp, MAX_CENTROIDS),
        _ => survivors,
    };
    // Optional final stage: greedy small-box centroiding (mutually exclusive with
    // watershed, enforced by the CLI). Tiles streaks into small centroids rather
    // than collapsing them, preserving the mobility profile.
    let survivors = match box_centroid_params {
        Some(bp) if filtered_here => box_centroid(&survivors, bp),
        _ => survivors,
    };

    let num_peaks = survivors.len() as u64;
    let summed_intensities: u64 = survivors.iter().map(|&(_, _, it)| it as u64).sum();
    let max_intensity = survivors.iter().map(|&(_, _, it)| it).max().unwrap_or(0);
    let record = encode_frame_type2(num_scans, &survivors);

    Ok(ProcessedFrame {
        frame_id,
        record,
        raw_points,
        num_peaks,
        max_intensity,
        summed_intensities,
        raw_summed,
        is_ms1,
        cropped: false,
    })
}

/// Build the diaPASEF MS1 out-of-window gate: read the isolation windows, pad each
/// in physical units (`mz_pad` Da, `im_pad` 1/K0) using the run's calibration,
/// convert to integer `(scan, TOF index)` boxes, and assemble the per-scan lookup.
/// Returns `None` for ddaPASEF (no windows) so the gate is skipped.
fn build_dia_ms1_gate(
    in_tdf: &Path,
    p: &DiaMs1WindowParams,
    meta: &[tdf::FrameMeta],
) -> Result<Option<DiaMs1Gate>> {
    let boxes = tdf::read_dia_ms1_boxes(in_tdf)?;
    if boxes.is_empty() {
        return Ok(None);
    }
    let md = MetadataReader::new(in_tdf).map_err(|e| DnoiseError::Metadata(e.to_string()))?;
    let num_scans = meta.iter().map(|m| m.num_scans).max().unwrap_or(0);
    if num_scans == 0 {
        return Ok(None);
    }

    let tof_boxes: Vec<TofScanBox> = boxes
        .iter()
        .map(|b| {
            // m/z edges -> TOF indices (monotonic), padded by mz_pad Da on each side.
            let t0 = md.mz_converter.invert(b.mz_lo - p.mz_pad);
            let t1 = md.mz_converter.invert(b.mz_hi + p.mz_pad);
            let tof_lo = t0.min(t1).floor().max(0.0) as u32;
            let tof_hi = t1.max(t0).ceil().max(0.0) as u32;

            // Scan range -> 1/K0 (monotonic decreasing), padded by im_pad, back to
            // scans. Take min/max so the result is correct regardless of direction.
            let im0 = md.im_converter.convert(b.scan_begin);
            let im1 = md.im_converter.convert(b.scan_end);
            let s0 = md.im_converter.invert(im0.max(im1) + p.im_pad);
            let s1 = md.im_converter.invert(im0.min(im1) - p.im_pad);
            let scan_lo = s0.min(s1).floor().max(0.0) as u32;
            let scan_hi = (s0.max(s1).ceil().max(0.0) as u32).min(num_scans as u32 - 1);

            TofScanBox {
                scan_lo,
                scan_hi,
                tof_lo,
                tof_hi,
            }
        })
        .collect();

    Ok(DiaMs1Gate::build(&tof_boxes, num_scans))
}

/// Build the MS1 selection-polygon gate: read the run's IMS PolygonFilter
/// `(m/z, 1/K0)` vertices, convert them to per-scan TOF-index intervals via the
/// run calibration (padded by `mz_pad` Da / `im_pad` 1/K0), and assemble the
/// per-scan lookup. Returns `None` when the run stores no polygon so the gate is
/// skipped.
///
/// **ddaPASEF only.** In ddaPASEF the IMS PolygonFilter is a single ring bounding
/// the precursor-selection region. In diaPASEF the same property instead stores
/// several disjoint quads (the window-placement anchors), which are *not* a
/// selection region — and diaPASEF MS1 windowing is already handled by the
/// [`crate::dia_ms1`] gate. So the polygon gate is skipped on any run that defines
/// a diaPASEF window scheme, to avoid misreading those quads as one polygon.
fn build_polygon_gate(
    in_tdf: &Path,
    p: &Ms1PolygonParams,
    meta: &[tdf::FrameMeta],
) -> Result<Option<PolygonGate>> {
    if !tdf::read_dia_windows(in_tdf)?.is_empty() {
        return Ok(None); // diaPASEF: the polygon property is multi-component here.
    }
    let Some((mz, im)) = tdf::read_selection_polygon(in_tdf)? else {
        return Ok(None);
    };
    let md = MetadataReader::new(in_tdf).map_err(|e| DnoiseError::Metadata(e.to_string()))?;
    let num_scans = meta.iter().map(|m| m.num_scans).max().unwrap_or(0);
    if num_scans == 0 {
        return Ok(None);
    }
    Ok(PolygonGate::build(
        &mz,
        &im,
        num_scans,
        |s| md.im_converter.convert(s),
        |mz| md.mz_converter.invert(mz),
        p.mz_pad,
        p.im_pad,
    ))
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

/// Pack an absolute `(scan, tof)` into a u64 key (scan in the high 32 bits).
fn frame_key(scan: u32, tof: u32) -> u64 {
    ((scan as u64) << 32) | tof as u64
}

/// Combine an MS1-frame neighborhood into one summed `(scan, tof)` spectrum, run the
/// vertical + horizontal-halo filter on it, and return the surviving `(scan, tof)`
/// keys. Used by the `frame_half_width` path: the combined spectrum only informs the
/// keep/drop decision (raising signal-to-noise for persistent features); the caller
/// prunes the current frame's native points by these keys, never merging intensities.
/// Mirrors the ddaPASEF MS/MS combine-then-decide path (`msms.rs::combine_and_filter`).
fn neighborhood_keep_keys(
    num_scans: usize,
    window: &[&FlatFrame],
    params: &FilterParams,
    halo: Option<&HaloParams>,
) -> HashSet<u64> {
    let mut acc: HashMap<u64, u64> = HashMap::new();
    for f in window {
        for k in 0..f.len() {
            *acc.entry(frame_key(f.scan[k], f.tof[k])).or_insert(0) += f.intensity[k] as u64;
        }
    }
    let mut scan = Vec::with_capacity(acc.len());
    let mut tof = Vec::with_capacity(acc.len());
    let mut intensity = Vec::with_capacity(acc.len());
    for (&k, &sum) in &acc {
        scan.push((k >> 32) as u32);
        tof.push((k & 0xFFFF_FFFF) as u32);
        intensity.push(sum.min(u32::MAX as u64) as u32);
    }
    let combined = FlatFrame {
        frame_id: 0,
        num_scans,
        scan,
        tof,
        intensity,
    };
    let mut keep = filter_iterated(&combined, params);
    if let Some(hp) = halo {
        apply_halo(&combined, hp, &mut keep);
    }
    let mut out = HashSet::new();
    for ((&keep_k, &scan), &tof) in keep.iter().zip(&combined.scan).zip(&combined.tof) {
        if keep_k {
            out.insert(frame_key(scan, tof));
        }
    }
    out
}

/// Recursively copy `src` into `dst`, skipping a top-level entry named `skip_top`.
fn copy_dir_except(src: &Path, dst: &Path, skip_top: &str) -> Result<()> {
    fs::create_dir_all(dst)?;
    for entry in fs::read_dir(src)? {
        let entry = entry?;
        let name = entry.file_name();
        let from = entry.path();
        let to = dst.join(&name);
        if from.is_dir() {
            copy_dir_recursive(&from, &to)?;
        } else if name != skip_top {
            fs::copy(&from, &to)?;
        }
    }
    Ok(())
}

fn copy_dir_recursive(src: &Path, dst: &Path) -> Result<()> {
    fs::create_dir_all(dst)?;
    for entry in fs::read_dir(src)? {
        let entry = entry?;
        let from = entry.path();
        let to = dst.join(entry.file_name());
        if from.is_dir() {
            copy_dir_recursive(&from, &to)?;
        } else {
            fs::copy(&from, &to)?;
        }
    }
    Ok(())
}
