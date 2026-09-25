//! Orchestration: copy the source `.d`, rewrite `analysis.tdf_bin` with filtered
//! frames (re-encoded as type 2), and fix up the `analysis.tdf` SQLite database.

use crate::box_centroid::box_centroid;
use crate::codec::try_encode_frame_type2;
use crate::crop::CropGate;
use crate::dia_ms1::{DiaMs1Gate, TofScanBox};
use crate::dia_window::{filter_per_window, in_window_mask};
use crate::error::{DnoiseError, Result};
use crate::filter::filter_iterated;
use crate::frame::FlatFrame;
use crate::halo::horizontal_halo_keep_mask;
use crate::mobility::{self, MobilityScale, ScanToMobility};
use crate::msms::{MsmsKeep, build_msms_keep};
use crate::neighbor::NeighborIndex;
use crate::params::{
    CropParams, DiaMs1WindowParams, FilterParams, HaloParams, Ms1PolygonParams, MsmsFilterParams,
    Stages,
};
use crate::polygon::PolygonGate;
use crate::provenance::NeighborUsage;
use crate::smooth::box_average;
use crate::tdf::{self, DiaWindows, FrameUpdate, PrmWindows};
use crate::watershed::watershed_centroid;
use rayon::prelude::*;
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
#[derive(Debug, Clone, Copy, Default, serde::Serialize)]
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
    /// Input MS1 summed intensity, including RT-cropped frames.
    pub raw_ms1_summed_intensity: u64,
    /// Output MS1 summed intensity after all enabled stages.
    pub kept_ms1_summed_intensity: u64,
    /// Input MS/MS summed intensity, including RT-cropped frames.
    pub raw_msms_summed_intensity: u64,
    /// Output MS/MS summed intensity after all enabled stages.
    pub kept_msms_summed_intensity: u64,
    /// Actual MS1 temporal evidence (central frames are excluded).
    pub ms1_neighbor_usage: NeighborUsage,
    /// Actual PRM/DIA temporal evidence (central events are excluded).
    pub msms_neighbor_usage: NeighborUsage,
    /// True when this was a dry run (no output written).
    pub dry_run: bool,
    /// Actual input frame binary size, not a point-count estimate.
    pub input_binary_bytes: u64,
    /// Actual output frame binary size (zero for a dry run).
    pub output_binary_bytes: u64,
    /// Elapsed processing/validation time in seconds.
    pub elapsed_seconds: f64,
    /// Worker pool size actually used by this run.
    pub worker_threads: usize,
    /// Active geometry gates, after acquisition detection.
    pub active_gates: crate::provenance::ActiveGates,
    /// Multiple calibration segments were detected for physical gates or crops.
    pub multiple_calibrations: bool,
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
    /// Requires at least one crop bound; invalid configurations return an error.
    pub crop_only: bool,
    /// Dry-run frame sampling for a fast reduction estimate (`None` = all frames).
    /// Only honoured together with `dry_run`.
    pub sample: Option<SampleSpec>,
    /// Cooperative cancellation token. When set and flipped to `true`, the run
    /// stops before another frame starts (in-flight frames finish) and returns [`DnoiseError::Cancelled`];
    /// partial temporary output is removed automatically. `None` = never cancelled.
    pub cancel: Option<&'a AtomicBool>,
    /// Maximum frames encoded per batch. None preserves the default (2048).
    pub frame_batch_size: Option<usize>,
    /// Skip separate full-frame decoding passes on input and output. Structural
    /// checks, checked encoding, path protection, and staged installation remain.
    /// False by default. Damage detectable only by decoding may go undetected.
    pub skip_validation: bool,
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
    execute(input, output, params, stages, options, progress, false)
}

/// Replace an input only after successfully processing and validating a temporary copy.
/// A failed installation restores the original or reports its preserved backup path.
pub fn denoise_in_place<F: FnMut(Progress)>(
    input: &Path,
    params: &FilterParams,
    stages: &Stages,
    options: &RunOptions,
    progress: F,
) -> Result<DenoiseStats> {
    if options.dry_run {
        return Err(DnoiseError::InvalidInput(
            "in-place cannot be a dry run".into(),
        ));
    }
    execute(input, input, params, stages, options, progress, true)
}

#[allow(clippy::too_many_arguments)]
fn execute<F: FnMut(Progress)>(
    input: &Path,
    output: &Path,
    params: &FilterParams,
    stages: &Stages,
    options: &RunOptions,
    progress: F,
    in_place: bool,
) -> Result<DenoiseStats> {
    let started = std::time::Instant::now();
    crate::validation::parameters(params, stages, options)?;
    if !options.dry_run && !in_place {
        crate::output::check_disjoint(input, output)?;
    }
    if options.cancel.is_some_and(|c| c.load(Ordering::Relaxed)) {
        return Err(DnoiseError::Cancelled);
    }
    crate::validation::inspect_cancellable(input, !options.skip_validation, options.cancel)?;
    crate::provenance::read(input)?;
    if options.dry_run {
        let mut stats = run(input, output, params, stages, options, progress)?;
        stats.elapsed_seconds = started.elapsed().as_secs_f64();
        return Ok(stats);
    }
    let tx = crate::output::OutputTransaction::begin(
        input,
        output,
        options.force || in_place,
        in_place,
    )?;
    let mut stats = run(input, tx.path(), params, stages, options, progress)?;
    crate::validation::inspect_cancellable(tx.path(), !options.skip_validation, options.cancel)?;
    if options.cancel.is_some_and(|c| c.load(Ordering::Relaxed)) {
        return Err(DnoiseError::Cancelled);
    }
    stats.elapsed_seconds = started.elapsed().as_secs_f64();
    crate::provenance::write(input, tx.path(), params, stages, options, &stats)?;
    tx.commit()?;
    Ok(stats)
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
    denoise_with_options(input, output, params, stages, &opts, progress)
}

/// The writer and streaming reader use the same acquisition policy. Targeted,
/// mixed and unknown runs may use the MS1 filter, but never discovery geometry
/// or the whole-frame MS/MS fallback. Cropping remains an explicit operation.
fn acquisition_stages<'a>(
    acquisition: crate::Acquisition,
    stages: &Stages<'a>,
    crop_only: bool,
) -> Result<Stages<'a>> {
    let mut effective = *stages;
    if matches!(
        acquisition,
        crate::Acquisition::PrmPasef | crate::Acquisition::Mixed | crate::Acquisition::Unknown
    ) {
        if !crop_only && acquisition == crate::Acquisition::PrmPasef && stages.filter_all_frames {
            return Err(DnoiseError::InvalidInput(
                "--all-frames (filter_all_frames) is unsupported for prm-PASEF; use --denoise-msms for experimental per-event fragment filtering".into()
            ));
        }
        if !crop_only
            && acquisition != crate::Acquisition::PrmPasef
            && (stages.denoise_msms.is_some() || stages.filter_all_frames)
        {
            return Err(DnoiseError::InvalidInput(format!(
                "MS/MS denoising is not supported for {acquisition}; use MS1-only processing (disable denoise_msms and filter_all_frames)"
            )));
        }
        if crop_only || acquisition != crate::Acquisition::PrmPasef {
            effective.denoise_msms = None;
        }
        effective.filter_all_frames = false;
        effective.ms1_polygon = None;
        effective.dia_ms1 = None;
        effective.dia_window = None;
        effective.dda_window = None;
        effective.dia_per_window = false;
    }
    if acquisition == crate::Acquisition::DiaPasef && stages.neighbors.dia_radius > 0 {
        effective.dia_per_window = true;
    }
    Ok(effective)
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
        force: _,
        dry_run,
        crop,
        crop_only,
        sample,
        cancel,
        frame_batch_size,
        skip_validation: _,
    } = options;
    let in_tdf = input.join("analysis.tdf");
    let in_bin = input.join("analysis.tdf_bin");
    if !in_tdf.is_file() || !in_bin.is_file() {
        return Err(DnoiseError::NotADotD(input.to_path_buf()));
    }

    // Frame metadata (ordered by Id == timsrust index). Empty frames are handled
    // without timsrust, which cannot decode their absent payload.
    let meta = tdf::read_frame_meta(&in_tdf)?;

    let acquisition = tdf::inspect_acquisition(&in_tdf, &meta)?;
    let scheme = acquisition.kind;
    let reader = FrameReader::new(input).map_err(|e| DnoiseError::OpenFrames(e.to_string()))?;
    let n_frames = reader.len();
    let effective_stages = acquisition_stages(scheme, stages, crop_only)?;
    let stages = &effective_stages;
    // Unpack the stages this function builds gates from; the per-frame stages
    // (smoothing, centroiding, etc.) are forwarded to `process_frame` via `stages`.
    let &Stages {
        halo,
        denoise_msms,
        dia_window,
        dda_window,
        dia_per_window,
        dia_ms1,
        ms1_polygon,
        ..
    } = stages;

    let neighbors = if crop_only {
        None
    } else {
        NeighborIndex::build(&in_tdf, &meta, scheme, stages)?
    };

    // Destination is an owned temporary directory; the transaction installs it later.
    if !dry_run {
        copy_dir_except(input, output, "analysis.tdf_bin")?;
    }
    let n_ms1 = meta.iter().filter(|m| m.is_ms1()).count();
    let n_empty = meta.iter().filter(|m| m.num_peaks == 0).count();
    if n_ms1 == 0 && !crop_only && stages.denoise_msms.is_none() && !stages.filter_all_frames {
        warn!("no MS1 frames; MS1-only denoising will retain all points");
    }
    info!(input = %input.display(), output = %output.display(), "denoise: starting");
    info!(
        scheme = %scheme,
        frames = n_frames,
        ms1 = n_ms1,
        msms = n_frames - n_ms1,
        empty = n_empty,
        "denoise: frame inventory"
    );

    // MS/MS denoising splits by acquisition scheme, driven by the same
    // `denoise_msms` params:
    //   * ddaPASEF — each precursor is re-isolated across several frames, so we
    //     build per-precursor keep sets up front (PasefFrameMsMsInfo) and combine
    //     a precursor's fragment scans across frames before filtering.
    //   * diaPASEF — each isolation window is filtered independently by default;
    //     optional neighbors provide evidence across compatible observations.
    //   * prm-PASEF — recorded target boundaries are always enforced, including
    //     when nearby observations supply temporal evidence.
    let prm_windows = (scheme == crate::Acquisition::PrmPasef && denoise_msms.is_some())
        .then_some(acquisition.prm_windows);
    if prm_windows.is_some() {
        warn!(
            "experimental prm-PASEF MS/MS denoising: target boundaries enforced; validate downstream quantification"
        );
    }
    let (msms_keep, dia_msms) = match denoise_msms {
        Some(_) if prm_windows.is_some() => (None, None),
        Some(mp) => {
            let windows = tdf::read_pasef_msms(&in_tdf)?;
            if windows.is_empty() {
                info!("MS/MS denoise: diaPASEF whole-frame path (no PasefFrameMsMsInfo)");
                (None, Some(mp))
            } else {
                let keep = build_msms_keep(&reader, &meta, &windows, mp, halo, cancel)?;
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
    let dia_regions = if scheme == crate::Acquisition::DiaPasef
        && dia_per_window
        && (denoise_msms.is_some() || stages.filter_all_frames)
        && !crop_only
    {
        Some(tdf::dia::read(&in_tdf)?)
    } else {
        None
    };

    // ddaPASEF isolation-event intervals for the MS/MS out-of-window gate, in the
    // same per-frame shape as the diaPASEF scheme. Empty for diaPASEF (no
    // PasefFrameMsMsInfo), so the gate no-ops there — and on standard timsTOF
    // ddaPASEF files it is expected to remove nothing (the acquisition writes
    // MS/MS scans only inside scheduled isolation events); it runs as a guarantee.
    let dda_windows: Option<DiaWindows> = if dda_window.is_some() {
        let w = DiaWindows::from_pasef(&tdf::read_pasef_msms(&in_tdf)?);
        if w.is_empty() {
            debug!(
                "ddaPASEF window gate requested but no isolation events found (diaPASEF?) — skipped"
            );
            None
        } else {
            info!("ddaPASEF MS/MS out-of-window gate active");
            Some(w)
        }
    } else {
        None
    };
    let dda_windows_ref = dda_windows.as_ref();

    // diaPASEF MS1 out-of-window gate: build the padded `(scan, TOF)` lookup once
    // from the isolation windows + calibration. `None` for ddaPASEF (no windows).
    let dia_ms1_gate = match dia_ms1.filter(|_| !crop_only) {
        Some(mp) => build_dia_ms1_gate(&in_tdf, mp, &meta, stages.mobility_scale)?,
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
    let polygon_gate = match ms1_polygon.filter(|_| !crop_only) {
        Some(pp) => build_polygon_gate(&in_tdf, pp, &meta, stages.mobility_scale)?,
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
            require_single_calibration(
                &in_tdf,
                cp.mz_min.is_some() || cp.mz_max.is_some(),
                cp.im_min.is_some() || cp.im_max.is_some(),
                "physical crop (remove the m/z or mobility bounds)",
            )?;
            let md =
                MetadataReader::new(&in_tdf).map_err(|e| DnoiseError::Metadata(e.to_string()))?;
            let num_scans = meta.iter().map(|m| m.num_scans).max().unwrap_or(0);
            let im = if cp.im_min.is_some() || cp.im_max.is_some() {
                mobility::load(&in_tdf, stages.mobility_scale, md.im_converter)?
            } else {
                ScanToMobility::Linear(md.im_converter)
            };
            let g = CropGate::build(
                cp,
                num_scans,
                |mz| md.mz_converter.invert(mz),
                |k0| im.invert(k0),
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
    // gates derived above, the crop, plus compatible-observation indices. Bundling these
    // keeps `process_frame` to a handful of arguments.
    let ctx = FrameCtx {
        msms: msms_ref,
        prm_msms: denoise_msms.zip(prm_windows.as_ref()),
        dia_msms,
        dia_windows: dia_windows_ref,
        dia_regions: dia_regions.as_ref(),
        dda_windows: dda_windows_ref,
        dia_ms1: dia_ms1_ref,
        polygon: polygon_ref,
        crop: crop_ref,
        crop_only,
        rt_keep: &rt_keep,
        neighbors: neighbors.as_ref(),
        cancel,
    };

    let mut offset: u64 = header_len;
    let mut updates: Vec<FrameUpdate> = Vec::with_capacity(n_process);
    let mut raw_points: u64 = 0;
    let mut kept_points: u64 = 0;
    let mut raw_ms1: u64 = 0;
    let mut kept_ms1: u64 = 0;
    let mut raw_summed: u64 = 0;
    let mut kept_summed: u64 = 0;
    let mut raw_ms1_summed = 0;
    let mut kept_ms1_summed = 0;
    let mut ms1_neighbor_usage = NeighborUsage::default();
    let mut msms_neighbor_usage = NeighborUsage::default();
    let mut cropped_frames: usize = 0;
    let mut frames_done: usize = 0;

    for chunk in selected.chunks(frame_batch_size.unwrap_or(CHUNK)) {
        // Cooperative cancellation: check once per chunk (a real run's partial
        // output is incomplete, so the caller discards it on Cancelled).
        if let Some(c) = cancel {
            if c.load(Ordering::Relaxed) {
                return Err(DnoiseError::Cancelled);
            }
        }
        let processed: Vec<ProcessedFrame> = chunk
            .par_iter()
            .map(|&i| {
                if cancel.is_some_and(|c| c.load(Ordering::Relaxed)) {
                    return Err(DnoiseError::Cancelled);
                }
                process_frame(&reader, &meta, i, params, stages, &ctx)
            })
            .collect::<Result<_>>()?;

        for pf in processed {
            if cancel.is_some_and(|c| c.load(Ordering::Relaxed)) {
                return Err(DnoiseError::Cancelled);
            }
            raw_points += pf.raw_points;
            kept_points += pf.num_peaks;
            raw_summed += pf.raw_summed;
            kept_summed += pf.summed_intensities;
            if pf.is_ms1 {
                raw_ms1 += pf.raw_points;
                kept_ms1 += pf.num_peaks;
                raw_ms1_summed += pf.raw_summed;
                kept_ms1_summed += pf.summed_intensities;
                ms1_neighbor_usage.add(pf.neighbor_usage);
            } else {
                msms_neighbor_usage.add(pf.neighbor_usage);
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
        b.get_ref().sync_all()?;
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
        raw_ms1_summed_intensity: raw_ms1_summed,
        kept_ms1_summed_intensity: kept_ms1_summed,
        raw_msms_summed_intensity: raw_summed - raw_ms1_summed,
        kept_msms_summed_intensity: kept_summed - kept_ms1_summed,
        ms1_neighbor_usage,
        msms_neighbor_usage,
        dry_run,
        input_binary_bytes: fs::metadata(&in_bin)?.len(),
        output_binary_bytes: if dry_run {
            0
        } else {
            fs::metadata(output.join("analysis.tdf_bin"))?.len()
        },
        elapsed_seconds: 0.0,
        worker_threads: rayon::current_num_threads(),
        active_gates: crate::provenance::ActiveGates {
            ms1_neighbors: neighbors.is_some() && stages.frame_half_width > 0 && n_ms1 > 0,
            prm_neighbors: neighbors.is_some()
                && stages.neighbors.prm_radius > 0
                && scheme == crate::Acquisition::PrmPasef,
            dia_neighbors: neighbors.is_some()
                && stages.neighbors.dia_radius > 0
                && scheme == crate::Acquisition::DiaPasef,
            dia_scan_varying: dia_regions.as_ref().is_some_and(|r| r.has_scanning()),
            prm_per_event: prm_windows.is_some() && !crop_only,
            ms1_polygon: polygon_ref.is_some() && !crop_only,
            dia_ms1: dia_ms1_ref.is_some() && !crop_only,
            dia_window: dia_windows_ref.is_some()
                && dia_window.is_some()
                && !crop_only
                && (stages.filter_all_frames || denoise_msms.is_some()),
            dda_window: dda_windows_ref.is_some()
                && !crop_only
                && (stages.filter_all_frames || denoise_msms.is_some()),
            dia_per_window: dia_windows_ref.is_some()
                && dia_per_window
                && !crop_only
                && (stages.filter_all_frames || denoise_msms.is_some()),
        },
        multiple_calibrations: {
            let (mz, im) = tdf::count_calibration_segments(&in_tdf)?;
            mz > 1 || im > 1
        },
    })
}

struct ProcessedFrame {
    neighbor_usage: NeighborUsage,
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
/// once in [`run`], plus the crop and compatible-observation indices. Lets the
/// per-frame worker take the run's derived state as a single value.
pub struct FrameCtx<'a> {
    /// PRM event-local filter knobs and checked, unmerged isolation intervals.
    prm_msms: Option<(&'a MsmsFilterParams, &'a PrmWindows)>,
    /// ddaPASEF per-precursor keep sets (`None` unless MS/MS denoising on ddaPASEF).
    msms: Option<&'a MsmsKeep>,
    /// diaPASEF MS/MS filter knobs (`None` unless MS/MS denoising on diaPASEF).
    dia_msms: Option<&'a MsmsFilterParams>,
    /// diaPASEF isolation windows (`None` for ddaPASEF or when unused).
    dia_windows: Option<&'a DiaWindows>,
    /// Checked static windows and continuous scan-dependent DIA regions.
    dia_regions: Option<&'a tdf::dia::DiaRegions>,
    /// ddaPASEF isolation-event intervals (`None` for diaPASEF or when unused).
    dda_windows: Option<&'a DiaWindows>,
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
    /// Compatible event neighborhoods for temporal filtering evidence.
    neighbors: Option<&'a NeighborIndex>,
    /// Cooperative cancellation while decoding supporting observations.
    cancel: Option<&'a std::sync::atomic::AtomicBool>,
}

/// A frame after every denoising / crop stage has run, but *before* it is
/// re-encoded into a `.d` record. This is the unit the streaming API
/// ([`crate::RunContext`]) hands to in-process callers that want the surviving
/// points directly (e.g. a feature finder) instead of a rewritten `.d` on disk.
/// [`process_frame`] wraps this with the type-2 encoder to write the file.
pub struct DecodedFrame {
    /// Bruker frame `Id` (== timsrust frame index).
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

/// Thin file-writer wrapper over [`process_frame_decoded`]: decode the frame's
/// survivors, then encode them into a `.d` type-2 record plus the per-frame stats
/// the run's metadata fixup needs. Behaviour-identical to the pre-streaming
/// implementation (proven by the byte-identity test in `tests/`).
fn process_frame(
    reader: &FrameReader,
    meta: &[tdf::FrameMeta],
    i: usize,
    params: &FilterParams,
    stages: &Stages,
    ctx: &FrameCtx,
) -> Result<ProcessedFrame> {
    let d = process_frame_decoded(reader, meta, i, params, stages, ctx)?;
    let num_peaks = d.survivors.len() as u64;
    let summed_intensities: u64 = d.survivors.iter().map(|&(_, _, it)| it as u64).sum();
    let max_intensity = d.survivors.iter().map(|&(_, _, it)| it).max().unwrap_or(0);
    let record = try_encode_frame_type2(d.num_scans, &d.survivors)?;
    Ok(ProcessedFrame {
        neighbor_usage: d.neighbor_usage,
        frame_id: d.frame_id,
        record,
        raw_points: d.raw_points,
        num_peaks,
        max_intensity,
        summed_intensities,
        raw_summed: d.raw_summed,
        is_ms1: d.is_ms1,
        cropped: d.cropped,
    })
}

/// Decode one frame through every enabled denoising / crop stage and return its
/// surviving points, without touching the output `.d`. This is the shared core
/// behind both the file writer ([`process_frame`]) and the streaming API
/// ([`crate::RunContext::process`]) — a single implementation, so an in-process
/// caller can never drift from the standalone tool.
pub fn process_frame_decoded(
    reader: &FrameReader,
    meta: &[tdf::FrameMeta],
    i: usize,
    params: &FilterParams,
    stages: &Stages,
    ctx: &FrameCtx,
) -> Result<DecodedFrame> {
    // The stages that act per frame; the polygon/dia_ms1/denoise_msms knobs were
    // already consumed into the gates and keep sets held by `ctx`.
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
    } = stages;
    let &FrameCtx {
        prm_msms,
        msms,
        dia_msms,
        dia_windows,
        dia_regions,
        dda_windows,
        dia_ms1,
        polygon,
        crop,
        crop_only,
        rt_keep,
        neighbors,
        cancel,
    } = ctx;
    let meta_i = &meta[i];
    let is_ms1 = meta_i.is_ms1();
    // Empty input frames are not read through timsrust, which cannot decode an
    // absent payload. They are still WRITTEN the long way (an all-zero scan
    // table, ~24 compressed bytes) so the output stays readable by timsrust;
    // see `codec::encode_empty_frame_type2`.
    if meta_i.num_peaks == 0 {
        return Ok(DecodedFrame {
            frame_id: meta_i.id,
            is_ms1,
            num_scans: meta_i.num_scans,
            rt_seconds: meta_i.rt,
            survivors: Vec::new(),
            raw_points: 0,
            raw_summed: 0,
            neighbor_usage: NeighborUsage::default(),
            cropped: false,
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

    // Decode RT-excluded frames too: QC denominators must include removed intensity.
    if !rt_keep[i] {
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
        stages
            .denoise_msms
            .map(|p| p.as_filter_params())
            .unwrap_or(*params)
    };
    let neighborhood_mask = match neighbors.filter(|_| !crop_only) {
        Some(n) => n.keep_mask(
            reader,
            meta,
            i,
            &flat,
            &neighbor_params,
            halo,
            rt_keep,
            cancel,
        )?,
        None => None,
    };
    let (neighborhood_mask, neighbor_usage) = match neighborhood_mask {
        Some((mask, usage)) => (Some(mask), usage),
        None => (None, NeighborUsage::default()),
    };
    let to_filter: &FlatFrame = &flat;

    // diaPASEF isolation-window scan intervals for this frame (None for MS1, for
    // ddaPASEF, or when neither DIA feature is enabled). Drives both per-window
    // MS/MS filtering and the out-of-window gate below.
    let dia_iv = dia_regions
        .and_then(|r| r.intervals(meta_i.id))
        .or_else(|| dia_windows.and_then(|dw| dw.intervals(meta_i.id)));

    // MS1 frames: vertical filter then (optional) horizontal-halo on the survivors.
    // MS/MS frames: pruned by the precursor keep sets when MS/MS denoising is on,
    // otherwise re-encoded unchanged (or vertical-filtered if `filter_all_frames`).
    // In `crop_only` mode no denoising runs at all — every point survives to the
    // crop below, which is the sole filter.
    let mut keep = if crop_only {
        vec![true; to_filter.len()]
    } else if meta_i.is_ms1() {
        let mut keep = if let Some(mask) = neighborhood_mask {
            mask
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
            let mut mask = gate.keep_mask(&to_filter.scan, &to_filter.tof);
            if let Some(reach) = gate.overlap {
                mask = overlap_mask(to_filter, &keep, &mask, params, reach);
            }
            for (slot, in_win) in keep.iter_mut().zip(mask) {
                *slot &= in_win;
            }
        }
        // MS1 selection-polygon gate: drop surviving points outside the run's IMS
        // PolygonFilter region (never-selected precursor space). Also ANDed in.
        if let Some(gate) = polygon {
            let mut mask = gate.keep_mask(&to_filter.scan, &to_filter.tof);
            if let Some(reach) = gate.overlap {
                mask = overlap_mask(to_filter, &keep, &mask, params, reach);
            }
            for (slot, inside) in keep.iter_mut().zip(mask) {
                *slot &= inside;
            }
        }
        keep
    } else if let Some(mask) = neighborhood_mask {
        mask
    } else if let Some((mp, windows)) = prm_msms {
        let intervals = windows.get(&meta_i.id).ok_or_else(|| {
            DnoiseError::InvalidInput(format!(
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
        let dda_iv = dda_windows.and_then(|w| w.intervals(meta_i.id));
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
    if let Some(cg) = crop {
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

// Raw-coordinate filters and RT/intensity-only crops do not need these conversions.
fn require_single_calibration(path: &Path, mz: bool, im: bool, operation: &str) -> Result<()> {
    if !mz && !im {
        return Ok(());
    }
    let (mz_count, im_count) = tdf::count_calibration_segments(path)?;
    if (mz && mz_count > 1) || (im && im_count > 1) {
        return Err(DnoiseError::InvalidInput(format!(
            "{operation} cannot use run-level conversion with multiple calibration references (m/z: {mz_count}, mobility: {im_count}); per-calibration physical filtering is not supported"
        )));
    }
    Ok(())
}

/// Build the diaPASEF MS1 out-of-window gate: read the isolation windows, pad each
/// in physical units (`mz_pad` Da, `im_pad` 1/K0) using the run's calibration,
/// convert to integer `(scan, TOF index)` boxes, and assemble the per-scan lookup.
/// Returns `None` for ddaPASEF (no windows) so the gate is skipped.
fn build_dia_ms1_gate(
    in_tdf: &Path,
    p: &DiaMs1WindowParams,
    meta: &[tdf::FrameMeta],
    scale: MobilityScale,
) -> Result<Option<DiaMs1Gate>> {
    let boxes = tdf::read_dia_ms1_boxes(in_tdf)?;
    if boxes.is_empty() {
        return Ok(None);
    }
    require_single_calibration(
        in_tdf,
        true,
        true,
        "DIA MS1 gate (disable dia_ms1_window / --no-dia-ms1-window)",
    )?;
    let md = MetadataReader::new(in_tdf).map_err(|e| DnoiseError::Metadata(e.to_string()))?;
    let num_scans = meta.iter().map(|m| m.num_scans).max().unwrap_or(0);
    if num_scans == 0 {
        return Ok(None);
    }
    let im = mobility::load(in_tdf, scale, md.im_converter)?;

    let tof_boxes: Vec<TofScanBox> = boxes
        .iter()
        .map(|b| {
            // m/z edges -> TOF indices (monotonic), padded by mz_pad Da on each side.
            let t0 = md.mz_converter.invert(b.mz_lo - p.mz_pad);
            let t1 = md.mz_converter.invert(b.mz_hi + p.mz_pad);
            let tof_lo = t0.min(t1).floor().max(0.0) as u32;
            let tof_hi = t1.max(t0).ceil().max(0.0) as u32;

            let (scan_lo, scan_hi) =
                window_scans(b.scan_begin, b.scan_end, p.im_pad, &im, num_scans);

            TofScanBox {
                scan_lo,
                scan_hi,
                tof_lo,
                tof_hi,
            }
        })
        .collect();

    let reach = p
        .overlap
        .then(|| reach_in_scans(p.overlap_reach, |s| im.convert(s as f64), num_scans));
    Ok(DiaMs1Gate::build(&tof_boxes, num_scans).map(|mut g| {
        g.overlap = reach;
        g
    }))
}

/// Inclusive scan range of a `[scan_begin, scan_end)` isolation window, padded by
/// `im_pad` in 1/K0. With no pad the window's own scans are returned exactly;
/// otherwise the edges go scan -> 1/K0 -> padded -> scan, taking min/max so the
/// result is right for either conversion direction. A padded window is rounded
/// outward; the rounding tolerance keeps float noise in the round trip from
/// adding a scan at either edge.
fn window_scans(
    scan_begin: u32,
    scan_end: u32,
    im_pad: f64,
    im: &ScanToMobility,
    num_scans: usize,
) -> (u32, u32) {
    const EPS: f64 = 1e-3; // scans: far above round-trip noise, far below one scan
    let last = num_scans.saturating_sub(1) as u32;
    let first_in = scan_begin;
    let last_in = scan_end.saturating_sub(1).max(scan_begin);
    if im_pad <= 0.0 {
        return (first_in.min(last), last_in.min(last));
    }
    let im0 = im.convert(first_in.into());
    let im1 = im.convert(last_in.into());
    let s0 = im.invert(im0.max(im1) + im_pad);
    let s1 = im.invert(im0.min(im1) - im_pad);
    let lo = (s0.min(s1) + EPS).floor().max(0.0) as u32;
    let hi = ((s0.max(s1) - EPS).ceil().max(0.0) as u32).min(last);
    (lo.min(hi), hi)
}

/// Convert an overlap reach in 1/K0 to mobility scans using the run's mean
/// 1/K0-per-scan slope, rounded up. `0.0` (unlimited) maps to `u32::MAX`.
fn reach_in_scans(reach: f64, im_at_scan: impl Fn(u32) -> f64, num_scans: usize) -> u32 {
    if reach <= 0.0 {
        return u32::MAX;
    }
    let last = num_scans.saturating_sub(1).max(1) as u32;
    let per_scan = (im_at_scan(0) - im_at_scan(last)).abs() / f64::from(last);
    if per_scan.is_nan() || per_scan <= 0.0 {
        return u32::MAX;
    }
    (reach / per_scan).ceil().min(f64::from(u32::MAX - 1)) as u32
}

/// Feature-level gate mask: link surviving points with the streak filter's own
/// adjacency (column half-width, bridged gap + 1) and keep features that touch
/// the gate anywhere, up to `reach` scans beyond their inside points. See
/// [`crate::overlap`].
fn overlap_mask(
    frame: &FlatFrame,
    keep: &[bool],
    inside: &[bool],
    params: &FilterParams,
    reach: u32,
) -> Vec<bool> {
    crate::overlap::extend_to_features(
        &frame.scan,
        &frame.tof,
        keep,
        inside,
        params.mz_half_width,
        params.max_internal_gap as u32 + 1,
        reach,
    )
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
    scale: MobilityScale,
) -> Result<Option<PolygonGate>> {
    if !tdf::read_dia_windows(in_tdf)?.is_empty() {
        return Ok(None); // diaPASEF: the polygon property is multi-component here.
    }
    let Some((mz, im)) = tdf::read_selection_polygon(in_tdf)? else {
        return Ok(None);
    };
    require_single_calibration(
        in_tdf,
        true,
        true,
        "MS1 polygon gate (disable ms1_polygon / --no-ms1-polygon)",
    )?;
    let md = MetadataReader::new(in_tdf).map_err(|e| DnoiseError::Metadata(e.to_string()))?;
    let num_scans = meta.iter().map(|m| m.num_scans).max().unwrap_or(0);
    if num_scans == 0 {
        return Ok(None);
    }
    let k0 = mobility::load(in_tdf, scale, md.im_converter)?;
    let reach = p
        .overlap
        .then(|| reach_in_scans(p.overlap_reach, |s| k0.convert(s as f64), num_scans));
    let im_at_scan = |s: u32| k0.convert(s as f64);
    let mz_to_tof = |mz: f64| md.mz_converter.invert(mz);
    let Some(mut gate) = PolygonGate::build(
        &mz, &im, num_scans, im_at_scan, mz_to_tof, p.mz_pad, p.im_pad,
    ) else {
        return Ok(None);
    };
    // Padding may only add points: refuse to run with a gate that would drop
    // signal inside the instrument's own selection polygon.
    gate.check_contains_unpadded(&mz, &im, im_at_scan, mz_to_tof)
        .map_err(|e| {
            DnoiseError::InvalidInput(format!(
                "{}: MS1 polygon gate self-check failed (m/z pad {} Da, 1/K0 pad {}): the \
                 padded gate does not contain the selection polygon ({e}). Pads must be >= 0; \
                 otherwise this is a dnoise bug, please report it. Disable the gate with \
                 ms1_polygon = false / --no-ms1-polygon to proceed.",
                in_tdf.display(),
                p.mz_pad,
                p.im_pad
            ))
        })?;
    gate.overlap = reach;
    Ok(Some(gate))
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

/// Recursively copy `src` into `dst`, skipping a top-level entry named `skip_top`.
fn copy_dir_except(src: &Path, dst: &Path, skip_top: &str) -> Result<()> {
    fs::create_dir_all(dst)?;
    for entry in fs::read_dir(src)? {
        let entry = entry?;
        let name = entry.file_name();
        let from = entry.path();
        let to = dst.join(&name);
        if !entry.file_type()?.is_dir() && !entry.file_type()?.is_file() {
            return Err(DnoiseError::InvalidInput(format!(
                "non-regular files (including symlinks) inside acquisitions are unsupported: {}",
                from.display()
            )));
        }
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
        if !entry.file_type()?.is_dir() && !entry.file_type()?.is_file() {
            return Err(DnoiseError::InvalidInput(format!(
                "non-regular files (including symlinks) inside acquisitions are unsupported: {}",
                from.display()
            )));
        }
        if from.is_dir() {
            copy_dir_recursive(&from, &to)?;
        } else {
            fs::copy(&from, &to)?;
        }
    }
    Ok(())
}

/// The run's `.d` calibration, exposed so an in-process caller can turn the
/// integer `(scan, tof_idx)` of a [`DecodedFrame`]'s survivors into physical
/// `(1/K0, m/z)` without re-opening the metadata. m/z uses timsrust's
/// converter; 1/K0 uses the run's [`Stages::mobility_scale`] (Bruker's
/// acquisition calibration by default, see [`crate::mobility`]).
pub struct Calibration {
    tof2mz: timsrust::converters::Tof2MzConverter,
    scan2im: ScanToMobility,
}

impl Calibration {
    /// Convert a TOF index to m/z (Da).
    pub fn tof_to_mz(&self, tof: u32) -> f64 {
        self.tof2mz.convert(tof as f64)
    }

    /// Convert a scan index to ion mobility (`1/K0`).
    pub fn scan_to_im(&self, scan: u32) -> f64 {
        self.scan2im.convert(scan as f64)
    }
}

/// In-process streaming denoiser: open a `.d` once, build the run's gates, and
/// hand back each frame's surviving points via [`RunContext::process`] — the
/// exact stage code the standalone tool runs (both go through
/// [`process_frame_decoded`]), with no rewritten `.d` on disk.
///
/// The caller drives parallelism, e.g.
/// ```ignore
/// let ctx = RunContext::open(input, &params, &stages)?;
/// let frames: Vec<_> = (0..ctx.len())
///     .into_par_iter()
///     .filter(|&i| ctx.is_ms1(i))
///     .map(|i| ctx.process(i))
///     .collect::<Result<_>>()?;
/// ```
///
/// This is the non-crop, non-dry-run path (crop/RT-crop belong to the file
/// writer). The per-run gate wiring mirrors [`run`]; the [`streaming_matches_writer`]
/// test asserts the two never diverge on a real `.d`.
pub struct RunContext<'a> {
    reader: FrameReader,
    meta: Vec<tdf::FrameMeta>,
    params: FilterParams,
    stages: Stages<'a>,
    msms_keep: Option<MsmsKeep>,
    prm_windows: Option<PrmWindows>,
    dia_msms: Option<&'a MsmsFilterParams>,
    dia_windows: Option<DiaWindows>,
    dia_regions: Option<tdf::dia::DiaRegions>,
    dda_windows: Option<DiaWindows>,
    dia_ms1_gate: Option<DiaMs1Gate>,
    polygon_gate: Option<PolygonGate>,
    rt_keep: Vec<bool>,
    neighbors: Option<NeighborIndex>,
    calibration: Calibration,
    n_ms1: usize,
}

impl<'a> RunContext<'a> {
    /// Open the `.d` and build every per-run gate once. `params` and `stages`
    /// must outlive the context (the DIA MS/MS params are borrowed from `stages`).
    pub fn open(input: &Path, params: &FilterParams, stages: &'a Stages<'a>) -> Result<Self> {
        crate::validation::parameters(params, stages, &RunOptions::default())?;
        let in_tdf = input.join("analysis.tdf");
        let in_bin = input.join("analysis.tdf_bin");
        if !in_tdf.is_file() || !in_bin.is_file() {
            return Err(DnoiseError::NotADotD(input.to_path_buf()));
        }

        let meta = tdf::read_frame_meta(&in_tdf)?;
        let acquisition = tdf::inspect_acquisition(&in_tdf, &meta)?;
        let scheme = acquisition.kind;
        let reader = FrameReader::new(input).map_err(|e| DnoiseError::OpenFrames(e.to_string()))?;
        let effective_stages = acquisition_stages(scheme, stages, false)?;
        let stages = &effective_stages;
        let n_frames = meta.len();
        let n_ms1 = meta.iter().filter(|m| m.is_ms1()).count();

        let neighbors = NeighborIndex::build(&in_tdf, &meta, scheme, stages)?;

        // MS/MS denoise: ddaPASEF per-precursor keep sets vs diaPASEF whole-frame.
        let prm_windows = (scheme == crate::Acquisition::PrmPasef && stages.denoise_msms.is_some())
            .then_some(acquisition.prm_windows);
        let (msms_keep, dia_msms) = match stages.denoise_msms {
            Some(_) if prm_windows.is_some() => (None, None),
            Some(mp) => {
                let windows = tdf::read_pasef_msms(&in_tdf)?;
                if windows.is_empty() {
                    (None, Some(mp))
                } else {
                    let keep = build_msms_keep(&reader, &meta, &windows, mp, stages.halo, None)?;
                    (Some(keep), None)
                }
            }
            None => (None, None),
        };

        // diaPASEF isolation windows, shared by the out-of-window gate and per-window
        // MS/MS filtering. Empty for ddaPASEF, so both features no-op there.
        let dia_windows = if stages.dia_window.is_some() || stages.dia_per_window {
            let w = tdf::read_dia_windows(&in_tdf)?;
            if w.is_empty() { None } else { Some(w) }
        } else {
            None
        };

        let dia_regions = if scheme == crate::Acquisition::DiaPasef
            && stages.dia_per_window
            && (stages.denoise_msms.is_some() || stages.filter_all_frames)
        {
            Some(tdf::dia::read(&in_tdf)?)
        } else {
            None
        };

        // ddaPASEF isolation-event intervals for the MS/MS out-of-window gate
        // (mirrors `run`). Empty for diaPASEF, so the gate no-ops there.
        let dda_windows = if stages.dda_window.is_some() {
            let w = DiaWindows::from_pasef(&tdf::read_pasef_msms(&in_tdf)?);
            if w.is_empty() { None } else { Some(w) }
        } else {
            None
        };

        let dia_ms1_gate = match stages.dia_ms1 {
            Some(mp) => build_dia_ms1_gate(&in_tdf, mp, &meta, stages.mobility_scale)?,
            None => None,
        };
        let polygon_gate = match stages.ms1_polygon {
            Some(pp) => build_polygon_gate(&in_tdf, pp, &meta, stages.mobility_scale)?,
            None => None,
        };

        // No crop on the streaming path: every frame is in the RT window.
        let rt_keep = vec![true; n_frames];

        let md = MetadataReader::new(&in_tdf).map_err(|e| DnoiseError::Metadata(e.to_string()))?;
        // The accessor serves the run's scale; a multi-calibration run (which the
        // gates refuse) keeps the linear scale here, as before 0.4.0.
        // The gates above already refused a run they cannot calibrate, so a
        // failure here only affects this accessor: warn and serve the linear scale
        // rather than refuse a run no gate needs a calibration for.
        let (_, im_segments) = tdf::count_calibration_segments(&in_tdf)?;
        let scan2im = if im_segments > 1 {
            warn!(
                "run has {im_segments} mobility calibrations; Calibration::scan_to_im uses the linear scale"
            );
            ScanToMobility::Linear(md.im_converter)
        } else {
            mobility::load(&in_tdf, stages.mobility_scale, md.im_converter).unwrap_or_else(|e| {
                warn!("{e}; Calibration::scan_to_im uses the linear scale");
                ScanToMobility::Linear(md.im_converter)
            })
        };
        let calibration = Calibration {
            tof2mz: md.mz_converter,
            scan2im,
        };

        Ok(RunContext {
            reader,
            meta,
            params: *params,
            stages: *stages,
            msms_keep,
            prm_windows,
            dia_msms,
            dia_windows,
            dia_regions,
            dda_windows,
            dia_ms1_gate,
            polygon_gate,
            rt_keep,
            neighbors,
            calibration,
            n_ms1,
        })
    }

    /// Total frame count (MS1 + MS/MS), the valid range for [`Self::process`].
    pub fn len(&self) -> usize {
        self.meta.len()
    }

    /// True when the context holds no frames.
    pub fn is_empty(&self) -> bool {
        self.meta.is_empty()
    }

    /// Number of MS1 frames.
    pub fn ms1_frames(&self) -> usize {
        self.n_ms1
    }

    /// True when frame `i` is an MS1 frame.
    pub fn is_ms1(&self, i: usize) -> bool {
        self.meta[i].is_ms1()
    }

    /// The run's calibration, for converting survivor `(scan, tof)` to `(1/K0, m/z)`.
    pub fn calibration(&self) -> &Calibration {
        &self.calibration
    }

    /// Decode frame `i` through every enabled stage and return its survivors.
    /// Safe to call concurrently across frames (`&self`).
    pub fn process(&self, i: usize) -> Result<DecodedFrame> {
        let ctx = FrameCtx {
            msms: self.msms_keep.as_ref(),
            prm_msms: self.stages.denoise_msms.zip(self.prm_windows.as_ref()),
            dia_msms: self.dia_msms,
            dia_windows: self.dia_windows.as_ref(),
            dia_regions: self.dia_regions.as_ref(),
            dda_windows: self.dda_windows.as_ref(),
            dia_ms1: self.dia_ms1_gate.as_ref(),
            polygon: self.polygon_gate.as_ref(),
            crop: None,
            crop_only: false,
            rt_keep: &self.rt_keep,
            neighbors: self.neighbors.as_ref(),
            cancel: None,
        };
        process_frame_decoded(
            &self.reader,
            &self.meta,
            i,
            &self.params,
            &self.stages,
            &ctx,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn benchmark_calibration() -> ScanToMobility {
        // The TimsCalibration row used by the mobility.rs load test.
        ScanToMobility::Calibrated(
            mobility::TimsCalibrationModel::new(
                2,
                [
                    1.0,
                    935.0,
                    239.34640606518187,
                    102.96946384662338,
                    33.64485981308411,
                    1.0,
                    -0.026580764926972034,
                    171.42849749723894,
                    16.838457909616054,
                    1732.6649859338625,
                ],
            )
            .unwrap(),
        )
    }

    #[test]
    fn window_scans_without_pad_is_the_half_open_window_exactly() {
        let im = benchmark_calibration();
        for begin in 0..900u32 {
            assert_eq!(
                window_scans(begin, begin + 30, 0.0, &im, 936),
                (begin, begin + 29)
            );
        }
        assert_eq!(window_scans(920, 1000, 0.0, &im, 936), (920, 935));
    }

    #[test]
    fn window_scans_tiny_pad_adds_no_scan_from_float_noise() {
        // A pad far below one scan's 1/K0 width must not widen the window.
        let im = benchmark_calibration();
        for begin in 1..900u32 {
            assert_eq!(
                window_scans(begin, begin + 30, 1e-9, &im, 936),
                (begin, begin + 29)
            );
        }
    }

    #[test]
    fn window_scans_pad_widens_both_edges() {
        let im = benchmark_calibration();
        let (lo, hi) = window_scans(400, 430, 0.01, &im, 936);
        assert!(lo < 400 && hi > 429, "({lo}, {hi})");
        assert!((im.convert(f64::from(lo)) - im.convert(400.0)).abs() <= 0.01 + 1e-3);
    }

    #[test]
    fn frame_sampled_fraction_one_keeps_every_frame() {
        assert!((0..1000).all(|i| frame_sampled(i, 42, 1.0)));
    }

    #[test]
    fn frame_sampled_fraction_zero_keeps_essentially_none() {
        let kept = (0..1000).filter(|&i| frame_sampled(i, 42, 0.0)).count();
        assert!(kept <= 1, "expected ~0 kept at fraction 0.0, got {kept}");
    }

    #[test]
    fn frame_sampled_is_deterministic() {
        assert!((0..500).all(|i| frame_sampled(i, 7, 0.5) == frame_sampled(i, 7, 0.5)));
    }

    #[test]
    fn frame_sampled_roughly_matches_the_requested_fraction() {
        let n = 20_000;
        let kept = (0..n).filter(|&i| frame_sampled(i, 123, 0.25)).count();
        let frac = kept as f64 / n as f64;
        // Loose band around 0.25; the SplitMix hash is well-distributed.
        assert!((0.22..0.28).contains(&frac), "fraction {frac} out of band");
    }

    #[test]
    fn denoise_stats_default_is_zeroed() {
        let s = DenoiseStats::default();
        assert_eq!(s.frames, 0);
        assert_eq!(s.raw_points, 0);
        assert_eq!(s.kept_points, 0);
        assert!(!s.dry_run);
    }

    #[test]
    fn copy_dir_except_skips_named_top_level_entry_and_recurses() {
        let base = std::env::temp_dir().join(format!(
            "dnoise_copy_except_{}_{}",
            std::process::id(),
            "writer"
        ));
        let src = base.join("src");
        let dst = base.join("dst");
        let _ = fs::remove_dir_all(&base);
        fs::create_dir_all(src.join("sub")).unwrap();
        fs::write(src.join("keep.txt"), b"a").unwrap();
        fs::write(src.join("analysis.tdf_bin"), b"skip me").unwrap();
        fs::write(src.join("sub").join("nested.txt"), b"b").unwrap();

        copy_dir_except(&src, &dst, "analysis.tdf_bin").unwrap();

        assert!(dst.join("keep.txt").is_file());
        assert!(dst.join("sub").join("nested.txt").is_file());
        assert!(
            !dst.join("analysis.tdf_bin").exists(),
            "the skipped top-level entry must not be copied"
        );

        let _ = fs::remove_dir_all(&base);
    }
}
