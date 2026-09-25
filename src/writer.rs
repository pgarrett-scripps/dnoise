//! Orchestration: copy the source `.d`, rewrite `analysis.tdf_bin` with filtered
//! frames (re-encoded as type 2), and fix up the `analysis.tdf` SQLite database.

use crate::codec::try_encode_frame_type2;
use crate::crop::CropGate;
use crate::dia_ms1::DiaMs1Gate;
use crate::error::{DnoiseError, Result};
use crate::frame::FlatFrame;
use crate::mobility::{self, MobilityScale, ScanToMobility};
use crate::msms::build_msms_keep;
use crate::neighbor;
use crate::params::{CropParams, DiaMs1WindowParams, FilterParams, Ms1PolygonParams, Stages};
use crate::polygon::PolygonGate;
use crate::provenance::NeighborUsage;
use crate::tdf::{self, DiaWindows, FrameUpdate};
use crate::tsr::ConvertableDomain;
use crate::tsr::{FrameReader, MetadataReader};
use dnoise_core::{DecodedFrame, Denoiser};
use rayon::prelude::*;
use std::fs;
use std::io::{BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use tracing::{debug, info, warn};

/// Frames are read+filtered+encoded in parallel batches of this size, then the
/// batch is written sequentially (so offsets stay ordered) before the next.
/// Bounds peak memory to roughly this many encoded frames.
const CHUNK: usize = 2048;

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
    // The per-frame pipeline. It applies the acquisition policy to `stages`;
    // everything below is built from the effective stages it reports.
    let mut denoiser = Denoiser::new(params, stages, scheme, meta.clone(), crop_only)?;
    let effective_stages = *denoiser.stages();
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
        neighbor::build(&in_tdf, &meta, scheme, stages)?
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
    // Raw-unit parameters (TOF indices, scans) with their physical equivalents
    // for this run's calibration. Informational only; a run whose calibration
    // cannot be read here still runs (a gate that needs it fails on its own).
    if !crop_only {
        match crate::units::for_run(&in_tdf, params, stages) {
            Ok(eqs) => {
                for e in eqs {
                    info!("units: {}", e.line());
                }
            }
            Err(e) => debug!("units: physical equivalents unavailable ({e})"),
        }
    }

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

    // Hand the run state to the per-frame pipeline. Frames outside the RT crop
    // are emitted empty rather than deleted, so the frame axis stays valid.
    denoiser.set_neighbors(neighbors);
    denoiser.set_msms_keep(msms_keep);
    denoiser.set_prm_windows(prm_windows);
    denoiser.set_whole_frame_msms(dia_msms.is_some());
    denoiser.set_dia_windows(dia_windows);
    denoiser.set_dia_regions(dia_regions);
    denoiser.set_dda_windows(dda_windows);
    denoiser.set_dia_ms1_gate(dia_ms1_gate);
    denoiser.set_polygon_gate(polygon_gate);
    denoiser.set_crop_gate(crop_gate);
    denoiser.set_rt_crop(crop);
    let denoiser = denoiser;

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
                process_frame(&reader, &denoiser, i, cancel)
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
            ms1_neighbors: denoiser.neighbors().is_some()
                && stages.frame_half_width > 0
                && n_ms1 > 0,
            prm_neighbors: denoiser.neighbors().is_some()
                && stages.neighbors.prm_radius > 0
                && scheme == crate::Acquisition::PrmPasef,
            dia_neighbors: denoiser.neighbors().is_some()
                && stages.neighbors.dia_radius > 0
                && scheme == crate::Acquisition::DiaPasef,
            dia_scan_varying: denoiser.dia_regions().is_some_and(|r| r.has_scanning()),
            prm_per_event: denoiser.prm_windows().is_some() && !crop_only,
            ms1_polygon: denoiser.polygon_gate().is_some() && !crop_only,
            dia_ms1: denoiser.dia_ms1_gate().is_some() && !crop_only,
            dia_window: denoiser.dia_windows().is_some()
                && dia_window.is_some()
                && !crop_only
                && (stages.filter_all_frames || denoise_msms.is_some()),
            dda_window: denoiser.dda_windows().is_some()
                && !crop_only
                && (stages.filter_all_frames || denoise_msms.is_some()),
            dia_per_window: denoiser.dia_windows().is_some()
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

/// Read 0-based frame `j` for the per-frame pipeline.
fn read_flat(reader: &FrameReader, j: usize) -> dnoise_core::Result<FlatFrame> {
    reader
        .get(j)
        .map(|f| FlatFrame::from_frame(&f))
        .map_err(|e| dnoise_core::Error::FrameRead {
            index: j,
            message: e.to_string(),
        })
}

/// Thin file-writer wrapper over [`Denoiser::process`]: decode the frame's
/// survivors, then encode them into a `.d` type-2 record plus the per-frame stats
/// the run's metadata fixup needs. Behaviour-identical to the pre-streaming
/// implementation (proven by the byte-identity test in `tests/`).
fn process_frame(
    reader: &FrameReader,
    denoiser: &Denoiser,
    i: usize,
    cancel: Option<&AtomicBool>,
) -> Result<ProcessedFrame> {
    let d = denoiser.process(i, &|j| read_flat(reader, j), cancel)?;
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
/// in physical units (`mz_pad` Th, `im_pad` 1/K0) using the run's calibration,
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
    let mz_to_tof = |mz: f64| md.mz_converter.invert(mz);
    Ok(DiaMs1Gate::from_windows(
        &boxes, p, mz_to_tof, &im, num_scans,
    ))
}

/// Build the MS1 selection-polygon gate: read the run's IMS PolygonFilter
/// `(m/z, 1/K0)` vertices, convert them to per-scan TOF-index intervals via the
/// run calibration (padded by `mz_pad` Th / `im_pad` 1/K0), and assemble the
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
                "{}: MS1 polygon gate self-check failed (m/z pad {} Th, 1/K0 pad {}): the \
                 padded gate does not contain the selection polygon ({e}). Pads must be >= 0; \
                 otherwise this is a dnoise bug, please report it. Disable the gate with \
                 ms1_polygon = false / --no-ms1-polygon to proceed.",
                in_tdf.display(),
                p.mz_pad,
                p.im_pad
            ))
        })?;
    gate.overlap = p.overlap;
    Ok(Some(gate))
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
    tof2mz: crate::tsr::Tof2MzConverter,
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
/// [`Denoiser::process`]), with no rewritten `.d` on disk.
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
/// writer). The per-run gate wiring mirrors [`run`]; the `streaming_matches_writer`
/// test asserts the two never diverge on a real `.d`.
pub struct RunContext<'a> {
    reader: FrameReader,
    denoiser: Denoiser<'a>,
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
        let mut denoiser = Denoiser::new(params, stages, scheme, meta.clone(), false)?;
        let effective_stages = *denoiser.stages();
        let stages = &effective_stages;
        let n_ms1 = meta.iter().filter(|m| m.is_ms1()).count();

        let neighbors = neighbor::build(&in_tdf, &meta, scheme, stages)?;

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
        denoiser.set_neighbors(neighbors);
        denoiser.set_msms_keep(msms_keep);
        denoiser.set_prm_windows(prm_windows);
        denoiser.set_whole_frame_msms(dia_msms.is_some());
        denoiser.set_dia_windows(dia_windows);
        denoiser.set_dia_regions(dia_regions);
        denoiser.set_dda_windows(dda_windows);
        denoiser.set_dia_ms1_gate(dia_ms1_gate);
        denoiser.set_polygon_gate(polygon_gate);

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
            denoiser,
            calibration,
            n_ms1,
        })
    }

    /// Total frame count (MS1 + MS/MS), the valid range for [`Self::process`].
    pub fn len(&self) -> usize {
        self.denoiser.len()
    }

    /// True when the context holds no frames.
    pub fn is_empty(&self) -> bool {
        self.denoiser.is_empty()
    }

    /// Number of MS1 frames.
    pub fn ms1_frames(&self) -> usize {
        self.n_ms1
    }

    /// True when frame `i` is an MS1 frame.
    pub fn is_ms1(&self, i: usize) -> bool {
        self.denoiser.is_ms1(i)
    }

    /// The run's calibration, for converting survivor `(scan, tof)` to `(1/K0, m/z)`.
    pub fn calibration(&self) -> &Calibration {
        &self.calibration
    }

    /// Decode frame `i` through every enabled stage and return its survivors.
    /// Safe to call concurrently across frames (`&self`).
    pub fn process(&self, i: usize) -> Result<DecodedFrame> {
        Ok(self
            .denoiser
            .process(i, &|j| read_flat(&self.reader, j), None)?)
    }

    /// The run's per-frame pipeline, for callers that read frames themselves.
    pub fn denoiser(&self) -> &Denoiser<'a> {
        &self.denoiser
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
