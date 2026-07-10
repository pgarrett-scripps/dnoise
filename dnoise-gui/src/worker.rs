//! Background batch runner: denoises (or, in estimate mode, dry-runs) each queued
//! `.d` in turn on a worker thread, streaming progress, results, and log lines back
//! to the UI over a channel. The core denoiser is called in-process via
//! [`dnoise::denoise_with_options`] — no subprocess, no stdout parsing.

use crate::settings::Settings;
use dnoise::{
    DenoiseStats, DiaMs1WindowParams, DiaWindowParams, Ms1PolygonParams, Progress, RunOptions,
    SampleSpec, Stages, denoise_with_options,
};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::Sender;

/// Whether the batch writes output or just estimates the reduction.
#[derive(Clone, Copy)]
pub enum RunMode {
    /// Denoise for real, writing each output `.d`.
    Full,
    /// Dry-run a deterministic frame sample for a fast reduction estimate; writes
    /// nothing.
    Estimate { fraction: f64 },
}

impl RunMode {
    fn is_estimate(self) -> bool {
        matches!(self, RunMode::Estimate { .. })
    }
}

/// A message from the worker thread to the UI.
pub enum WorkerMsg {
    /// A line for the log pane.
    Log(String),
    /// Frame progress for the file currently being processed.
    Progress {
        file: usize,
        done: usize,
        total: usize,
    },
    /// A file finished a full (writing) run successfully.
    FileDone {
        file: usize,
        kept_pct: f64,
        out: PathBuf,
    },
    /// A file finished an estimate (dry-run) pass.
    Estimate { file: usize, kept_pct: f64 },
    /// A file failed; the batch continues with the next one.
    FileError { file: usize, error: String },
    /// The whole batch is done (or was cancelled).
    Finished,
}

/// Process every input in `inputs` sequentially in the given `mode`. Checks
/// `cancel` before each file so the UI's Cancel button stops the batch after the
/// current file completes.
pub fn run_batch(
    inputs: Vec<PathBuf>,
    settings: Settings,
    mode: RunMode,
    tx: Sender<WorkerMsg>,
    cancel: Arc<AtomicBool>,
) {
    let n = inputs.len();
    for (i, input) in inputs.iter().enumerate() {
        if cancel.load(Ordering::Relaxed) {
            let _ = tx.send(WorkerMsg::Log(
                "Cancelled — remaining files skipped.".to_string(),
            ));
            break;
        }
        if let Err(e) = process_one(i, n, input, &settings, mode, &tx) {
            let _ = tx.send(WorkerMsg::FileError { file: i, error: e });
        }
    }
    let _ = tx.send(WorkerMsg::Finished);
}

/// Run one file. Returns `Err(msg)` for a setup error (bad output path, etc.);
/// denoiser errors are reported the same way.
fn process_one(
    i: usize,
    n: usize,
    input: &Path,
    settings: &Settings,
    mode: RunMode,
    tx: &Sender<WorkerMsg>,
) -> Result<(), String> {
    // A dry-run estimate writes nothing, so it needs no output path or overwrite.
    let output = if mode.is_estimate() {
        PathBuf::from("dnoise-estimate-unused")
    } else {
        let out = settings.output_path(input)?;
        if same_folder(&out, input) {
            return Err("output path equals the input; change the folder or suffix".to_string());
        }
        out
    };

    let verb = if mode.is_estimate() {
        "estimate"
    } else {
        "denoise"
    };
    let _ = tx.send(WorkerMsg::Log(format!(
        "[{}/{}] {} {}",
        i + 1,
        n,
        verb,
        input.display()
    )));

    // Warnings (bad ppm/crop input) go to the log with a small indent.
    let mut logf = |s: String| {
        let _ = tx.send(WorkerMsg::Log(format!("  {s}")));
    };

    // Build the config for this file. Owned params live on the stack for the whole
    // call, which is all `Stages`' borrows and the crop reference need.
    let params = settings.filter_params(input, &mut logf);
    let crop = settings.crop_params(&mut logf);
    let gates = settings.gates(input);
    let halo = settings.halo_params();
    let ms1_polygon = Ms1PolygonParams::default();
    let dia_ms1 = DiaMs1WindowParams::default();
    let dia_window = DiaWindowParams::default();

    let stages = Stages {
        filter_all_frames: false,
        frame_half_width: 0,
        halo: halo.as_ref(),
        denoise_msms: None,
        smooth: None,
        watershed: None,
        box_centroid: None,
        dia_window: gates.dia_window.then_some(&dia_window),
        dia_per_window: false,
        dia_ms1: gates.dia_ms1.then_some(&dia_ms1),
        ms1_polygon: gates.ms1_polygon.then_some(&ms1_polygon),
    };

    let sample = match mode {
        RunMode::Estimate { fraction } => Some(SampleSpec { fraction, seed: 0 }),
        RunMode::Full => None,
    };
    let opts = RunOptions {
        force: settings.overwrite,
        dry_run: mode.is_estimate(),
        crop: (!crop.is_empty()).then_some(&crop),
        crop_only: settings.crop_only,
        sample,
    };

    let txp = tx.clone();
    let stats = denoise_with_options(input, &output, &params, &stages, &opts, |p: Progress| {
        let _ = txp.send(WorkerMsg::Progress {
            file: i,
            done: p.frames_done,
            total: p.frames_total,
        });
    })
    .map_err(|e| e.to_string())?;

    let kept_pct = if stats.raw_points > 0 {
        100.0 * stats.kept_points as f64 / stats.raw_points as f64
    } else {
        0.0
    };

    if mode.is_estimate() {
        let _ = tx.send(WorkerMsg::Estimate { file: i, kept_pct });
    } else {
        if settings.write_report {
            if let Err(e) = write_report(&output, input, &stats) {
                let _ = tx.send(WorkerMsg::Log(format!("  report write failed: {e}")));
            }
        }
        let _ = tx.send(WorkerMsg::FileDone {
            file: i,
            kept_pct,
            out: output,
        });
    }
    Ok(())
}

/// True when two paths resolve to the same folder (best effort: falls back to a
/// literal compare when either path does not yet exist, which is the common case
/// for the not-yet-created output).
fn same_folder(a: &Path, b: &Path) -> bool {
    match (a.canonicalize(), b.canonicalize()) {
        (Ok(x), Ok(y)) => x == y,
        _ => a == b,
    }
}

/// Write a small JSON report next to the output (`<output>.report.json`): the
/// acquisition scheme and the reduction statistics.
fn write_report(output: &Path, input: &Path, stats: &DenoiseStats) -> std::io::Result<()> {
    let scheme = dnoise::detect_acquisition(input)
        .map(|a| format!("{a:?}"))
        .unwrap_or_else(|_| "unknown".to_string());
    let pct = |k: u64, r: u64| {
        if r > 0 {
            (10_000.0 * k as f64 / r as f64).round() / 100.0
        } else {
            0.0
        }
    };
    let v = serde_json::json!({
        "input": input.display().to_string(),
        "output": output.display().to_string(),
        "acquisition": scheme,
        "stats": {
            "frames": stats.frames,
            "ms1_frames": stats.ms1_frames,
            "msms_frames": stats.msms_frames,
            "raw_points": stats.raw_points,
            "kept_points": stats.kept_points,
            "kept_pct": pct(stats.kept_points, stats.raw_points),
            "raw_ms1_points": stats.raw_ms1_points,
            "kept_ms1_points": stats.kept_ms1_points,
        },
    });
    let report_path = output.with_extension("report.json");
    std::fs::write(report_path, serde_json::to_string_pretty(&v)? + "\n")
}
