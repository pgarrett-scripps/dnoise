//! Background batch runner: denoises each queued `.d` in turn on a worker thread,
//! streaming progress, per-file results, and log lines back to the UI over a
//! channel. The core denoiser is called in-process — no subprocess, no stdout
//! parsing — via [`dnoise::denoise_with_options`].

use crate::settings::{Settings, resolve_gates};
use dnoise::{
    DenoiseStats, DiaMs1WindowParams, DiaWindowParams, FilterParams, HaloParams, Ms1PolygonParams,
    Progress, RunOptions, Stages, denoise_with_options,
};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::Sender;

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
    /// A file finished successfully.
    FileDone {
        file: usize,
        kept_pct: f64,
        out: PathBuf,
    },
    /// A file failed; the batch continues with the next one.
    FileError { file: usize, error: String },
    /// The whole batch is done (or was cancelled).
    Finished,
}

/// Denoise every input in `inputs` sequentially. Checks `cancel` before each file
/// so the UI's Cancel button stops the batch after the current file completes.
pub fn run_batch(
    inputs: Vec<PathBuf>,
    settings: Settings,
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

        let output = match settings.output_path(input) {
            Ok(o) => o,
            Err(e) => {
                let _ = tx.send(WorkerMsg::FileError { file: i, error: e });
                continue;
            }
        };
        // Never write over the input folder (a forced run would delete it first).
        if same_folder(&output, input) {
            let _ = tx.send(WorkerMsg::FileError {
                file: i,
                error: "output path equals the input; change the folder or suffix".to_string(),
            });
            continue;
        }
        let _ = tx.send(WorkerMsg::Log(format!(
            "[{}/{}] {}  ->  {}",
            i + 1,
            n,
            input.display(),
            output.display()
        )));

        // Build the filter config for this file. Params live on the stack for the
        // duration of the call, which is all `Stages`' borrows need.
        let gates = resolve_gates(settings.preset, input);
        let params = FilterParams::default();
        let halo = HaloParams::default();
        let ms1_polygon = Ms1PolygonParams::default();
        let dia_ms1 = DiaMs1WindowParams::default();
        let dia_window = DiaWindowParams::default();
        let stages = Stages {
            filter_all_frames: false,
            frame_half_width: 0,
            halo: Some(&halo),
            denoise_msms: None,
            smooth: None,
            watershed: None,
            box_centroid: None,
            dia_window: gates.dia_window.then_some(&dia_window),
            dia_per_window: false,
            dia_ms1: gates.dia_ms1.then_some(&dia_ms1),
            ms1_polygon: gates.ms1_polygon.then_some(&ms1_polygon),
        };
        let opts = RunOptions {
            force: settings.overwrite,
            dry_run: false,
            crop: None,
            crop_only: false,
            sample: None,
        };

        let txp = tx.clone();
        let result = denoise_with_options(input, &output, &params, &stages, &opts, |p: Progress| {
            let _ = txp.send(WorkerMsg::Progress {
                file: i,
                done: p.frames_done,
                total: p.frames_total,
            });
        });

        match result {
            Ok(stats) => {
                let kept_pct = if stats.raw_points > 0 {
                    100.0 * stats.kept_points as f64 / stats.raw_points as f64
                } else {
                    0.0
                };
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
            Err(e) => {
                let _ = tx.send(WorkerMsg::FileError {
                    file: i,
                    error: e.to_string(),
                });
            }
        }
    }
    let _ = tx.send(WorkerMsg::Finished);
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
