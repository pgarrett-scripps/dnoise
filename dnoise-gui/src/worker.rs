//! Background batch runner: denoises (or, in estimate mode, dry-runs) each queued
//! `.d` in turn on a worker thread, streaming progress, results, and log lines back
//! to the UI over a channel. The core denoiser is called in-process via
//! [`dnoise::denoise_with_options`] — no subprocess, no stdout parsing.

use crate::settings::Settings;
use dnoise::{Progress, RunOptions, SampleSpec, denoise_with_options};
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
/// `cancel` before each file and passes it into the writer for cancellation
/// between validation and processing frames.
pub fn run_batch(
    inputs: Vec<PathBuf>,
    settings: Settings,
    mode: RunMode,
    tx: Sender<WorkerMsg>,
    cancel: Arc<AtomicBool>,
) {
    let n = inputs.len();
    if !mode.is_estimate() {
        let jobs: Result<Vec<_>, String> = inputs
            .iter()
            .map(|input| {
                Ok(dnoise::batch::Job {
                    input: input.clone(),
                    output: settings.output_path(input)?,
                    config: settings.config_for(input)?,
                })
            })
            .collect();
        let result = jobs.and_then(|jobs| {
            dnoise::batch::preflight(&jobs, settings.overwrite).map_err(|e| e.to_string())
        });
        if let Err(error) = result {
            for file in 0..n {
                let _ = tx.send(WorkerMsg::FileError {
                    file,
                    error: format!("batch preflight: {error}"),
                });
            }
            let _ = tx.send(WorkerMsg::Finished);
            return;
        }
    }

    for (i, input) in inputs.iter().enumerate() {
        if cancel.load(Ordering::Relaxed) {
            let _ = tx.send(WorkerMsg::Log(
                "Cancelled — remaining files skipped.".to_string(),
            ));
            break;
        }
        if let Err(e) = process_one(i, n, input, &settings, mode, &tx, &cancel) {
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
    cancel: &AtomicBool,
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

    let config = settings.config_for(input)?;
    let resolved = config.resolve(input).map_err(|e| e.to_string())?;
    let params = resolved.filter;
    let crop = resolved.crop;
    let stages = resolved.stages();
    if dnoise::provenance::read(input)
        .map_err(|e| e.to_string())?
        .is_some()
    {
        let _ = tx.send(WorkerMsg::Log(
            "Warning: input already has dnoise processing history.".into(),
        ));
    }
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
        cancel: Some(cancel),
        frame_batch_size: config.frame_batch_size,
        skip_validation: config.skip_validation.unwrap_or(false),
    };

    let txp = tx.clone();
    let stats = match config.in_thread_pool(|| {
        denoise_with_options(input, &output, &params, &stages, &opts, |p: Progress| {
            let _ = txp.send(WorkerMsg::Progress {
                file: i,
                done: p.frames_done,
                total: p.frames_total,
            });
        })
    }) {
        Ok(s) => s,
        Err(e) => {
            if matches!(e, dnoise::DnoiseError::Cancelled) {
                let _ = tx.send(WorkerMsg::Log(
                    "Cancelled; existing files preserved.".into(),
                ));
                return Ok(());
            }
            return Err(e.to_string());
        }
    };

    let kept_pct = if stats.raw_points > 0 {
        100.0 * stats.kept_points as f64 / stats.raw_points as f64
    } else {
        0.0
    };

    for warning in dnoise::provenance::report(input, &params, &stages, &opts, &stats)["warnings"]
        .as_array()
        .into_iter()
        .flatten()
    {
        let _ = tx.send(WorkerMsg::Log(format!(
            "Warning: {}",
            warning.as_str().unwrap_or_default()
        )));
    }
    if mode.is_estimate() {
        let _ = tx.send(WorkerMsg::Estimate { file: i, kept_pct });
    } else {
        if settings.write_report {
            let report = dnoise::provenance::report(input, &params, &stages, &opts, &stats);
            if let Err(e) = write_report(&output, &report) {
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
fn write_report(output: &Path, report: &serde_json::Value) -> std::io::Result<()> {
    use std::io::Write;
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output.with_extension("report.json"))?;
    file.write_all((serde_json::to_string_pretty(report)? + "\n").as_bytes())
}
