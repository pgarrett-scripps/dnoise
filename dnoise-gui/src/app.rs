//! The eframe application: input queue, output options, and run controls, with a
//! log pane. Denoising happens on a worker thread ([`crate::worker`]); the UI
//! polls a channel for progress and results.

use crate::settings::{OutputMode, PresetChoice, Settings};
use crate::worker::{WorkerMsg, run_batch};
use eframe::egui;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{Receiver, channel};
use std::time::Duration;

/// One queued input `.d` plus its detected acquisition scheme (shown in the list).
struct QueuedInput {
    path: PathBuf,
    scheme: String,
}

/// Live progress of the current batch.
struct RunState {
    n_files: usize,
    current: usize,
    done_frames: usize,
    total_frames: usize,
}

/// The whole GUI state.
pub struct DnoiseApp {
    queue: Vec<QueuedInput>,
    path_entry: String,
    settings: Settings,
    run: Option<RunState>,
    rx: Option<Receiver<WorkerMsg>>,
    cancel: Arc<AtomicBool>,
    log: Vec<String>,
}

impl Default for DnoiseApp {
    fn default() -> Self {
        Self {
            queue: Vec::new(),
            path_entry: String::new(),
            settings: Settings::default(),
            run: None,
            rx: None,
            cancel: Arc::new(AtomicBool::new(false)),
            log: Vec::new(),
        }
    }
}

impl DnoiseApp {
    /// Add one folder to the queue, validating it looks like a `.d` and detecting
    /// its acquisition scheme for the list. Duplicates and non-`.d` folders are
    /// skipped with a log note.
    fn add_path(&mut self, p: PathBuf) {
        if !p.join("analysis.tdf").is_file() {
            self.log
                .push(format!("skipped {}: not a .d folder", p.display()));
            return;
        }
        if self.queue.iter().any(|q| q.path == p) {
            return;
        }
        let scheme = dnoise::detect_acquisition(&p)
            .map(|a| format!("{a:?}"))
            .unwrap_or_else(|_| "?".to_string());
        self.queue.push(QueuedInput { path: p, scheme });
    }

    /// Spawn the worker thread for the current queue + settings.
    fn start_run(&mut self, ctx: &egui::Context) {
        self.cancel = Arc::new(AtomicBool::new(false));
        let (tx, rx) = channel();
        self.rx = Some(rx);
        self.run = Some(RunState {
            n_files: self.queue.len(),
            current: 0,
            done_frames: 0,
            total_frames: 0,
        });
        let inputs: Vec<PathBuf> = self.queue.iter().map(|q| q.path.clone()).collect();
        let settings = self.settings.clone();
        let cancel = self.cancel.clone();
        let ctx = ctx.clone();
        self.log
            .push(format!("Started: {} file(s).", self.queue.len()));
        std::thread::spawn(move || {
            run_batch(inputs, settings, tx, cancel);
            ctx.request_repaint();
        });
    }

    /// Drain any pending worker messages into the UI state.
    fn poll_worker(&mut self) {
        let Some(rx) = &self.rx else { return };
        let mut done = false;
        while let Ok(msg) = rx.try_recv() {
            match msg {
                WorkerMsg::Log(s) => self.log.push(s),
                WorkerMsg::Progress { file, done, total } => {
                    if let Some(r) = &mut self.run {
                        r.current = file;
                        r.done_frames = done;
                        r.total_frames = total;
                    }
                }
                WorkerMsg::FileDone {
                    file,
                    kept_pct,
                    out,
                } => {
                    self.log.push(format!(
                        "  [{}] done: {kept_pct:.1}% kept  ->  {}",
                        file + 1,
                        out.display()
                    ));
                }
                WorkerMsg::FileError { file, error } => {
                    self.log.push(format!("  ERROR on file {}: {error}", file + 1));
                }
                WorkerMsg::Finished => {
                    self.log.push("Batch finished.".to_string());
                    done = true;
                }
            }
        }
        if done {
            self.run = None;
            self.rx = None;
        }
    }
}

impl eframe::App for DnoiseApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Ingest any folders dropped onto the window.
        let dropped: Vec<PathBuf> =
            ctx.input(|i| i.raw.dropped_files.iter().filter_map(|f| f.path.clone()).collect());
        for p in dropped {
            self.add_path(p);
        }

        self.poll_worker();
        let running = self.run.is_some();

        egui::TopBottomPanel::top("header").show(ctx, |ui| {
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                ui.heading("dnoise");
                ui.separator();
                ui.label("Acquisition preset:");
                ui.add_enabled_ui(!running, |ui| {
                    egui::ComboBox::from_id_salt("preset")
                        .selected_text(self.settings.preset.label())
                        .show_ui(ui, |ui| {
                            for p in PresetChoice::ALL {
                                ui.selectable_value(&mut self.settings.preset, p, p.label());
                            }
                        });
                });
                ui.weak("(Auto picks the right MS1 gate per file)");
            });
            ui.add_space(4.0);
        });

        egui::TopBottomPanel::bottom("log")
            .resizable(true)
            .min_height(130.0)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.strong("Log");
                    if ui.button("Clear").clicked() {
                        self.log.clear();
                    }
                });
                egui::ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .stick_to_bottom(true)
                    .show(ui, |ui| {
                        for line in &self.log {
                            ui.monospace(line);
                        }
                    });
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            // --- Input ---
            ui.group(|ui| {
                ui.strong("Input .d folders");
                ui.label("Drag folders onto the window, or paste a path below.");
                ui.add_enabled_ui(!running, |ui| {
                    ui.horizontal(|ui| {
                        let entry = egui::TextEdit::singleline(&mut self.path_entry)
                            .hint_text("/path/to/sample.d")
                            .desired_width(360.0);
                        ui.add(entry);
                        if ui.button("Add").clicked() {
                            let t = self.path_entry.trim().to_string();
                            if !t.is_empty() {
                                self.add_path(PathBuf::from(t));
                                self.path_entry.clear();
                            }
                        }
                        if ui
                            .add_enabled(!self.queue.is_empty(), egui::Button::new("Clear all"))
                            .clicked()
                        {
                            self.queue.clear();
                        }
                    });
                });

                let mut remove: Option<usize> = None;
                egui::ScrollArea::vertical()
                    .id_salt("queue")
                    .max_height(150.0)
                    .auto_shrink([false, true])
                    .show(ui, |ui| {
                        if self.queue.is_empty() {
                            ui.weak("no folders yet");
                        }
                        for (i, item) in self.queue.iter().enumerate() {
                            ui.horizontal(|ui| {
                                if ui
                                    .add_enabled(!running, egui::Button::new("✖").small())
                                    .clicked()
                                {
                                    remove = Some(i);
                                }
                                ui.monospace(
                                    item.path.file_name().unwrap_or_default().to_string_lossy(),
                                );
                                ui.weak(format!("[{}]", item.scheme));
                            });
                        }
                    });
                if let Some(i) = remove {
                    self.queue.remove(i);
                }
            });

            ui.add_space(6.0);

            // --- Output ---
            ui.group(|ui| {
                ui.strong("Output");
                ui.add_enabled_ui(!running, |ui| {
                    ui.radio_value(
                        &mut self.settings.output_mode,
                        OutputMode::Suffix,
                        "Next to each input, with a suffix",
                    );
                    if self.settings.output_mode == OutputMode::Suffix {
                        ui.horizontal(|ui| {
                            ui.label("    suffix:");
                            ui.add(
                                egui::TextEdit::singleline(&mut self.settings.suffix)
                                    .desired_width(120.0),
                            );
                        });
                    }
                    ui.radio_value(
                        &mut self.settings.output_mode,
                        OutputMode::NewDir,
                        "Into one output folder",
                    );
                    if self.settings.output_mode == OutputMode::NewDir {
                        ui.horizontal(|ui| {
                            ui.label("    folder:");
                            ui.add(
                                egui::TextEdit::singleline(&mut self.settings.output_dir)
                                    .hint_text("/path/to/output")
                                    .desired_width(360.0),
                            );
                        });
                    }
                    ui.checkbox(
                        &mut self.settings.overwrite,
                        "Overwrite if the output already exists",
                    );
                    ui.checkbox(
                        &mut self.settings.write_report,
                        "Write a JSON report next to each output",
                    );
                });
            });

            ui.add_space(8.0);

            // --- Run / Cancel + progress ---
            ui.horizontal(|ui| {
                if running {
                    if ui.button("■  Cancel").clicked() {
                        self.cancel.store(true, Ordering::Relaxed);
                        self.log
                            .push("Cancel requested — stops after the current file.".to_string());
                    }
                } else if ui
                    .add_enabled(!self.queue.is_empty(), egui::Button::new("▶  Run"))
                    .clicked()
                {
                    self.start_run(ctx);
                }
            });

            if let Some(r) = &self.run {
                let frac = if r.total_frames > 0 {
                    r.done_frames as f32 / r.total_frames as f32
                } else {
                    0.0
                };
                ui.add(egui::ProgressBar::new(frac).text(format!(
                    "file {} of {} — {}/{} frames",
                    r.current + 1,
                    r.n_files,
                    r.done_frames,
                    r.total_frames
                )));
            }
        });

        // While a batch runs, keep polling the channel ~10x/sec.
        if self.rx.is_some() {
            ctx.request_repaint_after(Duration::from_millis(100));
        }
    }
}
