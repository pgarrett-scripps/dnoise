//! The eframe application: input queue, output options, and run controls, with a
//! log pane. Denoising happens on a worker thread ([`crate::worker`]); the UI
//! polls a channel for progress and results.

use crate::settings::{OutputMode, PresetChoice, Settings};
use crate::worker::{RunMode, WorkerMsg, run_batch};
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
    /// Verb shown on the progress bar ("Denoising" / "Estimating").
    verb: &'static str,
}

/// The whole GUI state.
pub struct DnoiseApp {
    queue: Vec<QueuedInput>,
    path_entry: String,
    config_path: String,
    settings: Settings,
    run: Option<RunState>,
    rx: Option<Receiver<WorkerMsg>>,
    cancel: Arc<AtomicBool>,
    log: Vec<String>,
    completed: std::collections::HashMap<PathBuf, PathBuf>,
    active_inputs: Vec<PathBuf>,
    has_run: bool,
    preview: Option<crate::preview::Preview>,
}

impl Default for DnoiseApp {
    fn default() -> Self {
        Self {
            queue: Vec::new(),
            path_entry: String::new(),
            config_path: String::new(),
            settings: Settings::default(),
            run: None,
            rx: None,
            cancel: Arc::new(AtomicBool::new(false)),
            log: Vec::new(),
            completed: Default::default(),
            active_inputs: Vec::new(),
            has_run: false,
            preview: None,
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
        let scheme = match dnoise::detect_acquisition(&p) {
            Ok(a) => a.to_string(),
            Err(e) => {
                self.log.push(format!("skipped {}: {e}", p.display()));
                return;
            }
        };
        self.queue.push(QueuedInput { path: p, scheme });
    }

    /// Load a `dnoise.toml` into the advanced settings.
    fn load_config(&mut self, path: &std::path::Path) {
        match dnoise::config::Config::load(path) {
            Ok(c) => {
                self.settings.apply_config(&c);
                self.config_path = path.display().to_string();
                self.log
                    .push(format!("loaded settings from {}", path.display()));
            }
            Err(e) => self.log.push(format!("load failed: {e}")),
        }
    }

    /// Save the current advanced settings to a `dnoise.toml`.
    fn save_config(&mut self, path: &std::path::Path) {
        for item in &self.queue {
            if let Err(e) = dnoise::output::check_disjoint(&item.path, path) {
                self.log.push(e.to_string());
                return;
            }
        }
        match self.settings.checked_config().and_then(|c| c.save(path)) {
            Ok(()) => self
                .log
                .push(format!("saved settings to {}", path.display())),
            Err(e) => self.log.push(format!("save failed: {e}")),
        }
    }

    /// Spawn the worker thread for the current queue + settings in the given mode.
    fn start_run(&mut self, ctx: &egui::Context, mode: RunMode) {
        if let Err(e) = self.settings.checked_config() {
            self.log.push(format!("Invalid settings: {e}"));
            return;
        }
        if matches!(mode, RunMode::Full) {
            self.has_run = true;
            for item in &self.queue {
                self.completed.remove(&item.path);
            }
        }
        self.cancel = Arc::new(AtomicBool::new(false));
        let (tx, rx) = channel();
        self.rx = Some(rx);
        let verb = match mode {
            RunMode::Estimate { .. } => "Estimating",
            RunMode::Full => "Denoising",
        };
        self.run = Some(RunState {
            n_files: self.queue.len(),
            current: 0,
            done_frames: 0,
            total_frames: 0,
            verb,
        });
        let inputs: Vec<PathBuf> = self.queue.iter().map(|q| q.path.clone()).collect();
        self.active_inputs = inputs.clone();
        let settings = self.settings.clone();
        let cancel = self.cancel.clone();
        let ctx = ctx.clone();
        self.log
            .push(format!("{verb}: {} file(s)…", self.queue.len()));
        std::thread::spawn(move || {
            run_batch(inputs, settings, mode, tx, cancel);
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
                    if let Some(input) = self.active_inputs.get(file) {
                        self.completed.insert(input.clone(), out.clone());
                    }
                    self.log.push(format!(
                        "  [{}] done: {kept_pct:.1}% kept  ->  {}",
                        file + 1,
                        out.display()
                    ));
                }
                WorkerMsg::Estimate { file, kept_pct } => {
                    self.log.push(format!(
                        "  [{}] estimate: ~{kept_pct:.1}% kept (frame sample)",
                        file + 1
                    ));
                }
                WorkerMsg::FileError { file, error } => {
                    self.log
                        .push(format!("  ERROR on file {}: {error}", file + 1));
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

    /// The collapsible advanced-settings panel: every filter / halo / gate / crop /
    /// ppm knob. Defaults match the tuned CLI defaults, so leaving it closed
    /// reproduces the standard run.
    fn advanced_panel(&mut self, ui: &mut egui::Ui) {
        let s = &mut self.settings;
        ui.label("Loaded TOML options are preserved, including MS/MS, smoothing, and centroiding.");
        let msms = s.loaded_config.denoise_msms.get_or_insert(false);
        ui.checkbox(msms, "Also denoise MS/MS (can change IDs and quantities)")
            .on_hover_text("PRM uses experimental filtering within each isolation event. Validate downstream quantification before routine use.");
        ui.collapsing("Neighbor support (experimental)", |ui| {
            ui.label("1 uses the previous and next matching observation for filtering; output intensities stay native.");
            let mut ms1 = s.loaded_config.ms1_neighbor_radius.or(s.loaded_config.frame_half_width).unwrap_or(0);
            ui.horizontal(|ui| {
                ui.label("MS1 radius");
                if ui.add(egui::DragValue::new(&mut ms1).range(0..=100)).changed() {
                    s.loaded_config.frame_half_width = None;
                    s.loaded_config.ms1_neighbor_radius = Some(ms1);
                }
            });
            for (label, value) in [
                ("PRM target radius", &mut s.loaded_config.prm_neighbor_radius),
                ("DIA window radius", &mut s.loaded_config.dia_neighbor_radius),
            ] {
                ui.horizontal(|ui| {
                    ui.label(label);
                    ui.add(egui::DragValue::new(value.get_or_insert(0)).range(0..=100));
                });
            }
            ui.horizontal(|ui| {
                ui.label("Maximum RT distance (seconds)");
                ui.add(egui::DragValue::new(s.loaded_config.neighbor_max_rt_gap.get_or_insert(5.0)).range(0.01..=3600.0).speed(0.1));
            });
            ui.label("PRM/DIA radii require MS/MS denoising. Quantitative accuracy is unvalidated.");
        });
        let mut full_validation = !s.loaded_config.skip_validation.unwrap_or(false);
        if ui.checkbox(&mut full_validation, "Full input/output validation")
            .on_hover_text("Decode and check every frame before and after processing. Turning this off keeps structural and file-safety checks.")
            .changed()
        {
            s.loaded_config.skip_validation = Some(!full_validation);
        }
        let batch = s.loaded_config.frame_batch_size.get_or_insert(2048);
        ui.horizontal(|ui| {
            ui.label("Frames per batch:");
            ui.add(egui::DragValue::new(batch).range(1..=2048));
        });

        // Vertical filter (the core denoiser aggressiveness).
        ui.strong("Vertical filter");
        egui::Grid::new("vfilter").num_columns(2).show(ui, |ui| {
            ui.label("min feature length (scans)");
            ui.add(egui::DragValue::new(&mut s.min_feature_length).range(1..=100));
            ui.end_row();
            ui.label("max internal gap (scans)");
            ui.add(egui::DragValue::new(&mut s.max_internal_gap).range(0..=50));
            ui.end_row();
            ui.label("iterations");
            ui.add(egui::DragValue::new(&mut s.iterations).range(1..=10));
            ui.end_row();
            ui.label("min per-scan intensity");
            ui.add(egui::DragValue::new(&mut s.min_window_intensity).speed(1.0));
            ui.end_row();
            ui.label("min feature intensity");
            ui.add(egui::DragValue::new(&mut s.min_feature_intensity).speed(1.0));
            ui.end_row();
        });

        // m/z window: raw TOF half-width, or a ppm tolerance.
        ui.horizontal(|ui| {
            ui.checkbox(&mut s.use_ppm, "m/z window from ppm");
            if s.use_ppm {
                ui.add(
                    egui::DragValue::new(&mut s.mz_ppm)
                        .range(1.0..=100.0)
                        .suffix(" ppm"),
                );
            } else {
                ui.label("half-width");
                ui.add(egui::DragValue::new(&mut s.mz_half_width).range(1..=50));
                ui.weak("TOF idx");
            }
        });

        ui.separator();

        // Halo filter.
        ui.checkbox(&mut s.halo, "Horizontal-halo filter");
        if s.halo {
            egui::Grid::new("halo").num_columns(2).show(ui, |ui| {
                ui.label("peak fraction");
                ui.add(
                    egui::DragValue::new(&mut s.halo_peak_fraction)
                        .range(0.0..=1.0)
                        .speed(0.01),
                );
                ui.end_row();
                ui.label("m/z-idx half-width");
                ui.add(egui::DragValue::new(&mut s.halo_mz_idx_half_width).range(1..=500));
                ui.end_row();
                ui.label("scan half-width");
                ui.add(egui::DragValue::new(&mut s.halo_scan_half_width).range(0..=50));
                ui.end_row();
            });
        }

        ui.separator();

        // Gate override (else the preset decides).
        ui.checkbox(&mut s.follow_preset_gates, "Gates follow the preset");
        ui.add_enabled_ui(!s.follow_preset_gates, |ui| {
            ui.horizontal(|ui| {
                ui.checkbox(&mut s.ms1_polygon, "MS1 polygon");
                ui.checkbox(&mut s.dia_ms1_window, "DIA MS1 window");
                ui.checkbox(&mut s.dia_window, "DIA MS/MS window");
            });
        });

        ui.separator();

        // Region-of-interest crop.
        ui.strong("Crop / trim (blank = no bound)");
        let field = |ui: &mut egui::Ui, label: &str, v: &mut String, hint: &str| {
            ui.label(label);
            ui.add(
                egui::TextEdit::singleline(v)
                    .desired_width(80.0)
                    .hint_text(hint),
            );
        };
        egui::Grid::new("crop").num_columns(4).show(ui, |ui| {
            field(ui, "m/z min", &mut s.mz_min, "Da");
            field(ui, "m/z max", &mut s.mz_max, "Da");
            ui.end_row();
            field(ui, "1/K0 min", &mut s.im_min, "");
            field(ui, "1/K0 max", &mut s.im_max, "");
            ui.end_row();
            field(ui, "RT min", &mut s.rt_min, "min");
            field(ui, "RT max", &mut s.rt_max, "min");
            ui.end_row();
            field(ui, "intensity min", &mut s.min_intensity, "");
            field(ui, "intensity max", &mut s.max_intensity, "");
            ui.end_row();
        });
        ui.checkbox(&mut s.crop_only, "Crop only (skip denoising)");
    }
}

impl eframe::App for DnoiseApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let context = ui.ctx().clone();
        let ctx = &context;
        // Ingest dropped items: a `.toml` loads settings, anything else is queued as
        // an input `.d`.
        let dropped: Vec<PathBuf> = ctx.input(|i| {
            i.raw
                .dropped_files
                .iter()
                .map(|f| f.path().to_owned())
                .collect()
        });
        for p in dropped {
            if p.extension()
                .is_some_and(|e| e.eq_ignore_ascii_case("toml"))
            {
                self.load_config(&p);
            } else {
                self.add_path(p);
            }
        }

        self.poll_worker();
        let running = self.run.is_some();

        egui::Panel::top("header").show(ui, |ui| {
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

        egui::Panel::bottom("log")
            .resizable(true)
            .min_size(130.0)
            .show(ui, |ui| {
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

        egui::CentralPanel::default().show(ui, |ui| {
            egui::ScrollArea::vertical()
                .id_salt("main-scroll")
                .show(ui, |ui| {
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
                                if ui.button("Browse…").clicked()
                                    && let Some(paths) = rfd::FileDialog::new().pick_folders()
                                {
                                    for path in paths {
                                        self.add_path(path);
                                    }
                                }
                                if ui.button("Add").clicked() {
                                    let t = self.path_entry.trim().to_string();
                                    if !t.is_empty() {
                                        self.add_path(PathBuf::from(t));
                                        self.path_entry.clear();
                                    }
                                }
                                if ui
                                    .add_enabled(
                                        !self.queue.is_empty(),
                                        egui::Button::new("Clear all"),
                                    )
                                    .clicked()
                                {
                                    self.queue.clear();
                                }
                            });
                        });

                        let mut remove: Option<usize> = None;
                        let mut preview_pair = None;
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
                                            item.path
                                                .file_name()
                                                .unwrap_or_default()
                                                .to_string_lossy(),
                                        );
                                        ui.weak(format!("[{}]", item.scheme));
                                        if let Some(output) = self.completed.get(&item.path)
                                            && ui
                                                .add_enabled(!running, egui::Button::new("Inspect"))
                                                .clicked()
                                        {
                                            preview_pair =
                                                Some((item.path.clone(), output.clone()));
                                        }
                                    });
                                }
                            });
                        if let Some((input, output)) = preview_pair {
                            self.preview = Some(crate::preview::Preview::new(input, output));
                        }
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
                                            .desired_width(300.0),
                                    );
                                    if ui.button("Browse…").clicked()
                                        && let Some(path) = rfd::FileDialog::new().pick_folder()
                                    {
                                        self.settings.output_dir = path.display().to_string();
                                    }
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

                    ui.add_space(6.0);

                    // --- Advanced ---
                    ui.add_enabled_ui(!running, |ui| {
                        egui::CollapsingHeader::new("Advanced settings")
                            .default_open(false)
                            .show(ui, |ui| self.advanced_panel(ui));
                    });

                    // Save / load the advanced knobs as a CLI-compatible dnoise.toml.
                    ui.add_space(4.0);
                    ui.horizontal(|ui| {
                        ui.label("Settings file:");
                        ui.add(
                            egui::TextEdit::singleline(&mut self.config_path)
                                .hint_text("/path/to/dnoise.toml")
                                .desired_width(300.0),
                        );
                        let has = !self.config_path.trim().is_empty();
                        if ui
                            .add_enabled(has && !running, egui::Button::new("Save"))
                            .clicked()
                        {
                            let p = PathBuf::from(self.config_path.trim());
                            self.save_config(&p);
                        }
                        if ui
                            .add_enabled(has && !running, egui::Button::new("Load"))
                            .clicked()
                        {
                            let p = PathBuf::from(self.config_path.trim());
                            self.load_config(&p);
                        }
                    });
                    ui.weak("Tip: drag a .toml onto the window to load it.");

                    ui.add_space(8.0);

                    // --- Run / Estimate / Cancel + progress ---
                    ui.horizontal(|ui| {
                        if running {
                            if ui.button("■  Cancel").clicked() {
                                self.cancel.store(true, Ordering::Relaxed);
                                self.log.push("Cancel requested — stopping…".to_string());
                            }
                        } else {
                            let have = !self.queue.is_empty();
                            let retry = self
                                .queue
                                .iter()
                                .any(|q| !self.completed.contains_key(&q.path))
                                && self.has_run;
                            if ui
                                .add_enabled(retry, egui::Button::new("Retry unfinished"))
                                .clicked()
                            {
                                self.queue.retain(|q| !self.completed.contains_key(&q.path));
                                self.start_run(ctx, RunMode::Full);
                            }
                            if ui
                                .add_enabled(have, egui::Button::new("Save batch…"))
                                .clicked()
                                && let Some(path) = rfd::FileDialog::new()
                                    .set_file_name("dnoise-batch.json")
                                    .save_file()
                            {
                                let result = self
                                    .queue
                                    .iter()
                                    .map(|q| {
                                        Ok(dnoise::batch::Job {
                                            input: dnoise::output::resolved_path(&q.path)
                                                .map_err(|e| e.to_string())?,
                                            output: dnoise::output::resolved_path(
                                                &self.settings.output_path(&q.path)?,
                                            )
                                            .map_err(|e| e.to_string())?,
                                            config: self.settings.config_for(&q.path)?,
                                        })
                                    })
                                    .collect::<Result<Vec<_>, String>>()
                                    .and_then(|jobs| {
                                        for job in &jobs {
                                            dnoise::output::check_disjoint(&job.input, &path)
                                                .map_err(|e| e.to_string())?;
                                            dnoise::output::check_disjoint(&job.output, &path)
                                                .map_err(|e| e.to_string())?;
                                        }
                                        let manifest = dnoise::batch::Manifest {
                                            schema_version: 1,
                                            jobs,
                                        };
                                        let text = serde_json::to_string_pretty(&manifest)
                                            .map_err(|e| e.to_string())?;
                                        std::fs::write(&path, text).map_err(|e| e.to_string())
                                    });
                                match result {
                                    Ok(()) => self.log.push(format!("Saved {}", path.display())),
                                    Err(e) => self.log.push(e),
                                }
                            }

                            if ui.add_enabled(have, egui::Button::new("▶  Run")).clicked() {
                                self.start_run(ctx, RunMode::Full);
                            }
                            if ui
                                .add_enabled(have, egui::Button::new("⚡ Estimate reduction"))
                                .on_hover_text("Dry-run an 8% frame sample — no output written")
                                .clicked()
                            {
                                self.start_run(ctx, RunMode::Estimate { fraction: 0.08 });
                            }
                        }
                    });

                    if let Some(r) = &self.run {
                        let frac = if r.total_frames > 0 {
                            r.done_frames as f32 / r.total_frames as f32
                        } else {
                            0.0
                        };
                        ui.add(egui::ProgressBar::new(frac).text(format!(
                            "{} file {} of {} — {}/{} frames",
                            r.verb,
                            r.current + 1,
                            r.n_files,
                            r.done_frames,
                            r.total_frames
                        )));
                    }
                });
        });

        if let Some(preview) = &mut self.preview {
            preview.show(ctx);
        }
        // While a batch runs, keep polling the channel ~10x/sec.
        if self.rx.is_some() {
            ctx.request_repaint_after(Duration::from_millis(100));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eframe::App;
    #[test]
    fn main_view_renders_with_new_framework_without_a_native_window() {
        let context = egui::Context::default();
        let mut app = DnoiseApp::default();
        let mut frame = eframe::Frame::_new_kittest();
        let input = egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(900., 800.),
            )),
            ..Default::default()
        };
        let mut output = context.run_ui(input, |ui| app.ui(ui, &mut frame));
        assert!(!output.shapes.is_empty());
        output.textures_delta.clear();
    }
}
