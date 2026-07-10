//! Run settings shared between the UI and the worker thread, plus the small
//! preset -> gate and output-path helpers.

use dnoise::{Acquisition, CropParams, FilterParams, HaloParams};
use std::path::{Path, PathBuf};

/// Acquisition preset chosen in the UI. Mirrors the CLI `--preset`: it only
/// decides which MS1 gates are enabled (see [`resolve_gates`]).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum PresetChoice {
    /// Enable no gates.
    None,
    /// Detect ddaPASEF/diaPASEF per file and enable the matching gate(s).
    Auto,
    /// ddaPASEF: MS1 selection-polygon gate.
    Dda,
    /// diaPASEF: MS1 + MS/MS isolation-window gates.
    Dia,
}

impl PresetChoice {
    /// Every choice, in menu order.
    pub const ALL: [PresetChoice; 4] = [
        PresetChoice::Auto,
        PresetChoice::Dda,
        PresetChoice::Dia,
        PresetChoice::None,
    ];

    /// Human-readable menu label.
    pub fn label(self) -> &'static str {
        match self {
            PresetChoice::None => "None (no gates)",
            PresetChoice::Auto => "Auto-detect",
            PresetChoice::Dda => "ddaPASEF (polygon)",
            PresetChoice::Dia => "diaPASEF (windows)",
        }
    }
}

/// Which MS1 gates a resolved preset turns on.
pub struct Gates {
    pub ms1_polygon: bool,
    pub dia_ms1: bool,
    pub dia_window: bool,
}

/// Resolve a preset to concrete gate enables. `Auto` inspects the input `.d`
/// (ddaPASEF -> polygon, diaPASEF -> windows, anything else -> nothing), exactly
/// like the CLI's `--preset auto`.
pub fn resolve_gates(preset: PresetChoice, input: &Path) -> Gates {
    let effective = match preset {
        PresetChoice::Auto => match dnoise::detect_acquisition(input) {
            Ok(Acquisition::DdaPasef) => PresetChoice::Dda,
            Ok(Acquisition::DiaPasef) => PresetChoice::Dia,
            _ => PresetChoice::None,
        },
        other => other,
    };
    match effective {
        PresetChoice::Dda => Gates {
            ms1_polygon: true,
            dia_ms1: false,
            dia_window: false,
        },
        PresetChoice::Dia => Gates {
            ms1_polygon: false,
            dia_ms1: true,
            dia_window: true,
        },
        _ => Gates {
            ms1_polygon: false,
            dia_ms1: false,
            dia_window: false,
        },
    }
}

/// Where each output `.d` goes.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum OutputMode {
    /// All outputs land inside one chosen directory, keeping their folder names.
    NewDir,
    /// Each output sits next to its input with a suffix on the name.
    Suffix,
}

/// Everything a run needs, cloned onto the worker thread when Run is pressed.
/// Groups the output/preset choices (shown in the main view) with the advanced
/// filter, halo, gate, crop, and ppm knobs (shown in the collapsible panel).
#[derive(Clone)]
pub struct Settings {
    // --- Main view ---
    pub preset: PresetChoice,
    pub output_mode: OutputMode,
    pub output_dir: String,
    pub suffix: String,
    pub overwrite: bool,
    pub write_report: bool,

    // --- Advanced: vertical filter ---
    pub mz_half_width: u32,
    pub min_feature_length: usize,
    pub max_internal_gap: usize,
    pub iterations: usize,
    pub min_window_intensity: u64,
    pub min_feature_intensity: u64,

    // --- Advanced: ppm m/z window (overrides mz_half_width when on) ---
    pub use_ppm: bool,
    pub mz_ppm: f64,

    // --- Advanced: horizontal-halo filter ---
    pub halo: bool,
    pub halo_peak_fraction: f64,
    pub halo_mz_idx_half_width: u32,
    pub halo_scan_half_width: usize,

    // --- Advanced: gate override (else the preset decides) ---
    pub follow_preset_gates: bool,
    pub ms1_polygon: bool,
    pub dia_ms1_window: bool,
    pub dia_window: bool,

    // --- Advanced: region-of-interest crop (blank field = no bound) ---
    pub mz_min: String,
    pub mz_max: String,
    pub im_min: String,
    pub im_max: String,
    pub rt_min: String,
    pub rt_max: String,
    pub min_intensity: String,
    pub max_intensity: String,
    pub crop_only: bool,
}

impl Default for Settings {
    fn default() -> Self {
        let fp = FilterParams::default();
        let hp = HaloParams::default();
        Self {
            preset: PresetChoice::Auto,
            output_mode: OutputMode::Suffix,
            output_dir: String::new(),
            suffix: "_dnoise".to_string(),
            overwrite: false,
            write_report: false,

            mz_half_width: fp.mz_half_width,
            min_feature_length: fp.min_feature_length,
            max_internal_gap: fp.max_internal_gap,
            iterations: fp.num_iterations,
            min_window_intensity: fp.min_window_intensity,
            min_feature_intensity: fp.min_feature_intensity,

            use_ppm: false,
            mz_ppm: 20.0,

            halo: true,
            halo_peak_fraction: hp.peak_fraction,
            halo_mz_idx_half_width: hp.mz_idx_half_width,
            halo_scan_half_width: hp.scan_half_width,

            follow_preset_gates: true,
            ms1_polygon: false,
            dia_ms1_window: false,
            dia_window: false,

            mz_min: String::new(),
            mz_max: String::new(),
            im_min: String::new(),
            im_max: String::new(),
            rt_min: String::new(),
            rt_max: String::new(),
            min_intensity: String::new(),
            max_intensity: String::new(),
            crop_only: false,
        }
    }
}

/// Parse an optional numeric crop bound: blank = `None`; a bad value logs a warning
/// and is treated as `None` (forgiving, so one typo doesn't abort a batch).
fn parse_opt<T: std::str::FromStr>(s: &str, name: &str, log: &mut dyn FnMut(String)) -> Option<T> {
    let t = s.trim();
    if t.is_empty() {
        return None;
    }
    match t.parse::<T>() {
        Ok(v) => Some(v),
        Err(_) => {
            log(format!("ignoring invalid {name} = '{t}'"));
            None
        }
    }
}

impl Settings {
    /// The vertical-filter parameters, applying the ppm override when enabled
    /// (which needs the run calibration, hence `input`).
    pub fn filter_params(&self, input: &Path, log: &mut dyn FnMut(String)) -> FilterParams {
        let mut p = FilterParams {
            mz_half_width: self.mz_half_width,
            min_feature_length: self.min_feature_length,
            max_internal_gap: self.max_internal_gap,
            min_window_intensity: self.min_window_intensity,
            min_feature_intensity: self.min_feature_intensity,
            num_iterations: self.iterations,
        };
        if self.use_ppm {
            match dnoise::tof_half_width_for_ppm(input, self.mz_ppm, None) {
                Ok(hw) => p.mz_half_width = hw,
                Err(e) => log(format!(
                    "ppm conversion failed ({e}); using mz_half_width = {}",
                    p.mz_half_width
                )),
            }
        }
        p
    }

    /// The halo parameters, or `None` when the halo filter is disabled.
    pub fn halo_params(&self) -> Option<HaloParams> {
        self.halo.then_some(HaloParams {
            peak_fraction: self.halo_peak_fraction,
            mz_idx_half_width: self.halo_mz_idx_half_width,
            scan_half_width: self.halo_scan_half_width,
        })
    }

    /// Which MS1 gates to enable: from the preset (default) or the explicit
    /// override toggles.
    pub fn gates(&self, input: &Path) -> Gates {
        if self.follow_preset_gates {
            resolve_gates(self.preset, input)
        } else {
            Gates {
                ms1_polygon: self.ms1_polygon,
                dia_ms1: self.dia_ms1_window,
                dia_window: self.dia_window,
            }
        }
    }

    /// The parsed region-of-interest crop (blank fields become unset bounds).
    pub fn crop_params(&self, log: &mut dyn FnMut(String)) -> CropParams {
        CropParams {
            mz_min: parse_opt(&self.mz_min, "mz_min", log),
            mz_max: parse_opt(&self.mz_max, "mz_max", log),
            im_min: parse_opt(&self.im_min, "im_min", log),
            im_max: parse_opt(&self.im_max, "im_max", log),
            rt_min: parse_opt(&self.rt_min, "rt_min", log),
            rt_max: parse_opt(&self.rt_max, "rt_max", log),
            min_intensity: parse_opt(&self.min_intensity, "min_intensity", log),
            max_intensity: parse_opt(&self.max_intensity, "max_intensity", log),
        }
    }

    /// Compute the output `.d` path for one input, or a user-facing error string.
    pub fn output_path(&self, input: &Path) -> Result<PathBuf, String> {
        let name = input
            .file_name()
            .ok_or_else(|| "input has no folder name".to_string())?;
        match self.output_mode {
            OutputMode::NewDir => {
                let dir = self.output_dir.trim();
                if dir.is_empty() {
                    return Err("choose an output folder".to_string());
                }
                Ok(Path::new(dir).join(name))
            }
            OutputMode::Suffix => {
                let stem = input
                    .file_stem()
                    .ok_or_else(|| "input has no name".to_string())?
                    .to_string_lossy();
                let sfx = if self.suffix.trim().is_empty() {
                    "_dnoise"
                } else {
                    self.suffix.trim()
                };
                let parent = input.parent().unwrap_or_else(|| Path::new("."));
                Ok(parent.join(format!("{stem}{sfx}.d")))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn settings(mode: OutputMode) -> Settings {
        Settings {
            output_mode: mode,
            ..Settings::default()
        }
    }

    #[test]
    fn suffix_output_sits_next_to_input() {
        let s = settings(OutputMode::Suffix);
        let out = s.output_path(Path::new("/data/sample.d")).unwrap();
        assert_eq!(out, PathBuf::from("/data/sample_dnoise.d"));
    }

    #[test]
    fn empty_suffix_falls_back_to_default() {
        let mut s = settings(OutputMode::Suffix);
        s.suffix = "   ".to_string();
        let out = s.output_path(Path::new("/data/sample.d")).unwrap();
        assert_eq!(out, PathBuf::from("/data/sample_dnoise.d"));
    }

    #[test]
    fn newdir_output_keeps_folder_name() {
        let mut s = settings(OutputMode::NewDir);
        s.output_dir = "/out".to_string();
        let out = s.output_path(Path::new("/data/sample.d")).unwrap();
        assert_eq!(out, PathBuf::from("/out/sample.d"));
    }

    #[test]
    fn newdir_without_folder_errors() {
        let s = settings(OutputMode::NewDir);
        assert!(s.output_path(Path::new("/data/sample.d")).is_err());
    }

    #[test]
    fn filter_params_carry_advanced_knobs() {
        let mut s = Settings::default();
        s.min_feature_length = 9;
        s.iterations = 4;
        s.mz_half_width = 7;
        // use_ppm is false, so no calibration read — a dummy path is fine.
        let mut warnings = Vec::new();
        let p = s.filter_params(Path::new("/x.d"), &mut |w| warnings.push(w));
        assert_eq!(p.min_feature_length, 9);
        assert_eq!(p.num_iterations, 4);
        assert_eq!(p.mz_half_width, 7);
        assert!(warnings.is_empty());
    }

    #[test]
    fn crop_params_parse_blank_valid_and_invalid() {
        let mut s = Settings::default();
        s.mz_min = "400.5".to_string();
        s.mz_max = String::new(); // blank -> None
        s.min_intensity = "50".to_string();
        s.rt_min = "not-a-number".to_string(); // invalid -> None + warning
        let mut warnings = Vec::new();
        let c = s.crop_params(&mut |w| warnings.push(w));
        assert_eq!(c.mz_min, Some(400.5));
        assert_eq!(c.mz_max, None);
        assert_eq!(c.min_intensity, Some(50));
        assert_eq!(c.rt_min, None);
        assert_eq!(warnings.len(), 1);
        assert!(warnings[0].contains("rt_min"));
    }

    #[test]
    fn gate_override_bypasses_preset() {
        let mut s = Settings::default();
        s.preset = PresetChoice::None;
        s.follow_preset_gates = false;
        s.ms1_polygon = true;
        s.dia_window = true;
        let g = s.gates(Path::new("/x.d"));
        assert!(g.ms1_polygon && g.dia_window && !g.dia_ms1);
    }

    #[test]
    fn explicit_presets_map_to_gates() {
        // Non-Auto presets ignore the path, so a dummy is fine.
        let p = Path::new("/nonexistent.d");
        let dda = resolve_gates(PresetChoice::Dda, p);
        assert!(dda.ms1_polygon && !dda.dia_ms1 && !dda.dia_window);
        let dia = resolve_gates(PresetChoice::Dia, p);
        assert!(!dia.ms1_polygon && dia.dia_ms1 && dia.dia_window);
        let none = resolve_gates(PresetChoice::None, p);
        assert!(!none.ms1_polygon && !none.dia_ms1 && !none.dia_window);
    }
}
