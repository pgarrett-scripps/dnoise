//! Run settings shared between the UI and the worker thread, plus the small
//! preset -> gate and output-path helpers.

use dnoise::Acquisition;
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
#[derive(Clone)]
pub struct Settings {
    pub preset: PresetChoice,
    pub output_mode: OutputMode,
    pub output_dir: String,
    pub suffix: String,
    pub overwrite: bool,
    pub write_report: bool,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            preset: PresetChoice::Auto,
            output_mode: OutputMode::Suffix,
            output_dir: String::new(),
            suffix: "_dnoise".to_string(),
            overwrite: false,
            write_report: false,
        }
    }
}

impl Settings {
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
