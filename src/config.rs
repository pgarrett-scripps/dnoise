//! On-disk run configuration (`dnoise.toml`), shared by the CLI and the GUI.
//!
//! Every field is an `Option`; a missing key falls back to the caller's default
//! (the CLI resolves CLI flag > config value > built-in default). This is a plain
//! data record of the TOML schema — it does **not** apply any resolution or build
//! [`crate::FilterParams`] / [`crate::Stages`] itself; the front ends do that so
//! they can layer their own precedence and UI state on top.
//!
//! Available only with the `config` feature (which pulls in `serde` + `toml`); the
//! `cli` feature enables it, and the GUI opts in explicitly.

use serde::{Deserialize, Serialize};
use std::path::Path;

/// The full `dnoise.toml` schema. Unknown keys are rejected on load to catch
/// typos; `None` fields are omitted on save (TOML has no null), so a saved file
/// contains only the knobs that were actually set.
///
/// Each field maps one-to-one to the identically named CLI flag / config key; see
/// the README options table and `dnoise.toml` for the per-key meaning. Field docs
/// are omitted here to avoid restating that reference for all ~50 knobs.
#[allow(missing_docs)]
#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Config {
    // Vertical filter.
    pub mz_half_width: Option<u32>,
    pub min_feature_length: Option<usize>,
    pub max_internal_gap: Option<usize>,
    pub min_window_intensity: Option<u64>,
    pub min_feature_intensity: Option<u64>,
    pub iterations: Option<usize>,
    pub frame_half_width: Option<usize>,
    // Horizontal-halo filter.
    pub halo: Option<bool>,
    pub halo_peak_fraction: Option<f64>,
    pub halo_mz_idx_half_width: Option<u32>,
    pub halo_scan_half_width: Option<usize>,
    // ddaPASEF MS/MS denoising.
    pub denoise_msms: Option<bool>,
    pub msms_mz_half_width: Option<u32>,
    pub msms_min_feature_length: Option<usize>,
    pub msms_max_internal_gap: Option<usize>,
    pub msms_min_window_intensity: Option<u64>,
    pub msms_min_feature_intensity: Option<u64>,
    pub msms_iterations: Option<usize>,
    // Intensity smoothing.
    pub smooth: Option<bool>,
    pub smooth_mz_idx_half_width: Option<u32>,
    pub smooth_scan_half_width: Option<usize>,
    pub smooth_iterations: Option<usize>,
    // Watershed centroider.
    pub watershed: Option<bool>,
    pub watershed_box_scan: Option<u32>,
    pub watershed_box_mz_idx: Option<u32>,
    pub watershed_min_seed_intensity: Option<u64>,
    pub watershed_min_centroid_total: Option<u64>,
    pub watershed_max_tof_offset: Option<u32>,
    // Greedy box centroider.
    pub box_centroid: Option<bool>,
    pub box_centroid_mz_idx_half: Option<u32>,
    pub box_centroid_scan_half: Option<u32>,
    pub box_centroid_min_total: Option<u64>,
    // diaPASEF isolation-window features.
    pub dia_window: Option<bool>,
    pub dia_window_scan_pad: Option<u32>,
    pub dia_per_window: Option<bool>,
    // ddaPASEF MS/MS out-of-window gate.
    pub dda_window: Option<bool>,
    pub dda_window_scan_pad: Option<u32>,
    pub dia_ms1_window: Option<bool>,
    pub dia_ms1_mz_pad: Option<f64>,
    pub dia_ms1_im_pad: Option<f64>,
    // MS1 selection-polygon gate.
    pub ms1_polygon: Option<bool>,
    pub ms1_polygon_mz_pad: Option<f64>,
    pub ms1_polygon_im_pad: Option<f64>,
    // Region-of-interest crop.
    pub mz_min: Option<f64>,
    pub mz_max: Option<f64>,
    pub im_min: Option<f64>,
    pub im_max: Option<f64>,
    pub rt_min: Option<f64>,
    pub rt_max: Option<f64>,
    pub min_intensity: Option<u32>,
    pub max_intensity: Option<u32>,
    pub crop_only: Option<bool>,
    // ppm-based m/z window.
    pub mz_ppm: Option<f64>,
    pub mz_ppm_ref: Option<f64>,
    // Operational.
    pub all_frames: Option<bool>,
    pub threads: Option<usize>,
}

impl Config {
    /// Load and parse a TOML config file. Returns a human-readable error string on
    /// a read or parse failure (including unknown keys).
    pub fn load(path: &Path) -> Result<Self, String> {
        let text = std::fs::read_to_string(path)
            .map_err(|e| format!("reading config {}: {e}", path.display()))?;
        toml::from_str(&text).map_err(|e| format!("parsing config {}: {e}", path.display()))
    }

    /// Serialize to a TOML string (only the set keys appear).
    pub fn to_toml_string(&self) -> Result<String, String> {
        toml::to_string_pretty(self).map_err(|e| format!("serializing config: {e}"))
    }

    /// Serialize and write to `path`.
    pub fn save(&self, path: &Path) -> Result<(), String> {
        let text = self.to_toml_string()?;
        std::fs::write(path, text).map_err(|e| format!("writing config {}: {e}", path.display()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_through_toml() {
        let c = Config {
            mz_half_width: Some(4),
            min_feature_length: Some(6),
            ms1_polygon: Some(true),
            mz_min: Some(400.5),
            ..Config::default()
        };
        let s = c.to_toml_string().unwrap();
        // None fields are omitted.
        assert!(!s.contains("mz_max"));
        assert!(s.contains("mz_half_width = 4"));
        let back: Config = toml::from_str(&s).unwrap();
        assert_eq!(c, back);
    }

    #[test]
    fn unknown_key_is_rejected() {
        assert!(toml::from_str::<Config>("bogus_key = 3").is_err());
    }

    #[test]
    fn save_then_load_round_trips_through_a_file() {
        let path =
            std::env::temp_dir().join(format!("dnoise_cfg_roundtrip_{}.toml", std::process::id()));
        let _ = std::fs::remove_file(&path);
        let c = Config {
            mz_half_width: Some(9),
            halo: Some(true),
            im_max: Some(1.3),
            ..Config::default()
        };
        c.save(&path).unwrap();
        let back = Config::load(&path).unwrap();
        assert_eq!(c, back);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn load_missing_file_is_err() {
        let path =
            std::env::temp_dir().join(format!("dnoise_cfg_absent_{}.toml", std::process::id()));
        let _ = std::fs::remove_file(&path);
        assert!(Config::load(&path).is_err());
    }
}
