//! On-disk run configuration (`dnoise.toml`), shared by the CLI and the GUI.
//!
//! Missing keys use the shared resolver's defaults. Front ends merge their explicit
//! overrides into Config, then call `resolve` for a consistent executable recipe.
//! Available with the `config` feature; CLI and GUI both enable it.

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
    pub ms1_neighbor_radius: Option<usize>,
    pub prm_neighbor_radius: Option<usize>,
    pub dia_neighbor_radius: Option<usize>,
    pub neighbor_max_rt_gap: Option<f64>,
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
    pub dia_ms1_overlap: Option<bool>,
    /// Deprecated in 0.5.0 (reach is now unlimited): still accepted so 0.4.0
    /// configs load, but ignored with a warning and never written.
    #[serde(skip_serializing)]
    pub dia_ms1_overlap_reach: Option<f64>,
    // MS1 selection-polygon gate.
    pub ms1_polygon: Option<bool>,
    pub ms1_polygon_mz_pad: Option<f64>,
    pub ms1_polygon_im_pad: Option<f64>,
    pub ms1_polygon_overlap: Option<bool>,
    /// Deprecated in 0.5.0, like `dia_ms1_overlap_reach`.
    #[serde(skip_serializing)]
    pub ms1_polygon_overlap_reach: Option<f64>,
    /// `"calibrated"` (default) or `"linear"`.
    pub mobility_scale: Option<crate::mobility::MobilityScale>,
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
    pub frame_batch_size: Option<usize>,
    pub skip_validation: Option<bool>,
}

impl Config {
    /// Deprecated keys that are set: accepted for compatibility with older
    /// configs, ignored by [`Config::resolve`] (which warns about each).
    pub fn deprecated_keys(&self) -> Vec<&'static str> {
        let mut keys = Vec::new();
        if self.dia_ms1_overlap_reach.is_some() {
            keys.push("dia_ms1_overlap_reach");
        }
        if self.ms1_polygon_overlap_reach.is_some() {
            keys.push("ms1_polygon_overlap_reach");
        }
        keys
    }

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
    fn deprecated_overlap_reach_loads_but_is_ignored_and_not_saved() {
        // 0.4.0 configs set these; they must still load.
        let c: Config =
            toml::from_str("ms1_polygon_overlap_reach = 0.1\ndia_ms1_overlap_reach = 0.0").unwrap();
        assert_eq!(
            c.deprecated_keys(),
            vec!["dia_ms1_overlap_reach", "ms1_polygon_overlap_reach"]
        );
        let s = c.to_toml_string().unwrap();
        assert!(!s.contains("overlap_reach"));
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

use crate::*;
/// Owned, fully resolved settings shared by the CLI and desktop worker.
#[derive(Debug, Clone)]
pub struct ResolvedConfig {
    /// Resolved processing parameters.
    pub filter: FilterParams,
    /// Resolved processing parameters.
    pub halo: Option<HaloParams>,
    /// Resolved processing parameters.
    pub denoise_msms: Option<MsmsFilterParams>,
    /// Resolved processing parameters.
    pub smooth: Option<SmoothParams>,
    /// Resolved processing parameters.
    pub watershed: Option<WatershedParams>,
    /// Resolved processing parameters.
    pub box_centroid: Option<BoxCentroidParams>,
    /// Resolved processing parameters.
    pub dia_window: Option<DiaWindowParams>,
    /// Resolved processing parameters.
    pub dda_window: Option<DdaWindowParams>,
    /// Resolved processing parameters.
    pub dia_ms1: Option<DiaMs1WindowParams>,
    /// Resolved processing parameters.
    pub ms1_polygon: Option<Ms1PolygonParams>,
    /// Resolved scan ↔ 1/K0 scale.
    pub mobility_scale: crate::mobility::MobilityScale,
    /// Resolved processing option.
    pub crop: CropParams,
    /// Resolved processing option.
    pub crop_only: bool,
    /// Resolved processing option.
    pub all_frames: bool,
    /// Resolved processing option.
    pub frame_half_width: usize,
    /// Resolved temporal support settings.
    pub neighbors: NeighborParams,
    /// Resolved processing option.
    pub dia_per_window: bool,
}
impl ResolvedConfig {
    /// Borrow the optional stages for the library pipeline.
    pub fn stages(&self) -> Stages<'_> {
        Stages {
            filter_all_frames: self.all_frames,
            frame_half_width: self.frame_half_width,
            neighbors: self.neighbors,
            dia_per_window: self.dia_per_window,
            halo: self.halo.as_ref(),
            denoise_msms: self.denoise_msms.as_ref(),
            smooth: self.smooth.as_ref(),
            watershed: self.watershed.as_ref(),
            box_centroid: self.box_centroid.as_ref(),
            dia_window: self.dia_window.as_ref(),
            dda_window: self.dda_window.as_ref(),
            dia_ms1: self.dia_ms1.as_ref(),
            ms1_polygon: self.ms1_polygon.as_ref(),
            mobility_scale: self.mobility_scale,
        }
    }
}
impl Config {
    /// Resolve defaults and ppm conversion using the same rules in both front ends.
    pub fn resolve(&self, input: &Path) -> crate::Result<ResolvedConfig> {
        for key in self.deprecated_keys() {
            tracing::warn!(
                "{key} / --{} is deprecated and ignored: since 0.5.0 a feature kept by an \
                 MS1 gate is kept whole (unlimited reach)",
                key.replace('_', "-")
            );
        }
        if let (Some(old), Some(new)) = (self.frame_half_width, self.ms1_neighbor_radius)
            && old != new
        {
            return Err(DnoiseError::InvalidInput(
                "frame_half_width and ms1_neighbor_radius disagree; use only ms1_neighbor_radius"
                    .into(),
            ));
        }
        let msms = self.denoise_msms.unwrap_or(false) || self.all_frames.unwrap_or(false);
        let d = FilterParams::default();
        let filter = FilterParams {
            mz_half_width: self.mz_half_width.unwrap_or(d.mz_half_width),
            min_feature_length: self.min_feature_length.unwrap_or(d.min_feature_length),
            max_internal_gap: self.max_internal_gap.unwrap_or(d.max_internal_gap),
            min_window_intensity: self.min_window_intensity.unwrap_or(d.min_window_intensity),
            min_feature_intensity: self
                .min_feature_intensity
                .unwrap_or(d.min_feature_intensity),
            num_iterations: self.iterations.unwrap_or(d.num_iterations),
        };
        let d = HaloParams::default();
        let halo = HaloParams {
            peak_fraction: self.halo_peak_fraction.unwrap_or(d.peak_fraction),
            mz_idx_half_width: self.halo_mz_idx_half_width.unwrap_or(d.mz_idx_half_width),
            scan_half_width: self.halo_scan_half_width.unwrap_or(d.scan_half_width),
        };
        let halo = self.halo.unwrap_or(true).then_some(halo);
        let d = MsmsFilterParams::default();
        let denoise_msms = MsmsFilterParams {
            mz_half_width: self.msms_mz_half_width.unwrap_or(d.mz_half_width),
            min_feature_length: self.msms_min_feature_length.unwrap_or(d.min_feature_length),
            max_internal_gap: self.msms_max_internal_gap.unwrap_or(d.max_internal_gap),
            min_window_intensity: self
                .msms_min_window_intensity
                .unwrap_or(d.min_window_intensity),
            min_feature_intensity: self
                .msms_min_feature_intensity
                .unwrap_or(d.min_feature_intensity),
            num_iterations: self.msms_iterations.unwrap_or(d.num_iterations),
        };
        let denoise_msms = self.denoise_msms.unwrap_or(false).then_some(denoise_msms);
        let d = SmoothParams::default();
        let smooth = SmoothParams {
            mz_idx_half_width: self.smooth_mz_idx_half_width.unwrap_or(d.mz_idx_half_width),
            scan_half_width: self.smooth_scan_half_width.unwrap_or(d.scan_half_width),
            iterations: self.smooth_iterations.unwrap_or(d.iterations),
        };
        let smooth = self.smooth.unwrap_or(false).then_some(smooth);
        let d = WatershedParams::default();
        let watershed = WatershedParams {
            box_scan: self.watershed_box_scan.unwrap_or(d.box_scan),
            box_mz_idx: self.watershed_box_mz_idx.unwrap_or(d.box_mz_idx),
            min_seed_intensity: self
                .watershed_min_seed_intensity
                .unwrap_or(d.min_seed_intensity),
            min_centroid_total: self
                .watershed_min_centroid_total
                .unwrap_or(d.min_centroid_total),
            max_tof_offset: self.watershed_max_tof_offset.unwrap_or(d.max_tof_offset),
        };
        let watershed = self.watershed.unwrap_or(false).then_some(watershed);
        let d = BoxCentroidParams::default();
        let box_centroid = BoxCentroidParams {
            mz_idx_half_width: self.box_centroid_mz_idx_half.unwrap_or(d.mz_idx_half_width),
            scan_half_width: self.box_centroid_scan_half.unwrap_or(d.scan_half_width),
            min_centroid_total: self.box_centroid_min_total.unwrap_or(d.min_centroid_total),
        };
        let box_centroid = self.box_centroid.unwrap_or(false).then_some(box_centroid);
        let d = DiaWindowParams::default();
        let dia_window = DiaWindowParams {
            scan_pad: self.dia_window_scan_pad.unwrap_or(d.scan_pad),
        };
        let dia_window = self.dia_window.unwrap_or(msms).then_some(dia_window);
        let d = DdaWindowParams::default();
        let dda_window = DdaWindowParams {
            scan_pad: self.dda_window_scan_pad.unwrap_or(d.scan_pad),
        };
        let dda_window = self.dda_window.unwrap_or(msms).then_some(dda_window);
        let d = DiaMs1WindowParams::default();
        let dia_ms1 = DiaMs1WindowParams {
            mz_pad: self.dia_ms1_mz_pad.unwrap_or(d.mz_pad),
            im_pad: self.dia_ms1_im_pad.unwrap_or(d.im_pad),
            overlap: self.dia_ms1_overlap.unwrap_or(d.overlap),
        };
        let dia_ms1 = self.dia_ms1_window.unwrap_or(true).then_some(dia_ms1);
        let d = Ms1PolygonParams::default();
        let ms1_polygon = Ms1PolygonParams {
            mz_pad: self.ms1_polygon_mz_pad.unwrap_or(d.mz_pad),
            im_pad: self.ms1_polygon_im_pad.unwrap_or(d.im_pad),
            overlap: self.ms1_polygon_overlap.unwrap_or(d.overlap),
        };
        let ms1_polygon = self.ms1_polygon.unwrap_or(true).then_some(ms1_polygon);
        let mut filter = filter;
        if let Some(ppm) = self.mz_ppm {
            if !ppm.is_finite()
                || ppm <= 0.0
                || self.mz_ppm_ref.is_some_and(|m| !m.is_finite() || m <= 0.0)
            {
                return Err(DnoiseError::InvalidInput(
                    "ppm and reference m/z must be finite and positive".into(),
                ));
            }
            filter.mz_half_width = tof_half_width_for_ppm(input, ppm, self.mz_ppm_ref)?;
        }
        let resolved = ResolvedConfig {
            filter,
            halo,
            denoise_msms,
            smooth,
            watershed,
            box_centroid,
            dia_window,
            dda_window,
            dia_ms1,
            ms1_polygon,
            mobility_scale: self.mobility_scale.unwrap_or_default(),
            crop: CropParams {
                mz_min: self.mz_min,
                mz_max: self.mz_max,
                im_min: self.im_min,
                im_max: self.im_max,
                rt_min: self.rt_min,
                rt_max: self.rt_max,
                min_intensity: self.min_intensity,
                max_intensity: self.max_intensity,
            },
            crop_only: self.crop_only.unwrap_or(false),
            all_frames: self.all_frames.unwrap_or(false),
            frame_half_width: self
                .ms1_neighbor_radius
                .or(self.frame_half_width)
                .unwrap_or(0),
            neighbors: NeighborParams {
                prm_radius: self.prm_neighbor_radius.unwrap_or(0),
                dia_radius: self.dia_neighbor_radius.unwrap_or(0),
                max_rt_gap_seconds: self.neighbor_max_rt_gap.unwrap_or(5.0),
            },
            dia_per_window: self.dia_per_window.unwrap_or(msms),
        };
        crate::validation::parameters(
            &resolved.filter,
            &resolved.stages(),
            &RunOptions {
                crop: Some(&resolved.crop),
                crop_only: resolved.crop_only,
                frame_batch_size: self.frame_batch_size,
                ..Default::default()
            },
        )?;
        Ok(resolved)
    }
    /// Capture an executable TOML recipe from the parameters actually used.
    pub fn from_effective(filter: &FilterParams, stages: &Stages, options: &RunOptions) -> Self {
        let mut c = Self {
            frame_batch_size: options.frame_batch_size,
            skip_validation: Some(options.skip_validation),
            threads: Some(rayon::current_num_threads()),
            all_frames: Some(stages.filter_all_frames),
            ms1_neighbor_radius: Some(stages.frame_half_width),
            prm_neighbor_radius: Some(stages.neighbors.prm_radius),
            dia_neighbor_radius: Some(stages.neighbors.dia_radius),
            neighbor_max_rt_gap: Some(stages.neighbors.max_rt_gap_seconds),
            dia_per_window: Some(stages.dia_per_window),
            crop_only: Some(options.crop_only),
            mobility_scale: Some(stages.mobility_scale),
            ..Self::default()
        };
        {
            let p = filter;
            c.mz_half_width = Some(p.mz_half_width);
            c.min_feature_length = Some(p.min_feature_length);
            c.max_internal_gap = Some(p.max_internal_gap);
            c.min_window_intensity = Some(p.min_window_intensity);
            c.min_feature_intensity = Some(p.min_feature_intensity);
            c.iterations = Some(p.num_iterations);
        }
        c.halo = Some(stages.halo.is_some());
        if let Some(p) = stages.halo {
            c.halo_peak_fraction = Some(p.peak_fraction);
            c.halo_mz_idx_half_width = Some(p.mz_idx_half_width);
            c.halo_scan_half_width = Some(p.scan_half_width);
        }
        c.denoise_msms = Some(stages.denoise_msms.is_some());
        if let Some(p) = stages.denoise_msms {
            c.msms_mz_half_width = Some(p.mz_half_width);
            c.msms_min_feature_length = Some(p.min_feature_length);
            c.msms_max_internal_gap = Some(p.max_internal_gap);
            c.msms_min_window_intensity = Some(p.min_window_intensity);
            c.msms_min_feature_intensity = Some(p.min_feature_intensity);
            c.msms_iterations = Some(p.num_iterations);
        }
        c.smooth = Some(stages.smooth.is_some());
        if let Some(p) = stages.smooth {
            c.smooth_mz_idx_half_width = Some(p.mz_idx_half_width);
            c.smooth_scan_half_width = Some(p.scan_half_width);
            c.smooth_iterations = Some(p.iterations);
        }
        c.watershed = Some(stages.watershed.is_some());
        if let Some(p) = stages.watershed {
            c.watershed_box_scan = Some(p.box_scan);
            c.watershed_box_mz_idx = Some(p.box_mz_idx);
            c.watershed_min_seed_intensity = Some(p.min_seed_intensity);
            c.watershed_min_centroid_total = Some(p.min_centroid_total);
            c.watershed_max_tof_offset = Some(p.max_tof_offset);
        }
        c.box_centroid = Some(stages.box_centroid.is_some());
        if let Some(p) = stages.box_centroid {
            c.box_centroid_mz_idx_half = Some(p.mz_idx_half_width);
            c.box_centroid_scan_half = Some(p.scan_half_width);
            c.box_centroid_min_total = Some(p.min_centroid_total);
        }
        c.dia_window = Some(stages.dia_window.is_some());
        if let Some(p) = stages.dia_window {
            c.dia_window_scan_pad = Some(p.scan_pad);
        }
        c.dda_window = Some(stages.dda_window.is_some());
        if let Some(p) = stages.dda_window {
            c.dda_window_scan_pad = Some(p.scan_pad);
        }
        c.dia_ms1_window = Some(stages.dia_ms1.is_some());
        if let Some(p) = stages.dia_ms1 {
            c.dia_ms1_mz_pad = Some(p.mz_pad);
            c.dia_ms1_im_pad = Some(p.im_pad);
            c.dia_ms1_overlap = Some(p.overlap);
        }
        c.ms1_polygon = Some(stages.ms1_polygon.is_some());
        if let Some(p) = stages.ms1_polygon {
            c.ms1_polygon_mz_pad = Some(p.mz_pad);
            c.ms1_polygon_im_pad = Some(p.im_pad);
            c.ms1_polygon_overlap = Some(p.overlap);
        }
        if let Some(crop) = options.crop {
            c.mz_min = crop.mz_min;
            c.mz_max = crop.mz_max;
            c.im_min = crop.im_min;
            c.im_max = crop.im_max;
            c.rt_min = crop.rt_min;
            c.rt_max = crop.rt_max;
            c.min_intensity = crop.min_intensity;
            c.max_intensity = crop.max_intensity;
        }
        c
    }
}

impl Config {
    /// Execute work in the requested Rayon pool without changing process-global state.
    pub fn in_thread_pool<T: Send>(
        &self,
        work: impl FnOnce() -> crate::Result<T> + Send,
    ) -> crate::Result<T> {
        match self.threads {
            Some(threads) => rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| DnoiseError::InvalidInput(format!("cannot create worker pool: {e}")))?
                .install(work),
            None => work(),
        }
    }
}
