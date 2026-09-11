//! Processing history stored beside the vendor files, without custom TDF schema.
use crate::{DenoiseStats, DnoiseError, FilterParams, Result, RunOptions, Stages};
use serde_json::{Value, json};
use std::io::{Read, Write};
use std::path::Path;

/// The completion/history sidecar inside each processed `.d` folder.
pub const FILE_NAME: &str = "dnoise.provenance.json";
/// Executable settings for the latest operation (when the config feature is enabled).
pub const CONFIG_NAME: &str = "dnoise.config.toml";

/// Actual temporal evidence for nonempty central events in processed frames.
/// A neighbor is counted only if it supplies points inside the matching event;
/// the central observation is excluded. Counts include repeated uses across events.
#[derive(Debug, Clone, Copy, Default, serde::Serialize)]
pub struct NeighborUsage {
    /// Nonempty central events evaluated with temporal support enabled.
    pub events: u64,
    /// Evaluated events receiving no points from other observations.
    pub events_without_neighbors: u64,
    /// Total contributing neighbor observations across evaluated events.
    pub neighbors_used: u64,
    /// Largest contributing-neighbor count for one event.
    pub max_neighbors_used: u64,
}

impl NeighborUsage {
    pub(crate) fn add(&mut self, other: Self) {
        self.events += other.events;
        self.events_without_neighbors += other.events_without_neighbors;
        self.neighbors_used += other.neighbors_used;
        self.max_neighbors_used = self.max_neighbors_used.max(other.max_neighbors_used);
    }
}

/// Which requested geometry gates found geometry in this acquisition.
#[derive(Debug, Clone, Copy, Default, serde::Serialize)]
pub struct ActiveGates {
    /// Experimental temporal support enabled for MS1 filtering.
    pub ms1_neighbors: bool,
    /// Experimental temporal support enabled for PRM target filtering.
    pub prm_neighbors: bool,
    /// Experimental temporal support enabled for DIA window filtering.
    pub dia_neighbors: bool,
    /// Experimental independent prm-PASEF isolation-event filtering.
    pub prm_per_event: bool,
    /// DDA MS1 selection polygon.
    pub ms1_polygon: bool,
    /// DIA MS1 selection windows.
    pub dia_ms1: bool,
    /// DIA MS/MS scan windows.
    pub dia_window: bool,
    /// DDA MS/MS scan windows.
    pub dda_window: bool,
    /// Continuous scan-varying DIA regions were used for fragment filtering.
    pub dia_scan_varying: bool,
    /// Independent DIA MS/MS window filtering.
    pub dia_per_window: bool,
}

/// Read and validate existing dnoise history. Absence does not prove a file is raw:
/// old versions did not write markers, and users can remove ancillary files.
pub fn read(input: &Path) -> Result<Option<Value>> {
    let file = match std::fs::File::open(input.join(FILE_NAME)) {
        Ok(file) => file,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(e) => return Err(e.into()),
    };
    let mut bytes = Vec::new();
    file.take(16 * 1024 * 1024 + 1).read_to_end(&mut bytes)?;
    if bytes.len() > 16 * 1024 * 1024 {
        return Err(DnoiseError::InvalidInput(
            "provenance exceeds 16 MiB".into(),
        ));
    }
    let value: Value = serde_json::from_slice(&bytes)
        .map_err(|e| DnoiseError::InvalidInput(format!("invalid dnoise provenance: {e}")))?;
    if value["schema_version"] != 1
        || value["tool"] != "dnoise"
        || value["history"].as_array().is_none_or(|h| h.is_empty())
    {
        return Err(DnoiseError::InvalidInput(
            "unsupported or incomplete dnoise provenance".into(),
        ));
    }
    Ok(Some(value))
}

/// Shared report payload for CLI, GUI, and library users, including dry runs.
pub fn report(
    input: &Path,
    params: &FilterParams,
    stages: &Stages,
    options: &RunOptions,
    stats: &DenoiseStats,
) -> Value {
    let mut warnings = Vec::new();
    let acquisition = crate::detect_acquisition(input);
    if matches!(
        acquisition,
        Ok(crate::Acquisition::PrmPasef | crate::Acquisition::Mixed)
    ) {
        warnings.push("PRM/mixed acquisition: discovery acquisition gates disabled; explicit crops can change MS/MS data");
    }
    if stats.active_gates.prm_per_event {
        warnings.push("experimental prm-PASEF MS/MS denoising; target boundaries enforced; downstream quantitative fidelity is unvalidated");
    }
    if stats.active_gates.ms1_neighbors
        || stats.active_gates.prm_neighbors
        || stats.active_gates.dia_neighbors
    {
        warnings.push("experimental neighbor support: summed evidence changes keep decisions; native output intensities preserved unless smoothing/centroiding enabled; quantitative fidelity is unvalidated");
    }
    if stats.active_gates.dia_scan_varying {
        warnings.push("experimental scanning DIA filtering: continuous quadrupole steps grouped for mobility support; native scan-dependent geometry preserved; downstream identification and quantitative fidelity are unvalidated");
    }
    if stats.ms1_frames == 0
        && !options.crop_only
        && stages.denoise_msms.is_none()
        && !stages.filter_all_frames
    {
        warnings.push("no MS1 frames; MS1-only denoising has no points to filter");
    }
    if stats.multiple_calibrations {
        warnings.push("multiple calibration references detected; physical gates/crops needing these conversions are rejected");
    }
    if stages.ms1_polygon.is_some() && !stats.active_gates.ms1_polygon && !options.crop_only {
        warnings.push("MS1 polygon requested but no applicable polygon was active");
    }
    if stages.dia_ms1.is_some() && !stats.active_gates.dia_ms1 && !options.crop_only {
        warnings.push("DIA MS1 gate requested but no applicable windows were active");
    }
    let mut value = json!({
        "schema_version": 1, "tool": "dnoise", "software_version": env!("CARGO_PKG_VERSION"),
        "build_revision": env!("DNOISE_BUILD_REVISION"),
        "status": if stats.dry_run { "estimate" } else { "completed" },
        "operation": if options.crop_only { "crop_only" } else { "denoise" },
        "completed_unix_seconds": std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0),
        "input_name": input.file_name().map(|s|s.to_string_lossy()),
        "acquisition": acquisition.map(|a|format!("{a:?}")).unwrap_or_else(|_|"unknown".into()),
        "frame_batch_size": options.frame_batch_size.unwrap_or(2048),
        "validation": if options.skip_validation { "structural" } else { "full" },
        "config": {"filter": params, "stages": stages, "crop": options.crop, "crop_only": options.crop_only},
        "sample": options.sample.map(|s| json!({"fraction":s.fraction,"seed":s.seed})),
        "stats": stats, "warnings": warnings,
        "intensity_basis": "stored_integer",
    });
    #[cfg(feature = "config")]
    {
        value["effective_config"] = serde_json::to_value(crate::config::Config::from_effective(
            params, stages, options,
        ))
        .expect("validated config is serializable");
    }
    // Keep the binding mutable in library-only builds as well.
    #[cfg(feature = "config")]
    {
        value["effective_config"]["threads"] = json!(stats.worker_threads);
    }
    value["dry_run"] = json!(stats.dry_run);
    value["elapsed_seconds"] = json!(stats.elapsed_seconds);
    value["stats"]["kept_pct"] = json!(if stats.raw_points > 0 {
        100.0 * stats.kept_points as f64 / stats.raw_points as f64
    } else {
        0.0
    });
    value["stats"]["ms1_kept_pct"] = json!(if stats.raw_ms1_points > 0 {
        100.0 * stats.kept_ms1_points as f64 / stats.raw_ms1_points as f64
    } else {
        0.0
    });
    let percent = |kept: u64, raw: u64| {
        if raw == 0 {
            None
        } else {
            Some(100.0 * kept as f64 / raw as f64)
        }
    };
    value["stats"]["ms1_intensity_retained_pct"] = json!(percent(
        stats.kept_ms1_summed_intensity,
        stats.raw_ms1_summed_intensity
    ));
    value["stats"]["msms_intensity_retained_pct"] = json!(percent(
        stats.kept_msms_summed_intensity,
        stats.raw_msms_summed_intensity
    ));
    value
}

pub(crate) fn write(
    input: &Path,
    output: &Path,
    params: &FilterParams,
    stages: &Stages,
    options: &RunOptions,
    stats: &DenoiseStats,
) -> Result<()> {
    let mut record =
        read(input)?.unwrap_or_else(|| json!({"schema_version":1,"tool":"dnoise","history":[]}));
    record["history"]
        .as_array_mut()
        .unwrap()
        .push(report(input, params, stages, options, stats));
    let bytes =
        serde_json::to_vec_pretty(&record).map_err(|e| DnoiseError::InvalidInput(e.to_string()))?;
    let mut file = std::fs::File::create(output.join(FILE_NAME))?;
    file.write_all(&bytes)?;
    file.write_all(b"\n")?;
    file.sync_all()?;
    #[cfg(feature = "config")]
    crate::config::Config::from_effective(params, stages, options)
        .save(&output.join(CONFIG_NAME))
        .map_err(DnoiseError::InvalidInput)?;
    // A library-only run cannot export a new TOML recipe. Do not leave the
    // input's copied recipe looking like settings for this new operation.
    #[cfg(not(feature = "config"))]
    match std::fs::remove_file(output.join(CONFIG_NAME)) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => return Err(e.into()),
    }
    Ok(())
}
