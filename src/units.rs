//! Physical equivalents of the raw-unit parameters, for the run log and the
//! provenance record. The conversion lives in `dnoise_core::units` (re-exported
//! here); [`for_run`] reads the run's calibration from `analysis.tdf`.

pub use dnoise_core::units::{REFERENCE_MZ, RawAxis, RunScales, UnitEquivalent, equivalents};

use crate::mobility::{self, ScanToMobility};
use crate::params::{FilterParams, Stages};
use crate::tsr::MetadataReader;
use crate::{DnoiseError, Result};
use std::path::Path;

/// Read the run's calibration from `analysis.tdf` and compute [`equivalents`].
/// Uses the configured mobility scale, falling back to the linear scale when the
/// calibrated one cannot be read (e.g. several calibration rows).
pub fn for_run(
    tdf_path: &Path,
    params: &FilterParams,
    stages: &Stages,
) -> Result<Vec<UnitEquivalent>> {
    let md = MetadataReader::new(tdf_path).map_err(|e| DnoiseError::Metadata(e.to_string()))?;
    let conn = rusqlite::Connection::open_with_flags(
        tdf_path,
        rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
    )?;
    let num_scans: Option<i64> =
        conn.query_row("SELECT MAX(NumScans) FROM Frames", [], |r| r.get(0))?;
    let num_scans = num_scans.unwrap_or(0).max(0) as u32;
    if num_scans < 2 {
        return Err(DnoiseError::Metadata("no mobility scans".into()));
    }
    let im = mobility::load(tdf_path, stages.mobility_scale, md.im_converter)
        .unwrap_or(ScanToMobility::Linear(md.im_converter));
    let scales = RunScales {
        mz: md.mz_converter,
        im,
        num_scans,
    };
    Ok(equivalents(params, stages, &scales))
}
