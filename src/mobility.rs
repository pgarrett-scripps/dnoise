//! TIMS scan ↔ 1/K0 conversion: reading the run's `TimsCalibration` row.
//!
//! The conversion itself ([`MobilityScale`], [`TimsCalibrationModel`],
//! [`ScanToMobility`]) lives in `dnoise_core::mobility` and is re-exported here;
//! this module adds the `analysis.tdf` reader.

use std::path::Path;

pub use dnoise_core::mobility::{
    MobilityScale, SUPPORTED_MODEL_TYPE, ScanToMobility, TimsCalibrationModel,
};

use crate::error::{DnoiseError, Result};
use crate::tsr::Scan2ImConverter;

/// Read the conversion for `analysis.tdf`. `linear` is timsrust's converter for
/// the run, returned as is for [`MobilityScale::Linear`]. The calibrated scale
/// needs exactly one `TimsCalibration` row referenced by the frames; callers
/// already refuse multi-calibration runs before converting.
pub(crate) fn load(
    tdf_path: &Path,
    scale: MobilityScale,
    linear: Scan2ImConverter,
) -> Result<ScanToMobility> {
    if scale == MobilityScale::Linear {
        return Ok(ScanToMobility::Linear(linear));
    }
    let err = |e: String| {
        DnoiseError::Metadata(format!(
            "{}: cannot read the TimsCalibration mobility calibration ({e}); set \
             mobility_scale = \"linear\" (--linear-mobility) to use the uncalibrated scale",
            tdf_path.display()
        ))
    };
    let conn =
        rusqlite::Connection::open_with_flags(tdf_path, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
            .map_err(|e| err(e.to_string()))?;
    let ids: Vec<i64> = conn
        .prepare("SELECT DISTINCT TimsCalibration FROM Frames")
        .and_then(|mut s| s.query_map([], |r| r.get(0))?.collect())
        .map_err(|e| err(e.to_string()))?;
    let [id] = ids[..] else {
        return Err(err(format!(
            "frames reference {} calibration rows, expected one",
            ids.len()
        )));
    };
    let (model_type, c) = conn
        .query_row(
            "SELECT ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9 FROM TimsCalibration WHERE Id = ?1",
            [id],
            |r| {
                let mut c = [0.0; 10];
                for (i, v) in c.iter_mut().enumerate() {
                    *v = r.get(i + 1)?;
                }
                Ok((r.get::<_, i64>(0)?, c))
            },
        )
        .map_err(|e| err(e.to_string()))?;
    TimsCalibrationModel::new(model_type, c)
        .map(ScanToMobility::Calibrated)
        .map_err(err)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tsr::ConvertableDomain;
    use serde::Deserialize;

    #[derive(Deserialize)]
    struct Fixture {
        cases: Vec<Case>,
    }

    #[derive(Deserialize)]
    struct Case {
        run: String,
        calibration_id: i64,
        model_type: i64,
        coefficients: Vec<f64>,
        num_scans: u32,
        one_over_k0_acq_range: [f64; 2],
        scans: Vec<f64>,
        sdk_one_over_k0: Vec<f64>,
    }

    // Bruker SDK reference values, copied from koth (no raw data stored).
    fn fixture() -> Fixture {
        serde_json::from_str(include_str!("../tests/data/tims_calibration_sdk.json")).unwrap()
    }

    fn model(case: &Case) -> TimsCalibrationModel {
        TimsCalibrationModel::new(
            case.model_type,
            case.coefficients.clone().try_into().unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn matches_sdk_on_real_runs() {
        let fx = fixture();
        assert!(fx.cases.len() >= 5);
        for case in &fx.cases {
            let m = model(case);
            for (&s, &want) in case.scans.iter().zip(&case.sdk_one_over_k0) {
                let got = m.scan_to_inv_k0(s);
                assert!(
                    (got - want).abs() <= 1e-12,
                    "{} cal {} scan {s}: got {got} want {want}",
                    case.run,
                    case.calibration_id
                );
            }
        }
    }

    #[test]
    fn inverse_round_trips_on_all_three_pieces() {
        for case in &fixture().cases {
            let m = model(case);
            for &s in &case.scans {
                let back = m.inv_k0_to_scan(m.scan_to_inv_k0(s));
                assert!(
                    (back - s).abs() <= 1e-6 * s.abs().max(1.0),
                    "{} scan {s} -> {back}",
                    case.run
                );
            }
        }
    }

    #[test]
    fn linear_scale_is_timsrust_and_differs_from_calibrated() {
        let fx = fixture();
        let case = &fx.cases[0];
        let [lo, hi] = case.one_over_k0_acq_range;
        let lin = Scan2ImConverter::from_boundaries(lo, hi, case.num_scans);
        let l = ScanToMobility::Linear(lin);
        let c = ScanToMobility::Calibrated(model(case));
        for s in [0.0, 100.0, 468.25, 935.0] {
            assert_eq!(l.convert(s), lin.convert(s));
            assert_eq!(l.invert(l.convert(s)), lin.invert(lin.convert(s)));
        }
        // On the benchmark calibration the scales disagree by ~0.032 1/K0 at scan 0.
        let d = c.convert(0.0) - l.convert(0.0);
        assert!((d - 0.0318).abs() < 1e-3, "{d}");
        assert_eq!(c.scale(), MobilityScale::Calibrated);
        assert_eq!(l.scale(), MobilityScale::Linear);
    }

    #[test]
    fn refuses_unknown_model_type_and_bad_coefficients() {
        let c = [
            1.0, 935.0, 239.3, 103.0, 33.6, 1.0, -0.027, 171.4, 16.8, 1732.7,
        ];
        assert!(TimsCalibrationModel::new(1, c).is_err());
        assert!(TimsCalibrationModel::new(2, c).is_ok());
        let mut flat = c;
        flat[3] = flat[2];
        assert!(TimsCalibrationModel::new(2, flat).is_err());
    }

    #[test]
    fn load_reads_the_referenced_row_and_refuses_a_missing_table() {
        let dir = std::env::temp_dir().join(format!("dnoise_mobility_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("analysis.tdf");
        let _ = std::fs::remove_file(&path);
        let conn = rusqlite::Connection::open(&path).unwrap();
        conn.execute_batch(
            "CREATE TABLE Frames (Id INTEGER, TimsCalibration INTEGER);
            INSERT INTO Frames VALUES (1, 7), (2, 7);",
        )
        .unwrap();
        let lin = Scan2ImConverter::from_boundaries(0.64, 1.45, 936);
        assert!(load(&path, MobilityScale::Calibrated, lin).is_err());
        assert_eq!(
            load(&path, MobilityScale::Linear, lin).unwrap().scale(),
            MobilityScale::Linear
        );
        conn.execute_batch("CREATE TABLE TimsCalibration (Id INTEGER, ModelType INTEGER,
              C0 REAL, C1 REAL, C2 REAL, C3 REAL, C4 REAL, C5 REAL, C6 REAL, C7 REAL, C8 REAL, C9 REAL);
            INSERT INTO TimsCalibration VALUES (7, 2, 1, 935, 239.34640606518187, 102.96946384662338,
              33.64485981308411, 1, -0.026580764926972034, 171.42849749723894, 16.838457909616054,
              1732.6649859338625);").unwrap();
        let conv = load(&path, MobilityScale::Calibrated, lin).unwrap();
        assert!((conv.convert(0.0) - 1.481819026431424).abs() < 1e-12);
        // Windows cannot delete a file an open connection holds.
        drop(conn);
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
