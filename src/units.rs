//! Physical equivalents of the raw-unit parameters, for the run log and the
//! provenance record.
//!
//! The streak, halo, smoothing and centroiding parameters are set in raw
//! instrument units (TOF indices and mobility scans) and stay that way. What
//! they mean physically depends on the run's calibration: a TOF half-width is a
//! ppm window that narrows with m/z (timsTOF TOF is a sqrt law), and a scan
//! count is a 1/K0 span that depends on the mobility ramp. [`for_run`] converts
//! each active raw-unit parameter once per run so the log and provenance say
//! what the numbers meant for that run.

use crate::mobility::{self, ScanToMobility};
use crate::params::{FilterParams, Stages};
use crate::tsr::{ConvertableDomain, MetadataReader, Tof2MzConverter};
use crate::{DnoiseError, Result};
use std::path::Path;

/// Reference m/z values at which TOF-index parameters are reported in ppm.
pub const REFERENCE_MZ: [f64; 3] = [400.0, 800.0, 1200.0];

/// Axis a raw-unit parameter is measured on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum RawAxis {
    /// TOF index (m/z axis).
    TofIndex,
    /// Mobility scan.
    Scans,
}

/// One raw-unit parameter with its physical equivalent for a run.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct UnitEquivalent {
    /// Config key of the parameter (e.g. `mz_half_width`).
    pub param: &'static str,
    /// Its value in raw units.
    pub value: u64,
    /// The raw unit.
    pub axis: RawAxis,
    /// TOF parameters: `[m/z, ppm]` pairs, the ppm span of `value` TOF indices
    /// at each of [`REFERENCE_MZ`]. Empty for scan parameters.
    pub ppm_at_mz: Vec<[f64; 2]>,
    /// Scan parameters: the 1/K0 span of `value` scans at [`Self::at_inv_k0`].
    pub inv_k0: Option<f64>,
    /// Scan parameters: the 1/K0 at which the span was taken (the mid scan).
    pub at_inv_k0: Option<f64>,
}

impl UnitEquivalent {
    /// One human-readable log line.
    pub fn line(&self) -> String {
        match self.axis {
            RawAxis::TofIndex => {
                let ppm = self
                    .ppm_at_mz
                    .iter()
                    .map(|[mz, ppm]| format!("{ppm:.1} ppm at m/z {mz:.0}"))
                    .collect::<Vec<_>>()
                    .join(", ");
                format!("{} = {} TOF index = {ppm}", self.param, self.value)
            }
            RawAxis::Scans => format!(
                "{} = {} scans = {:.4} 1/K0 (at 1/K0 {:.3})",
                self.param,
                self.value,
                self.inv_k0.unwrap_or(f64::NAN),
                self.at_inv_k0.unwrap_or(f64::NAN)
            ),
        }
    }
}

/// The run's m/z and mobility conversions, as the gates use them.
pub struct RunScales {
    /// TOF index <-> m/z.
    pub mz: Tof2MzConverter,
    /// Scan <-> 1/K0 on the run's configured scale.
    pub im: ScanToMobility,
    /// Mobility scans per frame (`MAX(NumScans)`).
    pub num_scans: u32,
}

impl RunScales {
    /// ppm span of `n` TOF indices at `mz`.
    fn ppm(&self, n: u64, mz: f64) -> f64 {
        let t = self.mz.invert(mz);
        (self.mz.convert(t + n as f64) - mz) / mz * 1e6
    }

    /// 1/K0 at the mid scan and the 1/K0 span of `n` scans there.
    fn inv_k0(&self, n: u64) -> (f64, f64) {
        let mid = f64::from(self.num_scans) / 2.0;
        let per_scan = (self.im.convert(mid + 0.5) - self.im.convert(mid - 0.5)).abs();
        (self.im.convert(mid), per_scan * n as f64)
    }

    fn equivalent(&self, param: &'static str, value: u64, axis: RawAxis) -> UnitEquivalent {
        match axis {
            RawAxis::TofIndex => UnitEquivalent {
                param,
                value,
                axis,
                ppm_at_mz: REFERENCE_MZ.map(|mz| [mz, self.ppm(value, mz)]).to_vec(),
                inv_k0: None,
                at_inv_k0: None,
            },
            RawAxis::Scans => {
                let (at, span) = self.inv_k0(value);
                UnitEquivalent {
                    param,
                    value,
                    axis,
                    ppm_at_mz: Vec::new(),
                    inv_k0: Some(span),
                    at_inv_k0: Some(at),
                }
            }
        }
    }
}

/// Physical equivalents of every active raw-unit parameter, in pipeline order.
/// Intensity thresholds, iteration counts and physical-unit parameters (the MS1
/// gate pads, the crop) are not listed.
pub fn equivalents(params: &FilterParams, stages: &Stages, s: &RunScales) -> Vec<UnitEquivalent> {
    use RawAxis::{Scans, TofIndex};
    let as_u64 = |v: usize| v as u64;
    let mut out = vec![
        s.equivalent("mz_half_width", params.mz_half_width.into(), TofIndex),
        s.equivalent(
            "min_feature_length",
            as_u64(params.min_feature_length),
            Scans,
        ),
        s.equivalent("max_internal_gap", as_u64(params.max_internal_gap), Scans),
    ];
    if let Some(h) = stages.halo {
        out.push(s.equivalent(
            "halo_mz_idx_half_width",
            h.mz_idx_half_width.into(),
            TofIndex,
        ));
        out.push(s.equivalent("halo_scan_half_width", as_u64(h.scan_half_width), Scans));
    }
    if let Some(m) = stages.denoise_msms {
        out.push(s.equivalent("msms_mz_half_width", m.mz_half_width.into(), TofIndex));
        out.push(s.equivalent(
            "msms_min_feature_length",
            as_u64(m.min_feature_length),
            Scans,
        ));
        out.push(s.equivalent("msms_max_internal_gap", as_u64(m.max_internal_gap), Scans));
    }
    if let Some(w) = stages.dia_window {
        out.push(s.equivalent("dia_window_scan_pad", w.scan_pad.into(), Scans));
    }
    if let Some(w) = stages.dda_window {
        out.push(s.equivalent("dda_window_scan_pad", w.scan_pad.into(), Scans));
    }
    if let Some(p) = stages.smooth {
        out.push(s.equivalent(
            "smooth_mz_idx_half_width",
            p.mz_idx_half_width.into(),
            TofIndex,
        ));
        out.push(s.equivalent("smooth_scan_half_width", as_u64(p.scan_half_width), Scans));
    }
    if let Some(w) = stages.watershed {
        out.push(s.equivalent("watershed_box_scan", w.box_scan.into(), Scans));
        out.push(s.equivalent("watershed_box_mz_idx", w.box_mz_idx.into(), TofIndex));
        out.push(s.equivalent(
            "watershed_max_tof_offset",
            w.max_tof_offset.into(),
            TofIndex,
        ));
    }
    if let Some(b) = stages.box_centroid {
        out.push(s.equivalent(
            "box_centroid_mz_idx_half",
            b.mz_idx_half_width.into(),
            TofIndex,
        ));
        out.push(s.equivalent("box_centroid_scan_half", b.scan_half_width.into(), Scans));
    }
    out
}

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tsr::Scan2ImConverter;

    fn scales() -> RunScales {
        RunScales {
            // A typical timsTOF range: m/z 95..1705 over 400k TOF indices.
            mz: Tof2MzConverter::from_boundaries(95.0, 1705.0, 400_000),
            // 1/K0 1.45 at scan 0 down to 0.65 at scan 1000: 0.0008 per scan.
            im: ScanToMobility::Linear(Scan2ImConverter::from_boundaries(0.65, 1.45, 1000)),
            num_scans: 1000,
        }
    }

    #[test]
    fn tof_ppm_narrows_with_mz_as_the_sqrt_law() {
        let e = scales().equivalent("mz_half_width", 3, RawAxis::TofIndex);
        let ppm: Vec<f64> = e.ppm_at_mz.iter().map(|p| p[1]).collect();
        assert!(ppm[0] > ppm[1] && ppm[1] > ppm[2], "{ppm:?}");
        // dm/m = 2 dt slope / sqrt(m): ppm scales as 1/sqrt(m/z).
        assert!((ppm[0] / ppm[1] - 2f64.sqrt()).abs() < 1e-3, "{ppm:?}");
        assert!(e.line().starts_with("mz_half_width = 3 TOF index = "));
    }

    #[test]
    fn scan_span_uses_the_mobility_ramp() {
        let e = scales().equivalent("max_internal_gap", 2, RawAxis::Scans);
        assert!((e.inv_k0.unwrap() - 0.0016).abs() < 1e-9);
        assert!((e.at_inv_k0.unwrap() - 1.05).abs() < 1e-9);
        assert_eq!(
            e.line(),
            "max_internal_gap = 2 scans = 0.0016 1/K0 (at 1/K0 1.050)"
        );
    }

    #[test]
    fn only_active_stages_are_listed() {
        let p = FilterParams::default();
        let stages = Stages::default();
        let names: Vec<_> = equivalents(&p, &stages, &scales())
            .iter()
            .map(|e| e.param)
            .collect();
        assert_eq!(
            names,
            ["mz_half_width", "min_feature_length", "max_internal_gap"]
        );
        let halo = crate::params::HaloParams::default();
        let stages = Stages {
            halo: Some(&halo),
            ..Stages::default()
        };
        assert_eq!(equivalents(&p, &stages, &scales()).len(), 5);
    }
}
