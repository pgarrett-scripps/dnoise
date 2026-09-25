//! TIMS scan ↔ 1/K0 conversion on Bruker's acquisition-calibrated scale.
//!
//! timsrust 0.4 converts a scan index to inverse reduced mobility with a
//! straight line between `OneOverK0AcqRangeUpper` (scan 0) and
//! `OneOverK0AcqRangeLower` (the largest `NumScans`). Bruker's software (the
//! timsdata SDK, timsControl, DataAnalysis) uses the run's `TimsCalibration`
//! row instead, and the two differ by up to ~0.03 1/K0 on a 0.64-1.45
//! acquisition range. The selection polygon and the diaPASEF window pads are
//! defined on Bruker's scale, so the MS1 gates convert with this module by
//! default ([`MobilityScale::Calibrated`]).
//!
//! # The ModelType 2 calibration
//!
//! Ported from koth's `tims_calibration` module, which was verified against
//! `libtimsdata.so` on 304 runs (maximum error 1.4e-14 1/K0). For a scan `s`:
//!
//! ```text
//! V(s)      = C2 + (C2 - C3) / C1 * (C4 + C0 - s)     ramp voltage at scan s
//! g(V)      = V / (C6 * V + C7)                        1/K0 = 1 / (C6 + C7 / V)
//! 1/K0(s)   = g(V)                                     for C8 <= V <= C9
//!           = g(C8) + g'(C8) * (V - C8)                for V < C8
//!           = g(C9) + g'(C9) * (V - C9)                for V > C9
//! g'(V)     = C7 / (C6 * V + C7)^2
//! ```
//!
//! `V(s)` is linear and `g` is increasing, so the inverse is closed-form on
//! each of the three pieces. Any other `ModelType` is refused.

use crate::convert::{ConvertableDomain, Scan2ImConverter};
use serde::{Deserialize, Serialize};

/// Which 1/K0 scale the gates and the mobility crop use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MobilityScale {
    /// Bruker acquisition calibration (`TimsCalibration`), the scale the
    /// selection polygon and the diaPASEF windows are defined on.
    #[default]
    Calibrated,
    /// timsrust's straight line between the acquisition-range bounds: dnoise's
    /// only behaviour before 0.4.0.
    Linear,
}

/// The only calibration model type implemented (and the only one observed).
pub const SUPPORTED_MODEL_TYPE: i64 = 2;

/// One `TimsCalibration` row of `ModelType` 2.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TimsCalibrationModel {
    c: [f64; 10],
}

impl TimsCalibrationModel {
    /// Build from `ModelType` and `C0..C9`. Errors on an unsupported model type
    /// or on coefficients that would divide by zero or not be monotonic.
    pub fn new(model_type: i64, c: [f64; 10]) -> std::result::Result<Self, String> {
        if model_type != SUPPORTED_MODEL_TYPE {
            return Err(format!(
                "TimsCalibration ModelType {model_type} is not supported (only {SUPPORTED_MODEL_TYPE}); \
                 set mobility_scale = \"linear\" (--linear-mobility) to use the uncalibrated scale"
            ));
        }
        let m = Self { c };
        let ok = c.iter().all(|v| v.is_finite())
            && c[1] != 0.0
            && c[2] != c[3]
            && c[7] > 0.0
            && c[8] < c[9]
            && c[6] * c[8] + c[7] > 0.0
            && c[6] * c[9] + c[7] > 0.0;
        if !ok {
            return Err(format!("invalid TimsCalibration coefficients {c:?}"));
        }
        Ok(m)
    }

    #[inline]
    fn g(&self, v: f64) -> f64 {
        v / (self.c[6] * v + self.c[7])
    }

    #[inline]
    fn dg(&self, v: f64) -> f64 {
        let d = self.c[6] * v + self.c[7];
        self.c[7] / (d * d)
    }

    #[inline]
    fn voltage(&self, scan: f64) -> f64 {
        let c = &self.c;
        c[2] + (c[2] - c[3]) / c[1] * (c[4] + c[0] - scan)
    }

    #[inline]
    fn scan_at_voltage(&self, v: f64) -> f64 {
        let c = &self.c;
        c[4] + c[0] - (v - c[2]) * c[1] / (c[2] - c[3])
    }

    /// 1/K0 (V·s/cm²) at scan `scan` (fractional allowed), as the timsdata
    /// SDK's `tims_scannum_to_oneoverk0` computes it.
    #[inline]
    pub fn scan_to_inv_k0(&self, scan: f64) -> f64 {
        let v = self.voltage(scan);
        let (lo, hi) = (self.c[8], self.c[9]);
        if v < lo {
            self.g(lo) + self.dg(lo) * (v - lo)
        } else if v > hi {
            self.g(hi) + self.dg(hi) * (v - hi)
        } else {
            self.g(v)
        }
    }

    /// Fractional scan at which the calibration gives `inv_k0`: the exact
    /// inverse of [`Self::scan_to_inv_k0`].
    #[inline]
    pub fn inv_k0_to_scan(&self, inv_k0: f64) -> f64 {
        let (lo, hi) = (self.c[8], self.c[9]);
        let (k_lo, k_hi) = (self.g(lo), self.g(hi));
        let v = if inv_k0 < k_lo {
            lo + (inv_k0 - k_lo) / self.dg(lo)
        } else if inv_k0 > k_hi {
            hi + (inv_k0 - k_hi) / self.dg(hi)
        } else {
            inv_k0 * self.c[7] / (1.0 - inv_k0 * self.c[6])
        };
        self.scan_at_voltage(v)
    }
}

/// Run-level scan ↔ 1/K0 conversion on the chosen [`MobilityScale`], with the
/// same `convert`/`invert` shape as timsrust's converters.
#[derive(Debug, Clone, Copy)]
pub enum ScanToMobility {
    /// timsrust's linear scale.
    Linear(Scan2ImConverter),
    /// The run's `TimsCalibration` row.
    Calibrated(TimsCalibrationModel),
}

impl ScanToMobility {
    /// 1/K0 at (fractional) scan `scan`.
    #[inline]
    pub fn convert(&self, scan: f64) -> f64 {
        match self {
            Self::Linear(l) => l.convert(scan),
            Self::Calibrated(m) => m.scan_to_inv_k0(scan),
        }
    }

    /// Fractional scan at 1/K0 `inv_k0`.
    #[inline]
    pub fn invert(&self, inv_k0: f64) -> f64 {
        match self {
            Self::Linear(l) => l.invert(inv_k0),
            Self::Calibrated(m) => m.inv_k0_to_scan(inv_k0),
        }
    }

    /// Which scale this converter uses.
    pub fn scale(&self) -> MobilityScale {
        match self {
            Self::Linear(_) => MobilityScale::Linear,
            Self::Calibrated(_) => MobilityScale::Calibrated,
        }
    }
}
