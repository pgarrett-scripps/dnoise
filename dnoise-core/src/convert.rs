//! Calibration converters: TOF index <-> m/z and scan <-> 1/K0 on the linear
//! (acquisition-range) scale, with fractional inverses.
//!
//! These are timsrust 0.4.2's straight-line formulas, kept here as plain structs
//! so the gates can be built without any reader. Build them from the run's
//! `GlobalMetadata` (`MzAcqRangeLower/Upper`, `DigitizerNumSamples`,
//! `OneOverK0AcqRangeLower/Upper`, `MAX(NumScans)`); the `dnoise` crate's
//! `tsr::MetadataReader` does exactly that. For Bruker's calibrated mobility
//! scale see [`crate::mobility`].

/// Scalar conversion with a fractional inverse (timsrust 0.4's trait).
pub trait ConvertableDomain {
    /// Forward conversion (index -> physical value).
    fn convert<T: Into<f64> + Copy>(&self, value: T) -> f64;
    /// Inverse conversion (physical value -> fractional index).
    fn invert<T: Into<f64> + Copy>(&self, value: T) -> f64;
}

/// TOF index -> m/z: `mz = (sqrt(mz_min) + slope * tof)^2` (timsrust 0.4.2).
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Tof2MzConverter {
    tof_intercept: f64,
    tof_slope: f64,
}

impl Tof2MzConverter {
    /// Line between `mz_min` at TOF 0 and `mz_max` at `tof_max_index`.
    pub fn from_boundaries(mz_min: f64, mz_max: f64, tof_max_index: u32) -> Self {
        let tof_intercept = mz_min.sqrt();
        let tof_slope = (mz_max.sqrt() - tof_intercept) / tof_max_index as f64;
        Self {
            tof_intercept,
            tof_slope,
        }
    }
}

impl ConvertableDomain for Tof2MzConverter {
    fn convert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        let tof_index: f64 = value.into();
        (self.tof_intercept + self.tof_slope * tof_index).powi(2)
    }
    fn invert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        let mz_value: f64 = value.into();
        (mz_value.sqrt() - self.tof_intercept) / self.tof_slope
    }
}

/// Scan index -> 1/K0: straight line from `im_max` at scan 0 to `im_min` at
/// `max(NumScans)` (timsrust 0.4.2).
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Scan2ImConverter {
    scan_intercept: f64,
    scan_slope: f64,
}

impl Scan2ImConverter {
    /// Line between `im_max` at scan 0 and `im_min` at `scan_max_index`.
    pub fn from_boundaries(im_min: f64, im_max: f64, scan_max_index: u32) -> Self {
        let scan_intercept = im_max;
        let scan_slope = (im_min - scan_intercept) / scan_max_index as f64;
        Self {
            scan_intercept,
            scan_slope,
        }
    }
}

impl ConvertableDomain for Scan2ImConverter {
    fn convert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        let scan_index: f64 = value.into();
        self.scan_intercept + self.scan_slope * scan_index
    }
    fn invert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        let im_value: f64 = value.into();
        (im_value - self.scan_intercept) / self.scan_slope
    }
}
