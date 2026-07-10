//! Region-of-interest crop (a.k.a. trim): keep only points inside an axis-aligned
//! box in `(m/z, 1/K0, intensity)` space, plus a retention-time window handled at
//! the frame level by the writer.
//!
//! Unlike the denoising filters (which decide *signal vs noise*), the crop is a
//! blunt subset of the acquisition: it is how a user carves a smaller `.d` out of
//! a large one (a mass/mobility/RT region of interest, or an intensity floor) for
//! sharing, faster downstream searches, or building test fixtures. It composes
//! with the denoiser (applied as an extra AND on the keep mask) or can run on its
//! own via the CLI's `--crop-only`.
//!
//! Like [`crate::polygon`] / [`crate::dia_ms1`] the gate is calibration-free at
//! run time: the writer converts the physical `(m/z, 1/K0)` bounds into integer
//! `(TOF index, scan)` ranges once (using the run calibration), and each point is
//! then tested with a handful of integer comparisons. Retention time is compared
//! per frame against the raw `Frames.Time`, so RT-cropped frames are emitted empty
//! rather than deleted — the frame axis (and every table that references it) stays
//! structurally valid and Bruker-SDK-compatible.

use crate::params::CropParams;

/// Integer `(TOF index, scan, intensity)` bounds for the point-level crop, built
/// once from [`CropParams`] plus the run calibration. All bounds are inclusive; an
/// unset axis widens to the full representable range so it never rejects a point.
#[derive(Debug, Clone, Copy)]
pub struct CropGate {
    tof_lo: u32,
    tof_hi: u32,
    scan_lo: u32,
    scan_hi: u32,
    int_lo: u32,
    int_hi: u32,
    /// True when at least one point-level bound is narrower than the full range;
    /// lets the writer skip the per-point test entirely for an RT-only crop.
    active: bool,
}

impl CropGate {
    /// Build the point-level crop from physical bounds.
    ///
    /// * `mz_to_tof(mz)` — fractional TOF index of an m/z (monotonic increasing).
    /// * `k0_to_scan(k0)` — fractional mobility-scan index of a `1/K0` (monotonic
    ///   decreasing: larger `1/K0` maps to a smaller scan).
    /// * `num_scans` — clamps the upper scan bound to a real scan index.
    ///
    /// The retention-time bounds in `p` are not consumed here; the writer applies
    /// them per frame.
    pub fn build(
        p: &CropParams,
        num_scans: usize,
        mz_to_tof: impl Fn(f64) -> f64,
        k0_to_scan: impl Fn(f64) -> f64,
    ) -> Self {
        let clamp_u32 = |x: f64| x.max(0.0).min(u32::MAX as f64);
        // m/z increases with TOF index, so mz_min -> low index, mz_max -> high index.
        let tof_lo = p
            .mz_min
            .map(|mz| clamp_u32(mz_to_tof(mz).floor()) as u32)
            .unwrap_or(0);
        let tof_hi = p
            .mz_max
            .map(|mz| clamp_u32(mz_to_tof(mz).ceil()) as u32)
            .unwrap_or(u32::MAX);
        // 1/K0 decreases with scan, so the *max* mobility is the *low* scan bound.
        let scan_cap = num_scans.saturating_sub(1) as u32;
        let scan_lo = p
            .im_max
            .map(|k0| clamp_u32(k0_to_scan(k0).floor()) as u32)
            .unwrap_or(0);
        let scan_hi = p
            .im_min
            .map(|k0| (clamp_u32(k0_to_scan(k0).ceil()) as u32).min(scan_cap))
            .unwrap_or(u32::MAX);
        let int_lo = p.min_intensity.unwrap_or(0);
        let int_hi = p.max_intensity.unwrap_or(u32::MAX);

        let active = tof_lo > 0
            || tof_hi < u32::MAX
            || scan_lo > 0
            || scan_hi < u32::MAX
            || int_lo > 0
            || int_hi < u32::MAX;

        Self {
            tof_lo,
            tof_hi,
            scan_lo,
            scan_hi,
            int_lo,
            int_hi,
            active,
        }
    }

    /// True when the gate constrains at least one point-level axis (m/z, mobility,
    /// or intensity). An RT-only crop returns `false`, so the writer skips the
    /// per-point pass.
    pub fn is_active(&self) -> bool {
        self.active
    }

    /// Keep a point iff it lies inside every bound (all inclusive).
    #[inline]
    pub fn keep(&self, scan: u32, tof: u32, intensity: u32) -> bool {
        scan >= self.scan_lo
            && scan <= self.scan_hi
            && tof >= self.tof_lo
            && tof <= self.tof_hi
            && intensity >= self.int_lo
            && intensity <= self.int_hi
    }

    /// AND the crop into an existing keep mask over `frame`'s points.
    pub fn apply(&self, scan: &[u32], tof: &[u32], intensity: &[u32], keep: &mut [bool]) {
        if !self.active {
            return;
        }
        for i in 0..keep.len() {
            if keep[i] && !self.keep(scan[i], tof[i], intensity[i]) {
                keep[i] = false;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Identity converters: TOF index == m/z, scan == 1/K0 (still monotone
    /// decreasing is not required for these algebra checks). Lets the tests reason
    /// in gate-native integer units.
    fn gate(p: &CropParams) -> CropGate {
        CropGate::build(p, 1000, |mz| mz, |k0| k0)
    }

    #[test]
    fn unset_crop_is_inactive_and_keeps_everything() {
        let g = gate(&CropParams::default());
        assert!(!g.is_active());
        assert!(g.keep(0, 0, 0));
        assert!(g.keep(999, u32::MAX, u32::MAX));
    }

    #[test]
    fn mz_and_intensity_bounds_are_inclusive() {
        let p = CropParams {
            mz_min: Some(100.0),
            mz_max: Some(200.0),
            min_intensity: Some(10),
            max_intensity: Some(50),
            ..CropParams::default()
        };
        let g = gate(&p);
        assert!(g.is_active());
        assert!(g.keep(0, 100, 10)); // on both edges
        assert!(g.keep(0, 200, 50));
        assert!(!g.keep(0, 99, 30)); // below m/z
        assert!(!g.keep(0, 201, 30)); // above m/z
        assert!(!g.keep(0, 150, 9)); // below intensity
        assert!(!g.keep(0, 150, 51)); // above intensity
    }

    #[test]
    fn mobility_max_maps_to_low_scan() {
        // Real timsTOF calibration maps 1/K0 to scan monotonically *decreasing*, so
        // the higher mobility bound becomes the lower scan bound. Model that with
        // scan = 1000 - k0: im_max=300 -> scan_lo=700, im_min=100 -> scan_hi=900.
        let p = CropParams {
            im_min: Some(100.0),
            im_max: Some(300.0),
            ..CropParams::default()
        };
        let g = CropGate::build(&p, 1000, |mz| mz, |k0| 1000.0 - k0);
        assert!(g.keep(700, 0, 1));
        assert!(g.keep(900, 0, 1));
        assert!(!g.keep(699, 0, 1));
        assert!(!g.keep(901, 0, 1));
    }

    #[test]
    fn apply_ands_into_existing_mask() {
        let p = CropParams {
            mz_min: Some(10.0),
            ..CropParams::default()
        };
        let g = gate(&p);
        let scan = [0u32, 0, 0];
        let tof = [5u32, 20, 30];
        let inten = [1u32, 1, 1];
        let mut keep = [true, true, false];
        g.apply(&scan, &tof, &inten, &mut keep);
        assert_eq!(keep, [false, true, false]); // tof 5 cropped, 20 kept, 30 already off
    }
}
