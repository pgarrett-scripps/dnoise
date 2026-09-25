//! diaPASEF MS1 out-of-window gate.
//!
//! In diaPASEF the quadrupole tiles the precursor space into isolation windows,
//! each covering an m/z band over a mobility-scan interval. The union of all
//! windows is the set of precursors the method can ever fragment; an MS1 peak
//! outside every window is a precursor that is never isolated, so it can be
//! dropped from the MS1 survey scans.
//!
//! The gate itself is calibration-free: it consumes windows already expressed as
//! **padded integer `(scan, TOF index)` boxes**. [`DiaMs1Gate::from_windows`]
//! does the one-time m/z→TOF and 1/K0→scan conversion, including the
//! physical-unit padding, given the run's converters. It builds a per-scan list of
//! merged TOF intervals so each MS1 point is tested with a single binary search.

use crate::mobility::ScanToMobility;
use crate::params::DiaMs1WindowParams;
use crate::windows::DiaMs1Box;

/// A padded isolation window as an integer `(scan, TOF index)` box, inclusive on
/// both axes. Produced by the writer from a `DiaFrameMsMsWindows` row plus the
/// m/z / 1/K0 padding.
#[derive(Debug, Clone, Copy)]
pub struct TofScanBox {
    /// First mobility scan covered (inclusive).
    pub scan_lo: u32,
    /// Last mobility scan covered (inclusive).
    pub scan_hi: u32,
    /// Lowest TOF index covered (inclusive).
    pub tof_lo: u32,
    /// Highest TOF index covered (inclusive).
    pub tof_hi: u32,
}

/// Per-scan merged TOF intervals for the union of all isolation windows. A point
/// `(scan, tof)` is kept iff `tof` lands in one of `per_scan[scan]`'s intervals.
#[derive(Debug)]
pub struct DiaMs1Gate {
    /// `per_scan[s]` = sorted, non-overlapping `[tof_lo, tof_hi]` intervals
    /// (inclusive) covered at mobility scan `s`. Empty rows keep nothing.
    per_scan: Vec<Vec<(u32, u32)>>,
    /// Gate whole features by any overlap ([`crate::overlap`]) instead of points
    /// (`false` = point by point).
    pub overlap: bool,
}

impl DiaMs1Gate {
    /// Build the gate from padded integer boxes. `num_scans` sizes the per-scan
    /// table; boxes are clamped into `0..num_scans`. Returns `None` when there are
    /// no usable boxes (e.g. ddaPASEF), so callers can skip the gate entirely.
    pub fn build(boxes: &[TofScanBox], num_scans: usize) -> Option<Self> {
        if boxes.is_empty() || num_scans == 0 {
            return None;
        }
        let mut per_scan: Vec<Vec<(u32, u32)>> = vec![Vec::new(); num_scans];
        for b in boxes {
            if b.tof_hi < b.tof_lo {
                continue;
            }
            let lo = b.scan_lo.min(num_scans as u32 - 1);
            let hi = b.scan_hi.min(num_scans as u32 - 1);
            for s in lo..=hi {
                per_scan[s as usize].push((b.tof_lo, b.tof_hi));
            }
        }
        // Sort + merge each scan's intervals so membership is a binary search.
        for row in &mut per_scan {
            row.sort_unstable();
            let mut merged: Vec<(u32, u32)> = Vec::with_capacity(row.len());
            for &(lo, hi) in row.iter() {
                match merged.last_mut() {
                    // Touching or overlapping (allow a 1-index gap to coalesce
                    // adjacent windows): extend the previous interval.
                    Some(last) if lo <= last.1.saturating_add(1) => last.1 = last.1.max(hi),
                    _ => merged.push((lo, hi)),
                }
            }
            *row = merged;
        }
        Some(Self {
            per_scan,
            overlap: false,
        })
    }

    /// Inclusive TOF hull `(lo, hi)` of every interval the gate keeps, over all
    /// scans; `None` when the gate keeps nothing. No kept point lies outside it.
    pub fn tof_span(&self) -> Option<(u32, u32)> {
        let lo = self
            .per_scan
            .iter()
            .filter_map(|r| r.first())
            .map(|i| i.0)
            .min()?;
        let hi = self
            .per_scan
            .iter()
            .filter_map(|r| r.last())
            .map(|i| i.1)
            .max()?;
        Some((lo, hi))
    }

    /// True when the point lies inside some isolation window (and should be kept).
    pub fn contains(&self, scan: u32, tof: u32) -> bool {
        let Some(row) = self.per_scan.get(scan as usize) else {
            return false;
        };
        // First interval whose end is >= tof; the point is inside iff it also
        // clears that interval's start.
        let i = row.partition_point(|&(_, hi)| hi < tof);
        i < row.len() && tof >= row[i].0
    }

    /// Per-point keep mask (in input order): `false` for points outside every
    /// isolation window. `scan` and `tof` are parallel.
    pub fn keep_mask(&self, scan: &[u32], tof: &[u32]) -> Vec<bool> {
        scan.iter()
            .zip(tof)
            .map(|(&s, &t)| self.contains(s, t))
            .collect()
    }
}

impl DiaMs1Gate {
    /// Build the gate from diaPASEF isolation windows (`DiaFrameMsMsWindows`
    /// rows as m/z-by-scan boxes), padded by `p.mz_pad` Th and `p.im_pad` 1/K0
    /// using the run's calibration (`mz_to_tof`, `im`), with `p.overlap` set.
    /// `num_scans` is the largest `Frames.NumScans`. Returns `None` when there
    /// are no windows or no scans.
    pub fn from_windows(
        boxes: &[DiaMs1Box],
        p: &DiaMs1WindowParams,
        mz_to_tof: impl Fn(f64) -> f64,
        im: &ScanToMobility,
        num_scans: usize,
    ) -> Option<Self> {
        let tof_boxes: Vec<TofScanBox> = boxes
            .iter()
            .map(|b| padded_box(b, p, &mz_to_tof, im, num_scans))
            .collect();
        DiaMs1Gate::build(&tof_boxes, num_scans).map(|mut g| {
            g.overlap = p.overlap;
            g
        })
    }
}

/// One diaPASEF isolation window as a padded integer `(scan, TOF)` box: the m/z
/// band widened by `mz_pad` Th on each side (TOF edges rounded outward) and the
/// scan interval widened by `im_pad` 1/K0 ([`window_scans`]). A window is a
/// rectangle in `(m/z, scan)`, so the padded box is exactly the window's
/// Minkowski sum with the pad rectangle: every scan within `im_pad` of the window
/// gets the full padded m/z band.
pub fn padded_box(
    b: &DiaMs1Box,
    p: &DiaMs1WindowParams,
    mz_to_tof: impl Fn(f64) -> f64,
    im: &ScanToMobility,
    num_scans: usize,
) -> TofScanBox {
    // m/z edges -> TOF indices (monotonic), padded by mz_pad on each side.
    let t0 = mz_to_tof(b.mz_lo - p.mz_pad);
    let t1 = mz_to_tof(b.mz_hi + p.mz_pad);
    let tof_lo = t0.min(t1).floor().max(0.0) as u32;
    let tof_hi = t1.max(t0).ceil().max(0.0) as u32;
    let (scan_lo, scan_hi) = window_scans(b.scan_begin, b.scan_end, p.im_pad, im, num_scans);
    TofScanBox {
        scan_lo,
        scan_hi,
        tof_lo,
        tof_hi,
    }
}

/// Inclusive scan range of a `[scan_begin, scan_end)` isolation window, padded by
/// `im_pad` in 1/K0. With no pad the window's own scans are returned exactly;
/// otherwise the edges go scan -> 1/K0 -> padded -> scan, taking min/max so the
/// result is right for either conversion direction. A padded window is rounded
/// outward; the rounding tolerance keeps float noise in the round trip from
/// adding a scan at either edge.
pub fn window_scans(
    scan_begin: u32,
    scan_end: u32,
    im_pad: f64,
    im: &ScanToMobility,
    num_scans: usize,
) -> (u32, u32) {
    const EPS: f64 = 1e-3; // scans: far above round-trip noise, far below one scan
    let last = num_scans.saturating_sub(1) as u32;
    let first_in = scan_begin;
    let last_in = scan_end.saturating_sub(1).max(scan_begin);
    if im_pad <= 0.0 {
        return (first_in.min(last), last_in.min(last));
    }
    let im0 = im.convert(first_in.into());
    let im1 = im.convert(last_in.into());
    let s0 = im.invert(im0.max(im1) + im_pad);
    let s1 = im.invert(im0.min(im1) - im_pad);
    let lo = (s0.min(s1) + EPS).floor().max(0.0) as u32;
    let hi = ((s0.max(s1) - EPS).ceil().max(0.0) as u32).min(last);
    (lo.min(hi), hi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::convert::ConvertableDomain;

    fn benchmark_calibration() -> ScanToMobility {
        // The TimsCalibration row used by the mobility.rs load test.
        ScanToMobility::Calibrated(
            crate::mobility::TimsCalibrationModel::new(
                2,
                [
                    1.0,
                    935.0,
                    239.34640606518187,
                    102.96946384662338,
                    33.64485981308411,
                    1.0,
                    -0.026580764926972034,
                    171.42849749723894,
                    16.838457909616054,
                    1732.6649859338625,
                ],
            )
            .unwrap(),
        )
    }

    #[test]
    fn window_scans_without_pad_is_the_half_open_window_exactly() {
        let im = benchmark_calibration();
        for begin in 0..900u32 {
            assert_eq!(
                window_scans(begin, begin + 30, 0.0, &im, 936),
                (begin, begin + 29)
            );
        }
        assert_eq!(window_scans(920, 1000, 0.0, &im, 936), (920, 935));
    }

    #[test]
    fn window_scans_tiny_pad_adds_no_scan_from_float_noise() {
        // A pad far below one scan's 1/K0 width must not widen the window.
        let im = benchmark_calibration();
        for begin in 1..900u32 {
            assert_eq!(
                window_scans(begin, begin + 30, 1e-9, &im, 936),
                (begin, begin + 29)
            );
        }
    }

    #[test]
    fn window_scans_pad_widens_both_edges() {
        let im = benchmark_calibration();
        let (lo, hi) = window_scans(400, 430, 0.01, &im, 936);
        assert!(lo < 400 && hi > 429, "({lo}, {hi})");
        assert!((im.convert(f64::from(lo)) - im.convert(400.0)).abs() <= 0.01 + 1e-3);
    }

    #[test]
    fn padded_dia_ms1_windows_contain_the_window_and_every_point_within_the_pads() {
        // Containment for the DIA MS1 gate, as tests/polygon_props.rs does for the
        // polygon: every point of the unpadded window is kept, and so is every
        // point whose 1/K0 lies within im_pad and m/z within mz_pad of it.
        let im = benchmark_calibration();
        let mz = crate::convert::Tof2MzConverter::from_boundaries(95.0, 1705.0, 400_000);
        let mz_to_tof = |m: f64| mz.invert(m);
        let num_scans = 936;
        let windows = [
            (0u32, 60u32, 400.0, 425.0),
            (100, 180, 612.5, 637.5),
            (430, 431, 800.0, 801.0),
            (500, 700, 1000.0, 1100.0),
            (880, 936, 1200.0, 1225.0),
        ];
        for (mz_pad, im_pad) in [(0.0, 0.0), (3.0, 0.015), (0.5, 0.003), (5.0, 0.05)] {
            let p = DiaMs1WindowParams {
                mz_pad,
                im_pad,
                overlap: false,
            };
            for &(scan_begin, scan_end, mz_lo, mz_hi) in &windows {
                let b = DiaMs1Box {
                    scan_begin,
                    scan_end,
                    mz_lo,
                    mz_hi,
                };
                let bx = padded_box(&b, &p, mz_to_tof, &im, num_scans);
                let gate = DiaMs1Gate::build(&[bx], num_scans).unwrap();
                let (k_a, k_b) = (
                    im.convert(f64::from(scan_begin)),
                    im.convert(f64::from(scan_end - 1)),
                );
                let (k_lo, k_hi) = (k_a.min(k_b) - im_pad, k_a.max(k_b) + im_pad);
                let tof_lo = mz_to_tof(mz_lo - mz_pad).ceil() as u32;
                let tof_hi = mz_to_tof(mz_hi + mz_pad).floor() as u32;
                for s in 0..num_scans as u32 {
                    let k = im.convert(f64::from(s));
                    let in_window = (scan_begin..scan_end).contains(&s);
                    let in_pad = (k_lo..=k_hi).contains(&k);
                    if in_window || in_pad {
                        for t in [tof_lo, (tof_lo + tof_hi) / 2, tof_hi] {
                            assert!(
                                gate.contains(s, t),
                                "pads ({mz_pad}, {im_pad}), window {scan_begin}..{scan_end}: \
                                 scan {s} TOF {t} dropped"
                            );
                        }
                    }
                    // Over-inclusion is bounded by one scan / one TOF index.
                    if gate.contains(s, (tof_lo + tof_hi) / 2) {
                        let slack = (im.convert(f64::from(s.saturating_sub(1)))
                            - im.convert(f64::from(s + 1)))
                        .abs();
                        assert!(
                            in_window || (k_lo - slack..=k_hi + slack).contains(&k),
                            "scan {s} (1/K0 {k}) kept beyond the pad"
                        );
                    }
                }
                assert!(!gate.contains(bx.scan_lo, bx.tof_lo.saturating_sub(1)) || bx.tof_lo == 0);
                assert!(!gate.contains(bx.scan_lo, bx.tof_hi + 1));
                assert!(
                    tof_lo.saturating_sub(bx.tof_lo) <= 1 && bx.tof_hi.saturating_sub(tof_hi) <= 1
                );
            }
        }
    }

    fn boxed(scan_lo: u32, scan_hi: u32, tof_lo: u32, tof_hi: u32) -> TofScanBox {
        TofScanBox {
            scan_lo,
            scan_hi,
            tof_lo,
            tof_hi,
        }
    }

    #[test]
    fn empty_boxes_yields_no_gate() {
        assert!(DiaMs1Gate::build(&[], 100).is_none());
    }

    #[test]
    fn point_inside_window_is_kept_outside_is_dropped() {
        // One window: scans 10..=20, TOF 1000..=2000.
        let gate = DiaMs1Gate::build(&[boxed(10, 20, 1000, 2000)], 100).unwrap();
        assert!(gate.contains(15, 1500)); // dead center
        assert!(gate.contains(10, 1000)); // inclusive corners
        assert!(gate.contains(20, 2000));
        assert!(!gate.contains(15, 999)); // just left in TOF
        assert!(!gate.contains(15, 2001)); // just right in TOF
        assert!(!gate.contains(9, 1500)); // just below in scan
        assert!(!gate.contains(21, 1500)); // just above in scan
    }

    #[test]
    fn padding_is_already_baked_into_the_box() {
        // The writer pads before building, so a box that has been widened by the
        // equivalent of +500 TOF keeps a point 400 past the nominal edge.
        let nominal_hi = 2000;
        let padded = boxed(10, 20, 1000, nominal_hi + 500);
        let gate = DiaMs1Gate::build(&[padded], 100).unwrap();
        assert!(gate.contains(15, nominal_hi + 400)); // isotope just past the edge
        assert!(!gate.contains(15, nominal_hi + 600)); // beyond the pad
    }

    #[test]
    fn overlapping_windows_on_a_scan_merge() {
        // Two windows overlap at scan 15 in TOF; the gap between them is bridged.
        let gate = DiaMs1Gate::build(&[boxed(10, 20, 1000, 1500), boxed(12, 18, 1490, 2000)], 100)
            .unwrap();
        // A point in the seam between the two TOF ranges is kept (merged).
        assert!(gate.contains(15, 1495));
        assert!(gate.contains(15, 1000));
        assert!(gate.contains(15, 2000));
        // Outside both still dropped.
        assert!(!gate.contains(15, 2500));
    }

    #[test]
    fn keep_mask_matches_per_point_contains() {
        let gate = DiaMs1Gate::build(&[boxed(10, 20, 1000, 2000)], 100).unwrap();
        let scan = [15, 15, 9];
        let tof = [1500, 3000, 1500];
        assert_eq!(gate.keep_mask(&scan, &tof), vec![true, false, false]);
    }
}
