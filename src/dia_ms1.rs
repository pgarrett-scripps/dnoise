//! diaPASEF MS1 out-of-window gate.
//!
//! In diaPASEF the quadrupole tiles the precursor space into isolation windows,
//! each covering an m/z band over a mobility-scan interval. The union of all
//! windows is the set of precursors the method can ever fragment; an MS1 peak
//! outside every window is a precursor that is never isolated, so it can be
//! dropped from the MS1 survey scans.
//!
//! This module is calibration-free: it consumes windows already expressed as
//! **padded integer `(scan, TOF index)` boxes** (the writer does the one-time
//! m/z→TOF and 1/K0→scan conversion, including the physical-unit padding). It
//! builds a per-scan list of merged TOF intervals so each MS1 point is tested with
//! a single binary search.

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
        Some(Self { per_scan })
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

#[cfg(test)]
mod tests {
    use super::*;

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
