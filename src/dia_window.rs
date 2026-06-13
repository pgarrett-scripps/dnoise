//! diaPASEF isolation-window handling for MS/MS frames.
//!
//! In diaPASEF the quadrupole steps through a set of `(mobility, m/z)` isolation
//! windows per cycle; each window occupies a contiguous mobility-scan interval
//! `[begin, end)` (read from `DiaFrameMsMsWindows`, see [`crate::tdf`]). Two things
//! follow that the whole-frame MS1 filter cannot see:
//!
//! 1. **Out-of-window signal.** Points whose scan falls in no window were never
//!    isolated — [`in_window_mask`] gates them out.
//! 2. **Cross-window linking.** The vertical filter ([`crate::filter`]) builds
//!    mobility runs by scan adjacency, so a TOF column straddling a window
//!    boundary is wrongly fused into one feature even though the two windows
//!    isolated unrelated precursor m/z bands. [`filter_per_window`] runs the
//!    filter independently inside each window's scan slice to prevent this.

use crate::filter::filter_iterated;
use crate::frame::FlatFrame;
use crate::halo::horizontal_halo_keep_mask;
use crate::params::{FilterParams, HaloParams};

/// Keep mask (in `scans` order) selecting points whose scan lies inside any
/// `[begin - scan_pad, end + scan_pad)` interval. `intervals` must be sorted by
/// `begin` and non-overlapping (as produced by [`crate::tdf::read_dia_windows`]);
/// an empty list keeps nothing. `scan_pad` widens each window symmetrically to
/// tolerate signal a few scans past an isolation edge.
pub fn in_window_mask(scans: &[u32], intervals: &[(u32, u32)], scan_pad: u32) -> Vec<bool> {
    scans
        .iter()
        .map(|&s| {
            // First interval whose padded end is strictly after `s`; the point is
            // in-window iff it also clears that interval's padded begin.
            let i = intervals.partition_point(|&(_, end)| end.saturating_add(scan_pad) <= s);
            i < intervals.len() && s >= intervals[i].0.saturating_sub(scan_pad)
        })
        .collect()
}

/// Run the MS/MS vertical filter — and, when `halo` is `Some`, the horizontal-halo
/// filter on its survivors — independently within each isolation window's scan
/// interval, returning a keep mask in `frame`'s point order. Points outside every
/// window are dropped (no interval visits them), so this also performs the
/// [`in_window_mask`] gate. Scans are remapped to a local `0..end-begin` origin
/// inside each window so the filter's per-scan profile is tight and run adjacency
/// cannot reach across a boundary.
pub fn filter_per_window(
    frame: &FlatFrame,
    intervals: &[(u32, u32)],
    params: &FilterParams,
    halo: Option<&HaloParams>,
) -> Vec<bool> {
    let mut keep = vec![false; frame.len()];
    for &(begin, end) in intervals {
        if end <= begin {
            continue;
        }
        let local: Vec<usize> = (0..frame.len())
            .filter(|&i| frame.scan[i] >= begin && frame.scan[i] < end)
            .collect();
        if local.is_empty() {
            continue;
        }
        let num_scans = (end - begin) as usize;
        let sub = FlatFrame {
            frame_id: frame.frame_id,
            num_scans,
            scan: local.iter().map(|&i| frame.scan[i] - begin).collect(),
            tof: local.iter().map(|&i| frame.tof[i]).collect(),
            intensity: local.iter().map(|&i| frame.intensity[i]).collect(),
        };

        let mut mask = filter_iterated(&sub, params);
        if let Some(hp) = halo {
            apply_halo(&sub, hp, &mut mask);
        }
        for (j, &orig) in local.iter().enumerate() {
            if mask[j] {
                keep[orig] = true;
            }
        }
    }
    keep
}

/// Turn off `mask` for points the horizontal-halo filter removes (operates on the
/// currently-kept points only). Mirrors `writer::apply_halo`, duplicated here to
/// keep this module independent of the writer's private helper.
fn apply_halo(frame: &FlatFrame, hp: &HaloParams, mask: &mut [bool]) {
    let idx: Vec<usize> = (0..frame.len()).filter(|&i| mask[i]).collect();
    if idx.is_empty() {
        return;
    }
    let scan: Vec<u32> = idx.iter().map(|&i| frame.scan[i]).collect();
    let tof: Vec<u32> = idx.iter().map(|&i| frame.tof[i]).collect();
    let inten: Vec<u32> = idx.iter().map(|&i| frame.intensity[i]).collect();
    let hmask = horizontal_halo_keep_mask(&scan, &tof, &inten, frame.num_scans, hp);
    for (k, &i) in idx.iter().enumerate() {
        if !hmask[k] {
            mask[i] = false;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn in_window_mask_basic() {
        // Two windows [10,20) and [30,40); scans on edges and in the gap.
        let scans = [5, 10, 15, 19, 20, 25, 30, 39, 40];
        let iv = [(10u32, 20u32), (30, 40)];
        let keep = in_window_mask(&scans, &iv, 0);
        assert_eq!(
            keep,
            vec![false, true, true, true, false, false, true, true, false]
        );
    }

    #[test]
    fn in_window_mask_pad_extends_edges() {
        let scans = [8, 9, 20, 21];
        let iv = [(10u32, 20u32)];
        // pad=2: [8,22) effectively. 8 in, 9 in, 20 in, 21 in.
        assert_eq!(in_window_mask(&scans, &iv, 2), vec![true, true, true, true]);
        // pad=0: only 8 and 9 below begin, 20/21 at/after end -> all dropped except none.
        assert_eq!(
            in_window_mask(&scans, &iv, 0),
            vec![false, false, false, false]
        );
    }

    #[test]
    fn empty_intervals_keep_nothing() {
        assert_eq!(
            in_window_mask(&[1, 2, 3], &[], 0),
            vec![false, false, false]
        );
    }

    #[test]
    fn per_window_does_not_link_across_boundary() {
        // A single TOF column occupied continuously across the boundary between
        // window [0,5) and window [5,10). min_feature_length=4 means neither
        // window's 3-scan slice survives on its own, but the fused 6-scan run
        // would. Per-window filtering must NOT keep it.
        let scan: Vec<u32> = (2..8).collect(); // scans 2..7 -> 3 in each window
        let frame = FlatFrame {
            frame_id: 1,
            num_scans: 10,
            scan: scan.clone(),
            tof: vec![1000; scan.len()],
            intensity: vec![100; scan.len()],
        };
        let params = FilterParams {
            min_feature_length: 4,
            max_internal_gap: 0,
            num_iterations: 1,
            ..FilterParams::default()
        };
        let iv = [(0u32, 5u32), (5, 10)];
        let keep = filter_per_window(&frame, &iv, &params, None);
        assert!(
            keep.iter().all(|&k| !k),
            "boundary-straddling run must split"
        );

        // Whole-frame filtering, by contrast, fuses the run and keeps it.
        let whole = filter_iterated(&frame, &params);
        assert!(whole.iter().all(|&k| k));
    }

    #[test]
    fn per_window_keeps_self_contained_feature() {
        // A 5-scan run fully inside window [5,12) survives min_feature_length=4.
        let scan: Vec<u32> = (6..11).collect();
        let frame = FlatFrame {
            frame_id: 1,
            num_scans: 20,
            scan: scan.clone(),
            tof: vec![1000; scan.len()],
            intensity: vec![100; scan.len()],
        };
        let params = FilterParams {
            min_feature_length: 4,
            max_internal_gap: 0,
            num_iterations: 1,
            ..FilterParams::default()
        };
        let keep = filter_per_window(&frame, &[(5, 12)], &params, None);
        assert!(keep.iter().all(|&k| k));
    }
}
