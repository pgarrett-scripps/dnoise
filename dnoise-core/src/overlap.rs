//! Feature-level ("any overlap") mode for the MS1 acquisition gates.
//!
//! The selection-polygon gate ([`crate::polygon`]) and the diaPASEF MS1 window
//! gate ([`crate::dia_ms1`]) are point masks: each surviving MS1 point is kept or
//! dropped on its own `(scan, TOF)` position, so a feature that straddles a gate
//! edge is cut in two, and the pads exist only to soften that cut.
//!
//! In overlap mode the decision is made per feature instead. The points that
//! survived the streak filter are grouped into connected features (points within
//! `tof_tol` TOF indices of each other and at most `scan_gap` scans apart, the
//! same adjacency the streak filter bridges), and a feature is kept whole when
//! **any** of its points lies inside the gate. A precursor whose mobility spread
//! or apex falls partly outside the region is therefore kept intact, and a feature
//! lying wholly outside is dropped, with no fixed pad deciding where to cut.
//!
//! A kept feature is kept over its whole extent: there is no limit on how far it
//! may run beyond the mobility range of its inside points. (0.4.0 capped this at
//! an `overlap_reach` of 0.1 1/K0; 0.5.0 removed the cap.)

/// Extend a point-level gate mask to whole features.
///
/// * `scan`, `tof` — parallel per-point coordinates of the frame.
/// * `keep` — the current keep mask; only kept points form features.
/// * `inside` — the gate's point mask (`true` = inside the gate).
/// * `tof_tol` — TOF-index distance that links two points (the streak filter's
///   column half-width).
/// * `scan_gap` — scan distance that links two points (the streak filter's
///   largest bridged gap plus one).
///
/// Returns a mask that is `true` for kept points belonging to a feature with at
/// least one inside point. Points not in `keep` are always `false`.
pub fn extend_to_features(
    scan: &[u32],
    tof: &[u32],
    keep: &[bool],
    inside: &[bool],
    tof_tol: u32,
    scan_gap: u32,
) -> Vec<bool> {
    let n = scan.len();
    let mut out = vec![false; n];

    // Kept points ordered by (scan, tof), so each scan is a contiguous run sorted
    // by TOF and neighbours can be found by binary search.
    let mut idx: Vec<usize> = (0..n).filter(|&i| keep[i]).collect();
    if idx.is_empty() {
        return out;
    }
    idx.sort_unstable_by_key(|&i| (scan[i], tof[i]));
    let m = idx.len();

    // Start offset of each distinct scan in `idx`.
    let mut starts: Vec<(u32, usize)> = Vec::new();
    for (k, &i) in idx.iter().enumerate() {
        if starts.last().is_none_or(|&(s, _)| s != scan[i]) {
            starts.push((scan[i], k));
        }
    }
    let run_of = |si: usize| -> (usize, usize) {
        let lo = starts[si].1;
        let hi = starts.get(si + 1).map_or(m, |&(_, k)| k);
        (lo, hi)
    };

    let mut parent: Vec<usize> = (0..m).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    fn union(p: &mut [usize], a: usize, b: usize) {
        let (ra, rb) = (find(p, a), find(p, b));
        if ra != rb {
            p[ra] = rb;
        }
    }

    for si in 0..starts.len() {
        let (lo, hi) = run_of(si);
        let s = starts[si].0;
        // Same scan: TOF-sorted, so only the next point can be the nearest link.
        for k in lo..hi.saturating_sub(1) {
            if tof[idx[k + 1]] - tof[idx[k]] <= tof_tol {
                union(&mut parent, k, k + 1);
            }
        }
        // Later scans within the gap.
        for (sj, &(s2, _)) in starts.iter().enumerate().skip(si + 1) {
            if s2 - s > scan_gap {
                break;
            }
            let (lo2, hi2) = run_of(sj);
            let row: Vec<u32> = idx[lo2..hi2].iter().map(|&i| tof[i]).collect();
            for k in lo..hi {
                let t = tof[idx[k]];
                let a = row.partition_point(|&x| x < t.saturating_sub(tof_tol));
                let b = row.partition_point(|&x| x <= t.saturating_add(tof_tol));
                for j in a..b {
                    union(&mut parent, k, lo2 + j);
                }
            }
        }
    }

    // Features (roots) with at least one inside point.
    let mut hit = vec![false; m];
    for k in 0..m {
        if inside[idx[k]] {
            let r = find(&mut parent, k);
            hit[r] = true;
        }
    }
    for k in 0..m {
        let r = find(&mut parent, k);
        out[idx[k]] = hit[r];
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn feature_touching_the_gate_is_kept_whole() {
        // One streak at TOF 100 over scans 0..6; only scan 5 is inside the gate.
        let scan = vec![0, 1, 2, 3, 4, 5];
        let tof = vec![100; 6];
        let keep = vec![true; 6];
        let mut inside = vec![false; 6];
        inside[5] = true;
        assert_eq!(
            extend_to_features(&scan, &tof, &keep, &inside, 3, 3),
            vec![true; 6]
        );
    }

    #[test]
    fn feature_wholly_outside_is_dropped() {
        // Two streaks far apart in TOF; only the first touches the gate.
        let scan = vec![0, 1, 2, 0, 1, 2];
        let tof = vec![100, 101, 100, 500, 500, 501];
        let keep = vec![true; 6];
        let inside = vec![true, false, false, false, false, false];
        assert_eq!(
            extend_to_features(&scan, &tof, &keep, &inside, 3, 3),
            vec![true, true, true, false, false, false]
        );
    }

    #[test]
    fn gap_beyond_the_bridge_splits_features() {
        // Scans 0 and 10 at the same TOF are not linked with scan_gap 3.
        let scan = vec![0, 10];
        let tof = vec![100, 100];
        let keep = vec![true, true];
        let inside = vec![true, false];
        assert_eq!(
            extend_to_features(&scan, &tof, &keep, &inside, 3, 3),
            vec![true, false]
        );
    }

    #[test]
    fn a_long_feature_is_kept_whole_however_far_it_runs_past_the_gate() {
        // One streak over scans 0..200; only scan 199 is inside. No reach cap.
        let scan: Vec<u32> = (0..200).collect();
        let tof = vec![100; 200];
        let keep = vec![true; 200];
        let mut inside = vec![false; 200];
        inside[199] = true;
        assert_eq!(
            extend_to_features(&scan, &tof, &keep, &inside, 3, 3),
            vec![true; 200]
        );
    }

    #[test]
    fn removed_points_neither_link_nor_survive() {
        let scan = vec![0, 1, 2];
        let tof = vec![100, 100, 100];
        let keep = vec![true, false, true];
        let inside = vec![true, true, false];
        // Scan 0 and 2 are 2 apart: still linked with scan_gap 3.
        assert_eq!(
            extend_to_features(&scan, &tof, &keep, &inside, 3, 3),
            vec![true, false, true]
        );
        // With scan_gap 1 the removed middle point no longer bridges them.
        assert_eq!(
            extend_to_features(&scan, &tof, &keep, &inside, 3, 1),
            vec![true, false, false]
        );
    }
}
