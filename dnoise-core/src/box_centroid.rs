//! Greedy small-box centroider in integer `(scan, TOF index)` space.
//!
//! An optional *final* reduction stage. Unlike the watershed centroider — which
//! grows transitively and collapses a whole ion streak into a single point — this
//! consolidates points inside small, fixed, non-transitive boxes. A long mobility
//! streak is therefore tiled into many small centroids, so the streak's *shape*
//! (its mobility profile) survives at coarser resolution while the redundant
//! per-sample points within each box are merged. The intent: shrink point count
//! and preserve m/z (TOF) precision, treating ion-mobility precision as cheaper.
//!
//! Points are processed in descending intensity order. Each not-yet-consumed
//! point seeds a box `[tof ± mz_idx_half_width] x [scan ± scan_half_width]`; all
//! unconsumed points inside are summed into one intensity-weighted
//! `(scan, tof, intensity)` centroid (rounded to the integer grid) and marked
//! consumed. A box whose summed intensity is below `min_centroid_total` is
//! dropped (its points are still consumed) — that floor is the optional
//! denoising knob; at 0 the stage conserves total intensity exactly.

use crate::params::BoxCentroidParams;
use crate::watershed::BuildU64Hasher;
use std::collections::HashMap;

#[inline]
fn key(cs: u32, ct: u32) -> u64 {
    ((cs as u64) << 32) | (ct as u64)
}

/// Consolidate `points` (`(scan, tof, intensity)` triples, any order) into a
/// smaller set of intensity-weighted centroids via greedy small-box merging.
/// Returns centroids in arbitrary order (the type-2 encoder regroups by scan/TOF).
pub fn box_centroid(points: &[(u32, u32, u32)], p: &BoxCentroidParams) -> Vec<(u32, u32, u32)> {
    let n = points.len();
    if n == 0 {
        return Vec::new();
    }
    let cell_w_scan = p.scan_half_width.max(1);
    let cell_w_tof = p.mz_idx_half_width.max(1);

    // Brightest first; deterministic tiebreak on (scan, tof).
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_by(|&a, &b| {
        let (sa, ta, ia) = points[a as usize];
        let (sb, tb, ib) = points[b as usize];
        ib.cmp(&ia).then(sa.cmp(&sb)).then(ta.cmp(&tb))
    });

    // Bucket grid: cell side = box half-widths, so any in-box neighbour lives in
    // one of the 3x3 cells around the seed.
    let mut grid: HashMap<u64, Vec<u32>, BuildU64Hasher> = HashMap::with_hasher(BuildU64Hasher);
    for i in 0..n as u32 {
        let (s, t, _) = points[i as usize];
        grid.entry(key(s / cell_w_scan, t / cell_w_tof))
            .or_default()
            .push(i);
    }

    let mut consumed = vec![false; n];
    let mut out: Vec<(u32, u32, u32)> = Vec::new();

    for &seed in &order {
        let si = seed as usize;
        if consumed[si] {
            continue;
        }
        let (ss, st, _) = points[si];
        let cell_s = ss / cell_w_scan;
        let cell_t = st / cell_w_tof;

        let mut sum_int: u64 = 0;
        let mut sum_scan_w: f64 = 0.0;
        let mut sum_tof_w: f64 = 0.0;
        let mut members: Vec<u32> = Vec::new();

        for ds in -1i64..=1 {
            for dt in -1i64..=1 {
                let cs = match (cell_s as i64).checked_add(ds) {
                    Some(v) if v >= 0 => v as u32,
                    _ => continue,
                };
                let ct = match (cell_t as i64).checked_add(dt) {
                    Some(v) if v >= 0 => v as u32,
                    _ => continue,
                };
                if let Some(bucket) = grid.get(&key(cs, ct)) {
                    for &q in bucket {
                        let qi = q as usize;
                        if consumed[qi] {
                            continue;
                        }
                        let (qs, qt, qint) = points[qi];
                        if ss.abs_diff(qs) > p.scan_half_width
                            || st.abs_diff(qt) > p.mz_idx_half_width
                        {
                            continue;
                        }
                        members.push(q);
                        let w = qint as f64;
                        sum_int += qint as u64;
                        sum_scan_w += qs as f64 * w;
                        sum_tof_w += qt as f64 * w;
                    }
                }
            }
        }

        // Everything in the box is absorbed (consumed), even if the box is then
        // dropped by the intensity floor.
        for &q in &members {
            consumed[q as usize] = true;
        }
        if sum_int == 0 || sum_int < p.min_centroid_total {
            continue;
        }
        let scan = (sum_scan_w / sum_int as f64).round() as u32;
        let tof = (sum_tof_w / sum_int as f64).round() as u32;
        out.push((scan, tof, sum_int.min(u32::MAX as u64) as u32));
    }

    // Merge any centroids that round to the same (scan, tof) cell so the encoder
    // never sees duplicate coordinates within a scan.
    let mut merged: HashMap<u64, u64, BuildU64Hasher> = HashMap::with_hasher(BuildU64Hasher);
    let mut keys_in_order: Vec<u64> = Vec::with_capacity(out.len());
    for &(s, t, i) in &out {
        let k = key(s, t);
        if !merged.contains_key(&k) {
            keys_in_order.push(k);
        }
        *merged.entry(k).or_insert(0) += i as u64;
    }
    keys_in_order
        .into_iter()
        .map(|k| {
            let scan = (k >> 32) as u32;
            let tof = (k & 0xFFFF_FFFF) as u32;
            (scan, tof, merged[&k].min(u32::MAX as u64) as u32)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> BoxCentroidParams {
        BoxCentroidParams {
            mz_idx_half_width: 2,
            scan_half_width: 2,
            min_centroid_total: 0,
        }
    }

    fn total(pts: &[(u32, u32, u32)]) -> u64 {
        pts.iter().map(|&(_, _, i)| i as u64).sum()
    }

    #[test]
    fn empty_in_empty_out() {
        assert!(box_centroid(&[], &params()).is_empty());
    }

    #[test]
    fn colocated_points_merge_and_conserve_intensity() {
        let pts = vec![(10, 100, 10), (10, 100, 20), (11, 101, 5)];
        let out = box_centroid(&pts, &params());
        assert_eq!(out.len(), 1); // all within ±2
        assert_eq!(out[0].2, 35); // intensity conserved
    }

    #[test]
    fn long_streak_is_tiled_not_collapsed() {
        // A 15-scan vertical streak at tof=100. Watershed would make ONE point;
        // small ±2-scan boxes tile it into several centroids that preserve the
        // streak's extent, and total intensity is conserved.
        let pts: Vec<(u32, u32, u32)> = (0..15).map(|s| (10 + s, 100, 10)).collect();
        let out = box_centroid(&pts, &params());
        assert!(
            out.len() >= 3,
            "expected the streak tiled into several centroids, got {}",
            out.len()
        );
        assert!(out.len() < pts.len(), "should still reduce the point count");
        assert_eq!(total(&out), total(&pts)); // intensity conserved at floor 0
        assert!(out.iter().all(|&(_, t, _)| t == 100)); // m/z preserved
    }

    #[test]
    fn min_total_drops_weak_isolated_box() {
        let pts = vec![(10, 100, 3), (50, 500, 100)];
        let p = BoxCentroidParams {
            min_centroid_total: 10,
            ..params()
        };
        let out = box_centroid(&pts, &p);
        assert_eq!(out.len(), 1); // the weak isolated point is dropped
        assert_eq!(out[0], (50, 500, 100));
    }
}
