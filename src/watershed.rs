//! Watershed-segmentation centroider in integer `(scan, TOF index)` space.
//!
//! Unlike the vertical and halo filters — which produce a *keep mask* over the
//! existing points — the watershed centroider produces a smaller set of *new*
//! points: each output is one watershed group collapsed to a single
//! `(scan, tof, intensity)` triple. The position is the intensity-weighted mean
//! of the group's members, rounded to the integer grid; the intensity is the
//! group's summed intensity (saturated to `u32::MAX`).
//!
//! Points are processed in descending intensity order. For each point:
//!   * If the nearest already-assigned point inside the
//!     `(±box_scan, ±box_mz_idx)` box exists AND that group's seed sits within
//!     `±max_tof_offset` TOF indices of the candidate, join its group (Manhattan
//!     distance, tiebreak prefers the group whose seed had higher intensity —
//!     keeps the watershed boundary stable).
//!   * Else if `intensity >= min_seed_intensity`, open a new group.
//!   * Else drop the point as an orphan; it does NOT enter the grid, so later
//!     weaker points cannot attach to it.
//!
//! The `max_tof_offset` cap is measured against the group's *seed* (not the
//! nearest follower), so a long chain of followers cannot drag the group past
//! its real TOF extent. Groups whose summed intensity is below
//! `min_centroid_total` are dropped at emit time.
//!
//! Ported from `koth_rust`'s `watershed_centroid`, adapted to emit integer
//! `(scan, tof, intensity)` triples (rather than converted float centroids) so
//! the result can be re-encoded as a type-2 `tdf_bin` like every other stage.

use crate::params::WatershedParams;
use std::collections::HashMap;
use std::hash::{BuildHasher, Hasher};

/// Hasher for `u64` composite keys built from packed `(scan_cell, tof_cell)`
/// integer coordinates. SipHash is overkill here (and far too slow); the keys
/// are well-spread but their low bits can have poor entropy, so a SplitMix64
/// finalizer scrambles bits into bucket positions. Ported from `koth_rust`.
#[derive(Default, Clone, Copy)]
pub(crate) struct U64Hasher(u64);

impl Hasher for U64Hasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, _: &[u8]) {
        unreachable!("U64Hasher is only used with u64 keys")
    }
    fn write_u64(&mut self, n: u64) {
        // SplitMix64 finalizer.
        let mut x = n;
        x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        x ^= x >> 31;
        self.0 = x;
    }
}

#[derive(Default, Clone, Copy)]
pub(crate) struct BuildU64Hasher;

impl BuildHasher for BuildU64Hasher {
    type Hasher = U64Hasher;
    fn build_hasher(&self) -> U64Hasher {
        U64Hasher(0)
    }
}

/// Centroid `points` (`(scan, tof, intensity)` triples, any order) into a
/// smaller set of `(scan, tof, intensity)` centroids via watershed segmentation.
///
/// At most `max_centroids` groups are formed (a guard against pathological
/// frames); once the cap is hit, remaining points are dropped. Returns the
/// centroids in arbitrary order — the type-2 encoder regroups by scan/TOF.
pub fn watershed_centroid(
    points: &[(u32, u32, u32)],
    p: &WatershedParams,
    max_centroids: usize,
) -> Vec<(u32, u32, u32)> {
    let n = points.len();
    if n == 0 {
        return Vec::new();
    }

    let cell_w_scan = p.box_scan.max(1);
    let cell_w_tof = p.box_mz_idx.max(1);

    // Sort indices by (intensity desc, scan asc, tof asc) for deterministic
    // tiebreaks across runs.
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_by(|&a, &b| {
        let (sa, ta, ia) = points[a as usize];
        let (sb, tb, ib) = points[b as usize];
        ib.cmp(&ia).then(sa.cmp(&sb)).then(ta.cmp(&tb))
    });

    // Bucket grid: cell side = (box_scan, box_mz_idx), so any in-box neighbour
    // must live in one of the 3x3 cells around the query.
    let mut grid: HashMap<u64, Vec<u32>, BuildU64Hasher> = HashMap::with_hasher(BuildU64Hasher);
    let mut group_id: Vec<i32> = vec![-1; n];
    let mut seed_intensities: Vec<u32> = Vec::new();
    let mut seed_tofs: Vec<u32> = Vec::new();

    let key = |cs: u32, ct: u32| -> u64 { ((cs as u64) << 32) | (ct as u64) };

    for &p_idx in &order {
        if seed_intensities.len() >= max_centroids {
            break;
        }
        let pi = p_idx as usize;
        let (ps, pt, p_int) = points[pi];
        let cell_s = ps / cell_w_scan;
        let cell_t = pt / cell_w_tof;

        let mut best_dist: u32 = u32::MAX;
        let mut best_q: i32 = -1;
        let mut best_seed: u32 = 0;
        let mut found = false;

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
                    for &q_idx in bucket {
                        let qi = q_idx as usize;
                        let (qs, qt, _) = points[qi];
                        let ds_abs = ps.abs_diff(qs);
                        let dt_abs = pt.abs_diff(qt);
                        if ds_abs > p.box_scan || dt_abs > p.box_mz_idx {
                            continue;
                        }
                        // Seed-distance cap: refuse this candidate if joining its
                        // group would put us further than `max_tof_offset` from
                        // that group's seed.
                        let gid = group_id[qi] as usize;
                        if pt.abs_diff(seed_tofs[gid]) > p.max_tof_offset {
                            continue;
                        }
                        let dist = ds_abs + dt_abs;
                        let q_seed = seed_intensities[gid];
                        if !found || dist < best_dist || (dist == best_dist && q_seed > best_seed) {
                            found = true;
                            best_dist = dist;
                            best_q = q_idx as i32;
                            best_seed = q_seed;
                        }
                    }
                }
            }
        }

        if best_q >= 0 {
            // Join existing group.
            let g = group_id[best_q as usize];
            group_id[pi] = g;
            grid.entry(key(cell_s, cell_t)).or_default().push(p_idx);
        } else if p_int as u64 >= p.min_seed_intensity {
            // Promote to new seed.
            let g = seed_intensities.len() as i32;
            group_id[pi] = g;
            seed_intensities.push(p_int);
            seed_tofs.push(pt);
            grid.entry(key(cell_s, cell_t)).or_default().push(p_idx);
        }
        // else: orphan — drop without inserting into the grid.
    }

    // Aggregate centroids: intensity-weighted average position per group.
    let n_groups = seed_intensities.len();
    if n_groups == 0 {
        return Vec::new();
    }
    let mut sum_int = vec![0u64; n_groups];
    let mut sum_scan_w = vec![0.0f64; n_groups];
    let mut sum_tof_w = vec![0.0f64; n_groups];

    for i in 0..n {
        let g = group_id[i];
        if g < 0 {
            continue;
        }
        let g = g as usize;
        let (s, t, int) = points[i];
        let w = int as f64;
        sum_int[g] += int as u64;
        sum_scan_w[g] += s as f64 * w;
        sum_tof_w[g] += t as f64 * w;
    }

    // Emit one integer point per surviving group, merging any groups that round
    // to the same (scan, tof) cell so the encoder never sees duplicate coords
    // within a scan. Intensity saturates at u32::MAX.
    let mut merged: HashMap<u64, u64, BuildU64Hasher> = HashMap::with_hasher(BuildU64Hasher);
    let mut keys_in_order: Vec<u64> = Vec::with_capacity(n_groups);
    for g in 0..n_groups {
        if sum_int[g] < p.min_centroid_total {
            continue;
        }
        let total = sum_int[g] as f64;
        let scan = (sum_scan_w[g] / total).round() as u32;
        let tof = (sum_tof_w[g] / total).round() as u32;
        let k = ((scan as u64) << 32) | (tof as u64);
        let slot = merged.entry(k);
        if let std::collections::hash_map::Entry::Vacant(_) = slot {
            keys_in_order.push(k);
        }
        *merged.entry(k).or_insert(0) += sum_int[g];
    }

    keys_in_order
        .into_iter()
        .map(|k| {
            let scan = (k >> 32) as u32;
            let tof = (k & 0xFFFF_FFFF) as u32;
            let intensity = merged[&k].min(u32::MAX as u64) as u32;
            (scan, tof, intensity)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> WatershedParams {
        WatershedParams {
            box_scan: 5,
            box_mz_idx: 3,
            min_seed_intensity: 0,
            min_centroid_total: 0,
            max_tof_offset: 10,
        }
    }

    /// Count distinct centroids emitted (positions may collide, so this counts
    /// the emitted points, which is what downstream sees).
    fn n_out(points: &[(u32, u32, u32)], p: &WatershedParams) -> usize {
        watershed_centroid(points, p, usize::MAX).len()
    }

    #[test]
    fn empty_in_empty_out() {
        assert_eq!(watershed_centroid(&[], &params(), usize::MAX), Vec::new());
    }

    #[test]
    fn two_separated_peaks_emit_two_centroids() {
        // Two distinct islands, well outside each other's box.
        let pts = vec![
            (10, 100, 10),
            (10, 101, 5),
            (11, 100, 5),
            (40, 200, 9),
            (40, 201, 4),
            (41, 200, 4),
        ];
        assert_eq!(n_out(&pts, &params()), 2);
    }

    #[test]
    fn single_ridge_merges_via_followers() {
        // A 5-point ridge along the scan axis at tof=100. With box_scan=2 the
        // seed at scan=12 reaches scan 10..14 via followers — one group.
        let pts = vec![
            (10, 100, 4),
            (11, 100, 8),
            (12, 100, 10),
            (13, 100, 8),
            (14, 100, 4),
        ];
        let p = WatershedParams {
            box_scan: 2,
            box_mz_idx: 1,
            ..params()
        };
        let out = watershed_centroid(&pts, &p, usize::MAX);
        assert_eq!(out.len(), 1);
        // Symmetric ridge -> weighted mean lands back on the apex scan=12, tof=100.
        assert_eq!(out[0].0, 12);
        assert_eq!(out[0].1, 100);
        // Intensity is the group sum.
        assert_eq!(out[0].2, 4 + 8 + 10 + 8 + 4);
    }

    #[test]
    fn orphan_dropped_below_seed_floor() {
        // A single point below the seed floor with no neighbours: dropped.
        let pts = vec![(10, 100, 1)];
        let p = WatershedParams {
            min_seed_intensity: 5,
            ..params()
        };
        assert_eq!(n_out(&pts, &p), 0);
    }

    #[test]
    fn max_tof_offset_caps_follower_chain() {
        // A staircase along TOF: seed at tof=100 stepping +2. With box_mz_idx=2
        // every adjacent pair is in-box; without a seed cap the whole chain is
        // one group, with a tight cap it splits.
        let pts = vec![
            (10, 100, 10),
            (10, 102, 9),
            (10, 104, 8),
            (10, 106, 7),
            (10, 108, 6),
            (10, 110, 5),
        ];
        let loose = WatershedParams {
            box_scan: 1,
            box_mz_idx: 2,
            max_tof_offset: 100,
            ..params()
        };
        assert_eq!(n_out(&pts, &loose), 1);

        let tight = WatershedParams {
            box_scan: 1,
            box_mz_idx: 2,
            max_tof_offset: 4,
            ..params()
        };
        assert!(n_out(&pts, &tight) > 1);
    }

    #[test]
    fn min_centroid_total_drops_weak_group() {
        // One isolated point whose group sum is below the centroid-total floor.
        let pts = vec![(10, 100, 3)];
        let p = WatershedParams {
            min_centroid_total: 5,
            ..params()
        };
        assert_eq!(n_out(&pts, &p), 0);
    }

    #[test]
    fn colliding_groups_are_merged() {
        // Two points at the identical (scan, tof) collapse to one emitted point
        // with summed intensity — the encoder never sees duplicate coords.
        let pts = vec![(10, 100, 10), (10, 100, 20)];
        let out = watershed_centroid(&pts, &params(), usize::MAX);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0], (10, 100, 30));
    }
}
