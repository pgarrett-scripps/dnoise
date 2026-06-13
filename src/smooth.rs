//! Box-averaging smoother in integer `(scan, TOF index)` space.
//!
//! Rewrites each point's intensity with the mean intensity of all points inside
//! a `(±scan_half_width, ±mz_idx_half_width)` box around it (inclusive of the
//! point itself). Coordinates are never moved; only intensities change. The
//! mean is taken over the points that actually exist in the box (a sparse mean,
//! not a dense one that treats empty grid cells as zero), so an isolated point
//! keeps its own value.
//!
//! This runs after the halo filter and before the watershed centroider: the
//! watershed seeds and grows regions in descending intensity order, so noise-
//! driven local intensity spikes can split a single ion into several centroids.
//! Averaging flattens those spikes, making the seed ordering — and thus the
//! centroids — more stable. It is intensity-only and order-independent (each
//! pass reads a snapshot of the previous pass's intensities).

use crate::params::SmoothParams;

/// Return `points` with each intensity replaced by the box-mean of the points
/// in its `(±scan_half_width, ±mz_idx_half_width)` neighbourhood. `scan`/`tof`
/// are preserved; the result is in the same order as the input. Iterated
/// `iterations` times (each pass feeds the next).
pub fn box_average(
    points: &[(u32, u32, u32)],
    num_scans: usize,
    p: &SmoothParams,
) -> Vec<(u32, u32, u32)> {
    let n = points.len();
    if n == 0 || num_scans == 0 || p.iterations == 0 {
        return points.to_vec();
    }

    // Sort by (scan, tof) so each scan is a contiguous, TOF-sorted block. The
    // coordinates never change, so this layout is reused across all iterations.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| (points[i].0, points[i].1));
    let scan_s: Vec<u32> = order.iter().map(|&i| points[i].0).collect();
    let tof_s: Vec<u32> = order.iter().map(|&i| points[i].1).collect();
    let mut int_s: Vec<u64> = order.iter().map(|&i| points[i].2 as u64).collect();

    // Per-scan block boundaries into the sorted arrays.
    let mut block_start = vec![0usize; num_scans];
    let mut block_len = vec![0usize; num_scans];
    let mut k = 0;
    while k < n {
        let s = scan_s[k] as usize;
        let start = k;
        while k < n && scan_s[k] as usize == s {
            k += 1;
        }
        block_start[s] = start;
        block_len[s] = k - start;
    }

    let sw = p.scan_half_width as i64;
    let tw = p.mz_idx_half_width;
    for _ in 0..p.iterations {
        let mut next = vec![0u64; n];
        for i in 0..n {
            let si = scan_s[i] as i64;
            let ti = tof_s[i];
            let lo_scan = (si - sw).max(0) as usize;
            let hi_scan = (si + sw).min(num_scans as i64 - 1) as usize;
            let lo_tof = ti.saturating_sub(tw);
            let hi_tof = ti.saturating_add(tw);

            let mut sum = 0u64;
            let mut cnt = 0u64;
            for v in lo_scan..=hi_scan {
                let bl = block_len[v];
                if bl == 0 {
                    continue;
                }
                let bs = block_start[v];
                let be = bs + bl;
                // First index in this block with tof >= lo_tof.
                let mut left = bs;
                let mut right = be;
                while left < right {
                    let mid = (left + right) / 2;
                    if tof_s[mid] < lo_tof {
                        left = mid + 1;
                    } else {
                        right = mid;
                    }
                }
                let mut j = left;
                while j < be && tof_s[j] <= hi_tof {
                    sum += int_s[j];
                    cnt += 1;
                    j += 1;
                }
            }
            // Round-half-up mean; cnt >= 1 because the box always contains self.
            next[i] = if cnt > 0 {
                (sum + cnt / 2) / cnt
            } else {
                int_s[i]
            };
        }
        int_s = next;
    }

    // Write the smoothed intensities back into the original point order.
    let mut out = points.to_vec();
    for (sorted_k, &orig) in order.iter().enumerate() {
        out[orig].2 = int_s[sorted_k].min(u32::MAX as u64) as u32;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params(tw: u32, sw: usize, it: usize) -> SmoothParams {
        SmoothParams {
            mz_idx_half_width: tw,
            scan_half_width: sw,
            iterations: it,
        }
    }

    #[test]
    fn empty_and_noop() {
        assert_eq!(box_average(&[], 700, &params(1, 1, 1)), Vec::new());
        let pts = vec![(10, 100, 5)];
        // iterations = 0 is a no-op.
        assert_eq!(box_average(&pts, 700, &params(1, 1, 0)), pts);
    }

    #[test]
    fn isolated_point_keeps_value() {
        // No neighbour within the box -> mean over {self} = self.
        let pts = vec![(10, 100, 42), (10, 500, 7)];
        let out = box_average(&pts, 700, &params(1, 1, 1));
        assert_eq!(out[0].2, 42);
        assert_eq!(out[1].2, 7);
    }

    #[test]
    fn two_neighbours_average_to_their_mean() {
        // Same scan, adjacent TOF, within tw=1: both become mean(10, 20) = 15.
        let pts = vec![(10, 100, 10), (10, 101, 20)];
        let out = box_average(&pts, 700, &params(1, 0, 1));
        assert_eq!(out[0].2, 15);
        assert_eq!(out[1].2, 15);
    }

    #[test]
    fn coordinates_are_preserved() {
        let pts = vec![(10, 100, 10), (11, 100, 30), (10, 101, 50)];
        let out = box_average(&pts, 700, &params(1, 1, 1));
        let coords: Vec<(u32, u32)> = out.iter().map(|&(s, t, _)| (s, t)).collect();
        assert_eq!(coords, vec![(10, 100), (11, 100), (10, 101)]);
    }

    #[test]
    fn spike_is_flattened() {
        // A bright spike flanked by weak points: averaging pulls the spike down
        // and the neighbours up, so the local maximum is less dominant.
        let pts = vec![(10, 99, 2), (10, 100, 100), (10, 101, 2)];
        let out = box_average(&pts, 700, &params(1, 0, 1));
        // center: mean(2,100,2)=34.67 -> 35 ; edges: mean of their 2-point boxes.
        assert_eq!(out[1].2, 35);
        assert!(out[1].2 < 100 && out[0].2 > 2 && out[2].2 > 2);
    }

    #[test]
    fn window_respects_box_bounds() {
        // tw=1 so TOF 100 and 103 are NOT in each other's box: each averages
        // only with itself -> unchanged.
        let pts = vec![(10, 100, 10), (10, 103, 40)];
        let out = box_average(&pts, 700, &params(1, 0, 1));
        assert_eq!(out[0].2, 10);
        assert_eq!(out[1].2, 40);
    }
}
