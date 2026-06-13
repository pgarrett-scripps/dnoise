//! Horizontal (m/z-axis) halo filter: remove the weak halo flanking bright ions,
//! left/right only.
//!
//! High-intensity ions are flanked by a halo of weak peaks (from charge
//! interactions or detector effects) that do not resolve to precise m/z values.
//! A real ion forms a *vertical streak* along the ion-mobility axis (the same
//! TOF index across many scans), so this filter only removes peaks to the left
//! and right (different TOF index) of a bright neighbor and never above/below
//! (same TOF index).
//!
//! For each peak it computes a reference intensity — the maximum intensity in
//! the surrounding `±scan_half_width × ±mz_idx_half_width` box **excluding the
//! peak's own m/z column** — and drops the peak when its intensity falls below
//! `peak_fraction` of that reference. Excluding the own column is what keeps the
//! vertical streak from ever counting against a point. Operates entirely in
//! integer `(scan, TOF index)` space — no calibration. Ported from `tdfpy`'s
//! `HorizontalHaloFilter`.

use crate::params::HaloParams;

/// Keep-mask for the horizontal-halo filter. `scan` and `tof` are the integer
/// ion-mobility and TOF indices of each point and `intensity` its intensity;
/// all three are parallel and in the same order. Returns a per-point keep mask
/// in that order.
pub fn horizontal_halo_keep_mask(
    scan: &[u32],
    tof: &[u32],
    intensity: &[u32],
    num_scans: usize,
    p: &HaloParams,
) -> Vec<bool> {
    let n = intensity.len();
    let mut keep = vec![true; n];
    if n == 0 || num_scans == 0 {
        return keep;
    }

    // Sort by (scan, tof) so each scan is a contiguous, TOF-sorted block.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| (scan[a], tof[a]).cmp(&(scan[b], tof[b])));
    let scan_s: Vec<u32> = order.iter().map(|&i| scan[i]).collect();
    let tof_s: Vec<u32> = order.iter().map(|&i| tof[i]).collect();
    let int_s: Vec<u32> = order.iter().map(|&i| intensity[i]).collect();

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
    let mw = p.mz_idx_half_width;
    for i in 0..n {
        let si = scan_s[i] as i64;
        let mi = tof_s[i];
        let lo_scan = (si - sw).max(0) as usize;
        let hi_scan = (si + sw).min(num_scans as i64 - 1) as usize;
        let lo_mz = mi.saturating_sub(mw);
        let hi_mz = mi.saturating_add(mw);

        // Max intensity over the box, excluding the point's own TOF column.
        let mut best = 0u32;
        for v in lo_scan..=hi_scan {
            let bl = block_len[v];
            if bl == 0 {
                continue;
            }
            let bs = block_start[v];
            let be = bs + bl;
            // First index in this block with tof >= lo_mz.
            let mut left = bs;
            let mut right = be;
            while left < right {
                let mid = (left + right) / 2;
                if tof_s[mid] < lo_mz {
                    left = mid + 1;
                } else {
                    right = mid;
                }
            }
            let mut j = left;
            while j < be && tof_s[j] <= hi_mz {
                if tof_s[j] != mi && int_s[j] > best {
                    best = int_s[j];
                }
                j += 1;
            }
        }

        // Keep unless strictly below the fraction of the off-column max.
        let alive = int_s[i] as f64 >= p.peak_fraction * best as f64;
        keep[order[i]] = alive;
    }
    keep
}
