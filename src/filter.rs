//! Stage 1 of ALGORITHM.md: the iterative vertical-IM feature filter.
//!
//! A real ion produces a vertical streak in `(tof_index × scan)` space. The
//! filter walks each unique TOF index, sums the IM profile in a small TOF
//! window around it, and keeps only points belonging to long-enough,
//! intense-enough vertical runs. The whole pass is iterated: each pass operates
//! on the survivors of the previous one.

use crate::frame::FlatFrame;
use crate::params::FilterParams;

/// One pass of the vertical-IM filter. Returns a per-point keep mask in the
/// input point order.
pub fn filter_once(frame: &FlatFrame, p: &FilterParams) -> Vec<bool> {
    let n = frame.len();
    let mut keep = vec![false; n];
    if n == 0 {
        return keep;
    }
    let num_scans = frame.num_scans;

    // Sort point indices by TOF so window lookups and unique-TOF blocks are contiguous.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by_key(|&i| frame.tof[i]);
    let sorted_tof: Vec<u32> = order.iter().map(|&i| frame.tof[i]).collect();

    let w = p.mz_half_width;
    let mut profile = vec![0u64; num_scans]; // summed intensity per scan in the current window
    let mut touched: Vec<usize> = Vec::new(); // scans we incremented, for cheap reset
    // Reused across TOF blocks (cleared, not reallocated): the filter visits one
    // block per unique TOF index, so per-block allocation here is the dominant
    // malloc churn in the hot path.
    let mut occupied: Vec<usize> = Vec::new();
    let mut kept_spans: Vec<(usize, usize)> = Vec::new();

    let mut k = 0;
    while k < n {
        let c = sorted_tof[k];

        // Block of points sharing this exact TOF index (they see the same window).
        let block_lo = k;
        let mut block_hi = k + 1;
        while block_hi < n && sorted_tof[block_hi] == c {
            block_hi += 1;
        }

        // Window [c - w, c + w] located by binary search on the sorted TOF array.
        let lo_val = c.saturating_sub(w);
        let hi_val = c.saturating_add(w);
        let w_lo = sorted_tof.partition_point(|&t| t < lo_val);
        let w_hi = sorted_tof.partition_point(|&t| t <= hi_val);

        // Build the IM profile over the window.
        for &idx in &order[w_lo..w_hi] {
            let s = frame.scan[idx] as usize;
            if profile[s] == 0 {
                touched.push(s);
            }
            profile[s] += frame.intensity[idx] as u64;
        }

        // Occupied scans: positive AND clearing the per-scan floor. Sorted ascending.
        occupied.clear();
        occupied.extend(
            touched
                .iter()
                .copied()
                .filter(|&s| profile[s] > 0 && profile[s] >= p.min_window_intensity),
        );
        occupied.sort_unstable();

        // Gap-closed runs: break where the gap exceeds `max_internal_gap` empty scans.
        kept_spans.clear();
        if !occupied.is_empty() {
            let mut run_start = occupied[0];
            let mut prev = occupied[0];
            let close = |run_start: usize,
                         run_end: usize,
                         profile: &[u64],
                         spans: &mut Vec<(usize, usize)>| {
                // Length is the number of OCCUPIED scans in the run, not the
                // gap-inclusive end-to-end span: scans bridged by `max_internal_gap`
                // (profile == 0, or below the window floor) do not count, so
                // `min_feature_length` means "this many points actually seen". The
                // kept range stays [run_start, run_end]; only acceptance changes.
                let occupied_count = profile[run_start..=run_end]
                    .iter()
                    .filter(|&&v| v > 0 && v >= p.min_window_intensity)
                    .count();
                if occupied_count < p.min_feature_length {
                    return;
                }
                let total: u64 = profile[run_start..=run_end].iter().sum();
                if total >= p.min_feature_intensity {
                    spans.push((run_start, run_end));
                }
            };
            for &s in &occupied[1..] {
                if s - prev > p.max_internal_gap + 1 {
                    close(run_start, prev, &profile, &mut kept_spans);
                    run_start = s;
                }
                prev = s;
            }
            close(run_start, prev, &profile, &mut kept_spans);
        }

        // Keep points at TOF `c` whose scan falls inside a kept run's span.
        // `kept_spans` is ascending and non-overlapping, so locate the candidate
        // span by binary search (first span ending at or after `s`) instead of a
        // linear scan over every run.
        for &idx in &order[block_lo..block_hi] {
            let s = frame.scan[idx] as usize;
            let i = kept_spans.partition_point(|&(_, re)| re < s);
            if i < kept_spans.len() && s >= kept_spans[i].0 {
                keep[idx] = true;
            }
        }

        // Reset only the scans we touched.
        for &s in &touched {
            profile[s] = 0;
        }
        touched.clear();

        k = block_hi;
    }
    keep
}

/// Iterate [`filter_once`] over its own survivors `num_iterations` times,
/// composing the mask back into the original point order each pass.
pub fn filter_iterated(frame: &FlatFrame, p: &FilterParams) -> Vec<bool> {
    let n = frame.len();
    let mut cumulative = vec![true; n];
    if n == 0 || p.num_iterations == 0 {
        return cumulative;
    }
    for _ in 0..p.num_iterations {
        let active: Vec<usize> = (0..n).filter(|&i| cumulative[i]).collect();
        if active.is_empty() {
            break;
        }
        let sub = FlatFrame {
            frame_id: frame.frame_id,
            num_scans: frame.num_scans,
            scan: active.iter().map(|&i| frame.scan[i]).collect(),
            tof: active.iter().map(|&i| frame.tof[i]).collect(),
            intensity: active.iter().map(|&i| frame.intensity[i]).collect(),
        };
        let mask = filter_once(&sub, p);
        let mut next = vec![false; n];
        for (j, &orig) in active.iter().enumerate() {
            if mask[j] {
                next[orig] = true;
            }
        }
        let any = next.iter().any(|&b| b);
        cumulative = next;
        if !any {
            break;
        }
    }
    cumulative
}
