//! ddaPASEF MS/MS denoising (precursor-centric).
//!
//! In ddaPASEF a precursor's fragments occupy one ion-mobility scan window
//! (`PasefFrameMsMsInfo.ScanNumBegin..ScanNumEnd`) that is re-isolated across one
//! or more frames. This module combines all of a precursor's fragment scans into
//! a single spectrum (summing intensity at aligned `(scan, TOF)`), denoises that
//! combined spectrum with the vertical + halo filters, and produces a keep set of
//! `(scan, TOF)` that is then applied back to the precursor's points in every
//! individual frame.
//!
//! Public surface is the pure [`combine_and_filter`]; the [`MsmsKeep`] builder
//! and frame map-back are crate-internal plumbing used by [`crate::writer`].

use crate::error::{DnoiseError, Result};
use crate::filter::filter_iterated;
use crate::frame::FlatFrame;
use crate::halo::horizontal_halo_keep_mask;
use crate::params::{FilterParams, HaloParams, MsmsFilterParams};
use crate::tdf::{FrameMeta, PasefWindow};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use timsrust::readers::FrameReader;

/// Pack an absolute `(scan, tof)` into a single key.
#[inline]
fn key(scan: u32, tof: u32) -> u64 {
    ((scan as u64) << 32) | tof as u64
}

/// Combine a precursor's fragment points and return the kept `(scan, tof)` keys.
///
/// `points` are `(scan, tof, intensity)` in absolute scan coordinates; they are
/// summed per `(scan, tof)` into one combined spectrum (scans remapped to a local
/// `0..num_local_scans` window starting at `scan_begin`), filtered with the
/// vertical filter `params` and, when `halo` is `Some`, the horizontal-halo
/// filter on the survivors. Returned keys are in absolute scan coordinates.
pub fn combine_and_filter(
    points: &[(u32, u32, u32)],
    scan_begin: u32,
    num_local_scans: usize,
    params: &FilterParams,
    halo: Option<&HaloParams>,
) -> HashSet<u64> {
    let mut out = HashSet::new();
    if points.is_empty() || num_local_scans == 0 {
        return out;
    }

    // Sum intensities per local (scan', tof).
    let mut acc: HashMap<u64, u64> = HashMap::new();
    for &(s, t, i) in points {
        let sl = s.saturating_sub(scan_begin);
        *acc.entry(key(sl, t)).or_insert(0) += i as u64;
    }
    let n = acc.len();
    let mut scan = Vec::with_capacity(n);
    let mut tof = Vec::with_capacity(n);
    let mut intensity = Vec::with_capacity(n);
    for (&k, &sum) in &acc {
        scan.push((k >> 32) as u32);
        tof.push((k & 0xFFFF_FFFF) as u32);
        intensity.push(sum.min(u32::MAX as u64) as u32);
    }
    let frame = FlatFrame {
        frame_id: 0,
        num_scans: num_local_scans,
        scan,
        tof,
        intensity,
    };

    let mut keep = filter_iterated(&frame, params);
    if let Some(hp) = halo {
        // Halo on the vertical survivors only (mirrors the MS1 path).
        let idx: Vec<usize> = (0..frame.len()).filter(|&i| keep[i]).collect();
        if !idx.is_empty() {
            let s: Vec<u32> = idx.iter().map(|&i| frame.scan[i]).collect();
            let t: Vec<u32> = idx.iter().map(|&i| frame.tof[i]).collect();
            let it: Vec<u32> = idx.iter().map(|&i| frame.intensity[i]).collect();
            let hmask = horizontal_halo_keep_mask(&s, &t, &it, num_local_scans, hp);
            for (k, &i) in idx.iter().enumerate() {
                if !hmask[k] {
                    keep[i] = false;
                }
            }
        }
    }

    for (i, &k) in keep.iter().enumerate() {
        if k {
            out.insert(key(frame.scan[i] + scan_begin, frame.tof[i]));
        }
    }
    out
}

/// Per-precursor keep sets plus the per-frame isolation windows, used to prune
/// MS/MS frames at write time.
pub(crate) struct MsmsKeep {
    /// Kept `(scan, tof)` keys per precursor id (indexed by id; 0 unused).
    keep: Vec<HashSet<u64>>,
    /// Frame `Id` -> isolation windows `(scan_begin, scan_end, precursor)`, sorted by scan_begin.
    frame_windows: HashMap<usize, Vec<(u32, u32, u32)>>,
}

impl MsmsKeep {
    /// Keep mask for one MS/MS frame: a point is kept iff its `(scan, tof)` is in
    /// its precursor's keep set. Points outside every isolation window are kept.
    pub(crate) fn keep_mask(&self, frame: &FlatFrame, frame_id: usize) -> Vec<bool> {
        let windows = match self.frame_windows.get(&frame_id) {
            Some(w) => w,
            None => return vec![true; frame.len()],
        };
        let mut keep = vec![true; frame.len()];
        for (i, slot) in keep.iter_mut().enumerate() {
            let s = frame.scan[i];
            if let Some(&(_, _, prec)) = windows.iter().find(|&&(sb, se, _)| s >= sb && s < se) {
                *slot = self.keep[prec as usize].contains(&key(s, frame.tof[i]));
            }
        }
        keep
    }
}

/// Build per-precursor keep sets: read each ddaPASEF MS/MS frame once, partition
/// its points to precursors, combine across frames, and filter each precursor's
/// combined spectrum.
pub(crate) fn build_msms_keep(
    reader: &FrameReader,
    meta: &[FrameMeta],
    windows: &[PasefWindow],
    params: &MsmsFilterParams,
    halo: Option<&HaloParams>,
) -> Result<MsmsKeep> {
    let max_prec = windows.iter().map(|w| w.precursor).max().unwrap_or(0) as usize;
    let mut sb0 = vec![u32::MAX; max_prec + 1];
    let mut se0 = vec![0u32; max_prec + 1];
    let mut frame_windows: HashMap<usize, Vec<(u32, u32, u32)>> = HashMap::new();
    for w in windows {
        let p = w.precursor as usize;
        sb0[p] = sb0[p].min(w.scan_begin);
        se0[p] = se0[p].max(w.scan_end);
        frame_windows
            .entry(w.frame)
            .or_default()
            .push((w.scan_begin, w.scan_end, w.precursor));
    }
    for v in frame_windows.values_mut() {
        v.sort_by_key(|&(sb, _, _)| sb);
    }

    // Accumulate a precursor's fragment points across all its frames.
    let mut raw: Vec<Vec<(u32, u32, u32)>> = vec![Vec::new(); max_prec + 1];
    for (i, m) in meta.iter().enumerate() {
        if m.num_peaks == 0 {
            continue;
        }
        let ws = match frame_windows.get(&m.id) {
            Some(w) => w,
            None => continue,
        };
        let frame = reader.get(i).map_err(|e| DnoiseError::FrameRead {
            index: i,
            message: e.to_string(),
        })?;
        let flat = FlatFrame::from_frame(&frame);
        for k in 0..flat.len() {
            let s = flat.scan[k];
            if let Some(&(_, _, prec)) = ws.iter().find(|&&(sb, se, _)| s >= sb && s < se) {
                raw[prec as usize].push((s, flat.tof[k], flat.intensity[k]));
            }
        }
    }

    let fp = params.as_filter_params();
    let keep: Vec<HashSet<u64>> = raw
        .into_par_iter()
        .enumerate()
        .map(|(p, pts)| {
            if p == 0 || pts.is_empty() {
                HashSet::new()
            } else {
                let num_local = se0[p].saturating_sub(sb0[p]) as usize;
                combine_and_filter(&pts, sb0[p], num_local, &fp, halo)
            }
        })
        .collect();

    Ok(MsmsKeep {
        keep,
        frame_windows,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vertical_params() -> FilterParams {
        FilterParams {
            mz_half_width: 0,
            min_feature_length: 3,
            max_internal_gap: 0,
            min_window_intensity: 0,
            min_feature_intensity: 0,
            num_iterations: 1,
        }
    }

    #[test]
    fn combine_and_filter_empty_input_is_empty() {
        assert!(combine_and_filter(&[], 0, 10, &vertical_params(), None).is_empty());
    }

    #[test]
    fn combine_and_filter_zero_local_scans_is_empty() {
        let pts = [(5, 100, 10)];
        assert!(combine_and_filter(&pts, 5, 0, &vertical_params(), None).is_empty());
    }

    #[test]
    fn combine_and_filter_keeps_a_vertical_streak_in_absolute_coords() {
        // One TOF column occupied across 10 consecutive scans starting at absolute
        // scan 10 -> a run of length 10 >= min_feature_length, so it survives.
        let pts: Vec<(u32, u32, u32)> = (10..20).map(|s| (s, 500, 100)).collect();
        let out = combine_and_filter(&pts, 10, 10, &vertical_params(), None);
        assert!(!out.is_empty());
        // Keys are packed in absolute scan coordinates.
        assert!(out.contains(&(((10u64) << 32) | 500)));
        assert!(out.contains(&(((19u64) << 32) | 500)));
    }

    #[test]
    fn combine_and_filter_drops_an_isolated_point() {
        // A single occupied scan is shorter than min_feature_length -> dropped.
        let pts = [(12, 500, 100)];
        assert!(combine_and_filter(&pts, 10, 10, &vertical_params(), None).is_empty());
    }

    #[test]
    fn keep_mask_passes_all_points_when_frame_has_no_windows() {
        let mk = MsmsKeep {
            keep: vec![HashSet::new()],
            frame_windows: HashMap::new(),
        };
        let frame = FlatFrame {
            frame_id: 3,
            num_scans: 5,
            scan: vec![0, 1, 2],
            tof: vec![10, 20, 30],
            intensity: vec![1, 1, 1],
        };
        assert_eq!(mk.keep_mask(&frame, 3), vec![true, true, true]);
    }

    #[test]
    fn keep_mask_prunes_in_window_points_by_precursor_keep_set() {
        // Frame 3 has one isolation window [scan 0, 2) belonging to precursor 1,
        // whose keep set contains only (scan 0, tof 10).
        let mut frame_windows = HashMap::new();
        frame_windows.insert(3usize, vec![(0u32, 2u32, 1u32)]);
        let mut prec1 = HashSet::new();
        prec1.insert(key(0, 10));
        let mk = MsmsKeep {
            keep: vec![HashSet::new(), prec1],
            frame_windows,
        };
        let frame = FlatFrame {
            frame_id: 3,
            num_scans: 5,
            scan: vec![0, 1, 3],
            tof: vec![10, 20, 30],
            intensity: vec![1, 1, 1],
        };
        // (0,10) in-window & kept; (1,20) in-window & not kept; (3,30) outside window.
        assert_eq!(mk.keep_mask(&frame, 3), vec![true, false, true]);
    }
}
