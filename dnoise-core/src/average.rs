//! Optional pre-filter smoothing: a centered moving average over the MS1-frame
//! subsequence, applied *before* the vertical-IM filter (and any downstream
//! centroider). This stage is completely decoupled from the filter — it just
//! produces ordinary [`FlatFrame`]s that the rest of the pipeline consumes
//! exactly as it consumes raw frames.
//!
//! For each MS1 frame we replace its points with the per-bin average of the
//! `2r+1` nearest MS1 frames `[F-r ..= F+r]` (clamped at the run ends). A real
//! ion recurs at the same `(scan, tof)` across the window so its average stays
//! near its native intensity, while noise lands on different bins each frame and
//! is divided down — raising signal-to-noise for the filter that follows.

use crate::frame::FlatFrame;
use std::collections::HashMap;

/// Centered moving average of the `window` frames into a new [`FlatFrame`].
///
/// Intensities are summed per exact `(scan, tof)` bin across the window and
/// divided by the number of frames present (`window.len()`), so the result
/// stays on the native single-frame intensity scale regardless of window size
/// (including the smaller windows at the start/end of the run). Bins that round
/// to zero are dropped. Point order is unspecified — the type-2 encoder
/// canonicalises by scan group and TOF, so output bytes are deterministic.
pub fn running_average(frame_id: usize, num_scans: usize, window: &[&FlatFrame]) -> FlatFrame {
    let n = window.len().max(1) as u64;

    // Sum intensities per (scan, tof) bin. Key packs scan into the high 32 bits.
    let mut acc: HashMap<u64, u64> = HashMap::new();
    for f in window {
        for i in 0..f.len() {
            let key = ((f.scan[i] as u64) << 32) | f.tof[i] as u64;
            *acc.entry(key).or_insert(0) += f.intensity[i] as u64;
        }
    }

    let mut scan = Vec::with_capacity(acc.len());
    let mut tof = Vec::with_capacity(acc.len());
    let mut intensity = Vec::with_capacity(acc.len());
    for (key, sum) in acc {
        // Round to nearest; an averaged bin is bounded by u32::MAX, but clamp defensively.
        let avg = ((sum + n / 2) / n).min(u32::MAX as u64) as u32;
        if avg == 0 {
            continue;
        }
        scan.push((key >> 32) as u32);
        tof.push((key & 0xFFFF_FFFF) as u32);
        intensity.push(avg);
    }

    FlatFrame {
        frame_id,
        num_scans,
        scan,
        tof,
        intensity,
    }
}
