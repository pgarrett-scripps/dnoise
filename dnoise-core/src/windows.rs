//! Plain per-run acquisition metadata the stages consume: frame metadata,
//! PASEF isolation events, diaPASEF window schemes and prm-PASEF event windows.
//!
//! These are the in-memory shapes of rows from `analysis.tdf` (`Frames`,
//! `PasefFrameMsMsInfo`, `DiaFrameMsMsInfo` / `DiaFrameMsMsWindows`,
//! `PrmFrameMsMsInfo`). The `dnoise` crate reads and checks them with SQLite;
//! an embedding reader can fill them from its own metadata.

use std::collections::HashMap;
use std::sync::Arc;

/// Minimal per-frame metadata from the `Frames` table, ordered by `Id` (the
/// position in this list is the 0-based frame index every stage uses).
#[derive(Debug, Clone, PartialEq)]
pub struct FrameMeta {
    /// `Frames.Id`.
    pub id: usize,
    /// `Frames.NumScans`.
    pub num_scans: usize,
    /// `Frames.NumPeaks`. Frames with zero peaks are never decoded.
    pub num_peaks: u64,
    /// Bruker `MsMsType`: 0 = MS1; 8 = ddaPASEF, 9 = diaPASEF, 10 = prm-PASEF.
    pub ms_ms_type: i64,
    /// Retention time in **seconds** (`Frames.Time`). Used by the RT crop and
    /// the neighbor retention-time limit.
    pub rt: f64,
}

impl FrameMeta {
    /// True for an MS1 frame (`MsMsType == 0`).
    pub fn is_ms1(&self) -> bool {
        self.ms_ms_type == 0
    }
}

/// One ddaPASEF MS/MS isolation event: in `frame`, scans `[scan_begin, scan_end)`
/// were isolated and fragmented for `precursor` (a `PasefFrameMsMsInfo` row).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PasefWindow {
    /// Frame `Id` (matches `Frames.Id`).
    pub frame: usize,
    /// First scan of the isolation window (inclusive).
    pub scan_begin: u32,
    /// Last scan of the isolation window (exclusive).
    pub scan_end: u32,
    /// Precursor `Id`.
    pub precursor: u32,
}

/// diaPASEF isolation-window scheme: every MS/MS frame belongs to one
/// `WindowGroup`, and each group defines a set of mobility-scan intervals
/// `[ScanNumBegin, ScanNumEnd)` over which the quadrupole isolated a precursor
/// m/z band. Signal outside every interval was never isolated, and signal in two
/// different intervals comes from unrelated isolation events — so the per-window
/// MS/MS filter ([`crate::dia_window`]) uses these intervals both to drop
/// out-of-window points and to filter each window independently (no cross-talk).
///
/// [`DiaWindows::from_pasef`] builds the same per-frame interval map from
/// ddaPASEF isolation events, so the out-of-window gate serves both acquisitions.
#[derive(Debug, Default, Clone)]
pub struct DiaWindows {
    /// `Frames.Id` -> sorted, non-overlapping `[scan_begin, scan_end)` intervals.
    frame_intervals: HashMap<usize, Vec<(u32, u32)>>,
}

impl DiaWindows {
    /// Wrap a ready `Frames.Id -> intervals` map. Each list must already be
    /// sorted by `begin` and non-overlapping (merge touching windows first), as
    /// [`crate::dia_window::in_window_mask`] requires.
    pub fn from_frame_intervals(frame_intervals: HashMap<usize, Vec<(u32, u32)>>) -> Self {
        Self { frame_intervals }
    }

    /// True when no window scheme was found (e.g. ddaPASEF data, where the
    /// `DiaFrameMsMs*` tables are absent or empty).
    pub fn is_empty(&self) -> bool {
        self.frame_intervals.is_empty()
    }

    /// Sorted, non-overlapping isolation-window scan intervals for one MS/MS
    /// frame, or `None` if the frame has no window-group entry.
    pub fn intervals(&self, frame_id: usize) -> Option<&[(u32, u32)]> {
        self.frame_intervals.get(&frame_id).map(Vec::as_slice)
    }

    /// The same per-frame interval map built from ddaPASEF isolation events
    /// ([`PasefWindow`]) instead of a diaPASEF window scheme, so the same
    /// out-of-window gate serves both acquisitions. Events are sorted per frame
    /// and merged where they touch or overlap, matching the contract
    /// [`crate::dia_window::in_window_mask`] requires. Empty input (diaPASEF,
    /// non-PASEF) yields an empty map, which callers treat as "no gate".
    pub fn from_pasef(windows: &[PasefWindow]) -> Self {
        let mut frame_intervals: HashMap<usize, Vec<(u32, u32)>> = HashMap::new();
        for w in windows {
            if w.scan_end > w.scan_begin {
                frame_intervals
                    .entry(w.frame)
                    .or_default()
                    .push((w.scan_begin, w.scan_end));
            }
        }
        for v in frame_intervals.values_mut() {
            v.sort_unstable();
            let mut merged: Vec<(u32, u32)> = Vec::with_capacity(v.len());
            for &(sb, se) in v.iter() {
                match merged.last_mut() {
                    Some(last) if sb <= last.1 => last.1 = last.1.max(se),
                    _ => merged.push((sb, se)),
                }
            }
            *v = merged;
        }
        DiaWindows { frame_intervals }
    }
}

/// One diaPASEF isolation window as a 2-D precursor-space box: an m/z band over a
/// mobility-scan interval. Used by the MS1 out-of-window gate ([`crate::dia_ms1`]),
/// which (unlike the MS/MS scan gate) needs the m/z extent too.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiaMs1Box {
    /// First mobility scan of the window (inclusive).
    pub scan_begin: u32,
    /// Last mobility scan of the window (exclusive, as stored).
    pub scan_end: u32,
    /// Low m/z edge (`IsolationMz - IsolationWidth/2`).
    pub mz_lo: f64,
    /// High m/z edge (`IsolationMz + IsolationWidth/2`).
    pub mz_hi: f64,
}

/// Checked prm-PASEF isolation events: `Frames.Id` -> unmerged
/// `[scan_begin, scan_end)` event intervals, sorted by begin.
pub type PrmWindows = HashMap<usize, Vec<(u32, u32)>>;

/// One `DiaFrameMsMsWindows` row, used to build [`DiaRegions`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiaIsolationWindow {
    /// First scan (inclusive).
    pub begin: u32,
    /// Last scan (exclusive).
    pub end: u32,
    /// `IsolationMz`.
    pub mz: f64,
    /// `IsolationWidth`.
    pub width: f64,
    /// `CollisionEnergy`.
    pub energy: f64,
}

impl DiaIsolationWindow {
    fn signature(&self) -> [u64; 5] {
        [
            self.begin as u64,
            self.end as u64,
            self.mz.to_bits(),
            self.width.to_bits(),
            self.energy.to_bits(),
        ]
    }

    fn joins(&self, next: &Self, direction: f64) -> bool {
        let delta = next.mz - self.mz;
        self.end - self.begin == 1
            && next.end - next.begin == 1
            && self.end == next.begin
            && (self.width - next.width).abs() <= 1e-9 * self.width.max(1.0)
            && delta != 0.0
            && (direction == 0.0 || delta.signum() == direction)
            && delta.abs() < (self.width + next.width) * 0.5
    }
}

/// One DIA filter region: a static window, or a run of consecutive one-scan,
/// overlapping quadrupole steps joined into one monotonic scanning region.
#[derive(Debug, Clone)]
pub struct DiaRegion {
    /// First scan (inclusive).
    pub begin: u32,
    /// Last scan (exclusive).
    pub end: u32,
    /// Complete scan-dependent geometry, used for exact neighbor matching.
    pub signature: Arc<[u64]>,
    /// True when the region joins several scan-dependent steps.
    pub scan_varying: bool,
}

/// Checked DIA filter regions per window group, plus the frame -> group map.
/// Static windows stay separate; consecutive one-scan, overlapping quadrupole
/// steps form a monotonic scanning region.
#[derive(Debug, Default, Clone)]
pub struct DiaRegions {
    /// `WindowGroup` -> regions, sorted by scan.
    pub groups: HashMap<i64, Vec<DiaRegion>>,
    /// `Frames.Id` -> `WindowGroup`.
    pub frames: HashMap<usize, i64>,
    intervals: HashMap<i64, Vec<(u32, u32)>>,
}

impl DiaRegions {
    /// Build the regions from each group's windows (sorted by `begin`,
    /// non-overlapping, physically valid: the caller checks the rows). Fill
    /// [`Self::frames`] afterwards.
    pub fn from_groups(grouped: HashMap<i64, Vec<DiaIsolationWindow>>) -> Self {
        let mut result = DiaRegions::default();
        for (group, windows) in grouped {
            let mut regions = Vec::new();
            let mut start = 0;
            while start < windows.len() {
                let mut end = start + 1;
                let mut direction = 0.0;
                while end < windows.len() && windows[end - 1].joins(&windows[end], direction) {
                    direction = (windows[end].mz - windows[end - 1].mz).signum();
                    end += 1;
                }
                let signature: Vec<_> = windows[start..end]
                    .iter()
                    .flat_map(DiaIsolationWindow::signature)
                    .collect();
                regions.push(DiaRegion {
                    begin: windows[start].begin,
                    end: windows[end - 1].end,
                    signature: signature.into(),
                    scan_varying: end - start > 1,
                });
                start = end;
            }
            result
                .intervals
                .insert(group, regions.iter().map(|r| (r.begin, r.end)).collect());
            result.groups.insert(group, regions);
        }
        result
    }

    /// Region scan intervals for frame `Frames.Id`, or `None` when it has no group.
    pub fn intervals(&self, frame: usize) -> Option<&[(u32, u32)]> {
        self.frames
            .get(&frame)
            .and_then(|g| self.intervals.get(g))
            .map(Vec::as_slice)
    }

    /// True when any region joins scan-dependent steps.
    pub fn has_scanning(&self) -> bool {
        self.groups.values().flatten().any(|r| r.scan_varying)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn w(begin: u32, end: u32, mz: f64) -> DiaIsolationWindow {
        DiaIsolationWindow {
            begin,
            end,
            mz,
            width: 25.0,
            energy: 30.0,
        }
    }

    #[test]
    fn static_windows_stay_separate_and_steps_join() {
        let mut grouped = HashMap::new();
        grouped.insert(1, vec![w(0, 10, 400.0), w(20, 30, 500.0)]);
        grouped.insert(2, vec![w(0, 1, 400.0), w(1, 2, 401.0), w(2, 3, 402.0)]);
        let mut r = DiaRegions::from_groups(grouped);
        r.frames.insert(5, 1);
        r.frames.insert(6, 2);
        assert_eq!(r.intervals(5), Some(&[(0, 10), (20, 30)][..]));
        assert_eq!(r.intervals(6), Some(&[(0, 3)][..]));
        assert_eq!(r.intervals(7), None);
        assert!(r.has_scanning());
    }

    #[test]
    fn from_pasef_merges_touching_events() {
        let ev = |frame, scan_begin, scan_end| PasefWindow {
            frame,
            scan_begin,
            scan_end,
            precursor: 1,
        };
        let d = DiaWindows::from_pasef(&[ev(2, 30, 40), ev(2, 10, 20), ev(2, 20, 25), ev(3, 5, 5)]);
        assert_eq!(d.intervals(2), Some(&[(10, 25), (30, 40)][..]));
        assert_eq!(d.intervals(3), None);
    }
}
