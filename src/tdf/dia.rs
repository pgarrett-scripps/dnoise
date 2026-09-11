//! Checked DIA filter regions. Static windows stay separate; consecutive
//! one-scan, overlapping quadrupole steps form a monotonic scanning region.
use crate::{DnoiseError, Result};
use rusqlite::{Connection, OpenFlags};
use std::{collections::HashMap, path::Path, sync::Arc};

pub(crate) struct DiaRegion {
    pub begin: u32,
    pub end: u32,
    // Complete scan-dependent geometry, used for exact neighbor matching.
    pub signature: Arc<[u64]>,
    pub scan_varying: bool,
}
#[derive(Default)]
pub(crate) struct DiaRegions {
    pub groups: HashMap<i64, Vec<DiaRegion>>,
    pub frames: HashMap<usize, i64>,
    intervals: HashMap<i64, Vec<(u32, u32)>>,
}
impl DiaRegions {
    pub fn intervals(&self, frame: usize) -> Option<&[(u32, u32)]> {
        self.frames
            .get(&frame)
            .and_then(|g| self.intervals.get(g))
            .map(Vec::as_slice)
    }
    pub fn has_scanning(&self) -> bool {
        self.groups.values().flatten().any(|r| r.scan_varying)
    }
}
struct Window {
    begin: u32,
    end: u32,
    mz: f64,
    width: f64,
    energy: f64,
}
impl Window {
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
fn invalid(s: impl Into<String>) -> DnoiseError {
    DnoiseError::InvalidInput(format!("invalid DIA geometry: {}", s.into()))
}

pub(crate) fn read(path: &Path) -> Result<DiaRegions> {
    let db = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY)?;
    let mut grouped = HashMap::<i64, Vec<Window>>::new();
    let mut stmt=db.prepare("SELECT WindowGroup,ScanNumBegin,ScanNumEnd,IsolationMz,IsolationWidth,CollisionEnergy FROM DiaFrameMsMsWindows ORDER BY WindowGroup,ScanNumBegin")?;
    let mut rows = stmt.query([])?;
    while let Some(row) = rows.next()? {
        let group: i64 = row.get(0)?;
        let begin: i64 = row.get(1)?;
        let end: i64 = row.get(2)?;
        let mz: f64 = row.get(3)?;
        let width: f64 = row.get(4)?;
        let energy: f64 = row.get(5)?;
        let v = grouped.entry(group).or_default();
        if begin < 0
            || end <= begin
            || end > u32::MAX as i64
            || v.last().is_some_and(|w| begin < (w.end as i64))
        {
            return Err(invalid(format!(
                "invalid or overlapping interval in window group {group}"
            )));
        }
        if !mz.is_finite()
            || mz <= 0.0
            || !width.is_finite()
            || width <= 0.0
            || !energy.is_finite()
            || energy < 0.0
        {
            return Err(invalid(format!("nonphysical window in group {group}")));
        }
        v.push(Window {
            begin: begin as u32,
            end: end as u32,
            mz,
            width,
            energy,
        });
    }
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
                .flat_map(Window::signature)
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
    let mut stmt=db.prepare("SELECT p.Frame,p.WindowGroup,f.MsMsType,f.NumScans FROM DiaFrameMsMsInfo p LEFT JOIN Frames f ON f.Id=p.Frame ORDER BY p.Frame")?;
    let mut rows = stmt.query([])?;
    while let Some(row) = rows.next()? {
        let frame: usize = row.get(0)?;
        let group: i64 = row.get(1)?;
        let kind: Option<i64> = row.get(2)?;
        let scans: Option<i64> = row.get(3)?;
        let regions = result
            .groups
            .get(&group)
            .ok_or_else(|| invalid(format!("frame {frame} references missing group {group}")))?;
        if frame == 0
            || kind != Some(9)
            || scans.is_none_or(|s| s < 0 || regions.iter().any(|r| r.end as i64 > s))
            || result.frames.insert(frame, group).is_some()
        {
            return Err(invalid(format!("invalid frame/group reference at {frame}")));
        }
    }
    let missing:bool=db.query_row("SELECT EXISTS(SELECT 1 FROM Frames f WHERE f.MsMsType=9 AND f.NumPeaks>0 AND NOT EXISTS(SELECT 1 FROM DiaFrameMsMsInfo p WHERE p.Frame=f.Id))",[],|r| r.get(0))?;
    if missing {
        return Err(invalid(
            "nonempty DIA frame lacks an isolation-window group",
        ));
    }
    Ok(result)
}
