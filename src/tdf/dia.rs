//! Checked DIA filter regions. Static windows stay separate; consecutive
//! one-scan, overlapping quadrupole steps form a monotonic scanning region.
use crate::{DnoiseError, Result};
use dnoise_core::windows::{DiaIsolationWindow as Window, DiaRegions};
use rusqlite::{Connection, OpenFlags};
use std::{collections::HashMap, path::Path};

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
    let mut result = DiaRegions::from_groups(grouped);
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
