//! Acquisition classification and checked prm-PASEF metadata.
use super::{FrameMeta, PrmWindows};
use crate::{Acquisition, DnoiseError, Result};
use rusqlite::{Connection, OpenFlags};
use std::collections::HashSet;
use std::path::Path;

fn invalid(message: impl Into<String>) -> DnoiseError {
    DnoiseError::InvalidInput(format!("invalid prm-PASEF metadata: {}", message.into()))
}

pub(crate) struct AcquisitionInfo {
    pub kind: Acquisition,
    pub prm_windows: PrmWindows,
}

pub(crate) fn detect(path: &Path, frames: &[FrameMeta]) -> Result<Acquisition> {
    Ok(inspect(path, frames)?.kind)
}

pub(crate) fn inspect(path: &Path, frames: &[FrameMeta]) -> Result<AcquisitionInfo> {
    let mut prm_windows = PrmWindows::new();
    let db = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY)?;
    let has_table = |name: &str| -> Result<bool> {
        Ok(db.query_row(
            "SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type='table' AND name=?1)",
            [name],
            |r| r.get(0),
        )?)
    };
    let has_prm = frames.iter().any(|f| f.ms_ms_type == 10);
    let has_events = has_table("PrmFrameMsMsInfo")?;
    let populated = has_events
        && db.query_row("SELECT EXISTS(SELECT 1 FROM PrmFrameMsMsInfo)", [], |r| {
            r.get::<_, bool>(0)
        })?;
    // Empty PRM tables occur in non-PRM files and must not activate PRM handling.
    if has_prm || populated {
        if !has_events || !has_table("PrmTargets")? {
            return Err(invalid("PrmFrameMsMsInfo and PrmTargets are required"));
        }
        let mut covered = HashSet::new();
        let mut ends = std::collections::HashMap::<i64, i64>::new();
        let mut stmt = db.prepare(
            "SELECT p.Frame, p.ScanNumBegin, p.ScanNumEnd, p.IsolationMz,
                    p.IsolationWidth, p.CollisionEnergy, p.Target,
                    f.MsMsType, f.NumScans, t.Id
             FROM PrmFrameMsMsInfo p
             LEFT JOIN Frames f ON f.Id=p.Frame
             LEFT JOIN PrmTargets t ON t.Id=p.Target
             ORDER BY p.Frame, p.ScanNumBegin",
        )?;
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            let frame: i64 = row.get(0)?;
            let begin: i64 = row.get(1)?;
            let end: i64 = row.get(2)?;
            let mz: f64 = row.get(3)?;
            let width: f64 = row.get(4)?;
            let energy: f64 = row.get(5)?;
            let target: i64 = row.get(6)?;
            let kind: Option<i64> = row.get(7)?;
            let scans: Option<i64> = row.get(8)?;
            let target_id: Option<i64> = row.get(9)?;
            if frame <= 0 || kind != Some(10) || target_id != Some(target) {
                return Err(invalid(format!(
                    "event at frame {frame} must reference a PRM frame (MsMsType 10) and an existing target"
                )));
            }
            if begin < 0
                || end <= begin
                || scans.is_none_or(|s| end > s || s > u32::MAX as i64)
                || ends.get(&frame).is_some_and(|&last| begin < last)
            {
                return Err(invalid(format!(
                    "invalid or overlapping scan interval [{begin}, {end}) at frame {frame}"
                )));
            }
            if !mz.is_finite()
                || mz <= 0.0
                || !width.is_finite()
                || width <= 0.0
                || !energy.is_finite()
                || energy < 0.0
            {
                return Err(invalid(format!(
                    "invalid isolation geometry at frame {frame}"
                )));
            }
            // Touching events remain distinct. Never merge target boundaries.
            ends.insert(frame, end);
            covered.insert(frame as usize);
            prm_windows
                .entry(frame as usize)
                .or_default()
                .push((begin as u32, end as u32));
        }
        for frame in frames {
            if frame.ms_ms_type == 10 && frame.num_peaks > 0 && !covered.contains(&frame.id) {
                return Err(invalid(format!(
                    "nonempty PRM frame {} has no isolation events",
                    frame.id
                )));
            }
        }
        for table in ["PasefFrameMsMsInfo", "DiaFrameMsMsInfo"] {
            if has_table(table)? {
                let conflict: bool = db.query_row(
                    &format!("SELECT EXISTS(SELECT 1 FROM {table} p JOIN Frames f ON f.Id=p.Frame WHERE f.MsMsType=10)"),
                    [], |r| r.get(0),
                )?;
                if conflict {
                    return Err(invalid(format!(
                        "PRM frames also have isolation events in {table}"
                    )));
                }
            }
        }
        if has_table("PrmFrameMeasurementMode")? {
            let invalid_frame: bool = db.query_row(
                "SELECT EXISTS(SELECT 1 FROM PrmFrameMeasurementMode p
                 LEFT JOIN Frames f ON f.Id=p.Frame WHERE f.Id IS NULL OR f.MsMsType!=10)",
                [],
                |r| r.get(0),
            )?;
            if invalid_frame {
                return Err(invalid(
                    "PrmFrameMeasurementMode references a missing or non-PRM frame",
                ));
            }
        }
    }
    let kinds: HashSet<_> = frames
        .iter()
        .filter(|f| !f.is_ms1())
        .map(|f| f.ms_ms_type)
        .collect();
    let kind = if kinds.len() > 1 {
        Acquisition::Mixed
    } else {
        match kinds.iter().next() {
            None => Acquisition::Ms1Only,
            Some(8) => Acquisition::DdaPasef,
            Some(9) => Acquisition::DiaPasef,
            Some(10) => Acquisition::PrmPasef,
            _ => Acquisition::Unknown,
        }
    };
    if kind == Acquisition::DiaPasef {
        // Validate references before timsrust sees them: orphan DIA frame IDs
        // otherwise reach unchecked indexing in the vendor-independent reader.
        super::dia::read(path)?;
    }
    Ok(AcquisitionInfo { kind, prm_windows })
}
