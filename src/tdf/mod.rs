//! timsTOF `.d` I/O: type-2 frame codec and `analysis.tdf` SQLite fixups.

pub mod encode;

use anyhow::{Context, Result};
use rusqlite::{Connection, OptionalExtension};
use std::path::Path;

/// Byte length of the leading header in `analysis.tdf_bin` that precedes the
/// first frame. Bruker reserves a (typically 64-byte, sometimes empty) block at
/// the start of the file; it equals the smallest `Frames.TimsId`. We copy it
/// verbatim and shift all rewritten offsets past it so the layout matches Bruker.
pub fn binary_header_len(tdf_path: &Path) -> Result<u64> {
    let conn = Connection::open(tdf_path)?;
    let min: Option<i64> = conn
        .query_row("SELECT MIN(TimsId) FROM Frames", [], |r| r.get(0))
        .optional()?;
    Ok(min.unwrap_or(0).max(0) as u64)
}

/// Minimal per-frame metadata read from the `Frames` table, ordered by `Id`
/// (which matches timsrust's frame index). Used to detect empty frames, which
/// timsrust cannot read.
pub struct FrameMeta {
    pub id: usize,
    pub num_scans: usize,
    pub num_peaks: u64,
    /// Bruker `MsMsType`: 0 = MS1, non-zero = MS/MS (8 = ddaPASEF, 9 = diaPASEF).
    pub ms_ms_type: i64,
}

impl FrameMeta {
    /// True for an MS1 frame (`MsMsType == 0`).
    pub fn is_ms1(&self) -> bool {
        self.ms_ms_type == 0
    }
}

/// Read `(Id, NumScans, NumPeaks, MsMsType)` for every frame, ordered by `Id`.
pub fn read_frame_meta(tdf_path: &Path) -> Result<Vec<FrameMeta>> {
    let conn =
        Connection::open(tdf_path).with_context(|| format!("open {}", tdf_path.display()))?;
    let mut stmt =
        conn.prepare("SELECT Id, NumScans, NumPeaks, MsMsType FROM Frames ORDER BY Id")?;
    let rows = stmt.query_map([], |r| {
        Ok(FrameMeta {
            id: r.get::<_, i64>(0)? as usize,
            num_scans: r.get::<_, i64>(1)? as usize,
            num_peaks: r.get::<_, i64>(2)? as u64,
            ms_ms_type: r.get::<_, i64>(3)?,
        })
    })?;
    Ok(rows.collect::<Result<_, _>>()?)
}

/// Per-frame values written back to the `Frames` table after filtering.
pub struct FrameUpdate {
    /// `Frames.Id` of this frame.
    pub frame_id: usize,
    /// New byte offset of the frame's record in the rewritten `analysis.tdf_bin`.
    pub tims_id: u64,
    /// Surviving peak count.
    pub num_peaks: u64,
    /// Max surviving intensity.
    pub max_intensity: u32,
    /// Sum of surviving intensities.
    pub summed_intensities: u64,
}

/// Apply all `Frames` updates and set `TimsCompressionType = 2` in one transaction.
pub fn update_metadata(tdf_path: &Path, updates: &[FrameUpdate]) -> Result<()> {
    let mut conn =
        Connection::open(tdf_path).with_context(|| format!("open {}", tdf_path.display()))?;
    let tx = conn.transaction()?;
    {
        let mut stmt = tx.prepare(
            "UPDATE Frames SET TimsId=?1, NumPeaks=?2, MaxIntensity=?3, SummedIntensities=?4 \
             WHERE Id=?5",
        )?;
        for u in updates {
            stmt.execute(rusqlite::params![
                u.tims_id as i64,
                u.num_peaks as i64,
                u.max_intensity as i64,
                u.summed_intensities as i64,
                u.frame_id as i64,
            ])?;
        }
    }
    tx.execute(
        "UPDATE GlobalMetadata SET Value='2' WHERE Key='TimsCompressionType'",
        [],
    )?;
    tx.commit()?;
    Ok(())
}
