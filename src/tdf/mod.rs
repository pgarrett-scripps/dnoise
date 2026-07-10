//! Internal `analysis.tdf` SQLite access and fixups (crate-private plumbing).

use crate::error::Result;
use rusqlite::{Connection, OptionalExtension};
use std::collections::HashMap;
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
    /// Retention time in **seconds** (`Frames.Time`). Used by the RT crop.
    pub rt: f64,
}

impl FrameMeta {
    /// True for an MS1 frame (`MsMsType == 0`).
    pub fn is_ms1(&self) -> bool {
        self.ms_ms_type == 0
    }
}

/// One ddaPASEF MS/MS isolation event: in `frame`, scans `[scan_begin, scan_end)`
/// were isolated and fragmented for `precursor`.
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

/// Read every `PasefFrameMsMsInfo` row (ddaPASEF). Empty for non-DDA data —
/// diaPASEF `.d` files omit the table entirely, which is treated as "no ddaPASEF
/// windows" (the caller then falls back to whole-frame MS/MS filtering).
pub fn read_pasef_msms(tdf_path: &Path) -> Result<Vec<PasefWindow>> {
    let conn = Connection::open(tdf_path)?;
    let has_table: bool = conn
        .query_row(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='PasefFrameMsMsInfo'",
            [],
            |_| Ok(true),
        )
        .optional()?
        .unwrap_or(false);
    if !has_table {
        return Ok(Vec::new());
    }
    let mut stmt = conn.prepare(
        "SELECT Frame, ScanNumBegin, ScanNumEnd, Precursor FROM PasefFrameMsMsInfo \
         ORDER BY Frame",
    )?;
    let rows = stmt.query_map([], |r| {
        Ok(PasefWindow {
            frame: r.get::<_, i64>(0)? as usize,
            scan_begin: r.get::<_, i64>(1)? as u32,
            scan_end: r.get::<_, i64>(2)? as u32,
            precursor: r.get::<_, i64>(3)? as u32,
        })
    })?;
    Ok(rows.collect::<std::result::Result<_, rusqlite::Error>>()?)
}

/// diaPASEF isolation-window scheme: every MS/MS frame belongs to one
/// `WindowGroup`, and each group defines a set of mobility-scan intervals
/// `[ScanNumBegin, ScanNumEnd)` over which the quadrupole isolated a precursor
/// m/z band. Signal outside every interval was never isolated, and signal in two
/// different intervals comes from unrelated isolation events — so the per-window
/// MS/MS filter ([`crate::dia_window`]) uses these intervals both to drop
/// out-of-window points and to filter each window independently (no cross-talk).
#[derive(Debug, Default)]
pub struct DiaWindows {
    /// `Frames.Id` -> sorted, non-overlapping `[scan_begin, scan_end)` intervals.
    frame_intervals: HashMap<usize, Vec<(u32, u32)>>,
}

impl DiaWindows {
    /// True when no diaPASEF window scheme was found (e.g. ddaPASEF data, where
    /// the `DiaFrameMsMs*` tables are absent or empty).
    pub fn is_empty(&self) -> bool {
        self.frame_intervals.is_empty()
    }

    /// Sorted, non-overlapping isolation-window scan intervals for one MS/MS
    /// frame, or `None` if the frame has no window-group entry.
    pub fn intervals(&self, frame_id: usize) -> Option<&[(u32, u32)]> {
        self.frame_intervals.get(&frame_id).map(Vec::as_slice)
    }
}

/// Read the diaPASEF isolation-window scheme by joining `DiaFrameMsMsInfo`
/// (`Frame` -> `WindowGroup`) with `DiaFrameMsMsWindows` (`WindowGroup` ->
/// `[ScanNumBegin, ScanNumEnd)`). Returns an empty [`DiaWindows`] when either
/// table is absent or empty (ddaPASEF, or non-PASEF data), which callers treat as
/// "no DIA windows" and fall back to whole-frame handling. Per-group intervals are
/// sorted by `ScanNumBegin` and merged so adjacent/overlapping windows coalesce.
pub fn read_dia_windows(tdf_path: &Path) -> Result<DiaWindows> {
    let conn = Connection::open(tdf_path)?;
    let has = |name: &str| -> Result<bool> {
        Ok(conn
            .query_row(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?1",
                [name],
                |_| Ok(true),
            )
            .optional()?
            .unwrap_or(false))
    };
    if !has("DiaFrameMsMsInfo")? || !has("DiaFrameMsMsWindows")? {
        return Ok(DiaWindows::default());
    }

    // WindowGroup -> sorted, merged [begin, end) intervals.
    let mut group_intervals: HashMap<i64, Vec<(u32, u32)>> = HashMap::new();
    {
        let mut stmt = conn.prepare(
            "SELECT WindowGroup, ScanNumBegin, ScanNumEnd FROM DiaFrameMsMsWindows \
             ORDER BY WindowGroup, ScanNumBegin",
        )?;
        let rows = stmt.query_map([], |r| {
            Ok((
                r.get::<_, i64>(0)?,
                r.get::<_, i64>(1)? as u32,
                r.get::<_, i64>(2)? as u32,
            ))
        })?;
        for row in rows {
            let (g, sb, se) = row?;
            if se <= sb {
                continue;
            }
            let v = group_intervals.entry(g).or_default();
            // Rows arrive sorted by ScanNumBegin, so we can merge against the last.
            match v.last_mut() {
                Some(last) if sb <= last.1 => last.1 = last.1.max(se),
                _ => v.push((sb, se)),
            }
        }
    }

    // Frame -> WindowGroup, resolved to the group's interval list.
    let mut frame_intervals: HashMap<usize, Vec<(u32, u32)>> = HashMap::new();
    {
        let mut stmt = conn.prepare("SELECT Frame, WindowGroup FROM DiaFrameMsMsInfo")?;
        let rows = stmt.query_map([], |r| {
            Ok((r.get::<_, i64>(0)? as usize, r.get::<_, i64>(1)?))
        })?;
        for row in rows {
            let (frame, group) = row?;
            if let Some(iv) = group_intervals.get(&group) {
                frame_intervals.insert(frame, iv.clone());
            }
        }
    }

    Ok(DiaWindows { frame_intervals })
}

/// One diaPASEF isolation window as a 2-D precursor-space box: an m/z band over a
/// mobility-scan interval. Used by the MS1 out-of-window gate ([`crate::dia_ms1`]),
/// which (unlike the MS/MS scan gate) needs the m/z extent too.
#[derive(Debug, Clone, Copy)]
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

/// Read the distinct diaPASEF isolation windows as 2-D `(scan, m/z)` boxes — the
/// union of the precursor-space tiling, independent of which frame samples each
/// window. Returns an empty vector when `DiaFrameMsMsWindows` is absent or carries
/// no m/z information (ddaPASEF or non-PASEF data), which the caller treats as "no
/// MS1 gate".
pub fn read_dia_ms1_boxes(tdf_path: &Path) -> Result<Vec<DiaMs1Box>> {
    let conn = Connection::open(tdf_path)?;
    let has_table: bool = conn
        .query_row(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='DiaFrameMsMsWindows'",
            [],
            |_| Ok(true),
        )
        .optional()?
        .unwrap_or(false);
    if !has_table {
        return Ok(Vec::new());
    }
    let mut stmt = conn.prepare(
        "SELECT DISTINCT ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth \
         FROM DiaFrameMsMsWindows \
         WHERE IsolationMz IS NOT NULL AND IsolationWidth IS NOT NULL",
    )?;
    let rows = stmt.query_map([], |r| {
        let scan_begin = r.get::<_, i64>(0)? as u32;
        let scan_end = r.get::<_, i64>(1)? as u32;
        let iso_mz = r.get::<_, f64>(2)?;
        let iso_width = r.get::<_, f64>(3)?;
        Ok(DiaMs1Box {
            scan_begin,
            scan_end,
            mz_lo: iso_mz - iso_width / 2.0,
            mz_hi: iso_mz + iso_width / 2.0,
        })
    })?;
    Ok(rows
        .collect::<std::result::Result<Vec<_>, rusqlite::Error>>()?
        .into_iter()
        .filter(|b| b.scan_end > b.scan_begin)
        .collect())
}

/// Read the ddaPASEF/PASEF MS1 selection polygon (the "IMS PolygonFilter") as
/// parallel `(m/z, 1/K0)` vertex arrays. Bruker stores the two arrays as
/// little-endian `f64` BLOBs in `GroupProperties`, keyed by the
/// `IMS_PolygonFilter_Mass` / `IMS_PolygonFilter_Mobility` property definitions
/// (resolved by permanent name, as the numeric IDs are not guaranteed stable).
///
/// Returns `None` when the run carries no usable polygon (property absent, no
/// stored value, fewer than 3 vertices, or mismatched array lengths), which the
/// caller treats as "no MS1 polygon gate".
pub fn read_selection_polygon(tdf_path: &Path) -> Result<Option<(Vec<f64>, Vec<f64>)>> {
    let conn = Connection::open(tdf_path)?;
    let has_tables: bool = conn
        .query_row(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='PropertyDefinitions'",
            [],
            |_| Ok(true),
        )
        .optional()?
        .unwrap_or(false);
    if !has_tables {
        return Ok(None);
    }

    let prop_id = |name: &str| -> Result<Option<i64>> {
        Ok(conn
            .query_row(
                "SELECT Id FROM PropertyDefinitions WHERE PermanentName=?1",
                [name],
                |r| r.get::<_, i64>(0),
            )
            .optional()?)
    };
    let (Some(mz_id), Some(im_id)) = (
        prop_id("IMS_PolygonFilter_Mass")?,
        prop_id("IMS_PolygonFilter_Mobility")?,
    ) else {
        return Ok(None);
    };

    // One polygon per run: take the first group that stores each array.
    let read_blob = |prop: i64| -> Result<Option<Vec<f64>>> {
        let blob: Option<Vec<u8>> = conn
            .query_row(
                "SELECT Value FROM GroupProperties WHERE Property=?1 LIMIT 1",
                [prop],
                |r| r.get::<_, Vec<u8>>(0),
            )
            .optional()?;
        Ok(blob.map(|b| {
            b.chunks_exact(8)
                .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
                .collect()
        }))
    };
    let (Some(mz), Some(im)) = (read_blob(mz_id)?, read_blob(im_id)?) else {
        return Ok(None);
    };
    if mz.len() < 3 || mz.len() != im.len() {
        return Ok(None);
    }
    Ok(Some((mz, im)))
}

/// Read `(Id, NumScans, NumPeaks, MsMsType, Time)` for every frame, ordered by `Id`.
pub fn read_frame_meta(tdf_path: &Path) -> Result<Vec<FrameMeta>> {
    let conn = Connection::open(tdf_path)?;
    let mut stmt =
        conn.prepare("SELECT Id, NumScans, NumPeaks, MsMsType, Time FROM Frames ORDER BY Id")?;
    let rows = stmt.query_map([], |r| {
        Ok(FrameMeta {
            id: r.get::<_, i64>(0)? as usize,
            num_scans: r.get::<_, i64>(1)? as usize,
            num_peaks: r.get::<_, i64>(2)? as u64,
            ms_ms_type: r.get::<_, i64>(3)?,
            rt: r.get::<_, f64>(4)?,
        })
    })?;
    Ok(rows.collect::<std::result::Result<_, rusqlite::Error>>()?)
}

/// Read the acquired m/z range `(lower, upper)` from `GlobalMetadata`
/// (`MzAcqRangeLower` / `MzAcqRangeUpper`). Used to pick a reference m/z for the
/// ppm-to-TOF-index conversion. Returns `None` if either key is absent or
/// unparseable.
pub fn read_mz_acq_range(tdf_path: &Path) -> Result<Option<(f64, f64)>> {
    let conn = Connection::open(tdf_path)?;
    let val = |key: &str| -> Result<Option<f64>> {
        Ok(conn
            .query_row(
                "SELECT Value FROM GlobalMetadata WHERE Key=?1",
                [key],
                |r| r.get::<_, String>(0),
            )
            .optional()?
            .and_then(|s| s.trim().parse::<f64>().ok()))
    };
    match (val("MzAcqRangeLower")?, val("MzAcqRangeUpper")?) {
        (Some(lo), Some(hi)) if hi > lo => Ok(Some((lo, hi))),
        _ => Ok(None),
    }
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
    let mut conn = Connection::open(tdf_path)?;
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a throwaway `.tdf` with the two diaPASEF window tables populated,
    /// returning its path (caller removes it).
    fn temp_dia_tdf(windows: &[(i64, u32, u32)], info: &[(usize, i64)]) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!(
            "dnoise_dia_test_{}_{:?}.tdf",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let conn = Connection::open(&path).unwrap();
        conn.execute_batch(
            "CREATE TABLE DiaFrameMsMsWindows (WindowGroup INTEGER, ScanNumBegin INTEGER, \
             ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, CollisionEnergy REAL);\
             CREATE TABLE DiaFrameMsMsInfo (Frame INTEGER, WindowGroup INTEGER);",
        )
        .unwrap();
        for &(g, sb, se) in windows {
            conn.execute(
                "INSERT INTO DiaFrameMsMsWindows (WindowGroup, ScanNumBegin, ScanNumEnd) \
                 VALUES (?1, ?2, ?3)",
                rusqlite::params![g, sb as i64, se as i64],
            )
            .unwrap();
        }
        for &(f, g) in info {
            conn.execute(
                "INSERT INTO DiaFrameMsMsInfo (Frame, WindowGroup) VALUES (?1, ?2)",
                rusqlite::params![f as i64, g],
            )
            .unwrap();
        }
        path
    }

    #[test]
    fn read_dia_windows_absent_tables_is_empty() {
        let path =
            std::env::temp_dir().join(format!("dnoise_dia_empty_{}.tdf", std::process::id()));
        let conn = Connection::open(&path).unwrap();
        conn.execute_batch("CREATE TABLE Frames (Id INTEGER);")
            .unwrap();
        drop(conn);
        let w = read_dia_windows(&path).unwrap();
        assert!(w.is_empty());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn read_dia_windows_maps_frames_and_merges_intervals() {
        // Group 1: contiguous windows that must merge into one [100,300).
        // Group 2: a gap that must stay two intervals.
        let windows = [
            (1i64, 100u32, 200u32),
            (1, 200, 300),
            (2, 100, 150),
            (2, 250, 400),
        ];
        // Two frames in group 1, one in group 2, one frame with no group entry.
        let info = [(10usize, 1i64), (11, 1), (12, 2)];
        let path = temp_dia_tdf(&windows, &info);

        let w = read_dia_windows(&path).unwrap();
        assert!(!w.is_empty());
        assert_eq!(w.intervals(10), Some(&[(100u32, 300u32)][..]));
        assert_eq!(w.intervals(11), Some(&[(100u32, 300u32)][..]));
        assert_eq!(
            w.intervals(12),
            Some(&[(100u32, 150u32), (250u32, 400u32)][..])
        );
        assert_eq!(w.intervals(99), None); // frame absent from DiaFrameMsMsInfo

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn read_dia_ms1_boxes_absent_table_is_empty() {
        let path =
            std::env::temp_dir().join(format!("dnoise_dia_ms1_empty_{}.tdf", std::process::id()));
        let conn = Connection::open(&path).unwrap();
        conn.execute_batch("CREATE TABLE Frames (Id INTEGER);")
            .unwrap();
        drop(conn);
        assert!(read_dia_ms1_boxes(&path).unwrap().is_empty());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn read_dia_ms1_boxes_reads_mz_band_and_dedupes() {
        let path = std::env::temp_dir().join(format!(
            "dnoise_dia_ms1_{}_{:?}.tdf",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let conn = Connection::open(&path).unwrap();
        conn.execute_batch(
            "CREATE TABLE DiaFrameMsMsWindows (WindowGroup INTEGER, ScanNumBegin INTEGER, \
             ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, CollisionEnergy REAL);",
        )
        .unwrap();
        // Same window appears in two groups (must dedupe to one box); a second,
        // distinct window; and a degenerate row (scan_end <= begin) that is dropped.
        let rows = [
            (1i64, 100u32, 200u32, 500.0f64, 50.0f64),
            (2, 100, 200, 500.0, 50.0), // identical -> deduped
            (1, 200, 300, 700.0, 25.0),
            (1, 50, 50, 400.0, 10.0), // empty scan range -> dropped
        ];
        for &(g, sb, se, mz, w) in &rows {
            conn.execute(
                "INSERT INTO DiaFrameMsMsWindows \
                 (WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth) \
                 VALUES (?1, ?2, ?3, ?4, ?5)",
                rusqlite::params![g, sb as i64, se as i64, mz, w],
            )
            .unwrap();
        }
        drop(conn);

        let mut boxes = read_dia_ms1_boxes(&path).unwrap();
        boxes.sort_by(|a, b| a.mz_lo.partial_cmp(&b.mz_lo).unwrap());
        assert_eq!(boxes.len(), 2, "expected 2 boxes after dedupe + drop");
        assert_eq!((boxes[0].scan_begin, boxes[0].scan_end), (100, 200));
        assert!((boxes[0].mz_lo - 475.0).abs() < 1e-9);
        assert!((boxes[0].mz_hi - 525.0).abs() < 1e-9);
        assert_eq!((boxes[1].scan_begin, boxes[1].scan_end), (200, 300));
        assert!((boxes[1].mz_lo - 687.5).abs() < 1e-9);
        assert!((boxes[1].mz_hi - 712.5).abs() < 1e-9);

        let _ = std::fs::remove_file(&path);
    }
}
