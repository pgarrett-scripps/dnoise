//! Read-only structural and decoded-content validation of type-2 acquisitions.
use crate::{DnoiseError, Result};
use rayon::prelude::*;
use rusqlite::{Connection, OpenFlags};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

// Keep compressed buffering and simultaneous decoded allocations bounded even
// on machines with large Rayon pools. An oversized record runs by itself.
const BATCH_BYTES: usize = 16 * 1024 * 1024;
const BATCH_FRAMES: usize = 256;
const DECODE_WORKERS: usize = 4;

struct Record {
    id: i64,
    offset: u64,
    length: usize,
    scans: usize,
    peaks: i64,
}

fn check_cancel(cancel: Option<&AtomicBool>) -> Result<()> {
    if cancel.is_some_and(|c| c.load(Ordering::Relaxed)) {
        return Err(DnoiseError::Cancelled);
    }
    Ok(())
}

/// Summary of a validated acquisition. SDK validation remains an independent check.
#[derive(Debug, Clone, Copy)]
pub struct Validation {
    /// Number of frames, including empty frames.
    pub frames: usize,
    /// Total points according to the validated metadata.
    pub points: u64,
    /// Actual frame binary size.
    pub binary_bytes: u64,
}

/// Check SQLite integrity, compression, frame offsets/counts, and optionally every
/// decoded record. No files or metadata are changed. Supports compression type 2.
/// Decoding uses bounded batches and up to four tasks in the current Rayon pool.
pub fn inspect(input: &Path, decode: bool) -> Result<Validation> {
    inspect_cancellable(input, decode, None)
}

pub(crate) fn inspect_cancellable(
    input: &Path,
    decode: bool,
    cancel: Option<&AtomicBool>,
) -> Result<Validation> {
    check_cancel(cancel)?;
    let tdf = input.join("analysis.tdf");
    let bin = input.join("analysis.tdf_bin");
    if !tdf.is_file() || !bin.is_file() {
        return Err(DnoiseError::NotADotD(input.to_owned()));
    }
    let db = Connection::open_with_flags(&tdf, OpenFlags::SQLITE_OPEN_READ_ONLY)?;
    let integrity: String = db.query_row("PRAGMA quick_check", [], |r| r.get(0))?;
    if integrity != "ok" {
        return Err(DnoiseError::InvalidInput(format!(
            "SQLite integrity: {integrity}"
        )));
    }
    let compression: String = db.query_row(
        "SELECT Value FROM GlobalMetadata WHERE Key='TimsCompressionType'",
        [],
        |r| r.get(0),
    )?;
    if compression != "2" {
        return Err(DnoiseError::InvalidInput(format!(
            "unsupported compression type {compression}; expected 2"
        )));
    }
    let mut file = File::open(bin)?;
    let size = file.metadata()?.len();
    let mut stmt =
        db.prepare("SELECT Id, TimsId, NumScans, NumPeaks, Time FROM Frames ORDER BY Id")?;
    let mut rows = stmt.query([])?;
    let mut frames = 0;
    let mut points = 0u64;
    let mut records = Vec::new();
    while let Some(row) = rows.next()? {
        check_cancel(cancel)?;
        let id: i64 = row.get(0)?;
        let offset: i64 = row.get(1)?;
        let scans: i64 = row.get(2)?;
        let peaks: i64 = row.get(3)?;
        let time: f64 = row.get(4)?;
        if id != frames as i64 + 1
            || offset < 0
            || scans < 0
            || scans > u32::MAX as i64
            || peaks < 0
            || !time.is_finite()
        {
            return Err(DnoiseError::InvalidInput(format!(
                "invalid frame metadata at frame {id}"
            )));
        }
        let offset = offset as u64;
        if offset.saturating_add(8) > size {
            return Err(DnoiseError::InvalidInput(format!(
                "invalid binary offset for frame {id}"
            )));
        }
        file.seek(SeekFrom::Start(offset))?;
        let mut header = [0u8; 8];
        file.read_exact(&mut header)?;
        let length = u32::from_le_bytes(header[..4].try_into().unwrap()) as u64;
        let header_scans = u32::from_le_bytes(header[4..].try_into().unwrap());
        if length < 8
            || offset.saturating_add(length) > size
            || header_scans as i64 != scans
            || (peaks > 0 && scans == 0)
        {
            return Err(DnoiseError::InvalidInput(format!(
                "invalid record header for frame {id}"
            )));
        }
        if decode && length > crate::codec::MAX_DECODE_BYTES as u64 {
            return Err(DnoiseError::InvalidInput(
                "frame exceeds validation size limit".into(),
            ));
        }
        records.push(Record {
            id,
            offset,
            length: length as usize,
            scans: scans as usize,
            peaks,
        });
        frames += 1;
        points = points
            .checked_add(peaks as u64)
            .ok_or_else(|| DnoiseError::InvalidInput("point count overflow".into()))?;
    }
    records.sort_unstable_by_key(|r| r.offset);
    if records
        .windows(2)
        .any(|w| w[0].offset + w[0].length as u64 > w[1].offset)
    {
        return Err(DnoiseError::InvalidInput(
            "overlapping frame records".into(),
        ));
    }
    if frames == 0 {
        return Err(DnoiseError::InvalidInput(
            "acquisition has no frames".into(),
        ));
    }
    crate::detect_acquisition(input)?;
    if decode {
        validate_records(&mut file, &records, cancel)?;
    }
    check_cancel(cancel)?;
    Ok(Validation {
        frames,
        points,
        binary_bytes: size,
    })
}

// One reader visits records in physical file order. Rayon handles only the
// independent decoding work, avoiding shared seek cursors or concurrent I/O.
fn validate_records(
    file: &mut File,
    records: &[Record],
    cancel: Option<&AtomicBool>,
) -> Result<()> {
    let workers = rayon::current_num_threads().min(DECODE_WORKERS);
    let mut remaining = records;
    while !remaining.is_empty() {
        let mut batch = Vec::new();
        let mut bytes = 0;
        for record in remaining.iter().take(BATCH_FRAMES) {
            if !batch.is_empty() && bytes + record.length > BATCH_BYTES {
                break;
            }
            check_cancel(cancel)?;
            file.seek(SeekFrom::Start(record.offset))?;
            let mut encoded = vec![0; record.length];
            file.read_exact(&mut encoded)?;
            bytes += record.length;
            batch.push((record, encoded));
        }
        // At most four tasks, regardless of the caller's pool size. Each task
        // drops its decoded frame before starting another; no decoded run is kept.
        batch
            .par_chunks(batch.len().div_ceil(workers))
            .try_for_each(|chunk| {
                for (record, encoded) in chunk {
                    check_cancel(cancel)?;
                    let (scans, decoded) = crate::codec::decode_frame_type2(encoded)?;
                    if scans != record.scans || decoded.len() as i64 != record.peaks {
                        return Err(DnoiseError::InvalidInput(format!(
                            "metadata/decoded counts disagree for frame {}",
                            record.id
                        )));
                    }
                }
                Ok(())
            })?;
        remaining = &remaining[batch.len()..];
    }
    Ok(())
}

/// Validate supplied processing parameters consistently for all front ends.
pub fn parameters(
    params: &crate::FilterParams,
    stages: &crate::Stages,
    options: &crate::RunOptions,
) -> Result<()> {
    if options.frame_batch_size == Some(0) {
        return Err(DnoiseError::InvalidInput(
            "frame batch size must be positive".into(),
        ));
    }
    let invalid = |message: &str| DnoiseError::InvalidInput(message.into());
    if !stages.neighbors.max_rt_gap_seconds.is_finite()
        || stages.neighbors.max_rt_gap_seconds <= 0.0
    {
        return Err(invalid(
            "neighbor_max_rt_gap must be finite and positive (seconds)",
        ));
    }
    if params.max_internal_gap == usize::MAX {
        return Err(invalid("max_internal_gap is too large"));
    }
    if stages.watershed.is_some() && stages.box_centroid.is_some() {
        return Err(invalid("choose only one centroider"));
    }
    if let Some(h) = stages.halo {
        if !h.peak_fraction.is_finite()
            || !(0.0..=1.0).contains(&h.peak_fraction)
            || h.scan_half_width > i64::MAX as usize
        {
            return Err(invalid(
                "halo fraction must be finite in [0,1], with a valid scan width",
            ));
        }
    }
    for (mz, im) in stages
        .ms1_polygon
        .map(|p| (p.mz_pad, p.im_pad))
        .into_iter()
        .chain(stages.dia_ms1.map(|p| (p.mz_pad, p.im_pad)))
    {
        if !mz.is_finite() || !im.is_finite() || mz < 0.0 || im < 0.0 {
            return Err(invalid("gate padding must be finite and nonnegative"));
        }
    }
    if options.crop_only && options.crop.is_none_or(|c| c.is_empty()) {
        return Err(invalid("crop-only requires at least one crop bound"));
    }
    if let Some(c) = options.crop {
        for (lo, hi) in [
            (c.mz_min, c.mz_max),
            (c.im_min, c.im_max),
            (c.rt_min, c.rt_max),
        ] {
            if lo.into_iter().chain(hi).any(|x| !x.is_finite() || x < 0.0)
                || lo.zip(hi).is_some_and(|(l, h)| l > h)
            {
                return Err(invalid(
                    "crop bounds must be finite, nonnegative, and ordered",
                ));
            }
        }
        if c.min_intensity
            .zip(c.max_intensity)
            .is_some_and(|(l, h)| l > h)
        {
            return Err(invalid("intensity crop bounds are reversed"));
        }
    }
    if let Some(s) = options.sample {
        if !options.dry_run || !s.fraction.is_finite() || s.fraction <= 0.0 || s.fraction > 1.0 {
            return Err(invalid("sampling requires dry-run and a fraction in (0,1]"));
        }
    }
    Ok(())
}

/// Read one zero-based frame with the checked codec (including empty frames).
/// Useful for inspection without running the denoising pipeline again.
pub fn read_frame(input: &Path, index: usize) -> Result<crate::codec::DecodedFrame> {
    let db =
        Connection::open_with_flags(input.join("analysis.tdf"), OpenFlags::SQLITE_OPEN_READ_ONLY)?;
    let offset: i64 = db.query_row(
        "SELECT TimsId FROM Frames ORDER BY Id LIMIT 1 OFFSET ?1",
        [index as i64],
        |r| r.get(0),
    )?;
    if offset < 0 {
        return Err(DnoiseError::InvalidInput("negative frame offset".into()));
    }
    let mut file = File::open(input.join("analysis.tdf_bin"))?;
    file.seek(SeekFrom::Start(offset as u64))?;
    let mut header = [0u8; 8];
    file.read_exact(&mut header)?;
    let length = u32::from_le_bytes(header[..4].try_into().unwrap()) as usize;
    if !(8..=crate::codec::MAX_DECODE_BYTES).contains(&length) {
        return Err(DnoiseError::InvalidInput("invalid frame size".into()));
    }
    let mut record = vec![0; length];
    record[..8].copy_from_slice(&header);
    file.read_exact(&mut record[8..])?;
    Ok(crate::codec::decode_frame_type2(&record)?)
}
