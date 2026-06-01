//! Type-2 `tdf_bin` frame codec.
//!
//! On-disk a type-2 frame record is:
//! `[u32 total_byte_count][u32 scan_count][ zstd( byte-transposed u32 array ) ]`.
//!
//! The decompressed payload is a logical array of `u32` values stored
//! *byte-transposed*: byte `j` of value `i` lives at offset `i + j*n`, where
//! `n` is the number of values. The logical array is laid out as:
//!
//! - `v[0]            = scan_count`
//! - `v[1 ..= scan_count-1] = 2 * peaks_in_scan[s]` for scans `0 ..= scan_count-2`
//!   (the last scan's count is implied by `peak_count`)
//! - `v[scan_count + 2*p]     = TOF delta for peak p` (cumulative within a scan; `tof = sum - 1`)
//! - `v[scan_count + 1 + 2*p] = intensity for peak p`
//!
//! This mirrors, inverted, timsrust's `read_scan_offsets` / `read_intensities`
//! / `read_tof_indices` and the `TdfBlob` byte layout. The [`decode_frame_type2`]
//! here is a faithful reimplementation of that read path, used to gate the
//! encoder via a round-trip test.

use anyhow::{Result, bail};

/// A decoded frame: scan count plus `(scan, tof, intensity)` points in scan/TOF order.
pub type DecodedFrame = (usize, Vec<(u32, u32, u32)>);

/// Two `u32` header fields precede the compressed payload.
const HEADER_BYTES: usize = 8;

/// Encode a frame into a complete type-2 `tdf_bin` record.
///
/// `points` are `(scan, tof, intensity)` triples in any order; they are grouped
/// by scan and sorted by TOF ascending internally.
pub fn encode_frame_type2(num_scans: usize, points: &[(u32, u32, u32)]) -> Vec<u8> {
    let peak_count = points.len();

    // Group points by scan via counting sort into a single flat buffer, avoiding
    // one heap-allocated Vec per scan (num_scans is ~700+, mostly empty).
    let mut counts = vec![0u32; num_scans];
    for &(scan, _, _) in points {
        counts[scan as usize] += 1;
    }
    let mut offsets = vec![0usize; num_scans + 1];
    for s in 0..num_scans {
        offsets[s + 1] = offsets[s] + counts[s] as usize;
    }
    let mut flat: Vec<(u32, u32)> = vec![(0, 0); peak_count];
    let mut cursor = offsets[..num_scans].to_vec(); // per-scan write head
    for &(scan, tof, intensity) in points {
        let s = scan as usize;
        flat[cursor[s]] = (tof, intensity);
        cursor[s] += 1;
    }
    // Sort each scan's run by TOF ascending.
    for s in 0..num_scans {
        flat[offsets[s]..offsets[s + 1]].sort_unstable_by_key(|&(tof, _)| tof);
    }

    let n = num_scans + 2 * peak_count;
    let mut values = vec![0u32; n];
    values[0] = num_scans as u32;

    // Per-scan sizes for scans 0 ..= num_scans-2 (the read path stops one short).
    for s in 0..num_scans.saturating_sub(1) {
        values[s + 1] = counts[s] * 2;
    }

    // Interleaved TOF-delta / intensity in scan order, TOF ascending within a scan.
    let mut p = 0usize;
    for s in 0..num_scans {
        let mut current_sum: u32 = 0;
        for &(tof, intensity) in &flat[offsets[s]..offsets[s + 1]] {
            let delta = (tof + 1) - current_sum; // safe: run is sorted ascending
            current_sum = tof + 1;
            values[num_scans + 2 * p] = delta;
            values[num_scans + 1 + 2 * p] = intensity;
            p += 1;
        }
    }

    // Byte-transpose: byte j of value i -> position i + j*n.
    let mut transposed = vec![0u8; n * 4];
    for (i, &v) in values.iter().enumerate() {
        transposed[i] = (v & 0xFF) as u8;
        transposed[i + n] = ((v >> 8) & 0xFF) as u8;
        transposed[i + 2 * n] = ((v >> 16) & 0xFF) as u8;
        transposed[i + 3 * n] = ((v >> 24) & 0xFF) as u8;
    }

    // Use the one-shot (bulk) encoder, not the streaming `encode_all`: it records the
    // decompressed size in the zstd frame header (single-segment, FCS present). Bruker's
    // reader requires a known content size and rejects frames without it.
    let compressed = zstd::bulk::compress(&transposed, 1).expect("zstd encode is infallible here");
    let total_byte_count = (HEADER_BYTES + compressed.len()) as u32;
    let mut record = Vec::with_capacity(total_byte_count as usize);
    record.extend_from_slice(&total_byte_count.to_le_bytes());
    record.extend_from_slice(&(num_scans as u32).to_le_bytes());
    record.extend_from_slice(&compressed);
    record
}

/// Encode an empty (0-peak) frame exactly as Bruker stores it: an 8-byte record
/// `[u32 total_byte_count = 8][u32 scan_count]` with no compressed payload.
///
/// timsTOF DDA files contain empty MS/MS frames; timsrust errors trying to
/// zstd-decode their absent payload, so we never round-trip these through the
/// reader — we emit the canonical empty record directly.
pub fn encode_empty_frame_type2(num_scans: usize) -> Vec<u8> {
    let mut record = Vec::with_capacity(HEADER_BYTES);
    record.extend_from_slice(&(HEADER_BYTES as u32).to_le_bytes());
    record.extend_from_slice(&(num_scans as u32).to_le_bytes());
    record
}

/// Decode a type-2 record back to `(num_scans, points)`, where `points` are
/// `(scan, tof, intensity)` triples in scan/TOF order. Faithful reimplementation
/// of timsrust's read path; used for round-trip testing.
pub fn decode_frame_type2(record: &[u8]) -> Result<DecodedFrame> {
    if record.len() < HEADER_BYTES {
        bail!("record shorter than header");
    }
    let total_byte_count = u32::from_le_bytes(record[0..4].try_into().unwrap()) as usize;
    if total_byte_count > record.len() || total_byte_count < HEADER_BYTES {
        bail!("invalid total_byte_count {total_byte_count}");
    }
    let payload = &record[HEADER_BYTES..total_byte_count];
    let bytes = zstd::decode_all(payload)?;
    if !bytes.len().is_multiple_of(4) {
        bail!("decompressed blob not u32-aligned");
    }
    let n = bytes.len() / 4;
    let get = |i: usize| -> u32 {
        bytes[i] as u32
            | (bytes[i + n] as u32) << 8
            | (bytes[i + 2 * n] as u32) << 16
            | (bytes[i + 3 * n] as u32) << 24
    };

    let scan_count = get(0) as usize;
    if scan_count > n {
        bail!("scan_count {scan_count} exceeds blob length {n}");
    }
    let peak_count = (n - scan_count) / 2;

    let mut scan_offsets = Vec::with_capacity(scan_count + 1);
    scan_offsets.push(0usize);
    for scan_index in 0..scan_count.saturating_sub(1) {
        let size = (get(scan_index + 1) / 2) as usize;
        scan_offsets.push(scan_offsets[scan_index] + size);
    }
    scan_offsets.push(peak_count);

    let mut points = Vec::with_capacity(peak_count);
    for scan_index in 0..scan_count {
        let start = scan_offsets[scan_index];
        let end = scan_offsets[scan_index + 1];
        let mut current_sum: u32 = 0;
        for peak_index in start..end {
            current_sum += get(scan_count + 2 * peak_index);
            let tof = current_sum - 1;
            let intensity = get(scan_count + 1 + 2 * peak_index);
            points.push((scan_index as u32, tof, intensity));
        }
    }
    Ok((scan_count, points))
}
