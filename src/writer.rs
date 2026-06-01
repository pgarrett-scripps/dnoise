//! Orchestration: copy the source `.d`, rewrite `analysis.tdf_bin` with filtered
//! frames (re-encoded as type 2), and fix up the `analysis.tdf` SQLite database.

use crate::average::running_average;
use crate::filter::filter_iterated;
use crate::frame::FlatFrame;
use crate::params::FilterParams;
use crate::tdf::{self, FrameUpdate, encode::encode_frame_type2};
use anyhow::{Context, Result, bail};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::fs;
use std::io::{BufWriter, Read, Write};
use std::path::Path;
use timsrust::readers::FrameReader;

/// Frames are read+filtered+encoded in parallel batches of this size, then the
/// batch is written sequentially (so offsets stay ordered) before the next.
/// Bounds peak memory to roughly this many encoded frames.
const CHUNK: usize = 2048;

/// Summary returned to the CLI.
pub struct DenoiseStats {
    pub frames: usize,
    pub raw_points: u64,
    pub kept_points: u64,
}

/// Denoise `input` (.d) into a new `output` (.d).
///
/// `filter_all_frames`: when false (default), only MS1 frames are filtered and
/// MS/MS frames are re-encoded unchanged — the vertical-IM filter is an MS1
/// algorithm and strips most MS/MS fragment signal, collapsing DDA IDs. When
/// true, every frame is filtered.
///
/// `frame_half_width`: when > 0, each MS1 frame is replaced by the centered
/// running average of its `2*frame_half_width+1` MS1-frame neighborhood before
/// filtering (see [`crate::average`]). MS/MS frames are never averaged. 0 is the
/// default and reproduces the unsmoothed pipeline exactly.
pub fn denoise(
    input: &Path,
    output: &Path,
    params: &FilterParams,
    filter_all_frames: bool,
    frame_half_width: usize,
    force: bool,
) -> Result<DenoiseStats> {
    let in_tdf = input.join("analysis.tdf");
    let in_bin = input.join("analysis.tdf_bin");
    if !in_tdf.is_file() || !in_bin.is_file() {
        bail!(
            "{} is not a Bruker .d folder (missing analysis.tdf / analysis.tdf_bin)",
            input.display()
        );
    }
    if output.exists() {
        if force {
            fs::remove_dir_all(output)
                .with_context(|| format!("remove existing {}", output.display()))?;
        } else {
            bail!(
                "output {} already exists (use --force to overwrite)",
                output.display()
            );
        }
    }

    // Copy everything except the binary (we regenerate that), so calibration.sqlite,
    // analysis.tdf, etc. come along.
    copy_dir_except(input, output, "analysis.tdf_bin")?;

    let reader = FrameReader::new(input).map_err(|e| anyhow::anyhow!("open frames: {e}"))?;
    let n_frames = reader.len();
    // Frame metadata (ordered by Id == timsrust index). Empty frames are handled
    // without timsrust, which cannot decode their absent payload.
    let meta = tdf::read_frame_meta(&in_tdf)?;

    // MS1 subsequence used by the running-average pre-filter: `ms1_indices` maps a
    // position in the MS1-only stream to a global frame index; `ms1_pos` is the
    // reverse (None for MS/MS frames). The window skips interleaved MS/MS frames.
    let ms1_indices: Vec<usize> = (0..n_frames).filter(|&i| meta[i].is_ms1()).collect();
    let mut ms1_pos: Vec<Option<usize>> = vec![None; n_frames];
    for (p, &gi) in ms1_indices.iter().enumerate() {
        ms1_pos[gi] = Some(p);
    }

    let out_bin = output.join("analysis.tdf_bin");
    let mut bin = BufWriter::new(
        fs::File::create(&out_bin).with_context(|| format!("create {}", out_bin.display()))?,
    );

    // Preserve the leading header that precedes the first frame (Bruker reserves a
    // block at the start of the .tdf_bin), and start writing frames after it so the
    // new TimsId offsets land in the same layout the Bruker reader expects.
    let header_len = tdf::binary_header_len(&in_tdf)?;
    if header_len > 0 {
        let mut header = vec![0u8; header_len as usize];
        fs::File::open(&in_bin)
            .and_then(|mut f| f.read_exact(&mut header))
            .with_context(|| format!("read {header_len}-byte header from {}", in_bin.display()))?;
        bin.write_all(&header)?;
    }

    let pb = ProgressBar::new(n_frames as u64);
    pb.set_style(ProgressStyle::with_template("{bar:40} {pos}/{len} frames").unwrap());

    let mut offset: u64 = header_len;
    let mut updates: Vec<FrameUpdate> = Vec::with_capacity(n_frames);
    let mut raw_points: u64 = 0;
    let mut kept_points: u64 = 0;

    for start in (0..n_frames).step_by(CHUNK) {
        let end = (start + CHUNK).min(n_frames);
        let processed: Vec<ProcessedFrame> = (start..end)
            .into_par_iter()
            .map(|i| {
                process_frame(
                    &reader,
                    &meta,
                    i,
                    params,
                    filter_all_frames,
                    frame_half_width,
                    &ms1_indices,
                    &ms1_pos,
                )
            })
            .collect::<Result<_>>()?;

        for pf in processed {
            raw_points += pf.raw_points;
            kept_points += pf.num_peaks;
            bin.write_all(&pf.record)?;
            updates.push(FrameUpdate {
                frame_id: pf.frame_id,
                tims_id: offset,
                num_peaks: pf.num_peaks,
                max_intensity: pf.max_intensity,
                summed_intensities: pf.summed_intensities,
            });
            offset += pf.record.len() as u64;
            pb.inc(1);
        }
    }
    bin.flush()?;
    drop(bin);
    pb.finish_and_clear();

    tdf::update_metadata(&output.join("analysis.tdf"), &updates)?;

    Ok(DenoiseStats {
        frames: n_frames,
        raw_points,
        kept_points,
    })
}

struct ProcessedFrame {
    frame_id: usize,
    record: Vec<u8>,
    raw_points: u64,
    num_peaks: u64,
    max_intensity: u32,
    summed_intensities: u64,
}

#[allow(clippy::too_many_arguments)]
fn process_frame(
    reader: &FrameReader,
    meta: &[tdf::FrameMeta],
    i: usize,
    params: &FilterParams,
    filter_all_frames: bool,
    frame_half_width: usize,
    ms1_indices: &[usize],
    ms1_pos: &[Option<usize>],
) -> Result<ProcessedFrame> {
    let meta_i = &meta[i];
    // Empty frames: timsrust cannot decode their absent payload, so emit the
    // canonical empty record directly (Bruker stores these too).
    if meta_i.num_peaks == 0 {
        return Ok(ProcessedFrame {
            frame_id: meta_i.id,
            record: crate::tdf::encode::encode_empty_frame_type2(meta_i.num_scans),
            raw_points: 0,
            num_peaks: 0,
            max_intensity: 0,
            summed_intensities: 0,
        });
    }

    let frame = reader
        .get(i)
        .map_err(|e| anyhow::anyhow!("read frame {i}: {e}"))?;
    let flat = FlatFrame::from_frame(&frame);
    let raw_points = flat.len() as u64;
    let num_scans = flat.num_scans;
    let frame_id = flat.frame_id;

    // Optional pre-filter smoothing (decoupled from the filter): replace this MS1
    // frame with the centered running average of its MS1-frame neighborhood. The
    // filter then runs on the smoothed frame exactly as it would on a raw one.
    let smoothed;
    let to_filter: &FlatFrame = if frame_half_width > 0 && meta_i.is_ms1() {
        let p = ms1_pos[i].expect("MS1 frame must have an MS1-stream position");
        let lo = p.saturating_sub(frame_half_width);
        let hi = (p + frame_half_width).min(ms1_indices.len() - 1);
        let mut neighbors: Vec<FlatFrame> = Vec::with_capacity(hi - lo);
        for &gi in &ms1_indices[lo..=hi] {
            // The center frame is already flattened; empty neighbors contribute nothing.
            if gi == i || meta[gi].num_peaks == 0 {
                continue;
            }
            let nf = reader
                .get(gi)
                .map_err(|e| anyhow::anyhow!("read frame {gi}: {e}"))?;
            neighbors.push(FlatFrame::from_frame(&nf));
        }
        let mut window: Vec<&FlatFrame> = neighbors.iter().collect();
        window.push(&flat);
        smoothed = running_average(frame_id, num_scans, &window);
        &smoothed
    } else {
        &flat
    };

    // The vertical-IM filter is MS1-specific; MS/MS frames are re-encoded unchanged
    // unless `filter_all_frames` is set.
    let keep = if filter_all_frames || meta_i.is_ms1() {
        filter_iterated(to_filter, params)
    } else {
        vec![true; to_filter.len()]
    };
    let survivors = to_filter.survivors(&keep);

    let num_peaks = survivors.len() as u64;
    let summed_intensities: u64 = survivors.iter().map(|&(_, _, it)| it as u64).sum();
    let max_intensity = survivors.iter().map(|&(_, _, it)| it).max().unwrap_or(0);
    let record = encode_frame_type2(num_scans, &survivors);

    Ok(ProcessedFrame {
        frame_id,
        record,
        raw_points,
        num_peaks,
        max_intensity,
        summed_intensities,
    })
}

/// Recursively copy `src` into `dst`, skipping a top-level entry named `skip_top`.
fn copy_dir_except(src: &Path, dst: &Path, skip_top: &str) -> Result<()> {
    fs::create_dir_all(dst).with_context(|| format!("create {}", dst.display()))?;
    for entry in fs::read_dir(src)? {
        let entry = entry?;
        let name = entry.file_name();
        let from = entry.path();
        let to = dst.join(&name);
        if from.is_dir() {
            copy_dir_recursive(&from, &to)?;
        } else if name != skip_top {
            fs::copy(&from, &to).with_context(|| format!("copy {}", from.display()))?;
        }
    }
    Ok(())
}

fn copy_dir_recursive(src: &Path, dst: &Path) -> Result<()> {
    fs::create_dir_all(dst)?;
    for entry in fs::read_dir(src)? {
        let entry = entry?;
        let from = entry.path();
        let to = dst.join(entry.file_name());
        if from.is_dir() {
            copy_dir_recursive(&from, &to)?;
        } else {
            fs::copy(&from, &to).with_context(|| format!("copy {}", from.display()))?;
        }
    }
    Ok(())
}
