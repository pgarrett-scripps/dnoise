//! Orchestration: copy the source `.d`, rewrite `analysis.tdf_bin` with filtered
//! frames (re-encoded as type 2), and fix up the `analysis.tdf` SQLite database.

use crate::average::running_average;
use crate::box_centroid::box_centroid;
use crate::codec::encode_frame_type2;
use crate::dia_ms1::{DiaMs1Gate, TofScanBox};
use crate::dia_window::{filter_per_window, in_window_mask};
use crate::error::{DnoiseError, Result};
use crate::filter::filter_iterated;
use crate::frame::FlatFrame;
use crate::halo::horizontal_halo_keep_mask;
use crate::msms::{MsmsKeep, build_msms_keep};
use crate::params::{
    BoxCentroidParams, DiaMs1WindowParams, DiaWindowParams, FilterParams, HaloParams,
    MsmsFilterParams, SmoothParams, WatershedParams,
};
use crate::smooth::box_average;
use crate::tdf::{self, DiaWindows, FrameUpdate};
use crate::watershed::watershed_centroid;
use rayon::prelude::*;
use std::fs;
use std::io::{BufWriter, Read, Write};
use std::path::Path;
use timsrust::converters::ConvertableDomain;
use timsrust::readers::{FrameReader, MetadataReader};

/// Frames are read+filtered+encoded in parallel batches of this size, then the
/// batch is written sequentially (so offsets stay ordered) before the next.
/// Bounds peak memory to roughly this many encoded frames.
const CHUNK: usize = 2048;

/// Upper bound on watershed groups formed per frame — a guard against
/// pathological frames. Real MS1 frames centroid to far fewer than this.
const MAX_CENTROIDS: usize = 100_000;

/// Summary returned by [`denoise`].
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct DenoiseStats {
    /// Total frames processed (MS1 + MS/MS + empty).
    pub frames: usize,
    /// Total input points across all frames.
    pub raw_points: u64,
    /// Total points written after filtering.
    pub kept_points: u64,
}

/// Progress update passed to the callback of [`denoise_with_progress`].
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct Progress {
    /// Frames written so far.
    pub frames_done: usize,
    /// Total frames to process.
    pub frames_total: usize,
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
///
/// `halo`: when `Some`, the horizontal-halo filter ([`crate::halo`]) runs after
/// the vertical filter to remove the weak m/z halo flanking bright peaks.
/// `None` disables it.
///
/// `denoise_msms`: when `Some`, MS/MS frames are denoised with these filter
/// knobs instead of being passed through (`None` leaves them unchanged). The
/// acquisition scheme is auto-detected: **ddaPASEF** combines each precursor's
/// fragment scans across frames before filtering (see [`crate::msms`]);
/// **diaPASEF** (no `PasefFrameMsMsInfo`) runs the same filter on each whole
/// MS/MS frame as-is, since each isolation window is sampled once per cycle.
///
/// `smooth`: when `Some`, the box-averaging smoother ([`crate::smooth`]) runs on
/// each filtered frame's survivors — after the halo filter, before `watershed` —
/// rewriting intensities with their `(scan, TOF index)` box mean to stabilise
/// watershed seeding. Applied to the same frames as the vertical filter. `None`
/// (default) disables it.
///
/// `watershed`: when `Some`, the watershed centroider ([`crate::watershed`])
/// runs as a final stage on each filtered frame's survivors, collapsing them
/// into intensity-weighted centroids. Applied to the same frames the vertical
/// filter touches (MS1 always; MS/MS only when `filter_all_frames`). `None`
/// (default) disables it, leaving the filtered points as-is.
///
/// `box_centroid_params`: when `Some`, the greedy small-box centroider
/// ([`crate::box_centroid`]) runs as the final stage instead — consolidating
/// points within small fixed boxes (tiling streaks rather than collapsing them).
/// Mutually exclusive with `watershed`. Applied to the same frames. `None`
/// (default) disables it.
///
/// `dia_window`: when `Some` (diaPASEF only), MS/MS points whose mobility scan
/// falls outside every isolation window for their frame are dropped (the
/// [`crate::dia_window::in_window_mask`] gate, with the struct's `scan_pad`). Has
/// no effect on ddaPASEF data (no `DiaFrameMsMs*` tables). `None` disables it.
///
/// `dia_per_window`: when `true` (diaPASEF only), the MS/MS vertical/halo filter
/// is run independently inside each isolation window's scan slice
/// ([`crate::dia_window::filter_per_window`]) instead of over the whole frame,
/// preventing the filter from fusing a mobility run across a window boundary
/// (cross-talk between unrelated isolation events). Applies wherever the MS/MS
/// filter runs (`denoise_msms` on diaPASEF, or `filter_all_frames`); ignored on
/// ddaPASEF and when no MS/MS filter runs.
///
/// `dia_ms1`: when `Some` (diaPASEF only), **MS1** points that fall outside every
/// isolation window's `(m/z, mobility)` region are dropped ([`crate::dia_ms1`]).
/// The windows (`DiaFrameMsMsWindows`) are padded by the struct's `mz_pad` (Da) and
/// `im_pad` (1/K0) — converted to TOF indices / scans once via the run's
/// calibration — so a precursor near a window edge keeps its full isotopic
/// envelope. Has no effect on ddaPASEF data (no `DiaFrameMsMs*` tables). `None`
/// disables it.
///
/// This reports no progress; use [`denoise_with_progress`] to receive
/// [`Progress`] updates as frames are written.
#[allow(clippy::too_many_arguments)]
pub fn denoise(
    input: &Path,
    output: &Path,
    params: &FilterParams,
    filter_all_frames: bool,
    frame_half_width: usize,
    halo: Option<&HaloParams>,
    denoise_msms: Option<&MsmsFilterParams>,
    smooth: Option<&SmoothParams>,
    watershed: Option<&WatershedParams>,
    box_centroid_params: Option<&BoxCentroidParams>,
    dia_window: Option<&DiaWindowParams>,
    dia_per_window: bool,
    dia_ms1: Option<&DiaMs1WindowParams>,
    force: bool,
) -> Result<DenoiseStats> {
    denoise_with_progress(
        input,
        output,
        params,
        filter_all_frames,
        frame_half_width,
        halo,
        denoise_msms,
        smooth,
        watershed,
        box_centroid_params,
        dia_window,
        dia_per_window,
        dia_ms1,
        force,
        |_| {},
    )
}

/// Like [`denoise`], but invokes `progress` once before processing and again
/// after each frame is written, so callers (e.g. a CLI) can drive a progress bar
/// without the library depending on any UI crate.
#[allow(clippy::too_many_arguments)]
pub fn denoise_with_progress<F: FnMut(Progress)>(
    input: &Path,
    output: &Path,
    params: &FilterParams,
    filter_all_frames: bool,
    frame_half_width: usize,
    halo: Option<&HaloParams>,
    denoise_msms: Option<&MsmsFilterParams>,
    smooth: Option<&SmoothParams>,
    watershed: Option<&WatershedParams>,
    box_centroid_params: Option<&BoxCentroidParams>,
    dia_window: Option<&DiaWindowParams>,
    dia_per_window: bool,
    dia_ms1: Option<&DiaMs1WindowParams>,
    force: bool,
    mut progress: F,
) -> Result<DenoiseStats> {
    let in_tdf = input.join("analysis.tdf");
    let in_bin = input.join("analysis.tdf_bin");
    if !in_tdf.is_file() || !in_bin.is_file() {
        return Err(DnoiseError::NotADotD(input.to_path_buf()));
    }
    if output.exists() {
        if force {
            fs::remove_dir_all(output)?;
        } else {
            return Err(DnoiseError::OutputExists(output.to_path_buf()));
        }
    }

    // Copy everything except the binary (we regenerate that), so calibration.sqlite,
    // analysis.tdf, etc. come along.
    copy_dir_except(input, output, "analysis.tdf_bin")?;

    let reader = FrameReader::new(input).map_err(|e| DnoiseError::OpenFrames(e.to_string()))?;
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

    // MS/MS denoising splits by acquisition scheme, both driven by the same
    // `denoise_msms` params:
    //   * ddaPASEF — each precursor is re-isolated across several frames, so we
    //     build per-precursor keep sets up front (PasefFrameMsMsInfo) and combine
    //     a precursor's fragment scans across frames before filtering.
    //   * diaPASEF — has no PasefFrameMsMsInfo (each isolation window is sampled
    //     once per cycle, nothing to combine). We detect this by an empty
    //     ddaPASEF window table and instead run the same MS/MS filter on each
    //     whole MS/MS frame as-is (the `dia_msms` path in `process_frame`).
    let (msms_keep, dia_msms) = match denoise_msms {
        Some(mp) => {
            let windows = tdf::read_pasef_msms(&in_tdf)?;
            if windows.is_empty() {
                (None, Some(mp))
            } else {
                (
                    Some(build_msms_keep(&reader, &meta, &windows, mp, halo)?),
                    None,
                )
            }
        }
        None => (None, None),
    };
    let msms_ref = msms_keep.as_ref();

    // diaPASEF isolation windows, read once and shared. Needed by either the
    // out-of-window gate (`dia_window`) or per-window MS/MS filtering
    // (`dia_per_window`). Empty for ddaPASEF, so both features no-op there.
    let dia_windows: Option<DiaWindows> = if dia_window.is_some() || dia_per_window {
        let w = tdf::read_dia_windows(&in_tdf)?;
        if w.is_empty() { None } else { Some(w) }
    } else {
        None
    };
    let dia_windows_ref = dia_windows.as_ref();

    // diaPASEF MS1 out-of-window gate: build the padded `(scan, TOF)` lookup once
    // from the isolation windows + calibration. `None` for ddaPASEF (no windows).
    let dia_ms1_gate = match dia_ms1 {
        Some(mp) => build_dia_ms1_gate(&in_tdf, mp, &meta)?,
        None => None,
    };
    let dia_ms1_ref = dia_ms1_gate.as_ref();

    let out_bin = output.join("analysis.tdf_bin");
    let mut bin = BufWriter::new(fs::File::create(&out_bin)?);

    // Preserve the leading header that precedes the first frame (Bruker reserves a
    // block at the start of the .tdf_bin), and start writing frames after it so the
    // new TimsId offsets land in the same layout the Bruker reader expects.
    let header_len = tdf::binary_header_len(&in_tdf)?;
    if header_len > 0 {
        let mut header = vec![0u8; header_len as usize];
        fs::File::open(&in_bin).and_then(|mut f| f.read_exact(&mut header))?;
        bin.write_all(&header)?;
    }

    progress(Progress {
        frames_done: 0,
        frames_total: n_frames,
    });

    let mut offset: u64 = header_len;
    let mut updates: Vec<FrameUpdate> = Vec::with_capacity(n_frames);
    let mut raw_points: u64 = 0;
    let mut kept_points: u64 = 0;
    let mut frames_done: usize = 0;

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
                    halo,
                    msms_ref,
                    dia_msms,
                    smooth,
                    watershed,
                    box_centroid_params,
                    dia_windows_ref,
                    dia_window,
                    dia_per_window,
                    dia_ms1_ref,
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
            frames_done += 1;
            progress(Progress {
                frames_done,
                frames_total: n_frames,
            });
        }
    }
    bin.flush()?;
    drop(bin);

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
    halo: Option<&HaloParams>,
    msms: Option<&MsmsKeep>,
    dia_msms: Option<&MsmsFilterParams>,
    smooth: Option<&SmoothParams>,
    watershed: Option<&WatershedParams>,
    box_centroid_params: Option<&BoxCentroidParams>,
    dia_windows: Option<&DiaWindows>,
    dia_window: Option<&DiaWindowParams>,
    dia_per_window: bool,
    dia_ms1: Option<&DiaMs1Gate>,
    ms1_indices: &[usize],
    ms1_pos: &[Option<usize>],
) -> Result<ProcessedFrame> {
    let meta_i = &meta[i];
    // Empty frames: timsrust cannot decode their absent payload, so emit the
    // canonical empty record directly (Bruker stores these too).
    if meta_i.num_peaks == 0 {
        return Ok(ProcessedFrame {
            frame_id: meta_i.id,
            record: crate::codec::encode_empty_frame_type2(meta_i.num_scans),
            raw_points: 0,
            num_peaks: 0,
            max_intensity: 0,
            summed_intensities: 0,
        });
    }

    let frame = reader.get(i).map_err(|e| DnoiseError::FrameRead {
        index: i,
        message: e.to_string(),
    })?;
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
            let nf = reader.get(gi).map_err(|e| DnoiseError::FrameRead {
                index: gi,
                message: e.to_string(),
            })?;
            neighbors.push(FlatFrame::from_frame(&nf));
        }
        let mut window: Vec<&FlatFrame> = neighbors.iter().collect();
        window.push(&flat);
        smoothed = running_average(frame_id, num_scans, &window);
        &smoothed
    } else {
        &flat
    };

    // diaPASEF isolation-window scan intervals for this frame (None for MS1, for
    // ddaPASEF, or when neither DIA feature is enabled). Drives both per-window
    // MS/MS filtering and the out-of-window gate below.
    let dia_iv = dia_windows.and_then(|dw| dw.intervals(meta_i.id));

    // MS1 frames: vertical filter then (optional) horizontal-halo on the survivors.
    // MS/MS frames: pruned by the precursor keep sets when MS/MS denoising is on,
    // otherwise re-encoded unchanged (or vertical-filtered if `filter_all_frames`).
    let mut keep = if meta_i.is_ms1() {
        let mut keep = filter_iterated(to_filter, params);
        if let Some(hp) = halo {
            apply_halo(to_filter, hp, &mut keep);
        }
        // diaPASEF MS1 out-of-window gate: drop surviving points whose (scan, TOF)
        // is in no padded isolation window. Composes as an AND on the keep mask.
        if let Some(gate) = dia_ms1 {
            for (slot, in_win) in keep
                .iter_mut()
                .zip(gate.keep_mask(&to_filter.scan, &to_filter.tof))
            {
                *slot &= in_win;
            }
        }
        keep
    } else if let Some(mk) = msms {
        mk.keep_mask(to_filter, meta_i.id)
    } else if let Some(mp) = dia_msms {
        // diaPASEF MS/MS: run the same MS/MS filter on each whole frame. With
        // `dia_per_window`, filter each isolation window's scan slice on its own
        // instead, so a mobility run cannot be fused across a window boundary
        // (cross-talk between unrelated isolation events). No cross-frame combine
        // (each DIA window is sampled once per cycle); the `msms_*` knobs apply via
        // FilterParams just like the ddaPASEF path.
        let fp = mp.as_filter_params();
        match dia_iv {
            Some(iv) if dia_per_window => filter_per_window(to_filter, iv, &fp, halo),
            _ => {
                let mut keep = filter_iterated(to_filter, &fp);
                if let Some(hp) = halo {
                    apply_halo(to_filter, hp, &mut keep);
                }
                keep
            }
        }
    } else if filter_all_frames {
        match dia_iv {
            Some(iv) if dia_per_window => filter_per_window(to_filter, iv, params, halo),
            _ => {
                let mut keep = filter_iterated(to_filter, params);
                if let Some(hp) = halo {
                    apply_halo(to_filter, hp, &mut keep);
                }
                keep
            }
        }
    } else {
        vec![true; to_filter.len()]
    };

    // Out-of-window gate: drop any MS/MS point whose scan falls outside every
    // isolation window for this frame. Independent of the streak filter, so it
    // also trims mobility-edge noise when no MS/MS filtering runs. (Per-window
    // filtering already excludes these points, making this a no-op there.)
    if let (Some(dp), Some(iv)) = (dia_window, dia_iv) {
        let mask = in_window_mask(&to_filter.scan, iv, dp.scan_pad);
        for (slot, keep_pt) in keep.iter_mut().zip(mask) {
            *slot &= keep_pt;
        }
    }
    let survivors = to_filter.survivors(&keep);

    // Optional intensity smoothing then watershed centroiding, both applied only
    // to frames the vertical filter actually processed — MS1 always, MS/MS only
    // under `filter_all_frames` (and never on the separate per-precursor MS/MS-
    // denoise path, which has its own keep logic). Smoothing runs first so the
    // watershed seeds on the stabilised intensities.
    let filtered_here =
        meta_i.is_ms1() || dia_msms.is_some() || (msms.is_none() && filter_all_frames);
    let survivors = match smooth {
        Some(sp) if filtered_here => box_average(&survivors, num_scans, sp),
        _ => survivors,
    };
    let survivors = match watershed {
        Some(wp) if filtered_here => watershed_centroid(&survivors, wp, MAX_CENTROIDS),
        _ => survivors,
    };
    // Optional final stage: greedy small-box centroiding (mutually exclusive with
    // watershed, enforced by the CLI). Tiles streaks into small centroids rather
    // than collapsing them, preserving the mobility profile.
    let survivors = match box_centroid_params {
        Some(bp) if filtered_here => box_centroid(&survivors, bp),
        _ => survivors,
    };

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

/// Build the diaPASEF MS1 out-of-window gate: read the isolation windows, pad each
/// in physical units (`mz_pad` Da, `im_pad` 1/K0) using the run's calibration,
/// convert to integer `(scan, TOF index)` boxes, and assemble the per-scan lookup.
/// Returns `None` for ddaPASEF (no windows) so the gate is skipped.
fn build_dia_ms1_gate(
    in_tdf: &Path,
    p: &DiaMs1WindowParams,
    meta: &[tdf::FrameMeta],
) -> Result<Option<DiaMs1Gate>> {
    let boxes = tdf::read_dia_ms1_boxes(in_tdf)?;
    if boxes.is_empty() {
        return Ok(None);
    }
    let md = MetadataReader::new(in_tdf).map_err(|e| DnoiseError::Metadata(e.to_string()))?;
    let num_scans = meta.iter().map(|m| m.num_scans).max().unwrap_or(0);
    if num_scans == 0 {
        return Ok(None);
    }

    let tof_boxes: Vec<TofScanBox> = boxes
        .iter()
        .map(|b| {
            // m/z edges -> TOF indices (monotonic), padded by mz_pad Da on each side.
            let t0 = md.mz_converter.invert(b.mz_lo - p.mz_pad);
            let t1 = md.mz_converter.invert(b.mz_hi + p.mz_pad);
            let tof_lo = t0.min(t1).floor().max(0.0) as u32;
            let tof_hi = t1.max(t0).ceil().max(0.0) as u32;

            // Scan range -> 1/K0 (monotonic decreasing), padded by im_pad, back to
            // scans. Take min/max so the result is correct regardless of direction.
            let im0 = md.im_converter.convert(b.scan_begin);
            let im1 = md.im_converter.convert(b.scan_end);
            let s0 = md.im_converter.invert(im0.max(im1) + p.im_pad);
            let s1 = md.im_converter.invert(im0.min(im1) - p.im_pad);
            let scan_lo = s0.min(s1).floor().max(0.0) as u32;
            let scan_hi = (s0.max(s1).ceil().max(0.0) as u32).min(num_scans as u32 - 1);

            TofScanBox {
                scan_lo,
                scan_hi,
                tof_lo,
                tof_hi,
            }
        })
        .collect();

    Ok(DiaMs1Gate::build(&tof_boxes, num_scans))
}

/// Run the horizontal-halo filter on the currently-kept points of `frame` and
/// turn off `keep` for any the filter removes. Operates in integer
/// `(scan, TOF index)` space — no calibration needed.
fn apply_halo(frame: &FlatFrame, hp: &HaloParams, keep: &mut [bool]) {
    let idx: Vec<usize> = (0..frame.len()).filter(|&i| keep[i]).collect();
    if idx.is_empty() {
        return;
    }
    let scan: Vec<u32> = idx.iter().map(|&i| frame.scan[i]).collect();
    let tof: Vec<u32> = idx.iter().map(|&i| frame.tof[i]).collect();
    let inten: Vec<u32> = idx.iter().map(|&i| frame.intensity[i]).collect();

    let hmask = horizontal_halo_keep_mask(&scan, &tof, &inten, frame.num_scans, hp);
    for (k, &i) in idx.iter().enumerate() {
        if !hmask[k] {
            keep[i] = false;
        }
    }
}

/// Recursively copy `src` into `dst`, skipping a top-level entry named `skip_top`.
fn copy_dir_except(src: &Path, dst: &Path, skip_top: &str) -> Result<()> {
    fs::create_dir_all(dst)?;
    for entry in fs::read_dir(src)? {
        let entry = entry?;
        let name = entry.file_name();
        let from = entry.path();
        let to = dst.join(&name);
        if from.is_dir() {
            copy_dir_recursive(&from, &to)?;
        } else if name != skip_top {
            fs::copy(&from, &to)?;
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
            fs::copy(&from, &to)?;
        }
    }
    Ok(())
}
