//! Search-free precursor MS1 areas, before and after denoising.
//!
//! Every ddaPASEF run carries the instrument's own precursor picks in the
//! `Precursors` table of `analysis.tdf` (monoisotopic m/z, charge, mobility scan,
//! parent MS1 frame). This example seeds from that table and integrates the same
//! window in two arms of the same run -- the original `.d` and a dnoise output --
//! so the two areas differ only by what the filter removed. No search engine is
//! involved, which is what lets the paper's cross-platform SI section run on any
//! public ddaPASEF deposit.
//!
//! `--mz-offset` shifts the whole isotope ladder by a fixed number of Da, which
//! turns the same seed list into a decoy control: a window of identical shape in
//! the same frames and mobility scans, at an m/z where this precursor has no
//! isotope. Comparing decoy to target retention separates "the filter removed the
//! peptide" from "the filter removed background that fell inside the box".
//!
//! The window is the monoisotopic peak plus the next `--isotopes - 1` isotopes,
//! each `--ppm` wide, over scans `ScanNumber +/- --scan-pad`, over the parent MS1
//! frame and its `--frame-pad` MS1 neighbours on each side. Both arms share the
//! calibration tables (dnoise copies them), so a window computed once in TOF-index
//! and scan space is identical in both.
//!
//! Usage:
//!   cargo run --release --example precursor_area -- \
//!       <original.d> <denoised.d> <out.csv> [--ppm 20] [--scan-pad 3] [--frame-pad 2] \
//!       [--isotopes 3] [--mz-offset 0.0]
//!
//! Output columns: precursor_id, mz, charge, scan, parent_frame, instrument_intensity,
//! area_original, area_denoised, points_original, points_denoised.

use dnoise::tsr::ConvertableDomain;
use dnoise::tsr::{FrameReader, MetadataReader};
use rayon::prelude::*;
use rusqlite::Connection;
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

const ISOTOPE_SPACING: f64 = 1.003_355;

struct Precursor {
    id: i64,
    mz: f64,
    charge: i64,
    scan: i64,
    parent: i64,
    intensity: f64,
}

struct Opts {
    ppm: f64,
    mz_offset: f64,
    scan_pad: i64,
    frame_pad: usize,
    isotopes: usize,
}

fn read_precursors(tdf: &Path) -> Result<Vec<Precursor>, Box<dyn std::error::Error>> {
    let conn = Connection::open(tdf)?;
    let mut stmt = conn.prepare(
        "SELECT Id, COALESCE(MonoisotopicMz, LargestPeakMz), COALESCE(Charge, 0), \
         ScanNumber, Parent, COALESCE(Intensity, 0) FROM Precursors ORDER BY Id",
    )?;
    let rows = stmt.query_map([], |r| {
        Ok(Precursor {
            id: r.get(0)?,
            mz: r.get(1)?,
            charge: r.get(2)?,
            scan: r.get::<_, f64>(3)?.round() as i64,
            parent: r.get(4)?,
            intensity: r.get(5)?,
        })
    })?;
    Ok(rows.collect::<Result<Vec<_>, _>>()?)
}

/// MS1 frame Ids in acquisition order, and the check that Frames.Id is the
/// contiguous 1-based sequence timsrust indexes by (index = Id - 1).
fn read_ms1_frames(tdf: &Path) -> Result<Vec<i64>, Box<dyn std::error::Error>> {
    let conn = Connection::open(tdf)?;
    let (lo, hi, n): (i64, i64, i64) =
        conn.query_row("SELECT MIN(Id), MAX(Id), COUNT(*) FROM Frames", [], |r| {
            Ok((r.get(0)?, r.get(1)?, r.get(2)?))
        })?;
    if lo != 1 || hi != n {
        return Err(format!("Frames.Id is not contiguous 1..{n} (min {lo}, max {hi})").into());
    }
    let mut stmt = conn.prepare("SELECT Id FROM Frames WHERE MsMsType = 0 ORDER BY Id")?;
    let ids = stmt.query_map([], |r| r.get::<_, i64>(0))?;
    Ok(ids.collect::<Result<Vec<_>, _>>()?)
}

/// Per precursor: the TOF-index ranges (one per isotope) and scan range to sum.
struct Window {
    tof: Vec<(u32, u32)>,
    scan_lo: usize,
    scan_hi: usize,
}

fn windows(precursors: &[Precursor], meta: &dnoise::tsr::Metadata, o: &Opts) -> Vec<Window> {
    precursors
        .iter()
        .map(|p| {
            let z = if p.charge > 0 { p.charge as f64 } else { 1.0 };
            let tof = (0..o.isotopes)
                .map(|k| {
                    let mz = p.mz + o.mz_offset + k as f64 * ISOTOPE_SPACING / z;
                    let lo = meta.mz_converter.invert(mz * (1.0 - o.ppm * 1e-6));
                    let hi = meta.mz_converter.invert(mz * (1.0 + o.ppm * 1e-6));
                    (lo.floor().max(0.0) as u32, hi.ceil().max(0.0) as u32)
                })
                .collect();
            Window {
                tof,
                scan_lo: (p.scan - o.scan_pad).max(0) as usize,
                scan_hi: (p.scan + o.scan_pad).max(0) as usize,
            }
        })
        .collect()
}

/// Sum every point of `frame` that falls in a precursor's window, for each
/// precursor listed against this frame. Returns (precursor index, area, points).
fn integrate(
    frame: &dnoise::tsr::Frame,
    wanted: &[usize],
    win: &[Window],
) -> Vec<(usize, f64, u64)> {
    let n_scans = frame.scan_offsets.len().saturating_sub(1);
    wanted
        .iter()
        .map(|&pi| {
            let w = &win[pi];
            let (mut area, mut points) = (0.0f64, 0u64);
            let hi = w.scan_hi.min(n_scans.saturating_sub(1));
            for scan in w.scan_lo..=hi {
                let (a, b) = (frame.scan_offsets[scan], frame.scan_offsets[scan + 1]);
                for i in a..b {
                    let t = frame.tof_indices[i];
                    if w.tof.iter().any(|&(lo, hi)| t >= lo && t <= hi) {
                        area += frame.intensities[i] as f64;
                        points += 1;
                    }
                }
            }
            (pi, area, points)
        })
        .collect()
}

fn areas(
    d: &Path,
    by_frame: &HashMap<usize, Vec<usize>>,
    win: &[Window],
    n: usize,
) -> Result<(Vec<f64>, Vec<u64>), Box<dyn std::error::Error>> {
    let reader = FrameReader::new(d).map_err(|e| format!("{}: {e}", d.display()))?;
    let mut frames: Vec<usize> = by_frame.keys().copied().collect();
    frames.sort_unstable();
    let partial: Vec<Vec<(usize, f64, u64)>> = frames
        .par_iter()
        .map(|&idx| match reader.get(idx) {
            Ok(frame) => integrate(&frame, &by_frame[&idx], win),
            Err(e) => {
                eprintln!("warning: frame index {idx}: {e}");
                Vec::new()
            }
        })
        .collect();
    let mut area = vec![0.0; n];
    let mut points = vec![0u64; n];
    for part in partial {
        for (pi, a, p) in part {
            area[pi] += a;
            points[pi] += p;
        }
    }
    Ok((area, points))
}

fn parse_args() -> Result<(PathBuf, PathBuf, PathBuf, Opts), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let usage = "usage: precursor_area <original.d> <denoised.d> <out.csv> \
                 [--ppm 20] [--scan-pad 3] [--frame-pad 2] [--isotopes 3] [--mz-offset 0.0]";
    if args.len() < 3 {
        return Err(usage.into());
    }
    let mut o = Opts {
        ppm: 20.0,
        mz_offset: 0.0,
        scan_pad: 3,
        frame_pad: 2,
        isotopes: 3,
    };
    let mut i = 3;
    while i < args.len() {
        let v = args
            .get(i + 1)
            .ok_or_else(|| format!("{} needs a value", args[i]))?;
        match args[i].as_str() {
            "--ppm" => o.ppm = v.parse().map_err(|_| usage)?,
            "--mz-offset" => o.mz_offset = v.parse().map_err(|_| usage)?,
            "--scan-pad" => o.scan_pad = v.parse().map_err(|_| usage)?,
            "--frame-pad" => o.frame_pad = v.parse().map_err(|_| usage)?,
            "--isotopes" => o.isotopes = v.parse().map_err(|_| usage)?,
            other => return Err(format!("unknown option {other}\n{usage}")),
        }
        i += 2;
    }
    if o.isotopes == 0 {
        return Err("--isotopes must be at least 1".into());
    }
    Ok((
        args[0].clone().into(),
        args[1].clone().into(),
        args[2].clone().into(),
        o,
    ))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (original, denoised, out, o) = parse_args()?;
    let tdf = original.join("analysis.tdf");
    let meta = MetadataReader::new(&original).map_err(|e| format!("metadata: {e}"))?;
    let precursors = read_precursors(&tdf)?;
    let ms1 = read_ms1_frames(&tdf)?;
    if precursors.is_empty() {
        return Err("Precursors table is empty: not a ddaPASEF run?".into());
    }
    // Position of each MS1 frame Id in the MS1-only sequence, so "frame-pad
    // neighbours" means neighbouring survey scans, not neighbouring frame Ids.
    let ms1_pos: HashMap<i64, usize> = ms1.iter().enumerate().map(|(i, &id)| (id, i)).collect();
    let win = windows(&precursors, &meta, &o);

    let mut by_frame: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut orphans = 0usize;
    for (pi, p) in precursors.iter().enumerate() {
        let Some(&pos) = ms1_pos.get(&p.parent) else {
            orphans += 1;
            continue;
        };
        let lo = pos.saturating_sub(o.frame_pad);
        let hi = (pos + o.frame_pad).min(ms1.len() - 1);
        for &frame in &ms1[lo..=hi] {
            by_frame.entry((frame - 1) as usize).or_default().push(pi);
        }
    }
    if orphans > 0 {
        eprintln!("warning: {orphans} precursors whose Parent is not an MS1 frame were skipped");
    }
    eprintln!(
        "{} precursors, {} MS1 frames to read per arm (ppm {}, scans +/-{}, frames +/-{}, {} isotopes)",
        precursors.len(),
        by_frame.len(),
        o.ppm,
        o.scan_pad,
        o.frame_pad,
        o.isotopes
    );

    let (a_orig, p_orig) = areas(&original, &by_frame, &win, precursors.len())?;
    let (a_den, p_den) = areas(&denoised, &by_frame, &win, precursors.len())?;

    let mut w = std::io::BufWriter::new(std::fs::File::create(&out)?);
    writeln!(
        w,
        "precursor_id,mz,charge,scan,parent_frame,instrument_intensity,\
         area_original,area_denoised,points_original,points_denoised"
    )?;
    for (i, p) in precursors.iter().enumerate() {
        writeln!(
            w,
            "{},{:.5},{},{},{},{},{},{},{},{}",
            p.id,
            p.mz,
            p.charge,
            p.scan,
            p.parent,
            p.intensity,
            a_orig[i],
            a_den[i],
            p_orig[i],
            p_den[i]
        )?;
    }
    let kept = a_den
        .iter()
        .zip(&a_orig)
        .filter(|(d, o)| **d > 0.0 && **o > 0.0)
        .count();
    let total: f64 = a_orig.iter().sum();
    let total_den: f64 = a_den.iter().sum();
    eprintln!(
        "wrote {} rows to {}: {} with area in both arms, summed area retained {:.1}%",
        precursors.len(),
        out.display(),
        kept,
        100.0 * total_den / total.max(1.0)
    );
    Ok(())
}
