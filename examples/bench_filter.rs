//! Dev/maintainer tool, not a usage example for the library API.
//!
//! CPU-bound profiling harness for the per-frame denoise work.
//!
//! Loads real MS1 frames from a `.d` folder once (untimed), then loops the
//! per-frame CPU work — `filter_iterated` + `encode_frame_type2` — many times so
//! a sampling profiler (cargo flamegraph) gets dense samples on the actual hot
//! path, without file I/O dominating the profile.
//!
//! Usage:
//!   cargo flamegraph --profile profiling --example bench_filter -- <input.d> [max_frames] [reps]

use dnoise::FilterParams;
use dnoise::codec::encode_frame_type2;
use dnoise::filter::filter_iterated;
use dnoise::frame::FlatFrame;
use rusqlite::Connection;
use std::path::PathBuf;
use std::time::Instant;
use timsrust::readers::FrameReader;

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let input = PathBuf::from(
        args.next()
            .expect("usage: bench_filter <input.d> [max_frames] [reps]"),
    );
    let max_frames: usize = args.next().map(|s| s.parse().unwrap()).unwrap_or(1000);
    let reps: usize = args.next().map(|s| s.parse().unwrap()).unwrap_or(20);

    // Match benchmark/config/dnoise.toml.
    let params = FilterParams {
        num_iterations: 2,
        ..FilterParams::default()
    };

    let _ = reps;
    let reader = FrameReader::new(&input).map_err(|e| anyhow::anyhow!("open frames: {e}"))?;

    // Non-empty MS1 frame indices (0-based, == timsrust index), read straight from
    // the Frames table since the library's SQLite plumbing is crate-private.
    let conn = Connection::open(input.join("analysis.tdf"))?;
    let mut stmt = conn.prepare("SELECT NumPeaks, MsMsType FROM Frames ORDER BY Id")?;
    let ms1: Vec<usize> = stmt
        .query_map([], |r| Ok((r.get::<_, i64>(0)?, r.get::<_, i64>(1)?)))?
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .enumerate()
        .filter(|(_, (num_peaks, ms_ms_type))| *num_peaks > 0 && *ms_ms_type == 0)
        .map(|(i, _)| i)
        .take(max_frames)
        .collect();

    let mut t_read = 0.0;
    let mut t_filter = 0.0;
    let mut t_encode = 0.0;
    let mut total_points: u64 = 0;
    let mut kept_total: u64 = 0;
    let mut bytes_total: u64 = 0;

    for &i in &ms1 {
        let t = Instant::now();
        let frame = reader
            .get(i)
            .map_err(|e| anyhow::anyhow!("read frame {i}: {e}"))?;
        let flat = FlatFrame::from_frame(&frame);
        t_read += t.elapsed().as_secs_f64();

        let t = Instant::now();
        let keep = filter_iterated(&flat, &params);
        t_filter += t.elapsed().as_secs_f64();

        let t = Instant::now();
        let survivors = flat.survivors(&keep);
        let record = encode_frame_type2(flat.num_scans, &survivors);
        t_encode += t.elapsed().as_secs_f64();

        total_points += flat.len() as u64;
        kept_total += survivors.len() as u64;
        bytes_total += record.len() as u64;
    }

    let total = t_read + t_filter + t_encode;
    eprintln!(
        "{} MS1 frames, {:.2}M points, kept {kept_total}, {bytes_total} encoded bytes",
        ms1.len(),
        total_points as f64 / 1e6,
    );
    eprintln!("  read   : {t_read:.3}s ({:.1}%)", 100.0 * t_read / total);
    eprintln!(
        "  filter : {t_filter:.3}s ({:.1}%)",
        100.0 * t_filter / total
    );
    eprintln!(
        "  encode : {t_encode:.3}s ({:.1}%)",
        100.0 * t_encode / total
    );
    eprintln!("  total  : {total:.3}s");
    Ok(())
}
