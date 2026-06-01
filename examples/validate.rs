//! Validate a (denoised) .d folder: re-read every frame with timsrust and confirm
//! the per-frame peak count matches `Frames.NumPeaks` in analysis.tdf.
//!
//! Empty frames (NumPeaks == 0) are expected: timsrust cannot decode their absent
//! payload, so a read error on a frame the DB marks empty is treated as OK.
//!
//! Usage: cargo run --example validate -- <PATH.d>

use rusqlite::Connection;
use std::path::PathBuf;
use timsrust::readers::FrameReader;

fn main() -> anyhow::Result<()> {
    let path: PathBuf = std::env::args()
        .nth(1)
        .expect("usage: validate <PATH.d>")
        .into();

    // (Id, NumPeaks) ordered by Id == timsrust frame index.
    let conn = Connection::open(path.join("analysis.tdf"))?;
    let mut stmt = conn.prepare("SELECT Id, NumPeaks FROM Frames ORDER BY Id")?;
    let meta: Vec<(usize, i64)> = stmt
        .query_map([], |r| {
            Ok((r.get::<_, i64>(0)? as usize, r.get::<_, i64>(1)?))
        })?
        .collect::<Result<_, _>>()?;
    let comp: String = conn.query_row(
        "SELECT Value FROM GlobalMetadata WHERE Key='TimsCompressionType'",
        [],
        |r| r.get(0),
    )?;

    let reader = FrameReader::new(&path).map_err(|e| anyhow::anyhow!("{e}"))?;
    let n = reader.len();

    let mut total = 0u64;
    let mut empty = 0u64;
    let mut mismatches = 0u64;
    for i in 0..n {
        let (id, db_peaks) = meta.get(i).copied().unwrap_or((0, -1));
        match reader.get(i) {
            Ok(frame) => {
                let peaks = frame.tof_indices.len();
                total += peaks as u64;
                if db_peaks >= 0 && db_peaks as usize != peaks {
                    mismatches += 1;
                    if mismatches <= 5 {
                        eprintln!("frame {i} (id {id}): binary {peaks} != db {db_peaks}");
                    }
                }
            }
            Err(e) => {
                // A read failure is acceptable only for a genuinely empty frame.
                if db_peaks == 0 {
                    empty += 1;
                } else {
                    return Err(anyhow::anyhow!(
                        "frame {i} (id {id}, db peaks {db_peaks}): {e}"
                    ));
                }
            }
        }
    }
    println!(
        "validate: {n} frames re-read OK ({empty} empty), {total} total points, compression={comp}, {mismatches} db/bin mismatches"
    );
    if mismatches > 0 {
        std::process::exit(1);
    }
    Ok(())
}
