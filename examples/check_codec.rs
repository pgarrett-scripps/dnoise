//! Dev/maintainer tool, not a usage example for the library API.
//!
//! Cross-check our type-2 codec against real Bruker bytes: read each frame's raw
//! record straight from analysis.tdf_bin, decode it with our [`decode_frame_type2`],
//! and confirm the points match what timsrust returns for the same frame.
//!
//! This proves the codec matches Bruker's on-disk format, not merely itself.
//!
//! Usage: cargo run --release --example check_codec -- <PATH.d> [num_frames]

use dnoise::codec::decode_frame_type2;
use dnoise::frame::FlatFrame;
use rusqlite::Connection;
use std::collections::HashMap;
use std::io::{Read, Seek, SeekFrom};
use std::path::PathBuf;
use timsrust::readers::FrameReader;

fn canonical(mut v: Vec<(u32, u32, u32)>) -> Vec<(u32, u32, u32)> {
    v.sort_unstable_by_key(|&(s, t, _)| (s, t));
    v
}

fn main() -> anyhow::Result<()> {
    let path: PathBuf = std::env::args()
        .nth(1)
        .expect("usage: check_codec <PATH.d> [n]")
        .into();
    let limit: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(20);

    let conn = Connection::open(path.join("analysis.tdf"))?;
    let mut stmt = conn.prepare("SELECT Id, TimsId FROM Frames")?;
    let tims_id: HashMap<usize, u64> = stmt
        .query_map([], |r| {
            Ok((r.get::<_, i64>(0)? as usize, r.get::<_, i64>(1)? as u64))
        })?
        .collect::<Result<_, _>>()?;

    let reader = FrameReader::new(&path).map_err(|e| anyhow::anyhow!("{e}"))?;
    let mut bin = std::fs::File::open(path.join("analysis.tdf_bin"))?;

    let n = reader.len().min(limit);
    let mut ok = 0usize;
    for i in 0..n {
        let frame = reader
            .get(i)
            .map_err(|e| anyhow::anyhow!("frame {i}: {e}"))?;
        let offset = tims_id[&frame.index];

        // Read the raw record: u32 byte_count at the offset, then byte_count bytes.
        bin.seek(SeekFrom::Start(offset))?;
        let mut len_buf = [0u8; 4];
        bin.read_exact(&mut len_buf)?;
        let byte_count = u32::from_le_bytes(len_buf) as usize;
        let mut record = vec![0u8; byte_count];
        bin.seek(SeekFrom::Start(offset))?;
        bin.read_exact(&mut record)?;

        let (scan_count, ours) = decode_frame_type2(&record)?;
        let flat = FlatFrame::from_frame(&frame);
        let theirs: Vec<(u32, u32, u32)> = (0..flat.len())
            .map(|p| (flat.scan[p], flat.tof[p], flat.intensity[p]))
            .collect();

        assert_eq!(scan_count, flat.num_scans, "frame {i}: scan count mismatch");
        assert_eq!(
            canonical(ours),
            canonical(theirs),
            "frame {i}: point mismatch"
        );
        ok += 1;
    }
    println!("check_codec: {ok}/{n} frames — our codec matches Bruker bytes exactly");
    Ok(())
}
