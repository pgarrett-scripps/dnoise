//! Print one digest over every frame `tsr::FrameReader` decodes from each `.d`,
//! so a reader change can be checked for byte-identical output on real runs.
//!
//! The digest covers each frame's Id, scan offsets, TOF indices and
//! intensities, in Id order. Frames that decode to no points, or that the
//! reader refuses, are counted and left out of the digest: readers differ on
//! whether an empty frame is an error, and dnoise never decodes one. The
//! hasher is std's SipHash with fixed keys, so compare digests built with the
//! same toolchain.
//!
//! Usage: cargo run --release --example frame_digest -- <PATH.d>...

use dnoise::tsr::FrameReader;
use std::hash::{DefaultHasher, Hash, Hasher};

fn main() -> anyhow::Result<()> {
    for path in std::env::args().skip(1) {
        let reader = match FrameReader::new(&path) {
            Ok(r) => r,
            Err(e) => {
                println!("{path}\topen failed: {e}");
                continue;
            }
        };
        let mut h = DefaultHasher::new();
        let mut points = 0usize;
        let mut refused = Vec::new();
        let mut empty = 0usize;
        for i in 0..reader.len() {
            match reader.get(i) {
                Ok(f) if f.tof_indices.is_empty() => empty += 1,
                Ok(f) => {
                    f.index.hash(&mut h);
                    f.scan_offsets.hash(&mut h);
                    f.tof_indices.hash(&mut h);
                    f.intensities.hash(&mut h);
                    points += f.tof_indices.len();
                }
                Err(_) => refused.push(i + 1),
            }
        }
        println!(
            "{path}\tframes={}\tpoints={points}\tempty={empty}\trefused={} (first {:?})\tdigest={:016x}",
            reader.len(),
            refused.len(),
            &refused[..refused.len().min(5)],
            h.finish()
        );
    }
    Ok(())
}
