//! Dump one frame's points as `(m/z, 1/K0, intensity)` CSV, using timsrust's
//! calibration converters. Used to render the before/after and concept figures
//! for the paper.
//!
//! Usage: cargo run --release --example dump_frame -- <PATH.d> <frame_index0> [out.csv]
//! `frame_index0` is the 0-based timsrust frame index (= Frames.Id - 1).

use dnoise::frame::FlatFrame;
use std::io::Write;
use std::path::PathBuf;
use timsrust::converters::ConvertableDomain;
use timsrust::readers::{FrameReader, MetadataReader};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let input = PathBuf::from(
        args.next()
            .expect("usage: dump_frame <PATH.d> <frame_index0> [out.csv]"),
    );
    let idx: usize = args.next().expect("frame_index0").parse()?;
    let out = args.next().unwrap_or_else(|| "frame.csv".to_string());

    let meta = MetadataReader::new(&input).map_err(|e| format!("metadata: {e}"))?;
    let reader = FrameReader::new(&input).map_err(|e| format!("frames: {e}"))?;
    let frame = reader
        .get(idx)
        .map_err(|e| format!("get frame {idx}: {e}"))?;
    let flat = FlatFrame::from_frame(&frame);

    let mut w = std::io::BufWriter::new(std::fs::File::create(&out)?);
    writeln!(w, "mz,one_over_k0,intensity,scan,tof")?;
    for i in 0..flat.len() {
        let mz = meta.mz_converter.convert(flat.tof[i]);
        let im = meta.im_converter.convert(flat.scan[i]);
        writeln!(
            w,
            "{mz:.5},{im:.5},{},{},{}",
            flat.intensity[i], flat.scan[i], flat.tof[i]
        )?;
    }
    eprintln!("wrote {} points to {out}", flat.len());
    Ok(())
}
