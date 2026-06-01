//! Compare the zstd frame header Bruker writes vs what our two encode options write.
//! Usage: cargo run --release --example zstd_probe -- <PATH.d>

use std::io::{Read, Seek, SeekFrom};
use std::path::PathBuf;

fn describe(tag: &str, frame: &[u8]) {
    let magic = &frame[0..4];
    let desc = frame[4];
    let fcs_flag = (desc >> 6) & 3;
    let single_segment = (desc >> 5) & 1;
    let checksum = (desc >> 2) & 1;
    let dict_flag = desc & 3;
    println!(
        "{tag}: len={} magic={:02x?} desc=0x{desc:02x} FCS_flag={fcs_flag} single_segment={single_segment} checksum={checksum} dict={dict_flag} head={:02x?}",
        frame.len(),
        magic,
        &frame[..frame.len().min(10)]
    );
}

fn main() -> anyhow::Result<()> {
    let path: PathBuf = std::env::args()
        .nth(1)
        .expect("usage: zstd_probe <PATH.d>")
        .into();

    // First frame offset = MIN(TimsId).
    let conn = rusqlite::Connection::open(path.join("analysis.tdf"))?;
    let off: i64 = conn.query_row("SELECT MIN(TimsId) FROM Frames", [], |r| r.get(0))?;

    let mut f = std::fs::File::open(path.join("analysis.tdf_bin"))?;
    f.seek(SeekFrom::Start(off as u64))?;
    let mut hdr = [0u8; 8];
    f.read_exact(&mut hdr)?;
    let bin_size = u32::from_le_bytes(hdr[0..4].try_into().unwrap()) as usize;
    let mut payload = vec![0u8; bin_size - 8];
    f.read_exact(&mut payload)?;

    describe("bruker  ", &payload);

    let raw = zstd::decode_all(&payload[..])?;
    println!("decompressed transposed-blob len = {}", raw.len());

    let streamed = zstd::encode_all(&raw[..], 1)?;
    describe("encode_all", &streamed);

    let bulk = zstd::bulk::compress(&raw, 1)?;
    describe("bulk::comp", &bulk);

    Ok(())
}
