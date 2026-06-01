//! Round-trip the type-2 codec: `decode(encode(points))` must recover the input
//! point set (after grouping by scan / sorting by TOF). This is the gate that the
//! encoder matches the layout timsrust reads.

use dnoise::tdf::encode::{decode_frame_type2, encode_frame_type2};

/// Group/sort points the way the codec canonicalizes them, for comparison.
fn canonical(points: &[(u32, u32, u32)]) -> Vec<(u32, u32, u32)> {
    let mut v = points.to_vec();
    v.sort_unstable_by_key(|&(scan, tof, _)| (scan, tof));
    v
}

#[test]
fn roundtrip_typical_frame() {
    let num_scans = 700;
    let mut points = Vec::new();
    // A vertical streak at tof 12345 across scans 100..130.
    for s in 100..130u32 {
        points.push((s, 12345, 50 + s));
    }
    // Some scattered peaks with multiple TOFs per scan (ascending and not).
    points.push((5, 999, 7));
    points.push((5, 200, 3));
    points.push((5, 50000, 11));
    points.push((699, 1, 1)); // last scan, smallest tof
    points.push((0, 4_000_000, 42)); // large tof exercises the high byte

    let record = encode_frame_type2(num_scans, &points);
    let (decoded_scans, decoded) = decode_frame_type2(&record).expect("decode");

    assert_eq!(decoded_scans, num_scans);
    assert_eq!(canonical(&decoded), canonical(&points));
}

#[test]
fn roundtrip_empty_frame() {
    let record = encode_frame_type2(512, &[]);
    let (scans, pts) = decode_frame_type2(&record).expect("decode");
    assert_eq!(scans, 512);
    assert!(pts.is_empty());
}

#[test]
fn roundtrip_single_point() {
    let record = encode_frame_type2(100, &[(42, 7, 9)]);
    let (scans, pts) = decode_frame_type2(&record).expect("decode");
    assert_eq!(scans, 100);
    assert_eq!(pts, vec![(42, 7, 9)]);
}

#[test]
fn record_length_matches_header() {
    let record = encode_frame_type2(300, &[(1, 10, 5), (1, 20, 6), (2, 30, 7)]);
    let total = u32::from_le_bytes(record[0..4].try_into().unwrap()) as usize;
    assert_eq!(
        total,
        record.len(),
        "header total_byte_count must equal record length"
    );
    let scan_count = u32::from_le_bytes(record[4..8].try_into().unwrap());
    assert_eq!(scan_count, 300);
}
