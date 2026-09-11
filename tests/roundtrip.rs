//! Round-trip the type-2 codec: `decode(encode(points))` must recover the input
//! point set (after grouping by scan / sorting by TOF). This is the gate that the
//! encoder matches the layout timsrust reads.

use dnoise::codec::{
    decode_frame_type2, encode_empty_frame_type2, encode_frame_type2, try_encode_frame_type2,
};

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

/// `encode_empty_frame_type2` is a DIFFERENT function from the one
/// `roundtrip_empty_frame` above exercises, and until cargo-mutants pointed it
/// out nothing called it: replacing its whole body with `vec![]`, `vec![0]` or
/// `vec![1]` left the suite green. It is no longer on the write path (see
/// `empty_frames_are_written_the_long_way` below) -- it documents the shape a
/// READER must tolerate, and its output cannot round-trip through timsrust,
/// which zstd-decodes an absent payload and fails. Assert the exact bytes the
/// format specifies.
#[test]
fn empty_frame_record_is_the_canonical_eight_bytes() {
    let record = encode_empty_frame_type2(512);
    assert_eq!(record.len(), 8, "an empty record is header-only");
    assert_eq!(
        u32::from_le_bytes(record[0..4].try_into().unwrap()),
        8,
        "total_byte_count must be the 8-byte header itself"
    );
    assert_eq!(
        u32::from_le_bytes(record[4..8].try_into().unwrap()),
        512,
        "scan count must survive verbatim"
    );
}

/// REGRESSION. An empty frame must be written the LONG way -- header plus a
/// zstd payload encoding an all-zero scan table -- not as the header-only
/// 8-byte record. timsrust slices everything after the header and zstd-decodes
/// it unconditionally, so a payload-less frame errors with "Decompression
/// fails", and one unreadable frame makes the entire `.d` unreadable to every
/// downstream consumer. dnoise's own reader tolerates both shapes, which is
/// exactly why this needs a test: nothing in this repo notices the difference.
///
/// MS/MS denoising empties frames routinely (354 of them across the 18-file
/// dda_5min benchmark arm), so this is the common case, not an edge case.
#[test]
fn empty_frames_are_written_the_long_way() {
    for num_scans in [1, 16, 512, 936] {
        let record = try_encode_frame_type2(num_scans, &[]).expect("encode");
        assert!(
            record.len() > 8,
            "num_scans={num_scans}: empty frame was written header-only, which \
             timsrust cannot decode"
        );
        assert_eq!(
            &record[8..12],
            &[0x28, 0xb5, 0x2f, 0xfd],
            "num_scans={num_scans}: payload must be a real zstd frame"
        );
        let (scans, pts) = decode_frame_type2(&record).expect("decode");
        assert_eq!(scans, num_scans);
        assert!(pts.is_empty());
    }
}

/// The one shape with no long form: zero scans means there is no scan table to
/// compress, so the header-only record is all there is.
#[test]
fn zero_scan_frame_is_header_only() {
    assert_eq!(encode_frame_type2(0, &[]).len(), 8);
}

/// The byte-transpose splits each u32 across four planes, but every value in the
/// other tests fits in three bytes (the largest is a TOF of 4,000,000, under
/// 2^24). That left the top plane always zero, so mutating `>> 24` to `<< 24` --
/// on either the encode or the decode side -- changed nothing observable.
#[test]
fn roundtrip_values_spanning_all_four_bytes() {
    let big = u32::MAX / 2; // 0x7FFF_FFFF: every byte plane non-zero
    let points = vec![(0, 1, big), (1, big, 1), (2, 0x0100_0000, 0x00FF_0000)];
    let record = encode_frame_type2(16, &points);
    let (scans, decoded) = decode_frame_type2(&record).expect("decode");
    assert_eq!(scans, 16);
    assert_eq!(canonical(&decoded), canonical(&points));
}

/// The header-validation branches. Nothing fed `decode_frame_type2` a malformed
/// record, so every comparison guarding these two errors could be mutated freely
/// (`<` to `<=`, `>` to `<`, `||` to `&&`) without a test noticing. They are the
/// checks that stand between a truncated or corrupt `tdf_bin` and an out-of-
/// bounds slice, which makes them worth pinning even though the happy path is
/// what the paper measures.
#[test]
fn malformed_records_are_rejected() {
    use dnoise::error::DecodeError;

    let short = [0u8; 7];
    assert!(
        matches!(decode_frame_type2(&short), Err(DecodeError::ShortRecord)),
        "a record below the 8-byte header must be ShortRecord"
    );

    assert_eq!(
        decode_frame_type2(&encode_empty_frame_type2(4)).unwrap(),
        (4, vec![])
    );

    // Declares more bytes than the record actually holds.
    let mut overlong = encode_empty_frame_type2(4);
    overlong[0..4].copy_from_slice(&99u32.to_le_bytes());
    assert!(
        matches!(
            decode_frame_type2(&overlong),
            Err(DecodeError::InvalidByteCount(99))
        ),
        "total_byte_count past the end of the record must be rejected"
    );

    // Declares fewer bytes than the fixed header.
    let mut undersized = encode_empty_frame_type2(4);
    undersized[0..4].copy_from_slice(&4u32.to_le_bytes());
    assert!(
        matches!(
            decode_frame_type2(&undersized),
            Err(DecodeError::InvalidByteCount(4))
        ),
        "total_byte_count below the header size must be rejected"
    );
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

#[test]
fn checked_encoder_rejects_bad_coordinates_and_handles_zero_scans() {
    use dnoise::codec::try_encode_frame_type2;
    assert!(try_encode_frame_type2(1, &[(1, 2, 3)]).is_err());
    assert!(try_encode_frame_type2(1, &[(0, u32::MAX, 3)]).is_err());
    assert_eq!(
        decode_frame_type2(&try_encode_frame_type2(0, &[]).unwrap()).unwrap(),
        (0, vec![])
    );
}

#[test]
fn corrupt_headers_and_random_payloads_never_panic() {
    let mut seed = 42u64;
    for size in 8usize..256 {
        let mut bytes = vec![0; size];
        for byte in &mut bytes {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            *byte = (seed >> 32) as u8;
        }
        bytes[..4].copy_from_slice(&(size as u32).to_le_bytes());
        let _ = decode_frame_type2(&bytes);
    }
    let mut frame = encode_frame_type2(3, &[(1, 20, 7)]);
    frame[4..8].copy_from_slice(&2u32.to_le_bytes());
    assert!(decode_frame_type2(&frame).is_err());
}
