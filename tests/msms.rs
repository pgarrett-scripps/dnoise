//! Behavioral tests for the MS/MS combine-and-filter core.

use dnoise::FilterParams;
use dnoise::msms::combine_and_filter;

fn key(scan: u32, tof: u32) -> u64 {
    ((scan as u64) << 32) | tof as u64
}

#[test]
fn keeps_fragment_streak_drops_isolated() {
    // Precursor window starts at scan 100, 25 scans wide. A fragment streak at
    // TOF 1000 spans scans 100..106, observed in two frames (summed). Plus one
    // isolated noise hit.
    let mut pts = Vec::new();
    for _frame in 0..2 {
        for s in 100..106u32 {
            pts.push((s, 1000, 50));
        }
    }
    pts.push((110, 2000, 40)); // isolated noise

    let params = FilterParams {
        min_feature_length: 3,
        ..Default::default()
    };
    let keep = combine_and_filter(&pts, 100, 25, &params, None);

    for s in 100..106u32 {
        assert!(
            keep.contains(&key(s, 1000)),
            "streak scan {s} kept (absolute)"
        );
    }
    assert!(!keep.contains(&key(110, 2000)), "isolated noise dropped");
}

#[test]
fn empty_points_empty_set() {
    let keep = combine_and_filter(&[], 0, 25, &FilterParams::default(), None);
    assert!(keep.is_empty());
}

#[test]
fn keys_are_absolute_scan_coordinates() {
    // A streak well inside the window maps back to absolute scans (offset by sb0).
    let mut pts = Vec::new();
    for s in 200..206u32 {
        pts.push((s, 500, 100));
    }
    let params = FilterParams {
        min_feature_length: 3,
        ..Default::default()
    };
    let keep = combine_and_filter(&pts, 195, 25, &params, None);
    assert!(keep.contains(&key(200, 500)));
    assert!(
        !keep.contains(&key(5, 500)),
        "must not be local coordinates"
    );
}
