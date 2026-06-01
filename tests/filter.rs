//! Behavioral tests for the vertical-IM feature filter on synthetic frames.

use dnoise::FilterParams;
use dnoise::filter::{filter_iterated, filter_once};
use dnoise::frame::FlatFrame;

/// Build a FlatFrame from explicit `(scan, tof, intensity)` points.
fn frame(num_scans: usize, points: &[(u32, u32, u32)]) -> FlatFrame {
    FlatFrame {
        frame_id: 1,
        num_scans,
        scan: points.iter().map(|&(s, _, _)| s).collect(),
        tof: points.iter().map(|&(_, t, _)| t).collect(),
        intensity: points.iter().map(|&(_, _, i)| i).collect(),
    }
}

#[test]
fn keeps_long_vertical_streak_drops_isolated() {
    let mut pts = Vec::new();
    // Vertical streak at tof 1000 across 8 consecutive scans (>= min_feature_length 5).
    for s in 10..18u32 {
        pts.push((s, 1000, 100));
    }
    // Isolated single hits far away and short.
    pts.push((50, 2000, 100));
    pts.push((60, 2000, 100)); // only 2 scans apart -> two separate length-1 runs

    let f = frame(700, &pts);
    let params = FilterParams::default();
    let keep = filter_once(&f, &params);

    // Streak points kept (the first 8 points are the streak).
    assert!(
        keep[..8].iter().all(|&k| k),
        "all streak points should be kept"
    );
    // Isolated points dropped (run length 1 < 5).
    assert!(!keep[8]);
    assert!(!keep[9]);
}

#[test]
fn gap_closing_bridges_small_gaps() {
    // tof 1000: scans 10,11,12, gap at 13, 14,15,16 -> with max_internal_gap=1 it's one run of span 7.
    let pts: Vec<(u32, u32, u32)> = [10u32, 11, 12, 14, 15, 16]
        .iter()
        .map(|&s| (s, 1000, 100))
        .collect();
    let f = frame(700, &pts);
    let params = FilterParams {
        min_feature_length: 5,
        max_internal_gap: 1,
        ..Default::default()
    };
    let keep = filter_once(&f, &params);
    assert!(
        keep.iter().all(|&k| k),
        "all points in the bridged run should survive"
    );
}

#[test]
fn gap_too_large_splits_into_short_runs() {
    // Same points but gap of 3 empty scans (13,14,15 missing) splits into two runs of 3 -> both < 5.
    let pts: Vec<(u32, u32, u32)> = [10u32, 11, 12, 16, 17, 18]
        .iter()
        .map(|&s| (s, 1000, 100))
        .collect();
    let f = frame(700, &pts);
    let params = FilterParams {
        min_feature_length: 5,
        max_internal_gap: 1,
        ..Default::default()
    };
    let keep = filter_once(&f, &params);
    assert!(
        keep.iter().all(|&k| !k),
        "split short runs should all be dropped"
    );
}

#[test]
fn window_aggregates_neighboring_tofs() {
    // Two adjacent TOF columns (1000 and 1001) each 3 scans long: individually < 5,
    // but within mz_half_width=2 they share a window and together span >= 5 scans.
    let mut pts = Vec::new();
    for s in 10..13u32 {
        pts.push((s, 1000, 100));
    }
    for s in 13..16u32 {
        pts.push((s, 1001, 100));
    }
    let f = frame(700, &pts);
    let params = FilterParams {
        mz_half_width: 2,
        min_feature_length: 5,
        ..Default::default()
    };
    let keep = filter_once(&f, &params);
    assert!(
        keep.iter().all(|&k| k),
        "neighboring TOFs should aggregate into one feature"
    );
}

#[test]
fn iteration_is_monotonic() {
    let mut pts = Vec::new();
    for s in 10..18u32 {
        pts.push((s, 1000, 100));
    }
    pts.push((50, 2000, 100));
    let f = frame(700, &pts);

    let p1 = FilterParams {
        num_iterations: 1,
        ..Default::default()
    };
    let p3 = FilterParams {
        num_iterations: 3,
        ..Default::default()
    };
    let k1 = filter_iterated(&f, &p1);
    let k3 = filter_iterated(&f, &p3);

    let c1 = k1.iter().filter(|&&b| b).count();
    let c3 = k3.iter().filter(|&&b| b).count();
    assert!(c3 <= c1, "more iterations must not keep more points");
    // The stable streak survives any number of passes.
    assert!(k3[0..8].iter().all(|&b| b));
}
