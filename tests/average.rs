//! Behavioral tests for the running-average pre-filter smoothing pass.

use dnoise::average::running_average;
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

/// Look up the averaged intensity at a `(scan, tof)` bin, if present.
fn at(f: &FlatFrame, scan: u32, tof: u32) -> Option<u32> {
    (0..f.len())
        .find(|&i| f.scan[i] == scan && f.tof[i] == tof)
        .map(|i| f.intensity[i])
}

#[test]
fn recurring_peak_keeps_native_intensity() {
    // The same (scan, tof) bin at intensity 100 in all three windowed frames
    // averages back to 100 — averaging preserves a persistent ion's scale.
    let f0 = frame(700, &[(10, 1000, 100)]);
    let f1 = frame(700, &[(10, 1000, 100)]);
    let f2 = frame(700, &[(10, 1000, 100)]);
    let avg = running_average(1, 700, &[&f0, &f1, &f2]);
    assert_eq!(at(&avg, 10, 1000), Some(100));
}

#[test]
fn single_frame_noise_is_divided_down() {
    // A bin present in only one of three frames is divided by 3 (300 -> 100),
    // while bins present in none stay absent.
    let f0 = frame(700, &[(20, 2000, 300)]);
    let f1 = frame(700, &[]);
    let f2 = frame(700, &[]);
    let avg = running_average(1, 700, &[&f0, &f1, &f2]);
    assert_eq!(at(&avg, 20, 2000), Some(100));
    assert_eq!(at(&avg, 21, 2000), None);
}

#[test]
fn rounds_to_zero_drops_bin() {
    // A lone count of 1 spread across 3 frames rounds to 0 (1/3 -> 0) and is dropped.
    let f0 = frame(700, &[(5, 500, 1)]);
    let f1 = frame(700, &[]);
    let f2 = frame(700, &[]);
    let avg = running_average(1, 700, &[&f0, &f1, &f2]);
    assert_eq!(at(&avg, 5, 500), None);
    assert!(avg.is_empty());
}

#[test]
fn clamped_edge_window_divides_by_actual_count() {
    // At a run edge the window is smaller; dividing by the actual frame count
    // (2 here) keeps intensities on the native scale rather than under-counting.
    let f0 = frame(700, &[(10, 1000, 100)]);
    let f1 = frame(700, &[(10, 1000, 100)]);
    let avg = running_average(1, 700, &[&f0, &f1]);
    assert_eq!(at(&avg, 10, 1000), Some(100));
}

#[test]
fn distinct_bins_are_preserved() {
    // Non-overlapping bins across frames all survive at their divided intensities.
    let f0 = frame(700, &[(10, 1000, 90)]);
    let f1 = frame(700, &[(11, 1001, 90), (10, 1000, 90)]);
    let avg = running_average(1, 700, &[&f0, &f1]);
    // (10,1000): (90+90)/2 = 90; (11,1001): 90/2 = 45.
    assert_eq!(at(&avg, 10, 1000), Some(90));
    assert_eq!(at(&avg, 11, 1001), Some(45));
}
