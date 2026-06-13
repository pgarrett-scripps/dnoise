//! Behavioral tests for the horizontal-halo filter.

use dnoise::HaloParams;
use dnoise::halo::horizontal_halo_keep_mask;

fn params() -> HaloParams {
    HaloParams {
        peak_fraction: 0.1,
        mz_idx_half_width: 100,
        scan_half_width: 2,
    }
}

#[test]
fn weak_left_right_halo_removed() {
    // Bright peak at (scan 10, tof 1000); a weak point one scan over at a nearby
    // (different) TOF index is left/right halo -> removed (5 < 0.1 * 10000).
    let scan = vec![10u32, 11];
    let tof = vec![1000u32, 1050];
    let intensity = vec![10_000u32, 5];
    let keep = horizontal_halo_keep_mask(&scan, &tof, &intensity, 700, &params());
    assert!(keep[0], "bright peak kept");
    assert!(!keep[1], "weak off-column halo removed");
}

#[test]
fn vertical_streak_never_used_against_a_point() {
    // A weak point directly above/below a bright one at the SAME TOF index is in
    // the same column, which is excluded from the reference -> it is kept.
    let scan = vec![10u32, 11];
    let tof = vec![1000u32, 1000];
    let intensity = vec![10_000u32, 5];
    let keep = horizontal_halo_keep_mask(&scan, &tof, &intensity, 700, &params());
    assert!(keep.iter().all(|&k| k), "same-column (vertical) point kept");
}

#[test]
fn far_in_mz_survives() {
    // A weak point beyond mz_idx_half_width (150 > 100) has no off-column
    // neighbor in its box -> kept.
    let scan = vec![10u32, 10];
    let tof = vec![1000u32, 1150];
    let intensity = vec![10_000u32, 5];
    let keep = horizontal_halo_keep_mask(&scan, &tof, &intensity, 700, &params());
    assert!(keep.iter().all(|&k| k), "out-of-box weak peak kept");
}

#[test]
fn far_in_mobility_survives() {
    // Different TOF, but beyond scan_half_width (5 > 2): outside the box -> kept.
    let scan = vec![10u32, 15];
    let tof = vec![1000u32, 1050];
    let intensity = vec![10_000u32, 5];
    let keep = horizontal_halo_keep_mask(&scan, &tof, &intensity, 700, &params());
    assert!(keep.iter().all(|&k| k), "point beyond scan window kept");
}

#[test]
fn comparable_neighbor_kept() {
    // A near off-column neighbor above the 10% reference survives.
    let scan = vec![10u32, 11];
    let tof = vec![1000u32, 1050];
    let intensity = vec![10_000u32, 8_000];
    let keep = horizontal_halo_keep_mask(&scan, &tof, &intensity, 700, &params());
    assert!(keep.iter().all(|&k| k), "comparable real ion kept");
}

#[test]
fn empty_is_empty() {
    let keep = horizontal_halo_keep_mask(&[], &[], &[], 700, &params());
    assert!(keep.is_empty());
}
