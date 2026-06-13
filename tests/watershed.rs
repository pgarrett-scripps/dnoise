//! Behavioral tests for the watershed centroider via the public API.

use dnoise::{WatershedParams, watershed_centroid};

fn params() -> WatershedParams {
    WatershedParams::default()
}

#[test]
fn dense_blob_collapses_to_one_centroid() {
    // A 5x5 block of raw points around (scan 50, tof 1000) is one ion: the
    // watershed should collapse it to a single centroid near the centre.
    let mut pts = Vec::new();
    for s in 48..=52u32 {
        for t in 998..=1002u32 {
            pts.push((s, t, 100u32));
        }
    }
    let out = watershed_centroid(&pts, &params(), usize::MAX);
    assert_eq!(out.len(), 1, "uniform blob -> one centroid");
    let (scan, tof, intensity) = out[0];
    assert_eq!(scan, 50, "centroid sits at the blob's scan centre");
    assert_eq!(tof, 1000, "centroid sits at the blob's tof centre");
    assert_eq!(intensity, 25 * 100, "intensity is the group sum");
}

#[test]
fn many_points_reduce_to_a_small_fraction() {
    // A handful of well-separated ion blobs scattered across the frame plus a
    // dusting of isolated single points: the blobs centroid, and with a seed
    // floor the lone points drop, leaving far fewer points than went in.
    let mut pts = Vec::new();
    // Ten dense blobs, each 9 points, spaced well outside each other's box.
    for b in 0..10u32 {
        let s0 = 20 + b * 40;
        let t0 = 500 + b * 400;
        for ds in 0..3u32 {
            for dt in 0..3u32 {
                pts.push((s0 + ds, t0 + dt, 200));
            }
        }
    }
    // Plus 500 weak, isolated noise points below the seed floor.
    for k in 0..500u32 {
        pts.push((5 + (k % 3), 100 + k * 7, 1));
    }
    let n_in = pts.len();

    let p = WatershedParams {
        min_seed_intensity: 50, // kills the weak isolated noise
        ..params()
    };
    let out = watershed_centroid(&pts, &p, usize::MAX);
    assert_eq!(out.len(), 10, "ten blobs -> ten centroids, noise dropped");
    assert!(
        (out.len() as f64) < 0.05 * n_in as f64,
        "centroiding shrinks the cloud to a tiny fraction: {} -> {}",
        n_in,
        out.len(),
    );
}

#[test]
fn intensity_is_conserved_within_a_group() {
    // The summed intensity of a single group's centroid equals the sum of its
    // members (no double counting, no loss).
    let pts = vec![(10, 100, 7), (10, 101, 11), (11, 100, 13), (11, 101, 17)];
    let out = watershed_centroid(&pts, &params(), usize::MAX);
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].2, 7 + 11 + 13 + 17);
}
