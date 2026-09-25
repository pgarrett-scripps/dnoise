//! Property tests for the MS1 selection-polygon gate ([`dnoise::PolygonGate`]):
//! containment, full padding (Minkowski sum with the pad rectangle) and
//! pad-0 equality with the plain polygon, each checked on a dense
//! `(scan, TOF index)` grid for synthetic shapes that stress the scan-line
//! logic (concave, thin spikes, notches, axis-aligned edges, holes).

#[path = "common/polygon_ref.rs"]
mod polygon_ref;

use dnoise::PolygonGate;
use polygon_ref::{Grid, Violations, check_gate};
use std::sync::OnceLock;

// Grid: 171 scans from 1/K0 1.45 down to 0.60 (step 0.005, decreasing like a
// real timsTOF), and a sqrt-law TOF axis (~0.2-0.4 Th per index over 100-1750).
const NUM_SCANS: u32 = 171;
fn im_at_scan(s: u32) -> f64 {
    1.45 - 0.005 * s as f64
}
fn mz_to_tof(mz: f64) -> f64 {
    (mz.sqrt() - 10.0) * 150.0
}
fn tof_to_mz(t: f64) -> f64 {
    let r = t / 150.0 + 10.0;
    r * r
}

/// `(m/z pad Da, 1/K0 pad)` values checked for every shape.
const PADS: [(f64, f64); 6] = [
    (0.0, 0.0),
    (2.0, 0.0),
    (0.0, 0.015),
    (2.0, 0.01),
    (3.0, 0.02),
    (5.0, 0.05),
];

type Shape = (&'static str, Vec<f64>, Vec<f64>);

fn shapes() -> Vec<Shape> {
    // (name, m/z vertices, 1/K0 vertices)
    let mut v: Vec<Shape> = vec![
        (
            "convex_hexagon",
            vec![400.0, 900.0, 1400.0, 1500.0, 1000.0, 500.0],
            vec![0.80, 0.70, 0.90, 1.20, 1.35, 1.10],
        ),
        (
            // U: notch cut down from the top between m/z 700 and 1100.
            "concave_u",
            vec![300.0, 1500.0, 1500.0, 1100.0, 1100.0, 700.0, 700.0, 300.0],
            vec![0.70, 0.70, 1.40, 1.40, 0.95, 0.95, 1.40, 1.40],
        ),
        (
            // Body plus a sideways spike to m/z 1500 only 0.003 1/K0 thick,
            // lying wholly between the scan lines at 1.000 and 1.005.
            "thin_spike_mz",
            vec![300.0, 700.0, 700.0, 1500.0, 700.0, 700.0, 300.0],
            vec![0.80, 0.80, 1.001, 1.0025, 1.004, 1.20, 1.20],
        ),
        (
            // Body plus an upward spike 0.5 Th wide (about 2 TOF indices)
            // reaching 1/K0 1.40, tip between scan lines.
            "thin_spike_im",
            vec![300.0, 1300.0, 1300.0, 900.5, 900.25, 900.0, 300.0],
            vec![0.70, 0.70, 0.90, 0.90, 1.4023, 0.90, 0.90],
        ),
        (
            // V-notch cut in from the right, apex at m/z 900.
            "v_notch",
            vec![300.0, 1500.0, 1500.0, 900.0, 1500.0, 1500.0, 300.0],
            vec![0.70, 0.70, 1.00, 1.051, 1.10, 1.40, 1.40],
        ),
        (
            // Staircase: only vertical and horizontal edges, some horizontal
            // edges exactly on a scan line (1.10 = scan 70, 0.90 = scan 110).
            "staircase",
            vec![300.0, 1500.0, 1500.0, 1100.0, 1100.0, 700.0, 700.0, 300.0],
            vec![
                0.70,
                0.70,
                im_at_scan(110),
                im_at_scan(110),
                im_at_scan(70),
                im_at_scan(70),
                1.30,
                1.30,
            ],
        ),
        (
            "triangle",
            vec![200.0, 1600.0, 700.0],
            vec![0.65, 0.90, 1.40],
        ),
        (
            // Realistic ddaPASEF charge band: diagonal strip from low m/z / low
            // 1/K0 to high m/z / high 1/K0, cutting off singly charged ions.
            "charge_band",
            vec![250.0, 450.0, 1700.0, 1700.0, 1350.0, 250.0],
            vec![0.62, 0.62, 1.12, 1.45, 1.45, 0.80],
        ),
    ];
    // Two rings in one vertex list, each closed by repeating its first vertex
    // (a square with a square hole). The even-odd rule treats the zero-width
    // bridge between the rings as cancelling.
    v.push((
        "square_with_hole",
        vec![
            300.0, 1500.0, 1500.0, 300.0, 300.0, 700.0, 1100.0, 1100.0, 700.0, 700.0,
        ],
        vec![0.70, 0.70, 1.40, 1.40, 0.70, 0.90, 0.90, 1.20, 1.20, 0.90],
    ));
    v
}

struct CaseResult {
    shape: &'static str,
    mz_pad: f64,
    im_pad: f64,
    v: Violations,
}

/// Every shape x pad, computed once and shared by the tests below.
fn results() -> &'static Vec<CaseResult> {
    static R: OnceLock<Vec<CaseResult>> = OnceLock::new();
    R.get_or_init(|| {
        let grid = Grid {
            num_scans: NUM_SCANS,
            im_at_scan: &im_at_scan,
            tof_to_mz: &tof_to_mz,
            tof_lo: 0,
            tof_hi: mz_to_tof(1750.0) as u32,
            tof_step: 1,
        };
        let mut out = Vec::new();
        for (name, mz, im) in shapes() {
            for &(mz_pad, im_pad) in &PADS {
                let gate = PolygonGate::build(
                    &mz,
                    &im,
                    NUM_SCANS as usize,
                    im_at_scan,
                    mz_to_tof,
                    mz_pad,
                    im_pad,
                )
                .expect("gate");
                let v = check_gate(&gate, &mz, &im, &grid, mz_pad, im_pad);
                out.push(CaseResult {
                    shape: name,
                    mz_pad,
                    im_pad,
                    v,
                });
            }
        }
        out
    })
}

fn report(pick: impl Fn(&Violations) -> usize, what: &str) {
    let mut failed = Vec::new();
    for r in results() {
        let n = pick(&r.v);
        if n > 0 {
            failed.push(format!(
                "{:>18} pad ({:.1} Th, {:.3} 1/K0): {n} of {} points; e.g. {:?}",
                r.shape, r.mz_pad, r.im_pad, r.v.checked, r.v.examples
            ));
        }
    }
    assert!(failed.is_empty(), "{what} violated:\n{}", failed.join("\n"));
}

#[test]
fn containment_every_shape_every_pad() {
    report(|v| v.containment, "containment");
}

#[test]
fn full_padding_every_shape_every_pad() {
    report(|v| v.full_padding, "full padding (Minkowski sum)");
}

#[test]
fn no_over_inclusion_every_shape_every_pad() {
    report(|v| v.over_inclusion, "over-inclusion bound");
}

/// Pad 0 equals the plain polygon to within one scan / TOF index: containment
/// (inside => kept) plus the over-inclusion bound (kept => within one step).
#[test]
fn pad_zero_equals_polygon() {
    for r in results()
        .iter()
        .filter(|r| r.mz_pad == 0.0 && r.im_pad == 0.0)
    {
        assert_eq!(
            r.v.containment + r.v.over_inclusion,
            0,
            "{}: pad-0 gate differs from polygon: {:?}",
            r.shape,
            r.v.examples
        );
    }
}

/// Padding only ever adds points: every pad-0 kept point stays kept.
#[test]
fn padded_gate_is_superset_of_unpadded() {
    let grid_hi = mz_to_tof(1750.0) as u32;
    for (name, mz, im) in shapes() {
        let base = PolygonGate::build(
            &mz,
            &im,
            NUM_SCANS as usize,
            im_at_scan,
            mz_to_tof,
            0.0,
            0.0,
        )
        .unwrap();
        for &(mz_pad, im_pad) in &PADS[1..] {
            let g = PolygonGate::build(
                &mz,
                &im,
                NUM_SCANS as usize,
                im_at_scan,
                mz_to_tof,
                mz_pad,
                im_pad,
            )
            .unwrap();
            for s in 0..NUM_SCANS {
                for t in 0..=grid_hi {
                    if base.contains(s, t) {
                        assert!(
                            g.contains(s, t),
                            "{name} pad ({mz_pad},{im_pad}) lost s{s} t{t}"
                        );
                    }
                }
            }
        }
    }
}
