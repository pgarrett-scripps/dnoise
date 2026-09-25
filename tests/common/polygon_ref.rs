//! Brute-force reference geometry for checking [`dnoise::PolygonGate`] on a
//! grid of `(scan, TOF index)` points. Deliberately independent of the gate's
//! own scan-line code: point-in-polygon by ray crossing, and "does this
//! rectangle meet the polygon" by clipping every edge against the rectangle.
//!
//! Shared by `tests/polygon_props.rs` (synthetic shapes) and
//! `examples/polygon_check.rs` (real run polygons) via `#[path]`.
#![allow(dead_code)]

use dnoise::PolygonGate;

/// Even-odd point-in-polygon, half-open in y (same convention as the gate).
pub fn point_in_polygon(mz: &[f64], im: &[f64], x: f64, y: f64) -> bool {
    let n = mz.len();
    let mut inside = false;
    for i in 0..n {
        let j = (i + 1) % n;
        let (yi, yj) = (im[i], im[j]);
        if (yi > y) != (yj > y) {
            let xc = mz[i] + (y - yi) / (yj - yi) * (mz[j] - mz[i]);
            if x < xc {
                inside = !inside;
            }
        }
    }
    inside
}

/// True when segment `a`-`b` meets the closed rectangle `[x0,x1] x [y0,y1]`
/// (Liang-Barsky clip; a zero-width or zero-height rectangle is allowed).
fn segment_meets_rect(a: (f64, f64), b: (f64, f64), x0: f64, x1: f64, y0: f64, y1: f64) -> bool {
    let (dx, dy) = (b.0 - a.0, b.1 - a.1);
    let (mut t0, mut t1) = (0.0f64, 1.0f64);
    for (p, q) in [
        (-dx, a.0 - x0),
        (dx, x1 - a.0),
        (-dy, a.1 - y0),
        (dy, y1 - a.1),
    ] {
        if p == 0.0 {
            if q < 0.0 {
                return false;
            }
        } else {
            let r = q / p;
            if p < 0.0 {
                if r > t1 {
                    return false;
                }
                t0 = t0.max(r);
            } else {
                if r < t0 {
                    return false;
                }
                t1 = t1.min(r);
            }
        }
    }
    t0 <= t1
}

/// True when the closed rectangle `[x0,x1] x [y0,y1]` meets the polygon, i.e.
/// the rectangle's centre lies in the Minkowski sum of the polygon with the
/// rectangle's half-extents.
pub fn rect_meets_polygon(mz: &[f64], im: &[f64], x0: f64, x1: f64, y0: f64, y1: f64) -> bool {
    let n = mz.len();
    for i in 0..n {
        let j = (i + 1) % n;
        if segment_meets_rect((mz[i], im[i]), (mz[j], im[j]), x0, x1, y0, y1) {
            return true;
        }
    }
    // No edge touches the rectangle: it is wholly inside or wholly outside.
    point_in_polygon(mz, im, 0.5 * (x0 + x1), 0.5 * (y0 + y1))
}

/// The grid the gate is checked on: every scan, and TOF indices
/// `tof_lo, tof_lo + tof_step, ...` up to `tof_hi`.
pub struct Grid<'a> {
    pub num_scans: u32,
    pub im_at_scan: &'a dyn Fn(u32) -> f64,
    pub tof_to_mz: &'a dyn Fn(f64) -> f64,
    pub tof_lo: u32,
    pub tof_hi: u32,
    pub tof_step: u32,
}

/// Violation counts for one gate. See [`check_gate`].
#[derive(Debug, Default, Clone)]
pub struct Violations {
    /// Grid points checked.
    pub checked: usize,
    /// Point inside the original polygon but outside the gate.
    pub containment: usize,
    /// Point whose padded box meets the polygon, but neither it nor any grid
    /// neighbour (one scan / one TOF index away) is in the gate.
    pub full_padding: usize,
    /// Point inside the gate although its padded box grown by one scan and two
    /// TOF steps does not meet the polygon (gate too generous).
    pub over_inclusion: usize,
    /// First few violations, for the failure message.
    pub examples: Vec<String>,
}

impl Violations {
    pub fn any(&self) -> bool {
        self.containment + self.full_padding + self.over_inclusion > 0
    }
}

/// True when the gate holds `(scan, tof)` or a grid neighbour one scan and/or
/// one TOF index away.
fn near_gate(gate: &PolygonGate, s: u32, t: u32) -> bool {
    (s.saturating_sub(1)..=s + 1)
        .any(|ss| (t.saturating_sub(1)..=t + 1).any(|tt| gate.contains(ss, tt)))
}

/// Check `gate` (built from polygon `mz`/`im` with `mz_pad` Da / `im_pad` 1/K0)
/// against the brute-force reference on every point of `grid`:
///
/// * **containment** — inside polygon => inside gate (exact, no tolerance);
/// * **full padding** — the point's `±mz_pad x ±im_pad` box meets the polygon
///   => the gate holds the point or a neighbour one scan / one TOF index away
///   (discretization tolerance, e.g. an edge lying exactly on a scan line);
/// * **over-inclusion** — inside gate => the box grown by one scan and two TOF
///   steps meets the polygon. With both pads 0 this plus containment is "the
///   gate equals the polygon to within one scan / TOF index".
pub fn check_gate(
    gate: &PolygonGate,
    mz: &[f64],
    im: &[f64],
    grid: &Grid,
    mz_pad: f64,
    im_pad: f64,
) -> Violations {
    let mut v = Violations::default();
    let (bx0, bx1) = mz
        .iter()
        .fold((f64::MAX, f64::MIN), |a, &x| (a.0.min(x), a.1.max(x)));
    let (by0, by1) = im
        .iter()
        .fold((f64::MAX, f64::MIN), |a, &y| (a.0.min(y), a.1.max(y)));
    let note = |v: &mut Violations, kind: &str, s: u32, t: u32, x: f64, y: f64| {
        if v.examples.len() < 5 {
            v.examples.push(format!(
                "{kind}: scan {s} tof {t} (m/z {x:.4}, 1/K0 {y:.5})"
            ));
        }
    };
    for s in 0..grid.num_scans {
        let y = (grid.im_at_scan)(s);
        // One scan step (the larger of the two neighbours).
        let dy = [
            s.checked_sub(1),
            Some(s + 1).filter(|&n| n < grid.num_scans),
        ]
        .into_iter()
        .flatten()
        .map(|n| ((grid.im_at_scan)(n) - y).abs())
        .fold(0.0f64, f64::max);
        let mut t = grid.tof_lo;
        while t <= grid.tof_hi {
            let x = (grid.tof_to_mz)(t as f64);
            let dx = (grid.tof_to_mz)(t as f64 + 1.0) - x;
            let inside_gate = gate.contains(s, t);
            v.checked += 1;
            // Quick reject: grown box misses the polygon's bounding box.
            let (gx, gy) = (mz_pad + 2.0 * dx, im_pad + dy);
            let far = x + gx < bx0 || x - gx > bx1 || y + gy < by0 || y - gy > by1;
            if far {
                if inside_gate {
                    v.over_inclusion += 1;
                    note(&mut v, "over-inclusion", s, t, x, y);
                }
            } else if inside_gate {
                if !rect_meets_polygon(mz, im, x - gx, x + gx, y - gy, y + gy) {
                    v.over_inclusion += 1;
                    note(&mut v, "over-inclusion", s, t, x, y);
                }
            } else {
                if point_in_polygon(mz, im, x, y) {
                    v.containment += 1;
                    note(&mut v, "containment", s, t, x, y);
                }
                // Exact padded box; a miss is excused when the gate holds a
                // neighbouring grid point (one scan / one TOF index away).
                if rect_meets_polygon(mz, im, x - mz_pad, x + mz_pad, y - im_pad, y + im_pad)
                    && !near_gate(gate, s, t)
                {
                    v.full_padding += 1;
                    note(&mut v, "full-padding", s, t, x, y);
                }
            }
            t += grid.tof_step;
        }
    }
    v
}
