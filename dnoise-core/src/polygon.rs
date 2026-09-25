//! ddaPASEF/PASEF MS1 selection-polygon gate.
//!
//! timsTOF PASEF acquisition methods restrict precursor selection to a polygon in
//! the `(m/z, 1/K0)` plane — the "IMS PolygonFilter", stored in `analysis.tdf`.
//! MS1 signal outside this polygon sits in a region the instrument never schedules
//! for fragmentation: in ddaPASEF it can never become a precursor, so it is noise
//! as far as identification is concerned and can be dropped from the survey scans.
//! This gate reproduces that selection region as a hard mask on MS1 points.
//!
//! Like [`crate::dia_ms1`] it is calibration-free at run time: the writer converts
//! the polygon vertices `(m/z, 1/K0)` into a per-scan list of TOF-index intervals
//! once (using the run calibration plus any padding), and each MS1 point is then
//! tested with a single binary search.

/// Per-scan merged TOF intervals for the interior of the selection polygon. A
/// point `(scan, tof)` is kept iff `tof` lands in one of `per_scan[scan]`'s
/// intervals.
#[derive(Debug)]
pub struct PolygonGate {
    /// `per_scan[s]` = sorted, non-overlapping `[tof_lo, tof_hi]` intervals
    /// (inclusive) of the polygon interior at mobility scan `s`. Empty rows keep
    /// nothing (the polygon does not cover that mobility).
    per_scan: Vec<Vec<(u32, u32)>>,
    /// Gate whole features by any overlap ([`crate::overlap`]) instead of points
    /// (`false` = point by point).
    pub overlap: bool,
}

impl PolygonGate {
    /// Build the gate from polygon vertices in `(m/z, 1/K0)` space.
    ///
    /// * `mz` / `im` — parallel vertex coordinates (the ring is closed
    ///   automatically; the last vertex connects back to the first).
    /// * `num_scans` — sizes the per-scan table.
    /// * `im_at_scan(s)` — the `1/K0` of mobility scan `s` (run calibration).
    /// * `mz_to_tof(mz)` — the fractional TOF index of an m/z; must be monotonic
    ///   increasing (it is, for the timsTOF √-law calibration).
    /// * `mz_pad` widens each in-polygon m/z interval by this many Th per side and
    ///   `im_pad` widens the test by this much `1/K0` per side, so a precursor near
    ///   an edge keeps its isotopic envelope / mobility spread. `0.0`/`0.0`
    ///   reproduces the literal polygon. With both pads the gate is the polygon's
    ///   Minkowski sum with the `±mz_pad x ±im_pad` rectangle, sampled at each
    ///   scan's 1/K0 (exact for any ring shape, including thin spikes and sharp
    ///   vertices between scans); a non-positive `im_pad` means no mobility pad.
    ///
    /// Returns `None` when the polygon is degenerate (< 3 vertices, mismatched
    /// lengths), `num_scans == 0`, or the polygon covers no scan's mobility (so the
    /// caller skips the gate rather than dropping every MS1 point).
    pub fn build(
        mz: &[f64],
        im: &[f64],
        num_scans: usize,
        im_at_scan: impl Fn(u32) -> f64,
        mz_to_tof: impl Fn(f64) -> f64,
        mz_pad: f64,
        im_pad: f64,
    ) -> Option<Self> {
        let n = mz.len();
        if n < 3 || im.len() != n || num_scans == 0 {
            return None;
        }

        // 1/K0 levels where the band projection can change shape (vertices and
        // edge crossings); only needed when there is a mobility pad.
        let critical = if im_pad > 0.0 {
            critical_levels(mz, im)
        } else {
            Vec::new()
        };

        let mut per_scan: Vec<Vec<(u32, u32)>> = vec![Vec::new(); num_scans];
        for (s, slot) in per_scan.iter_mut().enumerate() {
            let y0 = im_at_scan(s as u32);

            // With no im_pad: the single scan-line at y0. With an im_pad: the
            // exact m/z projection of the polygon inside the whole band
            // [y0 - im_pad, y0 + im_pad], so a point within `im_pad` of the
            // polygon (in mobility) is kept however thin or pointed the polygon
            // is there. (Sampling a few lines in the band misses spikes and
            // sharp vertices lying between them.)
            let mut spans = if im_pad > 0.0 {
                band_spans(mz, im, &critical, y0 - im_pad, y0 + im_pad)
            } else {
                scanline_spans(mz, im, y0)
            };
            if spans.is_empty() {
                continue;
            }

            // Merge in m/z space after padding each span; converting the merged,
            // ordered spans through the monotone mz->tof keeps them ordered.
            spans.sort_by(|a, b| a.0.total_cmp(&b.0));
            let mut tof_iv: Vec<(u32, u32)> = Vec::new();
            let mut cur = (spans[0].0 - mz_pad, spans[0].1 + mz_pad);
            for &(lo, hi) in &spans[1..] {
                let (lo, hi) = (lo - mz_pad, hi + mz_pad);
                if lo <= cur.1 {
                    cur.1 = cur.1.max(hi);
                } else {
                    push_tof_interval(&mut tof_iv, cur, &mz_to_tof);
                    cur = (lo, hi);
                }
            }
            push_tof_interval(&mut tof_iv, cur, &mz_to_tof);

            // Coalesce touching/adjacent TOF intervals (rounding can make merged
            // m/z spans abut) so membership is a clean binary search.
            tof_iv.sort_unstable();
            let mut merged: Vec<(u32, u32)> = Vec::with_capacity(tof_iv.len());
            for (lo, hi) in tof_iv {
                match merged.last_mut() {
                    Some(last) if lo <= last.1.saturating_add(1) => last.1 = last.1.max(hi),
                    _ => merged.push((lo, hi)),
                }
            }
            *slot = merged;
        }

        if per_scan.iter().all(|row| row.is_empty()) {
            return None;
        }
        Some(Self {
            per_scan,
            overlap: false,
        })
    }

    /// Build-time self-check: every scan's *unpadded* polygon interval (the
    /// literal scan-line at that scan's 1/K0, converted to TOF exactly as
    /// [`Self::build`] does) must lie inside one of this gate's intervals for the
    /// same scan. Padding may only ever add points, so a failure means the gate
    /// would drop signal inside the instrument's own selection polygon — a bug,
    /// or a negative pad. Pass the same polygon and converters used to build.
    /// Costs one extra unpadded build (a few thousand scan-lines).
    pub fn check_contains_unpadded(
        &self,
        mz: &[f64],
        im: &[f64],
        im_at_scan: impl Fn(u32) -> f64,
        mz_to_tof: impl Fn(f64) -> f64,
    ) -> Result<(), String> {
        let Some(base) = Self::build(mz, im, self.per_scan.len(), im_at_scan, mz_to_tof, 0.0, 0.0)
        else {
            return Ok(()); // the polygon covers no scan: nothing to contain
        };
        for (s, (inner, outer)) in base.per_scan.iter().zip(&self.per_scan).enumerate() {
            for &(lo, hi) in inner {
                let i = outer.partition_point(|&(_, o_hi)| o_hi < lo);
                if !(i < outer.len() && outer[i].0 <= lo && hi <= outer[i].1) {
                    return Err(format!(
                        "scan {s}: polygon TOF interval [{lo}, {hi}] is not inside the padded \
                         gate's intervals {outer:?}"
                    ));
                }
            }
        }
        Ok(())
    }

    /// Inclusive TOF hull `(lo, hi)` of every interval the gate keeps, over all
    /// scans; `None` when the gate keeps nothing. No kept point lies outside it.
    pub fn tof_span(&self) -> Option<(u32, u32)> {
        let lo = self
            .per_scan
            .iter()
            .filter_map(|r| r.first())
            .map(|i| i.0)
            .min()?;
        let hi = self
            .per_scan
            .iter()
            .filter_map(|r| r.last())
            .map(|i| i.1)
            .max()?;
        Some((lo, hi))
    }

    /// True when the point lies inside the polygon (and should be kept).
    pub fn contains(&self, scan: u32, tof: u32) -> bool {
        let Some(row) = self.per_scan.get(scan as usize) else {
            return false;
        };
        let i = row.partition_point(|&(_, hi)| hi < tof);
        i < row.len() && tof >= row[i].0
    }

    /// Per-point keep mask (in input order): `false` for points outside the
    /// polygon. `scan` and `tof` are parallel.
    pub fn keep_mask(&self, scan: &[u32], tof: &[u32]) -> Vec<bool> {
        scan.iter()
            .zip(tof)
            .map(|(&s, &t)| self.contains(s, t))
            .collect()
    }
}

/// Convert an `(m/z_lo, m/z_hi)` span to an inclusive TOF interval via the
/// monotone `mz_to_tof`, widening to whole indices (floor lo / ceil hi) so a point
/// exactly on the edge is kept. Clamps negative m/z (after padding) to 0.
fn push_tof_interval(
    out: &mut Vec<(u32, u32)>,
    (lo, hi): (f64, f64),
    mz_to_tof: &impl Fn(f64) -> f64,
) {
    let t_lo = mz_to_tof(lo.max(0.0)).floor().max(0.0) as u32;
    let t_hi = mz_to_tof(hi.max(0.0)).ceil().max(0.0) as u32;
    if t_hi >= t_lo {
        out.push((t_lo, t_hi));
    }
}

/// m/z spans where the horizontal line at `y` (in 1/K0) lies inside the polygon,
/// found by the even-odd ray-crossing rule over all edges. Returns sorted,
/// pairwise `(lo, hi)` m/z spans.
fn scanline_spans(mz: &[f64], im: &[f64], y: f64) -> Vec<(f64, f64)> {
    let n = mz.len();
    let mut xs: Vec<f64> = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        let (yi, yj) = (im[i], im[j]);
        // Half-open crossing test: counts each edge once and is robust at vertices.
        if (yi > y) != (yj > y) {
            let t = (y - yi) / (yj - yi);
            xs.push(mz[i] + t * (mz[j] - mz[i]));
        }
    }
    xs.sort_by(f64::total_cmp);
    xs.chunks_exact(2).map(|c| (c[0], c[1])).collect()
}

/// Sorted, deduplicated 1/K0 levels at which the set of polygon edges crossing
/// a horizontal line, or their left-to-right order, can change: every vertex,
/// plus every proper crossing of two non-adjacent edges (only a self-intersecting
/// ring has those). Between two consecutive levels each inside span's ends move
/// linearly with 1/K0.
fn critical_levels(mz: &[f64], im: &[f64]) -> Vec<f64> {
    let n = mz.len();
    let mut ys: Vec<f64> = im.to_vec();
    for i in 0..n {
        let (ax, ay, bx, by) = (mz[i], im[i], mz[(i + 1) % n], im[(i + 1) % n]);
        for k in i + 2..n {
            if i == 0 && k == n - 1 {
                continue; // adjacent through the ring closure
            }
            let (cx, cy, dx, dy) = (mz[k], im[k], mz[(k + 1) % n], im[(k + 1) % n]);
            let den = (bx - ax) * (dy - cy) - (by - ay) * (dx - cx);
            if den == 0.0 {
                continue; // parallel or collinear: order along y cannot swap
            }
            let t = ((cx - ax) * (dy - cy) - (cy - ay) * (dx - cx)) / den;
            let u = ((cx - ax) * (by - ay) - (cy - ay) * (bx - ax)) / den;
            if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
                ys.push(ay + t * (by - ay));
            }
        }
    }
    ys.retain(|y| y.is_finite());
    ys.sort_by(f64::total_cmp);
    ys.dedup();
    ys
}

/// m/z spans covering the polygon's intersection with the closed band
/// `a <= 1/K0 <= b` (`a < b`), projected onto m/z: exactly the m/z values `x`
/// for which some point `(x, y)` with `y` in the band lies in the (closed)
/// polygon. `critical` is [`critical_levels`] of the same polygon.
///
/// The band is cut at the critical levels into slabs. Inside a slab the same
/// edges cross every horizontal line in the same order, so each inside span's
/// ends move linearly, and the union of that span over the slab is simply the
/// hull of its extents at the slab's two boundary lines. Spans may overlap; the
/// caller merges them.
fn band_spans(mz: &[f64], im: &[f64], critical: &[f64], a: f64, b: f64) -> Vec<(f64, f64)> {
    let n = mz.len();
    let first = critical.partition_point(|&y| y <= a);
    let mut levels = Vec::with_capacity(critical.len().min(8) + 2);
    levels.push(a);
    levels.extend(critical[first..].iter().copied().take_while(|&y| y < b));
    levels.push(b);

    let mut spans = Vec::new();
    // Per crossing edge: (m/z at slab middle, m/z at slab bottom, m/z at top).
    let mut xs: Vec<(f64, f64, f64)> = Vec::new();
    for w in levels.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        if hi <= lo {
            continue;
        }
        let mid = 0.5 * (lo + hi);
        xs.clear();
        for i in 0..n {
            let j = (i + 1) % n;
            let (yi, yj) = (im[i], im[j]);
            // Same half-open rule as `scanline_spans`. No vertex lies strictly
            // inside the slab, so an edge crossing `mid` spans all of [lo, hi].
            if (yi > mid) != (yj > mid) {
                let at = |y: f64| {
                    let t = ((y - yi) / (yj - yi)).clamp(0.0, 1.0);
                    mz[i] + t * (mz[j] - mz[i])
                };
                xs.push((at(mid), at(lo), at(hi)));
            }
        }
        xs.sort_by(|p, q| p.0.total_cmp(&q.0));
        for c in xs.chunks_exact(2) {
            spans.push((c[0].1.min(c[0].2), c[1].1.max(c[1].2)));
        }
    }
    spans
}

#[cfg(test)]
mod tests {
    use super::*;

    // Identity converters: TOF index == m/z, scan == 1/K0. Lets the tests reason
    // in polygon coordinates directly.
    fn id_im(s: u32) -> f64 {
        s as f64
    }
    fn id_tof(mz: f64) -> f64 {
        mz
    }

    /// A 100x100 axis-aligned square with corners (10,10)-(90,90) in (m/z, 1/K0).
    fn square() -> (Vec<f64>, Vec<f64>) {
        (vec![10.0, 90.0, 90.0, 10.0], vec![10.0, 10.0, 90.0, 90.0])
    }

    #[test]
    fn degenerate_polygon_yields_no_gate() {
        assert!(
            PolygonGate::build(&[0.0, 1.0], &[0.0, 1.0], 100, id_im, id_tof, 0.0, 0.0).is_none()
        );
    }

    #[test]
    fn square_keeps_inside_drops_outside() {
        let (mz, im) = square();
        let gate = PolygonGate::build(&mz, &im, 100, id_im, id_tof, 0.0, 0.0).unwrap();
        assert!(gate.contains(50, 50)); // center
        assert!(gate.contains(10, 10)); // lower-left corner (boundary kept)
        assert!(gate.contains(89, 89)); // interior just inside the upper-right
        assert!(!gate.contains(50, 9)); // left of the m/z span
        assert!(!gate.contains(50, 91)); // right of the m/z span
        assert!(!gate.contains(5, 50)); // below the polygon's mobility range
        assert!(!gate.contains(95, 50)); // above it
    }

    #[test]
    fn mz_pad_widens_the_kept_span() {
        let (mz, im) = square();
        let gate = PolygonGate::build(&mz, &im, 100, id_im, id_tof, 5.0, 0.0).unwrap();
        assert!(gate.contains(50, 7)); // 3 past the nominal left edge, within pad
        assert!(gate.contains(50, 94)); // within the right pad
        assert!(!gate.contains(50, 4)); // beyond the pad
    }

    #[test]
    fn im_pad_widens_the_mobility_band() {
        let (mz, im) = square();
        let gate = PolygonGate::build(&mz, &im, 100, id_im, id_tof, 0.0, 3.0).unwrap();
        assert!(gate.contains(8, 50)); // 2 below the nominal mobility edge, within pad
        assert!(gate.contains(92, 50)); // within the upper pad
        assert!(!gate.contains(5, 50)); // beyond the pad
    }

    #[test]
    fn concave_polygon_gives_two_spans_on_a_scan() {
        // A "U": full base, with a notch cut from the top center (m/z 30..70) down
        // to 1/K0 30, so a scan line through the arms cuts two disjoint m/z spans.
        let mz = vec![0.0, 100.0, 100.0, 70.0, 70.0, 30.0, 30.0, 0.0];
        let im = vec![0.0, 0.0, 100.0, 100.0, 30.0, 30.0, 100.0, 0.0];
        let gate = PolygonGate::build(&mz, &im, 101, id_im, id_tof, 0.0, 0.0).unwrap();
        // Through the arms (scan 50): left and right spans in, notch out.
        assert!(gate.contains(50, 15)); // left arm
        assert!(gate.contains(50, 85)); // right arm
        assert!(!gate.contains(50, 50)); // in the notch
        // Below the notch floor (scan 10) the interior is one solid span.
        assert!(gate.contains(10, 50));
    }

    #[test]
    fn im_pad_catches_a_sliver_between_scan_lines() {
        // Body plus a sideways spike to m/z 90 that lies wholly between scan
        // lines 50 and 51. Sampling y0 - pad, y0, y0 + pad (the pre-0.5 rule)
        // misses it from scan 45; the band projection keeps it.
        let mz = vec![0.0, 20.0, 20.0, 90.0, 20.0, 20.0, 0.0];
        let im = vec![0.0, 0.0, 50.2, 50.5, 50.8, 100.0, 100.0];
        let gate = PolygonGate::build(&mz, &im, 101, id_im, id_tof, 0.0, 6.0).unwrap();
        for s in 45..=56 {
            assert!(gate.contains(s, 80), "scan {s}");
        }
        assert!(!gate.contains(44, 80));
        assert!(!gate.contains(57, 80));
    }

    #[test]
    fn band_keeps_a_sharp_vertex_between_sampled_lines() {
        // Triangle apex at m/z 50, 1/K0 50.5. With pad 3 the apex is within
        // pad of scans 48..=53 and must be kept there at m/z 50.
        let mz = vec![0.0, 100.0, 50.0];
        let im = vec![0.0, 0.0, 50.5];
        let gate = PolygonGate::build(&mz, &im, 101, id_im, id_tof, 0.0, 3.0).unwrap();
        for s in 48..=53 {
            assert!(gate.contains(s, 50), "scan {s}");
        }
        assert!(!gate.contains(54, 50));
    }

    #[test]
    fn self_check_accepts_padded_and_rejects_shrunk_gates() {
        let mz = vec![0.0, 100.0, 50.0];
        let im = vec![0.0, 0.0, 50.5];
        for (mzp, imp) in [(0.0, 0.0), (2.0, 0.0), (0.0, 3.0), (5.0, 5.0)] {
            let g = PolygonGate::build(&mz, &im, 101, id_im, id_tof, mzp, imp).unwrap();
            assert!(g.check_contains_unpadded(&mz, &im, id_im, id_tof).is_ok());
        }
        let shrunk = PolygonGate::build(&mz, &im, 101, id_im, id_tof, -3.0, 0.0).unwrap();
        assert!(
            shrunk
                .check_contains_unpadded(&mz, &im, id_im, id_tof)
                .is_err()
        );
    }

    #[test]
    fn keep_mask_matches_contains() {
        let (mz, im) = square();
        let gate = PolygonGate::build(&mz, &im, 100, id_im, id_tof, 0.0, 0.0).unwrap();
        let scan = [50, 50, 5];
        let tof = [50, 95, 50];
        assert_eq!(gate.keep_mask(&scan, &tof), vec![true, false, false]);
    }
}
