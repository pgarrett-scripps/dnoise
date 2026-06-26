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
    /// * `mz_pad` widens each in-polygon m/z interval by this many Da per side and
    ///   `im_pad` widens the test by this much `1/K0` per side, so a precursor near
    ///   an edge keeps its isotopic envelope / mobility spread. `0.0`/`0.0`
    ///   reproduces the literal polygon.
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

        let mut per_scan: Vec<Vec<(u32, u32)>> = vec![Vec::new(); num_scans];
        for (s, slot) in per_scan.iter_mut().enumerate() {
            let y0 = im_at_scan(s as u32);

            // Union the inside m/z intervals across the padded 1/K0 band, so a
            // point within `im_pad` of the polygon (in mobility) is kept. With no
            // im_pad this is a single scan-line at y0.
            let mut spans: Vec<(f64, f64)> = Vec::new();
            for y in [y0 - im_pad, y0, y0 + im_pad] {
                spans.extend(scanline_spans(mz, im, y));
                if im_pad == 0.0 {
                    break;
                }
            }
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
        Some(Self { per_scan })
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
    fn keep_mask_matches_contains() {
        let (mz, im) = square();
        let gate = PolygonGate::build(&mz, &im, 100, id_im, id_tof, 0.0, 0.0).unwrap();
        let scan = [50, 50, 5];
        let tof = [50, 95, 50];
        assert_eq!(gate.keep_mask(&scan, &tof), vec![true, false, false]);
    }
}
