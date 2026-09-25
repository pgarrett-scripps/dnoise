//! Check the MS1 selection-polygon gate on real runs: read each run's IMS
//! PolygonFilter and calibration, describe the polygon's shape, build the gate
//! at several pads and run the brute-force containment / full-padding /
//! over-inclusion checks from `tests/common/polygon_ref.rs` on a
//! `(scan, TOF index)` grid, plus the gate's own build-time self-check.
//!
//! Usage: cargo run --release --example polygon_check -- <RUN.d>...
//! Prints one tab-separated row per run and pad.

#[path = "../tests/common/polygon_ref.rs"]
mod polygon_ref;

use dnoise::PolygonGate;
use dnoise::mobility::TimsCalibrationModel;
use dnoise::tsr::ConvertableDomain;
use dnoise::tsr::MetadataReader;
use polygon_ref::{Grid, check_gate};
use rusqlite::{Connection, OpenFlags, OptionalExtension};
use std::fmt::Write;
use std::path::Path;

const PADS: [(f64, f64); 5] = [
    (0.0, 0.0),
    (2.0, 0.01),
    (3.0, 0.02),
    (3.0, 0.0),
    (5.0, 0.05),
];

type Poly = (Vec<f64>, Vec<f64>);

/// Every stored polygon (one per `GroupProperties` group), as dnoise decodes it.
fn read_polygons(conn: &Connection) -> rusqlite::Result<Vec<Poly>> {
    let id = |name: &str| -> rusqlite::Result<Option<i64>> {
        conn.query_row(
            "SELECT Id FROM PropertyDefinitions WHERE PermanentName=?1",
            [name],
            |r| r.get(0),
        )
        .optional()
    };
    let (Some(mz_id), Some(im_id)) = (
        id("IMS_PolygonFilter_Mass")?,
        id("IMS_PolygonFilter_Mobility")?,
    ) else {
        return Ok(Vec::new());
    };
    let decode = |b: Vec<u8>| -> Vec<f64> {
        b.chunks_exact(8)
            .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
            .collect()
    };
    let mut st = conn.prepare(
        "SELECT a.Value, b.Value FROM GroupProperties a JOIN GroupProperties b \
         ON a.PropertyGroup = b.PropertyGroup WHERE a.Property=?1 AND b.Property=?2 \
         ORDER BY a.PropertyGroup",
    )?;
    let rows = st.query_map([mz_id, im_id], |r| {
        Ok((
            decode(r.get::<_, Vec<u8>>(0)?),
            decode(r.get::<_, Vec<u8>>(1)?),
        ))
    })?;
    rows.collect()
}

fn cross(o: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
}

/// Shape summary: (convex?, reflex vertex count, repeated vertices, proper
/// self-intersections, signed area in Th x 1/K0).
fn describe(mz: &[f64], im: &[f64]) -> (bool, usize, usize, usize, f64) {
    let n = mz.len();
    let p = |i: usize| (mz[i % n], im[i % n]);
    let area: f64 = (0..n)
        .map(|i| p(i).0 * p(i + 1).1 - p(i + 1).0 * p(i).1)
        .sum::<f64>()
        / 2.0;
    let sign = area.signum();
    let mut reflex = 0;
    for i in 0..n {
        let c = cross(p(i + n - 1), p(i), p(i + 1));
        if c * sign < -1e-12 {
            reflex += 1;
        }
    }
    let mut repeated = 0;
    for i in 0..n {
        if (0..i).any(|k| p(k) == p(i)) {
            repeated += 1;
        }
    }
    let mut crossings = 0;
    for i in 0..n {
        for k in i + 2..n {
            if i == 0 && k == n - 1 {
                continue;
            }
            let (a, b, c, d) = (p(i), p(i + 1), p(k), p(k + 1));
            let (d1, d2) = (cross(c, d, a), cross(c, d, b));
            let (d3, d4) = (cross(a, b, c), cross(a, b, d));
            if d1 * d2 < 0.0 && d3 * d4 < 0.0 {
                crossings += 1;
            }
        }
    }
    (
        reflex == 0 && crossings == 0 && repeated == 0,
        reflex,
        repeated,
        crossings,
        area,
    )
}

fn run(path: &Path, out: &mut String) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let name = path.file_name().unwrap().to_string_lossy().to_string();
    let tdf = path.join("analysis.tdf");
    let conn = Connection::open_with_flags(&tdf, OpenFlags::SQLITE_OPEN_READ_ONLY)?;
    let polys = read_polygons(&conn)?;
    let dia: i64 = conn
        .query_row(
            "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='DiaFrameMsMsWindows'",
            [],
            |r| r.get(0),
        )
        .unwrap_or(0);
    let dia = dia > 0
        && conn
            .query_row("SELECT COUNT(*) FROM DiaFrameMsMsWindows", [], |r| {
                r.get::<_, i64>(0)
            })
            .unwrap_or(0)
            > 0;
    let Some((mz, im)) = polys.first().cloned() else {
        writeln!(out, "{name}\tno polygon")?;
        return Ok(());
    };
    let distinct = {
        let mut d: Vec<&Poly> = Vec::new();
        for p in &polys {
            if !d.contains(&p) {
                d.push(p);
            }
        }
        d.len()
    };
    let (convex, reflex, repeated, crossings, area) = describe(&mz, &im);
    let verts: Vec<String> = mz
        .iter()
        .zip(&im)
        .map(|(x, y)| format!("({x:.1},{y:.3})"))
        .collect();
    writeln!(
        out,
        "# {name}: acq={} groups={} distinct={} vertices={} convex={convex} reflex={reflex} \
         repeated={repeated} self_crossings={crossings} area={area:.2} poly=[{}]",
        if dia { "dia" } else { "dda" },
        polys.len(),
        distinct,
        mz.len(),
        verts.join(" ")
    )?;

    let md = MetadataReader::new(&tdf).map_err(|e| format!("metadata: {e}"))?;
    let num_scans: i64 = conn.query_row("SELECT MAX(NumScans) FROM Frames", [], |r| r.get(0))?;
    let (model_type, c) = conn.query_row(
        "SELECT ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9 FROM TimsCalibration \
         WHERE Id = (SELECT TimsCalibration FROM Frames LIMIT 1)",
        [],
        |r| {
            let mut c = [0.0; 10];
            for (i, v) in c.iter_mut().enumerate() {
                *v = r.get(i + 1)?;
            }
            Ok((r.get::<_, i64>(0)?, c))
        },
    )?;
    let k0 = TimsCalibrationModel::new(model_type, c)?;
    let im_at = |s: u32| k0.scan_to_inv_k0(s as f64);
    let mz_to_tof = |m: f64| md.mz_converter.invert(m);
    let tof_to_mz = |t: f64| md.mz_converter.convert(t);

    let (x0, x1) = mz
        .iter()
        .fold((f64::MAX, f64::MIN), |a, &x| (a.0.min(x), a.1.max(x)));
    for (mz_pad, im_pad) in PADS {
        let Some(gate) = PolygonGate::build(
            &mz,
            &im,
            num_scans as usize,
            im_at,
            mz_to_tof,
            mz_pad,
            im_pad,
        ) else {
            writeln!(out, "{name}\t{mz_pad}\t{im_pad}\tno gate")?;
            continue;
        };
        let self_check = gate
            .check_contains_unpadded(&mz, &im, im_at, mz_to_tof)
            .map_or_else(|e| format!("FAIL {e}"), |()| "ok".to_string());
        let lo = mz_to_tof((x0 - mz_pad - 5.0).max(1.0)).max(0.0) as u32;
        let hi = mz_to_tof(x1 + mz_pad + 5.0) as u32;
        let grid = Grid {
            num_scans: num_scans as u32,
            im_at_scan: &im_at,
            tof_to_mz: &tof_to_mz,
            tof_lo: lo,
            tof_hi: hi,
            // ~6000 TOF samples per scan; every index near the edges is not
            // needed to catch a shape-level failure.
            tof_step: ((hi - lo) / 6000).max(1),
        };
        let v = check_gate(&gate, &mz, &im, &grid, mz_pad, im_pad);
        writeln!(
            out,
            "{name}\t{mz_pad}\t{im_pad}\tchecked={}\tcontainment={}\tfull_padding={}\tover_inclusion={}\tself_check={self_check}\t{:?}",
            v.checked, v.containment, v.full_padding, v.over_inclusion, v.examples
        )?;
    }
    Ok(())
}

fn main() {
    use rayon::prelude::*;
    let args: Vec<String> = std::env::args().skip(1).collect();
    let rows: Vec<String> = args
        .par_iter()
        .map(|arg| {
            let mut out = String::new();
            if let Err(e) = run(Path::new(arg), &mut out) {
                let _ = writeln!(out, "{arg}\terror: {e}");
            }
            out
        })
        .collect();
    for r in rows {
        print!("{r}");
    }
}
