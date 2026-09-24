#![cfg(feature = "config")]
mod common;
use dnoise::{Acquisition, RunContext, RunOptions, config::Config, denoise_with_options};
use rusqlite::{Connection, types::Value};
use std::path::Path;

fn rows(path: &Path, query: &str) -> Vec<Vec<Value>> {
    let db = Connection::open(path.join("analysis.tdf")).unwrap();
    let mut stmt = db.prepare(query).unwrap();
    let cols = stmt.column_count();
    stmt.query_map([], |r| (0..cols).map(|i| r.get(i)).collect())
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap()
}

// Rebuild a tiny fixture after replacing selected frames' native points.
type Points = Vec<(u32, u32, u32)>;
fn replace_points(input: &Path, changes: &[(usize, Points)]) {
    use std::io::Write;
    let mut frames: Vec<_> = (0..4)
        .map(|i| dnoise::validation::read_frame(input, i).unwrap())
        .collect();
    for (index, points) in changes {
        frames[*index].1 = points.clone();
    }
    let db = Connection::open(input.join("analysis.tdf")).unwrap();
    let mut file = std::fs::File::create(input.join("analysis.tdf_bin")).unwrap();
    file.write_all(&[0; 64]).unwrap();
    let mut offset = 64;
    for (index, (scans, points)) in frames.iter().enumerate() {
        let record = dnoise::codec::encode_frame_type2(*scans, points);
        file.write_all(&record).unwrap();
        db.execute("UPDATE Frames SET TimsId=?1,NumPeaks=?2,MaxIntensity=?3,SummedIntensities=?4 WHERE Id=?5",
            rusqlite::params![offset as i64,points.len() as i64,points.iter().map(|p|p.2).max().unwrap_or(0),points.iter().map(|p|p.2 as u64).sum::<u64>() as i64,(index+1) as i64]).unwrap();
        offset += record.len();
    }
}

#[test]
fn prm_msms_respects_touching_events_and_matches_streaming_and_dry_run() {
    let root = tempfile::tempdir().unwrap();
    let input = common::prm_fixture(root.path(), "input.d");
    // Each target has four scans: together they would incorrectly qualify at five.
    for (minimum, expected) in [(3, 8), (5, 0)] {
        let config = Config {
            denoise_msms: Some(true),
            msms_min_feature_length: Some(minimum),
            halo: Some(false),
            ..Default::default()
        };
        let resolved = config.resolve(&input).unwrap();
        let stages = resolved.stages();
        let ctx = RunContext::open(&input, &resolved.filter, &stages).unwrap();
        let output = root.path().join(format!("out-{minimum}.d"));
        let stats = denoise_with_options(
            &input,
            &output,
            &resolved.filter,
            &stages,
            &RunOptions {
                frame_batch_size: Some(1),
                ..Default::default()
            },
            |_| {},
        )
        .unwrap();
        assert!(stats.active_gates.prm_per_event);
        assert_eq!(
            dnoise::validation::read_frame(&output, 1).unwrap().1.len(),
            expected
        );
        for index in 0..4 {
            let actual = dnoise::validation::read_frame(&output, index).unwrap().1;
            assert_eq!(actual, ctx.process(index).unwrap().survivors);
            let raw = dnoise::validation::read_frame(&input, index).unwrap().1;
            assert!(
                actual.iter().all(|p| raw.contains(p)),
                "native coordinates/intensities changed"
            );
        }
        let estimate = denoise_with_options(
            &input,
            &root.path().join("unused"),
            &resolved.filter,
            &stages,
            &RunOptions {
                dry_run: true,
                ..Default::default()
            },
            |_| {},
        )
        .unwrap();
        assert_eq!(
            (estimate.raw_points, estimate.kept_points),
            (stats.raw_points, stats.kept_points)
        );
        for query in [
            "SELECT * FROM PrmTargets ORDER BY Id",
            "SELECT * FROM PrmFrameMsMsInfo ORDER BY Frame,ScanNumBegin",
            "SELECT * FROM PrmFrameMeasurementMode ORDER BY Frame",
        ] {
            assert_eq!(rows(&input, query), rows(&output, query));
        }
        assert_eq!(
            dnoise::validation::inspect(&output, true).unwrap().frames,
            4
        );
    }
}

#[test]
fn prm_repeated_targets_do_not_pool_across_frames() {
    let root = tempfile::tempdir().unwrap();
    let input = common::prm_fixture(root.path(), "input.d");
    Connection::open(input.join("analysis.tdf")).unwrap().execute_batch(
        "UPDATE Frames SET MsMsType=10 WHERE Id=3; INSERT INTO PrmFrameMsMsInfo VALUES (3,3,7,500,2,30,1);"
    ).unwrap();
    replace_points(
        &input,
        &[
            (1, vec![(3, 10000, 100), (4, 10000, 100)]),
            (2, vec![(5, 10000, 100), (6, 10000, 100)]),
        ],
    );
    let resolved = Config {
        denoise_msms: Some(true),
        halo: Some(false),
        ..Default::default()
    }
    .resolve(&input)
    .unwrap();
    let stages = resolved.stages();
    let ctx = RunContext::open(&input, &resolved.filter, &stages).unwrap();
    assert!(ctx.process(1).unwrap().survivors.is_empty());
    assert!(ctx.process(2).unwrap().survivors.is_empty());
}

#[test]
fn prm_halo_and_postprocessing_do_not_mix_targets() {
    let root = tempfile::tempdir().unwrap();
    let input = common::prm_fixture(root.path(), "input.d");
    // Dim event next to a bright event within the halo box, at different TOFs.
    replace_points(
        &input,
        &[(
            1,
            (3..11)
                .map(|s| {
                    (
                        s,
                        if s < 7 { 10000 } else { 10001 },
                        if s < 7 { 10 } else { 1000 },
                    )
                })
                .collect(),
        )],
    );
    let resolved = Config {
        denoise_msms: Some(true),
        ..Default::default()
    }
    .resolve(&input)
    .unwrap();
    let stages = resolved.stages();
    assert_eq!(
        RunContext::open(&input, &resolved.filter, &stages)
            .unwrap()
            .process(1)
            .unwrap()
            .survivors
            .len(),
        8
    );
    // Both events share a TOF; smoothing/centroiding must still remain separate.
    replace_points(
        &input,
        &[(
            1,
            (3..11)
                .map(|s| (s, 10000, if s < 7 { 10 } else { 1000 }))
                .collect(),
        )],
    );
    for config in [
        Config {
            smooth: Some(true),
            smooth_scan_half_width: Some(20),
            ..Default::default()
        },
        Config {
            box_centroid: Some(true),
            box_centroid_scan_half: Some(20),
            ..Default::default()
        },
        Config {
            watershed: Some(true),
            watershed_box_scan: Some(20),
            ..Default::default()
        },
    ] {
        let resolved = Config {
            denoise_msms: Some(true),
            halo: Some(false),
            ..config
        }
        .resolve(&input)
        .unwrap();
        let stages = resolved.stages();
        let ctx = RunContext::open(&input, &resolved.filter, &stages).unwrap();
        let points = ctx.process(1).unwrap().survivors;
        let sums: Vec<u64> = [(3, 7), (7, 11)]
            .iter()
            .map(|&(b, e)| {
                points
                    .iter()
                    .filter(|p| p.0 >= b && p.0 < e)
                    .map(|p| p.2 as u64)
                    .sum()
            })
            .collect();
        assert_eq!(sums, vec![40, 4000]);
        if stages.smooth.is_some() {
            assert_eq!(points.len(), 8);
        }
    }
}

#[test]
fn prm_default_preserves_fragments_metadata_and_streaming_parity() {
    let root = tempfile::tempdir().unwrap();
    let input = common::prm_fixture(root.path(), "input.d");
    assert_eq!(
        dnoise::detect_acquisition(&input).unwrap(),
        Acquisition::PrmPasef
    );
    // A real polygon that excludes the synthetic MS1 signal must be ignored.
    let db = Connection::open(input.join("analysis.tdf")).unwrap();
    db.execute_batch("CREATE TABLE PropertyDefinitions (Id INTEGER, PermanentName TEXT);
        INSERT INTO PropertyDefinitions VALUES (1,'IMS_PolygonFilter_Mass'),(2,'IMS_PolygonFilter_Mobility');
        CREATE TABLE GroupProperties (Property INTEGER, Value BLOB);").unwrap();
    for (id, vertices) in [
        (1, [1500.0_f64, 1600.0, 1600.0, 1500.0]),
        (2, [0.6, 0.6, 0.7, 0.7]),
    ] {
        let blob: Vec<u8> = vertices.iter().flat_map(|v| v.to_le_bytes()).collect();
        db.execute(
            "INSERT INTO GroupProperties VALUES (?1,?2)",
            rusqlite::params![id, blob],
        )
        .unwrap();
    }
    let resolved = Config::default().resolve(&input).unwrap();
    let stages = resolved.stages();
    let ctx = RunContext::open(&input, &resolved.filter, &stages).unwrap();
    let output = root.path().join("output.d");
    let stats = denoise_with_options(
        &input,
        &output,
        &resolved.filter,
        &stages,
        &RunOptions::default(),
        |_| {},
    )
    .unwrap();
    assert!(stats.kept_ms1_points > 0 && stats.kept_ms1_points < stats.raw_ms1_points);
    assert!(
        !stats.active_gates.ms1_polygon
            && !stats.active_gates.dia_ms1
            && !stats.active_gates.dia_window
            && !stats.active_gates.dda_window
    );
    for i in 0..4 {
        let actual = dnoise::validation::read_frame(&output, i).unwrap();
        assert_eq!(actual.1, ctx.process(i).unwrap().survivors);
        if i % 2 == 1 {
            assert_eq!(actual, dnoise::validation::read_frame(&input, i).unwrap());
        }
    }
    for query in [
        "SELECT * FROM PrmTargets ORDER BY Id",
        "SELECT * FROM PrmFrameMsMsInfo ORDER BY Frame,ScanNumBegin",
        "SELECT * FROM PrmFrameMeasurementMode ORDER BY Frame",
        "SELECT * FROM GroupProperties ORDER BY Property",
        "SELECT Id,Time,ScanMode,MsMsType,NumScans,AccumulationTime,MzCalibration,TimsCalibration FROM Frames ORDER BY Id",
    ] {
        assert_eq!(rows(&input, query), rows(&output, query));
    }
    let report = dnoise::provenance::read(&output).unwrap().unwrap();
    assert_eq!(report["history"][0]["acquisition"], "PrmPasef");
    assert_eq!(
        dnoise::validation::inspect(&output, true).unwrap().frames,
        4
    );
}

#[test]
fn prm_all_frames_is_rejected_without_replacing_output() {
    let root = tempfile::tempdir().unwrap();
    let input = common::prm_fixture(root.path(), "input.d");
    let output = root.path().join("output.d");
    std::fs::create_dir(&output).unwrap();
    std::fs::write(output.join("keep"), b"existing output").unwrap();
    for config in [
        Config {
            all_frames: Some(true),
            ..Default::default()
        },
        Config {
            all_frames: Some(true),
            denoise_msms: Some(true),
            ..Default::default()
        },
    ] {
        let resolved = config.resolve(&input).unwrap();
        let stages = resolved.stages();
        for dry_run in [false, true] {
            let err = denoise_with_options(
                &input,
                &output,
                &resolved.filter,
                &stages,
                &RunOptions {
                    force: true,
                    dry_run,
                    skip_validation: true,
                    ..Default::default()
                },
                |_| {},
            )
            .unwrap_err();
            assert!(
                err.to_string()
                    .contains("--all-frames (filter_all_frames) is unsupported for prm-PASEF")
            );
            assert_eq!(
                std::fs::read(output.join("keep")).unwrap(),
                b"existing output"
            );
        }
        let err = RunContext::open(&input, &resolved.filter, &stages)
            .err()
            .unwrap();
        assert!(
            err.to_string()
                .contains("--all-frames (filter_all_frames) is unsupported for prm-PASEF")
        );
    }
}

#[test]
fn malformed_prm_metadata_is_rejected_even_without_decoding() {
    for sql in [
        "DROP TABLE PrmTargets",
        "DROP TABLE PrmFrameMsMsInfo",
        "DELETE FROM PrmTargets WHERE Id=1",
        "UPDATE PrmFrameMsMsInfo SET Frame=99 WHERE Frame=2",
        "UPDATE PrmFrameMsMsInfo SET Frame=1 WHERE Frame=2",
        "UPDATE PrmFrameMsMsInfo SET ScanNumBegin=-1 WHERE Frame=2 AND ScanNumBegin=3",
        "UPDATE PrmFrameMsMsInfo SET ScanNumEnd=17 WHERE Frame=2",
        "UPDATE PrmFrameMsMsInfo SET ScanNumEnd=ScanNumBegin WHERE Frame=2",
        "UPDATE PrmFrameMsMsInfo SET ScanNumEnd=8 WHERE Frame=2 AND ScanNumBegin=3",
        "UPDATE PrmFrameMsMsInfo SET IsolationWidth=-1 WHERE Frame=2",
        "UPDATE PrmFrameMsMsInfo SET IsolationMz=1e999 WHERE Frame=2",
        "UPDATE PrmFrameMsMsInfo SET CollisionEnergy=-1 WHERE Frame=2",
        "DELETE FROM PrmFrameMsMsInfo WHERE Frame=2",
        "INSERT INTO DiaFrameMsMsInfo VALUES (2,1)",
        "INSERT INTO PrmFrameMeasurementMode VALUES (99,NULL)",
    ] {
        let root = tempfile::tempdir().unwrap();
        let input = common::prm_fixture(root.path(), "input.d");
        let db = Connection::open(input.join("analysis.tdf")).unwrap();
        db.pragma_update(None, "foreign_keys", false).unwrap();
        db.execute_batch(sql).unwrap();
        assert!(dnoise::detect_acquisition(&input).is_err(), "{sql}");
        assert!(dnoise::validation::inspect(&input, false).is_err(), "{sql}");
        let params = dnoise::FilterParams::default();
        let stages = dnoise::Stages::default();
        assert!(RunContext::open(&input, &params, &stages).is_err(), "{sql}");
    }
}

#[test]
fn empty_prm_tables_do_not_change_dda_or_dia_detection() {
    for dia in [false, true] {
        let root = tempfile::tempdir().unwrap();
        let input = common::fixture(root.path(), "input.d", dia);
        Connection::open(input.join("analysis.tdf")).unwrap().execute_batch(
            "CREATE TABLE PrmFrameMsMsInfo (Frame INTEGER); CREATE TABLE PrmTargets (Id INTEGER);"
        ).unwrap();
        assert_eq!(
            dnoise::detect_acquisition(&input).unwrap(),
            if dia {
                Acquisition::DiaPasef
            } else {
                Acquisition::DdaPasef
            }
        );
    }
}

#[test]
fn mixed_runs_disable_geometry_and_reject_fragment_denoising() {
    let root = tempfile::tempdir().unwrap();
    let input = common::prm_fixture(root.path(), "input.d");
    Connection::open(input.join("analysis.tdf")).unwrap().execute_batch("UPDATE Frames SET MsMsType=8 WHERE Id=3;
        CREATE TABLE PasefFrameMsMsInfo (Frame INTEGER, ScanNumBegin INTEGER, ScanNumEnd INTEGER, Precursor INTEGER);
        INSERT INTO PasefFrameMsMsInfo VALUES (3,3,11,1);").unwrap();
    assert_eq!(
        dnoise::detect_acquisition(&input).unwrap(),
        Acquisition::Mixed
    );
    let resolved = Config::default().resolve(&input).unwrap();
    let output = root.path().join("output.d");
    denoise_with_options(
        &input,
        &output,
        &resolved.filter,
        &resolved.stages(),
        &RunOptions::default(),
        |_| {},
    )
    .unwrap();
    for i in 1..4 {
        assert_eq!(
            dnoise::validation::read_frame(&input, i).unwrap(),
            dnoise::validation::read_frame(&output, i).unwrap()
        );
    }
    let resolved = Config {
        denoise_msms: Some(true),
        ..Default::default()
    }
    .resolve(&input)
    .unwrap();
    let stages = resolved.stages();
    assert!(
        RunContext::open(&input, &resolved.filter, &stages)
            .err()
            .unwrap()
            .to_string()
            .contains("mixed acquisition")
    );
}

#[test]
fn prm_only_run_reports_no_ms1_opportunity_and_cropping_is_explicit() {
    let root = tempfile::tempdir().unwrap();
    let input = common::prm_fixture(root.path(), "input.d");
    Connection::open(input.join("analysis.tdf"))
        .unwrap()
        .execute_batch(
            "UPDATE Frames SET MsMsType=10;
        INSERT INTO PrmFrameMsMsInfo VALUES (1,3,11,500,2,30,1),(3,3,11,500,2,30,1);",
        )
        .unwrap();
    let resolved = Config::default().resolve(&input).unwrap();
    let stages = resolved.stages();
    let options = RunOptions {
        dry_run: true,
        ..Default::default()
    };
    let stats = denoise_with_options(
        &input,
        &root.path().join("unused"),
        &resolved.filter,
        &stages,
        &options,
        |_| {},
    )
    .unwrap();
    assert_eq!(stats.raw_points, stats.kept_points);
    assert_eq!(stats.ms1_frames, 0);
    assert!(
        dnoise::provenance::report(&input, &resolved.filter, &stages, &options, &stats)["warnings"]
            .to_string()
            .contains("no MS1 frames")
    );
    let crop = dnoise::CropParams {
        rt_min: Some(25.0),
        ..Default::default()
    };
    let output = root.path().join("cropped.d");
    denoise_with_options(
        &input,
        &output,
        &resolved.filter,
        &stages,
        &RunOptions {
            crop: Some(&crop),
            crop_only: true,
            ..Default::default()
        },
        |_| {},
    )
    .unwrap();
    assert!(
        dnoise::validation::read_frame(&output, 1)
            .unwrap()
            .1
            .is_empty()
    );
    assert_eq!(
        rows(
            &input,
            "SELECT * FROM PrmFrameMsMsInfo ORDER BY Frame,ScanNumBegin"
        ),
        rows(
            &output,
            "SELECT * FROM PrmFrameMsMsInfo ORDER BY Frame,ScanNumBegin"
        )
    );
}
