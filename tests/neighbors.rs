#![cfg(feature = "config")]
mod common;
use dnoise::{RunContext, RunOptions, config::Config, denoise_with_options};
use rusqlite::{Connection, params};
use std::io::Write;
use std::path::{Path, PathBuf};

type Points = Vec<(u32, u32, u32)>;
fn fixture(root: &Path, kind: &str) -> PathBuf {
    let path = if kind == "prm" {
        common::prm_fixture(root, "input.d")
    } else {
        common::fixture(root, "input.d", kind == "dia")
    };
    let db = Connection::open(path.join("analysis.tdf")).unwrap();
    db.pragma_update(None, "foreign_keys", false).unwrap();
    db.execute_batch("DELETE FROM Frames").unwrap();
    if kind == "prm" {
        db.execute_batch("DELETE FROM PrmFrameMsMsInfo; DELETE FROM PrmFrameMeasurementMode;")
            .unwrap();
    }
    if kind == "dia" {
        db.execute_batch("DELETE FROM DiaFrameMsMsInfo; DELETE FROM DiaFrameMsMsWindows;
            INSERT INTO DiaFrameMsMsWindows VALUES (1,3,7,500,2,30),(1,7,11,600,2,30),(2,3,7,700,2,30);").unwrap();
    }
    let mut file = std::fs::File::create(path.join("analysis.tdf_bin")).unwrap();
    file.write_all(&[0; 64]).unwrap();
    let mut offset = 64;
    for id in 1..=7 {
        let odd = id % 2 == 1;
        let points = if id == 7 {
            Vec::new()
        } else if odd {
            vec![(3 + (id - 1) / 2, 10000, 10)]
        } else {
            (3..7).map(|s| (s, 10000, 1000)).collect::<Points>()
        };
        let record = dnoise::codec::encode_frame_type2(16, &points);
        let typ = match kind {
            "prm" => 10,
            "dia" => 9,
            _ => {
                if odd {
                    0
                } else {
                    8
                }
            }
        };
        db.execute(
            "INSERT INTO Frames VALUES (?1,?1,8,?2,?3,?4,?5,16,?6,100,1,1)",
            params![
                id,
                typ,
                offset,
                points.iter().map(|p| p.2).max().unwrap_or(0),
                points.iter().map(|p| p.2 as u64).sum::<u64>(),
                points.len()
            ],
        )
        .unwrap();
        if kind == "prm" {
            db.execute(
                "INSERT INTO PrmFrameMsMsInfo VALUES (?1,3,7,500,2,30,?2)",
                params![id, if odd { 1 } else { 1000001 }],
            )
            .unwrap();
            db.execute(
                "INSERT INTO PrmFrameMsMsInfo VALUES (?1,7,11,600,2,30,1000001)",
                [id],
            )
            .unwrap();
        }
        if kind == "dia" {
            db.execute(
                "INSERT INTO DiaFrameMsMsInfo VALUES (?1,?2)",
                params![id, if odd { 1 } else { 2 }],
            )
            .unwrap();
        }
        file.write_all(&record).unwrap();
        offset += record.len();
    }
    path
}
fn config(radius: usize) -> Config {
    Config {
        ms1_neighbor_radius: Some(radius),
        prm_neighbor_radius: Some(radius),
        dia_neighbor_radius: Some(radius),
        denoise_msms: Some(true),
        halo: Some(false),
        mz_half_width: Some(0),
        msms_mz_half_width: Some(0),
        min_feature_length: Some(3),
        msms_min_feature_length: Some(3),
        max_internal_gap: Some(0),
        msms_max_internal_gap: Some(0),
        iterations: Some(1),
        msms_iterations: Some(1),
        dia_ms1_window: Some(false),
        ms1_polygon: Some(false),
        ..Default::default()
    }
}
fn run(input: &Path, config: &Config, index: usize) -> Points {
    let resolved = config.resolve(input).unwrap();
    let stages = resolved.stages();
    RunContext::open(input, &resolved.filter, &stages)
        .unwrap()
        .process(index)
        .unwrap()
        .survivors
}
fn sql(input: &Path, sql: &str) {
    Connection::open(input.join("analysis.tdf"))
        .unwrap()
        .execute_batch(sql)
        .unwrap();
}

#[test]
fn radius_uses_matching_observations_on_both_sides_and_native_points_only() {
    for kind in ["ms1", "prm", "dia"] {
        let root = tempfile::tempdir().unwrap();
        let input = fixture(root.path(), kind);
        let mut cfg = config(0);
        if kind == "ms1" {
            cfg.denoise_msms = Some(false);
        }
        assert!(run(&input, &cfg, 2).is_empty(), "{kind}");
        cfg.ms1_neighbor_radius = Some(1);
        cfg.prm_neighbor_radius = Some(1);
        cfg.dia_neighbor_radius = Some(1);
        assert_eq!(
            run(&input, &cfg, 2),
            vec![(4, 10000, 10)],
            "{kind}: previous and next support, no imported points or intensities"
        );
        assert!(
            run(&input, &cfg, 0).is_empty(),
            "{kind}: edge has insufficient support"
        );
        assert!(
            run(&input, &cfg, 4).is_empty(),
            "{kind}: empty observations consume a neighbor slot"
        );
        assert!(
            run(&input, &cfg, 6).is_empty(),
            "{kind}: empty frame stays empty"
        );
        cfg.ms1_neighbor_radius = Some(2);
        cfg.prm_neighbor_radius = Some(2);
        cfg.dia_neighbor_radius = Some(2);
        assert_eq!(
            run(&input, &cfg, 0),
            vec![(3, 10000, 10)],
            "{kind}: radius two"
        );
        cfg.neighbor_max_rt_gap = Some(1.0);
        assert!(run(&input, &cfg, 2).is_empty(), "{kind}: RT bound");
    }
}

#[test]
fn compatibility_and_calibration_boundaries_prevent_pooling() {
    for kind in ["ms1", "prm", "dia"] {
        let mut mutations = vec![
            "UPDATE Frames SET MzCalibration=2 WHERE Id=3",
            "UPDATE Frames SET TimsCalibration=2 WHERE Id=2", // even an intervening segment blocks support
            "UPDATE Frames SET ScanMode=99 WHERE Id=3",
            "UPDATE Frames SET AccumulationTime=200 WHERE Id=3",
        ];
        if kind == "prm" {
            mutations.extend([
                "INSERT INTO PrmTargets SELECT 99,NULL,Time,OneOverK0,MonoisotopicMz,Charge,Description FROM PrmTargets WHERE Id=1; UPDATE PrmFrameMsMsInfo SET Target=99 WHERE Frame=3 AND ScanNumBegin=3",
                "UPDATE PrmFrameMsMsInfo SET CollisionEnergy=31 WHERE Frame=3",
                "UPDATE PrmFrameMsMsInfo SET IsolationWidth=4 WHERE Frame=3",
                "UPDATE PrmFrameMsMsInfo SET IsolationMz=501 WHERE Frame=3",
                "UPDATE PrmFrameMsMsInfo SET ScanNumEnd=6 WHERE Frame=3 AND ScanNumBegin=3",
                "INSERT INTO PrmFrameMeasurementMode VALUES (3,'different')",
            ]);
        }
        if kind == "dia" {
            mutations.push("INSERT INTO DiaFrameMsMsWindows VALUES (3,3,7,500,2,31); UPDATE DiaFrameMsMsInfo SET WindowGroup=3 WHERE Frame=3");
        }
        for change in mutations {
            let root = tempfile::tempdir().unwrap();
            let input = fixture(root.path(), kind);
            let mut cfg = config(1);
            if kind == "ms1" {
                cfg.denoise_msms = Some(false);
            }
            sql(&input, change);
            assert!(run(&input, &cfg, 2).is_empty(), "{kind}: {change}");
        }
    }
}

#[test]
fn writer_streaming_dry_run_and_recipes_agree_for_all_three_workflows() {
    for kind in ["ms1", "prm", "dia"] {
        let root = tempfile::tempdir().unwrap();
        let input = fixture(root.path(), kind);
        let mut cfg = config(1);
        if kind == "ms1" {
            cfg.denoise_msms = Some(false);
        }
        let resolved = cfg.resolve(&input).unwrap();
        let stages = resolved.stages();
        let ctx = RunContext::open(&input, &resolved.filter, &stages).unwrap();
        let mut total = None;
        for (batch, threads) in [(1, 1), (4, 2)] {
            let output = root.path().join(format!("out-{batch}.d"));
            let stats = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    denoise_with_options(
                        &input,
                        &output,
                        &resolved.filter,
                        &stages,
                        &RunOptions {
                            frame_batch_size: Some(batch),
                            ..Default::default()
                        },
                        |_| {},
                    )
                    .unwrap()
                });
            for i in 0..ctx.len() {
                let native = dnoise::validation::read_frame(&input, i).unwrap().1;
                let written = dnoise::validation::read_frame(&output, i).unwrap().1;
                assert_eq!(written, ctx.process(i).unwrap().survivors);
                assert!(written.iter().all(|p| native.contains(p)));
            }
            let dry = denoise_with_options(
                &input,
                &root.path().join("unused"),
                &resolved.filter,
                &stages,
                &RunOptions {
                    dry_run: true,
                    frame_batch_size: Some(batch),
                    ..Default::default()
                },
                |_| {},
            )
            .unwrap();
            assert_eq!(dry.kept_points, stats.kept_points);
            if let Some(t) = total {
                assert_eq!(stats.kept_points, t);
            }
            total = Some(stats.kept_points);
            let recipe = Config::load(&output.join(dnoise::provenance::CONFIG_NAME)).unwrap();
            assert_eq!(recipe.ms1_neighbor_radius, Some(1));
            assert_eq!(recipe.prm_neighbor_radius, Some(1));
            assert_eq!(recipe.dia_neighbor_radius, Some(1));
            assert_eq!(run(&input, &recipe, 2), ctx.process(2).unwrap().survivors);
            assert_eq!(
                dnoise::validation::inspect(&output, true).unwrap().frames,
                7
            );
        }
    }
}

#[test]
fn cropping_excludes_neighbor_evidence_and_crop_only_bypasses_it() {
    for kind in ["ms1", "prm", "dia"] {
        let root = tempfile::tempdir().unwrap();
        let input = fixture(root.path(), kind);
        let mut cfg = config(1);
        if kind == "ms1" {
            cfg.denoise_msms = Some(false);
        }
        let resolved = cfg.resolve(&input).unwrap();
        let stages = resolved.stages();
        let crop = dnoise::CropParams {
            rt_min: Some(2.0 / 60.0),
            rt_max: Some(4.0 / 60.0),
            ..Default::default()
        };
        for crop_only in [false, true] {
            let output = root.path().join(format!("cropped-{crop_only}.d"));
            denoise_with_options(
                &input,
                &output,
                &resolved.filter,
                &stages,
                &RunOptions {
                    crop: Some(&crop),
                    crop_only,
                    ..Default::default()
                },
                |_| {},
            )
            .unwrap();
            let points = dnoise::validation::read_frame(&output, 2).unwrap().1;
            assert_eq!(
                points,
                if crop_only {
                    vec![(4, 10000, 10)]
                } else {
                    vec![]
                }
            );
        }
    }
}

#[test]
fn invalid_settings_and_malformed_dia_events_are_rejected() {
    let root = tempfile::tempdir().unwrap();
    let input = fixture(root.path(), "dia");
    for gap in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            Config {
                neighbor_max_rt_gap: Some(gap),
                ..config(1)
            }
            .resolve(&input)
            .is_err()
        );
        let stages = dnoise::Stages {
            neighbors: dnoise::NeighborParams {
                max_rt_gap_seconds: gap,
                ..Default::default()
            },
            ..Default::default()
        };
        assert!(RunContext::open(&input, &dnoise::FilterParams::default(), &stages).is_err());
    }
    let cfg = Config {
        denoise_msms: Some(false),
        ..config(1)
    }
    .resolve(&input)
    .unwrap();
    assert!(RunContext::open(&input, &cfg.filter, &cfg.stages()).is_err());
    sql(
        &input,
        "INSERT INTO DiaFrameMsMsWindows VALUES (1,6,9,700,2,30)",
    );
    let cfg = config(1).resolve(&input).unwrap();
    assert!(RunContext::open(&input, &cfg.filter, &cfg.stages()).is_err());
}

fn replace_points(input: &Path, changes: &[(usize, Points)]) {
    let mut frames: Vec<_> = (0..7)
        .map(|i| dnoise::validation::read_frame(input, i).unwrap())
        .collect();
    for (i, points) in changes {
        frames[*i].1 = points.clone();
    }
    let db = Connection::open(input.join("analysis.tdf")).unwrap();
    let mut file = std::fs::File::create(input.join("analysis.tdf_bin")).unwrap();
    file.write_all(&[0; 64]).unwrap();
    let mut offset = 64;
    for (i, (scans, points)) in frames.iter().enumerate() {
        let record = dnoise::codec::encode_frame_type2(*scans, points);
        file.write_all(&record).unwrap();
        db.execute("UPDATE Frames SET TimsId=?1,NumPeaks=?2,MaxIntensity=?3,SummedIntensities=?4 WHERE Id=?5",params![offset,points.len(),points.iter().map(|p|p.2).max().unwrap_or(0),points.iter().map(|p|p.2 as u64).sum::<u64>(),i+1]).unwrap();
        offset += record.len();
    }
}

#[test]
fn summed_intensity_support_never_rewrites_native_intensities() {
    for kind in ["ms1", "prm", "dia"] {
        let root = tempfile::tempdir().unwrap();
        let input = fixture(root.path(), kind);
        replace_points(
            &input,
            &[
                (0, vec![(3, 10000, 2)]),
                (2, vec![(3, 10000, 2)]),
                (4, vec![(3, 10000, 2)]),
            ],
        );
        let mut cfg = config(1);
        if kind == "ms1" {
            cfg.denoise_msms = Some(false);
        }
        cfg.min_feature_length = Some(1);
        cfg.msms_min_feature_length = Some(1);
        cfg.min_window_intensity = Some(5);
        cfg.msms_min_window_intensity = Some(5);
        assert_eq!(run(&input, &cfg, 2), vec![(3, 10000, 2)]);
        assert!(run(&input, &cfg, 0).is_empty()); // only two observations = four counts
    }
}

#[test]
fn halo_and_postprocessing_stay_inside_adjacent_events() {
    for kind in ["prm", "dia"] {
        let root = tempfile::tempdir().unwrap();
        let input = fixture(root.path(), kind);
        let mut points = vec![(4, 10000, 10)];
        points.extend((7..11).map(|s| (s, 10001, 1000)));
        replace_points(&input, &[(2, points)]);
        let mut cfg = config(1);
        cfg.halo = Some(true);
        cfg.halo_scan_half_width = Some(10);
        // With neighboring evidence the dim target qualifies, but would be
        // removed by the adjacent target's halo if either event were mixed.
        assert!(run(&input, &cfg, 2).contains(&(4, 10000, 10)));
        for centroid in ["none", "box", "watershed"] {
            cfg.smooth = Some(true);
            cfg.smooth_scan_half_width = Some(10);
            cfg.box_centroid = Some(centroid == "box");
            cfg.box_centroid_scan_half = Some(10);
            cfg.watershed = Some(centroid == "watershed");
            cfg.watershed_box_scan = Some(10);
            assert!(
                run(&input, &cfg, 2).contains(&(4, 10000, 10)),
                "{kind}: {centroid}"
            );
        }
    }
}

#[test]
fn dia_all_frames_can_use_neighbors_and_still_separates_events() {
    let root = tempfile::tempdir().unwrap();
    let input = fixture(root.path(), "dia");
    let cfg = Config {
        denoise_msms: Some(false),
        all_frames: Some(true),
        dia_per_window: Some(false),
        ..config(1)
    };
    assert_eq!(run(&input, &cfg, 2), vec![(4, 10000, 10)]);
    // Touching windows cannot jointly form a minimum-length feature.
    replace_points(
        &input,
        &[(
            2,
            vec![
                (5, 12000, 10),
                (6, 12000, 10),
                (7, 12000, 10),
                (8, 12000, 10),
            ],
        )],
    );
    assert!(run(&input, &cfg, 2).is_empty());
}

fn scanning_windows(input: &Path) {
    let db = Connection::open(input.join("analysis.tdf")).unwrap();
    db.execute("DELETE FROM DiaFrameMsMsWindows", []).unwrap();
    for group in [1, 2] {
        for scan in 0..16 {
            db.execute(
                "INSERT INTO DiaFrameMsMsWindows VALUES (?1,?2,?3,?4,25,?5)",
                params![
                    group,
                    scan,
                    scan + 1,
                    500.0 + 200.0 * (group - 1) as f64 - scan as f64,
                    30.0 - scan as f64 * 0.02
                ],
            )
            .unwrap();
        }
    }
}

#[test]
fn scanning_neighbors_match_entire_trajectory_and_preserve_native_points() {
    let root = tempfile::tempdir().unwrap();
    let input = fixture(root.path(), "dia");
    scanning_windows(&input);
    assert_eq!(run(&input, &config(1), 2), vec![(4, 10000, 10)]);
    assert!(run(&input, &config(0), 2).is_empty());
    // Same interval, endpoints, and most of the window shape; just one internal
    // m/z or collision-energy change must prevent support from that other shape.
    for field in ["IsolationMz", "CollisionEnergy"] {
        sql(
            &input,
            "DELETE FROM DiaFrameMsMsWindows WHERE WindowGroup=3; INSERT INTO DiaFrameMsMsWindows SELECT 3,ScanNumBegin,ScanNumEnd,IsolationMz,IsolationWidth,CollisionEnergy FROM DiaFrameMsMsWindows WHERE WindowGroup=1; UPDATE DiaFrameMsMsInfo SET WindowGroup=3 WHERE Frame=3;",
        );
        sql(
            &input,
            &format!(
                "UPDATE DiaFrameMsMsWindows SET {field}={field}+0.1 WHERE WindowGroup=3 AND ScanNumBegin=8"
            ),
        );
        assert!(run(&input, &config(1), 2).is_empty(), "{field}");
    }
}

#[test]
fn scanning_single_frame_filter_uses_continuous_steps_but_not_jumps() {
    let root = tempfile::tempdir().unwrap();
    let input = fixture(root.path(), "dia");
    scanning_windows(&input);
    replace_points(
        &input,
        &[(2, vec![(4, 10000, 10), (5, 10000, 10), (6, 10000, 10)])],
    );
    // Three native occupied scans qualify despite each vendor metadata row being
    // just one scan long. This is the synchro/midia encoding used by real files.
    assert_eq!(run(&input, &config(0), 2).len(), 3);
    sql(
        &input,
        "UPDATE DiaFrameMsMsWindows SET IsolationMz=IsolationMz-100 WHERE WindowGroup=1 AND ScanNumBegin>=6",
    );
    assert!(run(&input, &config(0), 2).is_empty());
}

#[test]
fn slice_static_boundaries_are_preserved_even_without_neighbors() {
    let root = tempfile::tempdir().unwrap();
    let input = fixture(root.path(), "dia");
    replace_points(
        &input,
        &[(
            2,
            vec![
                (5, 10000, 10),
                (6, 10000, 10),
                (7, 10000, 10),
                (8, 10000, 10),
            ],
        )],
    );
    assert!(run(&input, &config(0), 2).is_empty());
    // Explicit whole-frame override keeps its legacy meaning.
    let cfg = Config {
        dia_per_window: Some(false),
        ..config(0)
    };
    assert_eq!(run(&input, &cfg, 2).len(), 4);
}

#[test]
fn scanning_writer_streaming_and_dry_run_agree_and_preserve_geometry() {
    let root = tempfile::tempdir().unwrap();
    let input = fixture(root.path(), "dia");
    scanning_windows(&input);
    let cfg = config(1).resolve(&input).unwrap();
    let stages = cfg.stages();
    let ctx = RunContext::open(&input, &cfg.filter, &stages).unwrap();
    let output = root.path().join("output.d");
    let stats = denoise_with_options(
        &input,
        &output,
        &cfg.filter,
        &stages,
        &RunOptions {
            frame_batch_size: Some(1),
            ..Default::default()
        },
        |_| {},
    )
    .unwrap();
    assert!(stats.active_gates.dia_scan_varying);
    for i in 0..7 {
        assert_eq!(
            dnoise::validation::read_frame(&output, i).unwrap().1,
            ctx.process(i).unwrap().survivors
        );
    }
    let dry = denoise_with_options(
        &input,
        &root.path().join("unused"),
        &cfg.filter,
        &stages,
        &RunOptions {
            dry_run: true,
            ..Default::default()
        },
        |_| {},
    )
    .unwrap();
    assert_eq!(stats.kept_points, dry.kept_points);
    let rows = |p: &Path| {
        let db = Connection::open(p.join("analysis.tdf")).unwrap();
        let mut s=db.prepare("SELECT WindowGroup,ScanNumBegin,ScanNumEnd,IsolationMz,IsolationWidth,CollisionEnergy FROM DiaFrameMsMsWindows ORDER BY WindowGroup,ScanNumBegin").unwrap();
        s.query_map([], |r| {
            (0..6)
                .map(|i| r.get::<_, rusqlite::types::Value>(i))
                .collect::<Result<Vec<_>, _>>()
        })
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap()
    };
    assert_eq!(rows(&input), rows(&output));
}

#[test]
fn malformed_dia_geometry_is_rejected_without_neighbors() {
    for change in [
        "UPDATE DiaFrameMsMsWindows SET ScanNumBegin=-1 WHERE WindowGroup=1 AND ScanNumBegin=3",
        "UPDATE DiaFrameMsMsWindows SET ScanNumEnd=17 WHERE WindowGroup=1 AND ScanNumBegin=7",
        "UPDATE DiaFrameMsMsWindows SET CollisionEnergy=-1 WHERE WindowGroup=1",
        "UPDATE DiaFrameMsMsInfo SET WindowGroup=99 WHERE Frame=3",
        "DELETE FROM DiaFrameMsMsInfo WHERE Frame=3",
        "INSERT INTO DiaFrameMsMsInfo VALUES (999,1)",
    ] {
        let root = tempfile::tempdir().unwrap();
        let input = fixture(root.path(), "dia");
        sql(&input, change);
        let cfg = config(0).resolve(&input).unwrap();
        assert!(
            RunContext::open(&input, &cfg.filter, &cfg.stages()).is_err(),
            "{change}"
        );
    }
}

#[test]
fn neighbor_usage_counts_actual_evidence_excluding_central_empty_and_rt_gaps() {
    for kind in ["ms1", "prm", "dia"] {
        let root = tempfile::tempdir().unwrap();
        let input = fixture(root.path(), kind);
        let mut cfg = config(1);
        if kind == "ms1" {
            cfg.denoise_msms = Some(false);
        }
        for (gap, expected) in [(5.0, 2), (1.0, 0)] {
            cfg.neighbor_max_rt_gap = Some(gap);
            let resolved = cfg.resolve(&input).unwrap();
            let stages = resolved.stages();
            let ctx = RunContext::open(&input, &resolved.filter, &stages).unwrap();
            let usage = ctx.process(2).unwrap().neighbor_usage;
            assert_eq!(usage.events, 1, "{kind}");
            assert_eq!(usage.neighbors_used, expected, "{kind}");
            assert_eq!(usage.max_neighbors_used, expected);
            assert_eq!(usage.events_without_neighbors, u64::from(expected == 0));
            assert_eq!(ctx.process(6).unwrap().neighbor_usage.events, 0);
            let stats = denoise_with_options(
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
            let usage = if kind == "ms1" {
                stats.ms1_neighbor_usage
            } else {
                stats.msms_neighbor_usage
            };
            let total: u64 = (0..ctx.len())
                .map(|i| ctx.process(i).unwrap().neighbor_usage.neighbors_used)
                .sum();
            assert_eq!(usage.neighbors_used, total);
        }
    }
}
