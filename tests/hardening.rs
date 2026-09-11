#![cfg(feature = "config")]
mod common;
use dnoise::{CropParams, RunContext, RunOptions, config::Config, denoise_with_options};
use rusqlite::Connection;

#[test]
fn calibration_guards_cover_gates_crops_streaming_and_protect_output() {
    for column in ["MzCalibration", "TimsCalibration"] {
        let root = tempfile::tempdir().unwrap();
        let input = common::fixture(root.path(), "input.d", true);
        Connection::open(input.join("analysis.tdf"))
            .unwrap()
            .execute(&format!("UPDATE Frames SET {column}=2 WHERE Id=3"), [])
            .unwrap();
        let config = Config::default().resolve(&input).unwrap();
        let stages = config.stages();
        assert!(
            RunContext::open(&input, &config.filter, &stages)
                .err()
                .unwrap()
                .to_string()
                .contains("multiple calibration")
        );
        let output = root.path().join("output.d");
        std::fs::create_dir(&output).unwrap();
        std::fs::write(output.join("sentinel"), b"original").unwrap();
        for dry_run in [false, true] {
            let error = denoise_with_options(
                &input,
                &output,
                &config.filter,
                &stages,
                &RunOptions {
                    force: true,
                    dry_run,
                    ..Default::default()
                },
                |_| {},
            )
            .unwrap_err();
            assert!(error.to_string().contains("multiple calibration"));
            assert_eq!(std::fs::read(output.join("sentinel")).unwrap(), b"original");
        }
        let config = Config {
            dia_ms1_window: Some(false),
            ..Default::default()
        }
        .resolve(&input)
        .unwrap();
        let stages = config.stages();
        RunContext::open(&input, &config.filter, &stages).unwrap();
        for (crop, blocked) in [
            (
                CropParams {
                    mz_min: Some(100.0),
                    ..Default::default()
                },
                column == "MzCalibration",
            ),
            (
                CropParams {
                    im_min: Some(0.5),
                    ..Default::default()
                },
                column == "TimsCalibration",
            ),
            (
                CropParams {
                    rt_min: Some(0.25),
                    ..Default::default()
                },
                false,
            ),
            (
                CropParams {
                    min_intensity: Some(5),
                    ..Default::default()
                },
                false,
            ),
        ] {
            // Crop-only must not build even explicitly requested acquisition gates.
            let default_config = Config::default().resolve(&input).unwrap();
            let result = denoise_with_options(
                &input,
                &output,
                &config.filter,
                &default_config.stages(),
                &RunOptions {
                    dry_run: true,
                    crop_only: true,
                    crop: Some(&crop),
                    ..Default::default()
                },
                |_| {},
            );
            assert_eq!(result.is_err(), blocked, "{column}: {crop:?}");
            if let Err(error) = result {
                assert!(error.to_string().contains("multiple calibration"));
            }
        }
    }
}

#[test]
fn polygon_gate_rejects_multiple_calibrations_only_when_geometry_exists() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    let db = Connection::open(input.join("analysis.tdf")).unwrap();
    db.execute_batch("UPDATE Frames SET TimsCalibration=2 WHERE Id=3;
        CREATE TABLE PropertyDefinitions (Id INTEGER, PermanentName TEXT);
        CREATE TABLE GroupProperties (Property INTEGER, Value BLOB);
        INSERT INTO PropertyDefinitions VALUES (1,'IMS_PolygonFilter_Mass'),(2,'IMS_PolygonFilter_Mobility');").unwrap();
    let cfg = Config::default().resolve(&input).unwrap();
    RunContext::open(&input, &cfg.filter, &cfg.stages()).unwrap();
    for (id, values) in [
        (1, [100.0_f64, 1700.0, 1700.0, 100.0]),
        (2, [0.5_f64, 0.5, 1.5, 1.5]),
    ] {
        let bytes: Vec<_> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
        db.execute(
            "INSERT INTO GroupProperties VALUES (?1,?2)",
            rusqlite::params![id, bytes],
        )
        .unwrap();
    }
    assert!(
        RunContext::open(&input, &cfg.filter, &cfg.stages())
            .err()
            .unwrap()
            .to_string()
            .contains("MS1 polygon gate")
    );
}

#[test]
fn intensity_qc_includes_rt_removed_frames_and_matches_native_output() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", true);
    let cfg = Config {
        dia_ms1_window: Some(false),
        ..Default::default()
    }
    .resolve(&input)
    .unwrap();
    let crop = CropParams {
        rt_min: Some(0.4),
        ..Default::default()
    };
    let output = root.path().join("output.d");
    let options = RunOptions {
        crop: Some(&crop),
        crop_only: true,
        ..Default::default()
    };
    let stats = denoise_with_options(
        &input,
        &output,
        &cfg.filter,
        &cfg.stages(),
        &options,
        |_| {},
    )
    .unwrap();
    assert_eq!(stats.raw_ms1_summed_intensity, 1602);
    assert_eq!(stats.raw_msms_summed_intensity, 801);
    assert_eq!(stats.kept_ms1_summed_intensity, 801);
    assert_eq!(stats.kept_msms_summed_intensity, 0);
    assert_eq!(stats.raw_summed_intensity, 2403);
    let report = dnoise::provenance::report(&input, &cfg.filter, &cfg.stages(), &options, &stats);
    assert_eq!(report["stats"]["ms1_intensity_retained_pct"], 50.0);
    assert_eq!(report["stats"]["msms_intensity_retained_pct"], 0.0);
    assert_eq!(stats.ms1_neighbor_usage.events, 0);
    let empty = dnoise::DenoiseStats::default();
    let report = dnoise::provenance::report(&input, &cfg.filter, &cfg.stages(), &options, &empty);
    assert!(report["stats"]["ms1_intensity_retained_pct"].is_null());
    let dry = denoise_with_options(
        &input,
        &output,
        &cfg.filter,
        &cfg.stages(),
        &RunOptions {
            dry_run: true,
            ..options
        },
        |_| {},
    )
    .unwrap();
    assert_eq!(dry.raw_summed_intensity, stats.raw_summed_intensity);
    assert_eq!(dry.kept_summed_intensity, stats.kept_summed_intensity);
}
