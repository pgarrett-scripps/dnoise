#![cfg(feature = "config")]
mod common;
use dnoise::{RunOptions, config::Config, denoise_with_options};

#[test]
fn optional_validation_preserves_valid_output_and_records_the_mode() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    let config = Config::default().resolve(&input).unwrap();
    let mut binaries = Vec::new();
    for skip_validation in [false, true] {
        let output = root.path().join(format!("skip-{skip_validation}.d"));
        denoise_with_options(
            &input,
            &output,
            &config.filter,
            &config.stages(),
            &RunOptions {
                skip_validation,
                ..Default::default()
            },
            |_| {},
        )
        .unwrap();
        let history = dnoise::provenance::read(&output).unwrap().unwrap();
        assert_eq!(
            history["history"][0]["validation"],
            if skip_validation {
                "structural"
            } else {
                "full"
            }
        );
        let recipe = Config::load(&output.join(dnoise::provenance::CONFIG_NAME)).unwrap();
        assert_eq!(recipe.skip_validation, Some(skip_validation));
        binaries.push(std::fs::read(output.join("analysis.tdf_bin")).unwrap());
    }
    assert_eq!(binaries[0], binaries[1]);
}

#[test]
fn skipping_full_validation_still_checks_structure_and_protects_files() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    let output = root.path().join("output.d");
    std::fs::create_dir(&output).unwrap();
    std::fs::write(output.join("keep"), b"original").unwrap();
    let config = Config::default().resolve(&input).unwrap();
    let db = rusqlite::Connection::open(input.join("analysis.tdf")).unwrap();
    // Full validation catches a payload in a frame incorrectly marked empty.
    // Structural mode trusts that metadata, demonstrating the actual tradeoff.
    db.execute("UPDATE Frames SET NumPeaks=0 WHERE Id=1", [])
        .unwrap();
    let run = |skip_validation| {
        denoise_with_options(
            &input,
            &output,
            &config.filter,
            &config.stages(),
            &RunOptions {
                force: true,
                skip_validation,
                ..Default::default()
            },
            |_| {},
        )
    };
    assert!(
        run(false)
            .unwrap_err()
            .to_string()
            .contains("counts disagree")
    );
    assert_eq!(std::fs::read(output.join("keep")).unwrap(), b"original");
    run(true).unwrap();
    assert!(
        dnoise::validation::read_frame(&output, 0)
            .unwrap()
            .1
            .is_empty()
    );
    let before = std::fs::read(output.join("analysis.tdf_bin")).unwrap();
    db.execute("UPDATE Frames SET TimsId=-1 WHERE Id=1", [])
        .unwrap();
    assert!(run(true).is_err());
    assert_eq!(
        std::fs::read(output.join("analysis.tdf_bin")).unwrap(),
        before
    );
    let input_before = std::fs::read(input.join("analysis.tdf_bin")).unwrap();
    assert!(
        denoise_with_options(
            &input,
            &input,
            &config.filter,
            &config.stages(),
            &RunOptions {
                force: true,
                skip_validation: true,
                ..Default::default()
            },
            |_| {}
        )
        .unwrap_err()
        .to_string()
        .contains("overlap")
    );
    assert_eq!(
        std::fs::read(input.join("analysis.tdf_bin")).unwrap(),
        input_before
    );
}

#[test]
fn dda_and_dia_outputs_validate_preserve_msms_and_match_streaming() {
    for dia in [false, true] {
        let root = tempfile::tempdir().unwrap();
        let input = common::fixture(root.path(), "input.d", dia);
        let output = root.path().join("output.d");
        let resolved = Config::default().resolve(&input).unwrap();
        let stages = resolved.stages();
        let stats = denoise_with_options(
            &input,
            &output,
            &resolved.filter,
            &stages,
            &RunOptions::default(),
            |_| {},
        )
        .unwrap();
        assert_eq!(
            dnoise::validation::inspect(&output, true).unwrap().points,
            stats.kept_points
        );
        let context = dnoise::RunContext::open(&input, &resolved.filter, &stages).unwrap();
        let reader = timsrust::readers::FrameReader::new(&output).unwrap();
        let raw = timsrust::readers::FrameReader::new(&input).unwrap();
        for i in 0..3 {
            let expected = context.process(i).unwrap().survivors;
            let actual = dnoise::FlatFrame::from_frame(&reader.get(i).unwrap());
            let mut actual: Vec<_> = (0..actual.len())
                .map(|j| (actual.scan[j], actual.tof[j], actual.intensity[j]))
                .collect();
            let mut expected = expected;
            actual.sort_unstable();
            expected.sort_unstable();
            assert_eq!(actual, expected);
            if i == 1 {
                assert_eq!(
                    reader.get(i).unwrap().intensities,
                    raw.get(i).unwrap().intensities
                );
                assert_eq!(
                    reader.get(i).unwrap().tof_indices,
                    raw.get(i).unwrap().tof_indices
                );
            }
        }
    }
}

#[test]
fn provenance_records_exact_recipe_history_and_crop_operation() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    let first = root.path().join("first.d");
    let second = root.path().join("second.d");
    let config = Config {
        halo_peak_fraction: Some(0.37),
        ms1_polygon_mz_pad: Some(7.5),
        ..Default::default()
    };
    let resolved = config.resolve(&input).unwrap();
    denoise_with_options(
        &input,
        &first,
        &resolved.filter,
        &resolved.stages(),
        &RunOptions::default(),
        |_| {},
    )
    .unwrap();
    let provenance = dnoise::provenance::read(&first).unwrap().unwrap();
    assert_eq!(
        provenance["history"][0]["effective_config"]["halo_peak_fraction"],
        0.37
    );
    assert_eq!(provenance["history"][0]["status"], "completed");
    let replay = Config::load(&first.join(dnoise::provenance::CONFIG_NAME)).unwrap();
    assert_eq!(replay.ms1_polygon_mz_pad, Some(7.5));
    let crop = dnoise::CropParams {
        rt_max: Some(0.4),
        ..Default::default()
    };
    denoise_with_options(
        &first,
        &second,
        &resolved.filter,
        &resolved.stages(),
        &RunOptions {
            crop: Some(&crop),
            crop_only: true,
            ..Default::default()
        },
        |_| {},
    )
    .unwrap();
    let history = dnoise::provenance::read(&second).unwrap().unwrap();
    assert_eq!(history["history"].as_array().unwrap().len(), 2);
    assert_eq!(history["history"][1]["operation"], "crop_only");
    let out = dnoise::validation::inspect(&second, true).unwrap();
    assert_eq!(out.frames, 4);
}

#[test]
fn dry_run_writes_no_marker_and_corrupt_record_is_rejected() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    let output = root.path().join("output.d");
    let config = Config::default().resolve(&input).unwrap();
    denoise_with_options(
        &input,
        &output,
        &config.filter,
        &config.stages(),
        &RunOptions {
            dry_run: true,
            ..Default::default()
        },
        |_| {},
    )
    .unwrap();
    assert!(!output.exists());
    assert!(dnoise::provenance::read(&input).unwrap().is_none());
    let db = rusqlite::Connection::open(input.join("analysis.tdf")).unwrap();
    db.execute("UPDATE Frames SET NumPeaks=999 WHERE Id=1", [])
        .unwrap();
    assert!(dnoise::validation::inspect(&input, true).is_err());
}

#[test]
fn unsorted_physical_offsets_are_valid_but_overlaps_are_not() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    let db = rusqlite::Connection::open(input.join("analysis.tdf")).unwrap();
    let first: i64 = db
        .query_row("SELECT TimsId FROM Frames WHERE Id=1", [], |r| r.get(0))
        .unwrap();
    let second: i64 = db
        .query_row("SELECT TimsId FROM Frames WHERE Id=2", [], |r| r.get(0))
        .unwrap();
    db.execute(
        "UPDATE Frames SET TimsId=CASE Id WHEN 1 THEN ?1 WHEN 2 THEN ?2 ELSE TimsId END",
        [second, first],
    )
    .unwrap();
    assert!(dnoise::validation::inspect(&input, true).is_ok());
    db.execute("UPDATE Frames SET TimsId=?1 WHERE Id=2", [second])
        .unwrap();
    assert!(dnoise::validation::inspect(&input, false).is_err());
}

#[test]
fn validation_agrees_across_pools_and_rejects_corruption_in_later_batches() {
    use std::io::{Seek, SeekFrom, Write};
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    let db = rusqlite::Connection::open(input.join("analysis.tdf")).unwrap();
    let record = dnoise::codec::encode_frame_type2(16, &[(3, 10000, 100)]);
    let path = input.join("analysis.tdf_bin");
    let mut binary = std::fs::OpenOptions::new()
        .append(true)
        .open(&path)
        .unwrap();
    let mut offset = binary.metadata().unwrap().len();
    for id in 5..=600 {
        db.execute(
            "INSERT INTO Frames SELECT ?1,Time,ScanMode,MsMsType,?2,MaxIntensity,100,NumScans,1,AccumulationTime,MzCalibration,TimsCalibration FROM Frames WHERE Id=1",
            rusqlite::params![id, offset],
        ).unwrap();
        binary.write_all(&record).unwrap();
        offset += record.len() as u64;
    }
    drop(binary);
    let before = std::fs::read(&path).unwrap();
    let pools: Vec<_> = [1, 4]
        .into_iter()
        .map(|threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
        })
        .collect();
    for pool in &pools {
        let valid = pool
            .install(|| dnoise::validation::inspect(&input, true))
            .unwrap();
        assert_eq!(
            (valid.frames, valid.points, valid.binary_bytes),
            (600, 623, offset)
        );
    }
    assert_eq!(std::fs::read(&path).unwrap(), before);
    db.execute("UPDATE Frames SET NumPeaks=2 WHERE Id=600", [])
        .unwrap();
    for pool in &pools {
        let error = pool
            .install(|| dnoise::validation::inspect(&input, true))
            .unwrap_err();
        assert!(error.to_string().contains("counts disagree for frame 600"));
    }
    db.execute("UPDATE Frames SET NumPeaks=1 WHERE Id=600", [])
        .unwrap();
    let mut binary = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
    binary
        .seek(SeekFrom::Start(offset - record.len() as u64 + 8))
        .unwrap();
    binary.write_all(&[0; 4]).unwrap();
    for pool in &pools {
        assert!(
            pool.install(|| dnoise::validation::inspect(&input, false))
                .is_ok()
        );
        assert!(matches!(
            pool.install(|| dnoise::validation::inspect(&input, true)),
            Err(dnoise::DnoiseError::Decode(_))
        ));
    }
}

#[test]
fn batching_does_not_change_output_and_invalid_crop_is_rejected() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    let resolved = Config::default().resolve(&input).unwrap();
    let mut outputs = Vec::new();
    for batch_size in [1, 64, 2048] {
        let output = root.path().join(format!("batch-{batch_size}.d"));
        denoise_with_options(
            &input,
            &output,
            &resolved.filter,
            &resolved.stages(),
            &RunOptions {
                frame_batch_size: Some(batch_size),
                ..Default::default()
            },
            |_| {},
        )
        .unwrap();
        outputs.push(std::fs::read(output.join("analysis.tdf_bin")).unwrap());
    }
    assert_eq!(outputs[0], outputs[1]);
    assert_eq!(outputs[1], outputs[2]);
    let config = Config {
        mz_min: Some(f64::NAN),
        ..Default::default()
    };
    assert!(config.resolve(&input).is_err());
}
