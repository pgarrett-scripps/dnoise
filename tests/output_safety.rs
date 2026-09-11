mod common;
use dnoise::{FilterParams, RunOptions, Stages, denoise, denoise_in_place, denoise_with_options};
use std::sync::atomic::{AtomicBool, Ordering};

#[test]
fn identical_and_nested_outputs_leave_input_unchanged() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    let before = std::fs::read(input.join("analysis.tdf_bin")).unwrap();
    for output in [
        input.clone(),
        input.join("nested.d"),
        root.path().to_owned(),
    ] {
        assert!(
            denoise(
                &input,
                &output,
                &FilterParams::default(),
                &Stages::default(),
                true
            )
            .is_err()
        );
        assert_eq!(
            std::fs::read(input.join("analysis.tdf_bin")).unwrap(),
            before
        );
    }
}

#[test]
fn invalid_input_preserves_existing_output() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    std::fs::write(input.join("analysis.tdf"), b"invalid").unwrap();
    let out = root.path().join("output.d");
    std::fs::create_dir(&out).unwrap();
    std::fs::write(out.join("keep"), b"original").unwrap();
    assert!(
        denoise(
            &input,
            &out,
            &FilterParams::default(),
            &Stages::default(),
            true
        )
        .is_err()
    );
    assert_eq!(std::fs::read(out.join("keep")).unwrap(), b"original");
}

#[test]
fn cancellation_preserves_existing_output() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    let out = root.path().join("output.d");
    std::fs::create_dir(&out).unwrap();
    std::fs::write(out.join("keep"), b"original").unwrap();
    let cancel = AtomicBool::new(false);
    let options = RunOptions {
        force: true,
        cancel: Some(&cancel),
        ..Default::default()
    };
    let result = denoise_with_options(
        &input,
        &out,
        &FilterParams::default(),
        &Stages::default(),
        &options,
        |p| {
            if p.frames_done > 0 {
                cancel.store(true, Ordering::Relaxed);
            }
        },
    );
    assert!(matches!(result, Err(dnoise::DnoiseError::Cancelled)));
    assert_eq!(std::fs::read(out.join("keep")).unwrap(), b"original");
}

#[test]
fn in_place_matches_separate_output_without_fixed_backup_names() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    std::fs::write(input.join("dnoise.config.toml"), b"iterations = 999\n").unwrap();
    let out = root.path().join("output.d");
    let old = root.path().join("input.d.dnoise-old");
    std::fs::create_dir(&old).unwrap();
    std::fs::write(old.join("keep"), b"unrelated").unwrap();
    denoise(
        &input,
        &out,
        &FilterParams::default(),
        &Stages::default(),
        false,
    )
    .unwrap();
    denoise_in_place(
        &input,
        &FilterParams::default(),
        &Stages::default(),
        &RunOptions::default(),
        |_| {},
    )
    .unwrap();
    assert_eq!(
        std::fs::read(input.join("analysis.tdf_bin")).unwrap(),
        std::fs::read(out.join("analysis.tdf_bin")).unwrap()
    );
    assert_eq!(std::fs::read(old.join("keep")).unwrap(), b"unrelated");
    for path in [&input, &out] {
        #[cfg(feature = "config")]
        assert!(
            !std::fs::read_to_string(path.join("dnoise.config.toml"))
                .unwrap()
                .contains("999")
        );
        #[cfg(not(feature = "config"))]
        assert!(!path.join("dnoise.config.toml").exists());
    }
}

#[cfg(unix)]
#[test]
fn copy_failure_discards_staging_without_touching_existing_output() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    let external = root.path().join("external");
    std::fs::write(&external, b"external original").unwrap();
    std::os::unix::fs::symlink(&external, input.join("linked-file")).unwrap();
    let output = root.path().join("output.d");
    std::fs::create_dir(&output).unwrap();
    std::fs::write(output.join("keep"), b"old output").unwrap();
    assert!(
        denoise(
            &input,
            &output,
            &FilterParams::default(),
            &Stages::default(),
            true
        )
        .is_err()
    );
    assert_eq!(std::fs::read(output.join("keep")).unwrap(), b"old output");
    assert_eq!(std::fs::read(external).unwrap(), b"external original");
    assert!(!std::fs::read_dir(root.path()).unwrap().flatten().any(|e| {
        e.file_name()
            .to_string_lossy()
            .starts_with(".dnoise-stage-")
    }));
}
