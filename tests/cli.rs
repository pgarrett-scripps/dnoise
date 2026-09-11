#![cfg(feature = "cli")]
mod common;
use std::process::Command;
fn cli() -> Command {
    Command::new(env!("CARGO_BIN_EXE_dnoise"))
}

#[test]
fn prm_cli_detects_preserves_and_supports_opt_in_fragment_filtering() {
    let root = tempfile::tempdir().unwrap();
    let input = common::prm_fixture(root.path(), "input.d");
    let validation = cli().arg("validate").arg(&input).output().unwrap();
    assert!(validation.status.success());
    assert!(String::from_utf8_lossy(&validation.stdout).contains("prm-PASEF"));
    let output = root.path().join("output.d");
    let result = cli().arg(&input).arg(&output).output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        dnoise::validation::read_frame(&input, 1).unwrap(),
        dnoise::validation::read_frame(&output, 1).unwrap()
    );
    let filtered = root.path().join("filtered.d");
    let result = cli()
        .arg(&input)
        .arg(&filtered)
        .arg("--denoise-msms")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        dnoise::validation::read_frame(&filtered, 1)
            .unwrap()
            .1
            .len(),
        8
    );
    let report = dnoise::provenance::read(&filtered).unwrap().unwrap();
    assert_eq!(
        report["history"][0]["stats"]["active_gates"]["prm_per_event"],
        true
    );
    {
        let flag = "--all-frames";
        let rejected = root.path().join("rejected.d");
        let result = cli().arg(&input).arg(&rejected).arg(flag).output().unwrap();
        assert!(!result.status.success());
        assert!(
            String::from_utf8_lossy(&result.stderr)
                .contains("--all-frames (filter_all_frames) is unsupported for prm-PASEF")
        );
        assert!(!rejected.exists());
    }
}

#[test]
fn validation_flags_override_configs_in_single_and_batch_runs() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    for batch in [false, true] {
        for (index, (configured_skip, flag, expected_skip)) in [
            (false, None, false),
            (true, None, true),
            (false, Some("--skip-validation"), true),
            (true, Some("--validate"), false),
        ]
        .into_iter()
        .enumerate()
        {
            let output = root.path().join(format!("output-{batch}-{index}.d"));
            let mut command = cli();
            if batch {
                let manifest = root.path().join("batch.json");
                std::fs::write(&manifest, serde_json::to_vec(&serde_json::json!({
                    "schema_version":1, "jobs":[{"input":input,"output":output,"config":{"skip_validation":configured_skip}}]
                })).unwrap()).unwrap();
                command.arg("batch").arg(manifest);
            } else {
                let config = root.path().join("settings.toml");
                std::fs::write(&config, format!("skip_validation = {configured_skip}\n")).unwrap();
                command.arg(&input).arg(&output).arg("--config").arg(config);
            }
            if let Some(flag) = flag {
                command.arg(flag);
            }
            let result = command.output().unwrap();
            assert!(
                result.status.success(),
                "{}",
                String::from_utf8_lossy(&result.stderr)
            );
            let history = dnoise::provenance::read(&output).unwrap().unwrap();
            assert_eq!(
                history["history"][0]["validation"],
                if expected_skip { "structural" } else { "full" }
            );
            assert_eq!(
                history["history"][0]["effective_config"]["skip_validation"],
                expected_skip
            );
        }
    }
    assert!(
        !cli()
            .arg(&input)
            .args(["--dry-run", "--skip-validation", "--validate"])
            .output()
            .unwrap()
            .status
            .success()
    );
}

#[test]
fn cli_flags_override_config_and_metadata_is_queryable() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    let output = root.path().join("output.d");
    let config = root.path().join("settings.toml");
    std::fs::write(&config, "halo = true\nmz_half_width = 8\n").unwrap();
    let result = cli()
        .args([input.as_os_str(), output.as_os_str()])
        .arg("--config")
        .arg(config)
        .args(["--no-halo", "--mz-half-width", "3", "--threads", "2"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let recipe =
        dnoise::config::Config::load(&output.join(dnoise::provenance::CONFIG_NAME)).unwrap();
    assert_eq!(recipe.halo, Some(false));
    assert_eq!(recipe.threads, Some(2));
    assert_eq!(recipe.mz_half_width, Some(3));
    assert!(
        cli()
            .arg("validate")
            .arg(&output)
            .output()
            .unwrap()
            .status
            .success()
    );
    let metadata = cli().arg("metadata").arg(&output).output().unwrap();
    let value: serde_json::Value = serde_json::from_slice(&metadata.stdout).unwrap();
    assert_eq!(value["history"][0]["status"], "completed");
}

#[test]
fn batch_paths_are_relative_to_manifest_and_collisions_write_nothing() {
    let root = tempfile::tempdir().unwrap();
    common::fixture(root.path(), "input.d", true);
    let manifest = root.path().join("batch.json");
    let job = serde_json::json!({"input":"input.d","output":"output.d"});
    std::fs::write(
        &manifest,
        serde_json::to_vec(
            &serde_json::json!({"schema_version":1,"jobs":[job.clone(),job.clone()]}),
        )
        .unwrap(),
    )
    .unwrap();
    assert!(
        !cli()
            .arg("batch")
            .arg(&manifest)
            .output()
            .unwrap()
            .status
            .success()
    );
    assert!(!root.path().join("output.d").exists());
    std::fs::write(
        &manifest,
        serde_json::to_vec(&serde_json::json!({"schema_version":1,"jobs":[job]})).unwrap(),
    )
    .unwrap();
    let result = cli()
        .current_dir(std::env::temp_dir())
        .arg("batch")
        .arg(&manifest)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let report: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(report["failures"], 0);
    assert!(
        root.path()
            .join("output.d/dnoise.provenance.json")
            .is_file()
    );
}

#[test]
fn report_cannot_overwrite_input_even_in_dry_run() {
    let root = tempfile::tempdir().unwrap();
    let input = common::fixture(root.path(), "input.d", false);
    let tdf = input.join("analysis.tdf");
    let original = std::fs::read(&tdf).unwrap();
    let result = cli()
        .arg(&input)
        .args(["--dry-run", "--report"])
        .arg(&tdf)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(&tdf).unwrap(), original);
}

#[test]
fn neighbor_cli_overrides_aliases_and_exports_replayable_settings() {
    let root = tempfile::tempdir().unwrap();
    let input = common::prm_fixture(root.path(), "input.d");
    let cfg = root.path().join("settings.toml");
    std::fs::write(&cfg,"frame_half_width = 2\nprm_neighbor_radius = 2\ndia_neighbor_radius = 2\nneighbor_max_rt_gap = 1.0\n").unwrap();
    for flag in ["--ms1-neighbor-radius", "--frame-half-width"] {
        let output = root.path().join(format!("{flag}.d"));
        let result = cli()
            .arg(&input)
            .arg(&output)
            .arg("--config")
            .arg(&cfg)
            .args([
                "--denoise-msms",
                flag,
                "1",
                "--prm-neighbor-radius",
                "1",
                "--dia-neighbor-radius",
                "0",
                "--neighbor-max-rt-gap",
                "30",
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let recipe =
            dnoise::config::Config::load(&output.join(dnoise::provenance::CONFIG_NAME)).unwrap();
        assert_eq!(recipe.frame_half_width, None);
        assert_eq!(recipe.ms1_neighbor_radius, Some(1));
        assert_eq!(recipe.prm_neighbor_radius, Some(1));
        assert_eq!(recipe.dia_neighbor_radius, Some(0));
        assert_eq!(recipe.neighbor_max_rt_gap, Some(30.0));
        let report = dnoise::provenance::read(&output).unwrap().unwrap();
        assert_eq!(
            report["history"][0]["stats"]["active_gates"]["prm_neighbors"],
            true
        );
        assert_eq!(
            report["history"][0]["stats"]["active_gates"]["dia_neighbors"],
            false
        );
    }
    let invalid = cli()
        .arg(&input)
        .args(["--dry-run", "--prm-neighbor-radius", "1"])
        .output()
        .unwrap();
    assert!(!invalid.status.success());
    assert!(String::from_utf8_lossy(&invalid.stderr).contains("requires MS/MS denoising"));
}
