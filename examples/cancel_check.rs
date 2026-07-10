//! Smoke-check the cooperative cancellation hook: run `denoise` with a token that
//! is already set, and confirm it returns `DnoiseError::Cancelled` without
//! producing a complete output.
//!
//! Usage: cargo run --example cancel_check -- <PATH.d>

use dnoise::{DnoiseError, FilterParams, RunOptions, Stages, denoise_with_options};
use std::path::PathBuf;
use std::sync::atomic::AtomicBool;

fn main() {
    let input: PathBuf = std::env::args()
        .nth(1)
        .expect("usage: cancel_check <PATH.d>")
        .into();
    let output = std::env::temp_dir().join("dnoise-cancel-check.d");
    let _ = std::fs::remove_dir_all(&output);

    let cancel = AtomicBool::new(true); // already cancelled
    let opts = RunOptions {
        force: true,
        cancel: Some(&cancel),
        ..RunOptions::default()
    };

    let result = denoise_with_options(
        &input,
        &output,
        &FilterParams::default(),
        &Stages::default(),
        &opts,
        |_| {},
    );

    match result {
        Err(DnoiseError::Cancelled) => println!("OK: returned Cancelled as expected"),
        Err(e) => {
            eprintln!("FAIL: expected Cancelled, got {e}");
            std::process::exit(1);
        }
        Ok(_) => {
            eprintln!("FAIL: expected Cancelled, run completed");
            std::process::exit(1);
        }
    }
    let _ = std::fs::remove_dir_all(&output);
}
