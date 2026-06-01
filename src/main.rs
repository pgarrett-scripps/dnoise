//! dnoise CLI.

use anyhow::{Context, Result};
use clap::Parser;
use dnoise::FilterParams;
use serde::Deserialize;
use std::path::{Path, PathBuf};

/// Denoise a Bruker timsTOF .d folder via the iterative vertical-IM feature filter.
#[derive(Parser)]
#[command(name = "dnoise", version, about)]
struct Cli {
    /// Input Bruker .d folder.
    input: PathBuf,
    /// Output .d folder (created; must not exist unless --force).
    output: PathBuf,

    /// TOML config file with filter parameters. Explicit CLI flags override its values.
    #[arg(long, short = 'c', value_name = "FILE")]
    config: Option<PathBuf>,

    /// Column half-width in TOF indices.
    #[arg(long)]
    mz_half_width: Option<u32>,
    /// Minimum total span (scans) of a kept feature.
    #[arg(long)]
    min_feature_length: Option<usize>,
    /// Max empty scans tolerated inside a feature.
    #[arg(long)]
    max_internal_gap: Option<usize>,
    /// Per-scan summed-intensity floor for occupancy.
    #[arg(long)]
    min_window_intensity: Option<u64>,
    /// Total summed-intensity floor for a kept feature.
    #[arg(long)]
    min_feature_intensity: Option<u64>,
    /// Filter passes (each re-applies to the previous survivors).
    #[arg(long)]
    iterations: Option<usize>,

    /// Pre-filter smoothing: replace each MS1 frame with the centered running
    /// average of its `2r+1` MS1-frame neighborhood before filtering (0 = off).
    #[arg(long)]
    frame_half_width: Option<usize>,

    /// Filter MS/MS frames too. By default only MS1 frames are filtered (the
    /// vertical-IM filter is MS1-specific and strips most MS/MS fragment signal).
    #[arg(long)]
    all_frames: bool,
    /// Worker threads (default: all cores).
    #[arg(long)]
    threads: Option<usize>,
    /// Overwrite the output folder if it already exists.
    #[arg(long)]
    force: bool,
}

/// On-disk config. Every field is optional; missing keys fall back to CLI flags,
/// then to [`FilterParams::default`]. Unknown keys are rejected to catch typos.
#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct FileConfig {
    mz_half_width: Option<u32>,
    min_feature_length: Option<usize>,
    max_internal_gap: Option<usize>,
    min_window_intensity: Option<u64>,
    min_feature_intensity: Option<u64>,
    iterations: Option<usize>,
    frame_half_width: Option<usize>,
    all_frames: Option<bool>,
    threads: Option<usize>,
}

impl FileConfig {
    fn load(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("reading config file {}", path.display()))?;
        toml::from_str(&text).with_context(|| format!("parsing config file {}", path.display()))
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let cfg = match &cli.config {
        Some(path) => FileConfig::load(path)?,
        None => FileConfig::default(),
    };

    // Precedence for each knob: explicit CLI flag > config file > built-in default.
    let defaults = FilterParams::default();
    let params = FilterParams {
        mz_half_width: cli
            .mz_half_width
            .or(cfg.mz_half_width)
            .unwrap_or(defaults.mz_half_width),
        min_feature_length: cli
            .min_feature_length
            .or(cfg.min_feature_length)
            .unwrap_or(defaults.min_feature_length),
        max_internal_gap: cli
            .max_internal_gap
            .or(cfg.max_internal_gap)
            .unwrap_or(defaults.max_internal_gap),
        min_window_intensity: cli
            .min_window_intensity
            .or(cfg.min_window_intensity)
            .unwrap_or(defaults.min_window_intensity),
        min_feature_intensity: cli
            .min_feature_intensity
            .or(cfg.min_feature_intensity)
            .unwrap_or(defaults.min_feature_intensity),
        num_iterations: cli
            .iterations
            .or(cfg.iterations)
            .unwrap_or(defaults.num_iterations),
    };

    // Pre-filter smoothing radius (decoupled from the filter knobs above): explicit
    // CLI flag > config file > 0 (off).
    let frame_half_width = cli.frame_half_width.or(cfg.frame_half_width).unwrap_or(0);

    // Boolean/operational flags: the CLI flag can only turn `all_frames` on; when
    // absent it falls back to the config value (default false).
    let all_frames = cli.all_frames || cfg.all_frames.unwrap_or(false);
    let threads = cli.threads.or(cfg.threads);

    if let Some(t) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .ok();
    }

    let stats = dnoise::denoise(
        &cli.input,
        &cli.output,
        &params,
        all_frames,
        frame_half_width,
        cli.force,
    )?;
    let pct = if stats.raw_points > 0 {
        100.0 * stats.kept_points as f64 / stats.raw_points as f64
    } else {
        0.0
    };
    println!(
        "dnoise: {} frames, {} -> {} points kept ({:.1}%)",
        stats.frames, stats.raw_points, stats.kept_points, pct
    );
    Ok(())
}
