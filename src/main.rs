//! dnoise CLI.

use anyhow::{Context, Result};
use clap::Parser;
use dnoise::{
    BoxCentroidParams, DiaMs1WindowParams, DiaWindowParams, FilterParams, HaloParams,
    Ms1PolygonParams, MsmsFilterParams, SmoothParams, Stages, WatershedParams,
};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use serde::Deserialize;
use std::io::IsTerminal;
use std::path::{Path, PathBuf};
use std::time::Instant;
use tracing::{info, warn};

/// Denoise a Bruker timsTOF .d folder via the iterative vertical-IM feature filter.
#[derive(Parser)]
#[command(name = "dnoise", version, about)]
struct Cli {
    /// Input Bruker .d folder.
    input: PathBuf,
    /// Output .d folder (created; must not exist unless --force). Omit when using
    /// --in-place.
    output: Option<PathBuf>,

    /// TOML config file with filter parameters. Explicit CLI flags override its values.
    #[arg(long, short = 'c', value_name = "FILE")]
    config: Option<PathBuf>,

    /// Column half-width in TOF indices.
    #[arg(long)]
    mz_half_width: Option<u32>,
    /// Minimum number of occupied scans in a kept feature (bridged gaps not counted).
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

    /// Disable the horizontal-halo filter (on by default), which removes the weak
    /// m/z halo flanking bright ions (left/right) after the vertical filter.
    #[arg(long)]
    no_halo: bool,
    /// Halo: drop a peak below this fraction of its off-column box-max reference.
    #[arg(long)]
    halo_peak_fraction: Option<f64>,
    /// Halo: reference-box half-width along TOF index.
    #[arg(long)]
    halo_mz_idx_half_width: Option<u32>,
    /// Halo: reference-box half-width along ion-mobility scan.
    #[arg(long)]
    halo_scan_half_width: Option<usize>,

    /// Denoise ddaPASEF MS/MS frames precursor-by-precursor (off by default):
    /// combine each precursor's fragment scans across frames, filter the combined
    /// spectrum, and prune the individual scans. Changes MS/MS spectra (and IDs).
    #[arg(long)]
    denoise_msms: bool,
    /// MS/MS filter: column half-width in TOF indices.
    #[arg(long)]
    msms_mz_half_width: Option<u32>,
    /// MS/MS filter: minimum number of occupied scans in a kept run.
    #[arg(long)]
    msms_min_feature_length: Option<usize>,
    /// MS/MS filter: max empty scans tolerated inside a feature.
    #[arg(long)]
    msms_max_internal_gap: Option<usize>,
    /// MS/MS filter: per-scan summed-intensity floor.
    #[arg(long)]
    msms_min_window_intensity: Option<u64>,
    /// MS/MS filter: total summed-intensity floor for a kept feature.
    #[arg(long)]
    msms_min_feature_intensity: Option<u64>,
    /// MS/MS filter: passes over the combined spectrum's survivors.
    #[arg(long)]
    msms_iterations: Option<usize>,

    /// Box-average point intensities after the halo filter (off by default):
    /// each surviving point's intensity is replaced by the mean over its
    /// (scan, TOF-index) box. Stabilises the watershed centroider against noise.
    #[arg(long)]
    smooth: bool,
    /// Smoothing: averaging-box half-width along TOF index.
    #[arg(long)]
    smooth_mz_idx_half_width: Option<u32>,
    /// Smoothing: averaging-box half-width along ion-mobility scan.
    #[arg(long)]
    smooth_scan_half_width: Option<usize>,
    /// Smoothing: passes over the smoother's own output.
    #[arg(long)]
    smooth_iterations: Option<usize>,

    /// Centroid each filtered frame's survivors with the watershed centroider
    /// as a final stage (off by default). Lossy: collapses groups of raw points
    /// into intensity-weighted centroids, typically shrinking the point count to
    /// a small fraction. Applied to the same frames the vertical filter touches.
    #[arg(long)]
    watershed: bool,
    /// Watershed: nearest-neighbour reach along the ion-mobility scan axis.
    #[arg(long)]
    watershed_box_scan: Option<u32>,
    /// Watershed: nearest-neighbour reach along the TOF-index axis.
    #[arg(long)]
    watershed_box_mz_idx: Option<u32>,
    /// Watershed: minimum intensity for a point to open a new group.
    #[arg(long)]
    watershed_min_seed_intensity: Option<u64>,
    /// Watershed: drop groups whose summed intensity is below this.
    #[arg(long)]
    watershed_min_centroid_total: Option<u64>,
    /// Watershed: max follower distance from the group seed, in TOF indices.
    #[arg(long)]
    watershed_max_tof_offset: Option<u32>,

    /// Final stage: greedy small-box centroiding (off by default). Consolidates
    /// points within small fixed (scan, TOF-index) boxes into intensity-weighted
    /// centroids — tiling mobility streaks rather than collapsing them (cf.
    /// --watershed). Mutually exclusive with --watershed.
    #[arg(long)]
    box_centroid: bool,
    /// Box-centroid: box half-width along TOF index (m/z); keep tight.
    #[arg(long)]
    box_centroid_mz_idx_half: Option<u32>,
    /// Box-centroid: box half-width along ion-mobility scan.
    #[arg(long)]
    box_centroid_scan_half: Option<u32>,
    /// Box-centroid: drop boxes whose summed intensity is below this.
    #[arg(long)]
    box_centroid_min_total: Option<u64>,

    /// diaPASEF only: drop MS/MS points whose mobility scan falls outside every
    /// isolation window for their frame (out-of-window noise). No effect on
    /// ddaPASEF. Independent of the MS/MS streak filter.
    #[arg(long)]
    dia_window: bool,
    /// diaPASEF gate: scans of leniency added to each side of every isolation
    /// window before a point is treated as out-of-window.
    #[arg(long)]
    dia_window_scan_pad: Option<u32>,
    /// diaPASEF only: when the MS/MS filter runs (--denoise-msms or --all-frames),
    /// filter each isolation window's scan slice independently instead of the whole
    /// frame, so a mobility run cannot be fused across a window boundary. No effect
    /// on ddaPASEF.
    #[arg(long)]
    dia_per_window: bool,

    /// diaPASEF only: drop MS1 points whose (m/z, mobility) falls outside every
    /// isolation window (precursors that are never fragmented). Windows are padded
    /// per --dia-ms1-mz-pad / --dia-ms1-im-pad so edge precursors keep their full
    /// isotopic envelope. No effect on ddaPASEF.
    #[arg(long)]
    dia_ms1_window: bool,
    /// diaPASEF MS1 gate: m/z leniency added to each side of every window, in Da.
    #[arg(long)]
    dia_ms1_mz_pad: Option<f64>,
    /// diaPASEF MS1 gate: ion-mobility leniency added to each side, in 1/K0.
    #[arg(long)]
    dia_ms1_im_pad: Option<f64>,

    /// Drop MS1 points outside the run's ddaPASEF/PASEF selection polygon (the IMS
    /// PolygonFilter stored in analysis.tdf) — signal in never-selected precursor
    /// space. Auto-detected: no-op when the run stores no polygon.
    #[arg(long)]
    ms1_polygon: bool,
    /// MS1 polygon gate: m/z leniency added to each side, in Da (keeps an edge
    /// precursor's isotopic envelope).
    #[arg(long)]
    ms1_polygon_mz_pad: Option<f64>,
    /// MS1 polygon gate: ion-mobility leniency added to each side, in 1/K0.
    #[arg(long)]
    ms1_polygon_im_pad: Option<f64>,

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
    /// Denoise the input folder in place: write to a temporary sibling folder and,
    /// on success, atomically replace the input with it. Omit the OUTPUT argument.
    #[arg(long, conflicts_with = "output")]
    in_place: bool,

    /// Increase log verbosity: -v adds debug detail, -vv adds trace. Logs go to
    /// stderr; stdout carries only the final result line. Overridden by RUST_LOG.
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,
    /// Quiet: log only warnings and errors. Overridden by RUST_LOG.
    #[arg(short, long, conflicts_with = "verbose")]
    quiet: bool,
}

/// Install the stderr `tracing` subscriber. `RUST_LOG` (if set) wins; otherwise the
/// level comes from `-v`/`-q`: quiet=warn, default=info, -v=debug, -vv=trace. Only
/// the `dnoise` crate is set to that level (dependencies stay at `warn`) so the
/// output is the pipeline's own narration, not noise from libraries.
fn init_logging(verbose: u8, quiet: bool) {
    use tracing_subscriber::{EnvFilter, fmt};
    let level = if quiet {
        "warn"
    } else {
        match verbose {
            0 => "info",
            1 => "debug",
            _ => "trace",
        }
    };
    let filter = EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| EnvFilter::new(format!("dnoise={level},warn")));
    fmt()
        .with_env_filter(filter)
        .with_writer(std::io::stderr)
        .with_target(true)
        .init();
}

/// Build a sibling path by appending `suffix` to `path`'s full name (so an
/// `input.d` folder yields e.g. `input.d.dnoise-tmp`).
fn sibling_with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
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
    halo: Option<bool>,
    halo_peak_fraction: Option<f64>,
    halo_mz_idx_half_width: Option<u32>,
    halo_scan_half_width: Option<usize>,
    denoise_msms: Option<bool>,
    msms_mz_half_width: Option<u32>,
    msms_min_feature_length: Option<usize>,
    msms_max_internal_gap: Option<usize>,
    msms_min_window_intensity: Option<u64>,
    msms_min_feature_intensity: Option<u64>,
    msms_iterations: Option<usize>,
    smooth: Option<bool>,
    smooth_mz_idx_half_width: Option<u32>,
    smooth_scan_half_width: Option<usize>,
    smooth_iterations: Option<usize>,
    watershed: Option<bool>,
    watershed_box_scan: Option<u32>,
    watershed_box_mz_idx: Option<u32>,
    watershed_min_seed_intensity: Option<u64>,
    watershed_min_centroid_total: Option<u64>,
    watershed_max_tof_offset: Option<u32>,
    box_centroid: Option<bool>,
    box_centroid_mz_idx_half: Option<u32>,
    box_centroid_scan_half: Option<u32>,
    box_centroid_min_total: Option<u64>,
    dia_window: Option<bool>,
    dia_window_scan_pad: Option<u32>,
    dia_per_window: Option<bool>,
    dia_ms1_window: Option<bool>,
    dia_ms1_mz_pad: Option<f64>,
    dia_ms1_im_pad: Option<f64>,
    ms1_polygon: Option<bool>,
    ms1_polygon_mz_pad: Option<f64>,
    ms1_polygon_im_pad: Option<f64>,
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

/// Resolve one knob with explicit CLI flag > config-file value > built-in default
/// precedence. `$cli` and `$cfg` are `Option`s; `$default` is the fallback value.
macro_rules! pick {
    ($cli:expr, $cfg:expr, $default:expr $(,)?) => {
        $cli.or($cfg).unwrap_or($default)
    };
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    init_logging(cli.verbose, cli.quiet);

    let cfg = match &cli.config {
        Some(path) => FileConfig::load(path)?,
        None => FileConfig::default(),
    };

    // Each knob below resolves via `pick!`: CLI flag > config file > built-in default.
    let d = FilterParams::default();
    let params = FilterParams {
        mz_half_width: pick!(cli.mz_half_width, cfg.mz_half_width, d.mz_half_width),
        min_feature_length: pick!(
            cli.min_feature_length,
            cfg.min_feature_length,
            d.min_feature_length
        ),
        max_internal_gap: pick!(
            cli.max_internal_gap,
            cfg.max_internal_gap,
            d.max_internal_gap
        ),
        min_window_intensity: pick!(
            cli.min_window_intensity,
            cfg.min_window_intensity,
            d.min_window_intensity
        ),
        min_feature_intensity: pick!(
            cli.min_feature_intensity,
            cfg.min_feature_intensity,
            d.min_feature_intensity
        ),
        num_iterations: pick!(cli.iterations, cfg.iterations, d.num_iterations),
    };

    // Pre-filter smoothing radius (decoupled from the filter knobs above): explicit
    // CLI flag > config file > 0 (off).
    let frame_half_width = cli.frame_half_width.or(cfg.frame_half_width).unwrap_or(0);

    // Horizontal-halo filter: on by default; `--no-halo` (or `halo = false` in the
    // config) disables it. Its knobs follow the same CLI > config > default.
    let halo_enabled = if cli.no_halo {
        false
    } else {
        cfg.halo.unwrap_or(true)
    };
    let d = HaloParams::default();
    let halo = HaloParams {
        peak_fraction: pick!(
            cli.halo_peak_fraction,
            cfg.halo_peak_fraction,
            d.peak_fraction
        ),
        mz_idx_half_width: pick!(
            cli.halo_mz_idx_half_width,
            cfg.halo_mz_idx_half_width,
            d.mz_idx_half_width
        ),
        scan_half_width: pick!(
            cli.halo_scan_half_width,
            cfg.halo_scan_half_width,
            d.scan_half_width
        ),
    };

    // MS/MS denoising (ddaPASEF): off unless --denoise-msms or `denoise_msms = true`.
    let msms_enabled = cli.denoise_msms || cfg.denoise_msms.unwrap_or(false);
    let d = MsmsFilterParams::default();
    let msms = MsmsFilterParams {
        mz_half_width: pick!(
            cli.msms_mz_half_width,
            cfg.msms_mz_half_width,
            d.mz_half_width
        ),
        min_feature_length: pick!(
            cli.msms_min_feature_length,
            cfg.msms_min_feature_length,
            d.min_feature_length
        ),
        max_internal_gap: pick!(
            cli.msms_max_internal_gap,
            cfg.msms_max_internal_gap,
            d.max_internal_gap
        ),
        min_window_intensity: pick!(
            cli.msms_min_window_intensity,
            cfg.msms_min_window_intensity,
            d.min_window_intensity
        ),
        min_feature_intensity: pick!(
            cli.msms_min_feature_intensity,
            cfg.msms_min_feature_intensity,
            d.min_feature_intensity
        ),
        num_iterations: pick!(cli.msms_iterations, cfg.msms_iterations, d.num_iterations),
    };

    // Box-averaging smoother (post-halo, pre-watershed): off unless --smooth or
    // `smooth = true`. CLI > config > default for its knobs.
    let smooth_enabled = cli.smooth || cfg.smooth.unwrap_or(false);
    let d = SmoothParams::default();
    let smooth = SmoothParams {
        mz_idx_half_width: pick!(
            cli.smooth_mz_idx_half_width,
            cfg.smooth_mz_idx_half_width,
            d.mz_idx_half_width
        ),
        scan_half_width: pick!(
            cli.smooth_scan_half_width,
            cfg.smooth_scan_half_width,
            d.scan_half_width
        ),
        iterations: pick!(cli.smooth_iterations, cfg.smooth_iterations, d.iterations),
    };

    // Watershed centroider (final stage): off unless --watershed or `watershed =
    // true`. Its knobs follow the same CLI > config > default precedence.
    let watershed_enabled = cli.watershed || cfg.watershed.unwrap_or(false);
    let d = WatershedParams::default();
    let watershed = WatershedParams {
        box_scan: pick!(cli.watershed_box_scan, cfg.watershed_box_scan, d.box_scan),
        box_mz_idx: pick!(
            cli.watershed_box_mz_idx,
            cfg.watershed_box_mz_idx,
            d.box_mz_idx
        ),
        min_seed_intensity: pick!(
            cli.watershed_min_seed_intensity,
            cfg.watershed_min_seed_intensity,
            d.min_seed_intensity
        ),
        min_centroid_total: pick!(
            cli.watershed_min_centroid_total,
            cfg.watershed_min_centroid_total,
            d.min_centroid_total
        ),
        max_tof_offset: pick!(
            cli.watershed_max_tof_offset,
            cfg.watershed_max_tof_offset,
            d.max_tof_offset
        ),
    };

    // Greedy small-box centroider (alternative final stage): off unless
    // --box-centroid or `box_centroid = true`. Mutually exclusive with watershed.
    let box_centroid_enabled = cli.box_centroid || cfg.box_centroid.unwrap_or(false);
    if box_centroid_enabled && watershed_enabled {
        anyhow::bail!(
            "--box-centroid and --watershed are mutually exclusive (both are terminal centroiders)"
        );
    }
    let d = BoxCentroidParams::default();
    let box_centroid = BoxCentroidParams {
        mz_idx_half_width: pick!(
            cli.box_centroid_mz_idx_half,
            cfg.box_centroid_mz_idx_half,
            d.mz_idx_half_width
        ),
        scan_half_width: pick!(
            cli.box_centroid_scan_half,
            cfg.box_centroid_scan_half,
            d.scan_half_width
        ),
        min_centroid_total: pick!(
            cli.box_centroid_min_total,
            cfg.box_centroid_min_total,
            d.min_centroid_total
        ),
    };

    // diaPASEF isolation-window features (both off by default, diaPASEF-only).
    // `dia_window` gates out-of-window MS/MS points; `dia_per_window` makes the
    // MS/MS filter run window-by-window. Both share the DiaFrameMsMs* tables.
    let dia_window_enabled = cli.dia_window || cfg.dia_window.unwrap_or(false);
    let d = DiaWindowParams::default();
    let dia_window = DiaWindowParams {
        scan_pad: pick!(cli.dia_window_scan_pad, cfg.dia_window_scan_pad, d.scan_pad),
    };
    let dia_per_window = cli.dia_per_window || cfg.dia_per_window.unwrap_or(false);

    // diaPASEF MS1 out-of-window gate (off by default, diaPASEF-only): drop MS1
    // points outside every isolation window's (m/z, mobility) region, padded in
    // physical units. CLI > config > default for the pads.
    let dia_ms1_enabled = cli.dia_ms1_window || cfg.dia_ms1_window.unwrap_or(false);
    let d = DiaMs1WindowParams::default();
    let dia_ms1 = DiaMs1WindowParams {
        mz_pad: pick!(cli.dia_ms1_mz_pad, cfg.dia_ms1_mz_pad, d.mz_pad),
        im_pad: pick!(cli.dia_ms1_im_pad, cfg.dia_ms1_im_pad, d.im_pad),
    };

    // MS1 selection-polygon gate (off by default): drop MS1 points outside the
    // run's IMS PolygonFilter. CLI > config > default for the pads. Auto-detects
    // polygon presence (no-op otherwise).
    let ms1_polygon_enabled = cli.ms1_polygon || cfg.ms1_polygon.unwrap_or(false);
    let d = Ms1PolygonParams::default();
    let ms1_polygon = Ms1PolygonParams {
        mz_pad: pick!(cli.ms1_polygon_mz_pad, cfg.ms1_polygon_mz_pad, d.mz_pad),
        im_pad: pick!(cli.ms1_polygon_im_pad, cfg.ms1_polygon_im_pad, d.im_pad),
    };

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

    let stages = Stages {
        filter_all_frames: all_frames,
        frame_half_width,
        halo: halo_enabled.then_some(&halo),
        denoise_msms: msms_enabled.then_some(&msms),
        smooth: smooth_enabled.then_some(&smooth),
        watershed: watershed_enabled.then_some(&watershed),
        box_centroid: box_centroid_enabled.then_some(&box_centroid),
        dia_window: dia_window_enabled.then_some(&dia_window),
        dia_per_window,
        dia_ms1: dia_ms1_enabled.then_some(&dia_ms1),
        ms1_polygon: ms1_polygon_enabled.then_some(&ms1_polygon),
    };

    // Echo the effective configuration so a caller can confirm exactly which knobs
    // and stages this run used (after config-file + CLI + default resolution) — the
    // single most useful thing for reproducing or debugging a run.
    info!(
        mz_half_width = params.mz_half_width,
        min_feature_length = params.min_feature_length,
        max_internal_gap = params.max_internal_gap,
        min_window_intensity = params.min_window_intensity,
        min_feature_intensity = params.min_feature_intensity,
        iterations = params.num_iterations,
        "config: vertical filter"
    );
    info!(
        all_frames,
        frame_half_width,
        halo = halo_enabled,
        denoise_msms = msms_enabled,
        smooth = smooth_enabled,
        watershed = watershed_enabled,
        box_centroid = box_centroid_enabled,
        dia_window = dia_window_enabled,
        dia_per_window,
        dia_ms1 = dia_ms1_enabled,
        ms1_polygon = ms1_polygon_enabled,
        "config: enabled stages"
    );

    // Resolve where to write. In-place mode writes to a temp sibling folder and
    // swaps it over the input on success; otherwise the OUTPUT argument is used.
    let out_path = if cli.in_place {
        sibling_with_suffix(&cli.input, ".dnoise-tmp")
    } else {
        cli.output
            .clone()
            .context("provide an OUTPUT folder or use --in-place")?
    };
    // The temp folder is ours to clobber, so force-overwrite any stale leftover.
    let force = cli.force || cli.in_place;

    // Progress rendering adapts to the output: an interactive bar when stderr is a
    // terminal, otherwise periodic log lines at ~10% steps so piped/captured output
    // (e.g. an agent reading the logs) stays clean and line-oriented instead of
    // filling with carriage-return bar redraws.
    let interactive = std::io::stderr().is_terminal();
    let pb = ProgressBar::new(0);
    if interactive {
        pb.set_style(ProgressStyle::with_template("{bar:40} {pos}/{len} frames").unwrap());
    } else {
        pb.set_draw_target(ProgressDrawTarget::hidden());
    }
    let start = Instant::now();
    let mut last_decile = 0u64;
    let result = dnoise::denoise_with_progress(&cli.input, &out_path, &params, &stages, force, |p| {
        pb.set_length(p.frames_total as u64);
        pb.set_position(p.frames_done as u64);
        if !interactive && p.frames_total > 0 {
            let pct = 100 * p.frames_done as u64 / p.frames_total as u64;
            let decile = pct - pct % 10;
            if decile > last_decile {
                last_decile = decile;
                info!(
                    frames_done = p.frames_done,
                    frames_total = p.frames_total,
                    pct,
                    "denoise: progress"
                );
            }
        }
    });
    pb.finish_and_clear();
    // In-place mode owns the temp sibling folder, so clean up a partial one on
    // failure rather than leaving it behind for the next run to clobber.
    let stats = match result {
        Ok(stats) => stats,
        Err(e) => {
            if cli.in_place {
                std::fs::remove_dir_all(&out_path).ok();
            }
            return Err(e.into());
        }
    };

    // In-place swap: move the original aside, install the denoised folder under the
    // input's name, then drop the backup. On a failed install, restore the original.
    if cli.in_place {
        let backup = sibling_with_suffix(&cli.input, ".dnoise-old");
        if backup.exists() {
            std::fs::remove_dir_all(&backup)
                .with_context(|| format!("removing stale backup {}", backup.display()))?;
        }
        std::fs::rename(&cli.input, &backup)
            .with_context(|| format!("moving original {} aside", cli.input.display()))?;
        if let Err(e) = std::fs::rename(&out_path, &cli.input) {
            std::fs::rename(&backup, &cli.input).ok();
            return Err(anyhow::Error::new(e).context(format!(
                "installing denoised folder at {}; original restored",
                cli.input.display()
            )));
        }
        // The swap already succeeded, so the denoised folder is correctly installed
        // at the input path. Failing to remove the backup is not a failure of the
        // operation — warn and leave it for manual cleanup instead of reporting the
        // whole run as failed (which would misleadingly imply the data is bad).
        if let Err(e) = std::fs::remove_dir_all(&backup) {
            warn!(
                input = %cli.input.display(),
                backup = %backup.display(),
                error = %e,
                "denoised folder installed, but the backup could not be removed (leftover on disk)"
            );
        }
    }
    let elapsed = start.elapsed();
    let pct = if stats.raw_points > 0 {
        100.0 * stats.kept_points as f64 / stats.raw_points as f64
    } else {
        0.0
    };
    // Canonical result line on stdout (logs went to stderr): easy to grep/parse.
    println!(
        "dnoise: {} frames, {} -> {} points kept ({:.1}%) in {:.1}s",
        stats.frames,
        stats.raw_points,
        stats.kept_points,
        pct,
        elapsed.as_secs_f64()
    );
    Ok(())
}
