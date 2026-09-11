//! dnoise CLI.

mod commands;

use anyhow::{Context, Result};
use clap::Parser;
use dnoise::{RunOptions, SampleSpec};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use std::io::IsTerminal;
use std::path::PathBuf;
use std::time::Instant;
use tracing::info;

/// Denoise a Bruker timsTOF .d folder via the iterative vertical-IM feature filter.
#[derive(Parser)]
#[command(
    name = "dnoise",
    version,
    about,
    after_help = "Commands: dnoise validate INPUT.d | dnoise metadata INPUT.d | dnoise batch MANIFEST.json"
)]
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

    /// MS1 support: previous/next compatible MS1 observations for the filtering
    /// decision only; native output intensities are preserved (0 = off).
    #[arg(long = "ms1-neighbor-radius", visible_alias = "frame-half-width")]
    frame_half_width: Option<usize>,
    /// PRM support: previous/next matching target observations (0 = off).
    /// Requires --denoise-msms on PRM runs; ignored on other acquisition types.
    #[arg(long)]
    prm_neighbor_radius: Option<usize>,
    /// DIA support: previous/next matching isolation-window observations (0 = off).
    /// Requires --denoise-msms or --all-frames on DIA runs; forces event boundaries.
    #[arg(long)]
    dia_neighbor_radius: Option<usize>,
    /// Maximum distance from the current frame to a supporting neighbor, seconds
    /// (default 5). Applies to MS1, PRM and DIA; must be finite and positive.
    #[arg(long)]
    neighbor_max_rt_gap: Option<f64>,

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

    /// Denoise MS/MS (off by default): DDA precursor pooling, DIA window filtering,
    /// or experimental PRM filtering within each isolation event. Changes spectra
    /// and quantitative results. Unsupported for mixed or unknown acquisitions.
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
    /// ddaPASEF. Independent of the MS/MS streak filter. On by default whenever
    /// MS/MS frames are filtered (--denoise-msms / --all-frames); accepting this
    /// flag explicitly is a harmless no-op.
    #[arg(long)]
    dia_window: bool,
    /// Disable the diaPASEF MS/MS out-of-window gate (otherwise on whenever MS/MS
    /// frames are filtered).
    #[arg(long)]
    no_dia_window: bool,
    /// diaPASEF gate: scans of leniency added to each side of every isolation
    /// window before a point is treated as out-of-window.
    #[arg(long)]
    dia_window_scan_pad: Option<u32>,
    /// diaPASEF only: when the MS/MS filter runs (--denoise-msms or --all-frames),
    /// filter each isolation window's scan slice independently instead of the whole
    /// frame, so a mobility run cannot be fused across a window boundary. On by
    /// default whenever MS/MS frames are filtered; accepting this flag explicitly is
    /// a harmless no-op. No effect on ddaPASEF.
    #[arg(long)]
    dia_per_window: bool,
    /// Disable diaPASEF per-window MS/MS filtering, reverting to whole-frame (the
    /// MS/MS filter runs across the whole frame instead of window-by-window).
    #[arg(long)]
    no_dia_per_window: bool,

    /// ddaPASEF only: drop MS/MS points whose mobility scan falls outside every
    /// PasefFrameMsMsInfo isolation event for their frame. Standard timsTOF
    /// ddaPASEF acquisitions record no such points, so this is a guarantee
    /// rather than a reduction. No effect on diaPASEF. On by default whenever
    /// MS/MS frames are filtered (--denoise-msms / --all-frames); accepting this
    /// flag explicitly is a harmless no-op.
    #[arg(long)]
    dda_window: bool,
    /// Disable the ddaPASEF MS/MS out-of-window gate (otherwise on whenever
    /// MS/MS frames are filtered).
    #[arg(long)]
    no_dda_window: bool,
    /// ddaPASEF gate: scans of leniency added to each side of every isolation
    /// event before a point is treated as out-of-window.
    #[arg(long)]
    dda_window_scan_pad: Option<u32>,

    /// diaPASEF only: drop MS1 points whose (m/z, mobility) falls outside every
    /// isolation window (precursors that are never fragmented). Windows are padded
    /// per --dia-ms1-mz-pad / --dia-ms1-im-pad so edge precursors keep their full
    /// isotopic envelope. No effect on ddaPASEF. On by default; passing this flag
    /// explicitly is a harmless no-op.
    #[arg(long)]
    dia_ms1_window: bool,
    /// Disable the diaPASEF MS1 out-of-window gate (otherwise on by default).
    #[arg(long)]
    no_dia_ms1_window: bool,
    /// diaPASEF MS1 gate: m/z leniency added to each side of every window, in Da.
    #[arg(long)]
    dia_ms1_mz_pad: Option<f64>,
    /// diaPASEF MS1 gate: ion-mobility leniency added to each side, in 1/K0.
    #[arg(long)]
    dia_ms1_im_pad: Option<f64>,

    /// Drop MS1 points outside the run's ddaPASEF/PASEF selection polygon (the IMS
    /// PolygonFilter stored in analysis.tdf) — signal in never-selected precursor
    /// space. On by default; auto-detected: a no-op when the run stores no polygon
    /// or defines a diaPASEF window scheme. Passing this flag explicitly is a
    /// harmless no-op.
    #[arg(long)]
    ms1_polygon: bool,
    /// Disable the MS1 selection-polygon gate (otherwise on by default for
    /// ddaPASEF/PASEF).
    #[arg(long)]
    no_ms1_polygon: bool,
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

    /// Crop: keep only points at or above this m/z (Da). Applies to all frames.
    #[arg(long, value_name = "MZ")]
    mz_min: Option<f64>,
    /// Crop: keep only points at or below this m/z (Da). Applies to all frames.
    #[arg(long, value_name = "MZ")]
    mz_max: Option<f64>,
    /// Crop: keep only points at or above this ion mobility (1/K0).
    #[arg(long, value_name = "K0")]
    im_min: Option<f64>,
    /// Crop: keep only points at or below this ion mobility (1/K0).
    #[arg(long, value_name = "K0")]
    im_max: Option<f64>,
    /// Crop: keep only frames at or after this retention time (minutes); earlier
    /// frames are emitted empty (never deleted, so the frame axis stays valid).
    #[arg(long, value_name = "MIN")]
    rt_min: Option<f64>,
    /// Crop: keep only frames at or before this retention time (minutes).
    #[arg(long, value_name = "MIN")]
    rt_max: Option<f64>,
    /// Crop: drop points below this intensity.
    #[arg(long, value_name = "N")]
    min_intensity: Option<u32>,
    /// Crop: drop points above this intensity.
    #[arg(long, value_name = "N")]
    max_intensity: Option<u32>,
    /// Apply only the crop (--mz-*/--im-*/--rt-*/--*-intensity) and skip all
    /// denoising, so the output is a raw subset of the input. Requires a crop bound.
    #[arg(long)]
    crop_only: bool,

    /// Set the vertical filter's m/z window from a mass tolerance in ppm rather than
    /// raw TOF indices, converted at a reference m/z via the run calibration.
    /// Overrides --mz-half-width when set.
    #[arg(long, value_name = "PPM")]
    mz_ppm: Option<f64>,
    /// Reference m/z (Da) for --mz-ppm. Default: midpoint of the acquired m/z range.
    #[arg(long, value_name = "MZ")]
    mz_ppm_ref: Option<f64>,

    /// Estimate the reduction without writing any output: prints the stats (and the
    /// --report JSON if given) and leaves the output folder untouched.
    #[arg(long)]
    dry_run: bool,
    /// With --dry-run, process only this fraction (0 < f <= 1) of frames, chosen
    /// deterministically, for a fast estimate. Ignored without --dry-run.
    #[arg(long, value_name = "FRACTION")]
    sample: Option<f64>,
    /// Seed for --sample frame selection (deterministic; default 0).
    #[arg(long, value_name = "N", default_value_t = 0)]
    sample_seed: u64,
    /// Write a JSON run report (effective config + reduction stats) to this file.
    #[arg(long, value_name = "FILE")]
    report: Option<PathBuf>,

    /// Worker threads (default: all cores).
    #[arg(long)]
    threads: Option<usize>,
    /// Maximum frames per processing batch (default 2048); smaller uses less memory.
    #[arg(long)]
    frame_batch_size: Option<usize>,
    /// Skip full input/output decoding checks; structural and file-safety checks remain.
    #[arg(long, alias = "no-validation", conflicts_with = "validate")]
    skip_validation: bool,
    /// Enable full validation (default), overriding skip_validation in the config.
    #[arg(long, conflicts_with = "skip_validation")]
    validate: bool,
    /// Overwrite the output folder if it already exists.
    #[arg(long)]
    force: bool,
    /// Denoise the input folder in place: write to a temporary sibling folder and,
    /// on success, replace the input with recovery on installation failure. Omit the OUTPUT argument.
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

/// Map the `-v` (repeatable) / `-q` flags to the `dnoise` log level used when
/// `RUST_LOG` is unset: `quiet` wins and forces `warn`; otherwise `-v` steps up
/// info -> debug -> trace (saturating).
fn verbosity_level(verbose: u8, quiet: bool) -> &'static str {
    if quiet {
        "warn"
    } else {
        match verbose {
            0 => "info",
            1 => "debug",
            _ => "trace",
        }
    }
}

/// Install the stderr `tracing` subscriber. `RUST_LOG` (if set) wins; otherwise the
/// level comes from `-v`/`-q`: quiet=warn, default=info, -v=debug, -vv=trace. Only
/// the `dnoise` crate is set to that level (dependencies stay at `warn`) so the
/// output is the pipeline's own narration, not noise from libraries.
fn init_logging(verbose: u8, quiet: bool) {
    use tracing_subscriber::{EnvFilter, fmt};
    let level = verbosity_level(verbose, quiet);
    let filter = EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| EnvFilter::new(format!("dnoise={level},warn")));
    fmt()
        .with_env_filter(filter)
        .with_writer(std::io::stderr)
        .with_target(true)
        .init();
}

use dnoise::config::Config as FileConfig;
fn main() -> Result<()> {
    if let Some(result) = commands::dispatch() {
        return result;
    }
    let cli = Cli::parse();
    init_logging(cli.verbose, cli.quiet);

    let mut cfg = match &cli.config {
        Some(path) => FileConfig::load(path).map_err(anyhow::Error::msg)?,
        None => FileConfig::default(),
    };
    if let Some(value) = cli.mz_half_width {
        cfg.mz_half_width = Some(value);
    }
    if let Some(value) = cli.min_feature_length {
        cfg.min_feature_length = Some(value);
    }
    if let Some(value) = cli.max_internal_gap {
        cfg.max_internal_gap = Some(value);
    }
    if let Some(value) = cli.min_window_intensity {
        cfg.min_window_intensity = Some(value);
    }
    if let Some(value) = cli.min_feature_intensity {
        cfg.min_feature_intensity = Some(value);
    }
    if let Some(value) = cli.iterations {
        cfg.iterations = Some(value);
    }
    if let Some(value) = cli.frame_half_width {
        cfg.frame_half_width = None;
        cfg.ms1_neighbor_radius = Some(value);
    }
    if let Some(value) = cli.prm_neighbor_radius {
        cfg.prm_neighbor_radius = Some(value);
    }
    if let Some(value) = cli.dia_neighbor_radius {
        cfg.dia_neighbor_radius = Some(value);
    }
    if let Some(value) = cli.neighbor_max_rt_gap {
        cfg.neighbor_max_rt_gap = Some(value);
    }
    if let Some(value) = cli.halo_peak_fraction {
        cfg.halo_peak_fraction = Some(value);
    }
    if let Some(value) = cli.halo_mz_idx_half_width {
        cfg.halo_mz_idx_half_width = Some(value);
    }
    if let Some(value) = cli.halo_scan_half_width {
        cfg.halo_scan_half_width = Some(value);
    }
    if cli.denoise_msms {
        cfg.denoise_msms = Some(true);
    }
    if let Some(value) = cli.msms_mz_half_width {
        cfg.msms_mz_half_width = Some(value);
    }
    if let Some(value) = cli.msms_min_feature_length {
        cfg.msms_min_feature_length = Some(value);
    }
    if let Some(value) = cli.msms_max_internal_gap {
        cfg.msms_max_internal_gap = Some(value);
    }
    if let Some(value) = cli.msms_min_window_intensity {
        cfg.msms_min_window_intensity = Some(value);
    }
    if let Some(value) = cli.msms_min_feature_intensity {
        cfg.msms_min_feature_intensity = Some(value);
    }
    if let Some(value) = cli.msms_iterations {
        cfg.msms_iterations = Some(value);
    }
    if cli.smooth {
        cfg.smooth = Some(true);
    }
    if let Some(value) = cli.smooth_mz_idx_half_width {
        cfg.smooth_mz_idx_half_width = Some(value);
    }
    if let Some(value) = cli.smooth_scan_half_width {
        cfg.smooth_scan_half_width = Some(value);
    }
    if let Some(value) = cli.smooth_iterations {
        cfg.smooth_iterations = Some(value);
    }
    if cli.watershed {
        cfg.watershed = Some(true);
    }
    if let Some(value) = cli.watershed_box_scan {
        cfg.watershed_box_scan = Some(value);
    }
    if let Some(value) = cli.watershed_box_mz_idx {
        cfg.watershed_box_mz_idx = Some(value);
    }
    if let Some(value) = cli.watershed_min_seed_intensity {
        cfg.watershed_min_seed_intensity = Some(value);
    }
    if let Some(value) = cli.watershed_min_centroid_total {
        cfg.watershed_min_centroid_total = Some(value);
    }
    if let Some(value) = cli.watershed_max_tof_offset {
        cfg.watershed_max_tof_offset = Some(value);
    }
    if cli.box_centroid {
        cfg.box_centroid = Some(true);
    }
    if let Some(value) = cli.box_centroid_mz_idx_half {
        cfg.box_centroid_mz_idx_half = Some(value);
    }
    if let Some(value) = cli.box_centroid_scan_half {
        cfg.box_centroid_scan_half = Some(value);
    }
    if let Some(value) = cli.box_centroid_min_total {
        cfg.box_centroid_min_total = Some(value);
    }
    if cli.dia_window {
        cfg.dia_window = Some(true);
    }
    if cli.no_dia_window {
        cfg.dia_window = Some(false);
    }
    if let Some(value) = cli.dia_window_scan_pad {
        cfg.dia_window_scan_pad = Some(value);
    }
    if cli.dia_per_window {
        cfg.dia_per_window = Some(true);
    }
    if cli.no_dia_per_window {
        cfg.dia_per_window = Some(false);
    }
    if cli.dda_window {
        cfg.dda_window = Some(true);
    }
    if cli.no_dda_window {
        cfg.dda_window = Some(false);
    }
    if let Some(value) = cli.dda_window_scan_pad {
        cfg.dda_window_scan_pad = Some(value);
    }
    if cli.dia_ms1_window {
        cfg.dia_ms1_window = Some(true);
    }
    if cli.no_dia_ms1_window {
        cfg.dia_ms1_window = Some(false);
    }
    if let Some(value) = cli.dia_ms1_mz_pad {
        cfg.dia_ms1_mz_pad = Some(value);
    }
    if let Some(value) = cli.dia_ms1_im_pad {
        cfg.dia_ms1_im_pad = Some(value);
    }
    if cli.ms1_polygon {
        cfg.ms1_polygon = Some(true);
    }
    if cli.no_ms1_polygon {
        cfg.ms1_polygon = Some(false);
    }
    if let Some(value) = cli.ms1_polygon_mz_pad {
        cfg.ms1_polygon_mz_pad = Some(value);
    }
    if let Some(value) = cli.ms1_polygon_im_pad {
        cfg.ms1_polygon_im_pad = Some(value);
    }
    if let Some(value) = cli.mz_min {
        cfg.mz_min = Some(value);
    }
    if let Some(value) = cli.mz_max {
        cfg.mz_max = Some(value);
    }
    if let Some(value) = cli.im_min {
        cfg.im_min = Some(value);
    }
    if let Some(value) = cli.im_max {
        cfg.im_max = Some(value);
    }
    if let Some(value) = cli.rt_min {
        cfg.rt_min = Some(value);
    }
    if let Some(value) = cli.rt_max {
        cfg.rt_max = Some(value);
    }
    if let Some(value) = cli.min_intensity {
        cfg.min_intensity = Some(value);
    }
    if let Some(value) = cli.max_intensity {
        cfg.max_intensity = Some(value);
    }
    if cli.crop_only {
        cfg.crop_only = Some(true);
    }
    if let Some(value) = cli.mz_ppm {
        cfg.mz_ppm = Some(value);
    }
    if let Some(value) = cli.mz_ppm_ref {
        cfg.mz_ppm_ref = Some(value);
    }
    if cli.all_frames {
        cfg.all_frames = Some(true);
    }
    if let Some(value) = cli.threads {
        cfg.threads = Some(value);
    }
    if let Some(value) = cli.frame_batch_size {
        cfg.frame_batch_size = Some(value);
    }
    if cli.skip_validation {
        cfg.skip_validation = Some(true);
    } else if cli.validate {
        cfg.skip_validation = Some(false);
    }
    if cli.no_halo {
        cfg.halo = Some(false);
    }
    let resolved = cfg.resolve(&cli.input)?;
    let params = resolved.filter;
    let crop = resolved.crop;
    let crop_only = resolved.crop_only;
    let stages = resolved.stages();
    // Dry-run / sampling / report.
    let dry_run = cli.dry_run;
    if dry_run && cli.in_place {
        anyhow::bail!("--dry-run writes nothing, so it cannot be combined with --in-place");
    }
    let sample = match cli.sample {
        Some(f) => {
            if !dry_run {
                anyhow::bail!("--sample only applies with --dry-run");
            }
            if !(f > 0.0 && f <= 1.0) {
                anyhow::bail!("--sample must be a fraction in (0, 1], got {f}");
            }
            Some(SampleSpec {
                fraction: f,
                seed: cli.sample_seed,
            })
        }
        None => None,
    };

    // Resolve where to write. In-place mode writes to a temp sibling folder and
    // swaps it over the input on success; otherwise the OUTPUT argument is used. A
    // dry run writes nothing, so OUTPUT is optional and the path is only a
    // placeholder the pipeline never touches.
    let out_path = if cli.in_place {
        cli.input.clone()
    } else if dry_run {
        cli.output
            .clone()
            .unwrap_or_else(|| PathBuf::from("dnoise-dry-run-unused"))
    } else {
        cli.output
            .clone()
            .context("provide an OUTPUT folder or use --in-place")?
    };
    // Explicit in-place mode authorizes replacing the original after validation.
    let force = cli.force || cli.in_place;

    let options = RunOptions {
        force,
        dry_run,
        crop: (!crop.is_empty()).then_some(&crop),
        crop_only,
        sample,
        cancel: None,
        frame_batch_size: cfg.frame_batch_size,
        skip_validation: cfg.skip_validation.unwrap_or(false),
    };

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
    if let Some(report) = &cli.report {
        dnoise::output::check_disjoint(&cli.input, report)?;
        if !dry_run {
            dnoise::output::check_disjoint(&out_path, report)?;
        }
        if report.exists() {
            anyhow::bail!("report already exists: {}", report.display());
        }
    }
    if dnoise::provenance::read(&cli.input)?.is_some() {
        tracing::warn!("input has dnoise history; this run processes already modified data");
    }
    let start = Instant::now();
    let mut last_decile = 0u64;
    let progress = |p: dnoise::Progress| {
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
    };
    let stats = cfg.in_thread_pool(|| {
        if cli.in_place {
            dnoise::denoise_in_place(&cli.input, &params, &stages, &options, progress)
        } else {
            dnoise::denoise_with_options(
                &cli.input, &out_path, &params, &stages, &options, progress,
            )
        }
    })?;
    pb.finish_and_clear();
    let elapsed = start.elapsed();
    let pct = if stats.raw_points > 0 {
        100.0 * stats.kept_points as f64 / stats.raw_points as f64
    } else {
        0.0
    };

    // Optional JSON report: the effective config plus the reduction stats, for
    // parameter sweeps and provenance. Written for both real and dry runs.
    if let Some(path) = &cli.report {
        let report = dnoise::provenance::report(&cli.input, &params, &stages, &options, &stats);
        use std::io::Write;
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)?;
        file.write_all((serde_json::to_string_pretty(&report)? + "\n").as_bytes())?;
        info!(report = %path.display(), "wrote run report");
    }

    // Canonical result line on stdout (logs went to stderr): easy to grep/parse. A
    // dry run is flagged so the line is not mistaken for a written output.
    let tag = if stats.dry_run { " [dry-run]" } else { "" };
    let sampled = if stats.processed_frames != stats.frames {
        format!(" ({} sampled)", stats.processed_frames)
    } else {
        String::new()
    };
    println!(
        "dnoise:{} {} frames{}, {} -> {} points kept ({:.1}%) in {:.1}s",
        tag,
        stats.frames,
        sampled,
        stats.raw_points,
        stats.kept_points,
        pct,
        elapsed.as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::verbosity_level;

    #[test]
    fn verbosity_maps_v_flags_to_levels() {
        assert_eq!(verbosity_level(0, false), "info");
        assert_eq!(verbosity_level(1, false), "debug");
        assert_eq!(verbosity_level(2, false), "trace");
        assert_eq!(verbosity_level(9, false), "trace"); // saturates
    }

    #[test]
    fn quiet_overrides_verbose() {
        assert_eq!(verbosity_level(0, true), "warn");
        assert_eq!(verbosity_level(3, true), "warn");
    }
}
