# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **diaPASEF isolation-window MS/MS filtering** (`dia_window` module,
  `DiaWindowParams`, `in_window_mask`, `filter_per_window`; `tdf::read_dia_windows`
  joining `DiaFrameMsMsInfo` + `DiaFrameMsMsWindows`). Two opt-in, diaPASEF-only
  features sharing one window reader: `--dia-window` drops MS/MS points whose
  mobility scan falls outside every isolation window for their frame
  (out-of-window noise; `--dia-window-scan-pad` adds edge leniency), and
  `--dia-per-window` runs the MS/MS filter independently within each window's scan
  slice so a mobility run is never fused across a window boundary (no cross-talk
  between the unrelated precursor m/z bands adjacent windows isolate). Both are
  no-ops on ddaPASEF. `denoise`/`denoise_with_progress` gain
  `dia_window: Option<&DiaWindowParams>` and `dia_per_window: bool` arguments.
- ddaPASEF **MS/MS denoising** (`msms` module, `MsmsFilterParams`,
  `combine_and_filter`): for each precursor, combine its fragment scans across the
  frames it was isolated in, denoise the combined spectrum (vertical + halo), and
  prune the individual scans to the surviving `(scan, TOF)`. Opt-in via
  `--denoise-msms` (off by default) with separate `--msms-*` knobs;
  `denoise`/`denoise_with_progress` gain a `denoise_msms: Option<&MsmsFilterParams>`
  argument. Unlike MS1 denoising, this modifies MS/MS spectra and identifications.
- Horizontal-halo filter (`halo` module, `HaloParams`,
  `horizontal_halo_keep_mask`): removes the weak m/z halo flanking bright ions
  (left/right only) by comparing each peak to the max of its surrounding box
  excluding its own TOF column, in integer `(scan, TOF index)` space. Runs after
  the vertical filter, **on by default**; disable with `--no-halo`. Ported from
  `tdfpy`'s `HorizontalHaloFilter`. `denoise`/`denoise_with_progress` gain a
  `halo: Option<&HaloParams>` argument.
- Public library API with two tiers: high-level `denoise` /
  `denoise_with_progress`, and low-level `FlatFrame`, `filter_once`,
  `filter_iterated`, `running_average`, and the type-2 `codec` module.
- Typed error types `DnoiseError` and `DecodeError` (replacing `anyhow` in the
  library; `anyhow` is now used only by the CLI binary).
- `denoise_with_progress` with a `Progress` callback, so the library no longer
  depends on any terminal-UI crate.
- `cli` cargo feature (enabled by default) gating the binary and its
  `clap`/`indicatif`/`serde`/`toml`/`anyhow` dependencies. Build with
  `--no-default-features` for a library-only dependency tree.
- **Experimental** `--frame-half-width` / `frame_half_width`: pre-average each
  MS1 frame over its `2r+1` MS1-frame neighborhood before filtering. See the
  README for its current limitations.

### Changed
- The SQLite `tdf` plumbing is now crate-private; the type-2 codec moved to the
  public `dnoise::codec` module (`dnoise::tdf::encode::*` → `dnoise::codec::*`).
