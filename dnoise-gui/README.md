# dnoise-gui

A point-and-click desktop front end for [`dnoise`](../), for users who would
rather not touch the command line. It links the `dnoise` **library** directly
(no subprocess), so denoising runs in-process with live progress.

Built with [`egui`/`eframe`](https://github.com/emilk/egui) — one self-contained
native binary per OS, no web runtime.

The queue identifies prm-PASEF acquisitions automatically. PRM preserves fragments
by default; “Also denoise MS/MS” enables experimental filtering within each
isolation event. Target/event metadata is checked before a run enters the queue. See
[targeted support](../docs/targeted.md) for behavior and validation limits.

## Run / build

Building the GUI requires Rust **1.95+** (the CLI/library minimum remains 1.85).

```bash
cargo run -p dnoise-gui              # from the workspace root
cargo build -p dnoise-gui --release # release binary at target/release/dnoise-gui
```

## What it does

- **Input:** use *Browse…*, drag `.d` folders onto the window, or paste a path.
  Each folder's acquisition scheme (ddaPASEF / diaPASEF / MS1-only) is detected and
  shown in the list.
- **Preset:** pick *Auto-detect* (chooses the right MS1 gate per file), *ddaPASEF*,
  *diaPASEF*, or *None*.
- **Output:** write next to each input with a suffix (`_dnoise`) or into one chosen
  folder; optional overwrite and per-output JSON report.
- **Validation:** enabled by default; uncheck *Full input/output validation* in
  advanced settings to skip extra decoding passes. Structural checks remain, and
  the selected mode is recorded inside the output.
- **Advanced settings** (collapsible): every filter, halo, gate, crop, and ppm knob,
  defaulting to the tuned CLI defaults — leave it closed for a standard run, open it
  to tune aggressiveness, add a region-of-interest crop, switch the m/z window to a
  ppm tolerance, or override the preset's gates.
- **Estimate reduction:** dry-runs an 8% frame sample and reports "~X% kept" per file
  in seconds, without writing anything — tune, estimate, repeat, then Run.
- **Settings file:** Save/Load the advanced knobs as a `dnoise.toml` that the CLI
  reads too (same schema); or drag a `.toml` onto the window to load it.
- **Run:** processes the queue on a worker thread with a per-file progress bar and a
  log pane. *Cancel* discards only the current temporary output, preserving existing files.
- **Recovery:** output collisions are checked before a batch starts; *Retry unfinished*
  reruns unfinished queue items in the current session.
- **Inspect:** after completion, compare before/after frames and removed intensity
  using identical axes and color scales.
- **Save batch…:** export the queue and per-file settings for `dnoise batch jobs.json`.

Everything runs through the same `dnoise::denoise_with_options` path as the CLI.

## Downloads

Pre-built binaries for Linux, macOS, and Windows are attached to each tagged
[release](https://github.com/pgarrett-scripps/dnoise/releases) (built by
`.github/workflows/release.yml`). Each archive contains both `dnoise` (CLI) and
`dnoise-gui`. The binaries are currently **unsigned**, so first launch may need an
OS "allow anyway" step (Windows SmartScreen / macOS Gatekeeper).

## Processing records

Completed outputs include `dnoise.provenance.json` and `dnoise.config.toml` inside
the `.d` folder. Loaded TOML settings use the same resolver as the CLI, including
options beyond the visible controls. Invalid crop text blocks a run.
See [provenance](../docs/provenance.md) and [batch workflows](../docs/batch.md).

Signed installers remain future work. Native dialogs and platform-specific
launch behavior still need testing on each release platform.

Advanced settings include **Neighbor support (experimental)**: separate MS1, PRM
target, and DIA window radii, plus the maximum RT distance in seconds. Radius 1
uses the previous/next compatible observation for filtering decisions only.
PRM/DIA require MS/MS denoising. Defaults are off; native output intensity does
not establish quantitative accuracy. See [details](../docs/neighbors.md).
