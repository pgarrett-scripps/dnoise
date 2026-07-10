# dnoise-gui

A point-and-click desktop front end for [`dnoise`](../), for users who would
rather not touch the command line. It links the `dnoise` **library** directly
(no subprocess), so denoising runs in-process with live progress.

Built with [`egui`/`eframe`](https://github.com/emilk/egui) — one self-contained
native binary per OS, no web runtime.

## Run / build

```bash
cargo run -p dnoise-gui              # from the workspace root
cargo build -p dnoise-gui --release # release binary at target/release/dnoise-gui
```

## What it does

- **Input:** drag `.d` folders onto the window, or paste a path and click *Add*.
  Each folder's acquisition scheme (ddaPASEF / diaPASEF / MS1-only) is detected and
  shown in the list.
- **Preset:** pick *Auto-detect* (chooses the right MS1 gate per file), *ddaPASEF*,
  *diaPASEF*, or *None*.
- **Output:** write next to each input with a suffix (`_dnoise`) or into one chosen
  folder; optional overwrite and per-output JSON report.
- **Advanced settings** (collapsible): every filter, halo, gate, crop, and ppm knob,
  defaulting to the tuned CLI defaults — leave it closed for a standard run, open it
  to tune aggressiveness, add a region-of-interest crop, switch the m/z window to a
  ppm tolerance, or override the preset's gates.
- **Estimate reduction:** dry-runs an 8% frame sample and reports "~X% kept" per file
  in seconds, without writing anything — tune, estimate, repeat, then Run.
- **Settings file:** Save/Load the advanced knobs as a `dnoise.toml` that the CLI
  reads too (same schema); or drag a `.toml` onto the window to load it.
- **Run:** processes the queue on a worker thread with a per-file progress bar and a
  log pane. *Cancel* stops the current file mid-way and removes its partial output.

Everything runs through the same `dnoise::denoise_with_options` path as the CLI.

## Downloads

Pre-built binaries for Linux, macOS, and Windows are attached to each tagged
[release](https://github.com/pgarrett-scripps/dnoise/releases) (built by
`.github/workflows/release.yml`). Each archive contains both `dnoise` (CLI) and
`dnoise-gui`. The binaries are currently **unsigned**, so first launch may need an
OS "allow anyway" step (Windows SmartScreen / macOS Gatekeeper).

## Roadmap (not yet implemented)

- **Native OS file pickers** (behind a platform-gated `rfd` dependency; drag-and-drop
  and path entry cover input for now).
- A before/after **frame preview** heatmap.
- **Signed** installers: Windows `.exe`, macOS `.app`/dmg, Linux AppImage.
