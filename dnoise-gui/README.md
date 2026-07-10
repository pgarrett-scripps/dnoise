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

## What it does (Phase 1 MVP)

- **Input:** drag `.d` folders onto the window, or paste a path and click *Add*.
  Each folder's acquisition scheme (ddaPASEF / diaPASEF / MS1-only) is detected and
  shown in the list.
- **Preset:** pick *Auto-detect* (chooses the right MS1 gate per file), *ddaPASEF*,
  *diaPASEF*, or *None*.
- **Output:** write next to each input with a suffix (`_dnoise`) or into one chosen
  folder; optional overwrite and per-output JSON report.
- **Run:** processes the queue on a worker thread with a per-file progress bar and a
  log pane. *Cancel* stops after the current file.

The MVP uses the tuned default filter parameters plus the preset's gate. Everything
runs through the same `dnoise::denoise_with_options` path as the CLI.

## Roadmap (not yet implemented)

- **Advanced settings** panel exposing every filter/gate/crop/ppm knob, with
  TOML config load & save (shared with the CLI).
- **Estimate reduction** button (dry-run + frame sample) for instant feedback
  before a full run.
- **Native OS file pickers** (behind a platform-gated `rfd` dependency; drag-and-drop
  and path entry cover input for now).
- **Mid-run cancellation** (needs a cooperative abort hook in the library) and a
  before/after **frame preview** heatmap.
- **Packaging:** signed Windows `.exe`, macOS `.app`/dmg, Linux AppImage via CI.
