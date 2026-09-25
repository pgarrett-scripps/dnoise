# dnoise-core

The in-memory denoising stages of [dnoise](https://crates.io/crates/dnoise), for
Bruker timsTOF frames, in pure Rust.

`dnoise-core` does no I/O. Its only dependencies are `serde` and `thiserror`: no
`timsrust`, no SQLite, no zstd, no rayon, and no crate with native code or a
`links` key. CI checks this (`.github/scripts/check-core-deps.sh`). A program
that already reads timsTOF data, such as a search engine with its own
`timsrust` and `rusqlite` versions, can embed it without dependency clashes.

The `dnoise` crate is the I/O layer: it reads `.d` folders, builds the run
state from `analysis.tdf`, calls this crate for every frame, and writes a new
`.d`. The CLI and the library go through the same per-frame call, so an embedder
gets the same output as `dnoise` for the same frames and run state.

## Usage

Build one `Denoiser` per run, then denoise frames one at a time:

```rust
use dnoise_core::{Acquisition, Denoiser, FilterParams, FlatFrame, FrameMeta, Stages};

// Frame metadata for the whole run, in `Frames.Id` order (here, one MS1 frame).
let meta = vec![FrameMeta { id: 1, num_scans: 700, num_peaks: 9, ms_ms_type: 0, rt: 0.5 }];
let params = FilterParams::default();
let stages = Stages::default();
let denoiser = Denoiser::new(&params, &stages, Acquisition::Ms1Only, meta, false)?;

// Decoded points of frame 0: a streak at TOF 1000 over 8 scans, and one noise point.
let mut scan: Vec<u32> = (10..18).collect();
scan.push(300);
let mut tof = vec![1000; 8];
tof.push(50_000);
let frame = FlatFrame { frame_id: 1, num_scans: 700, scan, tof, intensity: vec![100; 9] };

let out = denoiser.denoise_ms1(0, &frame)?;
assert_eq!(out.survivors.len(), 8); // (scan, tof, intensity) of the streak
# Ok::<(), dnoise_core::Error>(())
```

A frame read with `timsrust` converts with `FlatFrame::from_csr` (or implement
`CsrFrame` for the reader's frame type and use `FlatFrame::from_frame`).

## Run state and cross-frame state

`Denoiser::new` takes the run's frame metadata and the parameters. It applies
the acquisition policy (for example, no discovery gates on prm-PASEF), and
the resulting `Denoiser::stages()` say which optional state the run needs. Set
that state once, before the first frame:

| Setter | Contents | Source in a `.d` |
|---|---|---|
| `set_msms_keep` | ddaPASEF per-precursor keep sets | every MS/MS frame, via `msms::MsmsKeepBuilder` |
| `set_prm_windows` | prm-PASEF isolation events | `PrmFrameMsMsInfo` |
| `set_whole_frame_msms` | diaPASEF MS/MS filtering on | no `PasefFrameMsMsInfo` rows |
| `set_dia_windows`, `set_dia_regions` | diaPASEF isolation windows | `DiaFrameMsMsInfo`, `DiaFrameMsMsWindows` |
| `set_dda_windows` | ddaPASEF event intervals | `PasefFrameMsMsInfo` |
| `set_dia_ms1_gate`, `set_polygon_gate` | MS1 gates | windows or IMS polygon, plus calibration |
| `set_crop_gate`, `set_rt_crop` | region-of-interest crop | crop bounds, plus calibration |
| `set_neighbors` | temporal neighbor index | frame compatibility metadata |

After that the denoiser is read-only, so frames can be processed in any order
and on any thread. Two stages depend on other frames:

- **ddaPASEF MS/MS denoising** pools a precursor's scans over every frame that
  isolated it. `MsmsKeepBuilder` does this once, up front, over all MS/MS
  frames; each per-frame call only looks up the result.
- **Temporal neighbor support** (off by default) sums up to `radius`
  compatible frames either side of the frame. Use `Denoiser::process` with a
  frame source for this; `denoise_ms1` and `denoise_msms` return an error when a
  neighbor frame would be needed.

The order in which stages run on a frame is internal to `Denoiser` and may
change between releases. Callers never sequence stages themselves.

## License

MIT
