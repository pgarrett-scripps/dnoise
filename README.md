# dnoise

Denoise Bruker timsTOF `.d` folders by reducing the raw 3-D data
(`scan × TOF-index × intensity`) with the **iterative vertical-IM feature
filter** (Stage 1 of [ALGORITHM.md](ALGORITHM.md)). Centroiding and smoothing
are intentionally not implemented.

A real ion forms a vertical streak in `(TOF-index × scan)` space. The filter
walks each TOF index, sums the ion-mobility profile in a small TOF window
around it, and keeps only points belonging to long-enough, intense-enough
vertical runs. The pass is iterated: each pass re-filters the survivors of the
previous one.

## Usage

```bash
dnoise <INPUT.d> <OUTPUT.d> [options]
```

The source folder is never modified; a new `.d` is written with a rewritten
`analysis.tdf_bin` (re-encoded as compression **type 2**) and an updated
`analysis.tdf` (`Frames.TimsId/NumPeaks/MaxIntensity/SummedIntensities` and
`GlobalMetadata.TimsCompressionType`). The leading reserved header that Bruker
places before the first frame (the smallest `Frames.TimsId`, often 64 bytes) is
copied verbatim and all rewritten offsets are shifted past it, so the output is
byte-layout-compatible with the Bruker SDK / `timsdata` DLL.

| Option | Default | Meaning |
|---|---|---|
| `--mz-half-width` | 2 | Column half-width in TOF indices (`[c-w, c+w]`). |
| `--min-feature-length` | 5 | Minimum total span (scans) of a kept feature. |
| `--max-internal-gap` | 1 | Max empty scans tolerated inside a feature. |
| `--min-window-intensity` | 0 | Per-scan summed-intensity floor for occupancy. |
| `--min-feature-intensity` | 0 | Total summed-intensity floor for a kept feature. |
| `--iterations` | 1 | Filter passes (each re-applies to prior survivors). |
| `--all-frames` | off | Also filter MS/MS frames (default: MS1 only). |
| `--threads` | all cores | Worker threads. |
| `--config` / `-c` | — | Load parameters from a TOML file (see below). |
| `--force` | off | Overwrite an existing output folder. |

Filtering runs in parallel across frames (rayon); frames are written in order
so binary offsets stay consistent.

### Config file

Instead of (or alongside) flags, parameters can come from a TOML file:

```bash
dnoise <INPUT.d> <OUTPUT.d> --config dnoise.toml
```

Every key is optional and uses the same name as the flag with underscores. See
[dnoise.toml](dnoise.toml) for a fully-commented example. Precedence is
**explicit CLI flag > config file > built-in default**, so a config sets the
baseline and individual flags override it for one run. Unknown keys are rejected
to catch typos.

```toml
mz_half_width = 3
min_feature_length = 7
iterations = 2
all_frames = false
# threads = 8
```

## Validate output

```bash
cargo run --release --example validate -- <PATH.d>
```

Re-reads every frame with timsrust and checks each frame's decoded peak count
against `Frames.NumPeaks`.

To verify the type-2 codec against *real Bruker bytes* (decode raw frames
straight from `analysis.tdf_bin` and compare to timsrust):

```bash
cargo run --release --example check_codec -- <PATH.d> [num_frames]
```

## Limitation: input compression type

dnoise reads **compression type 2** (and uncompressed) input and always writes
type 2. It cannot yet read **type 3** (zstd + bitshuffle) `.d` files: timsrust's
`timscompress` feature depends on a `timscompress` crate that is only a stub on
crates.io, so the decoder does not build. Type-3 support is blocked on that
crate being published.

## Tests

```bash
cargo test     # codec round-trip + filter behavior
cargo clippy --all-targets
```

`tests/roundtrip.rs` gates the type-2 encoder by decoding what it encodes;
`tests/filter.rs` checks streak retention, gap-closing, window aggregation, and
iteration monotonicity.
