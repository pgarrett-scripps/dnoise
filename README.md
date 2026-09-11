# dnoise

[![CI](https://github.com/pgarrett-scripps/dnoise/actions/workflows/ci.yml/badge.svg)](https://github.com/pgarrett-scripps/dnoise/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/dnoise.svg)](https://crates.io/crates/dnoise)
[![Docs.rs](https://docs.rs/dnoise/badge.svg)](https://docs.rs/dnoise)
[![License](https://img.shields.io/crates/l/dnoise.svg)](#license)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21959649.svg)](https://doi.org/10.5281/zenodo.21959649)

Denoise Bruker timsTOF `.d` folders. Real ions form vertical streaks along the
ion-mobility axis. Chemical and electronic noise is short, isolated, or
scattered. dnoise keeps the streaks and drops the rest, writing a cleaned `.d`
that stays drop-in compatible with the Bruker SDK and existing search tools.

![A real timsTOF MS1 frame: kept signal forms streaks, discarded noise is scattered](docs/graphical_abstract.png)

Across 72 ddaPASEF + diaPASEF benchmark runs: **35-53% smaller** native
binaries, LFQ accuracy preserved, and at most a 2.3% change in
identifications. All of those runs used the same default parameters with no
per-run tuning. They come from one instrument and two gradients, though, so
other instruments and sample types are untested. Validate on your own data
before committing: `--dry-run` reports the reduction without writing
anything. MS1-only mode (the default) preserves MS/MS spectra; validate downstream
identification and quantification results for your workflow.

## Install

```bash
cargo install dnoise
```

Or download a prebuilt binary (CLI + GUI, Linux/macOS/Windows) from the
[releases page](https://github.com/pgarrett-scripps/dnoise/releases).

## Usage

The defaults are the configuration benchmarked in the paper, so no flags are
needed:

```bash
dnoise input.d output.d
```

```text
INFO dnoise::writer: denoise: frame inventory scheme="ddaPASEF" frames=8639 ms1=786 msms=7853
INFO dnoise::writer: MS1 selection-polygon gate active
INFO dnoise::writer: denoise: complete frames=8639 raw_points=300509979 kept_points=110092035 kept_pct=36.64
```

By default only **MS1** frames are filtered, so MS/MS spectra are
untouched. Acquisition-aware gates detect whether the
run is ddaPASEF or diaPASEF and apply the matching geometry automatically. On
runs where a gate's geometry is absent it is a silent no-op.

Useful variations:

| Command | What it does |
|---|---|
| `dnoise in.d --in-place` | Replace input after validation, with recovery on installation failure. |
| `dnoise in.d out.d --dry-run` | Report the reduction without writing anything. |
| `dnoise in.d out.d --denoise-msms` | Also denoise MS/MS spectra (changes IDs, re-search to measure). |
| `dnoise in.d out.d --config my.toml` | Load parameters from a TOML file ([example](dnoise.toml)). |
| `dnoise in.d out.d --report run.json` | Write effective config + reduction stats as JSON. |
| `dnoise in.d out.d --skip-validation` | Skip full input/output decoding checks; keep structural and file-safety checks. |
| `dnoise validate out.d` | Check metadata and decode every output frame. |
| `dnoise metadata out.d` | Read the processing history stored inside the folder. |
| `dnoise batch jobs.json` | Process a portable batch manifest. |

Every completed output includes **`dnoise.provenance.json`** (version, exact
settings, statistics, and processing history) and **`dnoise.config.toml`** (a
reusable recipe). These files travel with the `.d` folder. See
[processing metadata](docs/provenance.md) and [batch workflows](docs/batch.md).

Every knob (filter parameters, per-gate control, region-of-interest cropping,
smoothing and centroiding stages, logging) is documented in the
**[full reference](docs/reference.md)** and in `dnoise --help`. The method
itself is described in [ALGORITHM.md](ALGORITHM.md).

## Library usage

dnoise is also a Rust library ([docs.rs](https://docs.rs/dnoise)). Depend on it
without the CLI's dependencies:

```toml
[dependencies]
dnoise = { version = "0.1", default-features = false }
```

```rust,no_run
use dnoise::{FilterParams, Stages, denoise};
use std::path::Path;

let stats = denoise(
    Path::new("input.d"),
    Path::new("output.d"),
    &FilterParams::default(),
    &Stages::default(), // optional stages (halo, gates, smoothing, centroiders); all off
    false,              // don't overwrite an existing output
)?;
println!("{} -> {} points", stats.raw_points, stats.kept_points);
# Ok::<(), dnoise::DnoiseError>(())
```

A lower-level API exposes the filter on in-memory frames (`FlatFrame`,
`filter_iterated`) and the type-2 codec directly. See
[docs.rs](https://docs.rs/dnoise) and [docs/reference.md](docs/reference.md).

## Neighbor support (experimental)

Use nearby matching observations to support weak signals during filtering:

```bash
dnoise input.d output.d --denoise-msms \
  --prm-neighbor-radius 1 --dia-neighbor-radius 1 --ms1-neighbor-radius 1
```

Each radius defaults to **0** (off). `1` uses the previous and next compatible
observation, within **5 seconds** of the current frame. PRM matches targets; DIA
matches isolation windows; MS1 skips fragment frames. The combined spectrum only
informs filtering: output retains native points and intensities unless optional
smoothing/centroiding is enabled. Quantitative accuracy remains unvalidated.
See [neighbor options, boundaries, and validation](docs/neighbors.md).

## Compatibility

Type-2 synchro-PASEF, midia-PASEF, and Slice-PASEF examples are supported through
the DIA path, including experimental MS/MS and neighbor filtering. See
[downloaded examples, geometry handling, and validation](docs/acquisition-examples.md).

prm-PASEF type-2 `.d` files support automatic detection and checked target/event
metadata. By default, only MS1 is denoised and PRM fragments are preserved unless
explicitly cropped. **`--denoise-msms` enables experimental PRM fragment filtering
within each isolation event**; validate quantitative results before routine use.
Discovery acquisition gates stay disabled. `--all-frames` is rejected for PRM;
mixed/unknown acquisitions reject all fragment denoising. See
[targeted proteomics support and validation](docs/targeted.md).

dnoise reads compression **type 2** `.d` input and always
writes type 2, byte-layout compatible with the Bruker SDK / `timsdata` DLL.
Validate any output with
`cargo run --release --example validate -- <PATH.d>`.

## Reproducing the paper

The accompanying manuscript evaluates the immutable
[dnoise v0.1.0 release](https://github.com/pgarrett-scripps/dnoise/releases/tag/v0.1.0),
which is also archived on Zenodo under
[10.5281/zenodo.21959649](https://doi.org/10.5281/zenodo.21959649). The raw
benchmark `.d` files are public on PRIDE
([PXD070049](https://www.ebi.ac.uk/pride/archive/projects/PXD070049)); the
manuscript and Supporting Information report the complete dnoise parameters
and downstream search settings.

## Citing this work

If you use dnoise in your research, please cite it. Machine-readable metadata is
in [CITATION.cff](CITATION.cff) (GitHub's "Cite this repository" button reads it),
and each tagged release is archived on Zenodo.

> Garrett, P., Diedrich, J. K., & Yates, J. R. III. dnoise (version 0.1.0) [Software].
> Zenodo. https://doi.org/10.5281/zenodo.21959649

The accompanying paper has been submitted to the *Journal of
the American Society for Mass Spectrometry*. Its preprint and journal citation
will be added here when available.

Planned maintenance, usability improvements, and research extensions are in the
[project roadmap](ROADMAP.md).

## License

Licensed under the [MIT License](LICENSE).
