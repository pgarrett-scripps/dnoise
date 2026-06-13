# dnoise DDA benchmark

Measures how `dnoise` affects downstream **DDA** results — number of
identifications and **label-free quantification** — using a dataset with known
ground-truth ratios.

## Dataset

[PXD070049](https://www.ebi.ac.uk/pride/archive/projects/PXD070049) — *"LFQ
Benchmark Dataset - Generation Beta"*. HYE 3-species hybrid (ratios from the
dataset SDRF):

| Mixture (`Condition`) | Human | Yeast (*S. cerevisiae*) | *E. coli* |
|---|---|---|---|
| A | 65% | 30% | 5% |
| B | 65% | 15% | 20% |
| C | 65% | 3% | 32% |

→ expected log2: A/B = Human **0** / Yeast **+1** / E. coli **−2**; the A/C and
B/C pairs widen the dynamic range (yeast up to **+3.3**, ecoli down to **−2.7**).

We use the **timsTOF Ultra 2 ddaPASEF** runs (`LFQ_Ultra2_PASEF_*`; "PASEF" =
DDA, as opposed to `diaPASEF`), conditions A/B/C in six replicates → **18 runs
per gradient**. The deposit offers each mode at two gradients, **5 min and 15
min** (nothing longer); both are dnoise-readable `TimsCompressionType: 2` `.d`
(~0.75 GB each at 5 min, ~1.6 GB at 15 min). The search DB is the authors' own
`uniprotkb_proteome_HYE_UniversalContaminants.fasta`, re-tagged by species
([scripts/02_fetch_fasta.py](scripts/02_fetch_fasta.py)).

### Datasets (the `DATASET` namespace)

Everything is namespaced by a `DATASET` so multiple acquisition sets coexist
under `data/<dataset>/` and `results/<dataset>/` without clobbering. Known
datasets (see [scripts/01_list_files.py](scripts/01_list_files.py)):

| `DATASET` | Acquisition group | Engine |
|---|---|---|
| `dda_5min` *(default)* | `LFQ_Ultra2_PASEF_5min_50ng` | Sage |
| `dda_15min` | `LFQ_Ultra2_PASEF_15min_50ng` | Sage |
| `dia_5min` | `LFQ_Ultra2_diaPASEF_5min_50ng` | DIA-NN *(stub)* |
| `dia_15min` | `LFQ_Ultra2_diaPASEF_15min_50ng` | DIA-NN *(stub)* |

## Design

Three arms searched with an **identical** Sage config + FASTA, all parameterized
by [config/dnoise.toml](config/dnoise.toml) (pinned to dnoise defaults; edit it
for a sweep, or point `DNOISE_CONFIG=...` at an alternate file):

- **original** (`results/<dataset>/original`) — the raw `.d`
- **denoised** / MS1 (`.../denoised`) — `dnoise` on **MS1 frames only** (default)
- **msms** / MS1+MS/MS (`.../msms`) — adds per-precursor MS/MS denoising
  (`--denoise-msms`); changes fragment spectra, so IDs move on this arm

The MS1-only arm leaves MS/MS untouched, so its IDs are **identical** to original
— it isolates the quant effect. The msms arm trades a few % of IDs for cleaner
fragment spectra. Engine: **Sage** (reads Bruker `.d` directly, LFQ with
cross-run matching). FASTA headers are species-tagged (`HUMAN_`, `YEAST_`,
`ECOLI_`, `CONT_`) so proteins split cleanly by species.

Each denoised `.d` carries a `.dnoise_stamp` (sha256 of the config + dnoise
version + arm). [scripts/04_denoise.sh](scripts/04_denoise.sh) rebuilds an output
whenever its stamp doesn't match the current config/binary, instead of blindly
skipping anything that already exists — so a config or algorithm change can't
leave a stale `.d` behind. Force a full rebuild with `DENOISE_FORCE=1`.

Note: timsTOF DDA files contain empty MS/MS frames (0 peaks) that timsrust cannot
decode; dnoise handles these directly (emitting Bruker's canonical empty record).

## Run

```bash
cd benchmark
just fasta                      # build hybrid.fasta (shared across datasets)

# default dataset (dda_5min):
just list                       # write config/files.dda_5min.tsv (cheap, no download)
just download                   # ~13 GB zipped (18 files)
just denoise                    # dnoise each raw .d (MS1 + MS/MS arms)
just search                     # Sage x3 (SAGE_BATCH=2 by default; lower if RAM-limited)
just analyze && just datasize   # metrics, plots, data-reduction figure
# or, end to end:  just all

# a second gradient — same recipes, different namespace:
just dataset=dda_15min all      # ~29 GB zipped (18 files)
just compare                    # overlay dda_5min vs dda_15min -> results/compare/
```

Requires `uv`, `sage` on PATH, and the built `dnoise` binary
(`cargo build --release` in the repo root). `data/` and `results/` are gitignored.

## Outputs (`results/<dataset>/analysis/`)

- `summary.csv` — the three arms, side by side:
  - IDs @ 1% FDR: PSMs, peptides, protein groups (+ per-species peptide counts)
  - LFQ accuracy: per-species median log2(A/B), **bias** vs expected, MAD spread
  - LFQ precision: median protein CV within condition
  - completeness: # proteins quantified (≥2/3 reps in both conditions)
- `accuracy.csv` — observed vs expected median log2 for A/B, A/C, B/C × species
- `data_reduction.{png,csv}` — frame-binary size + peak counts (MS1/MS2) per arm
- `lfq_accuracy.png`, `lfq_ratio_violins.png`, `lfq_cv.png`, `id_counts.png`

Cross-dataset comparison (`results/compare/`, via `just compare`):
- `gradient_compare.{png,csv}` — data reduction, IDs, quantified proteins, and
  LFQ precision for each arm, side by side across the compared datasets.

## LFQ requires a Sage with ion-mobility quant

timsTOF MS1 LFQ needs **Sage 0.15.0-beta.1** (where ion-mobility LFQ was added).
[scripts/05_search.sh](scripts/05_search.sh) auto-fetches that Linux x86_64
prebuilt into `tools/` and uses it; override with `SAGE_BIN=/path/to/sage` (e.g.
on other platforms). It runs LFQ **directly on the `.d`** — no mzML conversion or
external quant tool. (`06_analyze.py` reads Sage's `lfq.tsv` straight from each
arm's `results/` dir.)

Version notes:
- **0.14.6 / stable (≤ 0.14.7)** reads the `.d` and searches fine, but its
  MS1-feature scoring collapses the mobility dimension and reports "0 target MS1
  peaks at 5% FDR" → empty `lfq.tsv` → `n_quantified = 0`.
- **0.15.0-beta.2** *prebuilt* panics resolving the `.d` path (unfinished
  `with_path` in timsrust 0.4.2) — do not use it.

ID counts are unaffected by Sage version. (For DDA, denoise leaves MS/MS
untouched, so IDs are identical across arms regardless — the quant columns are
the only place denoise can move.)

## Sanity gate (read before trusting the comparison)

On the **original** arm the three species' median log2(A/B) must land near
**0 / +1 / −2** as separated clouds. If not, the search/FASTA/quant setup is wrong
and the denoise comparison is meaningless — fix that first.
