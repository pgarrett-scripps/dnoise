# Repeatable acquisition compatibility checks

Run the downloaded example datasets through MS1-only, MS1 + MS/MS, and neighbor
recipes, then compare every output frame with the original through the Bruker SDK:

```bash
python3 -m pip install -r examples/requirements-validation.txt
cargo build --release --bin dnoise --locked
python3 examples/validate_acquisitions.py \
  --binary target/release/dnoise \
  --sdk /path/to/libtimsdata.so \
  --output-root test-data/acquisitions/validation-run-1
```

Requires Python 3.10+, NumPy, zstandard, and an independently installed Bruker SDK library
(`libtimsdata.so` or `timsdata.dll`). The SDK is not bundled or downloaded. The
output root must be new; the command never overwrites existing runs. Run from the
repository root, or supply absolute paths. To check just one method, add
`--method synchro-PASEF` (repeat `--method` to select several).

The default manifest is `test-data/acquisition-examples.json`; native file sizes
and SHA-256 checksums are verified before processing. Paths in a custom
`--manifest FILE` are relative to that manifest. Use the same schema to register
other type-2 DDA, DIA, or PRM examples; mixed/unsupported MS/MS types are rejected.
No datasets are downloaded by this command. A recorded checksum mismatch stops
processing rather than silently using different data.

## Recipes

The raw input is the reference and is never rewritten. The three generated
outputs use CLI defaults without smoothing, centroiding or cropping:

| Recipe | Processing |
| --- | --- |
| `ms1` | Default MS1 denoising; all native fragment points preserved |
| `msms` | MS1 and MS/MS denoising, with all temporal radii zero |
| `neighbors` | MS1 and MS/MS denoising, MS1 radius 1, DIA or PRM radius 1, default 5-second RT bound |

DDA has no MS/MS temporal-radius option, so its neighbor recipe adds only MS1
neighbors. Four threads and a 32-frame batch are explicit runner defaults to
limit queued output memory; adjust `--threads` and `--frame-batch-size` if needed.
Full input/output validation remains enabled. Single-calibration inputs use the
same filters as ordinary CLI runs. Multi-calibration inputs requesting physical
gates stop with the calibration error; the runner does not silently disable them.

## Evidence and failure behavior

`validation.json` records the status, binary/SDK/manifest SHA-256, input and output
native-file checksums, exact commands, complete CLI reports, Python/dependency versions, verifier checksum, and SDK-derived
point/intensity totals. Each recipe has a separate log and report beside its
output. The final status is `passed` only after:

- All frames decode through the SDK, preserving native retained coordinates and
  intensities; default MS1-only fragment spectra are identical.
- SQLite schema and non-writer metadata, including isolation trajectories,
  calibrations, frame identities and timing, are unchanged.
- Reported MS1/MS/MS intensity sums and total point counts agree exactly with a
  separate Python decoder of the stored integer payload. SDK and native point
  coordinates also agree, and retained intensities match input in both readers.
- Input native files, executable and SDK still match their original checksums.

Errors return a nonzero exit code. Once processing starts, failures preserve the
logs and a failed summary for diagnosis; there is no automatic overwrite/resume.
Preflight failures create no output root. Neighbor usage counts describe evidence
used by the filter, not whether a particular neighbor caused a point to survive.

These are file compatibility and reporting checks. They do not evaluate peptide
identification, fragment-ratio bias, quantitative precision or biological accuracy.
The report's intensity basis is `stored_integer`. SDK intensity sums are retained
separately as `sdk_totals`: the inspected synchro SDK reads differed slightly
from the stored integers, consistent with accumulation-time normalization and
integer rounding in the diagnostic frames. The runner does not assume SDK sums
and stored-integer sums are identical or relax the exact native-subset checks.

The earlier one-off results remain in `test-data/acquisition-verification.json`.
The old neighbor recipe there used DIA radius 1 with MS1 radius 0, so its MS1
retention is not directly comparable with this runner's neighbor recipe.

Verifier unit tests (no SDK required):

```bash
python3 -m unittest discover -s examples -p 'test_validate_acquisitions.py'
```

## Verified development run (2026-09-10)

All three downloaded examples passed all three recipes, with 42,444 SDK frame
reads and matching independent native-coordinate/retained-intensity checks.
The new report sums matched the separately decoded stored integers exactly.
Input native files were unchanged. Compact evidence, executable/SDK/verifier
hashes, effective configurations and per-recipe QC are recorded in
[`test-data/hardening-verification.json`](../test-data/hardening-verification.json).

The local workspace passed 170 Rust tests (two ignored), five verifier tests,
Clippy with warnings denied, formatting, the Rust 1.85 core check and 75
library-only tests. CI now includes the verifier's SDK-free regression tests;
full real-file checks require the separately supplied SDK and datasets.
