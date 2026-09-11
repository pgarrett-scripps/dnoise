# Maintenance review: 2026-09-09

Historical review of the code **before implementation**. The subsequent changes
and checks are recorded in [implementation validation](implementation-validation-2026-09-09.md).

**Recommendation: ship a file-safety patch first, then strengthen verification
and reporting before expanding the algorithm.** The code has a useful modular
core and passing tests, but the file lifecycle has confirmed destructive edge
cases that normal tests do not cover.

Reviewed local commit `71b87bacf32eb4581878455361c8cb6986aec74f` and compared with
the local paper tag `v0.1.0` (`fa750c92e2daed4ad119c7f142c0510be4d8c385`).
There are no changes in `src`, `dnoise-gui/src`, `Cargo.lock`, or the workflows
between those two revisions. The adjacent paper repository was read for the
manuscript's limitations and cross-platform results.

The [roadmap](../ROADMAP.md) translates these findings into proposed releases.
No production code was changed as part of this review.

## Verified state

| Check | Result |
|---|---|
| `cargo test --workspace --locked --offline` | 116 passed, including one doctest; two ignored |
| `cargo fmt --all --check` | Passed |
| `cargo clippy --workspace --all-targets --locked --offline -- -D warnings` | Passed |
| `cargo +1.85.0 check --all-targets --locked --offline` | Passed for the root library/CLI package, not the GUI |
| Library-only check and `cargo test --no-default-features --lib --locked --offline` | Passed; 70 library tests |
| `cargo build --bin dnoise --locked --offline` | Passed; executable used in disposable-file probes |

The ignored integration test requires a real `.d` fixture. The other ignored
item is a `RunContext` documentation example. No real acquisition was modified,
and no full scientific benchmark, native GUI interaction, independent SDK
validation, or current RustSec vulnerability scan was run during this review.

GitHub read-only checks returned no open issues and successful recent CI runs.
The most recent returned [CI run](https://github.com/pgarrett-scripps/dnoise/actions/runs/31915768135)
was for `2cc968f`, not the reviewed local HEAD. The
[v0.1.0 release](https://github.com/pgarrett-scripps/dnoise/releases/tag/v0.1.0)
has Linux x86_64, macOS ARM64, and Windows x86_64 archives with API-reported
SHA-256 digests. Their contents and correspondence to the archive were not
independently downloaded or checked here.

## 1. File lifecycle: confirmed first priority

**Reproduced:** the library deletes an existing destination under `force` before
opening the input reader. The existence of the two input filenames is the only
validation preceding deletion. See [writer.rs](../src/writer.rs), `run`, around
lines 222–244.

Two probes used freshly created disposable directories under `/tmp` and the
current built CLI:

| Probe | Observed result |
|---|---|
| Identical input/output paths, two dummy input files, `--force` | Exit 1; both original input files and a marker were gone |
| Invalid input contents, separate pre-existing output, `--force` | Exit 1; the old output marker was gone; input remained |

Related risks established by inspection, not crash/race reproduction:

1. `copy_dir_except` creates the destination before enumerating the source.
   Output nested inside input can therefore be copied recursively into itself.
   There is no shared overlap/alias guard in the library.
2. In-place mode uses predictable `.dnoise-tmp` and `.dnoise-old` siblings and
   removes existing siblings. A failed restore is ignored while the error text
   says the original was restored. The two-rename sequence is not an atomic
   crash-safe directory exchange. See [main.rs](../src/main.rs), lines 718–809.
3. GUI cancellation cleanup tests the shared cancellation flag after any error
   and removes the output path. It does not establish that the failed operation
   created that output. See [worker.rs](../dnoise-gui/src/worker.rs), lines 165–182.

Put path validation and ownership-aware staging in the library so all callers
receive the same guarantees. Add failure-injection tests before releasing it.

## 2. Verification and release gaps

The normal suite has good focused coverage of filtering, gate geometry,
centroiding, configuration, and codec round trips. Full writer/streaming parity
is ignored in ordinary CI. Its current `Stages::default()` also disables the
halo and acquisition gates, so it does not cover CLI-default wiring even when
enabled. See [streaming.rs](../tests/streaming.rs).

The [validation example](../examples/validate.rs) rereads using timsrust and
checks peak counts. This is useful, but not an independent Bruker SDK oracle.
Make both the automated fixture check and independent release validation explicit.

The [MSRV job](../.github/workflows/ci.yml) installs 1.85.0 but invokes plain
`cargo`; the repository pins 1.97.1. The upstream action sets the rustup default,
while the repository toolchain file takes precedence over that default.
Consequently, the job does not reliably test its named compiler. Use explicit
`cargo +1.85.0` and record `rustc --version`. The direct 1.85.0 check passed in
this review. See the [rustup precedence documentation](https://rust-lang.github.io/rustup/overrides.html)
and [action implementation](https://github.com/dtolnay/rust-toolchain/blob/master/action.yml).

Release jobs build across three platforms but do not themselves run the test
suite or depend on a testing job. Add test gates and use `--locked` to require the
resolved dependencies to match the lockfile. See [Cargo's build documentation](https://doc.rust-lang.org/cargo/commands/cargo-build.html).

## 3. Public API, configuration, and reports

Temporary Rust probes linked against the built library confirmed that
`encode_frame_type2(0, &[])` panics and decoding a canonical empty record emitted
by `encode_empty_frame_type2(3)` returns an error. The production writer has
separate empty-frame handling, so these probes do not establish a failure on
ordinary empty acquisition frames. They establish public codec edge cases.
See [codec.rs](../src/codec.rs).

Inspection also found unbounded `zstd::decode_all`, unchecked scan indexing in
the encoder, and incomplete decoded offset/header validation. Introduce checked
APIs and bounded decoding with focused malformed-input tests.

CLI JSON reports include vertical filter parameters but mostly booleans for
optional stages. Their halo thresholds, gate padding, and other active stage
parameters cannot be reconstructed. Some fields use raw CLI values rather than
resolved configuration. GUI reports contain statistics without effective
configuration or software version. See `build_report` in
[main.rs](../src/main.rs) and `write_report` in
[worker.rs](../dnoise-gui/src/worker.rs).

The GUI shares the TOML schema but applies only a subset of it. Its worker fixes
MS/MS denoising, averaging, smoothing, and centroiding off, and uses default gate
padding. Invalid numeric crop text becomes a warning and an omitted bound.
Resolve settings in a shared layer, reject invalid supplied bounds, and explain
unsupported settings on load. See [settings.rs](../dnoise-gui/src/settings.rs).

## 4. Performance and scientific scope

The writer parallelizes 2,048 frames per chunk and checks cancellation only at
chunk boundaries. DDA MS/MS processing first accumulates raw points across the
run in `build_msms_keep`. Bounded chunk output therefore does not imply bounded
memory for the whole pipeline. Measure both paths before choosing batch sizes
or replacing data structures. See [writer.rs](../src/writer.rs) and
[msms.rs](../src/msms.rs).

The manuscript already goes beyond the main 72-run benchmark with a search-free
cross-platform precursor-area check. Its local `paper/si/platforms_agreement.typ`
reports only 55.7–60.9% seeded precursor retention for the two 250 pg Ultra runs,
with 8.1–15.5% retention in the low-abundance quartile. Those are measurements
from the adjacent paper artifacts, not new validation performed here.

The manuscript's limitations section explicitly distinguishes that search-free
check from replicated identification/quantification validation. The highest
value research extension is therefore low-input fidelity with held-out,
searched data, rather than simply more compression or another platform table.

Multi-segment calibration currently produces a warning for acquisition gates,
which still use run-level conversion. Physical-unit crops also use a run-level
converter. Address that limitation before relying on more adaptive physical-unit
thresholds across instruments.

## 5. Maintenance decisions

Keep the paper release immutable and label changes to its algorithm/defaults
explicitly. Add a single reproduction manifest and a current compatibility
matrix. Update the README's submission status now and add the final journal
citation after acceptance.

Do dependency upgrades in reviewed groups after the safety patch and fixture
coverage. No claim is made here that a listed dependency is vulnerable or that
a particular newer version is required. A current lockfile audit remains work
to do; [RustSec](https://rustsec.org/) provides the advisory database and tooling.

Preserve the small Rust core and existing in-process GUI/library architecture.
Shared configuration, reporting, validation, and file lifecycle code offer
more immediate value than a rewrite, cloud service, or additional algorithms.
