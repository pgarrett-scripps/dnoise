# Implementation validation: 2026-09-09

Local, unreleased **0.2.0-dev** implementation of roadmap steps 1–3 and the batch
manifest extension from step 5. The paper’s `v0.1.0` tag is unchanged. No new
release or remote CI run was created as part of this work.

## Automated checks

| Check | Result |
|---|---|
| `cargo test --workspace --locked --offline` | 138 passed; two ignored |
| `cargo fmt --all --check` | Passed |
| `cargo clippy --workspace --all-targets --locked --offline -- -D warnings` | Passed |
| `cargo +1.85.0 check --all-targets --locked --offline` | Passed for CLI/library |
| `cargo test --no-default-features --lib --locked --offline` | 75 passed |
| `cargo audit` against the updated lockfile | Zero known vulnerabilities; one unmaintained warning |

Coverage includes owned temporary output, path aliases and overlaps, cancellation,
failed installation/restoration, ancillary copy failures, malformed codec input,
synthetic DDA/DIA round trips, streaming agreement, physically unordered frame
records, provenance history, replay settings, CLI overrides, and batch collisions.
Eleven GUI tests include shared settings, preview calculations, and headless UI
rendering. The ignored items are the opt-in real-file streaming parity test and a
`RunContext` documentation example; synthetic streaming parity runs normally.

The remaining advisory is `RUSTSEC-2024-0436` for transitive `paste` 1.0.15
(unmaintained). It is not suppressed. The GUI framework update raises its minimum
Rust version to 1.95; the root package retains 1.85.

## Real acquisitions and independent reads

Built the archived paper tag `v0.1.0`
(`fa750c92e2daed4ad119c7f142c0510be4d8c385`) in a separate temporary tree and ran
both binaries with their default MS1 settings. Original acquisitions were read
without modification. These are compatibility checks, not new identification or
quantification benchmarks.

| Acquisition | Frames | Input points | Retained points | Output binary bytes |
|---|---:|---:|---:|---:|
| DDA: Ultra2 PASEF, 5 min, 50 ng, condition B, replicate 1 | 8,638 | 301,495,029 | 118,461,309 | 348,012,138 |
| DIA: Ultra2 diaPASEF, 5 min, 50 ng, condition A, replicate 6 | 8,635 | 1,321,633,266 | 763,824,272 | 1,989,868,817 |

For each acquisition, the paper and development `analysis.tdf_bin` files had the
same SHA-256. The entire folder is intentionally different because it now has
sidecars; this assertion is about the frame binary.

```text
DDA 74474dbf180f47c74a788210687f339d658c846c00323d88436a34ed2efd85e6
DIA 75a49f6f6f26f2ae5ab81ba60ab66442ae7a7ba50cb8c5510633e8942a84a78f
```

Bruker's `libtimsdata.so` independently opened both development folders, including
the sidecars, and read every frame with `tims_read_scans_v2`. Decoded peak counts
matched `Frames.NumPeaks` for all 17,273 frames. Repeat this with an available SDK:

```sh
python examples/validate_sdk.py --sdk /path/to/libtimsdata.so output.d
```

Development logs and disposable outputs for these checks are under
`/tmp/dnoise-validation-4qjydp_b`; they are local evidence, not archived fixtures.
The final DDA binary was rechecked after the cancellation, copying, and
thread/report changes and retained the same hash. The DIA check preceded those
operational changes; the filtering algorithm did not change afterward.

## Performance and release limits

The timings below precede parallel validation. See the subsequent
[Rayon comparison](validation-performance.md) for current measurements.

A DDA batch-size probe produced identical binaries with sizes 2,048 and 64.
Peak RSS fell from about 1.29 GiB to 0.99 GiB, but elapsed time rose from 9.78 to
30.72 seconds. These probes overlapped build activity and preceded enabling full
input decoding during validation: they are not controlled release benchmarks.
The default remains 2,048; smaller batches are an explicit memory tradeoff.

Sequential final DDA runs, without a concurrent build, took **13.22 seconds**
for development versus **7.95 seconds** for the paper binary. Peak RSS was
1,398,976 versus 1,417,164 KiB. This is a single local comparison with a warm
filesystem cache, not a replicated benchmark or a speed improvement.

The later DIA run, with full input/output validation enabled, took 42.9 seconds
versus 21.2 seconds for the paper binary. This was one local run, not a general
performance estimate. Validation adds substantial work and needs a controlled
benchmark before release. Cancellation is checked between validation/processing
frames and during DDA preparation; in-flight frame work and file I/O may finish.
A release-build DDA probe set the cancellation flag 100 ms into input validation
and, separately, 100 ms after processing began. Cancellation-to-return latency
was **30.142 ms** and **45.525 ms**, respectively, including cleanup. Both runs
preserved a pre-existing output sentinel. The probe source and log are alongside
the temporary validation artifacts; timings do not cover every storage device,
acquisition, or MS/MS preparation workload.

Linux builds and local checks pass. Native dialogs, interactive desktop behavior,
and macOS/Windows builds still need release smoke checks; configured CI has not
run remotely for these changes. The preview shows binned intensity changes and
does not establish scientific fidelity. Calibration-segment handling and low-input
method development remain the deferred research track.
