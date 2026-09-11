# Parallel validation

Full-frame validation now decodes records using **up to four tasks in the current
Rayon pool**. Structural checks still verify SQLite integrity, offsets, headers,
and overlaps; every payload is still decoded and checked against metadata.
This applies to normal processing, dry runs, and `dnoise validate`.
Full validation is enabled by default; normal runs can opt out with
[`--skip-validation`](reference.md#optional-full-validation). The measurements
below compare serial and parallel **full** validation, not skipped validation.

One reader loads records in physical file order. Each batch holds at most 256
records and normally at most 16 MiB of compressed data; a larger supported
record runs alone. Each decoder releases its frame before starting another.
The existing 256 MiB per-record/decompressed-size limits remain in force.
The four-task cap limits simultaneous decode allocations on large machines;
a caller's smaller pool is respected. Cancellation is checked between reads
and decoded frames. In-flight frame decoding may finish before return.

## Local comparison — September 2026

Compared release binaries immediately before and after the validation change,
using the same DDA and DIA acquisitions as the
[implementation checks](implementation-validation-2026-09-09.md). Both binaries
perform full input/output validation. Linux, 20 logical CPUs, default Rayon pool,
and a warm filesystem cache; benchmark commands ran sequentially without builds.

| Check | Serial validation | Rayon validation | Difference |
|---|---:|---:|---:|
| DDA validation alone, median of 3 | 4.03 s | 1.58 s | 2.55× faster |
| DIA validation alone, median of 3 | 12.79 s | 5.71 s | 2.24× faster |
| Complete DDA run, one pair | 13.31 s | 10.21 s | 23% less elapsed time |
| Complete DIA run, one pair | 42.73 s | 29.49 s | 31% less elapsed time |

Standalone checks alternated serial/parallel order across repetitions. Complete
runs were serial followed by parallel. These are local measurements on two
acquisitions, not guarantees for other CPUs, disks, data, or cold caches. The
paper release performs fewer checks and remains a separate baseline.

Parallel decoding trades memory for speed. Median standalone peak RSS rose from
**24 to 120 MiB** for DDA and **39 to 203 MiB** for DIA. Whole-run peak RSS changed
from 1,350,888 to 1,389,456 KiB for DDA and 4,459,220 to 4,517,864 KiB for DIA.
Compressed-buffer limits are not a total-process memory cap.

Both complete output binaries retained the previously verified paper hashes:

```text
DDA 74474dbf180f47c74a788210687f339d658c846c00323d88436a34ed2efd85e6
DIA 75a49f6f6f26f2ae5ab81ba60ab66442ae7a7ba50cb8c5510633e8942a84a78f
```

The regression suite validates a 600-frame fixture with one- and four-thread
pools, including decoded-count mismatches and corrupt payloads in later batches.
Existing overlap, empty-frame, cancellation, and output-preservation tests remain.
All 139 workspace tests pass, along with Clippy, formatting, the Rust 1.85 core
check, and 75 library-only tests. The real DDA cancellation probe returned in
52.697 ms during validation and 53.922 ms during processing; both preserved the
existing output. These are individual probes, not maximum-latency guarantees.
Benchmark commands, raw timings, executable hashes, source snapshot, and outputs
are local artifacts in `/tmp/dnoise-rayon-zyxge35l`; they are not archived fixtures.
