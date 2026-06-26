#!/usr/bin/env python3
"""DIA whole-frame MS/MS denoising sweep on dia_15min Condition A.

Companion to 40_msms_sweep.py, for diaPASEF. On DIA the MS1 survey is already
saturated by the `--dia-ms1-window` gate (~91% of MS1 peaks removed) and MS2
dominates the data, so the compression headroom is the whole-frame MS2 denoiser.
We sweep its two knobs, `msms_max_internal_gap` x `msms_min_feature_length`, with
`--denoise-msms --dia-ms1-window --dia-window` (the DIA msms arm), scoring
compression (MS1/MS2/bytes) and quantification against the raw runs' 2nd-pass
empirical library.

SEARCH: single pass against results/dia_15min/original/report-lib.parquet, NO
`--reanalyse` (~20 min/arm vs ~2 h with reanalyse). This is also the right design
for a sweep: every arm is searched identically against one fixed library, so the
only thing that varies is the denoised data (no cross-run MBR coupling). It does
mean absolute counts differ from the main DIA benchmark (which uses `--reanalyse`).

Disk-safe + idempotent like 40_*: denoise -> record reduction.json -> single-pass
search -> delete the (large) denoised .d. An arm whose report.pg_matrix.tsv exists
is skipped, so the sweep resumes after an interruption (and reuses a pre-seeded
g8_l3 from the earlier one-off run).
"""

from __future__ import annotations

import json
import os
import shutil
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COND = os.environ.get("SWEEP_COND", "A")

RAW = ROOT / "data" / "dia_15min" / "raw"
SWEEP_DATA = ROOT / "data" / "dia_15min" / "msms_sweep"        # temp denoised arms
SWEEP_RES = ROOT / "results" / "dia_15min" / "msms_sweep"      # kept DIA-NN outputs
DNOISE = (ROOT / ".." / "target" / "release" / "dnoise").resolve()
CONFIG = ROOT / "config" / "dnoise.toml"
DIANN = Path(os.environ.get("DIANN_BIN",
             "/home/patrick-garrett/tools/diann-2.2.0/diann-linux"))
LIB = ROOT / "results" / "dia_15min" / "original" / "report-lib.parquet"
DIANN_FASTA = ROOT / "data" / "fasta" / "hybrid_diann.fasta"
THREADS = str(os.cpu_count())

GAPS = [3, 5, 8]
LENS = [3, 5]
GRID = [(g, l) for g in GAPS for l in LENS]


def peaks(d: Path) -> tuple[int, int]:
    c = sqlite3.connect(d / "analysis.tdf")
    m1 = c.execute("SELECT COALESCE(SUM(NumPeaks),0) FROM Frames WHERE MsMsType=0").fetchone()[0]
    m2 = c.execute("SELECT COALESCE(SUM(NumPeaks),0) FROM Frames WHERE MsMsType!=0").fetchone()[0]
    c.close()
    return int(m1), int(m2)


def bin_bytes(d: Path) -> int:
    return (d / "analysis.tdf_bin").stat().st_size


def run_diann(files: list[Path], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = str(DIANN.parent) + ":" + env.get("LD_LIBRARY_PATH", "")
    args = [str(DIANN)]
    for f in files:
        args += ["--f", str(f)]
    args += ["--lib", str(LIB), "--fasta", str(DIANN_FASTA), "--matrices",
             "--qvalue", "0.01", "--threads", THREADS,
             "--out", str(outdir / "report.parquet")]
    with open(outdir / "diann.log", "w") as log:
        subprocess.run(args, check=True, stdout=log, stderr=subprocess.STDOUT, env=env)


def main() -> int:
    for p in (DNOISE, DIANN, LIB, DIANN_FASTA):
        if not Path(p).exists():
            print(f"missing: {p}", file=sys.stderr)
            return 1
    files = sorted(RAW.glob(f"LFQ_Ultra2_diaPASEF_15min_50ng_Condition_{COND}_REP*.d"))
    if len(files) != 6:
        print(f"expected 6 replicates for condition {COND}, found {len(files)}", file=sys.stderr)
        return 1

    raw_ms1 = raw_ms2 = raw_bytes = 0
    for f in files:
        m1, m2 = peaks(f)
        raw_ms1 += m1
        raw_ms2 += m2
        raw_bytes += bin_bytes(f)
    print(f"raw (6 reps): MS1 {raw_ms1/1e9:.2f}G  MS2 {raw_ms2/1e9:.2f}G  bytes {raw_bytes/1e9:.1f}GB")
    print(f"grid: {len(GRID)} arms (gap {GAPS} x len {LENS})")

    for g, l in GRID:
        arm = f"g{g}_l{l}"
        res = SWEEP_RES / arm
        if (res / "report.pg_matrix.tsv").exists():
            print(f"skip {arm} (done)")
            continue

        t0 = time.time()
        den = SWEEP_DATA / arm
        if den.exists():
            shutil.rmtree(den)
        den.mkdir(parents=True, exist_ok=True)
        dfiles = []
        for f in files:
            out = den / f.name
            subprocess.run(
                [str(DNOISE), str(f), str(out), "--config", str(CONFIG),
                 "--denoise-msms", "--dia-ms1-window", "--dia-window",
                 "--msms-max-internal-gap", str(g), "--msms-min-feature-length", str(l)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            dfiles.append(out)
        km1 = km2 = kb = 0
        for d in dfiles:
            a, b = peaks(d)
            km1 += a
            km2 += b
            kb += bin_bytes(d)
        t_den = time.time() - t0

        res.mkdir(parents=True, exist_ok=True)
        (res / "reduction.json").write_text(json.dumps({
            "arm": arm, "gap": g, "len": l,
            "ms1_kept": km1, "ms1_raw": raw_ms1,
            "ms2_kept": km2, "ms2_raw": raw_ms2,
            "bytes_kept": kb, "bytes_raw": raw_bytes}))

        run_diann(dfiles, res)
        shutil.rmtree(den)
        print(f"{arm}: MS2 -{100*(1-km2/raw_ms2):.1f}%  bytes -{100*(1-kb/raw_bytes):.1f}%  "
              f"denoise {t_den:.0f}s  total {time.time()-t0:.0f}s")

    if SWEEP_DATA.exists() and not any(SWEEP_DATA.iterdir()):
        SWEEP_DATA.rmdir()
    print("dia msms sweep complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
