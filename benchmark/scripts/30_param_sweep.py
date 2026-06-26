#!/usr/bin/env python3
"""Default-parameter sweep on one dda_15min sample with replicates.

Picks one HYE condition (default A, 6 replicates) and sweeps the two main
streak-filter knobs - max_internal_gap x min_feature_length - re-denoising the
replicates from raw for each setting, then running Sage LFQ on the replicate
set. The companion 31_sweep_analyze.py scores each arm by quantified peptides /
proteins / precision vs. MS1 data reduction so we can pick a defensible default.

Disk-safe: each arm is denoised, searched, and then its (large) denoised .d
folders are deleted, keeping only the small Sage TSV outputs. Peak extra disk is
~6 GB (one arm of 6 x ~0.9 GB) regardless of grid size.

Idempotent: an arm whose results/.../sweep/<arm>/lfq.tsv already exists is
skipped, so the sweep resumes after an interruption.

Outputs (kept):
  results/dda_15min/sweep/original/          Sage outputs for the raw arm
  results/dda_15min/sweep/g<gap>_l<len>/     Sage outputs + reduction.json per arm
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
BATCH = os.environ.get("SAGE_BATCH", "3")

RAW = ROOT / "data" / "dda_15min" / "raw"
SWEEP_DATA = ROOT / "data" / "dda_15min" / "sweep"          # temp denoised arms
SWEEP_RES = ROOT / "results" / "dda_15min" / "sweep"        # kept Sage outputs
DNOISE = (ROOT / ".." / "target" / "release" / "dnoise").resolve()
CONFIG = ROOT / "config" / "dnoise.toml"
FASTA = ROOT / "data" / "fasta" / "hybrid.fasta"
SAGE_CFG = ROOT / "config" / "sage.json"
SAGE = (ROOT / "tools" / "sage-0.15.0-beta.1"
        / "sage-v0.15.0-beta.1-x86_64-unknown-linux-gnu" / "sage")

# 3 x 5 rectangular grid over the two main knobs (halo + iterations=2 from config),
# plus two extra long-streak arms at gap 3 (Figure S11) to probe larger
# min_feature_length than the base grid.
GAPS = [1, 2, 3]
LENS = [3, 4, 5, 6, 7]
EXTRA = [(3, 9), (3, 11)]
GRID = [(g, l) for g in GAPS for l in LENS] + EXTRA


def ms1_peaks(d: Path) -> int:
    c = sqlite3.connect(d / "analysis.tdf")
    v = c.execute(
        "SELECT COALESCE(SUM(NumPeaks),0) FROM Frames WHERE MsMsType=0"
    ).fetchone()[0]
    c.close()
    return int(v)


def run_sage(arm_files: list[Path], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [str(SAGE), str(SAGE_CFG), "-f", str(FASTA), "-o", str(outdir),
           "--batch-size", BATCH] + [str(f) for f in arm_files]
    with open(outdir / "sage_run.log", "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)


def main() -> int:
    for p in (DNOISE, SAGE, SAGE_CFG, FASTA):
        if not Path(p).exists():
            print(f"missing: {p}", file=sys.stderr)
            return 1

    files = sorted(RAW.glob(f"LFQ_Ultra2_PASEF_15min_50ng_Condition_{COND}_REP*.d"))
    if len(files) != 6:
        print(f"expected 6 replicates for condition {COND}, found {len(files)}",
              file=sys.stderr)
        return 1
    print(f"condition {COND}: {len(files)} replicates")
    print(f"grid: {len(GRID)} arms (gap {GAPS} x len {LENS}) + original")

    raw_ms1 = sum(ms1_peaks(f) for f in files)
    print(f"raw MS1 peaks (6 reps): {raw_ms1/1e6:.1f}M")

    # Original (raw) baseline.
    orig = SWEEP_RES / "original"
    if (orig / "lfq.tsv").exists():
        print("skip original (done)")
    else:
        t = time.time()
        run_sage(files, orig)
        (orig / "reduction.json").write_text(json.dumps(
            {"arm": "original", "gap": None, "len": None,
             "ms1_kept": raw_ms1, "ms1_raw": raw_ms1}))
        print(f"original searched in {time.time()-t:.0f}s")

    for g, l in GRID:
        arm = f"g{g}_l{l}"
        res = SWEEP_RES / arm
        if (res / "lfq.tsv").exists():
            print(f"skip {arm} (done)")
            continue

        t0 = time.time()
        den_dir = SWEEP_DATA / arm
        if den_dir.exists():
            shutil.rmtree(den_dir)
        den_dir.mkdir(parents=True, exist_ok=True)
        den_files = []
        for f in files:
            out = den_dir / f.name
            subprocess.run(
                [str(DNOISE), str(f), str(out), "--config", str(CONFIG),
                 "--ms1-polygon",  # DDA benchmark default gate (see 04_denoise.sh)
                 "--max-internal-gap", str(g), "--min-feature-length", str(l)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            den_files.append(out)
        kept_ms1 = sum(ms1_peaks(d) for d in den_files)
        t_den = time.time() - t0

        res.mkdir(parents=True, exist_ok=True)
        (res / "reduction.json").write_text(json.dumps(
            {"arm": arm, "gap": g, "len": l,
             "ms1_kept": kept_ms1, "ms1_raw": raw_ms1}))

        run_sage(den_files, res)
        shutil.rmtree(den_dir)  # reclaim disk before the next arm
        red = 100 * (1 - kept_ms1 / raw_ms1)
        print(f"{arm}: MS1 -{red:.1f}%  denoise {t_den:.0f}s  "
              f"total {time.time()-t0:.0f}s")

    if SWEEP_DATA.exists() and not any(SWEEP_DATA.iterdir()):
        SWEEP_DATA.rmdir()
    print("sweep complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
