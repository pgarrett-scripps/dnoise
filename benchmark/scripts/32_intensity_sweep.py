#!/usr/bin/env python3
"""Intensity-threshold control arms for the default-parameter sweep.

Same one-condition replicate set as 30_param_sweep.py, but instead of the streak
filter each arm is a pure per-point MS1 intensity threshold (keep iff intensity
>= T; config/dnoise.intensity.15min.toml collapses dnoise to this). T values are
calibrated (scripts/_calib_intensity.py) to remove ~20/40/60/80% of MS1 points,
so the control sits at matched data reduction to the streak arms and 31 can plot
streak vs threshold on a common reduction axis.

Same disk-safe denoise -> search -> delete loop and per-arm idempotency as 30.

Outputs (kept):
  results/dda_15min/sweep/int<pct>/   Sage outputs + reduction.json (kind=intensity)
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
SWEEP_DATA = ROOT / "data" / "dda_15min" / "sweep"
SWEEP_RES = ROOT / "results" / "dda_15min" / "sweep"
DNOISE = (ROOT / ".." / "target" / "release" / "dnoise").resolve()
ICONFIG = ROOT / "config" / "dnoise.intensity.15min.toml"
FASTA = ROOT / "data" / "fasta" / "hybrid.fasta"
SAGE_CFG = ROOT / "config" / "sage.json"
SAGE = (ROOT / "tools" / "sage-0.15.0-beta.1"
        / "sage-v0.15.0-beta.1-x86_64-unknown-linux-gnu" / "sage")

# (label, T, nominal target%) from scripts/_calib_intensity.py.
ARMS = [("int20", 41, 20), ("int40", 54, 40), ("int60", 72, 60), ("int80", 104, 80)]


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
    for p in (DNOISE, SAGE, SAGE_CFG, FASTA, ICONFIG):
        if not Path(p).exists():
            print(f"missing: {p}", file=sys.stderr)
            return 1

    files = sorted(RAW.glob(f"LFQ_Ultra2_PASEF_15min_50ng_Condition_{COND}_REP*.d"))
    if len(files) != 6:
        print(f"expected 6 replicates for condition {COND}, found {len(files)}",
              file=sys.stderr)
        return 1
    raw_ms1 = sum(ms1_peaks(f) for f in files)
    print(f"condition {COND}: {len(files)} reps, raw MS1 {raw_ms1/1e6:.1f}M")
    print(f"intensity arms: {[(a, t) for a, t, _ in ARMS]}")

    for arm, T, tgt in ARMS:
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
                [str(DNOISE), str(f), str(out), "--config", str(ICONFIG),
                 "--ms1-polygon",  # match the streak arms' DDA polygon pre-gate
                 "--min-window-intensity", str(T)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            den_files.append(out)
        kept_ms1 = sum(ms1_peaks(d) for d in den_files)
        res.mkdir(parents=True, exist_ok=True)
        (res / "reduction.json").write_text(json.dumps(
            {"arm": arm, "kind": "intensity", "T": T, "target_pct": tgt,
             "gap": None, "len": None, "ms1_kept": kept_ms1, "ms1_raw": raw_ms1}))
        run_sage(den_files, res)
        shutil.rmtree(den_dir)
        red = 100 * (1 - kept_ms1 / raw_ms1)
        print(f"{arm} (T={T}): MS1 -{red:.1f}%  total {time.time()-t0:.0f}s")

    print("intensity sweep complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
