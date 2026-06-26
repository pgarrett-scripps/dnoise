#!/usr/bin/env python3
"""MS1 selection-polygon gate: with vs. without, on one dda_15min sample.

The ddaPASEF benchmark applies `--ms1-polygon` (drop MS1 points outside the
instrument's precursor-selection polygon) on top of the MS1 streak+halo filter.
This isolates what that gate alone adds: extra MS1 / byte compression, and
whether it costs any identifications or quantification.

Three arms on one HYE condition (default A, 6 replicates), MS1-only denoising
(no `--denoise-msms`, so MS/MS frames are byte-identical between the two denoised
arms and any ID difference is purely the polygon's MS1 effect):

  original : raw data, no denoising (reused from the MS1 sweep if present).
  nopoly   : dnoise --config dnoise.toml                 (streak + halo, no gate)
  poly     : dnoise --config dnoise.toml --ms1-polygon   (+ selection-polygon gate)

Disk-safe and idempotent like 30/40_*: each arm is denoised (recording MS1-peak
and tdf_bin-byte reduction), searched, then its denoised .d folders deleted. An
arm whose results/.../polygon_compare/<arm>/lfq.tsv exists is skipped.

Outputs (kept):
  results/dda_15min/polygon_compare/{original,nopoly,poly}/  Sage + reduction.json
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
CMP_DATA = ROOT / "data" / "dda_15min" / "polygon_compare"      # temp denoised arms
CMP_RES = ROOT / "results" / "dda_15min" / "polygon_compare"    # kept Sage outputs
MS1_SWEEP_RES = ROOT / "results" / "dda_15min" / "sweep"        # reuse its raw arm
DNOISE = (ROOT / ".." / "target" / "release" / "dnoise").resolve()
CONFIG = ROOT / "config" / "dnoise.toml"
FASTA = ROOT / "data" / "fasta" / "hybrid.fasta"
SAGE_CFG = ROOT / "config" / "sage.json"
SAGE = (ROOT / "tools" / "sage-0.15.0-beta.1"
        / "sage-v0.15.0-beta.1-x86_64-unknown-linux-gnu" / "sage")

# arm -> extra dnoise flags on top of `--config dnoise.toml` (both MS1-only).
ARMS = {"nopoly": [], "poly": ["--ms1-polygon"]}


def ms1_peaks(d: Path) -> int:
    c = sqlite3.connect(d / "analysis.tdf")
    v = c.execute(
        "SELECT COALESCE(SUM(NumPeaks),0) FROM Frames WHERE MsMsType=0").fetchone()[0]
    c.close()
    return int(v)


def bin_bytes(d: Path) -> int:
    return (d / "analysis.tdf_bin").stat().st_size


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

    raw_ms1 = sum(ms1_peaks(f) for f in files)
    raw_bytes = sum(bin_bytes(f) for f in files)
    print(f"raw: MS1 {raw_ms1/1e6:.1f}M peaks, tdf_bin {raw_bytes/1e9:.1f} GB")

    # Original (raw) baseline: reuse the identical raw search from the MS1 sweep.
    orig = CMP_RES / "original"
    if (orig / "lfq.tsv").exists():
        print("skip original (done)")
    else:
        src = MS1_SWEEP_RES / "original"
        orig.mkdir(parents=True, exist_ok=True)
        if (src / "lfq.tsv").exists() and (src / "results.sage.tsv").exists():
            for name in ("results.sage.tsv", "lfq.tsv"):
                shutil.copy(src / name, orig / name)
            print(f"reused original Sage outputs from {src}")
        else:
            t = time.time()
            run_sage(files, orig)
            print(f"original searched in {time.time()-t:.0f}s")
        (orig / "reduction.json").write_text(json.dumps(
            {"arm": "original", "polygon": False,
             "ms1_kept": raw_ms1, "ms1_raw": raw_ms1,
             "bytes_kept": raw_bytes, "bytes_raw": raw_bytes}))

    for arm, flags in ARMS.items():
        res = CMP_RES / arm
        if (res / "lfq.tsv").exists():
            print(f"skip {arm} (done)")
            continue

        t0 = time.time()
        den_dir = CMP_DATA / arm
        if den_dir.exists():
            shutil.rmtree(den_dir)
        den_dir.mkdir(parents=True, exist_ok=True)
        den_files = []
        for f in files:
            out = den_dir / f.name
            subprocess.run(
                [str(DNOISE), str(f), str(out), "--config", str(CONFIG), *flags],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            den_files.append(out)
        kept_ms1 = sum(ms1_peaks(d) for d in den_files)
        kept_bytes = sum(bin_bytes(d) for d in den_files)
        t_den = time.time() - t0

        res.mkdir(parents=True, exist_ok=True)
        (res / "reduction.json").write_text(json.dumps(
            {"arm": arm, "polygon": "--ms1-polygon" in flags,
             "ms1_kept": kept_ms1, "ms1_raw": raw_ms1,
             "bytes_kept": kept_bytes, "bytes_raw": raw_bytes}))

        run_sage(den_files, res)
        shutil.rmtree(den_dir)  # reclaim disk before the next arm
        print(f"{arm}: MS1 -{100*(1-kept_ms1/raw_ms1):.1f}%  "
              f"bytes -{100*(1-kept_bytes/raw_bytes):.1f}%  "
              f"denoise {t_den:.0f}s  total {time.time()-t0:.0f}s")

    if CMP_DATA.exists() and not any(CMP_DATA.iterdir()):
        CMP_DATA.rmdir()
    print("polygon comparison complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
