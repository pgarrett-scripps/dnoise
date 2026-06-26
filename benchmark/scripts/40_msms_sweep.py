#!/usr/bin/env python3
"""MS/MS-denoising parameter sweep on one dda_15min sample with replicates.

Companion to 30_param_sweep.py, but for the ddaPASEF MS/MS spectral denoiser
(`--denoise-msms`) instead of the MS1 streak filter. MS1 denoising leaves DDA
identifications unchanged (it only touches MS1 frames), so that sweep was scored
on quant. MS/MS denoising rewrites the fragment spectra (it combines each
precursor's re-isolated scans, filters the combined spectrum, and prunes the
individual scans), so it directly changes Sage identifications - which is the
question here: can we compress/clean MS/MS without losing IDs?

Picks one HYE condition (default A, 6 replicates) and sweeps the two MS/MS knobs,
`msms_max_internal_gap` x `msms_min_feature_length`, re-denoising the replicates
from raw for each setting with the benchmark MS1 pipeline held fixed
(`--config dnoise.toml --ms1-polygon`) plus `--denoise-msms`. MS1 denoising is
ID-neutral and contributes nothing to MS2 data reduction, so the ID and MS2-byte
deltas vs. the raw `original` arm isolate the MS/MS denoiser. 41_msms_sweep_analyze.py
scores each arm (IDs + MS2 reduction + quant/fidelity).

Disk-safe and idempotent, exactly like 30_param_sweep.py: each arm is denoised,
searched, then its (large) denoised .d folders are deleted, keeping only the
small Sage TSVs. An arm whose results/.../msms_sweep/<arm>/lfq.tsv already exists
is skipped, so the sweep resumes after an interruption.

Outputs (kept):
  results/dda_15min/msms_sweep/original/      Sage outputs for the raw arm
  results/dda_15min/msms_sweep/g<gap>_l<len>/ Sage outputs + reduction.json per arm
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
SWEEP_DATA = ROOT / "data" / "dda_15min" / "msms_sweep"      # temp denoised arms
SWEEP_RES = ROOT / "results" / "dda_15min" / "msms_sweep"    # kept Sage outputs
MS1_SWEEP_RES = ROOT / "results" / "dda_15min" / "sweep"     # reuse its raw arm
DNOISE = (ROOT / ".." / "target" / "release" / "dnoise").resolve()
CONFIG = ROOT / "config" / "dnoise.toml"
FASTA = ROOT / "data" / "fasta" / "hybrid.fasta"
SAGE_CFG = ROOT / "config" / "sage.json"
SAGE = (ROOT / "tools" / "sage-0.15.0-beta.1"
        / "sage-v0.15.0-beta.1-x86_64-unknown-linux-gnu" / "sage")

# 3 x 4 grid over the two MS/MS knobs (msms_mz_half_width=3, msms_iterations=1
# from config). MS/MS isolation windows are short (~25 scans, combined across
# re-isolations), so feature lengths stay small.
GAPS = [2, 5, 8]
LENS = [1, 2, 3, 4]
GRID = [(g, l) for g in GAPS for l in LENS]


def peak_counts(d: Path) -> tuple[int, int]:
    """(MS1 peaks, MS/MS peaks) for a .d folder."""
    c = sqlite3.connect(d / "analysis.tdf")
    ms1 = c.execute(
        "SELECT COALESCE(SUM(NumPeaks),0) FROM Frames WHERE MsMsType=0").fetchone()[0]
    ms2 = c.execute(
        "SELECT COALESCE(SUM(NumPeaks),0) FROM Frames WHERE MsMsType!=0").fetchone()[0]
    c.close()
    return int(ms1), int(ms2)


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
    print(f"grid: {len(GRID)} arms (msms gap {GAPS} x len {LENS}) + original")

    raw_ms1, raw_ms2 = 0, 0
    for f in files:
        m1, m2 = peak_counts(f)
        raw_ms1 += m1
        raw_ms2 += m2
    print(f"raw peaks (6 reps): MS1 {raw_ms1/1e6:.1f}M  MS/MS {raw_ms2/1e6:.1f}M")

    # Original (raw) baseline. Reuse the identical raw search from the MS1 sweep if
    # present (same 6 files, same Sage config/FASTA) to save one search; else run.
    orig = SWEEP_RES / "original"
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
            {"arm": "original", "kind": "msms", "gap": None, "len": None,
             "ms1_kept": raw_ms1, "ms1_raw": raw_ms1,
             "ms2_kept": raw_ms2, "ms2_raw": raw_ms2}))

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
                 "--ms1-polygon",       # DDA benchmark MS1 gate (held fixed)
                 "--denoise-msms",      # the feature under test
                 "--msms-max-internal-gap", str(g),
                 "--msms-min-feature-length", str(l)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            den_files.append(out)
        kept_ms1, kept_ms2 = 0, 0
        for d in den_files:
            m1, m2 = peak_counts(d)
            kept_ms1 += m1
            kept_ms2 += m2
        t_den = time.time() - t0

        res.mkdir(parents=True, exist_ok=True)
        (res / "reduction.json").write_text(json.dumps(
            {"arm": arm, "kind": "msms", "gap": g, "len": l,
             "ms1_kept": kept_ms1, "ms1_raw": raw_ms1,
             "ms2_kept": kept_ms2, "ms2_raw": raw_ms2}))

        run_sage(den_files, res)
        shutil.rmtree(den_dir)  # reclaim disk before the next arm
        red2 = 100 * (1 - kept_ms2 / raw_ms2)
        print(f"{arm}: MS/MS -{red2:.1f}%  denoise {t_den:.0f}s  "
              f"total {time.time()-t0:.0f}s")

    if SWEEP_DATA.exists() and not any(SWEEP_DATA.iterdir()):
        SWEEP_DATA.rmdir()
    print("msms sweep complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
