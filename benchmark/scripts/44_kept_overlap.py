#!/usr/bin/env python3
"""Streak filter vs. matched intensity threshold: how different are the points they
keep, and why?

At matched MS1-point removal the two filters keep nearly the same *number* of points
but largely different *sets*. This script quantifies the overlap and characterizes
the disjoint sets mechanistically, on the native (scan, TOF-index) grid:

  - overlap        Jaccard, streak-only%, intensity-only% per run, aggregated
  - intensity      median intensity of shared / streak-only / intensity-only points
  - structure      median vertical-run length (consecutive occupied mobility scans
                   in the same m/z window, the streak prior) per partition

Run-length uses the default streak geometry (mz_half_width=3, max_internal_gap=2);
it is an independent structural descriptor of each kept point, computed from the raw
frame, not a re-run of the filter.

Reads the cached arms produced by 04_denoise.sh:
  data/<ds>/raw, data/<ds>/denoised (streak), data/<ds>/denoised_intensity.
Usage: uv run scripts/44_kept_overlap.py
"""

from __future__ import annotations

import csv
import sqlite3
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DUMP = (ROOT / ".." / "target" / "release" / "examples" / "dump_frame").resolve()
MZ_HALF, MAX_GAP = 3, 2          # default streak geometry (config/dnoise.toml)
N_FRAMES = 18                    # MS1 frames sampled per run (kept small for footprint)
RUNS = ["Condition_A_REP1", "Condition_B_REP1"]
DATASETS = {
    "dda_5min": "LFQ_Ultra2_PASEF_5min_50ng_{run}.d",
    "dda_15min": "LFQ_Ultra2_PASEF_15min_50ng_{run}.d",
}


def dump(d: Path, fid: int) -> dict[tuple[int, int], int]:
    """Return {(scan, tof): intensity} for 0-based frame fid-1."""
    out = f"/tmp/ovl_{fid}.csv"
    subprocess.run([str(DUMP), str(d), str(fid - 1), out],
                   check=True, stderr=subprocess.DEVNULL)
    pts = {}
    with open(out) as fh:
        r = csv.reader(fh)
        next(r)
        for mz, k0, inten, scan, tof in r:
            pts[(int(scan), int(tof))] = int(inten)
    Path(out).unlink(missing_ok=True)
    return pts


def run_lengths(raw: dict[tuple[int, int], int],
                query: set[tuple[int, int]]) -> list[int]:
    """Vertical-run length of each query point: occupied mobility scans (any raw
    point within +/-MZ_HALF TOF indices) in the bridged run containing its scan."""
    col2scans: dict[int, set[int]] = {}
    for (s, t) in raw:
        col2scans.setdefault(t, set()).add(s)
    col_runlen: dict[int, dict[int, int]] = {}

    def lengths_for_col(c: int) -> dict[int, int]:
        occ = set()
        for t in range(c - MZ_HALF, c + MZ_HALF + 1):
            if t in col2scans:
                occ |= col2scans[t]
        scans = sorted(occ)
        out: dict[int, int] = {}
        run = [scans[0]]
        for prev, cur in zip(scans, scans[1:]):
            if cur - prev - 1 <= MAX_GAP:
                run.append(cur)
            else:
                for s in run:
                    out[s] = len(run)
                run = [cur]
        for s in run:
            out[s] = len(run)
        return out

    res = []
    for (s, c) in query:
        if c not in col_runlen:
            col_runlen[c] = lengths_for_col(c)
        res.append(col_runlen[c][s])
    return res


def analyze(ds: str, pattern: str) -> dict:
    base = ROOT / "data" / ds
    agg = {"jac": [], "s_only": [], "i_only": [],
           "med_shared": [], "med_sonly": [], "med_ionly": [],
           "rl_shared": [], "rl_sonly": [], "rl_ionly": [],
           "kept_s": [], "kept_i": []}
    for run in RUNS:
        f = pattern.format(run=run)
        raw_d, streak_d, int_d = (base / "raw" / f, base / "denoised" / f,
                                  base / "denoised_intensity" / f)
        if not raw_d.is_dir():
            print(f"  skip {ds}/{run}: no raw")
            return {}
        con = sqlite3.connect(raw_d / "analysis.tdf")
        ids = [r[0] for r in con.execute(
            "SELECT Id FROM Frames WHERE MsMsType=0 ORDER BY Id")]
        con.close()
        sel = [ids[i] for i in np.linspace(0, len(ids) - 1, N_FRAMES).astype(int)]
        nr = ns = ni = inter = union = 0
        ints = {"shared": [], "sonly": [], "ionly": []}
        rls = {"shared": [], "sonly": [], "ionly": []}
        for fid in sel:
            raw = dump(raw_d, fid)
            S = set(dump(streak_d, fid)); I = set(dump(int_d, fid))
            shared, sonly, ionly = S & I, S - I, I - S
            nr += len(raw); ns += len(S); ni += len(I)
            inter += len(shared); union += len(S | I)
            for part, pts in (("shared", shared), ("sonly", sonly), ("ionly", ionly)):
                ints[part].extend(raw[p] for p in pts)
                rls[part].extend(run_lengths(raw, pts))
        agg["jac"].append(100 * inter / union)
        agg["s_only"].append(100 * (ns - inter) / ns)
        agg["i_only"].append(100 * (ni - inter) / ni)
        agg["kept_s"].append(100 * ns / nr)
        agg["kept_i"].append(100 * ni / nr)
        agg["med_shared"].append(np.median(ints["shared"]))
        agg["med_sonly"].append(np.median(ints["sonly"]))
        agg["med_ionly"].append(np.median(ints["ionly"]))
        agg["rl_shared"].append(np.median(rls["shared"]))
        agg["rl_sonly"].append(np.median(rls["sonly"]))
        agg["rl_ionly"].append(np.median(rls["ionly"]))
        print(f"  {ds}/{run}: Jaccard {agg['jac'][-1]:.1f}%  "
              f"streak-only {agg['s_only'][-1]:.1f}%  int-only {agg['i_only'][-1]:.1f}%")
    return {k: (float(np.mean(v)), float(np.std(v))) for k, v in agg.items()}


def main() -> int:
    for stale in Path("/tmp").glob("ovl_*.csv"):   # clear any leftovers
        stale.unlink(missing_ok=True)
    for ds, pat in DATASETS.items():
        print(f"== {ds} ==")
        r = analyze(ds, pat)
        if not r:
            continue
        m = lambda k: f"{r[k][0]:.1f} +/- {r[k][1]:.1f}"
        print(f"  kept: streak {m('kept_s')}%  intensity {m('kept_i')}%")
        print(f"  Jaccard {m('jac')}%  | streak-only {m('s_only')}%  | int-only {m('i_only')}%")
        print(f"  median intensity: shared {m('med_shared')}  streak-only {m('med_sonly')}  int-only {m('med_ionly')}")
        print(f"  median run length: shared {m('rl_shared')}  streak-only {m('rl_sonly')}  int-only {m('rl_ionly')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
