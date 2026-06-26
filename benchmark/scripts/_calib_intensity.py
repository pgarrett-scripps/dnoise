#!/usr/bin/env python3
"""Calibrate per-point MS1 intensity thresholds T for target removal fractions.

The intensity-threshold control keeps a point iff its own intensity >= T (see
config/dnoise.intensity.toml). To compare it against the streak filter at matched
data reduction we need T values that remove ~20/40/60/80% of MS1 points. This
samples MS1 frames from one replicate, builds the intensity distribution, and
prints the smallest integer T whose removed fraction reaches each target.
"""

from __future__ import annotations

import sqlite3
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REP = ROOT / "data/dda_15min/raw/LFQ_Ultra2_PASEF_15min_50ng_Condition_A_REP1.d"
DUMP = (ROOT / ".." / "target" / "release" / "examples" / "dump_frame").resolve()
TARGETS = [0.2, 0.4, 0.6, 0.8]
N_SAMPLE = 24

con = sqlite3.connect(REP / "analysis.tdf")
ids = [r[0] for r in con.execute(
    "SELECT Id FROM Frames WHERE MsMsType=0 ORDER BY Id")]
con.close()
sel = [ids[i] for i in np.linspace(0, len(ids) - 1, N_SAMPLE).astype(int)]

chunks = []
for fid in sel:
    out = f"/tmp/cf_{fid}.csv"
    subprocess.run([str(DUMP), str(REP), str(fid - 1), out],
                   check=True, stderr=subprocess.DEVNULL)
    chunks.append(np.loadtxt(out, delimiter=",", skiprows=1, usecols=2))
    Path(out).unlink(missing_ok=True)

v = np.concatenate(chunks)
v.sort()
print(f"sampled {len(sel)} MS1 frames, {len(v)} points; "
      f"intensity min={int(v.min())} median={int(np.median(v))} max={int(v.max())}")
uniq = np.unique(v)
print("target  T   actual_removed")
for tgt in TARGETS:
    T = next((int(t) for t in uniq if (v < t).mean() >= tgt), int(uniq[-1]) + 1)
    print(f"{tgt:.2f}   {T:>4d}  {(v < T).mean():.3f}")
