#!/usr/bin/env python3
"""Compare denoisers on the SAME frame + precursor as Figure 1/2.

Uses the identical frame selection and peak-centering as 08_frame_figures.py
(most intense MS1 frame of the RAW .d; zoom centered on its most intense peak),
then dumps that frame from RAW and each denoised .d and plots full-frame + zoom
side by side — so dnoise and Bruker Minesweeper can be compared on identical data.

Default arms (same Condition_A_REP1 used in Fig 1):
  raw  |  dnoise_ms1  |  bruker_minesweeper_mf30  (~matched MS1 retention)

Usage: uv run scripts/20_frame_compare.py [raw.d label1=den1.d label2=den2.d ...]
"""

from __future__ import annotations

import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
DUMP = REPO / "target" / "release" / "examples" / "dump_frame"
OUT = ROOT / "results" / "single_compare"
F = "LFQ_Ultra2_PASEF_5min_50ng_Condition_A_REP1.d"

MZ_HALF = 5.0      # zoom 10 m/z wide  (matches 08_frame_figures.py)
K0_HALF = 0.075    # zoom 0.15 1/K0 wide
CMAP = "viridis"

# (label, .d path) — raw first.
RAW_D = ROOT / "data/dda_5min/raw" / F
ARMS = [
    ("dnoise_ms1", ROOT / "data/dda_5min/denoised" / F),
    ("bruker_minesweeper_mf30", ROOT / "data/dda_5min/bruker_demo/ms_mf30" / F),
]


def most_intense_ms1_index(d: Path) -> int:
    c = sqlite3.connect(d / "analysis.tdf")
    ids = [r[0] for r in c.execute("SELECT Id FROM Frames ORDER BY Id")]
    best = c.execute(
        "SELECT Id FROM Frames WHERE MsMsType=0 ORDER BY SummedIntensities DESC LIMIT 1"
    ).fetchone()[0]
    c.close()
    return ids.index(best)


def dump_frame(d: Path, idx: int, out: Path) -> pd.DataFrame:
    subprocess.run([str(DUMP), str(d), str(idx), str(out)], check=True, stderr=subprocess.DEVNULL)
    return pd.read_csv(out)


def main() -> int:
    if not DUMP.is_file():
        print(f"dump_frame not built: {DUMP}\n  cargo build --release --example dump_frame")
        return 1
    OUT.mkdir(parents=True, exist_ok=True)

    idx = most_intense_ms1_index(RAW_D)
    cols = [("raw", RAW_D)] + ARMS
    print(f"frame: 0-based index {idx}  (most intense MS1 frame of raw)")

    with tempfile.TemporaryDirectory() as tmp:
        frames = []
        for label, d in cols:
            df = dump_frame(d, idx, Path(tmp) / f"{label}.csv")
            frames.append((label, df))

    raw = frames[0][1]
    peak = raw.loc[raw["intensity"].idxmax()]
    mz0, k0_0 = float(peak["mz"]), float(peak["one_over_k0"])
    print(f"precursor (zoom center): m/z={mz0:.3f}, 1/K0={k0_0:.4f}, "
          f"intensity={int(peak['intensity'])}")

    vmax = float(raw["intensity"].max())
    norm = LogNorm(vmin=1, vmax=vmax)
    mz_edges = np.linspace(raw["mz"].min(), raw["mz"].max(), 400)
    k0_edges = np.linspace(raw["one_over_k0"].min(), raw["one_over_k0"].max(), 300)

    n = len(cols)
    fig, ax = plt.subplots(2, n, figsize=(5 * n, 9))
    im = sc = None
    for j, (label, df) in enumerate(frames):
        # full frame
        im = ax[0, j].hist2d(df["mz"], df["one_over_k0"], bins=[mz_edges, k0_edges],
                             weights=df["intensity"], norm=norm, cmap=CMAP)[3]
        ax[0, j].add_patch(Rectangle((mz0 - MZ_HALF, k0_0 - K0_HALF),
                                     2 * MZ_HALF, 2 * K0_HALF, fill=False,
                                     edgecolor="red", lw=1.2))
        ax[0, j].set_title(f"{label}  ({len(df):,} pts)")
        ax[0, j].set_xlabel("m/z"); ax[0, j].set_ylabel("ion mobility (1/K0)")
        # zoom
        sub = df[df["mz"].between(mz0 - MZ_HALF, mz0 + MZ_HALF)
                 & df["one_over_k0"].between(k0_0 - K0_HALF, k0_0 + K0_HALF)]
        sc = ax[1, j].scatter(sub["mz"], sub["one_over_k0"], c=sub["intensity"],
                              norm=norm, cmap=CMAP, s=9, edgecolors="none")
        ax[1, j].set_xlim(mz0 - MZ_HALF, mz0 + MZ_HALF)
        ax[1, j].set_ylim(k0_0 - K0_HALF, k0_0 + K0_HALF)
        ax[1, j].set_title(f"{label} — zoom ({len(sub):,} pts)")
        ax[1, j].set_xlabel("m/z"); ax[1, j].set_ylabel("ion mobility (1/K0)")

    fig.colorbar(im, ax=ax[0, :].tolist(), label="summed intensity", shrink=0.8)
    fig.colorbar(sc, ax=ax[1, :].tolist(), label="intensity", shrink=0.8)
    fig.suptitle(f"Same MS1 frame/precursor — denoiser comparison "
                 f"(frame idx {idx})", fontsize=13)
    out = OUT / "frame_compare_A_REP1.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
