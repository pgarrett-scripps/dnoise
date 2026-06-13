#!/usr/bin/env python3
"""Figure 1 (concept) and Figure 2 (before/after) from one MS1 frame.

Consumes two CSVs (mz, one_over_k0, intensity) produced by the dnoise
`dump_frame` example for the same frame of a raw and a denoised .d. Plots use
x = m/z, y = ion mobility (1/K0), so a real ion appears as a *vertical streak*
(constant m/z across many mobility scans), matching the filter's name.

Usage: uv run scripts/08_frame_figures.py [raw.csv] [denoised.csv]
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
OUT = ROOT / "results" / "analysis"
RAW_D = Path(sys.argv[1]) if len(sys.argv) > 1 else (
    ROOT / "data/raw/LFQ_Ultra2_PASEF_5min_50ng_Condition_A_REP1.d"
)
DEN_D = Path(sys.argv[2]) if len(sys.argv) > 2 else (
    ROOT / "data/denoised/LFQ_Ultra2_PASEF_5min_50ng_Condition_A_REP1.d"
)

MZ_HALF = 5.0      # zoom is 10 m/z wide
K0_HALF = 0.075    # zoom is 0.15 1/K0 wide
CMAP = "viridis"


def most_intense_ms1_index(d: Path) -> int:
    """0-based timsrust index of the MS1 frame with the largest summed intensity."""
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


def frame_hist(ax, df, mz_edges, k0_edges, norm):
    h = ax.hist2d(
        df["mz"], df["one_over_k0"],
        bins=[mz_edges, k0_edges], weights=df["intensity"],
        norm=norm, cmap=CMAP,
    )
    ax.set_xlabel("m/z")
    ax.set_ylabel("ion mobility (1/K0)")
    return h[3]


def zoom_scatter(ax, df, mz0, k0_0, norm):
    sub = df[
        (df["mz"].between(mz0 - MZ_HALF, mz0 + MZ_HALF))
        & (df["one_over_k0"].between(k0_0 - K0_HALF, k0_0 + K0_HALF))
    ]
    sc = ax.scatter(
        sub["mz"], sub["one_over_k0"], c=sub["intensity"],
        norm=norm, cmap=CMAP, s=9, edgecolors="none",
    )
    ax.set_xlim(mz0 - MZ_HALF, mz0 + MZ_HALF)
    ax.set_ylim(k0_0 - K0_HALF, k0_0 + K0_HALF)
    ax.set_xlabel("m/z")
    ax.set_ylabel("ion mobility (1/K0)")
    return sc, sub


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    if not DUMP.is_file():
        print(f"dump_frame not built at {DUMP}\n  cargo build --release --example dump_frame")
        return 1
    idx = most_intense_ms1_index(RAW_D)
    print(f"most intense MS1 frame: 0-based index {idx}")
    with tempfile.TemporaryDirectory() as tmp:
        raw = dump_frame(RAW_D, idx, Path(tmp) / "raw.csv")
        den = dump_frame(DEN_D, idx, Path(tmp) / "den.csv")

    # Center the zoom on the most intense peak in the raw frame.
    peak = raw.loc[raw["intensity"].idxmax()]
    mz0, k0_0 = float(peak["mz"]), float(peak["one_over_k0"])
    print(f"peak: m/z={mz0:.3f}, 1/K0={k0_0:.4f}, intensity={int(peak['intensity'])}")

    vmax = float(raw["intensity"].max())
    full_norm = LogNorm(vmin=1, vmax=vmax)
    pt_norm = LogNorm(vmin=1, vmax=vmax)

    mz_edges = np.linspace(raw["mz"].min(), raw["mz"].max(), 400)
    k0_edges = np.linspace(raw["one_over_k0"].min(), raw["one_over_k0"].max(), 300)

    # ---- Figure 2: before/after (full frame + zoom) ----
    fig, ax = plt.subplots(2, 2, figsize=(11, 9))
    for col, (df, label) in enumerate([(raw, "raw"), (den, "denoised")]):
        im = frame_hist(ax[0, col], df, mz_edges, k0_edges, full_norm)
        ax[0, col].add_patch(Rectangle(
            (mz0 - MZ_HALF, k0_0 - K0_HALF), 2 * MZ_HALF, 2 * K0_HALF,
            fill=False, edgecolor="red", lw=1.2,
        ))
        ax[0, col].set_title(f"{label}  ({len(df):,} points)")
        sc, sub = zoom_scatter(ax[1, col], df, mz0, k0_0, pt_norm)
        ax[1, col].set_title(f"{label} — zoom ({len(sub):,} points)")
    fig.colorbar(im, ax=ax[0, :].tolist(), label="summed intensity", shrink=0.8)
    fig.colorbar(sc, ax=ax[1, :].tolist(), label="intensity", shrink=0.8)
    fig.suptitle("MS1 frame before/after dnoise (most intense frame)", fontsize=13)
    fig.savefig(OUT / "fig2_beforeafter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"wrote {OUT}/fig2_beforeafter.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
