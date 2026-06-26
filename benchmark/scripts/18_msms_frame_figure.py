#!/usr/bin/env python3
"""MS/MS frame before/after dnoise, for the SI (DDA and DIA).

Picks the most intense MS/MS frame of a run, dumps it from the raw and the
MS/MS-denoised (`denoised_msms`) `.d` with the dnoise `dump_frame` example, and
plots before/after as point density over m/z x ion mobility (1/K0) plus a zoom
on the most intense fragment peak.

Works for ddaPASEF (MS/MS denoised per precursor across re-isolations) and
diaPASEF (MS/MS denoised whole-frame); the DATASET selects which.

Usage: uv run scripts/18_msms_frame_figure.py [dataset] [run_basename]
       dataset      : dda_5min (default) | dia_5min | ...
       run_basename : a .d folder name under data/<dataset>/raw (default: first
                      Condition_A_REP1 match, else first sorted)
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

# Keep all lettering >= 4.5 pt once the figure is downscaled to the text width.
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "figure.titlesize": 15,
})
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
DUMP = REPO / "target" / "release" / "examples" / "dump_frame"
MZ_HALF = 5.0      # zoom 10 m/z wide
K0_HALF = 0.075    # zoom 0.15 1/K0 wide
CMAP = "viridis"


def acq_label(dataset: str) -> str:
    return "diaPASEF" if "dia" in dataset else "ddaPASEF"


def pick_run(raw_dir: Path, want: str | None) -> Path:
    runs = sorted(raw_dir.glob("*.d"))
    if not runs:
        raise SystemExit(f"no raw .d in {raw_dir}")
    if want:
        for r in runs:
            if r.name == want or r.stem == want:
                return r
    for r in runs:
        if "Condition_A_REP1" in r.name:
            return r
    return runs[0]


def most_intense_msms_index(d: Path) -> tuple[int, int]:
    """(0-based timsrust index, Id) of the MS/MS frame with the largest summed intensity."""
    c = sqlite3.connect(d / "analysis.tdf")
    best = c.execute(
        "SELECT Id FROM Frames WHERE MsMsType!=0 ORDER BY SummedIntensities DESC LIMIT 1"
    ).fetchone()[0]
    idx = c.execute("SELECT COUNT(*) FROM Frames WHERE Id<?", (best,)).fetchone()[0]
    c.close()
    return idx, best


def dump_frame(d: Path, idx: int, out: Path) -> pd.DataFrame:
    subprocess.run([str(DUMP), str(d), str(idx), str(out)], check=True, stderr=subprocess.DEVNULL)
    return pd.read_csv(out)


def frame_hist(ax, df, mz_edges, k0_edges, norm):
    h = ax.hist2d(df["mz"], df["one_over_k0"], bins=[mz_edges, k0_edges],
                  weights=df["intensity"], norm=norm, cmap=CMAP)
    ax.set_xlabel("m/z")
    ax.set_ylabel("ion mobility (1/K0)")
    return h[3]


def zoom_scatter(ax, df, mz0, k0_0, norm):
    sub = df[(df["mz"].between(mz0 - MZ_HALF, mz0 + MZ_HALF))
             & (df["one_over_k0"].between(k0_0 - K0_HALF, k0_0 + K0_HALF))]
    sc = ax.scatter(sub["mz"], sub["one_over_k0"], c=sub["intensity"],
                    norm=norm, cmap=CMAP, s=9, edgecolors="none")
    ax.set_xlim(mz0 - MZ_HALF, mz0 + MZ_HALF)
    ax.set_ylim(k0_0 - K0_HALF, k0_0 + K0_HALF)
    ax.set_xlabel("m/z")
    ax.set_ylabel("ion mobility (1/K0)")
    return sc, sub


def main() -> int:
    dataset = sys.argv[1] if len(sys.argv) > 1 else "dda_5min"
    want = sys.argv[2] if len(sys.argv) > 2 else None
    if not DUMP.is_file():
        print(f"dump_frame not built at {DUMP}\n  cargo build --release --example dump_frame")
        return 1

    raw_dir = ROOT / "data" / dataset / "raw"
    msms_dir = ROOT / "data" / dataset / "denoised_msms"
    out_dir = ROOT / "results" / dataset / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    run = pick_run(raw_dir, want)
    raw_d = run
    msms_d = msms_dir / run.name
    if not msms_d.is_dir():
        print(f"no MS/MS-denoised counterpart: {msms_d}\n  run denoise with --denoise-msms")
        return 1

    idx, fid = most_intense_msms_index(raw_d)
    print(f"[{dataset}] {run.name}: most intense MS/MS frame Id={fid} (0-based idx {idx})")
    with tempfile.TemporaryDirectory() as tmp:
        raw = dump_frame(raw_d, idx, Path(tmp) / "raw.csv")
        den = dump_frame(msms_d, idx, Path(tmp) / "den.csv")

    peak = raw.loc[raw["intensity"].idxmax()]
    mz0, k0_0 = float(peak["mz"]), float(peak["one_over_k0"])
    print(f"  raw {len(raw):,} pts -> denoised {len(den):,} pts ({100*len(den)/max(len(raw),1):.0f}% kept)")

    vmax = float(raw["intensity"].max())
    norm = LogNorm(vmin=1, vmax=vmax)
    mz_edges = np.linspace(raw["mz"].min(), raw["mz"].max(), 400)
    k0_edges = np.linspace(raw["one_over_k0"].min(), raw["one_over_k0"].max(), 300)

    fig, ax = plt.subplots(2, 2, figsize=(11, 9))
    for col, (df, label) in enumerate([(raw, "raw"), (den, "denoised")]):
        im = frame_hist(ax[0, col], df, mz_edges, k0_edges, norm)
        ax[0, col].add_patch(Rectangle((mz0 - MZ_HALF, k0_0 - K0_HALF), 2 * MZ_HALF, 2 * K0_HALF,
                                       fill=False, edgecolor="red", lw=1.2))
        ax[0, col].set_title(f"{label}  ({len(df):,} points)")
        sc, sub = zoom_scatter(ax[1, col], df, mz0, k0_0, norm)
        ax[1, col].set_title(f"{label}, zoom ({len(sub):,} points)")
    fig.colorbar(im, ax=ax[0, :].tolist(), label="summed intensity", shrink=0.8)
    fig.colorbar(sc, ax=ax[1, :].tolist(), label="intensity", shrink=0.8)
    fig.suptitle(f"{acq_label(dataset)} MS/MS frame before/after dnoise "
                 f"(most intense MS/MS frame, {dataset})", fontsize=13)
    out = out_dir / "msms_beforeafter.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
