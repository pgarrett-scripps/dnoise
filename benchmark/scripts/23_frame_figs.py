#!/usr/bin/env python3
"""Frame-level figures on the most-intense MS1 frame of Condition_A_REP1.

Produces, into paper/figures/:
  fig2_frame_stages.png : raw -> vertical filter -> +halo (the two-stage MS1 filter)
  fig5_frames.png       : raw -> watershed -> box (the two centroiders)

Same frame and zoom as the original before/after figure (most intense MS1 frame
of the raw run; zoom centered on its most intense peak). Dumps each arm's frame
with the `dump_frame` example and plots full-frame density (top) + zoom scatter
(bottom) per column.
"""

from __future__ import annotations

import sqlite3
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

# Keep all lettering >= 4.5 pt after the figure is downscaled to the text width.
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "figure.titlesize": 15,
})

ROOT = Path(__file__).resolve().parents[1]
FIGS = ROOT.parent / "paper" / "figures"
DUMP = ROOT.parent / "target" / "release" / "examples" / "dump_frame"
DATA = ROOT / "data" / "dda_5min"
F = "LFQ_Ultra2_PASEF_5min_50ng_Condition_A_REP1.d"
MZ_HALF, K0_HALF, CMAP = 5.0, 0.075, "viridis"


def most_intense_ms1_index(d: Path) -> int:
    c = sqlite3.connect(d / "analysis.tdf")
    ids = [r[0] for r in c.execute("SELECT Id FROM Frames ORDER BY Id")]
    best = c.execute("SELECT Id FROM Frames WHERE MsMsType=0 "
                     "ORDER BY SummedIntensities DESC LIMIT 1").fetchone()[0]
    c.close()
    return ids.index(best)


def dump(d: Path, idx: int, out: Path) -> pd.DataFrame:
    subprocess.run([str(DUMP), str(d), str(idx), str(out)], check=True, stderr=subprocess.DEVNULL)
    return pd.read_csv(out)


def frame_fig(cols: list[tuple[str, Path]], out: Path, title: str) -> None:
    raw_d = cols[0][1]
    idx = most_intense_ms1_index(raw_d)
    with tempfile.TemporaryDirectory() as tmp:
        frames = [(lab, dump(d, idx, Path(tmp) / f"{lab}.csv")) for lab, d in cols]
    raw = frames[0][1]
    peak = raw.loc[raw["intensity"].idxmax()]
    mz0, k0 = float(peak["mz"]), float(peak["one_over_k0"])
    vmax = float(raw["intensity"].max())
    norm = LogNorm(vmin=1, vmax=vmax)
    mz_e = np.linspace(raw["mz"].min(), raw["mz"].max(), 400)
    k0_e = np.linspace(raw["one_over_k0"].min(), raw["one_over_k0"].max(), 300)

    n = len(cols)
    # Extra hspace so the top row's "m/z" xlabel does not collide with the
    # bottom row's titles.
    fig, ax = plt.subplots(2, n, figsize=(4.6 * n, 9.2),
                           gridspec_kw={"hspace": 0.45})
    im = sc = None
    for j, (lab, df) in enumerate(frames):
        im = ax[0, j].hist2d(df["mz"], df["one_over_k0"], bins=[mz_e, k0_e],
                             weights=df["intensity"], norm=norm, cmap=CMAP)[3]
        ax[0, j].add_patch(Rectangle((mz0 - MZ_HALF, k0 - K0_HALF), 2 * MZ_HALF, 2 * K0_HALF,
                                     fill=False, edgecolor="red", lw=1.2))
        ax[0, j].set_title(f"{lab}  ({len(df):,} pts)")
        ax[0, j].set_xlabel("m/z"); ax[0, j].set_ylabel("ion mobility (1/K0)")
        sub = df[df["mz"].between(mz0 - MZ_HALF, mz0 + MZ_HALF)
                 & df["one_over_k0"].between(k0 - K0_HALF, k0 + K0_HALF)]
        sc = ax[1, j].scatter(sub["mz"], sub["one_over_k0"], c=sub["intensity"],
                              norm=norm, cmap=CMAP, s=9, edgecolors="none")
        ax[1, j].set_xlim(mz0 - MZ_HALF, mz0 + MZ_HALF); ax[1, j].set_ylim(k0 - K0_HALF, k0 + K0_HALF)
        ax[1, j].set_title(f"{lab}, zoom ({len(sub):,} pts)")
        ax[1, j].set_xlabel("m/z"); ax[1, j].set_ylabel("ion mobility (1/K0)")
    fig.colorbar(im, ax=ax[0, :].tolist(), label="summed intensity", shrink=0.8)
    fig.colorbar(sc, ax=ax[1, :].tolist(), label="intensity", shrink=0.8)
    fig.suptitle(title, fontsize=13)
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    if not DUMP.is_file():
        print(f"dump_frame not built: {DUMP}\n  cargo build --release --example dump_frame")
        return 1
    frame_fig(
        [("raw", DATA / "raw" / F),
         ("vertical filter", DATA / "_fig_vertonly.d"),
         ("+ halo filter", DATA / "denoised" / F)],  # benchmark default (halo_peak_fraction 0.15)
        FIGS / "fig2_frame_stages.png",
        "MS1 frame through the two-stage filter (most intense MS1 frame)",
    )
    frame_fig(
        [("raw", DATA / "raw" / F),
         ("watershed", DATA / "watershed" / F),
         ("box", DATA / "box" / F)],
        FIGS / "fig5_frames.png",
        "MS1 frame after the two centroiders (watershed collapses; box tiles)",
    )
    frame_fig(
        [("raw", DATA / "raw" / F),
         ("streak filter", DATA / "denoised" / F),
         ("intensity threshold", DATA / "denoised_intensity" / F)],
        FIGS / "si_intensity_frame.png",
        "MS1 frame: streak filter vs. a matched strict intensity threshold "
        "(equal MS1-point removal)",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
