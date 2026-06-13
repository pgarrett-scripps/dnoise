#!/usr/bin/env python3
"""Figure 3: compression parameter-interaction grid on a single file.

Sweeps the vertical filter's `min_feature_length` (x-axis) against the pairing of
halo on/off and `max_internal_gap` (y-axis), with all other parameters at the
benchmark defaults (`mz_half_width=3`, `iterations=2`). For every cell it runs the
dnoise binary on one raw .d and measures the resulting data reduction. Two
heatmaps share the grid: % MS1 peaks removed and % analysis.tdf_bin bytes
removed. No search/quant — this characterizes compression only.

The grid shows how the knobs interact: vertical colour bands mean a knob is inert
given the others; gradients along both axes mean the knobs are complementary; the
halo on/off rows show the halo's (roughly constant) contribution on top of the
vertical filter.

Usage: uv run scripts/09_compression_sweep.py [raw.d]
"""

from __future__ import annotations

import os
import shutil
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
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
DNOISE = REPO / "target" / "release" / "dnoise"
# Use the benchmark config as the base (mz_half_width, halo params, iterations)
# so the grid reflects the real pipeline; the sweep only overrides the knobs on
# the two axes (min_feature_length, max_internal_gap) and the halo on/off flag.
CONFIG = ROOT / "config" / "dnoise.toml"
DATASET = os.environ.get("DATASET", "dda_5min")
OUT = ROOT / "results" / DATASET / "analysis"
RAW = Path(sys.argv[1]) if len(sys.argv) > 1 else (
    ROOT / "data" / DATASET / "raw" / "LFQ_Ultra2_PASEF_5min_50ng_Condition_A_REP1.d"
)

LENGTHS = [3, 5, 7, 11]          # x-axis: min_feature_length
GAPS = [0, 1, 2]                 # max_internal_gap, paired with halo on the y-axis
HALOS = [("off", ["--no-halo"]), ("on", [])]
ITERATIONS = 2                   # held at the benchmark default
DEFAULT = {"length": 5, "gap": 1, "halo": "on"}  # boxed reference cell


def measure(d: Path) -> tuple[int, int]:
    """(MS1 peaks, tdf_bin bytes) for a .d folder."""
    c = sqlite3.connect(d / "analysis.tdf")
    ms1 = c.execute("SELECT COALESCE(SUM(NumPeaks),0) FROM Frames WHERE MsMsType=0").fetchone()[0]
    c.close()
    return int(ms1), (d / "analysis.tdf_bin").stat().st_size


def heatmap(ax, grid, row_labels, title, cbar_label, fig, default_rc):
    im = ax.imshow(grid, cmap="viridis", aspect="auto", vmin=0, vmax=100, origin="upper")
    ax.set_xticks(range(len(LENGTHS)))
    ax.set_xticklabels(LENGTHS)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels)
    ax.set_xlabel("min_feature_length (streak span, scans)")
    # Annotate each cell; pick text colour for contrast against viridis.
    for r in range(grid.shape[0]):
        for cc in range(grid.shape[1]):
            v = grid[r, cc]
            ax.text(cc, r, f"{v:.0f}", ha="center", va="center", fontsize=8,
                    color="white" if v < 55 else "black")
    # Divider between the halo-off and halo-on blocks.
    ax.axhline(len(GAPS) - 0.5, color="white", lw=2)
    # Outline the default cell.
    dr, dc = default_rc
    ax.add_patch(Rectangle((dc - 0.5, dr - 0.5), 1, 1, fill=False, edgecolor="red", lw=2.2))
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=cbar_label, shrink=0.85)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    if not DNOISE.is_file():
        print(f"dnoise binary not found at {DNOISE} (run: cargo build --release)")
        return 1
    if not RAW.is_dir():
        print(f"raw .d not found: {RAW}")
        return 1

    base_ms1, base_bytes = measure(RAW)

    # Row order: halo off (gaps 0,1,2) then halo on (gaps 0,1,2).
    row_keys = [(h, g) for (h, _) in HALOS for g in GAPS]
    row_labels = [f"halo {h} · gap {g}" for (h, g) in row_keys]
    ms1_grid = np.full((len(row_keys), len(LENGTHS)), np.nan)
    bytes_grid = np.full((len(row_keys), len(LENGTHS)), np.nan)

    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for (halo_name, halo_flags) in HALOS:
            for g in GAPS:
                for L in LENGTHS:
                    out = Path(tmp) / f"h{halo_name}_g{g}_L{L}.d"
                    subprocess.run(
                        [str(DNOISE), str(RAW), str(out),
                         "--config", str(CONFIG),
                         "--min-feature-length", str(L),
                         "--max-internal-gap", str(g),
                         "--iterations", str(ITERATIONS),
                         *halo_flags, "--force"],
                        check=True, stdout=subprocess.DEVNULL,
                    )
                    ms1, nbytes = measure(out)
                    shutil.rmtree(out, ignore_errors=True)
                    ms1_pct = 100 * (1 - ms1 / base_ms1)
                    bytes_pct = 100 * (1 - nbytes / base_bytes)
                    ri = row_keys.index((halo_name, g))
                    ci = LENGTHS.index(L)
                    ms1_grid[ri, ci] = ms1_pct
                    bytes_grid[ri, ci] = bytes_pct
                    rows.append({
                        "halo": halo_name, "max_internal_gap": g, "min_feature_length": L,
                        "ms1_removed_pct": ms1_pct, "bytes_removed_pct": bytes_pct,
                    })
                    print(f"halo {halo_name:3s} gap {g} len {L:2d}: "
                          f"MS1 -{ms1_pct:.1f}%  bytes -{bytes_pct:.1f}%")

    pd.DataFrame(rows).to_csv(OUT / "compression_sweep.csv", index=False)

    default_rc = (row_keys.index((DEFAULT["halo"], DEFAULT["gap"])),
                  LENGTHS.index(DEFAULT["length"]))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    heatmap(axes[0], ms1_grid, row_labels, "MS1 peaks removed", "% MS1 removed", fig, default_rc)
    heatmap(axes[1], bytes_grid, [""] * len(row_labels), "tdf_bin bytes removed",
            "% bytes removed", fig, default_rc)
    axes[0].set_ylabel("halo state · max_internal_gap")
    fig.suptitle("dnoise compression: filter parameter interactions (single 5-min run; "
                 "red = default)", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "compression_sweep.png", dpi=150, bbox_inches="tight")
    if (REPO / "paper" / "figures").is_dir():
        fig.savefig(REPO / "paper" / "figures" / "compression_sweep.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}/compression_sweep.png and compression_sweep.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
