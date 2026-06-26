#!/usr/bin/env python3
"""Figure 3: compression parameter-interaction grid on a single file.

Top row (two heatmaps): sweeps the vertical filter's `min_feature_length`
(x-axis) against the pairing of halo on/off and `max_internal_gap` (y-axis), with
all other parameters at the benchmark defaults (`mz_half_width=3`,
`iterations=2`). For every cell it runs the dnoise binary on one raw .d and
measures the resulting data reduction: % MS1 peaks removed and % analysis.tdf_bin
bytes removed.

Bottom panel: holds every other parameter at the benchmark default and sweeps the
streak filter's intensity threshold `min_feature_intensity` (the total
summed-intensity floor a vertical run must clear to be kept), plotting % MS1 peaks
and % bytes removed against it. The default is 0 (no intensity floor), so this
panel shows how an added intensity gate raises removal beyond the geometric knobs.

No search/quant — this characterizes compression only. The grid shows how the
knobs interact: vertical colour bands mean a knob is inert given the others;
gradients along both axes mean the knobs are complementary; the halo on/off rows
show the halo's (roughly constant) contribution on top of the vertical filter.

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

# Keep all lettering >= 4.5 pt once the figure is downscaled to the text width.
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "figure.titlesize": 15,
})

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

LENGTHS = [3, 5, 7, 9, 11]       # x-axis: min_feature_length
GAPS = [0, 1, 2, 3]              # max_internal_gap, paired with halo on the y-axis
HALOS = [("off", ["--no-halo"]), ("on", [])]
ITERATIONS = 2                   # held at the benchmark default
DEFAULT = {"length": 5, "gap": 2, "halo": "on"}  # boxed reference cell
# Bottom panel: streak filter's per-feature intensity floor, swept with every
# other knob at the benchmark default. 0 is the benchmark default (no floor).
INTENSITIES = [0, 50, 100, 200, 500, 1000, 2000, 5000]


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
            ax.text(cc, r, f"{v:.0f}", ha="center", va="center", fontsize=11,
                    color="white" if v < 55 else "black")
    # Divider between the halo-off and halo-on blocks.
    ax.axhline(len(GAPS) - 0.5, color="white", lw=2)
    # Outline the default cell.
    dr, dc = default_rc
    ax.add_patch(Rectangle((dc - 0.5, dr - 0.5), 1, 1, fill=False, edgecolor="red", lw=2.2))
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=cbar_label, shrink=0.85)


def line_panel(ax, intensities, ms1_pct, bytes_pct, default_intensity):
    """Removal vs. the streak filter's min_feature_intensity floor."""
    xs = range(len(intensities))
    ax.plot(xs, ms1_pct, "o-", color="#1f77b4", lw=2, ms=7, label="MS1 peaks removed")
    ax.plot(xs, bytes_pct, "s-", color="#d62728", lw=2, ms=7, label="tdf_bin bytes removed")
    ax.set_xticks(list(xs))
    ax.set_xticklabels([str(i) for i in intensities])
    ax.set_xlabel("min_feature_intensity (streak total-intensity floor)")
    ax.set_ylabel("% removed")
    ax.set_ylim(0, 100)
    ax.grid(True, axis="y", alpha=0.3)
    # Mark the benchmark default (no intensity floor).
    di = intensities.index(default_intensity)
    ax.axvline(di, color="red", ls="--", lw=2, label="default (no floor)")
    ax.legend(fontsize=10, loc="center right")
    ax.set_title("Streak intensity floor (other knobs at default)")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    if not DNOISE.is_file():
        print(f"dnoise binary not found at {DNOISE} (run: cargo build --release)")
        return 1
    if not RAW.is_dir():
        print(f"raw .d not found: {RAW}")
        return 1

    # Row order: halo off (gaps 0,1,2) then halo on (gaps 0,1,2).
    row_keys = [(h, g) for (h, _) in HALOS for g in GAPS]
    row_labels = [f"halo {h} · gap {g}" for (h, g) in row_keys]
    ms1_grid = np.full((len(row_keys), len(LENGTHS)), np.nan)
    bytes_grid = np.full((len(row_keys), len(LENGTHS)), np.nan)

    # The sweep runs dnoise once per cell, which is slow; set SWEEP_REPLOT=1 to
    # rebuild the figure from the cached compression_sweep.csv (e.g. to restyle it
    # without re-running the sweep).
    csv_path = OUT / "compression_sweep.csv"
    int_csv_path = OUT / "compression_intensity_sweep.csv"
    int_ms1 = [np.nan] * len(INTENSITIES)   # bottom panel: MS1 removal vs. floor
    int_bytes = [np.nan] * len(INTENSITIES)
    if os.environ.get("SWEEP_REPLOT") and csv_path.is_file() and int_csv_path.is_file():
        df = pd.read_csv(csv_path)
        for r in df.itertuples():
            ri = row_keys.index((r.halo, r.max_internal_gap))
            ci = LENGTHS.index(r.min_feature_length)
            ms1_grid[ri, ci] = r.ms1_removed_pct
            bytes_grid[ri, ci] = r.bytes_removed_pct
        di = pd.read_csv(int_csv_path)
        for r in di.itertuples():
            ii = INTENSITIES.index(r.min_feature_intensity)
            int_ms1[ii] = r.ms1_removed_pct
            int_bytes[ii] = r.bytes_removed_pct
        print(f"replotting from {csv_path} and {int_csv_path} (SWEEP_REPLOT set; sweep skipped)")
    else:
        base_ms1, base_bytes = measure(RAW)
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
            # Bottom panel: sweep the streak intensity floor, all else at default.
            int_rows = []
            for ii, I in enumerate(INTENSITIES):
                out = Path(tmp) / f"int_{I}.d"
                subprocess.run(
                    [str(DNOISE), str(RAW), str(out),
                     "--config", str(CONFIG),
                     "--min-feature-length", str(DEFAULT["length"]),
                     "--max-internal-gap", str(DEFAULT["gap"]),
                     "--iterations", str(ITERATIONS),
                     "--min-feature-intensity", str(I),
                     "--force"],
                    check=True, stdout=subprocess.DEVNULL,
                )
                ms1, nbytes = measure(out)
                shutil.rmtree(out, ignore_errors=True)
                ms1_pct = 100 * (1 - ms1 / base_ms1)
                bytes_pct = 100 * (1 - nbytes / base_bytes)
                int_ms1[ii] = ms1_pct
                int_bytes[ii] = bytes_pct
                int_rows.append({
                    "min_feature_intensity": I,
                    "ms1_removed_pct": ms1_pct, "bytes_removed_pct": bytes_pct,
                })
                print(f"min_feature_intensity {I:5d}: MS1 -{ms1_pct:.1f}%  bytes -{bytes_pct:.1f}%")
        pd.DataFrame(rows).to_csv(csv_path, index=False)
        pd.DataFrame(int_rows).to_csv(int_csv_path, index=False)

    default_rc = (row_keys.index((DEFAULT["halo"], DEFAULT["gap"])),
                  LENGTHS.index(DEFAULT["length"]))

    fig, axd = plt.subplot_mosaic(
        [["ms1", "bytes"], ["intensity", "intensity"]],
        figsize=(14, 10),
    )
    heatmap(axd["ms1"], ms1_grid, row_labels, "MS1 peaks removed", "% MS1 removed", fig, default_rc)
    heatmap(axd["bytes"], bytes_grid, [""] * len(row_labels), "tdf_bin bytes removed",
            "% bytes removed", fig, default_rc)
    axd["ms1"].set_ylabel("halo state · max_internal_gap")
    line_panel(axd["intensity"], INTENSITIES, int_ms1, int_bytes, default_intensity=0)
    fig.suptitle("dnoise compression: filter parameter interactions (single 5-min run; "
                 "red = default)", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "compression_sweep.png", dpi=150, bbox_inches="tight")
    if (REPO / "paper" / "figures").is_dir():
        fig.savefig(REPO / "paper" / "figures" / "compression_sweep.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}/compression_sweep.png, compression_sweep.csv, "
          "and compression_intensity_sweep.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
