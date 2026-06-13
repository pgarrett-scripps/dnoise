#!/usr/bin/env python3
"""Merged compression figure (paper Fig 2) across gradients and arms.

Reads results/<dataset>/analysis/data_reduction.csv (written by
07_data_reduction.py) for each dataset and emits a single compression-only
figure to results/compare/:
  Panel A: on-disk frame binary (GB) per gradient, three arms (raw / MS1 /
           MS1+MS/MS), with % of raw labelled on the denoised bars.
  Panel B: peaks by MS level (billions), stacked MS1 + MS/MS, per gradient and
           arm, with total + % removed labelled.

Compression statistics only -- identification/LFQ metrics live in Fig 5
(15_id_figure.py / id_counts) and Table 2.

Usage: uv run scripts/14_compression_figure.py [dataset ...]
       (default: dda_5min dda_15min)
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "compare"
ARMS = ["raw", "MS1", "MS1+MS/MS"]
ARM_COLOR = {"raw": "#1f77b4", "MS1": "#d62728", "MS1+MS/MS": "#2ca02c"}
MS1_COLOR, MS2_COLOR = "#4c72b0", "#dd8452"
GRAD_LABEL = {"dda_5min": "5 min", "dda_15min": "15 min"}


def load(dataset: str) -> dict | None:
    p = ROOT / "results" / dataset / "analysis" / "data_reduction.csv"
    if not p.is_file():
        print(f"  skip {dataset}: no {p}")
        return None
    df = pd.read_csv(p)
    agg = {}
    for arm in ARMS:
        agg[arm] = {
            "bytes": df[f"bytes_{arm}"].sum(),
            "ms1": df[f"ms1_{arm}"].sum(),
            "ms2": df[f"ms2_{arm}"].sum(),
        }
    return {"dataset": dataset, "agg": agg, "nruns": len(df)}


def main() -> int:
    datasets = sys.argv[1:] or ["dda_5min", "dda_15min"]
    OUT.mkdir(parents=True, exist_ok=True)
    recs = [r for r in (load(d) for d in datasets) if r is not None]
    if not recs:
        print("no datasets with data_reduction.csv found")
        return 1

    labels = [GRAD_LABEL.get(r["dataset"], r["dataset"]) for r in recs]
    ng = len(recs)
    x = np.arange(ng)
    w = 0.8 / len(ARMS)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A: frame-binary size (GB), grouped by gradient, one bar per arm.
    for j, arm in enumerate(ARMS):
        gb = [r["agg"][arm]["bytes"] / 1e9 for r in recs]
        off = (-(len(ARMS) - 1) / 2 + j) * w
        bars = axA.bar(x + off, gb, w, label=arm, color=ARM_COLOR[arm])
        for k, r in enumerate(recs):
            raw_b = r["agg"]["raw"]["bytes"]
            pct = 100 * (1 - r["agg"][arm]["bytes"] / raw_b)
            txt = f"{gb[k]:.1f}" + (f"\n(−{pct:.0f}%)" if arm != "raw" else "")
            axA.text(x[k] + off, gb[k], txt, ha="center", va="bottom", fontsize=7)
    axA.set_xticks(x)
    axA.set_xticklabels(labels)
    axA.set_ylabel("analysis.tdf_bin (GB, summed over runs)")
    axA.set_title("On-disk frame binary")
    axA.legend()
    axA.set_ylim(0, max(r["agg"]["raw"]["bytes"] for r in recs) / 1e9 * 1.18)

    # Panel B: peaks (billions), stacked MS1 + MS/MS, grouped by gradient/arm.
    for j, arm in enumerate(ARMS):
        ms1 = [r["agg"][arm]["ms1"] / 1e9 for r in recs]
        ms2 = [r["agg"][arm]["ms2"] / 1e9 for r in recs]
        off = (-(len(ARMS) - 1) / 2 + j) * w
        axB.bar(x + off, ms1, w, color=MS1_COLOR)
        axB.bar(x + off, ms2, w, bottom=ms1, color=MS2_COLOR, hatch="//", edgecolor="white", linewidth=0)
        for k, r in enumerate(recs):
            tot = ms1[k] + ms2[k]
            raw_tot = (r["agg"]["raw"]["ms1"] + r["agg"]["raw"]["ms2"]) / 1e9
            pct = 100 * (1 - tot / raw_tot)
            lbl = f"{tot:.1f}" + (f"\n−{pct:.0f}%" if arm != "raw" else "")
            axB.text(x[k] + off, tot, lbl, ha="center", va="bottom", fontsize=6.5)
    axB.set_xticks(x)
    axB.set_xticklabels(labels)
    axB.set_ylabel("peaks (billions, summed over runs)")
    axB.set_title("Peaks by level  (bars within each gradient: raw, MS1, MS1+MS/MS)")
    axB.set_ylim(0, max((r["agg"]["raw"]["ms1"] + r["agg"]["raw"]["ms2"]) for r in recs) / 1e9 * 1.18)
    axB.legend(handles=[
        Patch(facecolor=MS1_COLOR, label="MS1 peaks"),
        Patch(facecolor=MS2_COLOR, hatch="//", label="MS/MS peaks"),
    ])

    fig.suptitle(f"dnoise data-volume reduction across gradients "
                 f"({' / '.join(f'{r['nruns']} runs at {GRAD_LABEL.get(r['dataset'], r['dataset'])}' for r in recs)})")
    fig.tight_layout()
    fig.savefig(OUT / "compression.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT}/compression.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
