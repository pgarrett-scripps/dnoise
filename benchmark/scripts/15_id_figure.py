#!/usr/bin/env python3
"""Identifications figure (paper Fig 5) across gradients and arms.

Reads results/<dataset>/analysis/summary.csv (written by 06_analyze.py) for each
dataset and emits one figure to results/compare/. One panel *per metric* (PSMs,
peptides, protein groups, quantified), each with its own y-axis so a high-count
metric (PSMs) does not flatten the low-count ones (protein groups, quantified).
Within each panel the x-axis is the gradient and bars are the three arms.

Usage: uv run scripts/15_id_figure.py [dataset ...]   (default: dda_5min dda_15min)
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "compare"
# summary.csv arm index -> display label / color
ARMS = ["original", "denoised", "msms"]
ARM_LABEL = {"original": "original", "denoised": "MS1", "msms": "MS1+MS/MS"}
ARM_COLOR = {"original": "#1f77b4", "denoised": "#d62728", "msms": "#2ca02c"}
METRICS = ["n_psm", "n_peptide", "n_protein", "n_quantified"]
METRIC_LABEL = ["PSMs", "peptides", "protein groups", "quantified"]
GRAD_LABEL = {"dda_5min": "5 min", "dda_15min": "15 min"}


def load(dataset: str) -> pd.DataFrame | None:
    p = ROOT / "results" / dataset / "analysis" / "summary.csv"
    if not p.is_file():
        print(f"  skip {dataset}: no {p}")
        return None
    return pd.read_csv(p, index_col=0)


def main() -> int:
    datasets = sys.argv[1:] or ["dda_5min", "dda_15min"]
    OUT.mkdir(parents=True, exist_ok=True)
    recs = [(d, s) for d in datasets if (s := load(d)) is not None]
    if not recs:
        print("no datasets with summary.csv found")
        return 1

    arms_present = [a for a in ARMS if any(a in s.index for _, s in recs)]
    labels = [GRAD_LABEL.get(d, d) for d, _ in recs]
    x = np.arange(len(recs))            # one group per gradient
    w = 0.8 / len(arms_present)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, metric, mlabel in zip(axes.ravel(), METRICS, METRIC_LABEL):
        for j, arm in enumerate(arms_present):
            vals = [float(s.loc[arm, metric]) if (arm in s.index and metric in s.columns) else 0
                    for _, s in recs]
            off = (-(len(arms_present) - 1) / 2 + j) * w
            bars = ax.bar(x + off, vals, w, label=ARM_LABEL[arm], color=ARM_COLOR[arm])
            ax.bar_label(bars, fmt="%d", fontsize=7, padding=2)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_title(mlabel)
        ax.set_ylabel("count @ 1% FDR")
        ax.set_ylim(0, max(ax.get_ylim()[1], 1) * 1.15)
        ax.margins(x=0.15)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Identifications and quantified proteins across gradients "
                 "(per-metric scale)", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "id_counts.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT}/id_counts.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
