#!/usr/bin/env python3
"""LFQ accuracy figure (paper Fig 5) across gradients and arms.

Reads results/<dataset>/analysis/accuracy.csv (written by 06_analyze.py) for each
dataset and emits one figure to results/compare/ with one panel per gradient:
observed median log2 ratio vs. expected, for all three condition pairs
(A/B, A/C, B/C) and species, for each arm. Shows accuracy across the full dynamic
range and that the arms overlay (denoising does not distort ratios) at both
gradients.

Usage: uv run scripts/16_lfq_accuracy_figure.py [dataset ...]
       (default: dda_5min dda_15min)
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "compare"
ARMS = ["original", "denoised", "msms"]
ARM_LABEL = {"original": "original", "denoised": "MS1", "msms": "MS1+MS/MS"}
ARM_COLOR = {"original": "#1f77b4", "denoised": "#d62728", "msms": "#2ca02c"}
SPECIES = ["HUMAN", "YEAST", "ECOLI"]
MARKER = {"HUMAN": "o", "YEAST": "s", "ECOLI": "^"}
GRAD_LABEL = {"dda_5min": "5 min", "dda_15min": "15 min"}
LIM = (-3.5, 4.0)


def load(dataset: str) -> pd.DataFrame | None:
    p = ROOT / "results" / dataset / "analysis" / "accuracy.csv"
    if not p.is_file():
        print(f"  skip {dataset}: no {p}")
        return None
    return pd.read_csv(p)


def main() -> int:
    datasets = sys.argv[1:] or ["dda_5min", "dda_15min"]
    OUT.mkdir(parents=True, exist_ok=True)
    recs = [(d, a) for d in datasets if (a := load(d)) is not None]
    if not recs:
        print("no datasets with accuracy.csv found")
        return 1

    fig, axes = plt.subplots(1, len(recs), figsize=(5.6 * len(recs), 5.6),
                             squeeze=False, sharex=True, sharey=True)
    arms_present = [a for a in ARMS if any(a in acc["arm"].values for _, acc in recs)]

    for ax, (dataset, acc) in zip(axes[0], recs):
        ax.plot(LIM, LIM, color="gray", ls="--", lw=1, zorder=0, label="ideal")
        for arm in arms_present:
            d = acc[acc["arm"] == arm]
            for sp in SPECIES:
                ds = d[d["species"] == sp]
                ax.scatter(ds["expected"], ds["observed"], color=ARM_COLOR[arm],
                           marker=MARKER[sp], s=50, alpha=0.8,
                           label=ARM_LABEL[arm] if sp == "HUMAN" else None)
        ax.set_xlim(*LIM)
        ax.set_ylim(*LIM)
        ax.set_xlabel("expected log2 ratio")
        ax.set_title(GRAD_LABEL.get(dataset, dataset))
    axes[0, 0].set_ylabel("observed median log2 ratio")
    axes[0, 0].legend(title="arm", loc="upper left", fontsize=8)

    fig.suptitle("LFQ accuracy across the dynamic range "
                 "(A/B, A/C, B/C;  ○ human  △ ecoli  □ yeast)", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "lfq_accuracy.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT}/lfq_accuracy.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
