#!/usr/bin/env python3
"""Compare original vs denoised arms: identifications and LFQ quant.

Reads Sage outputs from results/{original,denoised}/ (results.sage.tsv, lfq.tsv),
computes ID counts at 1% FDR and LFQ accuracy/precision against the known
hybrid-proteome ratios, and writes results/analysis/{summary.csv, *.png}.

Ground truth log2(A/B): HUMAN 0, YEAST +1, ECOLI -2 (Condition A 65/30/5,
B 65/15/20 w/w human/yeast/ecoli).
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from _metrics import (
    EXPECTED,
    SPECIES,
    SPECIES_COLOR,
    id_metrics,
    lfq_metrics,
    lfq_table,
    pair_accuracy,
)

ROOT = Path(__file__).resolve().parents[1]
DATASET = os.environ.get("DATASET", "dda_5min")
RESULTS = ROOT / "results" / DATASET
OUT = RESULTS / "analysis"
# Arms compared, in order; only those with a Sage result on disk are included.
# "msms" is the MS1+MS/MS-denoised arm (present only after a --denoise-msms run).
CANDIDATE_ARMS = ["original", "denoised", "intensity", "msms"]
ARMS = [a for a in CANDIDATE_ARMS if (RESULTS / a / "results.sage.tsv").is_file()]
ARM_COLOR = {"original": "#1f77b4", "denoised": "#d62728",
             "intensity": "#9467bd", "msms": "#2ca02c"}


# ---------- plots ----------

def plot_ratio_clouds(prots: dict[str, pd.DataFrame | None]) -> None:
    fig, axes = plt.subplots(1, len(ARMS), figsize=(11, 5), sharey=True)
    for ax, arm in zip(axes, ARMS):
        prot = prots[arm]
        ax.set_title(arm)
        ax.axhline(0, color="gray", lw=0.5)
        if prot is not None and not prot.empty:
            for i, sp in enumerate(SPECIES):
                vals = prot[prot["species"] == sp]["log2_ratio"].dropna().values
                if len(vals):
                    x = np.random.normal(i, 0.07, len(vals))
                    ax.scatter(x, vals, s=4, alpha=0.25, color=SPECIES_COLOR[sp])
                ax.hlines(EXPECTED[sp], i - 0.35, i + 0.35, color=SPECIES_COLOR[sp], lw=2)
        ax.set_xticks(range(len(SPECIES)))
        ax.set_xticklabels(SPECIES)
        ax.set_ylim(-5, 4)
    axes[0].set_ylabel("log2(A/B) protein ratio")
    fig.suptitle("LFQ ratio accuracy (bars = expected)")
    fig.tight_layout()
    fig.savefig(OUT / "lfq_ratio_clouds.png", dpi=150)
    plt.close(fig)


def plot_ratio_violins(prots: dict[str, pd.DataFrame | None]) -> None:
    """Per-species log2(A/B) distributions, one violin per arm per species,
    with the expected ratio drawn as a dashed line per species."""
    n = len(ARMS)
    width = 0.8 / n
    # Center the arms' violins within each species slot.
    offsets = [(-(n - 1) / 2 + j) * width for j in range(n)]
    fig, ax = plt.subplots(figsize=(8, 5))

    for j, arm in enumerate(ARMS):
        prot = prots[arm]
        if prot is None or prot.empty:
            continue
        data, pos = [], []
        for i, sp in enumerate(SPECIES):
            vals = prot[prot["species"] == sp]["log2_ratio"].dropna().values
            if len(vals) >= 2:  # violinplot needs a distribution
                data.append(vals)
                pos.append(i + offsets[j])
        if not data:
            continue
        vp = ax.violinplot(data, positions=pos, widths=width * 0.9, showmedians=True, showextrema=False)
        for body in vp["bodies"]:
            body.set_facecolor(ARM_COLOR[arm])
            body.set_alpha(0.55)
            body.set_edgecolor("black")
            body.set_linewidth(0.5)
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(1.2)

    for i, sp in enumerate(SPECIES):
        ax.hlines(EXPECTED[sp], i - 0.4, i + 0.4, color="gray", lw=2, ls="--", zorder=5)

    ax.set_xticks(range(len(SPECIES)))
    ax.set_xticklabels(SPECIES)
    ax.set_ylabel("log2(A/B) protein ratio")
    ax.set_ylim(-5, 4)
    ax.set_title("LFQ ratio distributions (dashed = expected)")
    ax.legend(handles=[Patch(facecolor=ARM_COLOR[a], alpha=0.55, edgecolor="black", label=a) for a in ARMS])
    fig.tight_layout()
    fig.savefig(OUT / "lfq_ratio_violins.png", dpi=150)
    plt.close(fig)


def plot_cv(prots: dict[str, pd.DataFrame | None]) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    for arm in ARMS:
        prot = prots[arm]
        if prot is None or prot.empty:
            continue
        cv = pd.concat([prot["cv_A"], prot["cv_B"]]).dropna()
        cv = cv[cv < 1.0]
        ax.hist(cv, bins=40, histtype="step", lw=2, label=f"{arm} (median {cv.median():.3f})")
    ax.set_xlabel("protein CV within condition")
    ax.set_ylabel("count")
    ax.set_title("LFQ precision")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "lfq_cv.png", dpi=150)
    plt.close(fig)


def plot_ids(summary: pd.DataFrame) -> None:
    metrics = ["n_psm", "n_peptide", "n_protein", "n_quantified"]
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(metrics))
    n = len(ARMS)
    w = 0.8 / n
    for i, arm in enumerate(ARMS):
        vals = [summary.loc[arm, m] for m in metrics]
        bars = ax.bar(x + (-(n - 1) / 2 + i) * w, vals, w, label=arm, color=ARM_COLOR.get(arm))
        ax.bar_label(bars, fmt="%d", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.set_ylabel("count @ 1% FDR")
    ax.set_title("Identifications & quantified proteins")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "id_counts.png", dpi=150)
    plt.close(fig)


def plot_accuracy(acc: pd.DataFrame) -> None:
    """Observed vs expected median log2 ratio across all pairs/species/arms —
    shows accuracy across the full dynamic range (−2.7 to +3.3)."""
    fig, ax = plt.subplots(figsize=(6, 6))
    lim = (-3.5, 4.0)
    ax.plot(lim, lim, color="gray", ls="--", lw=1, zorder=0, label="ideal")
    marker = {"HUMAN": "o", "YEAST": "s", "ECOLI": "^"}
    for arm in ARMS:
        d = acc[acc["arm"] == arm]
        for sp in SPECIES:
            ds = d[d["species"] == sp]
            ax.scatter(ds["expected"], ds["observed"], color=ARM_COLOR.get(arm),
                       marker=marker[sp], s=45, alpha=0.8,
                       label=arm if sp == "HUMAN" else None)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("expected log2 ratio")
    ax.set_ylabel("observed median log2 ratio")
    ax.set_title("LFQ accuracy across the dynamic range\n(A/B, A/C, B/C; ○ human △ ecoli □ yeast)")
    ax.legend(title="arm")
    fig.tight_layout()
    fig.savefig(OUT / "lfq_accuracy.png", dpi=150)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, prots, acc = {}, {}, []
    for arm in ARMS:
        prots[arm] = lfq_table(RESULTS / arm)
        rows[arm] = {**id_metrics(RESULTS / arm), **lfq_metrics(prots[arm], RESULTS / arm)}
        acc += pair_accuracy(arm, RESULTS / arm)

    summary = pd.DataFrame(rows).T
    # Order columns: counts first, then quant.
    summary = summary.reindex(sorted(summary.columns), axis=1)
    summary.to_csv(OUT / "summary.csv")

    acc_df = pd.DataFrame(acc)
    acc_df.to_csv(OUT / "accuracy.csv", index=False)

    plot_ratio_clouds(prots)
    plot_ratio_violins(prots)
    plot_cv(prots)
    plot_ids(summary)
    if not acc_df.empty:
        plot_accuracy(acc_df)

    pd.set_option("display.width", 200, "display.max_columns", 100)
    print(summary.T)
    if not acc_df.empty:
        print("\nLFQ accuracy across pairs (observed vs expected median log2):")
        print(acc_df.pivot_table(index=["pair", "species"], columns="arm",
                                 values="observed").round(2))
    print(f"\nWrote {OUT}/summary.csv, accuracy.csv, and plots.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
