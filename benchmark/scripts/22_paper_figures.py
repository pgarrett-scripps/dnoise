#!/usr/bin/env python3
"""Draft main-text figures from the computed dda_5min arms:

  Fig 5 (centroiders): original vs watershed vs box   -> paper/figures/fig5_centroiders.png
  Fig 6 (Bruker):      original vs dnoise vs bruker(default, matched)
                                                       -> paper/figures/fig6_bruker.png

Each: 3 panels — data reduction (MS1 kept % + binary %), quantified proteins
(1% FDR), and LFQ accuracy (observed vs expected median log2 across all pairs/
species). Uses the shared metric helpers; compression is summed over the 18 runs.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _metrics import SPECIES, id_metrics, lfq_metrics, lfq_table, pair_accuracy, stats

ARM3_COLOR = {"original": "#1f77b4", "denoised": "#d62728", "msms": "#2ca02c"}
ARM3_LABEL = {"original": "original", "denoised": "MS1", "msms": "MS1+MS/MS"}
# arm -> data .d subdir (results subdir is the arm name itself)
ARM3_DATA = {"original": "raw", "denoised": "denoised", "msms": "denoised_msms"}

ROOT = Path(__file__).resolve().parents[1]
FIGS = ROOT.parent / "paper" / "figures"
DDA = ROOT / "results" / "dda_5min"
DATA = ROOT / "data" / "dda_5min"

# label -> (data .d subdir, results subdir, color)
ARM = {
    "original":       ("raw",            "original",       "#1f77b4"),
    "dnoise":         ("denoised",       "denoised",       "#d62728"),
    "watershed":      ("watershed",      "watershed",      "#2ca02c"),
    "box":            ("box",            "box",            "#ff7f0e"),
    "bruker (default)": ("bruker_default", "bruker_default", "#9467bd"),
    "bruker (matched)": ("bruker",        "bruker",         "#8c564b"),
}
MARK = {"HUMAN": "o", "YEAST": "s", "ECOLI": "^"}


def compression(dsub: str) -> tuple[int, int]:
    ms1 = b = 0
    for d in sorted((DATA / dsub).glob("*.d")):
        a, _, bb = stats(d)
        ms1 += a
        b += bb
    return ms1, b


def collect(labels: list[str]) -> dict:
    o_ms1, o_bin = compression(ARM["original"][0])
    rows = {}
    for lab in labels:
        dsub, rsub, color = ARM[lab]
        ms1, b = compression(dsub)
        rdir = DDA / rsub
        lm = lfq_metrics(lfq_table(rdir))
        rows[lab] = {
            "ms1_pct": 100 * ms1 / o_ms1,
            "bin_pct": 100 * b / o_bin,
            "n_quant": lm.get("n_quantified", 0),
            "acc": pair_accuracy(lab, rdir),
            "color": color,
        }
    return rows


def make_fig(labels: list[str], title: str, out: Path) -> None:
    rows = collect(labels)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))

    # Panel 1: data reduction (MS1 kept %, binary %)
    x = np.arange(len(labels))
    w = 0.38
    ax[0].bar(x - w / 2, [rows[l]["ms1_pct"] for l in labels], w, label="MS1 peaks", color="#4c72b0")
    ax[0].bar(x + w / 2, [rows[l]["bin_pct"] for l in labels], w, label="binary size", color="#dd8452")
    ax[0].axhline(100, color="gray", lw=0.8, ls="--")
    ax[0].set_ylabel("% of original (18 runs)")
    ax[0].set_title("Data reduction")
    ax[0].set_xticks(x); ax[0].set_xticklabels(labels, rotation=25, ha="right")
    ax[0].legend(fontsize=8)

    # Panel 2: quantified proteins
    bars = ax[1].bar(x, [rows[l]["n_quant"] for l in labels],
                     color=[rows[l]["color"] for l in labels])
    ax[1].bar_label(bars, fmt="%d", fontsize=8)
    ax[1].set_ylabel("proteins quantified (1% FDR)")
    ax[1].set_title("Quantification coverage")
    ax[1].set_xticks(x); ax[1].set_xticklabels(labels, rotation=25, ha="right")
    lo = min(rows[l]["n_quant"] for l in labels)
    ax[1].set_ylim(lo * 0.95, max(rows[l]["n_quant"] for l in labels) * 1.02)

    # Panel 3: LFQ accuracy (observed vs expected, all pairs/species)
    lim = (-3.5, 4.0)
    ax[2].plot(lim, lim, color="gray", ls="--", lw=1, zorder=0)
    for lab in labels:
        for a in rows[lab]["acc"]:
            ax[2].scatter(a["expected"], a["observed"], color=rows[lab]["color"],
                          marker=MARK[a["species"]], s=40, alpha=0.8)
        ax[2].scatter([], [], color=rows[lab]["color"], label=lab)  # legend proxy
    ax[2].set_xlim(*lim); ax[2].set_ylim(*lim)
    ax[2].set_xlabel("expected log2 ratio"); ax[2].set_ylabel("observed median log2")
    ax[2].set_title("LFQ accuracy (○ human △ ecoli □ yeast)")
    ax[2].legend(fontsize=7)

    fig.suptitle(title, fontsize=13)
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def compress(dataset: str, dsub: str) -> tuple[int, int]:
    ms1 = b = 0
    for d in sorted((ROOT / "data" / dataset / dsub).glob("*.d")):
        a, _, bb = stats(d)
        ms1 += a
        b += bb
    return ms1, b


def _acc_scatter(ax, series, title="LFQ accuracy (○ human △ ecoli □ yeast)"):
    """series: list of (color, [accuracy-dicts])."""
    lim = (-3.5, 4.0)
    ax.plot(lim, lim, color="gray", ls="--", lw=1, zorder=0)
    for color, acc in series:
        for a in acc:
            ax.scatter(a["expected"], a["observed"], color=color,
                       marker=MARK[a["species"]], s=35, alpha=0.75)
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_xlabel("expected log2 ratio"); ax.set_ylabel("observed median log2")
    ax.set_title(title)


def _arm3_quad(binpct, quant, acc, arms, grads, title, out) -> None:
    """2x2 figure for a 3-arm, 2-gradient comparison:
      top:    data reduction (binary % of raw) + quantification coverage
      bottom: LFQ accuracy, a separate panel per gradient.
    binpct/quant: {grad: [per-arm]}; acc: {grad: [(arm_color, accuracy-dicts)]}."""
    x = np.arange(len(arms)); w = 0.38
    GC = {"5min": "#4c72b0", "15min": "#dd8452"}
    GL = {"5min": "5 min", "15min": "15 min"}
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    # [0,0] data reduction (binary % of raw), grouped by gradient
    a = ax[0, 0]
    for i, g in enumerate(grads):
        bars = a.bar(x + (i - 0.5) * w, binpct[g], w, label=GL[g], color=GC[g])
        a.bar_label(bars, fmt="%.0f", fontsize=7, padding=1)
    a.axhline(100, color="gray", lw=0.8, ls="--")
    a.set_ylabel("frame binary, % of raw"); a.set_title("Data reduction")
    a.set_xticks(x); a.set_xticklabels([ARM3_LABEL[m] for m in arms]); a.legend(fontsize=8)
    # [0,1] quantified proteins, grouped by gradient
    a = ax[0, 1]
    for i, g in enumerate(grads):
        bars = a.bar(x + (i - 0.5) * w, quant[g], w, label=GL[g], color=GC[g])
        a.bar_label(bars, fmt="%d", fontsize=7, padding=1)
    a.set_ylabel("proteins quantified (1% FDR)"); a.set_title("Quantification coverage")
    a.set_xticks(x); a.set_xticklabels([ARM3_LABEL[m] for m in arms]); a.legend(fontsize=8, loc="upper left")
    # [1,0],[1,1] LFQ accuracy, one panel per gradient (arms colored)
    for col, g in enumerate(grads):
        _acc_scatter(ax[1, col], acc[g], title=f"LFQ accuracy, {GL[g]} (○ human △ ecoli □ yeast)")
    for m in arms:
        ax[1, 0].scatter([], [], color=ARM3_COLOR[m], label=ARM3_LABEL[m])
    ax[1, 0].legend(fontsize=7)
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(); fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig); print(f"wrote {out}")


def make_dda() -> None:
    """Fig 3 — DDA, both gradients, three arms (Sage metrics)."""
    arms = ["original", "denoised", "msms"]
    grads = ["5min", "15min"]
    binpct = {g: [] for g in grads}; quant = {g: [] for g in grads}
    acc = {g: [] for g in grads}
    for g in grads:
        _, o_bin = compress(f"dda_{g}", "raw")
        for arm in arms:
            _, b = compress(f"dda_{g}", ARM3_DATA[arm])
            binpct[g].append(100 * b / o_bin)
            rdir = ROOT / "results" / f"dda_{g}" / arm
            quant[g].append(lfq_metrics(lfq_table(rdir)).get("n_quantified", 0))
            acc[g].append((ARM3_COLOR[arm], pair_accuracy(arm, rdir)))
    _arm3_quad(binpct, quant, acc, arms, grads,
               "DDA (ddaPASEF): original vs MS1 vs MS1+MS/MS denoising, 5 & 15 min",
               FIGS / "fig3_dda.png")


def make_dia() -> None:
    """Fig 4 — DIA (diaPASEF), both gradients, three arms (DIA-NN metrics)."""
    arms = ["original", "denoised", "msms"]
    grads = ["5min", "15min"]
    binpct = {g: [] for g in grads}; quant = {g: [] for g in grads}
    acc = {g: [] for g in grads}
    for g in grads:
        summ = pd.read_csv(ROOT / f"results/dia_{g}/analysis/summary.csv", index_col=0)
        accdf = pd.read_csv(ROOT / f"results/dia_{g}/analysis/accuracy.csv")
        _, o_bin = compress(f"dia_{g}", "raw")
        for arm in arms:
            _, b = compress(f"dia_{g}", ARM3_DATA[arm])
            binpct[g].append(100 * b / o_bin)
            quant[g].append(int(summ.loc[arm, "n_quantified"]) if arm in summ.index else 0)
            d = accdf[accdf["arm"] == arm]
            acc[g].append((ARM3_COLOR[arm],
                           [{"expected": r.expected, "observed": r.observed, "species": r.species}
                            for r in d.itertuples()]))
    _arm3_quad(binpct, quant, acc, arms, grads,
               "DIA (diaPASEF): original vs MS1 vs MS1+MS/MS denoising, 5 & 15 min",
               FIGS / "fig4_dia.png")


def main() -> int:
    make_dda()
    make_dia()
    make_fig(["original", "watershed", "box"],
             "Optional centroiding (MS1-only; identifications unchanged)",
             FIGS / "fig5_centroiders.png")
    make_fig(["original", "dnoise", "bruker (default)", "bruker (matched)"],
             "dnoise vs Bruker Minesweeper (MS1-only; identifications unchanged)",
             FIGS / "fig6_bruker.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
