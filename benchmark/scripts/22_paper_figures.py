#!/usr/bin/env python3
"""Draft main-text figures from the computed dda_5min / dia arms:

  Fig 3 (DDA):         original vs MS1 vs MS1+MS/MS    -> paper/figures/fig3_dda.png
  Fig 4 (DIA):         original vs MS1 vs MS1+MS/MS    -> paper/figures/fig4_dia.png
  Fig 5 (centroiders): original vs watershed vs box    -> paper/figures/fig5_centroiders.png

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

# Figures are placed at the full ~6.3 in text width; they are authored larger and
# downscaled, so fonts must be large enough that nothing drops below 4.5 pt at
# print. With the smallest text at 11 pt and the most-downscaled figure (~0.42x),
# 11 x 0.42 = 4.65 pt > 4.5 pt.
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "figure.titlesize": 15,
})

# Colorblind-safe qualitative palette (Wong, Nat. Methods 2011); avoids the
# red/green pairing. Species are additionally distinguished by marker shape.
ARM3_COLOR = {"original": "#0072B2", "denoised": "#E69F00", "msms": "#009E73"}
ARM3_LABEL = {"original": "original", "denoised": "MS1", "msms": "MS1+MS/MS"}
# arm -> data .d subdir (results subdir is the arm name itself)
ARM3_DATA = {"original": "raw", "denoised": "denoised", "msms": "denoised_msms"}

ROOT = Path(__file__).resolve().parents[1]
FIGS = ROOT.parent / "paper" / "figures"
DDA = ROOT / "results" / "dda_5min"
DATA = ROOT / "data" / "dda_5min"

# label -> (data .d subdir, results subdir, color)
ARM = {
    "original":       ("raw",            "original",       "#0072B2"),
    "dnoise":         ("denoised",       "denoised",       "#E69F00"),
    "watershed":      ("watershed",      "watershed",      "#009E73"),
    "box":            ("box",            "box",            "#D55E00"),
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
    fig, ax = plt.subplots(1, 3, figsize=(12, 4.2))

    # Panel 1: data reduction (MS1 kept %, binary %)
    x = np.arange(len(labels))
    w = 0.38
    ax[0].bar(x - w / 2, [rows[l]["ms1_pct"] for l in labels], w, label="MS1 peaks", color="#0072B2")
    ax[0].bar(x + w / 2, [rows[l]["bin_pct"] for l in labels], w, label="binary size", color="#E69F00")
    ax[0].axhline(100, color="gray", lw=0.8, ls="--")
    ax[0].set_ylabel("% of original (18 runs)")
    ax[0].set_title("Data reduction")
    ax[0].set_xticks(x); ax[0].set_xticklabels(labels, rotation=25, ha="right")
    ax[0].legend(fontsize=11)

    # Panel 2: quantified proteins
    bars = ax[1].bar(x, [rows[l]["n_quant"] for l in labels],
                     color=[rows[l]["color"] for l in labels])
    ax[1].bar_label(bars, fmt="%d", fontsize=11)
    ax[1].set_ylabel("proteins quantified (1% FDR)")
    ax[1].set_title("Quantification coverage")
    ax[1].set_xticks(x); ax[1].set_xticklabels(labels, rotation=25, ha="right")
    lo = min(rows[l]["n_quant"] for l in labels)
    ax[1].set_ylim(lo * 0.95, max(rows[l]["n_quant"] for l in labels) * 1.02)

    # Panel 3: LFQ accuracy, residual (observed - expected) vs expected
    xlim = (-3.5, 4.0)
    ylo, yhi = _resid_ylim([[(rows[l]["color"], rows[l]["acc"]) for l in labels]])
    ax[2].axhline(0, color="gray", ls="--", lw=1, zorder=0)
    for lab in labels:
        for a in rows[lab]["acc"]:
            ax[2].scatter(a["expected"], a["observed"] - a["expected"], color=rows[lab]["color"],
                          marker=MARK[a["species"]], s=60, alpha=0.8)
        ax[2].scatter([], [], color=rows[lab]["color"], label=lab)  # legend proxy
    ax[2].set_xlim(*xlim); ax[2].set_ylim(ylo, yhi)
    ax[2].set_xlabel("expected log2 ratio"); ax[2].set_ylabel("observed − expected (log2)")
    ax[2].set_title("LFQ accuracy (○ human △ ecoli □ yeast)")
    ax[2].legend(fontsize=11)

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


def _resid_ylim(series_list, pad_frac: float = 0.15, min_pad: float = 0.3):
    """series_list: list of `series` (each a list of (color, [accuracy-dicts])).
    Shared y-limits for a residual (observed - expected) scatter, computed
    across ALL series passed in so multiple panels (e.g. one per gradient)
    share one y-scale for a fair visual comparison."""
    resid = [a["observed"] - a["expected"]
             for series in series_list for _, acc in series for a in acc]
    if not resid:
        return (-1.0, 1.0)
    pad = max(min_pad, pad_frac * (max(resid) - min(resid)))
    return (min(resid) - pad, max(resid) + pad)


def _acc_scatter(ax, series, title="LFQ accuracy (○ human △ ecoli □ yeast)", ylim=None):
    """series: list of (color, [accuracy-dicts]). Plots the RESIDUAL
    (observed − expected) against expected, so the ideal line is horizontal
    at 0 instead of a 45° diagonal — arms that overlap tightly along the
    diagonal in an observed-vs-expected view separate out much more clearly
    once the shared slope is subtracted off. Pass `ylim` to share a y-scale
    across multiple panels; otherwise it is computed from this series alone."""
    xlim = (-3.5, 4.0)
    ylo, yhi = ylim if ylim is not None else _resid_ylim([series])
    ax.axhline(0, color="gray", ls="--", lw=1, zorder=0)
    for color, acc in series:
        for a in acc:
            ax.scatter(a["expected"], a["observed"] - a["expected"], color=color,
                       marker=MARK[a["species"]], s=60, alpha=0.75)
    ax.set_xlim(*xlim); ax.set_ylim(ylo, yhi)
    ax.set_xlabel("expected log2 ratio"); ax.set_ylabel("observed − expected (log2)")
    ax.set_title(title)


def _arm3_quad(binpct, quant, acc, arms, grads, title, out) -> None:
    """2x2 figure for a 3-arm, 2-gradient comparison:
      top:    data reduction (binary % of raw) + quantification coverage
      bottom: LFQ accuracy, a separate panel per gradient.
    binpct/quant: {grad: [per-arm]}; acc: {grad: [(arm_color, accuracy-dicts)]}."""
    x = np.arange(len(arms)); w = 0.38
    GC = {"5min": "#0072B2", "15min": "#E69F00"}
    GL = {"5min": "5 min", "15min": "15 min"}
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    # [0,0] data reduction (binary % of raw), grouped by gradient
    a = ax[0, 0]
    for i, g in enumerate(grads):
        bars = a.bar(x + (i - 0.5) * w, binpct[g], w, label=GL[g], color=GC[g])
        a.bar_label(bars, fmt="%.0f", fontsize=11, padding=1)
    a.axhline(100, color="gray", lw=0.8, ls="--")
    a.set_ylabel("frame binary, % of raw"); a.set_title("Data reduction")
    a.set_xticks(x); a.set_xticklabels([ARM3_LABEL[m] for m in arms]); a.legend(fontsize=11)
    # [0,1] quantified proteins, grouped by gradient
    a = ax[0, 1]
    for i, g in enumerate(grads):
        bars = a.bar(x + (i - 0.5) * w, quant[g], w, label=GL[g], color=GC[g])
        a.bar_label(bars, fmt="%d", fontsize=11, padding=1)
    a.set_ylabel("proteins quantified (1% FDR)"); a.set_title("Quantification coverage")
    a.set_xticks(x); a.set_xticklabels([ARM3_LABEL[m] for m in arms]); a.legend(fontsize=11, loc="upper left")
    # [1,0],[1,1] LFQ accuracy, one panel per gradient (arms colored), shared y-scale
    acc_ylim = _resid_ylim([acc[g] for g in grads])
    for col, g in enumerate(grads):
        _acc_scatter(ax[1, col], acc[g], title=f"LFQ accuracy, {GL[g]} (○ human △ ecoli □ yeast)",
                     ylim=acc_ylim)
    for m in arms:
        ax[1, 0].scatter([], [], color=ARM3_COLOR[m], label=ARM3_LABEL[m])
    ax[1, 0].legend(fontsize=11)
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
