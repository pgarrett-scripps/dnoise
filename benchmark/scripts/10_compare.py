#!/usr/bin/env python3
"""Compare dnoise across datasets (e.g. dda_5min vs dda_15min).

Reads the per-dataset analysis already written by 06_analyze.py and
07_data_reduction.py:
  results/<dataset>/analysis/summary.csv         (per-arm IDs + LFQ metrics)
  results/<dataset>/analysis/data_reduction.csv  (per-file peaks + binary bytes)
  results/<dataset>/analysis/accuracy.csv         (observed vs expected log2)

and writes a side-by-side comparison to results/compare/:
  gradient_compare.png  -- data reduction, IDs, LFQ precision, LFQ accuracy
  gradient_compare.csv  -- the tidy table behind the figure

Datasets to compare come from argv (default: dda_5min dda_15min). A dataset
whose analysis is missing is skipped with a warning, so this is safe to run
before the second arm has finished.
"""

from __future__ import annotations

import sys
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
    "legend.fontsize": 11,
    "figure.titlesize": 15,
})
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "compare"
ARMS = ["original", "denoised", "msms"]
ARM_LABEL = {"original": "raw", "denoised": "MS1", "msms": "MS1+MS/MS"}
# Colorblind-safe (Wong, Nat. Methods 2011); avoids the red/green pairing.
ARM_COLOR = {"original": "#0072B2", "denoised": "#E69F00", "msms": "#009E73"}


def load(dataset: str) -> dict | None:
    adir = ROOT / "results" / dataset / "analysis"
    summ_p, dr_p = adir / "summary.csv", adir / "data_reduction.csv"
    if not summ_p.is_file():
        print(f"  skip {dataset}: no {summ_p}")
        return None
    summary = pd.read_csv(summ_p, index_col=0)
    rec: dict = {"dataset": dataset, "summary": summary}
    if dr_p.is_file():
        dr = pd.read_csv(dr_p)
        agg = {}
        for arm in ARMS:
            label = ARM_LABEL[arm]
            bcol = f"bytes_{label}" if f"bytes_{label}" in dr.columns else None
            m1, m2 = f"ms1_{label}", f"ms2_{label}"
            if bcol:
                agg[arm] = {
                    "bytes": dr[bcol].sum(),
                    "peaks": (dr[m1].sum() + dr[m2].sum()) if m1 in dr else np.nan,
                }
        rec["reduction"] = agg
    acc_p = adir / "accuracy.csv"
    if acc_p.is_file():
        rec["accuracy"] = pd.read_csv(acc_p)
    return rec


def main() -> int:
    datasets = sys.argv[1:] or ["dda_5min", "dda_15min"]
    OUT.mkdir(parents=True, exist_ok=True)
    recs = [r for r in (load(d) for d in datasets) if r is not None]
    if len(recs) < 1:
        print("no datasets with analysis found")
        return 1
    present_arms = [a for a in ARMS
                    if any(a in r["summary"].index for r in recs)]

    # ---- tidy table ----
    rows = []
    for r in recs:
        ds, summ = r["dataset"], r["summary"]
        for arm in present_arms:
            if arm not in summ.index:
                continue
            row = {"dataset": ds, "arm": ARM_LABEL[arm]}
            for k in ("n_psm", "n_peptide", "n_protein", "n_quantified", "median_cv"):
                if k in summ.columns:
                    row[k] = summ.loc[arm, k]
            red = r.get("reduction", {}).get(arm)
            if red:
                raw = r["reduction"].get("original", {})
                if raw.get("bytes"):
                    row["binary_pct_of_raw"] = 100 * red["bytes"] / raw["bytes"]
                if raw.get("peaks") and not pd.isna(red.get("peaks", np.nan)):
                    row["peaks_pct_of_raw"] = 100 * red["peaks"] / raw["peaks"]
            rows.append(row)
    tidy = pd.DataFrame(rows)
    tidy.to_csv(OUT / "gradient_compare.csv", index=False)

    # ---- figure ----
    nds = len(recs)
    x = np.arange(nds)
    w = 0.8 / max(len(present_arms), 1)
    labels = [r["dataset"] for r in recs]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    def grouped(ax, valfn, title, ylabel, fmt="%.0f"):
        for j, arm in enumerate(present_arms):
            vals = [valfn(r, arm) for r in recs]
            off = (-(len(present_arms) - 1) / 2 + j) * w
            bars = ax.bar(x + off, vals, w, label=ARM_LABEL[arm], color=ARM_COLOR[arm])
            ax.bar_label(bars, fmt=fmt, fontsize=11)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=11)

    def sval(r, arm, col):
        s = r["summary"]
        return float(s.loc[arm, col]) if arm in s.index and col in s.columns else 0.0

    def binary_pct(r, arm):
        red = r.get("reduction", {})
        raw = red.get("original", {}).get("bytes")
        cur = red.get(arm, {}).get("bytes")
        return 100 * cur / raw if raw and cur else 0.0

    grouped(axes[0, 0], binary_pct, "Frame-binary size (% of raw)", "% of raw", "%.0f")
    grouped(axes[0, 1], lambda r, a: sval(r, a, "n_protein"),
            "Protein groups @ 1% FDR", "count")
    grouped(axes[1, 0], lambda r, a: sval(r, a, "n_quantified"),
            "Quantified proteins (LFQ)", "count")
    grouped(axes[1, 1], lambda r, a: sval(r, a, "median_cv"),
            "LFQ precision (median CV)", "median CV", "%.3f")

    fig.suptitle("dnoise across gradients: " + " vs ".join(labels), fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT / "gradient_compare.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    pd.set_option("display.width", 200, "display.max_columns", 50)
    print(tidy.to_string(index=False))
    print(f"\nWrote {OUT}/gradient_compare.png and gradient_compare.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
