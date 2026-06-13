#!/usr/bin/env python3
"""Analyze the single-sample denoiser comparison (scripts/18_single_compare.sh).

All arms are the SAME acquisition denoised different ways, searched in one Sage
LFQ run, so lfq.tsv has one intensity column per method for the same peptides.
We compare each method's per-peptide MS1 intensity against `original` WITHOUT
normalization (any difference is the denoising effect, not loading).

Outputs (results/single_compare/):
  summary.csv             -- per arm: IDs, peptides quantified, median log2 vs orig
  retention.csv           -- per arm: kept / lost / gained peptides vs original
  intensity_fold_change.png -- per-peptide log2(method/original) distributions
  id_counts.png           -- peptide/PSM IDs per arm (control: MS1-only == orig)
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "single_compare"
FDR = 0.01

# Arm order (original first = the reference); colors for plots.
ORDER = ["original", "dnoise_ms1", "dnoise_msms", "bruker_minesweeper",
         "bruker_minesweeper_mf30", "bruker_minesweeper_strong", "bruker_eh",
         "bruker_background"]
COLOR = {"original": "#1f77b4", "dnoise_ms1": "#d62728", "dnoise_msms": "#2ca02c",
         "bruker_minesweeper": "#ff7f0e", "bruker_minesweeper_mf30": "#e377c2",
         "bruker_minesweeper_strong": "#8c564b", "bruker_eh": "#9467bd",
         "bruker_background": "#17becf"}


def label_of(col: str) -> str:
    """Map an lfq.tsv run column (e.g. 'dnoise_ms1.d' or a path) to its arm label."""
    name = Path(col).name
    return name[:-2] if name.endswith(".d") else name


def load_lfq() -> tuple[pd.DataFrame, dict[str, str]]:
    df = pd.read_csv(OUT / "lfq.tsv", sep="\t")
    if "q_value" in df.columns:
        df = df[df["q_value"] <= FDR]
    dcols = [c for c in df.columns if c.endswith(".d")]
    col_for = {label_of(c): c for c in dcols}
    present = [a for a in ORDER if a in col_for]
    if "original" not in present:
        raise SystemExit("no 'original' column in lfq.tsv -- check staged names")
    # one row per peptide; intensities, 0 -> NaN (not quantified in that run)
    key = "peptide" if "peptide" in df.columns else df.columns[0]
    mat = df[[key] + [col_for[a] for a in present]].copy()
    mat.columns = [key] + present
    mat[present] = mat[present].replace(0, np.nan)
    mat = mat.groupby(key, as_index=False)[present].sum(min_count=1)
    return mat, {a: a for a in present}


def reduction() -> pd.DataFrame:
    """MS1/MS2 peak counts and tdf_bin size per arm, as % of original — read from
    the staged .d folders (symlinks created by 18_single_compare.sh)."""
    from _metrics import stats  # (ms1_peaks, ms2_peaks, tdf_bin bytes)
    base = OUT / "stage"
    rows = {}
    for a in ORDER:
        d = base / f"{a}.d"
        if not d.exists():
            continue
        ms1, ms2, nbytes = stats(d)
        rows[a] = {"ms1_peaks": ms1, "ms2_peaks": ms2, "bin_bytes": nbytes}
    df = pd.DataFrame(rows).T
    if "original" in df.index:
        o = df.loc["original"]
        df["ms1_kept_%"] = (100 * df["ms1_peaks"] / o["ms1_peaks"]).round(1)
        df["ms2_kept_%"] = (100 * df["ms2_peaks"] / o["ms2_peaks"]).round(1)
        df["bin_%_of_orig"] = (100 * df["bin_bytes"] / o["bin_bytes"]).round(1)
    return df


def id_counts() -> pd.DataFrame:
    """Peptide / PSM counts per arm at 1% FDR from results.sage.tsv."""
    df = pd.read_csv(OUT / "results.sage.tsv", sep="\t")
    tgt = df[df["label"] == 1].copy()
    tgt["arm"] = tgt["filename"].map(label_of)
    rows = {}
    for arm, g in tgt.groupby("arm"):
        psm = g[g["spectrum_q"] <= FDR]
        pep = g[g["peptide_q"] <= FDR].drop_duplicates("peptide")
        rows[arm] = {"n_psm": len(psm), "n_peptide": len(pep)}
    return pd.DataFrame(rows).T


def main() -> int:
    mat, _ = load_lfq()
    arms = [c for c in mat.columns if c in ORDER]
    methods = [a for a in arms if a != "original"]
    ref = mat["original"]
    ref_ok = ref.notna()

    # retention + fold-change vs original
    ret_rows, fold = {}, {}
    for m in methods:
        col = mat[m]
        both = ref_ok & col.notna()
        lost = int((ref_ok & col.isna()).sum())       # quant in orig, gone here
        gained = int((~ref_ok & col.notna()).sum())   # new here, absent in orig
        lr = np.log2(col[both] / ref[both])
        fold[m] = lr
        ret_rows[m] = {
            "n_quantified": int(col.notna().sum()),
            "n_vs_original": int(both.sum()),
            "n_lost_vs_original": lost,
            "n_gained_vs_original": gained,
            "median_log2_vs_orig": float(lr.median()) if len(lr) else np.nan,
            "iqr_log2_vs_orig": float(lr.quantile(0.75) - lr.quantile(0.25)) if len(lr) else np.nan,
        }
    ret = pd.DataFrame(ret_rows).T
    ret.loc["original"] = {"n_quantified": int(ref_ok.sum()), "n_vs_original": int(ref_ok.sum()),
                           "n_lost_vs_original": 0, "n_gained_vs_original": 0,
                           "median_log2_vs_orig": 0.0, "iqr_log2_vs_orig": 0.0}
    ret = ret.reindex([a for a in ORDER if a in ret.index])

    # IDs
    ids = id_counts().reindex([a for a in ORDER if a in id_counts().index])

    # Data reduction (peaks + binary size vs original)
    red = reduction()
    red.to_csv(OUT / "data_reduction.csv")
    red_cols = red[["ms1_kept_%", "ms2_kept_%", "bin_%_of_orig"]]

    summary = ret.join(ids, how="left").join(red_cols, how="left")
    OUT.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT / "summary.csv")
    ret.to_csv(OUT / "retention.csv")

    # ---- plot: per-peptide log2 fold-change vs original ----
    fig, ax = plt.subplots(figsize=(10, 5))
    data = [fold[m].clip(-4, 4).dropna().values for m in methods]  # clip for display
    pos = range(len(methods))
    vp = ax.violinplot([d if len(d) else [0] for d in data], positions=list(pos),
                       widths=0.8, showmedians=True, showextrema=False)
    for body, m in zip(vp["bodies"], methods):
        body.set_facecolor(COLOR.get(m, "#999")); body.set_alpha(0.6)
        body.set_edgecolor("black"); body.set_linewidth(0.5)
    vp["cmedians"].set_color("black")
    ax.axhline(0, color="gray", lw=1, ls="--")
    ax.set_xticks(list(pos))
    ax.set_xticklabels(methods, rotation=30, ha="right")
    ax.set_ylabel("log2(method / original) per peptide  [clipped ±4]")
    ax.set_title("MS1 intensity change per peptide vs original (same acquisition)\n"
                 "0 = intensity preserved; negative = signal eroded")
    for i, m in enumerate(methods):
        n = len(fold[m].dropna())
        ax.text(i, -4.3, f"n={n}\nlost {int(ret.loc[m,'n_lost_vs_original'])}",
                ha="center", va="top", fontsize=7)
    ax.set_ylim(-5, 4)
    fig.tight_layout()
    fig.savefig(OUT / "intensity_fold_change.png", dpi=150)
    plt.close(fig)

    # ---- plot: ID counts per arm ----
    if not ids.empty:
        fig, ax = plt.subplots(figsize=(9, 5))
        x = np.arange(len(ids.index))
        for j, metric in enumerate(["n_psm", "n_peptide"]):
            ax.bar(x + (j - 0.5) * 0.4, ids[metric], 0.4, label=metric)
        ax.set_xticks(x); ax.set_xticklabels(ids.index, rotation=30, ha="right")
        ax.set_ylabel("count @ 1% FDR"); ax.legend()
        ax.set_title("Identifications per arm (MS1-only denoisers should match original)")
        fig.tight_layout()
        fig.savefig(OUT / "id_counts.png", dpi=150)
        plt.close(fig)

    pd.set_option("display.width", 200, "display.max_columns", 100)
    print(summary)
    print(f"\nWrote {OUT}/summary.csv, retention.csv, intensity_fold_change.png, id_counts.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
