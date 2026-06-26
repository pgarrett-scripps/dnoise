#!/usr/bin/env python3
"""Score the default-parameter sweep (30_param_sweep.py).

For each arm (original + gap x len) reads the Sage outputs over the one-condition
replicate set and reports, against MS1 data reduction:
  - quantified peptides  : distinct peptides in lfq.tsv (q<=0.01) with a nonzero
                           LFQ intensity in >=1 replicate (the headline metric)
  - complete peptides    : quantified in ALL replicates (data completeness)
  - quantified proteins  : species-clean proteins with >=2 quantified peptides
  - median peptide CV    : run-to-run precision across the replicates
  - peptide / PSM IDs    : 1% FDR identifications (sanity; MS1 denoise shouldn't
                           move these much since IDs come from MS/MS)

Writes results/dda_15min/sweep/sweep_metrics.csv and sweep_quant.png, and prints
a table sorted by quantified peptides.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _metrics import FDR, file_columns, species_of

ROOT = Path(__file__).resolve().parents[1]
SWEEP = ROOT / "results" / "dda_15min" / "sweep"
MIN_PEP_PER_PROT = 2


def id_counts(arm: Path) -> tuple[int, int]:
    """(PSMs, distinct peptides) at 1% FDR from results.sage.tsv."""
    f = arm / "results.sage.tsv"
    if not f.exists():
        return 0, 0
    df = pd.read_csv(f, sep="\t")
    tgt = df[df["label"] == 1]
    psm = int((tgt["spectrum_q"] <= FDR).sum())
    pep = int(tgt.loc[tgt["peptide_q"] <= FDR, "peptide"].nunique())
    return psm, pep


def quant_metrics(arm: Path) -> dict:
    f = arm / "lfq.tsv"
    if not f.exists():
        return {}
    df = pd.read_csv(f, sep="\t")
    if "q_value" in df.columns:
        df = df[df["q_value"] <= FDR]
    fcols = file_columns(df)
    if not fcols or df.empty:
        return {}
    n = len(fcols)
    vals = df[fcols].to_numpy(dtype=float)
    nz = (vals > 0).sum(axis=1)
    df = df.assign(_nz=nz)

    quant_any = int(df.loc[df["_nz"] >= 1, "peptide"].nunique())
    quant_all = int(df.loc[df["_nz"] >= n, "peptide"].nunique())

    # Run-to-run precision: CV over nonzero replicate intensities per peptide row
    # (require >=3 observations so the CV is meaningful).
    cvs = []
    for row, k in zip(vals, nz):
        if k >= 3:
            v = row[row > 0]
            m = v.mean()
            if m > 0:
                cvs.append(v.std(ddof=1) / m)
    median_cv = float(np.median(cvs)) if cvs else np.nan

    # Quantified proteins: species-clean peptides, two-peptide rule.
    sub = df[df["_nz"] >= 1].copy()
    sub["species"] = sub["proteins"].map(species_of)
    sub = sub[sub["species"].notna()]
    sub["protein"] = sub["proteins"].map(lambda p: str(p).split(";")[0])
    npep = sub.groupby("protein")["peptide"].nunique()
    n_prot = int((npep >= MIN_PEP_PER_PROT).sum())
    clean_pep = int(sub.loc[sub["protein"].isin(npep[npep >= MIN_PEP_PER_PROT].index),
                            "peptide"].nunique())

    return {
        "quant_peptides": quant_any,
        "complete_peptides": quant_all,
        "quant_proteins": n_prot,
        "clean_quant_peptides": clean_pep,
        "median_cv": median_cv,
    }


def main() -> int:
    arms = sorted([d for d in SWEEP.iterdir() if d.is_dir()])
    rows = []
    for arm in arms:
        red = {}
        rj = arm / "reduction.json"
        if rj.exists():
            red = json.loads(rj.read_text())
        ms1_kept = red.get("ms1_kept")
        ms1_raw = red.get("ms1_raw")
        reduction = (100 * (1 - ms1_kept / ms1_raw)
                     if ms1_kept and ms1_raw else np.nan)
        psm, pep = id_counts(arm)
        rec = {
            "arm": arm.name,
            "kind": red.get("kind", "streak"),
            "T": red.get("T"),
            "gap": red.get("gap"),
            "len": red.get("len"),
            "ms1_reduction_pct": reduction,
            "psm_ids": psm,
            "peptide_ids": pep,
            **quant_metrics(arm),
        }
        rows.append(rec)

    df = pd.DataFrame(rows)
    if df.empty:
        print("no sweep arms found; run 30_param_sweep.py first")
        return 1
    df = df.sort_values("quant_peptides", ascending=False).reset_index(drop=True)
    out_csv = SWEEP / "sweep_metrics.csv"
    df.to_csv(out_csv, index=False)

    cols = ["arm", "ms1_reduction_pct", "quant_peptides", "complete_peptides",
            "quant_proteins", "median_cv", "peptide_ids", "psm_ids"]
    with pd.option_context("display.max_rows", None, "display.width", 160,
                           "display.float_format", lambda v: f"{v:.3f}"):
        print(df[cols].to_string(index=False))
    print(f"\nwrote {out_csv}")

    _figure(df)
    return 0


def _figure(df: pd.DataFrame) -> None:
    grid = df[(df["kind"] == "streak") & df["gap"].notna()].copy()
    grid["gap"] = grid["gap"].astype(int)
    grid["len"] = grid["len"].astype(int)
    intens = df[df["kind"] == "intensity"].sort_values("ms1_reduction_pct")
    orig = df[df["arm"] == "original"]
    o_red = float(orig["ms1_reduction_pct"].iloc[0]) if len(orig) else 0.0
    o_pep = int(orig["quant_peptides"].iloc[0]) if len(orig) else 0
    o_cv = float(orig["median_cv"].iloc[0]) if len(orig) else np.nan

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    cmap = {1: "#0072B2", 2: "#E69F00", 3: "#009E73"}

    # Panel A: quantified peptides vs MS1 reduction.
    ax = axes[0]
    for g, sub in grid.groupby("gap"):
        sub = sub.sort_values("ms1_reduction_pct")
        ax.plot(sub["ms1_reduction_pct"], sub["quant_peptides"], "-o",
                color=cmap[g], label=f"gap {g}")
        for _, r in sub.iterrows():
            ax.annotate(f"L{int(r['len'])}",
                        (r["ms1_reduction_pct"], r["quant_peptides"]),
                        fontsize=7, xytext=(3, 3), textcoords="offset points")
    if not intens.empty:
        ax.plot(intens["ms1_reduction_pct"], intens["quant_peptides"], "--s",
                color="#d62728", label="intensity threshold", zorder=4)
    ax.axhline(o_pep, ls="--", color="gray", lw=1)
    ax.scatter([o_red], [o_pep], marker="*", s=160, color="k", zorder=5,
               label="original")
    ax.set_xlabel("MS1 data reduction (%)")
    ax.set_ylabel("quantified peptides")
    ax.set_title("Quantified peptides vs reduction")
    ax.legend(fontsize=9)

    # Panel B: heatmap of quantified peptides over the gap x len grid.
    ax = axes[1]
    piv = grid.pivot(index="gap", columns="len", values="quant_peptides")
    im = ax.imshow(piv.values, aspect="auto", cmap="viridis", origin="lower")
    ax.set_xticks(range(len(piv.columns)))
    ax.set_xticklabels(piv.columns)
    ax.set_yticks(range(len(piv.index)))
    ax.set_yticklabels(piv.index)
    ax.set_xlabel("min_feature_length")
    ax.set_ylabel("max_internal_gap")
    ax.set_title("Quantified peptides")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{int(v)}", ha="center", va="center",
                        color="w", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046)

    # Panel C: median CV vs reduction (precision; lower is better).
    ax = axes[2]
    for g, sub in grid.groupby("gap"):
        sub = sub.sort_values("ms1_reduction_pct")
        ax.plot(sub["ms1_reduction_pct"], sub["median_cv"], "-o",
                color=cmap[g], label=f"gap {g}")
    if not intens.empty:
        ax.plot(intens["ms1_reduction_pct"], intens["median_cv"], "--s",
                color="#d62728", label="intensity threshold")
    if not np.isnan(o_cv):
        ax.axhline(o_cv, ls="--", color="gray", lw=1)
        ax.scatter([o_red], [o_cv], marker="*", s=160, color="k", zorder=5,
                   label="original")
    ax.set_xlabel("MS1 data reduction (%)")
    ax.set_ylabel("median peptide CV")
    ax.set_title("Quant precision vs reduction")
    ax.legend(fontsize=9)

    fig.suptitle("dda_15min Condition A (6 replicates): default-parameter sweep",
                 fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = SWEEP / "sweep_quant.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    raise SystemExit(main())
