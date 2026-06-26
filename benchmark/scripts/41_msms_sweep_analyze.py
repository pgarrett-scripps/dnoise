#!/usr/bin/env python3
"""Score the MS/MS-denoising sweep (40_msms_sweep.py).

For each arm (original + msms gap x len) reads the Sage outputs over the
one-condition replicate set and reports, vs. MS/MS data reduction:

  identifications : PSMs / peptides / proteins at 1% FDR (the primary metric -
                    MS/MS denoising rewrites fragment spectra, so it moves IDs).
  quant           : quantified peptides + run-to-run precision (median CV).
  fidelity        : per-peptide LFQ intensity vs. the raw arm, median-aligned
                    (constant scale removed) so it measures distortion, with the
                    global shift reported separately (cf. 33_sweep_fidelity.py).

Writes results/dda_15min/msms_sweep/msms_sweep_metrics.csv and msms_sweep.png.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _metrics import FDR, file_columns, id_metrics, norm_factors, species_of

ROOT = Path(__file__).resolve().parents[1]
SWEEP = ROOT / "results" / "dda_15min" / "msms_sweep"
MIN_PEP_PER_PROT = 2
LOG2_10PCT = np.log2(1.10)


def quant_metrics(arm: Path) -> dict:
    """Quantified peptides/proteins and run-to-run precision (mirrors 31)."""
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

    cvs = []
    for row, k in zip(vals, nz):
        if k >= 3:
            v = row[row > 0]
            m = v.mean()
            if m > 0:
                cvs.append(v.std(ddof=1) / m)
    median_cv = float(np.median(cvs)) if cvs else np.nan

    sub = df[df["_nz"] >= 1].copy()
    sub["species"] = sub["proteins"].map(species_of)
    sub = sub[sub["species"].notna()]
    sub["protein"] = sub["proteins"].map(lambda p: str(p).split(";")[0])
    npep = sub.groupby("protein")["peptide"].nunique()
    n_prot = int((npep >= MIN_PEP_PER_PROT).sum())
    return {"quant_peptides": quant_any, "quant_proteins": n_prot,
            "median_cv": median_cv}


def pep_means(arm: Path) -> pd.Series | None:
    """Per-peptide mean normalized LFQ intensity across the replicate set."""
    f = arm / "lfq.tsv"
    if not f.exists():
        return None
    df = pd.read_csv(f, sep="\t")
    if "q_value" in df.columns:
        df = df[df["q_value"] <= FDR]
    fcols = file_columns(df)
    if not fcols or df.empty:
        return None
    factors = norm_factors(df, fcols)
    vals = df[fcols].replace(0, np.nan).mul(pd.Series(factors), axis=1)
    s = pd.Series(vals.mean(axis=1).values, index=df["peptide"].values)
    s = s[s > 0]
    return s.groupby(level=0).mean()


def fidelity(arm: Path, raw: pd.Series | None) -> dict:
    """Median-aligned per-peptide intensity fidelity vs. raw (cf. 33)."""
    if raw is None:
        return {}
    s = pep_means(arm)
    if s is None:
        return {}
    j = pd.concat({"arm": s, "raw": raw}, axis=1).dropna()
    if j.empty:
        return {}
    lr = np.log2(j["arm"] / j["raw"])
    offset = float(lr.median())
    aligned = lr - offset
    return {
        "fid_offset_pct": float((2.0 ** offset - 1.0) * 100.0),
        "fid_med_abs_log2": float(aligned.abs().median()),
        "fid_within_10pct": float((aligned.abs() <= LOG2_10PCT).mean()),
    }


def main() -> int:
    if not SWEEP.is_dir():
        print("no msms_sweep results; run 40_msms_sweep.py first")
        return 1
    raw_means = pep_means(SWEEP / "original")

    rows = []
    for arm in sorted(d for d in SWEEP.iterdir() if d.is_dir()):
        rj = arm / "reduction.json"
        red = json.loads(rj.read_text()) if rj.exists() else {}
        ms2_kept, ms2_raw = red.get("ms2_kept"), red.get("ms2_raw")
        ms2_red = (100 * (1 - ms2_kept / ms2_raw)
                   if ms2_kept and ms2_raw else 0.0)
        if not (arm / "results.sage.tsv").exists():
            continue
        row = {"arm": arm.name, "gap": red.get("gap"), "len": red.get("len"),
               "ms2_reduction_pct": ms2_red}
        row.update(id_metrics(arm))
        row.update(quant_metrics(arm))
        if arm.name != "original":
            row.update(fidelity(arm, raw_means))
        rows.append(row)

    if not rows:
        print("no arms found; run 40_msms_sweep.py first")
        return 1

    df = pd.DataFrame(rows).sort_values("ms2_reduction_pct").reset_index(drop=True)
    out_csv = SWEEP / "msms_sweep_metrics.csv"
    df.to_csv(out_csv, index=False)
    cols = ["arm", "ms2_reduction_pct", "n_psm", "n_peptide", "n_protein",
            "quant_peptides", "median_cv", "fid_med_abs_log2", "fid_within_10pct"]
    with pd.option_context("display.width", 200,
                           "display.float_format", lambda v: f"{v:.4f}"):
        print(df[[c for c in cols if c in df.columns]].to_string(index=False))
    print(f"\nwrote {out_csv}")

    _figure(df)
    return 0


def _figure(df: pd.DataFrame) -> None:
    orig = df[df["arm"] == "original"]
    grid = df[df["arm"] != "original"].copy()
    grid["gap"] = grid["gap"].astype(int)
    grid["len"] = grid["len"].astype(int)
    cmap = {2: "#0072B2", 5: "#E69F00", 8: "#009E73"}

    def o(col):
        return float(orig[col].iloc[0]) if len(orig) and col in orig else np.nan

    fig, axd = plt.subplot_mosaic([["pep", "prot"], ["quant", "fid"]],
                                  figsize=(13, 9.5))

    def series(ax, ycol, ylab, title, oval=None):
        for g, sub in grid.groupby("gap"):
            sub = sub.sort_values("ms2_reduction_pct")
            ax.plot(sub["ms2_reduction_pct"], sub[ycol], "-o",
                    color=cmap[g], label=f"msms gap {g}")
            for _, r in sub.iterrows():
                ax.annotate(f"L{int(r['len'])}",
                            (r["ms2_reduction_pct"], r[ycol]),
                            textcoords="offset points", xytext=(4, 4), fontsize=8)
        if oval is not None and np.isfinite(oval):
            ax.scatter([0], [oval], marker="*", s=180, color="black",
                       zorder=5, label="original (raw)")
        ax.set_xlabel("MS/MS data reduction (%)")
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.legend(fontsize=8)

    series(axd["pep"], "n_peptide", "peptides (1% FDR)",
           "Peptide IDs vs MS/MS reduction", o("n_peptide"))
    series(axd["prot"], "quant_proteins", "quantified proteins (>=2 peptides)",
           "Quantified proteins vs MS/MS reduction", o("quant_proteins"))
    series(axd["quant"], "quant_peptides", "quantified peptides",
           "Quantified peptides vs MS/MS reduction", o("quant_peptides"))
    if "fid_med_abs_log2" in grid:
        series(axd["fid"], "fid_med_abs_log2",
               "median |log2(arm / raw)|, aligned",
               "LFQ fidelity vs MS/MS reduction (lower = truer)")

    fig.suptitle("dda_15min Condition A (6 reps): MS/MS-denoising parameter sweep",
                 fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = SWEEP / "msms_sweep.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")
    # Copy into the paper figures dir for the SI (cf. 34/43_*).
    paper_fig = ROOT.parent / "paper" / "figures" / "si_msms_sweep.png"
    if paper_fig.parent.is_dir():
        import shutil
        shutil.copy(out, paper_fig)
        print(f"copied -> {paper_fig}")


if __name__ == "__main__":
    raise SystemExit(main())
