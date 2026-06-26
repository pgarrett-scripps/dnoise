#!/usr/bin/env python3
"""Score the MS1 selection-polygon comparison (42_polygon_compare.py).

Reports, for original / nopoly / poly, what the `--ms1-polygon` gate adds:
  compression : % MS1 peaks removed and % tdf_bin bytes removed.
  IDs         : PSMs / peptides / proteins at 1% FDR (expected ~unchanged - the
                gate is MS1-only and DDA IDs come from MS/MS).
  quant       : quantified peptides/proteins + run-to-run precision (median CV).
  fidelity    : per-peptide LFQ intensity vs. raw, median-aligned (cf. 33).

Writes results/dda_15min/polygon_compare/polygon_compare_metrics.csv and
polygon_compare.png (paired nopoly-vs-poly bars).
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
CMP = ROOT / "results" / "dda_15min" / "polygon_compare"
MIN_PEP_PER_PROT = 2
LOG2_10PCT = np.log2(1.10)


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
    cvs = []
    for row, k in zip(vals, nz):
        if k >= 3:
            v = row[row > 0]
            if v.mean() > 0:
                cvs.append(v.std(ddof=1) / v.mean())
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
    if not CMP.is_dir():
        print("no polygon_compare results; run 42_polygon_compare.py first")
        return 1
    raw_means = pep_means(CMP / "original")

    order = ["original", "nopoly", "poly"]
    rows = []
    for name in order:
        arm = CMP / name
        if not (arm / "results.sage.tsv").exists():
            continue
        rj = arm / "reduction.json"
        red = json.loads(rj.read_text()) if rj.exists() else {}
        m1k, m1r = red.get("ms1_kept"), red.get("ms1_raw")
        bk, br = red.get("bytes_kept"), red.get("bytes_raw")
        row = {
            "arm": name,
            "ms1_reduction_pct": (100 * (1 - m1k / m1r)) if m1k and m1r else 0.0,
            "bytes_reduction_pct": (100 * (1 - bk / br)) if bk and br else 0.0,
        }
        row.update(id_metrics(arm))
        row.update(quant_metrics(arm))
        if name != "original":
            row.update(fidelity(arm, raw_means))
        rows.append(row)

    if not rows:
        print("no arms found; run 42_polygon_compare.py first")
        return 1
    df = pd.DataFrame(rows).set_index("arm").reindex([r["arm"] for r in rows])
    out_csv = CMP / "polygon_compare_metrics.csv"
    df.to_csv(out_csv)
    cols = ["ms1_reduction_pct", "bytes_reduction_pct", "n_psm", "n_peptide",
            "n_protein", "quant_peptides", "quant_proteins", "median_cv",
            "fid_med_abs_log2", "fid_within_10pct"]
    with pd.option_context("display.width", 200,
                           "display.float_format", lambda v: f"{v:.4f}"):
        print(df[[c for c in cols if c in df.columns]].to_string())

    # poly-vs-nopoly deltas (the gate's marginal effect).
    if {"nopoly", "poly"} <= set(df.index):
        d = df.loc["poly"] - df.loc["nopoly"]
        print("\npoly - nopoly:")
        print(f"  +{d['ms1_reduction_pct']:.1f} pts MS1 removed, "
              f"+{d['bytes_reduction_pct']:.1f} pts bytes removed")
        print(f"  peptides {d.get('n_peptide', 0):+.0f}, "
              f"proteins {d.get('n_protein', 0):+.0f}, "
              f"quant_peptides {d.get('quant_peptides', 0):+.0f}")
    print(f"\nwrote {out_csv}")

    _figure(df)
    return 0


def _figure(df: pd.DataFrame) -> None:
    arms = [a for a in ["nopoly", "poly"] if a in df.index]
    colors = {"nopoly": "#999999", "poly": "#0072B2"}
    pct = lambda v: f"{v:.1f}"
    cnt = lambda v: f"{v:,.0f}"
    fid = lambda v: f"{v:.3f}"
    panels = [
        ("ms1_reduction_pct", "% MS1 peaks removed", pct),
        ("bytes_reduction_pct", "% tdf_bin bytes removed", pct),
        ("n_peptide", "peptides (1% FDR)", cnt),
        ("n_protein", "proteins (1% FDR)", cnt),
        ("quant_peptides", "quantified peptides", cnt),
        ("fid_med_abs_log2", "LFQ distortion vs raw\n(median |log2|, aligned)", fid),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
    orig = df.loc["original"] if "original" in df.index else None
    for ax, (col, title, fmt) in zip(axes.flat, panels):
        if col not in df.columns:
            ax.set_visible(False)
            continue
        vals = [df.loc[a, col] for a in arms]
        ax.bar(arms, vals, color=[colors[a] for a in arms], width=0.6)
        for i, v in enumerate(vals):
            ax.text(i, v, fmt(v), ha="center", va="bottom", fontsize=10)
        # ID/quant panels: draw the raw original as a reference line.
        if orig is not None and col in ("n_peptide", "n_protein", "quant_peptides"):
            ov = orig.get(col)
            if ov is not None and np.isfinite(ov):
                ax.axhline(ov, color="black", ls="--", lw=1, label="original (raw)")
                ax.legend(fontsize=8)
        ax.set_title(title)
        ax.margins(y=0.18)
    fig.suptitle("dda_15min Condition A (6 reps): MS1 selection-polygon gate "
                 "(no-polygon vs. polygon, MS1-only denoising)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = CMP / "polygon_compare.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")
    # Copy into the paper figures dir for the SI (cf. 34_sweep_si_table.py).
    paper_fig = ROOT.parent / "paper" / "figures" / "si_polygon_compare.png"
    if paper_fig.parent.is_dir():
        import shutil
        shutil.copy(out, paper_fig)
        print(f"copied -> {paper_fig}")


if __name__ == "__main__":
    raise SystemExit(main())
