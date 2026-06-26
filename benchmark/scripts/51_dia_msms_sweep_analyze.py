#!/usr/bin/env python3
"""Score the DIA whole-frame MS/MS denoising sweep (50_dia_msms_sweep.py).

For each arm reports the compression-vs-quant trade-off on dia_15min Condition A:
  compression : % MS1 / MS2 peaks and % tdf_bin bytes removed (from reduction.json).
  quant       : protein groups and precursors quantified (>=1 of 6 reps) and their
                run-to-run median CV, from DIA-NN's report.pg/pr_matrix.tsv.

Anchors (single-pass, same library): the MS1-only arm (no MS2 denoising; the
len-0 reference) from results/dia_15min/quant_bench/denoised. The raw `original`
arm is shown for context but was searched with --reanalyse, so it is a soft
baseline only (annotated).

Writes results/dia_15min/msms_sweep/dia_msms_sweep_metrics.csv and
dia_msms_sweep.png (bytes removed vs proteins/precursors quantified).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SWEEP = ROOT / "results" / "dia_15min" / "msms_sweep"
QB = ROOT / "results" / "dia_15min" / "quant_bench"
ORIG = ROOT / "results" / "dia_15min" / "original"


def quant(matrix: Path) -> tuple[int, float]:
    """(# quantified rows with >=1 of 6 Condition-A reps, median CV over >=3)."""
    df = pd.read_csv(matrix, sep="\t")
    cols = sorted([c for c in df.columns if "Condition_A_REP" in c],
                  key=lambda c: int(re.search(r"REP(\d+)", c).group(1)))
    M = df[cols].apply(pd.to_numeric, errors="coerce").to_numpy()
    nz = (M > 0).sum(axis=1)
    n = int((nz >= 1).sum())
    cvs = [v[v > 0].std(ddof=1) / v[v > 0].mean()
           for v, c in zip(M, nz) if c >= 3 and v[v > 0].mean() > 0]
    return n, float(np.median(cvs)) if cvs else float("nan")


def quant_row(resdir: Path) -> dict:
    pg, pcv = quant(resdir / "report.pg_matrix.tsv")
    pr, rcv = quant(resdir / "report.pr_matrix.tsv")
    return {"proteins": pg, "protein_cv": pcv, "precursors": pr, "precursor_cv": rcv}


def main() -> int:
    rows = []
    # Raw anchor (context only; --reanalyse, soft baseline).
    if (ORIG / "report.pg_matrix.tsv").exists():
        rows.append({"arm": "original (raw)*", "gap": None, "len": None,
                     "ms1_pct": 0.0, "ms2_pct": 0.0, "bytes_pct": 0.0,
                     **quant_row(ORIG)})
    # MS1-only anchor (len 0; single-pass) -- compression from the on-disk arm.
    ms1_dir = QB / "denoised"
    if (ms1_dir / "report.pg_matrix.tsv").exists():
        # MS1-only compression: MS2 untouched (0%); MS1/bytes from data/.../denoised.
        rows.append({"arm": "MS1-only (L0)", "gap": 0, "len": 0,
                     "ms1_pct": 91.0, "ms2_pct": 0.0, "bytes_pct": 32.9,
                     **quant_row(ms1_dir)})

    for res in sorted(SWEEP.glob("g*_l*")):
        rj = res / "reduction.json"
        if not rj.exists() or not (res / "report.pg_matrix.tsv").exists():
            continue
        r = json.loads(rj.read_text())
        rows.append({
            "arm": res.name, "gap": r.get("gap"), "len": r.get("len"),
            "ms1_pct": 100 * (1 - r["ms1_kept"] / r["ms1_raw"]),
            "ms2_pct": 100 * (1 - r["ms2_kept"] / r["ms2_raw"]),
            "bytes_pct": 100 * (1 - r["bytes_kept"] / r["bytes_raw"]),
            **quant_row(res)})

    if not rows:
        print("no arms found; run 50_dia_msms_sweep.py first")
        return 1
    df = pd.DataFrame(rows)
    out_csv = SWEEP / "dia_msms_sweep_metrics.csv"
    df.to_csv(out_csv, index=False)
    with pd.option_context("display.width", 200, "display.float_format", lambda v: f"{v:.3f}"):
        print(df.to_string(index=False))
    print(f"\nwrote {out_csv}")

    _figure(df)
    return 0


def _figure(df: pd.DataFrame) -> None:
    sweep = df[df["gap"].fillna(-1) > 0].copy()
    ms1 = df[df["arm"] == "MS1-only (L0)"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    gapcol = {3: "#0072B2", 5: "#E69F00", 8: "#009E73"}
    for ax, ycol, ylab in ((axes[0], "proteins", "protein groups quantified"),
                           (axes[1], "precursors", "precursors quantified")):
        for g, sub in sweep.groupby("gap"):
            sub = sub.sort_values("bytes_pct")
            ax.plot(sub["bytes_pct"], sub[ycol], "-o", color=gapcol.get(int(g), "gray"),
                    label=f"msms gap {int(g)}")
            for _, r in sub.iterrows():
                ax.annotate(f"L{int(r['len'])}", (r["bytes_pct"], r[ycol]),
                            textcoords="offset points", xytext=(4, 4), fontsize=8)
        if len(ms1):
            ax.scatter(ms1["bytes_pct"], ms1[ycol], marker="*", s=180, color="black",
                       zorder=5, label="MS1-only (no MS2 denoise)")
        ax.set_xlabel("tdf_bin bytes removed (%)")
        ax.set_ylabel(ylab)
        ax.legend(fontsize=8)
    fig.suptitle("dia_15min Condition A (6 reps): whole-frame MS/MS denoising "
                 "compression vs. quant (single-pass, empirical lib)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = SWEEP / "dia_msms_sweep.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")
    paper_fig = ROOT.parent / "paper" / "figures" / "si_dia_msms_sweep.png"
    if paper_fig.parent.is_dir():
        import shutil
        shutil.copy(out, paper_fig)
        print(f"copied -> {paper_fig}")


if __name__ == "__main__":
    raise SystemExit(main())
