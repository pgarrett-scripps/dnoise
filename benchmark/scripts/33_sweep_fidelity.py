#!/usr/bin/env python3
"""Intensity fidelity of each sweep arm vs. the raw (original) arm.

Quantified-peptide *count* alone flatters a brightness cut at high removal (it
can quantify more peptides by keeping only clean, bright apexes), so it is a weak
quality metric. The discriminating question is whether an arm's reported LFQ
intensities stay TRUE to the raw intensities or get distorted. For every arm we
take the per-peptide mean intensity across the replicate set (total-intensity
normalized), join to the raw arm by peptide, and report:

  offset_pct    : global intensity shift vs raw (median log2(arm/raw), as %).
                  An arm can lower (or raise) every peptide by a constant; the
                  intensity threshold biases intensities down ~5% at high removal
                  by discarding dim streak members, the streak filter ~0%.
  med_abs_log2  : median |log2(arm / raw)| AFTER removing that global shift, so it
                  measures per-peptide distortion (relative-quant fidelity), not a
                  constant scale difference (0 = identical shape).
  within_10pct  : fraction of shared peptides within +/-10% of raw after alignment.
  n_shared      : peptides shared with raw (coverage of the comparison)

The global shift and the per-peptide distortion are reported separately: aligning
is the fair way to score relative quantification, but the shift itself is a real,
interpretable side effect worth keeping visible.

A structural streak filter should preserve fidelity (it keeps whole ion streaks,
intensities intact) far better than a per-point intensity threshold at matched
data reduction. Reads the arms produced by 30/32; writes sweep_fidelity.csv and
sweep_fidelity.png (streak vs intensity-threshold vs reduction).
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _metrics import FDR, file_columns, norm_factors

ROOT = Path(__file__).resolve().parents[1]
SWEEP = ROOT / "results" / "dda_15min" / "sweep"
LOG2_10PCT = np.log2(1.10)


def pep_means(arm_dir: Path) -> pd.Series | None:
    """Per-peptide mean normalized LFQ intensity across the replicate set."""
    f = arm_dir / "lfq.tsv"
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
    mean = vals.mean(axis=1)  # mean over replicates where quantified
    s = pd.Series(mean.values, index=df["peptide"].values)
    s = s[s > 0]
    return s.groupby(level=0).mean()


def main() -> int:
    raw = pep_means(SWEEP / "original")
    if raw is None:
        print("no original arm; run 30_param_sweep.py first")
        return 1

    rows = []
    for arm in sorted(d for d in SWEEP.iterdir() if d.is_dir()):
        if arm.name == "original":
            continue
        rj = arm / "reduction.json"
        red = json.loads(rj.read_text()) if rj.exists() else {}
        ms1_kept, ms1_raw = red.get("ms1_kept"), red.get("ms1_raw")
        reduction = (100 * (1 - ms1_kept / ms1_raw)
                     if ms1_kept and ms1_raw else np.nan)
        s = pep_means(arm)
        if s is None:
            continue
        j = pd.concat({"arm": s, "raw": raw}, axis=1).dropna()
        if j.empty:
            continue
        lr = np.log2(j["arm"] / j["raw"])
        # Global shift (median ratio) reported on its own; the distortion metrics
        # are measured on the median-aligned residual so a constant scale
        # difference between arms does not masquerade as per-peptide error.
        offset = float(lr.median())
        aligned = lr - offset
        rows.append({
            "arm": arm.name,
            "kind": red.get("kind", "streak"),
            "gap": red.get("gap"),
            "len": red.get("len"),
            "ms1_reduction_pct": reduction,
            "n_shared": len(j),
            "offset_pct": float((2.0 ** offset - 1.0) * 100.0),
            "med_abs_log2": float(aligned.abs().median()),
            "within_10pct": float((aligned.abs() <= LOG2_10PCT).mean()),
        })

    df = pd.DataFrame(rows).sort_values("med_abs_log2").reset_index(drop=True)
    out_csv = SWEEP / "sweep_fidelity.csv"
    df.to_csv(out_csv, index=False)
    with pd.option_context("display.width", 160,
                           "display.float_format", lambda v: f"{v:.4f}"):
        print(df[["arm", "ms1_reduction_pct", "offset_pct", "med_abs_log2",
                  "within_10pct", "n_shared"]].to_string(index=False))
    print(f"\nwrote {out_csv}")

    _figure(df)
    return 0


def _figure(df: pd.DataFrame) -> None:
    streak = df[df["kind"] == "streak"].copy()
    streak["gap"] = streak["gap"].astype(int)
    intens = df[df["kind"] == "intensity"].sort_values("ms1_reduction_pct")
    cmap = {1: "#0072B2", 2: "#E69F00", 3: "#009E73"}

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

    for ax, ycol, title, ylab in (
        (axes[0], "med_abs_log2",
         "Per-peptide distortion vs raw\n(median-aligned; lower = truer)",
         "median |log2(arm / raw)|, aligned"),
        (axes[1], "within_10pct",
         "Peptides within +/-10% of raw\n(median-aligned; higher = truer)",
         "fraction within +/-10%, aligned"),
        (axes[2], "offset_pct",
         "Global intensity shift vs raw\n(0 = matched scale)",
         "median(arm / raw) - 1  (%)"),
    ):
        for g, sub in streak.groupby("gap"):
            sub = sub.sort_values("ms1_reduction_pct")
            ax.plot(sub["ms1_reduction_pct"], sub[ycol], "-o",
                    color=cmap[g], label=f"streak gap {g}")
        if not intens.empty:
            ax.plot(intens["ms1_reduction_pct"], intens[ycol], "--s",
                    color="#d62728", label="intensity threshold")
        ax.set_xlabel("MS1 data reduction (%)")
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.legend(fontsize=8)
    axes[2].axhline(0, color="0.6", lw=1, ls=":", zorder=0)

    fig.suptitle("dda_15min Condition A (6 reps): LFQ intensity fidelity vs raw "
                 "(per-peptide distortion median-aligned; global shift shown separately)",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = SWEEP / "sweep_fidelity.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    raise SystemExit(main())
