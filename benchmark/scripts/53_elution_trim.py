#!/usr/bin/env python3
"""Future-directions estimate: additional file-size reduction from trimming the
dead-time flanks of the LC gradient (lead-in + column wash) that contain no
confident peptide identifications.

This is an *offline upper-bound*, not an implemented dnoise stage: it needs the
search results first (so it cannot run as a blind one-pass preprocessing step the
way dnoise does), and removing whole frames would require renumbering the `Frames`
table and dependent tables. We report it to size the opportunity.

Method, per run (DDA, both gradients):
  * Elution window from the `original` arm's confident peptide IDs (1% FDR target
    peptides for that file): strict = [min RT, max RT]; margin = strict +/- MARGIN.
  * Per-frame bytes from `Frames.TimsId` (byte offset into analysis.tdf_bin) on the
    *denoised* .d, so "trimmable" bytes are measured on top of dnoise's point
    removal (additional, not double-counted).
  * dnoise byte saving = 1 - denoised_bin / raw_bin (from file sizes on disk).

Writes results/<dataset>/analysis/elution_trim.csv and the SI figure
paper/figures/si_elution_trim.png (example run + stacked per-gradient summary).

Run:  cd benchmark && uv run scripts/53_elution_trim.py
"""

from __future__ import annotations

import os
import sqlite3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIGDIR = ROOT.parent / "paper" / "figures"
DATASETS = {"dda_5min": "5 min", "dda_15min": "15 min"}
MARGIN_MIN = 0.5  # safety pad (each side) so a deployable window never clips IDs.
FDR = 0.01
# Illustrative example run for panel A (must exist in dda_5min).
EXAMPLE = "LFQ_Ultra2_PASEF_5min_50ng_Condition_A_REP1"

plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
})
COLOR = {"dnoise": "#E69F00", "trim_strict": "#0072B2", "trim_margin": "#56B4E9"}


def id_rt_window(sage_tsv: Path, run: str) -> tuple[float, float] | None:
    """[min, max] RT (minutes) of 1% FDR target peptides for `run`, or None."""
    df = pd.read_csv(sage_tsv, sep="\t")
    df = df[df["label"] == 1]
    fcol = next(c for c in df.columns if "file" in c.lower())
    sub = df[df[fcol].astype(str).str.contains(run, regex=False)]
    qcol = "peptide_q" if "peptide_q" in df.columns else "spectrum_q"
    sub = sub[sub[qcol] <= FDR]
    if sub.empty:
        return None
    rt = sub["rt"].to_numpy()
    return float(rt.min()), float(rt.max())


def frame_bytes(dot_d: Path) -> pd.DataFrame:
    """Per-frame RT (minutes) and exact tdf_bin bytes via consecutive TimsId."""
    con = sqlite3.connect(dot_d / "analysis.tdf")
    fr = pd.read_sql(
        "SELECT Id, Time, MsMsType, TimsId FROM Frames ORDER BY TimsId", con
    )
    con.close()
    bin_size = os.path.getsize(dot_d / "analysis.tdf_bin")
    off = fr["TimsId"].to_numpy()
    fr["bytes"] = np.append(np.diff(off), bin_size - off[-1])
    fr["rt_min"] = fr["Time"] / 60.0
    return fr


def out_of_window_bytes(fr: pd.DataFrame, lo: float, hi: float) -> int:
    return int(fr.loc[(fr.rt_min < lo) | (fr.rt_min > hi), "bytes"].sum())


def analyze(dataset: str) -> pd.DataFrame:
    raw_dir = ROOT / "data" / dataset / "raw"
    den_dir = ROOT / "data" / dataset / "denoised"
    sage = ROOT / "results" / dataset / "original" / "results.sage.tsv"
    rows = []
    for raw_d in sorted(raw_dir.glob("*.d")):
        run = raw_d.stem
        den_d = den_dir / raw_d.name
        if not den_d.is_dir():
            continue
        win = id_rt_window(sage, run)
        if win is None:
            continue
        lo, hi = win
        raw_bin = os.path.getsize(raw_d / "analysis.tdf_bin")
        den_bin = os.path.getsize(den_d / "analysis.tdf_bin")
        fr_raw = frame_bytes(raw_d)
        fr_den = frame_bytes(den_d)
        rows.append({
            "dataset": dataset,
            "run": run,
            "gradient_min": float(fr_raw.rt_min.max()),
            "id_lo": lo,
            "id_hi": hi,
            "raw_bin": raw_bin,
            "den_bin": den_bin,
            # Dead-time flank bytes on the RAW file (gross opportunity)...
            "raw_strict_bytes": out_of_window_bytes(fr_raw, lo, hi),
            "raw_margin_bytes": out_of_window_bytes(fr_raw, lo - MARGIN_MIN, hi + MARGIN_MIN),
            # ...and on the DENOISED file (what is left to gain on top of dnoise).
            "den_strict_bytes": out_of_window_bytes(fr_den, lo, hi),
            "den_margin_bytes": out_of_window_bytes(fr_den, lo - MARGIN_MIN, hi + MARGIN_MIN),
        })
    df = pd.DataFrame(rows)
    # All reductions expressed as a fraction of the raw tdf_bin.
    df["dnoise_pct"] = 100 * (df.raw_bin - df.den_bin) / df.raw_bin
    df["flank_strict_pct"] = 100 * df.raw_strict_bytes / df.raw_bin
    df["flank_margin_pct"] = 100 * df.raw_margin_bytes / df.raw_bin
    df["addl_strict_pct"] = 100 * df.den_strict_bytes / df.raw_bin
    df["addl_margin_pct"] = 100 * df.den_margin_bytes / df.raw_bin
    return df


def panel_example(ax) -> None:
    """Per-frame bytes vs RT for one raw run, ID window shaded, flanks marked."""
    raw_d = ROOT / "data" / "dda_5min" / "raw" / f"{EXAMPLE}.d"
    sage = ROOT / "results" / "dda_5min" / "original" / "results.sage.tsv"
    lo, hi = id_rt_window(sage, EXAMPLE)
    fr = frame_bytes(raw_d)
    ms1 = fr[fr.MsMsType == 0]
    ax.fill_between(ms1.rt_min, ms1.bytes / 1e3, step="mid", color="#bbbbbb", lw=0)
    grad = float(fr.rt_min.max())
    ax.axvspan(0, lo, color=COLOR["trim_strict"], alpha=0.18)
    ax.axvspan(hi, grad, color=COLOR["trim_strict"], alpha=0.18, label="dead-time flanks")
    ax.axvline(lo, color=COLOR["trim_strict"], ls="--", lw=1.2)
    ax.axvline(hi, color=COLOR["trim_strict"], ls="--", lw=1.2,
               label="confident-ID RT window")
    ax.set_xlim(0, grad)
    ax.set_xlabel("retention time (min)")
    ax.set_ylabel("MS1 frame size (kB)")
    ax.set_title(f"(a) {EXAMPLE.split('_50ng_')[-1]}: signal vs ID window")
    ax.legend(loc="upper right", framealpha=0.9)


def panel_summary(ax, frames: dict[str, pd.DataFrame]) -> None:
    """Per gradient: dead-time flank fraction of the raw file vs what remains to
    gain on top of dnoise (strict window; mean +/- SD over runs)."""
    labels = list(DATASETS.values())
    x = np.arange(len(labels))
    w = 0.38

    def ms(col):
        return ([frames[d][col].mean() for d in DATASETS],
                [frames[d][col].std() for d in DATASETS])

    flank, flank_sd = ms("flank_strict_pct")
    addl, addl_sd = ms("addl_strict_pct")
    ax.bar(x - w / 2, flank, w, yerr=flank_sd, capsize=3, color=COLOR["trim_strict"],
           label="dead-time flanks (% of raw file)")
    ax.bar(x + w / 2, addl, w, yerr=addl_sd, capsize=3, color=COLOR["trim_margin"],
           label="remaining after dnoise (% of raw)")
    for xi, v in zip(x - w / 2, flank):
        ax.text(xi, v + 0.15, f"{v:.1f}%", ha="center", va="bottom", fontsize=9)
    for xi, v in zip(x + w / 2, addl):
        ax.text(xi, v + 0.15, f"{v:.1f}%", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("tdf_bin fraction (% of raw)")
    ax.set_title("(b) dead-time opportunity vs gain after dnoise")
    ax.legend(loc="upper right", framealpha=0.9)


def main() -> None:
    frames = {}
    for ds in DATASETS:
        df = analyze(ds)
        out = ROOT / "results" / ds / "analysis"
        out.mkdir(parents=True, exist_ok=True)
        df.to_csv(out / "elution_trim.csv", index=False)
        frames[ds] = df
        print(f"{ds}: n={len(df)}  dnoise={df.dnoise_pct.mean():.1f}%  "
              f"flank(strict)={df.flank_strict_pct.mean():.1f}%  "
              f"flank(margin)={df.flank_margin_pct.mean():.1f}%  "
              f"addl-after-dnoise(strict)={df.addl_strict_pct.mean():.1f}%  "
              f"(gradient {df.gradient_min.mean():.1f} min, "
              f"ID window {df.id_lo.mean():.2f}-{df.id_hi.mean():.2f} min)")

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.2))
    panel_example(axes[0])
    panel_summary(axes[1], frames)
    fig.tight_layout()
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / "si_elution_trim.png", dpi=150, bbox_inches="tight")
    print(f"wrote {FIGDIR / 'si_elution_trim.png'}")


if __name__ == "__main__":
    main()
