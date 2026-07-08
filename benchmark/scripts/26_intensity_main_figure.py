#!/usr/bin/env python3
"""Main-text figure: streak filter vs. a matched strict intensity threshold,
EXHAUSTIVE across both acquisition modes and both denoising levels.

Promotes the intensity-baseline control (SI Section S6) to a main-text figure
(paper/figures/fig_intensity.png), a 2-row x 3-column grid:
  row 1: ddaPASEF (MS1-only and MS1+MS/MS combined)
  row 2: diaPASEF (MS1-only and MS1+MS/MS combined)
Each row: (A) quantified proteins, (B) quantified peptides, (C) LFQ accuracy,
both gradients grouped within each bar panel, with 5 arms per gradient group
(original + {streak, threshold} x {MS1-only, MS1+MS/MS}).

ddaPASEF metrics come from Sage's lfq.tsv (via _metrics.py, two-peptide rule).
diaPASEF metrics come from DIA-NN's precomputed results/<ds>/analysis/{summary,
accuracy}.csv (written by 12_analyze_dia.py), which use the identical
two-peptide, >=2-replicate rule.

Also runs a paired permutation test (sign-flip, pure numpy) on the per-run
quantified-protein counts (streak vs. intensity, 18 runs per gradient) for
ddaPASEF at both levels, where per-run data is directly available, and prints
the observed per-run difference and p-value for citation in the text.
diaPASEF per-run counts are not extracted here (DIA-NN's per-run breakdown
would need a separate pr_matrix pass); the figure and SI Table S15 report
diaPASEF aggregates only.

Run: uv run scripts/26_intensity_main_figure.py
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
from _metrics import lfq_metrics, lfq_table

ROOT = Path(__file__).resolve().parents[1]
FIGS = ROOT.parent / "paper" / "figures"
MARK = {"HUMAN": "o", "YEAST": "s", "ECOLI": "^"}
GRADS = [("5 min", "5min"), ("15 min", "15min")]

# 5 arms per row: original (shared), streak/threshold at each of the two
# levels. Colorblind-safe (Wong, Nat. Methods 2011): blue/orange/vermillion
# for MS1-only (original/streak/threshold), bluish-green/reddish-purple for
# the MS1+MS/MS streak/threshold pair (matching 12_analyze_dia.py's palette).
ARMS = ["original", "streak (MS1)", "threshold (MS1)", "streak (MS1+MS/MS)", "threshold (MS1+MS/MS)"]
ARM_COLOR = {
    "original": "#0072B2",
    "streak (MS1)": "#E69F00", "threshold (MS1)": "#D55E00",
    "streak (MS1+MS/MS)": "#009E73", "threshold (MS1+MS/MS)": "#CC79A7",
}
ARM_SUBDIR = {
    "original": "original",
    "streak (MS1)": "denoised", "threshold (MS1)": "intensity",
    "streak (MS1+MS/MS)": "msms", "threshold (MS1+MS/MS)": "intensity_msms",
}

# (row title, dataset prefix, uses_dia_csv)
ROWS = [
    ("ddaPASEF", "dda", False),
    ("diaPASEF", "dia", True),
]


def perrun_protein_counts(arm_dir: Path) -> np.ndarray:
    import _metrics as M
    prot, fcols = M._rollup(arm_dir)
    return prot[fcols].notna().sum(axis=0).values


def paired_perm_p(a, b, iters: int = 100000, seed: int = 0) -> tuple[float, float]:
    """Two-sided paired sign-flip permutation test on the mean difference."""
    d = np.asarray(a, float) - np.asarray(b, float)
    obs = float(d.mean())
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1.0, 1.0], size=(iters, len(d)))
    null = (signs * np.abs(d)).mean(axis=1)
    p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (iters + 1)
    return obs, float(p)


def load_row_data(ds_prefix: str, uses_dia_csv: bool):
    """metrics[(grad_label, arm_label)] -> {n_quantified, n_quant_peptide};
    acc[(grad_label, arm_label)] -> list of {expected, observed, species} dicts."""
    metrics, acc = {}, {}
    for glab, g in GRADS:
        ds = f"{ds_prefix}_{g}"
        if uses_dia_csv:
            summ = pd.read_csv(ROOT / f"results/{ds}/analysis/summary.csv", index_col=0)
            accdf = pd.read_csv(ROOT / f"results/{ds}/analysis/accuracy.csv")
            for alab in ARMS:
                rsub = ARM_SUBDIR[alab]
                if rsub in summ.index:
                    metrics[(glab, alab)] = {
                        "n_quantified": int(summ.loc[rsub, "n_quantified"]),
                        "n_quant_peptide": int(summ.loc[rsub, "n_quant_peptide"]),
                    }
                else:
                    metrics[(glab, alab)] = {"n_quantified": 0, "n_quant_peptide": 0}
                d = accdf[accdf["arm"] == rsub]
                acc[(glab, alab)] = [{"expected": r.expected, "observed": r.observed, "species": r.species}
                                      for r in d.itertuples()]
        else:
            import _metrics as M
            for alab in ARMS:
                rsub = ARM_SUBDIR[alab]
                d = ROOT / "results" / ds / rsub
                prot = lfq_table(d)
                metrics[(glab, alab)] = lfq_metrics(prot, d)
                acc[(glab, alab)] = M.pair_accuracy(rsub, d)
    return metrics, acc


def main() -> int:
    print("Paired permutation test (per-run quantified proteins, streak vs intensity), ddaPASEF only:")
    for level, streak_sub, thresh_sub in [("MS1-only", "denoised", "intensity"), ("MS1+MS/MS", "msms", "intensity_msms")]:
        for glab, g in GRADS:
            ds = f"dda_{g}"
            s = perrun_protein_counts(ROOT / "results" / ds / streak_sub)
            i = perrun_protein_counts(ROOT / "results" / ds / thresh_sub)
            obs, p = paired_perm_p(s, i)
            print(f"  [{level}] {glab}: streak mean {s.mean():.1f}, intensity mean {i.mean():.1f}, "
                  f"delta={obs:+.1f}/run, n={len(s)} runs, p={p:.2e}")

    n_rows = len(ROWS)
    fig, ax = plt.subplots(n_rows, 3, figsize=(13, 4.6 * n_rows))
    plt.rcParams.update({
        "font.size": 12, "axes.titlesize": 13, "axes.labelsize": 12,
        "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 9,
    })
    x = np.arange(len(GRADS))
    xlim = (-3.5, 4.0)

    for ri, (row_title, ds_prefix, uses_dia_csv) in enumerate(ROWS):
        metrics, acc = load_row_data(ds_prefix, uses_dia_csv)
        w = 0.8 / len(ARMS)

        a0, a1, a2 = ax[ri, 0], ax[ri, 1], ax[ri, 2]
        for j, alab in enumerate(ARMS):
            off = (-(len(ARMS) - 1) / 2 + j) * w
            prot_vals = [metrics[(g, alab)]["n_quantified"] for g, _ in GRADS]
            pep_vals = [metrics[(g, alab)]["n_quant_peptide"] for g, _ in GRADS]
            b0 = a0.bar(x + off, prot_vals, w, label=alab, color=ARM_COLOR[alab])
            b1 = a1.bar(x + off, pep_vals, w, label=alab, color=ARM_COLOR[alab])
            a0.bar_label(b0, fmt="%d", fontsize=6, padding=1, rotation=90)
            a1.bar_label(b1, fmt="%d", fontsize=6, padding=1, rotation=90)
        a0.set_title("Proteins quantified" if ri == 0 else "")
        a1.set_title("Peptides quantified" if ri == 0 else "")
        a0.set_ylabel(f"{row_title}\nproteins (1% FDR, ≥2 pep.)")
        a1.set_ylabel("quantified peptides")
        a0.set_xticks(x); a0.set_xticklabels([g for g, _ in GRADS])
        a1.set_xticks(x); a1.set_xticklabels([g for g, _ in GRADS])
        a0.margins(y=0.15); a1.margins(y=0.15)
        if ri == 0:
            a0.legend(fontsize=8, loc="upper left")

        # Residual (observed - expected) vs expected: ideal line is horizontal
        # at 0 rather than a 45deg diagonal, so the 5 overlapping arms/levels
        # separate out far more clearly than an observed-vs-expected view.
        resid = [pt["observed"] - pt["expected"] for alab in ARMS for glab, _ in GRADS
                 for pt in acc[(glab, alab)]]
        pad = max(0.3, 0.15 * (max(resid) - min(resid))) if resid else 0.3
        ylo, yhi = (min(resid) - pad, max(resid) + pad) if resid else (-1, 1)
        a2.axhline(0, color="gray", ls="--", lw=1, zorder=0)
        for alab in ARMS:
            for glab, _ in GRADS:
                for pt in acc[(glab, alab)]:
                    a2.scatter(pt["expected"], pt["observed"] - pt["expected"], color=ARM_COLOR[alab],
                               marker=MARK[pt["species"]], s=40, alpha=0.75)
            a2.scatter([], [], color=ARM_COLOR[alab], label=alab)
        a2.set_xlim(*xlim); a2.set_ylim(ylo, yhi)
        a2.set_xlabel("expected log2 ratio" if ri == n_rows - 1 else "")
        a2.set_ylabel("observed − expected (log2)")
        a2.set_title("LFQ accuracy (○ human △ ecoli □ yeast)" if ri == 0 else "")

    fig.suptitle("Streak filter vs. a matched strict intensity threshold, "
                 "all conditions (equal data removal; ddaPASEF MS/MS untouched at MS1-only)",
                 fontsize=13)
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    out = FIGS / "fig_intensity.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
