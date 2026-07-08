#!/usr/bin/env python3
"""Main-text figure: streak filter vs. a matched strict intensity threshold,
EXHAUSTIVE across both acquisition modes and both denoising levels.

Promotes the intensity-baseline control (SI Section S6) to a main-text figure
(paper/figures/fig_intensity.png), a 4-row x 3-column grid:
  row 1: ddaPASEF, MS1-only      (original / streak / threshold)
  row 2: ddaPASEF, MS1+MS/MS     (original / streak_msms / threshold_msms)
  row 3: diaPASEF, MS1-only      (original / streak / threshold)
  row 4: diaPASEF, MS1+MS/MS     (original / streak_msms / threshold_msms)
Each row: (A) quantified proteins, (B) quantified peptides, (C) LFQ accuracy,
both gradients grouped within each bar panel.

ddaPASEF metrics come from Sage's lfq.tsv (via _metrics.py, two-peptide rule).
diaPASEF metrics come from DIA-NN's precomputed results/<ds>/analysis/{summary,
accuracy}.csv (written by 12_analyze_dia.py), which use the identical
two-peptide, >=2-replicate rule.

Also runs a paired permutation test (sign-flip, pure numpy) on the per-run
quantified-protein counts (streak vs. intensity, 18 runs per gradient) for the
two ddaPASEF rows, where per-run data is directly available, and prints the
observed per-run difference and p-value for citation in the text. diaPASEF
per-run counts are not extracted here (DIA-NN's per-run breakdown would need
a separate pr_matrix pass); the figure and SI Table S15 report diaPASEF
aggregates only.

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
ARM_COLOR = {"original": "#0072B2", "streak filter": "#E69F00", "intensity threshold": "#D55E00"}
GRADS = [("5 min", "5min"), ("15 min", "15min")]

# (row title, dataset prefix, {display label -> arm subdir}, uses_dia_csv)
ROWS = [
    ("ddaPASEF, MS1-only", "dda",
     {"original": "original", "streak filter": "denoised", "intensity threshold": "intensity"}, False),
    ("ddaPASEF, MS1+MS/MS", "dda",
     {"original": "original", "streak filter": "msms", "intensity threshold": "intensity_msms"}, False),
    ("diaPASEF, MS1-only", "dia",
     {"original": "original", "streak filter": "denoised", "intensity threshold": "intensity"}, True),
    ("diaPASEF, MS1+MS/MS", "dia",
     {"original": "original", "streak filter": "msms", "intensity threshold": "intensity_msms"}, True),
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


def load_row_data(ds_prefix: str, arm_map: dict, uses_dia_csv: bool):
    """metrics[(grad_label, arm_label)] -> {n_quantified, n_quant_peptide};
    acc[(grad_label, arm_label)] -> list of {expected, observed, species} dicts."""
    metrics, acc = {}, {}
    for glab, g in GRADS:
        ds = f"{ds_prefix}_{g}"
        if uses_dia_csv:
            summ = pd.read_csv(ROOT / f"results/{ds}/analysis/summary.csv", index_col=0)
            accdf = pd.read_csv(ROOT / f"results/{ds}/analysis/accuracy.csv")
            for alab, rsub in arm_map.items():
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
            for alab, rsub in arm_map.items():
                d = ROOT / "results" / ds / rsub
                prot = lfq_table(d)
                metrics[(glab, alab)] = lfq_metrics(prot, d)
                import _metrics as M
                acc[(glab, alab)] = M.pair_accuracy(alab, d)
    return metrics, acc


def main() -> int:
    print("Paired permutation test (per-run quantified proteins, streak vs intensity), ddaPASEF only:")
    for row_title, ds_prefix, arm_map, uses_dia_csv in ROWS:
        if uses_dia_csv:
            continue
        streak_sub = arm_map["streak filter"]
        thresh_sub = arm_map["intensity threshold"]
        for glab, g in GRADS:
            ds = f"{ds_prefix}_{g}"
            s = perrun_protein_counts(ROOT / "results" / ds / streak_sub)
            i = perrun_protein_counts(ROOT / "results" / ds / thresh_sub)
            obs, p = paired_perm_p(s, i)
            print(f"  [{row_title}] {glab}: streak mean {s.mean():.1f}, intensity mean {i.mean():.1f}, "
                  f"delta={obs:+.1f}/run, n={len(s)} runs, p={p:.2e}")

    n_rows = len(ROWS)
    fig, ax = plt.subplots(n_rows, 3, figsize=(13, 4.1 * n_rows))
    plt.rcParams.update({
        "font.size": 12, "axes.titlesize": 13, "axes.labelsize": 12,
        "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10,
    })
    x = np.arange(len(GRADS))
    lim = (-3.5, 4.0)

    for ri, (row_title, ds_prefix, arm_map, uses_dia_csv) in enumerate(ROWS):
        metrics, acc = load_row_data(ds_prefix, arm_map, uses_dia_csv)
        arms = list(arm_map)
        w = 0.8 / len(arms)

        a0, a1, a2 = ax[ri, 0], ax[ri, 1], ax[ri, 2]
        for j, alab in enumerate(arms):
            off = (-(len(arms) - 1) / 2 + j) * w
            prot_vals = [metrics[(g, alab)]["n_quantified"] for g, _ in GRADS]
            pep_vals = [metrics[(g, alab)]["n_quant_peptide"] for g, _ in GRADS]
            b0 = a0.bar(x + off, prot_vals, w, label=alab, color=ARM_COLOR[alab])
            b1 = a1.bar(x + off, pep_vals, w, label=alab, color=ARM_COLOR[alab])
            a0.bar_label(b0, fmt="%d", fontsize=7, padding=1)
            a1.bar_label(b1, fmt="%d", fontsize=7, padding=1)
        a0.set_title("Proteins quantified" if ri == 0 else "")
        a1.set_title("Peptides quantified" if ri == 0 else "")
        a0.set_ylabel(f"{row_title}\nproteins (1% FDR, ≥2 pep.)")
        a1.set_ylabel("quantified peptides")
        a0.set_xticks(x); a0.set_xticklabels([g for g, _ in GRADS])
        a1.set_xticks(x); a1.set_xticklabels([g for g, _ in GRADS])
        if ri == 0:
            a0.legend(fontsize=8)

        a2.plot(lim, lim, color="gray", ls="--", lw=1, zorder=0)
        for alab in arms:
            for glab, _ in GRADS:
                for pt in acc[(glab, alab)]:
                    a2.scatter(pt["expected"], pt["observed"], color=ARM_COLOR[alab],
                               marker=MARK[pt["species"]], s=40, alpha=0.75)
            a2.scatter([], [], color=ARM_COLOR[alab], label=alab)
        a2.set_xlim(*lim); a2.set_ylim(*lim)
        a2.set_xlabel("expected log2 ratio" if ri == n_rows - 1 else "")
        a2.set_ylabel("observed median log2")
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
