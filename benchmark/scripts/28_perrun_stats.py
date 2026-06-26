#!/usr/bin/env python3
"""Per-run paired statistics for the original vs MS1-denoised arms (reviewer ask):
substantiate the quantified-protein coverage change with a confidence interval and
test, and report the run-to-run SD that frames the "within run-to-run variance"
claim. Covers ddaPASEF (per-run quantified-protein counts) and diaPASEF (per-run
precursor and protein-group counts), both gradients.

For each gradient we pair the 18 runs (matched by .d basename), report each arm's
per-run mean +/- SD, the paired mean difference (denoised - original) with a
paired percentile bootstrap 95% CI (10,000 resamples of the 18 run-differences),
and a two-sided paired sign-flip permutation p-value (100,000 permutations, seed 0)
-- the same test family used for the intensity-threshold control.

Writes paper/si/perrun_stats.typ (a #table) and prints the numbers.
Usage: uv run scripts/28_perrun_stats.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _metrics as M  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SI = ROOT.parent / "paper" / "si"
BOOT, PERM, SEED = 10_000, 100_000, 0


def _basename_key(col: str) -> str:
    return Path(str(col)).name.replace(".d", "")


def dda_perrun(arm_dir: Path) -> pd.Series:
    """Per-run quantified-protein count (two-peptide rule), indexed by run basename."""
    prot, fcols = M._rollup(arm_dir)
    counts = prot[fcols].notna().sum(axis=0)
    counts.index = [_basename_key(c) for c in fcols]
    return counts


def dia_perrun(arm_dir: Path, matrix: str) -> pd.Series:
    """Per-run non-null count from a DIA-NN matrix (pr_matrix=precursors,
    pg_matrix=protein groups), indexed by run basename."""
    df = pd.read_csv(arm_dir / matrix, sep="\t")
    run_cols = [c for c in df.columns if str(c).endswith(".d")]
    counts = df[run_cols].notna().sum(axis=0)
    counts.index = [_basename_key(c) for c in run_cols]
    return counts


def paired(orig: pd.Series, den: pd.Series) -> dict:
    """Pair two per-run count Series by run; paired diff stats."""
    runs = sorted(set(orig.index) & set(den.index))
    o = orig.loc[runs].to_numpy(float)
    d = den.loc[runs].to_numpy(float)
    diff = d - o
    n = len(diff)
    obs = float(diff.mean())
    rng = np.random.default_rng(SEED)
    # paired bootstrap CI on the mean difference
    bidx = rng.integers(0, n, size=(BOOT, n))
    bmeans = diff[bidx].mean(axis=1)
    lo, hi = float(np.percentile(bmeans, 2.5)), float(np.percentile(bmeans, 97.5))
    # two-sided paired sign-flip permutation p
    signs = rng.choice([-1.0, 1.0], size=(PERM, n))
    null = (signs * np.abs(diff)).mean(axis=1)
    p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (PERM + 1)
    return {
        "n": n,
        "orig_mean": float(o.mean()), "orig_sd": float(o.std(ddof=1)),
        "den_mean": float(d.mean()), "den_sd": float(d.std(ddof=1)),
        "diff": obs, "ci_lo": lo, "ci_hi": hi, "p": float(p),
        "diff_min": float(diff.min()), "diff_max": float(diff.max()),
    }


def fmt_p(p: float) -> str:
    return f"{p:.2g}" if p > 1.0 / (PERM + 1) + 1e-12 else f"{1.0/(PERM+1):.1e}"


def main() -> int:
    rows = []  # (gradient, metric, stats)
    print("=== ddaPASEF: per-run quantified-protein counts (original vs MS1) ===")
    for g in ["dda_5min", "dda_15min"]:
        s = paired(dda_perrun(ROOT / "results" / g / "original"),
                   dda_perrun(ROOT / "results" / g / "denoised"))
        rows.append((g.replace("dda_", "").replace("min", " min"), "Proteins/run", s))
        print(f"  {g}: orig {s['orig_mean']:.1f}±{s['orig_sd']:.1f}, "
              f"MS1 {s['den_mean']:.1f}±{s['den_sd']:.1f}, "
              f"Δ={s['diff']:+.1f} [{s['ci_lo']:+.1f},{s['ci_hi']:+.1f}] "
              f"(range {s['diff_min']:+.0f}..{s['diff_max']:+.0f}), p={fmt_p(s['p'])}")

    print("\n=== diaPASEF: per-run precursor + protein-group counts (original vs MS1) ===")
    for g in ["dia_5min", "dia_15min"]:
        for label, mat in [("Precursors/run", "report.pr_matrix.tsv"),
                           ("Protein groups/run", "report.pg_matrix.tsv")]:
            s = paired(dia_perrun(ROOT / "results" / g / "original", mat),
                       dia_perrun(ROOT / "results" / g / "denoised", mat))
            rows.append((g.replace("dia_", "").replace("min", " min"), label, s))
            print(f"  {g} {label}: orig {s['orig_mean']:.0f}±{s['orig_sd']:.0f}, "
                  f"MS1 {s['den_mean']:.0f}±{s['den_sd']:.0f}, "
                  f"Δ={s['diff']:+.1f} [{s['ci_lo']:+.1f},{s['ci_hi']:+.1f}], p={fmt_p(s['p'])}")

    # ---- emit typst table ----
    SI.mkdir(parents=True, exist_ok=True)
    lines = [
        "#table(",
        "  columns: 7,",
        "  align: (left, left, center, center, center, center, center),",
        "  table.header(",
        "    [Gradient], [Metric], [Original (mean ± SD)], [MS1 (mean ± SD)],",
        "    [$Delta$ (MS1 #sym.minus orig.)], [95% CI of $Delta$], [perm. $p$],",
        "  ),",
    ]
    for grad, metric, s in rows:
        om = f"{s['orig_mean']:.0f} #sym.plus.minus {s['orig_sd']:.0f}"
        dm = f"{s['den_mean']:.0f} #sym.plus.minus {s['den_sd']:.0f}"
        d = f"{s['diff']:+.1f}".replace("-", "#sym.minus ")
        ci = f"[{s['ci_lo']:+.1f}, {s['ci_hi']:+.1f}]".replace("-", "#sym.minus ")
        lines.append(
            f"  [{grad}], [{metric}], [{om}], [{dm}], [{d}], [{ci}], [{fmt_p(s['p'])}],"
        )
    lines.append(")")
    (SI / "perrun_stats.typ").write_text("\n".join(lines) + "\n")
    print(f"\nwrote {SI/'perrun_stats.typ'} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
