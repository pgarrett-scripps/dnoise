#!/usr/bin/env python3
"""S6 table: streak filter vs. a matched strict intensity threshold.

Self-contained metrics for the Supporting-Information intensity-baseline
comparison (paper Section S6). Unlike 06_analyze.py / _metrics.py (which report
quant at >=2 replicates per condition), this script adds a *two-peptides-per-
protein* requirement for quant reporting and also reports the number of quantified
peptides (distinct sequences mapping to a reported protein). It is intentionally
kept separate so the main-text quant definition is unchanged.

For each arm it reads lfq.tsv, total-intensity-normalizes, keeps species-clean
peptides, requires a protein to be supported by >=2 distinct peptide sequences,
and (as in the main pipeline) requires >=2 replicates in each compared condition.
Run: uv run scripts/25_intensity_compare.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _metrics import (  # noqa: E402
    EXPECTED,
    FDR,
    SPECIES,
    condition_of,
    expected_log2,
    file_columns,
    norm_factors,
    species_of,
)

ROOT = Path(__file__).resolve().parents[1]
ARMS = {
    "original": "original",
    "denoised": "streak",
    "intensity": "intensity",
    "msms": "streak_msms",
    "intensity_msms": "intensity_msms",
}
MIN_PEPTIDES = 2  # two-peptide rule for quant reporting (S6 only)


def quant_metrics(arm_dir: Path) -> dict:
    """Quant metrics with the two-peptide-per-protein rule, headline pair A/B."""
    df = pd.read_csv(arm_dir / "lfq.tsv", sep="\t")
    if "q_value" in df.columns:
        df = df[df["q_value"] <= FDR]
    fcols = file_columns(df)
    if not fcols or df.empty:
        return {}

    factors = norm_factors(df, fcols)
    df = df.assign(species=df["proteins"].map(species_of))
    df = df[df["species"].notna()].copy()
    df["protein"] = df["proteins"].map(lambda p: str(p).split(";")[0])

    # Two-peptide rule: keep only proteins supported by >=2 distinct peptides.
    pep_per_prot = df.groupby("protein")["peptide"].nunique()
    keep_prot = set(pep_per_prot[pep_per_prot >= MIN_PEPTIDES].index)
    df = df[df["protein"].isin(keep_prot)]

    # Normalized per-file intensities, rolled up peptide -> (protein, species).
    vals = df[fcols].replace(0, np.nan).mul(pd.Series(factors), axis=1)
    rolled = pd.concat([df[["protein", "species"]], vals], axis=1)
    prot = rolled.groupby(["protein", "species"], as_index=False)[fcols].sum(min_count=1)

    a_cols = [c for c in fcols if condition_of(c) == "A"]
    b_cols = [c for c in fcols if condition_of(c) == "B"]
    a_ok = prot[a_cols].notna().sum(axis=1) >= 2
    b_ok = prot[b_cols].notna().sum(axis=1) >= 2
    prot = prot[a_ok & b_ok].copy()

    prot["mean_A"] = prot[a_cols].mean(axis=1)
    prot["mean_B"] = prot[b_cols].mean(axis=1)
    prot["log2_ratio"] = np.log2(prot["mean_A"] / prot["mean_B"])
    prot["cv_A"] = prot[a_cols].std(axis=1, ddof=1) / prot[a_cols].mean(axis=1)
    prot["cv_B"] = prot[b_cols].std(axis=1, ddof=1) / prot[b_cols].mean(axis=1)

    # Quantified peptides: distinct sequences mapping to a quantified protein.
    quant_proteins = set(prot["protein"])
    n_quant_pep = int(df[df["protein"].isin(quant_proteins)]["peptide"].nunique())

    cv = pd.concat([prot["cv_A"], prot["cv_B"]]).dropna()
    m = {
        "n_quant_protein": len(prot),
        "n_quant_peptide": n_quant_pep,
        "median_cv": float(cv.median()),
    }
    for sp in SPECIES:
        sub = prot[prot["species"] == sp]["log2_ratio"].dropna()
        m[f"n_{sp}"] = len(sub)
        m[f"log2_{sp}"] = float(sub.median()) if len(sub) else np.nan
        m[f"mad_{sp}"] = float((sub - sub.median()).abs().median()) if len(sub) else np.nan
    return m


def main() -> int:
    for ds in ("dda_5min", "dda_15min"):
        print(f"\n===== {ds} (two-peptide rule, A/B headline pair) =====")
        rows = {}
        for arm, label in ARMS.items():
            d = ROOT / "results" / ds / arm
            if not (d / "lfq.tsv").is_file():
                print(f"  [skip] no lfq.tsv for {arm}")
                continue
            rows[label] = quant_metrics(d)
        out = pd.DataFrame(rows).T
        exp = {f"log2_{sp}": EXPECTED[sp] for sp in SPECIES}
        pd.set_option("display.width", 200, "display.max_columns", 50)
        print(out.round(4).to_string())
        print(f"  (expected log2 A/B: {exp})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
