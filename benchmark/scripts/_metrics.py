#!/usr/bin/env python3
"""Shared ID/LFQ/compression metric helpers for the benchmark analysis scripts.

Imported by 06_analyze.py (3-arm comparison), 07_data_reduction.py (compression),
and 13_watershed_figure.py (Figure 6). Every metric function takes an explicit
arm results directory (`.../results/<dataset>/<arm>`) or a `.d` path, so callers
control which arms they read — there is no module-level dataset state here.

Ground-truth HYE mixing fractions and expected log2 ratios live here so all
scripts agree on the same truth.
"""

from __future__ import annotations

import re
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

SPECIES = ["HUMAN", "YEAST", "ECOLI"]
EXPECTED = {"HUMAN": 0.0, "YEAST": 1.0, "ECOLI": -2.0}  # A/B (headline pair)
# HYE mixing fractions per condition (from the dataset SDRF).
COND_FRAC = {
    "A": {"HUMAN": 0.65, "YEAST": 0.30, "ECOLI": 0.05},
    "B": {"HUMAN": 0.65, "YEAST": 0.15, "ECOLI": 0.20},
    "C": {"HUMAN": 0.65, "YEAST": 0.03, "ECOLI": 0.32},
}
PAIRS = [("A", "B"), ("A", "C"), ("B", "C")]
SPECIES_COLOR = {"HUMAN": "#444444", "YEAST": "#1f77b4", "ECOLI": "#d62728"}
FDR = 0.01


def species_of(proteins: str) -> str | None:
    """Species tag iff every (target) protein for the peptide is the same species.

    Excludes decoys, contaminants, and cross-species/shared peptides.
    """
    tags = set()
    for p in str(proteins).split(";"):
        p = p.strip()
        if not p or p.startswith("rev_"):
            continue
        tag = p.split("_", 1)[0]
        tags.add(tag)
    if len(tags) == 1 and (t := next(iter(tags))) in SPECIES:
        return t
    return None


# ---------- identifications ----------

def id_metrics(arm_dir: Path) -> dict:
    """PSM/peptide/protein counts at 1% FDR for one arm's `results.sage.tsv`."""
    df = pd.read_csv(arm_dir / "results.sage.tsv", sep="\t")
    tgt = df[df["label"] == 1].copy()
    psm = tgt[tgt["spectrum_q"] <= FDR]
    pep = tgt[tgt["peptide_q"] <= FDR].drop_duplicates("peptide")
    prot = tgt[tgt["protein_q"] <= FDR]

    prot_groups = set()
    for p in prot["proteins"]:
        prot_groups.update(x for x in str(p).split(";") if not x.startswith("rev_"))

    m = {"n_psm": len(psm), "n_peptide": len(pep), "n_protein": len(prot_groups)}
    pep = pep.assign(species=pep["proteins"].map(species_of))
    for sp in SPECIES:
        m[f"n_peptide_{sp}"] = int((pep["species"] == sp).sum())
    return m


# ---------- LFQ ----------

def file_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.endswith(".d")]


def condition_of(col: str) -> str | None:
    m = re.search(r"Condition_([ABC])", col)
    return m.group(1) if m else None


def norm_factors(df: pd.DataFrame, fcols: list[str]) -> dict[str, float]:
    """Total-intensity normalization factors: scale each run so its summed
    peptide intensity equals the cross-run mean. Corrects per-run loading /
    sensitivity differences which otherwise shift every species' log2 ratio by
    a constant offset."""
    totals = df[fcols].replace(0, np.nan).sum(axis=0, min_count=1)
    target = totals.mean()
    return {c: (target / t if t and t > 0 else 1.0) for c, t in totals.items()}


def _rollup(arm_dir: Path) -> tuple[pd.DataFrame | None, list[str]]:
    """Read an arm's lfq.tsv, total-intensity-normalize, keep species-clean
    peptides, and roll up to per-(protein, species) intensities per file."""
    df = pd.read_csv(arm_dir / "lfq.tsv", sep="\t")
    if "q_value" in df.columns:
        df = df[df["q_value"] <= FDR]
    fcols = file_columns(df)
    if not fcols or df.empty:
        return None, []
    factors = norm_factors(df, fcols)
    df = df.assign(species=df["proteins"].map(species_of))
    df = df[df["species"].notna()]
    df = df.assign(protein=df["proteins"].map(lambda p: str(p).split(";")[0]))
    vals = df[fcols].replace(0, np.nan).mul(pd.Series(factors), axis=1)
    df = pd.concat([df[["protein", "species"]], vals], axis=1)
    prot = df.groupby(["protein", "species"], as_index=False)[fcols].sum(min_count=1)
    return prot, fcols


def expected_log2(a: str, b: str, sp: str) -> float:
    return float(np.log2(COND_FRAC[a][sp] / COND_FRAC[b][sp]))


def pair_accuracy(arm: str, arm_dir: Path) -> list[dict]:
    """Observed median log2 ratio per (pair, species) for A/B, A/C, B/C, vs the
    SDRF-derived expected ratios. A protein contributes to a pair if quantified
    in >=2/3 replicates of both that pair's conditions. `arm` labels the rows."""
    prot, fcols = _rollup(arm_dir)
    if prot is None:
        return []
    means = {}
    for c in "ABC":
        cc = [x for x in fcols if condition_of(x) == c]
        if cc:
            ok = prot[cc].notna().sum(axis=1) >= 2
            means[c] = prot[cc].mean(axis=1).where(ok)
    out = []
    for a, b in PAIRS:
        if a not in means or b not in means:
            continue
        ratio = np.log2(means[a] / means[b])
        for sp in SPECIES:
            sub = ratio[prot["species"] == sp].dropna()
            if len(sub) == 0:
                continue
            out.append({
                "arm": arm, "pair": f"{a}/{b}", "species": sp,
                "expected": expected_log2(a, b, sp),
                "observed": float(sub.median()), "n": len(sub),
            })
    return out


def lfq_table(arm_dir: Path) -> pd.DataFrame | None:
    """Per-(protein, species) A vs B log2 ratio and within-condition CV."""
    prot, fcols = _rollup(arm_dir)
    if prot is None:
        return None

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
    return prot


def lfq_metrics(prot: pd.DataFrame | None) -> dict:
    if prot is None or prot.empty:
        return {"n_quantified": 0}
    m = {"n_quantified": len(prot)}
    cv = pd.concat([prot["cv_A"], prot["cv_B"]]).dropna()
    m["median_cv"] = float(cv.median())
    for sp in SPECIES:
        sub = prot[prot["species"] == sp]["log2_ratio"].dropna()
        m[f"n_{sp}"] = len(sub)
        m[f"median_log2_{sp}"] = float(sub.median()) if len(sub) else np.nan
        m[f"bias_{sp}"] = m[f"median_log2_{sp}"] - EXPECTED[sp] if len(sub) else np.nan
        m[f"mad_{sp}"] = float((sub - sub.median()).abs().median()) if len(sub) else np.nan
    return m


# ---------- compression ----------

def stats(d: Path) -> tuple[int, int, int]:
    """(MS1 peaks, MS2 peaks, tdf_bin bytes) for a .d folder."""
    c = sqlite3.connect(d / "analysis.tdf")
    ms1 = c.execute("SELECT COALESCE(SUM(NumPeaks),0) FROM Frames WHERE MsMsType=0").fetchone()[0]
    ms2 = c.execute("SELECT COALESCE(SUM(NumPeaks),0) FROM Frames WHERE MsMsType!=0").fetchone()[0]
    c.close()
    return int(ms1), int(ms2), (d / "analysis.tdf_bin").stat().st_size
