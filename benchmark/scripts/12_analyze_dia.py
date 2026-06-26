#!/usr/bin/env python3
"""Analyze DIA-NN results for a diaPASEF DATASET: IDs + LFQ vs the HYE truth.

Mirrors 06_analyze.py (Sage) but reads DIA-NN output and writes the SAME
results/<DATASET>/analysis/ schema (summary.csv, accuracy.csv, the four plots)
so scripts/10_compare.py picks it up unchanged.

Arms compared, in order (only those with DIA-NN output on disk are included):
  original  -> results/<DATASET>/original
  denoised  -> results/<DATASET>/denoised      (MS1-only dnoise)
  msms      -> results/<DATASET>/msms           (MS1 + whole-frame MS2 dnoise)

Inputs per arm (written by scripts/11_diann.sh with --matrices):
  report.pg_matrix.tsv  -- protein groups x runs, values = PG.MaxLFQ
  report.pr_matrix.tsv  -- precursors x runs (for ID counts)

SPECIES TAGGING: DIA-NN searches a stripped FASTA (standard >sp| headers), so the
SPECIES_ prefix the Sage path uses is gone. We rebuild an accession->species map
from the ORIGINAL prefixed hybrid.fasta (the authoritative tag) and map DIA-NN's
Protein.Group accessions back through it -- contaminants (CONT_) map to None and
are excluded, exactly like the Sage species_of rule.

Ground truth log2(A/B): HUMAN 0, YEAST +1, ECOLI -2 (COND_FRAC in _metrics.py).
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Keep all lettering >= 4.5 pt once the figure is downscaled to the text width.
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "figure.titlesize": 15,
})
import pandas as pd
from matplotlib.patches import Patch

from _metrics import (
    EXPECTED,
    MIN_PEPTIDES_PER_PROTEIN,
    PAIRS,
    SPECIES,
    condition_of,
    expected_log2,
    median_ci,
    norm_factors,
)

ROOT = Path(__file__).resolve().parents[1]
DATASET = os.environ.get("DATASET", "dia_5min")
RESULTS = ROOT / "results" / DATASET
OUT = RESULTS / "analysis"
FASTA = ROOT / "data" / "fasta" / "hybrid.fasta"  # the prefixed (Sage) FASTA

CANDIDATE_ARMS = ["original", "denoised", "msms"]
# Colorblind-safe (Wong, Nat. Methods 2011); avoids the red/green pairing.
ARM_COLOR = {"original": "#0072B2", "denoised": "#E69F00", "msms": "#009E73"}


# ---------- species map (accession -> HUMAN/YEAST/ECOLI, else None) ----------

def build_acc2species(fasta: Path) -> dict[str, str]:
    acc2sp: dict[str, str] = {}
    with fasta.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            head = line[1:]
            sp = head.split("_", 1)[0]
            parts = head.split("|")
            if len(parts) >= 2:
                acc2sp[parts[1]] = sp  # parts[1] = UniProt accession
    return acc2sp


def species_of_pg(protein_group: str, acc2sp: dict[str, str]) -> str | None:
    """Species iff every accession in the group maps to the same HYE species.

    Mirrors _metrics.species_of: mixed / contaminant / unknown -> None.
    Handles isoform suffixes (P12345-2 -> P12345).
    """
    tags = set()
    for acc in str(protein_group).split(";"):
        acc = acc.strip()
        if not acc:
            continue
        sp = acc2sp.get(acc) or acc2sp.get(acc.split("-", 1)[0])
        tags.add(sp)
    if len(tags) == 1 and (t := next(iter(tags))) in SPECIES:
        return t
    return None


# ---------- matrix loading + rollup ----------

def load_matrix(path: Path) -> tuple[pd.DataFrame, list[str]]:
    """Load a DIA-NN matrix. Run columns are the per-file quant columns, detected
    by the `.d` path suffix (robust to DIA-NN version differences in the leading
    metadata columns, e.g. 2.2.0 adds N.Sequences / drops Protein.Ids in pg)."""
    df = pd.read_csv(path, sep="\t")
    run_cols = [c for c in df.columns if c.endswith(".d")]
    return df, run_cols


def pg_peptides(arm_dir: Path) -> pd.DataFrame | None:
    """Distinct (protein group, peptide-sequence) pairs from the precursor matrix,
    for the two-peptide rule and the quantified-peptide count. DIA-NN matrices are
    already FDR-filtered, so every listed precursor is a quantified peptide ion."""
    pr = arm_dir / "report.pr_matrix.tsv"
    if not pr.is_file():
        return None
    df, _ = load_matrix(pr)
    if "Protein.Group" not in df.columns:
        return None
    seqcol = "Stripped.Sequence" if "Stripped.Sequence" in df.columns else "Modified.Sequence"
    sub = df[["Protein.Group", seqcol]].dropna()
    sub.columns = ["protein", "peptide"]
    return sub.astype(str).drop_duplicates()


def pg_table(arm_dir: Path, acc2sp: dict[str, str]) -> pd.DataFrame | None:
    """Per-protein-group, per-run MaxLFQ intensities, total-intensity normalized,
    tagged by species. One row per protein group (DIA-NN already rolled up).
    Protein groups are restricted to those supported by
    >= MIN_PEPTIDES_PER_PROTEIN distinct peptide sequences (two-peptide rule)."""
    p = arm_dir / "report.pg_matrix.tsv"
    if not p.is_file():
        return None
    df, run_cols = load_matrix(p)
    if not run_cols or df.empty:
        return None
    factors = norm_factors(df, run_cols)
    vals = df[run_cols].replace(0, np.nan).mul(pd.Series(factors), axis=1)
    out = pd.DataFrame({"protein": df["Protein.Group"].astype(str)})
    out["species"] = out["protein"].map(lambda g: species_of_pg(g, acc2sp))
    out = pd.concat([out, vals], axis=1)
    if MIN_PEPTIDES_PER_PROTEIN > 1:
        peps = pg_peptides(arm_dir)
        if peps is not None:
            counts = peps.groupby("protein")["peptide"].nunique()
            keep = set(counts[counts >= MIN_PEPTIDES_PER_PROTEIN].index)
            out = out[out["protein"].isin(keep)]
    return out


def quant_peptide_count(arm_dir: Path, p: pd.DataFrame | None) -> int:
    """Distinct quantified peptide sequences mapping to the quantified protein
    groups in `p` (the set reported by `lfq_frame` for this arm)."""
    if p is None or p.empty:
        return 0
    peps = pg_peptides(arm_dir)
    if peps is None:
        return 0
    return int(peps[peps["protein"].isin(set(p["protein"]))]["peptide"].nunique())


def cond_cols(run_cols: list[str], cond: str) -> list[str]:
    return [c for c in run_cols if condition_of(c) == cond]


# ---------- identifications ----------

def id_metrics(arm_dir: Path, acc2sp: dict[str, str]) -> dict:
    """Experiment-wide IDs from the FDR-filtered DIA-NN matrices: precursors,
    peptides (stripped sequences), and protein groups, plus per-species peptides."""
    m: dict = {"n_precursor": 0, "n_peptide": 0, "n_protein": 0}
    pr = arm_dir / "report.pr_matrix.tsv"
    if pr.is_file():
        df, _ = load_matrix(pr)
        m["n_precursor"] = int(df["Precursor.Id"].nunique()) if "Precursor.Id" in df else len(df)
        seqcol = "Stripped.Sequence" if "Stripped.Sequence" in df else "Modified.Sequence"
        pep = df.drop_duplicates(seqcol)
        m["n_peptide"] = int(pep[seqcol].nunique())
        pep = pep.assign(species=pep["Protein.Group"].map(lambda g: species_of_pg(g, acc2sp)))
        for sp in SPECIES:
            m[f"n_peptide_{sp}"] = int((pep["species"] == sp).sum())
    pg = arm_dir / "report.pg_matrix.tsv"
    if pg.is_file():
        df, _ = load_matrix(pg)
        m["n_protein"] = int(df["Protein.Group"].nunique())
    return m


# ---------- LFQ accuracy / precision ----------

def lfq_frame(prot: pd.DataFrame | None) -> pd.DataFrame | None:
    """Per-protein A-vs-B log2 ratio + within-condition CV (>=2/3 reps each)."""
    if prot is None or prot.empty:
        return None
    run_cols = [c for c in prot.columns if c not in ("protein", "species")]
    a_cols, b_cols = cond_cols(run_cols, "A"), cond_cols(run_cols, "B")
    if not a_cols or not b_cols:
        return None
    a_ok = prot[a_cols].notna().sum(axis=1) >= 2
    b_ok = prot[b_cols].notna().sum(axis=1) >= 2
    p = prot[a_ok & b_ok].copy()
    p["mean_A"] = p[a_cols].mean(axis=1)
    p["mean_B"] = p[b_cols].mean(axis=1)
    p["log2_ratio"] = np.log2(p["mean_A"] / p["mean_B"])
    p["cv_A"] = p[a_cols].std(axis=1, ddof=1) / p[a_cols].mean(axis=1)
    p["cv_B"] = p[b_cols].std(axis=1, ddof=1) / p[b_cols].mean(axis=1)
    return p


def lfq_metrics(p: pd.DataFrame | None, arm_dir: Path | None = None) -> dict:
    if p is None or p.empty:
        return {"n_quantified": 0, "n_quant_peptide": 0}
    m = {"n_quantified": len(p)}
    if arm_dir is not None:
        m["n_quant_peptide"] = quant_peptide_count(arm_dir, p)
    cv = pd.concat([p["cv_A"], p["cv_B"]]).dropna()
    m["median_cv"] = float(cv.median())
    for sp in SPECIES:
        sub = p[p["species"] == sp]["log2_ratio"].dropna()
        m[f"n_{sp}"] = len(sub)
        m[f"median_log2_{sp}"] = float(sub.median()) if len(sub) else np.nan
        m[f"bias_{sp}"] = (m[f"median_log2_{sp}"] - EXPECTED[sp]) if len(sub) else np.nan
        m[f"mad_{sp}"] = float((sub - sub.median()).abs().median()) if len(sub) else np.nan
    return m


def pair_accuracy(arm: str, prot: pd.DataFrame | None) -> list[dict]:
    """Observed median log2 ratio per (pair, species) for A/B, A/C, B/C."""
    if prot is None or prot.empty:
        return []
    run_cols = [c for c in prot.columns if c not in ("protein", "species")]
    means = {}
    for c in "ABC":
        cc = cond_cols(run_cols, c)
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
            lo, hi = median_ci(sub.values)
            out.append({"arm": arm, "pair": f"{a}/{b}", "species": sp,
                        "expected": expected_log2(a, b, sp),
                        "observed": float(sub.median()), "ci_lo": lo, "ci_hi": hi,
                        "n": len(sub)})
    return out


# ---------- plots (same filenames/shape as 06_analyze.py) ----------

def plot_ids(summary: pd.DataFrame, arms: list[str]) -> None:
    metrics = [m for m in ("n_precursor", "n_peptide", "n_protein", "n_quantified")
               if m in summary.columns]
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(metrics))
    n = len(arms)
    w = 0.8 / max(n, 1)
    for i, arm in enumerate(arms):
        vals = [summary.loc[arm, m] for m in metrics]
        bars = ax.bar(x + (-(n - 1) / 2 + i) * w, vals, w, label=arm, color=ARM_COLOR.get(arm))
        ax.bar_label(bars, fmt="%d", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.set_ylabel("count @ 1% FDR")
    ax.set_title(f"DIA-NN identifications & quantified proteins ({DATASET})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "id_counts.png", dpi=150)
    plt.close(fig)


def plot_violins(frames: dict[str, pd.DataFrame | None], arms: list[str]) -> None:
    n = len(arms)
    width = 0.8 / max(n, 1)
    offsets = [(-(n - 1) / 2 + j) * width for j in range(n)]
    fig, ax = plt.subplots(figsize=(8, 5))
    for j, arm in enumerate(arms):
        p = frames[arm]
        if p is None or p.empty:
            continue
        data, pos = [], []
        for i, sp in enumerate(SPECIES):
            vals = p[p["species"] == sp]["log2_ratio"].dropna().values
            if len(vals) >= 2:
                data.append(vals)
                pos.append(i + offsets[j])
        if not data:
            continue
        vp = ax.violinplot(data, positions=pos, widths=width * 0.9,
                           showmedians=True, showextrema=False)
        for body in vp["bodies"]:
            body.set_facecolor(ARM_COLOR[arm])
            body.set_alpha(0.55)
            body.set_edgecolor("black")
            body.set_linewidth(0.5)
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(1.2)
    for i, sp in enumerate(SPECIES):
        ax.hlines(EXPECTED[sp], i - 0.4, i + 0.4, color="gray", lw=2, ls="--", zorder=5)
    ax.set_xticks(range(len(SPECIES)))
    ax.set_xticklabels(SPECIES)
    ax.set_ylabel("log2(A/B) protein ratio")
    ax.set_ylim(-5, 4)
    ax.set_title(f"DIA LFQ ratio distributions (dashed = expected, {DATASET})")
    ax.legend(handles=[Patch(facecolor=ARM_COLOR[a], alpha=0.55, edgecolor="black", label=a)
                       for a in arms])
    fig.tight_layout()
    fig.savefig(OUT / "lfq_ratio_violins.png", dpi=150)
    plt.close(fig)


def plot_cv(frames: dict[str, pd.DataFrame | None], arms: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    for arm in arms:
        p = frames[arm]
        if p is None or p.empty:
            continue
        cv = pd.concat([p["cv_A"], p["cv_B"]]).dropna()
        cv = cv[cv < 1.0]
        ax.hist(cv, bins=40, histtype="step", lw=2, label=f"{arm} (median {cv.median():.3f})")
    ax.set_xlabel("protein CV within condition")
    ax.set_ylabel("count")
    ax.set_title(f"DIA LFQ precision ({DATASET})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "lfq_cv.png", dpi=150)
    plt.close(fig)


def plot_accuracy(acc: pd.DataFrame, arms: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    lim = (-3.5, 4.0)
    ax.plot(lim, lim, color="gray", ls="--", lw=1, zorder=0, label="ideal")
    marker = {"HUMAN": "o", "YEAST": "s", "ECOLI": "^"}
    for arm in arms:
        d = acc[acc["arm"] == arm]
        for sp in SPECIES:
            ds = d[d["species"] == sp]
            ax.scatter(ds["expected"], ds["observed"], color=ARM_COLOR.get(arm),
                       marker=marker[sp], s=45, alpha=0.8,
                       label=arm if sp == "HUMAN" else None)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("expected log2 ratio")
    ax.set_ylabel("observed median log2 ratio")
    ax.set_title(f"DIA LFQ accuracy across the dynamic range ({DATASET})\n"
                 "(A/B, A/C, B/C; o human ^ ecoli s yeast)")
    ax.legend(title="arm")
    fig.tight_layout()
    fig.savefig(OUT / "lfq_accuracy.png", dpi=150)
    plt.close(fig)


def main() -> int:
    arms = [a for a in CANDIDATE_ARMS if (RESULTS / a / "report.pg_matrix.tsv").is_file()]
    if not arms:
        print(f"no DIA-NN output found under {RESULTS}/(original|denoised|msms)/ "
              "-- run scripts/11_diann.sh first")
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    acc2sp = build_acc2species(FASTA)

    rows, frames, acc = {}, {}, []
    for arm in arms:
        prot = pg_table(RESULTS / arm, acc2sp)
        frames[arm] = lfq_frame(prot)
        rows[arm] = {**id_metrics(RESULTS / arm, acc2sp), **lfq_metrics(frames[arm], RESULTS / arm)}
        acc += pair_accuracy(arm, prot)

    summary = pd.DataFrame(rows).T.reindex(sorted({k for r in rows.values() for k in r}), axis=1)
    summary.to_csv(OUT / "summary.csv")
    acc_df = pd.DataFrame(acc)
    acc_df.to_csv(OUT / "accuracy.csv", index=False)

    plot_ids(summary, arms)
    plot_violins(frames, arms)
    plot_cv(frames, arms)
    if not acc_df.empty:
        plot_accuracy(acc_df, arms)

    pd.set_option("display.width", 200, "display.max_columns", 100)
    print(summary.T)
    if not acc_df.empty:
        print("\nDIA LFQ accuracy across pairs (observed vs expected median log2):")
        print(acc_df.pivot_table(index=["pair", "species"], columns="arm",
                                 values="observed").round(2))
    print(f"\nWrote {OUT}/summary.csv, accuracy.csv, and plots.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
