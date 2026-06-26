#!/usr/bin/env python3
"""MS1-based vs MS2-based diaPASEF quant, per arm -> SI table fragment.

The headline DIA quant (12_analyze_dia.py, @tab:dia-2grad) uses report.pg_matrix.tsv
= PG.MaxLFQ, built from MS2 *fragment* areas. dnoise's `denoised` arm modifies ONLY
MS1 frames, so an MS2-based number is largely insensitive to it -- it cannot show
what MS1 denoising does to the quantity it actually touches. This script rolls up
the MS1 precursor quantity (Ms1.Normalised, summed per protein group) from
report.parquet and recomputes ratio/CV alongside an identically-built MS2 rollup
(PG.MaxLFQ), so the ONLY thing that differs between the two columns is the
quantity: same FDR set, same two-peptide rule, same total-intensity normalization.

DIA-NN emits no MS1 MaxLFQ, so the MS1 rollup is a plain per-protein sum; this is
why MS1 CV is higher than MS2 MaxLFQ CV (a noisier estimator, not a denoising
effect). The comparison of interest is across arms within the MS1 column.

Writes paper/si/ms1_quant_dia.typ (mirrors @tab:dia-2grad layout) and prints the
full MS1-vs-MS2 table for both gradients.

Usage: uv run --with pyarrow python scripts/54_ms1_quant_dia.py [dataset ...]
       (default: dia_5min dia_15min)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from _metrics import EXPECTED, MIN_PEPTIDES_PER_PROTEIN, SPECIES, condition_of, norm_factors

ROOT = Path(__file__).resolve().parents[1]
SI = ROOT.parent / "paper" / "si"
FASTA = ROOT / "data" / "fasta" / "hybrid.fasta"
ARMS = ["original", "denoised", "msms"]


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
                acc2sp[parts[1]] = sp
    return acc2sp


def species_of_pg(pg: str, acc2sp: dict[str, str]) -> str | None:
    tags = set()
    for acc in str(pg).split(";"):
        acc = acc.strip()
        if not acc:
            continue
        tags.add(acc2sp.get(acc) or acc2sp.get(acc.split("-", 1)[0]))
    if len(tags) == 1 and (t := next(iter(tags))) in SPECIES:
        return t
    return None


def load_fdr(arm_dir: Path) -> pd.DataFrame:
    """FDR-filtered precursor rows (1% run + PG q-value) with rollup columns."""
    cols = ["Run", "Protein.Group", "Stripped.Sequence", "Decoy",
            "Q.Value", "PG.Q.Value", "Ms1.Normalised", "PG.MaxLFQ"]
    df = pd.read_parquet(arm_dir / "report.parquet", columns=cols)
    return df[(df["Decoy"] == 0) & (df["Q.Value"] <= 0.01) & (df["PG.Q.Value"] <= 0.01)].copy()


def protein_matrix(df: pd.DataFrame, acc2sp: dict, mode: str) -> pd.DataFrame:
    """protein x run matrix. 'ms1' -> sum of precursor Ms1.Normalised;
    'ms2' -> PG.MaxLFQ (DIA-NN's per-group value). Two-peptide rule on distinct
    stripped sequences (global), identical for both modes."""
    col, agg = ("Ms1.Normalised", "sum") if mode == "ms1" else ("PG.MaxLFQ", "max")
    g = df.groupby(["Protein.Group", "Run"])[col].agg(agg).unstack("Run").replace(0, np.nan)
    if MIN_PEPTIDES_PER_PROTEIN > 1:
        npep = df.groupby("Protein.Group")["Stripped.Sequence"].nunique()
        g = g.loc[g.index.intersection(npep[npep >= MIN_PEPTIDES_PER_PROTEIN].index)]
    out = g.reset_index().rename(columns={"Protein.Group": "protein"})
    out.insert(1, "species", out["protein"].map(lambda p: species_of_pg(p, acc2sp)))
    return out


def lfq_stats(mat: pd.DataFrame) -> dict:
    run_cols = [c for c in mat.columns if c not in ("protein", "species")]
    cond = {c: condition_of(c) for c in run_cols}
    a = [c for c in run_cols if cond[c] == "A"]
    b = [c for c in run_cols if cond[c] == "B"]
    factors = norm_factors(mat, run_cols)
    v = mat[run_cols].mul(pd.Series(factors), axis=1)
    p = mat[(v[a].notna().sum(axis=1) >= 2) & (v[b].notna().sum(axis=1) >= 2)].copy()
    va, vb = v.loc[p.index, a], v.loc[p.index, b]
    ratio = np.log2(va.mean(axis=1) / vb.mean(axis=1))
    cv = pd.concat([va.std(axis=1, ddof=1) / va.mean(axis=1),
                    vb.std(axis=1, ddof=1) / vb.mean(axis=1)]).dropna()
    m = {"n_quant": len(p), "median_cv": float(cv.median())}
    for sp in SPECIES:
        sub = ratio[(p["species"] == sp).to_numpy()].dropna()
        m[f"obs_{sp}"] = float(sub.median()) if len(sub) else np.nan
    return m


def fmt(x: float) -> str:
    return f"+{x:.2f}" if x >= 0 else f"{x:.2f}".replace("-", "#sym.minus ")


def write_si(per: dict[tuple[str, str], dict], datasets: list[str]) -> None:
    """One #table mirroring @tab:dia-2grad, but for the MS1-summed quantity."""
    g5, g15 = datasets[0], datasets[1]

    def cells(metric_key):
        return [per[(g, arm)][metric_key] for g in (g5, g15) for arm in ARMS]

    rows = [
        ("Quantified proteins", [f"{int(v):,}" for v in cells("n_quant")]),
        ("Median CV", [f"{v:.3f}" for v in cells("median_cv")]),
        ("$log_2(A/B)$ #sym.space human (exp. 0)", [fmt(v) for v in cells("obs_HUMAN")]),
        ("$log_2(A/B)$ #sym.space yeast (exp. +1)", [fmt(v) for v in cells("obs_YEAST")]),
        ("$log_2(A/B)$ #sym.space E. coli (exp. #sym.minus 2)", [fmt(v) for v in cells("obs_ECOLI")]),
    ]
    lines = [
        "// AUTO-GENERATED by scripts/54_ms1_quant_dia.py -- do not edit by hand.",
        "#table(",
        "  columns: 7,",
        "  align: (left, center, center, center, center, center, center),",
        "  table.header(",
        "    table.cell(rowspan: 2)[Metric (MS1 quantity)],",
        "    table.cell(colspan: 3)[5-minute], table.cell(colspan: 3)[15-minute],",
        "    [Orig.], [MS1], [+MS/MS], [Orig.], [MS1], [+MS/MS],",
        "  ),",
    ]
    for label, vals in rows:
        lines.append(f"  [{label}], " + ", ".join(f"[{v}]" for v in vals) + ",")
    lines.append(")")
    SI.mkdir(parents=True, exist_ok=True)
    (SI / "ms1_quant_dia.typ").write_text("\n".join(lines) + "\n")
    print(f"\n  wrote {SI / 'ms1_quant_dia.typ'}")


def main() -> None:
    datasets = sys.argv[1:] or ["dia_5min", "dia_15min"]
    acc2sp = build_acc2species(FASTA)
    per: dict[tuple[str, str], dict] = {}
    print("\n=== diaPASEF MS1-sum vs MS2-MaxLFQ protein quant (2-peptide rule) ===")
    print("observed median log2(A/B); exp: human 0, yeast +1, E.coli -2\n")
    for ds in datasets:
        res_dir = ROOT / "results" / ds
        for arm in ARMS:
            d = res_dir / arm
            if not (d / "report.parquet").is_file():
                continue
            df = load_fdr(d)
            m1 = lfq_stats(protein_matrix(df, acc2sp, "ms1"))
            m2 = lfq_stats(protein_matrix(df, acc2sp, "ms2"))
            per[(ds, arm)] = m1
            print(f"{ds:10s} {arm:9s} MS2  n={m2['n_quant']:6d} cv={m2['median_cv']:.3f} "
                  f"H={m2['obs_HUMAN']:+.2f} Y={m2['obs_YEAST']:+.2f} E={m2['obs_ECOLI']:+.2f}")
            print(f"{'':10s} {'':9s} MS1  n={m1['n_quant']:6d} cv={m1['median_cv']:.3f} "
                  f"H={m1['obs_HUMAN']:+.2f} Y={m1['obs_YEAST']:+.2f} E={m1['obs_ECOLI']:+.2f}")
    if len(datasets) == 2 and all((ds, "original") in per for ds in datasets):
        write_si(per, datasets)


if __name__ == "__main__":
    main()
