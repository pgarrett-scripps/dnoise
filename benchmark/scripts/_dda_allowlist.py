#!/usr/bin/env python3
"""Build a DDA-restricted FASTA for the DIA-NN search (a speed optimization).

Library-free DIA-NN over the full 31k-protein hybrid proteome is slow on CPU
(~14 min/file). Restricting the predicted library to proteins the *original*
(undenoised) DDA arm detected at a lenient protein FDR shrinks the search space
~3-5x with little practical loss: the restricted DB is applied IDENTICALLY to
every DIA arm, so the cross-arm "does denoise help?" comparison stays valid, and
a lenient FDR (default 10%) is inclusive enough that truly-present proteins are
rarely excluded. Contaminants are always kept regardless of detection.

Caveat (state it in the paper): this scopes the DIA result to the DDA-detectable
proteome, so it cannot measure denoise revealing proteins DDA never saw at all.
Use the full-proteome FASTA for a confirmatory final run if needed.

The allowlist is built from the DDA arm only (denoise-independent). Accessions
are the middle '|' field of Sage's protein IDs (HUMAN_sp|ACC|ENTRY) and of the
stripped DIA-NN FASTA headers (>sp|ACC|ENTRY), so they join directly.

Usage:
  uv run scripts/_dda_allowlist.py \
      --in-fasta data/fasta/hybrid_diann.fasta \
      --out-fasta data/fasta/hybrid_diann_dda10.fasta \
      --fdr 0.10 \
      --dda results/dda_15min/original/results.sage.tsv [more ...]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def accession(protein_id: str) -> str:
    """Middle '|' field if present (UniProt accession), else the whole id."""
    parts = protein_id.split("|")
    return parts[1] if len(parts) >= 2 else protein_id


def dda_allowlist(sage_tsvs: list[Path], fdr: float) -> set[str]:
    accs: set[str] = set()
    for tsv in sage_tsvs:
        df = pd.read_csv(tsv, sep="\t")
        tgt = df[(df["label"] == 1) & (df["protein_q"] <= fdr)]
        for proteins in tgt["proteins"].dropna():
            for pid in str(proteins).split(";"):
                pid = pid.strip()
                if not pid or pid.startswith("rev_"):
                    continue
                accs.add(accession(pid))
    return accs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-fasta", required=True, type=Path)
    ap.add_argument("--out-fasta", required=True, type=Path)
    ap.add_argument("--fdr", type=float, default=0.10)
    ap.add_argument("--dda", required=True, nargs="+", type=Path,
                    help="one or more original-arm Sage results.sage.tsv")
    a = ap.parse_args()

    allow = dda_allowlist(a.dda, a.fdr)
    print(f"DDA allowlist @ protein_q<={a.fdr}: {len(allow)} accessions "
          f"from {len(a.dda)} run(s)")

    kept = total = kept_cont = 0
    keep = False
    with a.in_fasta.open() as fin, a.out_fasta.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                total += 1
                head = line[1:]
                acc = accession(head.split()[0])
                # Entry name is the 3rd '|' field; contaminants carry a Cont_ acc.
                is_cont = acc.startswith("Cont_") or "Cont_" in head.split()[0]
                keep = (acc in allow) or is_cont
                if keep:
                    kept += 1
                    kept_cont += int(is_cont)
            if keep:
                fout.write(line)
    print(f"wrote {kept}/{total} proteins to {a.out_fasta} "
          f"({kept_cont} contaminants force-kept)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
