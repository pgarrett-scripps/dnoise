#!/usr/bin/env python3
"""Build the species-tagged hybrid FASTA for the PXD070049 HYE benchmark.

Downloads the authors' own search database from PRIDE
(`uniprotkb_proteome_HYE_UniversalContaminants.fasta`: Human + S. cerevisiae +
E. coli + universal contaminants) and prepends a species token to every header
so the analysis can split proteins by species. The token is derived from the
UniProt entry-name suffix (`..._HUMAN` / `_YEAST` / `_ECOLI`); everything else
(BOVIN, MOUSE, etc.) is tagged CONT. Sage adds decoys (`rev_`) itself.

Using the authors' DB (not a fresh UniProt pull) keeps us aligned with the
reference search and fixes the prior wrong-yeast-species bug.
"""

from __future__ import annotations

import re
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "fasta" / "raw"
OUT = ROOT / "data" / "fasta" / "hybrid.fasta"

SRC_NAME = "uniprotkb_proteome_HYE_UniversalContaminants.fasta"
SRC_URL = f"https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/02/PXD070049/{SRC_NAME}"

# UniProt entry-name organism code (header `>db|ACC|MNEMONIC_CODE`) -> species tag.
CODE_TO_TAG = {"HUMAN": "HUMAN", "YEAST": "YEAST", "ECOLI": "ECOLI"}
HEADER_CODE = re.compile(r"^>\S*_([A-Za-z0-9]+)\b")


def species_tag(header: str) -> str:
    """HUMAN/YEAST/ECOLI from the entry-name suffix; CONT for anything else."""
    m = HEADER_CODE.match(header.split()[0])
    return CODE_TO_TAG.get(m.group(1), "CONT") if m else "CONT"


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    src = RAW / SRC_NAME
    if not (src.exists() and src.stat().st_size > 0):
        print(f"downloading {SRC_NAME} ...")
        urllib.request.urlretrieve(SRC_URL, src)

    counts: dict[str, int] = {}
    with src.open() as fh, OUT.open("w") as out:
        for line in fh:
            if line.startswith(">"):
                tag = species_tag(line)
                counts[tag] = counts.get(tag, 0) + 1
                out.write(f">{tag}_{line[1:]}")
            elif line.strip():
                out.write(line)

    print("\nentries per species tag:")
    for tag in ("HUMAN", "YEAST", "ECOLI", "CONT"):
        if tag in counts:
            print(f"  {tag}: {counts[tag]}")
    print(f"\nWrote {OUT} ({OUT.stat().st_size/1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
