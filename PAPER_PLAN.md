# TODO — before JPR submission

## ⚠️ BLOCKER: DIA 15-min searched against the WRONG library (not yet fixed)
The two diaPASEF gradients used DIFFERENT DIA-NN libraries, so their absolute
counts are NOT comparable across gradients:
  - dia_5min  -> hybrid_diann_dda10  (restricted DDA-10%-FDR allowlist):
                 9,305 isoforms / 12,321 groups / 1.77M precursors  (432 MB speclib)
  - dia_15min -> hybrid_diann        (FULL ~31k proteome, no allowlist):
                 31,390 isoforms / 41,867 groups / 4.84M precursors (1.18 GB speclib)
Cause: /tmp/dia15_run.sh never set DIANN_FASTA, so 11_diann.sh defaulted to the
full FASTA. Within-gradient cross-arm (denoise) comparison is still VALID; only
the 5-min-vs-15-min absolute counts are on different bases.
Affected if not fixed: @fig:dia coverage panel (5min restricted vs 15min full),
its caption line "Counts are within the DDA-detectable proteome the search was
restricted to" (true for 5min, FALSE for 15min), and SI §S7 two-gradient prose.

FIX (recommended): re-run the 18 dia_15min runs x3 arms against the restricted
library so both gradients match and the existing SI rationale holds. Everything is
on disk:
    cd benchmark
    DIANN_FASTA=$PWD/data/fasta/hybrid_diann_dda10.fasta \
      DATASET=dia_15min just diann          # (or scripts/11_diann.sh)
    DATASET=dia_15min just analyze-dia       # 12_analyze_dia.py -> summary/accuracy.csv
    uv run scripts/22_paper_figures.py       # regen fig4_dia.png
Then finish SI §S7: per_run_dia_15min.typ table is ALREADY generated; expand the
DIA summary table (tab:dia-5min) to both gradients and write the IDs/LFQ prose.
(~5-8h compute; restricted lib is ~2.7x smaller than the full one already run.)

NOTE: S7 reduction prose + per-run tables (both gradients) are already updated and
correct — the data-volume numbers are library-independent. Only the DIA
IDs/quantification (summary table + prose) is blocked on the re-run above.

## Other pre-submission blockers (from review 2026-06-12)
- Front matter missing (JPR/ACS): co-authors, affiliation (Dept [TODO]), ORCID,
  Author Contributions, Funding, Conflict of Interest, Acknowledgements,
  Associated Content heading, designated TOC graphic.
- Citations: DIA-NN never cited (add Demichev 2020); timsrust + PXD070049 carry
  TODO notes in references.bib; intro "existing approaches" cites nothing; only 5
  refs total (thin for a full article).
- Bruker Minesweeper comparison: confirm permission/attribution (footnote flag).
- Two [TODO: exact CPU model] placeholders (paper Performance + SI S8 Hardware).
- 15-min CV framing: 0.059 -> 0.073 called "essentially unchanged" — soften.

---

overall the paper looks good.

figure 1 is great! fantastic visualization

figure2 also is good, but I want it to shwo compression for 5/15/DDA/DIA so eventaull it wuill show 4 different file/experimnt types.

figure 3 is ok.... but pretty weak. Will thik about how to make this better.

figure 4 shows identifications psm/peptide/protein/quantified. Eventaully i want this to shwo data for the 4 differetn file/experiemtns. Similar to figure 2. 

figure 5 - shows the lfq accuracy. this is good but should eventaully show for all 4 file/experiemnt types.

figure 6 - add in ms1 watershed centroider algorithm. This shoudl only be for dda 5 min. and should show the compression, IDs, and LFQ asccuracy. this one will only be run for the small testing ddatset.