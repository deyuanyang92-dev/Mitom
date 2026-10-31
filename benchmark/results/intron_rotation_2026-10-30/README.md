# Introns with the origin inside the gene: rotation test (2026-10-30)

Question: does MitoM find an intron when the record starts inside the gene that carries it (MJ20260820041 CHY
records: cox1 exon 2 at the start, exon 1 and the 1,142 bp group II intron at the end)?

Records: the 270 of intron_gold_2026-10-29 (accs.txt: 193 with 309 checked introns - Cnidaria, Annelida, Placozoa,
Porifera - and 76 without), fetched from NCBI (efetch). Every circular record with a cis-spliced intron gene (cox1
when it has one) is rotated (scripts/make_sets.py; rotations.tsv, rotations_v2.tsv):

| set | origin |
|---|---|
| orig | as in NCBI (bare FASTA, no topology label) |
| exon | middle of the gene's last exon |
| intron | middle of its last intron |
| near5 / near3 | 15 nt inside the intron's left / right genomic end (near3 = the CHY records) |

Run: benchmark/mito_bench_run.py --tools mitom_lo with MITOM_EXCLUDE_LEVEL=genus (the record's genus left out of the
intron-site library, its own proteins out of the core's references); scored with scripts/score.py (eval_introns.py
against the NCBI records, and every rotated set row by row against the unrotated one).

| origin | 171154d (before) | 738f84d | 66e0e09 | 2dae5b3 |
|---|---|---|---|---|
| exon / intron middle | 295 / 299 | 295 | 295 | 295 |
| 15 nt inside the left end | 190, 12 false introns | 295 | 295 | 295 |
| 15 nt inside the right end (CHY) | 232 | 249 | 295 | 295 |
| rows identical to the unrotated run | - | - | 300 / 300 each | 300 / 300 each |

- 738f84d: spliced re-check of a gene cut at the junction (exonerate protein2genome); window limited to ~6 kb of
  intron - the hexacoral nad5 intron (10-12 kb, holds other genes) failed (65 records linear).
- 66e0e09: the re-check over the whole rest of the circle on either side.
- 2dae5b3 (scores of 243f03a before the rebase): the part the spliced alignment adds must align alone on the far side
  (it joined contigs to chance matches), and a gene stopping short of the junction is checked too (placozoan nad5).
- The misses and > 6 nt placements are the same in every set (Trichoplax trans-spliced cox1, an Anemonia viridis
  224-nt NCBI 'intron', two Diadumene cox1 introns merged; most > 6 nt are NCBI placements the gold standard corrects).

Files: score_<commit>.txt (score_66e0e09.txt = run of 451d1a4, score_2dae5b3.txt = run of 243f03a: the same code
before the rebase onto the user's commits).
