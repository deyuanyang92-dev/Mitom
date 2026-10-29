# Introns: verified gold standard of 193 NCBI records and MitoM end to end (2026-10-29)

Every intron NCBI annotates in the 193 local NCBI records with introns (records.tsv: Cnidaria 172, Annelida 16,
Placozoa 3, Porifera 2; an INSDC copy of an NC_ record dropped; build_set.py) was checked, and MitoM was run on the
bare sequences (`mitom run`, current code) and scored against the checked placements.

## Gold standard (gold_README.md, introns.tsv, summary.tsv, unannotated.tsv, scripts/)

Per NCBI intron: (1) the spliced product (internal stops, start, stop, /translation); (2) the junction against 12
homolog proteins of other species (indel within +-5 aa); (3) splice signals - type by Rfam / group I subtype models,
the sliding window of placements with the identical spliced sequence, group I (exon ends U, intron ends omega-G) and
group II (GUGYG ... AY) rules; (4) the placement in related records (one site = junctions within 3 nt).

| phylum | A | B | consensus | D-corrected | D | E | in gold |
|---|---|---|---|---|---|---|---|
| Annelida | 5 | 1 | 9 | 8 | 0 | 1 | 24 |
| Cnidaria | 195 | 36 | 7 | 28 | 2 | 0 | 266 |
| Placozoa | 0 | 1 | 0 | 1 | 0 | 14 | 16 |
| Porifera | 2 | 1 | 0 | 0 | 0 | 0 | 3 |
| all | 202 | 39 | 16 | 37 | 2 | 15 | 309 |

A = NCBI's placement verified; B = same product, NCBI's placement slid off the canonical one (gold = canonical);
consensus = NCBI's product and junction clean, an alternative aligns better (gold moves only on strong evidence,
otherwise keeps NCBI's); D-corrected = NCBI's product or junction wrong (indel against the homologs, stop, a junction
side that does not align), gold = the unique correction; D = two NCBI intron features that are no introns (removed);
E = not fully checkable (placozoan trans-spliced cox1 pieces, the rRNA intron, the 1.5-kb fragment LC670763.1).
97 gold placements differ from NCBI's; NCBI's own placement equals the gold for 213 of 309. Rules derived from the
data: the hexacoral nad5 intron family ends exon ...T | intron ...A (137 of 146 records; no omega-G anywhere in the
sliding window), group II 5' end GYGYG (GCGCG in annelid cox1 introns). No NCBI record had an unannotated intron by
the criteria of unannotated.tsv. Built with the scripts in scripts/ from the local NCBI data (the gold GenBank
records, 7.8 MB, are in mitom_work/intron_gold_2026-10-29/gold/gb; MitoM's site library models/intron_sites.* is
built from them).

## MitoM end to end (e2e/)

`mitom run` on the bare records (mito_bench_run.py; 24 records at a time, 2 threads each), scored with
`benchmark/eval_introns.py` against the gold (e2e/score.txt):

| run | gold introns | exact | within 6 nt | more | missed | false introns (76 intron-free records) | NCBI ERROR |
|---|---|---|---|---|---|---|---|
| default (whole site library), 194 records | 309 | **309** | 0 | 0 | 0 | 0 | 0 |
| leave-one-out (the record left out of the core's references and of the site library: `MITOM_EXCLUDE_REF`), 270 records | 309 | **299** | 4 | 0 | 6 | 0 | 0 |

- Leave-one-out by class: A 202/202, B 39/39, D-corrected 37/37, consensus 13/16, E 8/15. Records not used while
  developing the intron code (155): 232 of 235 exact. Off: 3 hexacoral nad5 introns where the gold keeps NCBI's
  off-consensus placement (relatives give the consensus), 7 placozoan class E introns no other record carries.
- Leave-species-out / leave-genus-out (intron step on the same pipeline input, not end to end): 296 / 4 / 3 / 6.
- Two introns NCBI does not annotate are predicted, both probably real: KT989324.1 (Glycera unicornis) cox1
  6253..8609 - NCBI has only gene 5553..9433 + misc_feature 5553..6252 - and Trichoplax rrnL (default mode only; the
  site comes from Polyplacotoma).
- Against NCBI's raw annotation (311 introns): default 218 exact / 76 within 6 nt / 15 more / 2 missed (the two
  removed NCBI features); leave-one-out 208 / 80 / 15 / 8.
- Run time of the intron step: de novo (mito_intron.py) 10-20 s per record, site transfer (mito_intron_ref.py)
  2-5 s.
