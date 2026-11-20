# Short divergent exons (issue #15), 2026-11-20, cody

Question / 问题: can a homolog vote against a shuffled null (mito_intron_exon.py) recover the 19 bryozoan introns beside short exons that Exonerate does not align?

Data: 134 NCBI records (root `_realtest/intron_ncbi_2026-11-09`), scorer `benchmark/eval_introns.py`.

| run | exact | ≤6 nt | >6 nt | missed | fp | note |
|---|---|---|---|---|---|---|
| head (no search) | 148 | 47 | 31 | 22 | 3 | |
| short1/short2 (search on every chain) | 146 | 45 | 32 | 25 | 3 | broke cox1/cox2 chains that were right (false end exons); new exons ~40 nt off, internal stop |
| short3 (only genes with no model) | 148 | 47 | 31 | 22 | 3 | identical to head: no gain, no loss |

Conclusion / 结论: the locus search finds the right places (nad4, cob, cox2 of MW592987.1) but the exon boundaries and frame of the new exons are not refined, so the models fail the gates; extending existing chains adds false exons. Kept as `--short_exons` (off by default). Open: boundary refinement for the new exons; accept an extension only when it rates better than the original model.
Note: with near-identical homologs the shuffled null is weak (all copies hit the same random locus); real homologs diverge. Unit test uses divergent homologs.
