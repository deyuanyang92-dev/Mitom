# Intron gold standard: 193 NCBI records, 311 NCBI introns (2026-10-29)

Every intron NCBI annotates in the 193 local records with introns (`../ncbi_gb`, `../records.tsv`) was checked and
given a class; the gold placement is what MitoM's intron annotation has to reproduce. Local data only.

## Files

| file | content |
|---|---|
| `introns.tsv` | one row per NCBI intron: accession, gene, strand, `ncbi` and `gold` coordinates (1-based; `excluded` / `removed`), `class`, `verdict` (the A-D result of an E intron), `type`, the result of each check (`product`, `junction`, `ncbi_signals`, `window`, `best_alt`, `clade`, `homologs`) and the `reason` |
| `summary.tsv` | phylum x class; introns in the gold, removed, moved from NCBI's placement |
| `unannotated.tsv` | possible introns NCBI did not annotate in these records (criterion 5) |
| `gb/<acc>.gb` | the NCBI record with the gold introns (see "Scoring") |
| `selftest/` | `benchmark/eval_introns.py --gb gold/gb --pred gold/gb`: 309 introns, all exact, 0 fp |
| `scripts/` | `homologs.py` (homolog table), `verify.py` (evidence), `classify.py` (classes, gold placement), `build_gb.py`, `unannotated.py` |
| `work/` | `evidence.pkl` (all evidence), `gold.json`, `homologs.pkl`, the mito_intron run used for criterion 5 |

Rebuild: `python3 scripts/homologs.py work/homologs.pkl; python3 scripts/verify.py --no_cm; python3 scripts/verify.py
--cm_only; python3 scripts/classify.py; python3 scripts/build_gb.py; python3 scripts/unannotated.py`.

## Introns checked

The truth of `MitoM/benchmark/eval_introns.py`: every `intron` feature and every gap of >= 142 nt between consecutive
same-strand parts of a CDS join (measured across the origin when the join crosses it). This includes the cis introns
of the trans-spliced placozoan cox1 genes and the one rRNA intron (NC_041549.1 rrnL). Where an intron feature and the
CDS join disagree, the join (it defines the protein) is the NCBI placement; the disagreement is listed in `reason`.

## Checks (per intron)

1. **Product** - the CDS spliced as NCBI joins it, translated with the record's code: internal stops; start codon
   (or 5' partial); stop (complete, T / TA, or 3' partial); length; agreement with `/translation`.
2. **Junction homology** - the spliced protein aligned (semi-global, BLOSUM62, gap -11/-1) to 12 proteins of the
   same gene from other species (refpool_metazoa 15,803 RefSeq + the gold set), species without an intron in that gene
   first, then by shared lineage. Junction clean = fewer than half of the homologs have an insertion / deletion
   within +-5 aa of the junction; identity in the 15 aa before / after it.
3. **Splice signals** - type by Rfam RF00028 / RF00029 / group II D1-D4 + the 14 Selles Vidal 2025 group I subtype
   models (`mito_intron.intron_types`, all introns of a record in one cmsearch call, E <= 0.01; E <= 1e-5 fixes the
   type, weaker hits allow both rules) and NCBI's note. Sliding window = all placements giving the identical spliced
   sequence. Canonical: group I = exon 1 ends T (U) and the intron ends G (omega-G); group II = intron starts GYGYG
   (GUGYG; GCGCG as in the annelid cox1 introns) and ends AY; a group II placement with the 5' motif only is used
   when no fully canonical one exists. Alternative placements: the 5' and 3' ends moved independently by <= 60 nt,
   frame kept, canonical, no internal stop; scored against the same homologs.
4. **Clade consistency** - the junction as (residue, phase) in the Metridium senile (NC_000933.1) protein; introns
   of one gene and phylum whose junctions lie within 3 nt form a site; consensus = the most common (residue, phase)
   of the site's clean members (>= 3, >= half). When most consensus members lack exon-U + omega-G (the hexacoral nad5
   intron at residue 240: 137 of 146 end ...T | intron ...A), the consensus placement is the family's rule.
5. **Unannotated introns** - see `unannotated.tsv`.

## Classes

| class | meaning | gold |
|---|---|---|
| A | product ok, junction clean, NCBI placement canonical (group I / II rule, or the family rule at the clade consensus) | NCBI placement |
| B | same product, NCBI placement slid off the canonical (or clade-consensus) one: inside the window, or a synonymous difference | the canonical placement |
| C | product ok, junction clean, no canonical or consensus placement | NCBI placement (flagged) |
| consensus | NCBI's product and junction are clean, but NCBI's placement is not canonical / not the clade consensus and exactly one canonical or consensus placement aligns better (>= 80 % of the homologs) | the alternative only on strong evidence (better for >= 11 of 12 homologs and either a score gain >= 20 or a fully canonical placement of a confidently typed intron, E <= 1e-5, where NCBI's is not canonical); otherwise NCBI's placement, the alternative listed |
| D-corrected | NCBI's own product or junction is wrong (internal stop, insertion / deletion vs >= half of the homologs at the junction, or a junction side that does not align: identity < 0.3, >= 0.5 after the correction) and exactly one corrected placement is canonical (or the clade consensus), clean, valid and better for >= 2/3 of the homologs | the correction |
| D | NCBI error without a unique correction (`excluded`), or an intron feature that is no intron (`removed`: inside an exon, or part of a join gap that is evaluated on its own row) | none |
| E | not fully checkable: rRNA intron, the 1.5-kb fragment LC670763.1, introns of the trans-spliced placozoan cox1 genes (incl. the two trans-splicing half-introns of NC_041549.1); `verdict` holds the A-D result of what could be checked | verdict's placement |

## Scoring with benchmark/eval_introns.py

`python3 MitoM/benchmark/eval_introns.py --gb gold/gb --pred <MitoM runs> --label <name> -o <out>` scores against the
gold directly. In `gold/gb` each gold intron is an `intron` feature at the gold placement whose `/note` starts with
`gold_class=<class>`, and the host CDS / rRNA join is moved to the gold exon boundaries (`/translation` redone where the
protein changes). The evaluator (extended 2026-10-29, unchanged output on NCBI records apart from the new column)
writes the class to column `truth_class` and adds one summary row per class (`gold_class A` ...). An intron of unknown
placement would be a `misc_feature` whose `/note` starts with `eval_introns: ignore`: truth and predictions >= 50 %
inside it are not scored (class `ignored`); the current gold has none. Removed intron features are simply absent.
E introns are scored like the others; drop them from a summary by their `truth_class`.

## Result

| phylum | A | B | C | consensus | D-corrected | D | E | total | in gold |
|---|---|---|---|---|---|---|---|---|---|
| Annelida | 5 | 1 | 0 | 9 | 8 | 0 | 1 | 24 | 24 |
| Cnidaria | 195 | 36 | 0 | 7 | 28 | 2 | 0 | 268 | 266 |
| Placozoa | 0 | 1 | 0 | 0 | 1 | 0 | 14 | 16 | 16 |
| Porifera | 2 | 1 | 0 | 0 | 0 | 0 | 0 | 3 | 3 |
| all | 202 | 39 | 0 | 16 | 37 | 2 | 15 | 311 | 309 |

A = 65 by the group I / II rules + 137 hexacoral nad5 introns by the family rule. consensus: gold moved for 10 (strong
evidence), NCBI's placement kept for 6 (Decemunciger nad4 x3, NC_027526.1, NC_045931.1, NC_083272.1 nad5). E verdicts:
A 5, B 1, C 4, consensus 1, D-corrected 1, and 3 not checkable (rRNA intron, 2 half-introns). Clade consensus is
computed per phylum (the hexacoral nad5 rule is not applied to Trichoplax). 91 gold placements differ from NCBI's.

## Limits

- Junction homology cannot separate a true lineage-specific indel at a junction from a misplaced boundary when no
  canonical or consensus placement removes it; the 60-nt search window bounds the corrections.
- Corrections that change the protein rest on homolog alignment; low-identity junctions (Decemunciger nad4,
  identity 0.23 / 0.44; Decemunciger nad1, 5' group II motif only) are the least certain.
- The family rule (exon-U, intron-terminal A or C at the consensus junction) is derived from the 146 consensus
  members of the hexacoral nad5 site; no other site has enough members to set a family rule.
