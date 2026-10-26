# rRNA ends: how close can MitoM get to NCBI? (2026-10-26)

Question (user, 2026-10-26): on the 1,050-record NC_ set MitoM's 12S / 16S have both ends exactly as in NCBI for
51.6 % of the genes; can that reach 95 %? Short answer: not against NCBI, because NCBI's annotations of closely
related genomes disagree with each other about as often. What was measured, what was fixed, and which numbers can
honestly be reported is below. Data: `benchmark/sets/mito_bench_large_2026-10.tsv` (1,050 RefSeq records, 22 phyla,
root `mitom_work/bench_large`) and `mito_bench_annelida_2026-10.tsv` (669 annelids, `mitom_work/bench_annelida`);
`mitom_lo` = leave-out mode (the record's own proteins and every benchmark record left out of the references, rRNA
HMMs rebuilt without them) is the main result, `mitom_native` = shipped references.

## 1. The ceiling: NCBI disagrees with itself

`project.py` aligns each NCBI rRNA with the rRNAs of the 5 closest annotated relatives (other species, from a
RefSeq pool without the benchmark records) and projects their NCBI ends onto the record. If annotation were
consistent, a near-identical relative would put the ends at the same place.

| relative's rRNA identity | large: both ends agree | annelids: both ends agree |
|---|---|---|
| 98-100 % | 72 % (n = 308) | 70 % (n = 63) |
| 95-97 % | 64 % (n = 194) | 44 % (n = 50) |
| 90-95 % | 60 % (n = 441) | 42 % (n = 252) |
| < 85 % | 15 % (n = 639) | 7 % (n = 687) |

(`result_identity_bins.txt`.) Examples with identical rRNA sequences (proj_species.tsv, checked in the GenBank files):
Perinereis aibuhitensis NC_023943.1 / P. linea NC_063944.1 rrnL, both annotated 11854..13050 although P. linea has a
1-nt indel upstream, so the homologous ends differ by 1 nt at both ends (coordinates copied, not boundaries);
Sternaspis scutata MK548645.1 / S. liui NC_069059.1 rrnS, 5' end 2 nt and 3' end 7 nt apart; Artemia tibetiana
NC_021383.1 / A. sorgeloosi NC_072958.1 rrnL, 5 and 3 nt apart. A consistent annotator can match only one of each pair.
Picking, per gene, whichever of the 5 relatives matches NCBI exactly (an oracle that
looks at the truth) gives 57.9 % (large) and 40.1 % (annelids) (`result_projection_consensus.txt`). Any annotator,
human or program, that is consistent across genomes therefore cannot reproduce NCBI's exact ends for much more
than ~70 % of the genes. Most annelid differences are 1-3 nt: the submitters' rRNA overlaps the flanking tRNA by
1-2 nt (MITOS style) while MitoM abuts it (tRNA punctuation), and different submitters do it differently
(`conventions.py`, `result_refseq_conventions.txt`).

## 2. Things tried that did not help (kept for the record)

| idea | script | result |
|---|---|---|
| copy the relatives' end convention where they agree | `hybrid_sim.py` | +1.3 (large) / +5.6 (annelids) points exact at best; not adopted (MitoM would copy whichever convention the nearest submitters used) |
| fill to the flanking gene only where RefSeq of that phylum is contiguous | `r1sim.py` | no gain (dev 46.0 -> 45.7 %, test 50.7 -> 50.0 %) |
| cap the fill at the model columns the span is missing (+ slack 0-100 nt), optionally extrapolate | `r1cap.py` | exact: every variant <= the current rule on dev and test (e.g. cap50: dev 45.8 -> 44.7 %, test 50.8 -> 50.4 %); consensus ends within 3 nt mixed (cap50: dev 93.7 -> 93.8 %, test 93.9 -> 95.0 %) |
| extend ends without a neighbour within 300 nt by the missing model columns | `r1cap.py` (current+free_x) | dev 45.8 -> 44.9 %, test 50.8 -> 49.8 % |
| fill only to a tRNA (tRNA punctuation, Ojala et al. 1981), optionally only a tRNA on the rRNA's strand, + rRNA / protein neighbours, +x | `r1cap.py` (fill:...) | exact: dev 45.8 -> 43.4-45.3 %, test 50.8 -> 48.8-49.9 %, annelids 32.4 -> 32.1-33.1 %; consensus within 3: dev 93.7 -> 92.8-93.6 %, test 93.9 -> 93.3-94.2 %, annelids 95.7 -> 96.2 % at best |
| where NCBI puts an end, by nearest neighbour and gap (`result_endfill_*.txt`) | `endfill.py` | same-strand tRNA within 10 nt: NCBI at the tRNA 86 %; same-strand protein gene 31-100 nt away: at the gene 20 %, at the model end 20 %, elsewhere 60 % (large). Beyond small gaps NCBI follows neither convention, so no fill rule can match it much better |
| use MITOS2's covariance model when the HMM span is partial | (MITOS2 calls as proxy) | on the 569 partial spans MITOS2 is within 10 nt for 123, MitoM for 295 |
| RNA-seq as a transcript truth (5 records, ENA) | `run_one.sh`, `frags_from_bam.py`, `ends_from_rnaseq.py` | ordinary RNA-seq has no nucleotide-resolution ends (`result_rnaseq_summary.txt`) |

## 3. Fixes made

v26 (all records rerun):
1. `mito_rrna_hmm.py build`: 17 reference rRNAs are on the wrong strand in GenBank (Chaetognatha 5 of 7 per gene;
   the model matched the antisense and every chaetognath rRNA was called on the wrong strand); every reference is
   now checked with the MITOS2 covariance model and used reverse-complemented when only the other strand hits.
   References are sorted before sampling, so the library can be rebuilt identically.
2. `mito_core.annotate_rrnas`: an HMM span covering < 40 % of the model falls back to the covariance model
   (tunicates NC_009833.1, NC_013561.1).
3. `mito_curate.py` R1 handles rRNAs across the origin (echinoderm rrnL before cox1, NC_013881.1).

v27: duplicated rRNAs. The core kept only the best span per gene; `best_spans` already found a runner-up over the
same model columns elsewhere. `sim_copies.py` (leave-out libraries, all 1,719 records): 7 such copies, every one of
them NCBI's second copy (NC_044480.1 rrnS + rrnL, NC_008046.1 rrnS, NC_008047.1 rrnS + rrnL, NC_081951.1 rrnS +
rrnL), none elsewhere (`result_copies.txt`). A copy scoring >= 0.9 x the gene (COPY_FRAC, as for protein genes) is
now annotated by the core. The rRNA rescue step (`mito_rrna_hmm.py search`) still only reports a second copy
(docs/reviews/2026-10-21/A_fixes.md): its example, region a of the split Trichoplax 16S (NC_008151.3, 8766..9793),
scores 0.24 x region b with the shipped Placozoa model, so the core does not add it either.

v27 (2): an rRNA cut by the end of a linear sequence. The rescue step rejected every span over < 40 % of the model,
but at the end of a partial genome the gene runs on beyond the sequence, so the missing columns are expected.
`mito_rrna_hmm._cut_by_sequence_end`: a linear record's span that reaches the sequence end (within 30 nt) on every
side where model columns are missing is kept and marked partial there, and R1 keeps that < / > (it rebuilt the
location without it). Scans of every run's reports find one such span in all 1,719 records: DQ517436.1 (Scoloplos,
partial genome starting inside rrnL; NCBI <1..479, MitoZ and MITOS2 annotate it too): native now <1..479; the
leave-out span lacks 33 columns at its inner 3' end (> 30), so leave-out still misses it (the threshold was not
moved for one record). Only DQ517436.1 was rerun; no large-set or anchor record has a span the rule could change.

Not changed: split 16S of oysters (NC_016180.1, NC_022688.1, NC_036478.1: 5' and 3' pieces ~2.5-4.5 kb apart, only
the 3' piece is annotated); a rule that finds the 3 also has to reject 14 other far-apart fragments, too few
examples to set it. Tunicate rrnL (NC_012887.1, NC_024105.1, NC_026107.1) is missed by MitoM, MitoZ and MITOS2.

## 4. NCBI annotations that are wrong

The full MITOS2 covariance models (cmsearch default pipeline; independent of MitoM's HMMs) hit only the strand /
gene / position MitoM reports (`ncbi_errors_cm.py`, `result_ncbi_errors_cm.txt`):

| record | NCBI | covariance model (bits) |
|---|---|---|
| NC_028407.1 Cypridopsis vidua | rrnS on + | rrnS - 275, + 0 |
| NC_028085.1 Sinergasilus polycolpus | rrnL, rrnS on + | rrnL - 203, rrnS - 110, + 0 |
| NC_032061.1 Armillifer agkistrodontis | rrnL on + | - 362, + 0 |
| NC_071876.1 Ctenodiscus crispatus | rrnL on + | - 790, + 0 |
| NC_006162.1 Graptacme eborea | rrnS on + | - 347, + 0 |
| NC_005840.1 Siphonodentalium lobatum | rrnS on + | - 230, + 0 |
| NC_067052.1 Xenos yangi | /gene="rrnL" /product="12S ribosomal RNA" at 13985..14719 | rrnS 146, rrnL 0 |
| NC_028708.1 Granata lyrata | both rRNAs "16S ribosomal RNA"; the one at 12329..13506 | rrnS 415, rrnL 0 |
| NC_027943.1 (Mollusca) | rrnS 11235..12137, rrnL 13351..14667 | rrnS at 12838..13745 (422), rrnL at 14914..16240 (704) |
| NC_044703.1 Rhinebothrium reydai | /gene="rrnL" /product="12S ribosomal RNA" at 9610..10339 | rrnS 255, rrnL 0 |
| NC_068907.1 Sternaspis sendalli (annelid set) | both rRNAs /gene="rrn12"; the one at 12273..13590 (/product 16S) | rrnL 749, rrnS 0 |
| PP998669.1 Odontosyllis enopla (annelid set) | rrnL 8182..9304, rrnS 9386..9984 on + | rrnL - 628, rrnS - 366, + 0 |
| OQ807204.1 Teuthidodrilus samae (annelid set) | rrnS 15689..16071,1..432 on - (every other gene on +) | + 440, - 0 |

These 13 records (16 genes: 10 on the wrong strand, 4 with the wrong name, 2 in the wrong place; 10 records / 12
genes in the large set, 3 / 4 in the annelid set) count as MitoM errors in every table below (the truth is not
edited). MitoM (v27, leave-out and native) reports all 16 as the covariance model does.

## 5. Results (v27)

v27 = v26 + the two rules of section 3. `final_metrics.py` on `eval/all/features.tsv` (mito_bench_eval.py
--table2asn, 6 tools; report_v27.md in each root) with proj_species.tsv; runs in
mitom_work/bench_{large,annelida,anchors}/runs/mitom_{native,lo}, v26 kept as *_v26. NC_081951.1 (large) fails for
every tool in the benchmark harness (lift to the NCBI origin: its 1.6-kb duplicated region has no unique match) and
is left out everywhere. NCBI is the truth throughout; only the "NCBI errors scored by the covariance model" row scores
the 16 genes of section 4 against the covariance model instead (`ncbi_fixed_det.py`).

| leave-out (main result) | large, all 1,049 records MitoM finished | large, 918 records all four tools finished | annelids, 669 (all tools) |
|---|---|---|---|
| rRNA detected (right gene and strand, >= 50 % of the NCBI rRNA) | 98.0 % (2,052 / 2,093; v26 97.8 %) | 98.9 % (1,817 / 1,837) | 99.6 % (1,323 / 1,328) |
| same, NCBI errors scored by the covariance model | 98.6 % | 99.5 % | 99.9 % |
| consensus-supported ends within 3 nt | 94.0 % (n = 2,481) | 96.0 % (n = 2,311) | 95.7 % (n = 1,196) |
| consensus-supported ends exact | 84.8 % | 86.9 % | 66.1 % |
| both ends within 3 nt of NCBI | 65.3 % | 70.1 % | 70.1 % |
| both ends exact as NCBI | 48.3 % | 52.4 % | 32.5 % |

Native (shipped models): large, all records 98.5 % detected (99.1 % with the covariance model for NCBI's errors),
consensus ends within 3 nt 93.5 %, both exact 48.5 %; large, common records 98.9 %, 95.8 %, 52.6 %; annelids 99.7 %
(100.0 %), 95.7 %, 32.7 %. Same records, other tools (large common / annelids): MitoZ detected 96.7 / 99.7 %,
consensus ends within 3 nt 76.5 / 26.2 %, both exact 21.0 / 10.7 %; MITOS2 96.3 / 97.9 %, 64.8 / 27.3 %, 18.3 / 12.7 %.

Unchanged against v26: protein genes and tRNAs (every number of report_v27.md), rRNA false positives (large 25
leave-out / 23 native, annelids 4). Undetected v26 -> v27: large leave-out 46 -> 41 and native 36 -> 31 (copies of
NC_044480.1 x 2, NC_008046.1, NC_008047.1 x 2), annelid native 5 -> 4 (DQ517436.1).

Which numbers reach 95 %: detection everywhere (98.0-99.7 %); consensus-supported ends within 3 nt on the annelids
(95.7 %) and on the records every tool finished (96.0 %), not on all large-set records (94.0 %). The records MitoZ
cannot run are hard for MitoM too: Cnidaria (R1 fills to protein genes where NCBI leaves gaps: 42 of the 148 ends
missed by > 3 nt) and, in leave-out mode, phyla left without any reference, whose rRNAs fall back to MITOS2's
covariance model and are missed (Ctenophora 8 genes, Chaetognatha 2, Rotifera 1, Porifera 1; the shipped models
find all 12). Both ends exact as NCBI cannot reach 95 % for any annotator (section 1).

The 41 large-set rRNAs leave-out does not detect: 12 NCBI errors (section 4); the 12 lost to leave-out above;
tunicate rrnL / rrnS x 3 and Cnidaria NC_065161.1 (2.9-kb NCBI rrnL), missed by every tool; the oyster split 16S x 3;
divergent small phyla (Kinorhyncha 2, Xenacoelomorpha 2, Gnathostomulida, Tardigrada, Acanthocephala NC_022921.1,
Nematomorpha NC_085569.1 whose NCBI rrnL spans 12,037 nt); NC_015091.1 (rrnL region partly duplicated, every tool
places it where NCBI does not); NC_027510.1 (span without the first 238 model columns: 43 % of NCBI's 1,194-nt rrnS).
Annelids: the 4 NCBI errors and DQ517436.1 (section 3).

Experimental ends (section 6): v27 = v26. 16 'paper' ends: MitoM native 12 exact / 13 within 3 nt, leave-out 12 / 12,
NCBI 14 / 15, MitoZ 2 / 7, MITOS2 2 / 7.

Run time: report_v27 medians (leave-out, large: 107 s per record, v26 87 s) come from a busier machine (load 130-170
on 72 cores: the annelid runs were started beside the large ones); v27 adds no search (the copy is best_spans'
runner-up, computed before; the edge rule only changes a verdict).

## 6. Experimental ends (literature)

`literature_rrna_ends.md` surveys papers that mapped mature rRNA ends (S1 / primer extension / RACE / ESTs /
small-RNA-seq); `rrna_end_anchors.tsv` lists them, `benchmark/sets/mito_bench_anchors_2026-10.tsv` the 13 genomes,
`anchor_score.py` the 32 ends used ('paper' = position from the paper or the authors' clones; 'derived' = abstract
statement converted with the GenBank record, equal to NCBI by construction).

## Scripts

| script | what |
|---|---|
| `project.py` | project relatives' NCBI rRNA ends onto each record (-> proj_species.tsv) |
| `identity_bins.py`, `summ.py`, `consensus.py` | agreement by identity, projection summaries, consensus-supported ends |
| `final_metrics.py features.tsv proj_species.tsv TOOL...` | detection, exact / within 3 / within 10, consensus ends |
| `anchor_score.py ROOT TOOL...` | score against the literature ends |
| `r1sim.py`, `r1cap.py` | replay R1 (rRNA end rule) variants on finished runs (`result_r1cap_*.txt`, v26 leave-out runs) |
| `endfill.py` | where NCBI puts rRNA ends relative to the nearest neighbour |
| `ncbi_errors_cm.py` | full MITOS2 covariance-model check of the NCBI rRNAs MitoM contradicts |
| `hybrid_sim.py`, `conventions.py`, `endctx.py`, `deltas.py` | convention analyses |
| `sim_copies.py` | runner-up rRNA copies vs NCBI duplicates |
| `build_lib.py` | leave-out rRNA HMM libraries |
| `run_one.sh`, `frags_from_bam.py`, `ends_from_rnaseq.py`, `ena_find.py` | RNA-seq pilot |
