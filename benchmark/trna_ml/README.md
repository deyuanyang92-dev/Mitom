# MitoM-ML tRNA annotator: research scripts and status (2026-10-26)

Research pipeline behind `trna_ml.py` (the production module). Data, scans and trained
models are not in git (several GB); `scripts/fetch_invert.py` and `scripts/build_dataset.py`
rebuild the dataset from NCBI RefSeq.

## Method
- 220 clade x 22 family Infernal CMs (`cm_build.py`), built from the train split only.
- Scan (`scan.py`): circular sequence + 300 nt, a `--mid` rescue for families with no hit at E <= 1e-3,
  and an optional third pass `rescue_max.py` (cmsearch --max with the MITOS2 refseq89m pan CMs).
- Caller (`learn.py`): HistGradientBoosting over locus x family candidates; shuffled genomes
  (`make_shuffled.py`) are negatives with weight 0.1. Boundary model: `bounds.py`.
- Baselines: `run_baselines.py` (ARWEN, ARAGORN, tRNAscan-SE, MITOS2), `run_mitfi.py` (real MiTFi through `mito_core`).
- Scoring: `evaluate.py` (same strand, >= 50 % overlap, same amino acid; L1/L2, S1/S2 checked).

## Splits (`splits/`)
- `test_accs.txt`: 1,069 records not in refseq89m, with their families held out of training.
- `excl_accs.txt`: 417 records used to train the caller.
- `hard_accs.txt`: 22 records of divergent families (Tetranychidae, Mermithidae, Lithobiidae, Steganacaridae).
- `hardx_accs.txt`: 405 genomes of 138 divergent families (extra training data, not yet used in a result).
- `mitos_done_test.txt`: 246 test records also run through MITOS2.

## Results (measured)
| set | method | recall | exact ends | FP/record |
|---|---|---|---|---|
| test 1,069 | MitoM-ML (clx) | 98.5 % | 86.6 % | 0.17 |
| test 1,069 | own emulation of the MiTFi CMs (NOT real MiTFi) | 94.9 % | - | 0.19 |
| test 246 | MITOS2 | 96.6 % | 80.7 % | 0.34 |
| test 246 | MitoM-ML (clx) | 96.6 % | 81.5 % | 0.20 |
| test 246 | ARWEN / ARAGORN | 93.3 / 94.8 % | 43 % | 3.1 / 4.4 |
| hard 22 | real MiTFi | 49.6 % | 27.8 % | 1.82 |
| hard 22 | MitoM-ML (clx) | 45.2 % | 23.8 % | 1.64 |
| 150 shuffled genomes | MitoM-ML (clx) | - | - | 0 calls |

Honest reading:
- On ordinary genomes MitoM-ML equals MITOS2 in recall; it has fewer false positives (0.20 vs 0.34).
  Of 105 truth tRNAs both miss, 91 are also missed by ARWEN and ARAGORN; 38 of the 48 MitoM-ML
  "false positives" sit where MITOS2 also calls a tRNA. The remaining errors are likely RefSeq label
  problems (not checked by hand).
- On divergent families MitoM-ML is still WORSE than real MiTFi. Not solved.

## Open work
1. `run_final.sh`: add the `--max` rescue (scan/mx_*) and train with hardx; check whether hard recall beats MiTFi.
2. If it does, add `rescue_max` to `trna_ml.scan` and regenerate `models/trna_ml/*.tsv.gz`.
3. Real MiTFi on 50 shuffled genomes and the 246 test records (`run_mitfi.py`) was still running.
4. Miss analysis on hard: 128 misses had a right-family candidate the caller rejected, 95 only wrong-family
   candidates, 43 no candidate. The caller distrusts weak hits in genomes with few strong tRNAs.

## 2026-10-27: production path tested against MITOS2 on held-out data (second machine)

Sets: `sets/mito_bench_trna_test_2026-10.tsv` = the 1,057 usable records of `splits/test_accs.txt` (families unseen
by the clade CMs and the caller, not in refseq89m), split by sha1(family) into `splits/test_tune_dev.txt` (480,
tuning) and `test_tune_test.txt` (577, held out; md5 cannot split them: the test families were chosen by md5 % 4);
`sets/mito_bench_trna_vert_2026-10.tsv` = 188 vertebrates of 102 orders not in refseq89m nor in the large set
(`splits/vert_accs.txt`; the models never saw a vertebrate); the 501 GenBank-only records of the annelid set
(MitoM-ML trained on RefSeq only). The 22 `hard_accs.txt` records are NOT clean for the shipped models (their
families are in the clade CMs: MitoM-ML 86 % recall there vs 45 % in the table above), so they are left out.
Scoring: `mito_bench_eval.py` + `scripts/trna_metrics.py` (found = amino acid, Leu / Ser class by anticodon, strand,
>= 50 % of the NCBI tRNA; precision = found / all predicted tRNAs).

Scripts: `trna_only.py` (the core's tRNA step: MitoM-ML or MiTFi), `trna_dump.py` (MitoM-ML keeping every candidate),
`rescue_sim.py`, `ensemble_sim.py`, `replay_noclade.py`, `mitfi_ends.py`, `bounds_xfit.py`, `caller_xfit.py`,
`mitfi_variants.py`, `ext_tools.py`, `compare_runs.py`.

Methods (dev / held-out half, all records; recall, precision, both ends exact):

| method | recall | precision | both ends exact |
|---|---|---|---|
| ARWEN | 94.6 / 96.0 % | 85.5 / 86.5 % | 41.5 / 42.3 % |
| ARAGORN | 95.7 / 97.0 % | 80.7 / 82.0 % | 41.9 / 42.5 % |
| tRNAscan-SE 2 -O | 47.2 / 54.5 % | 98.8 / 99.0 % | 41.7 / 49.3 % |
| MiTFi as mito_core runs it (MITOS2's tRNA step) | 97.39 / 98.35 % | 98.82 / 98.85 % | 82.85 / 88.55 % |
| MiTFi glocal (cmsearch -g, -method g) | 97.6 / 98.4 % | 98.9 / 98.8 % | 82.8 / 88.5 % |
| MiTFi -onlycutoff / -overlap 50 / fast pass E 0.001 | 95.1 / 97.3 %; 97.3 %; 97.3 % | 99.0 / 98.9 %; 98.8 %; 98.8 % | 81.7 / 87.9 %; 82.8 %; 82.9 % |
| MitoM-ML as shipped | 97.86 / 98.70 % | 98.97 / 98.83 % | 83.56 / 88.67 % |
| + MiTFi calls of the families it did not call | 98.45 / 98.88 % | 98.88 / 98.79 % | 83.88 / 88.80 % |
| + no duplicate of a gene at the sequence start (final) | 98.45 / 98.88 % | 98.95 / 98.90 % | 83.87 / 88.78 % |

Tried and not taken (dev decided, held-out half confirmed): caller threshold 0.3-0.9 (0.5 best F1 on dev; 0.8 adds
0.07 points precision for 0.05 recall); votes of MitoM-ML / MiTFi / ARWEN / ARAGORN (>= 2 of 4: precision 88 %;
>= 3: recall 95-97 %); ARWEN + ARAGORN agreement for the families still missing (precision -0.9 points: they agree on
many loci without an NCBI tRNA; of the 239 NCBI tRNAs the final method misses only 26 are found by either);
dropping calls no other caller supports when p < 0.9 (+0.01 points); caller or boundary model retrained with the
other half (cross-fit) and neighbour-gap features (no gain); dropping tRNAs on a same-strand protein gene by > 40 nt
(66 correct, 9 wrong on dev); one call per family (NCBI has most duplicate families too).

Vertebrates: the learned ends never saw a vertebrate and are 1 nt off NCBI more often than MiTFi's; MiTFi's ends for
the tRNAs MiTFi also calls (`mito_core.mitfi_ends`, lineage with Vertebrata): 434 large-set chordates (tuning) both
ends exact 90.2 -> 91.4 % (MITOS2 91.1 %), 188 held-out vertebrates 91.3 -> 92.7 %. Without the Chordata clade CMs
(built from 28 invertebrate chordates) the ends were not better (89.5 %).

Bugs fixed: a tRNA at the sequence start was called twice (`trna_ml.drop_wrap_duplicates`: 21 duplicates in 1,057
records, 6 in 188 vertebrates); without scikit-learn the first MitoM-ML call stopped the core run (`mitom_ready` now
checks it, `install.sh` installs it); scikit-learn's OpenMP pool took every core (99 threads per record, load 160 on
72 cores; `threadpool_limits` to `--threads`). Full MitoM runs give the same tRNAs as the replay of the dumps
(`compare_runs.py`).

MiTFi variants on the dev half (480 records): default recall 97.39 %, precision 98.82 %; fast pass E 0.1: 97.47 /
98.71 %; with the clade CMs of the record's lineage given to MiTFi too: 98.18 / 96.14 %, 0.64 false tRNAs per
record instead of 0.08 (MiTFi's top-hit-per-model selection was made for 22 models; the learned caller is what makes
the clade models useful: MitoM-ML 97.86 / 98.97 % with the same models).

Final comparison with MITOS2 (runmitos full mode, refseq89m), records every tool finished (11 MITOS2 runs over 1 h
left out): recall / precision / both ends exact / false positives per record

| set | MitoM (final) | MITOS2 | MiTFi as MitoM runs it |
|---|---|---|---|
| unseen families, all 1,046 | 98.68 / 98.94 / 86.54 % / 0.076 | 96.59 / 98.79 / 84.63 % / 0.101 | 97.91 / 98.85 / 85.95 % / 0.086 |
| tuning half (478) | 98.44 / 98.95 / 83.86 % | 96.15 / 98.74 / 81.54 % | 97.39 / 98.82 / 82.84 % |
| held-out half (568) | 98.88 / 98.93 / 88.76 % | 96.95 / 98.82 / 87.20 % | 98.35 / 98.88 / 88.54 % |

By phylum (all 1,043 records every tool finished): Arthropoda 754 records 99.13 / 99.21 / 89.35 % vs MITOS2 96.81 /
99.20 / 87.35 %; Mollusca 147 98.24 / 98.21 / 79.14 % vs 97.60 / 98.08 / 78.53 %; Cnidaria 57 92.86 / 87.50 / 81.63 %
vs 90.82 / 68.99 / 79.59 % (MITOS2 0.65 false tRNAs per record); Nematoda 23 99.41 / 99.60 / 66.40 % vs 96.05 / 99.18
/ 61.66 %; Platyhelminthes 22 99.38 / 99.38 / 83.47 % vs 94.21 / 99.56 / 77.69 %; Echinodermata 20 96.83 / 96.83 /
80.95 % vs 96.15 / 96.80 / 79.37 %; Annelida 13 99.31 / 100 / 86.85 % vs 97.92 / 100 / 86.51 %; Acanthocephala 4
43.18 / 86.36 / 14.77 % vs 45.45 / 88.89 / 13.64 % (MitoM-ML alone 17 %); Kinorhyncha 1 68.2 % vs 68.2 % recall
(MitoM-ML alone 0 %). Full MitoM runs agree with the replay of the dumps (221 of 222 records the same tRNAs; 188 of
188 vertebrates).

| set | MitoM (final) | MITOS2 |
|---|---|---|
| 188 vertebrates (held out) | 98.43 / 98.55 / 92.71 % / 0.048 | 98.31 / 98.24 / 92.61 % / 0.080 |
| 501 GenBank-only annelids (species unseen by MitoM-ML) | 99.35 / 98.95 / 83.52 % / 0.154 | 97.61 / 99.01 / 81.47 % / 0.132 |

The annelids are the one set where MITOS2's precision is higher (by 0.06 points: 11 more calls without an NCBI tRNA
in 501 records); 58 of MitoM's 77 such calls are at loci where MITOS2 calls a tRNA too, and several of the rest lie
where NCBI's own protein genes run over a tRNA (OP537514.1 nad1 over trnI, MG428625.1 nd2 over trnS).

