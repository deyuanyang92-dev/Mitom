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
