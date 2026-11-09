#!/bin/bash
# ORFs of unknown function (mito_orf.py): library, the three test sets, MITOS2 on three Polydorini. Machine cody,
# mitotools environment. W = work data (not in git).
set -e
M=$(cd "$(dirname "$0")/../../.." && pwd)
W=${W:-/media/cody/18T_R1/yzz/_realtest/orf_2026-11-09}
P=${P:-python3}
$P $M/scripts/build_orf_library.py $W/refseq                       # 16,073 RefSeq records -> models/orf_library.fna.gz
# 1. Spionidae: the 46 NCBI records with >= 10 CDS (esearch Spionidae[Organism] AND mitochondrion[filter], 10-60 kb)
$P $M/mito_orf.py -i $W/eval/spionidae.gb -o $W/eval/spio_out --threads 16
# 2. 362 intron-free RefSeq records of 26 phyla (../intron_algo_2026-11-09/negatives_set.tsv; in the library, but
#    each record's own entries are left out)
$P $M/mito_orf.py -i $W/../intron_algo_2026-11-09/neg/gb -o $W/eval/neg_out --threads 32
# 3. 335 NCBI annelid records without any annotation, annotated by `mitom run` (../intron_annelid_2026-11-09/)
$P $M/mito_orf.py -i $W/eval/ann335_in -o $W/eval/ann335_out --threads 32
# 4. MITOS2 and MitoM on three Polydorini; NC_061377 alone with a relative given (--ref = mitom run --orf_ref)
$P $M/benchmark/mito_bench_run.py --set $W/mitos2_spio/set.tsv --root $W/mitos2_spio --split all --tools mitos2,mitom_native --jobs 6 --threads 2
$P $M/mito_orf.py -i $W/mitos2_spio/runs/mitom_native/NC_061377.1/raw/annotation/NC_061377.1.gb -o $W/mitos2_spio/nc061377_ref --ref $W/mitos2_spio/gb/OQ078742.1.gb
