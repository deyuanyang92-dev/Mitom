#!/bin/bash
# ORFs of unknown function (mito_orf.py): library, truth, the test sets, MITOS2 on three Polydorini. Machine cody,
# mitotools environment. W = work data (not in git).
set -e
H=$(cd "$(dirname "$0")" && pwd)
M=$(cd "$H/../../.." && pwd)
W=${W:-/media/cody/18T_R1/yzz/_realtest/orf_2026-11-09}
NEG=${NEG:-$W/../intron_algo_2026-11-09/neg/gb}       # 362 intron-free RefSeq records (build_negatives.py)
P=${P:-python3}
$P $M/scripts/build_orf_library.py $W/refseq                       # 16,073 RefSeq records -> models/orf_library.fna.gz
# 1. known ORFs: records without their non-standard CDS, each run alone (its own library entries left out)
$P $H/make_truth.py $W
mkdir -p $W/recall/out
ls $W/recall/in/*.gb | grep -v NC_137528 | xargs -P 40 -I{} sh -c "$P $M/mito_orf.py -i {} -o $W/recall/out/\$(basename {} .gb)"
$P $H/score_recall.py $W/recall out
$P $H/length_baseline.py $W $NEG > $H/length_baseline.tsv
# 2. intron-free, Spionidae (46 NCBI records with >= 10 CDS), 335 unannotated annelids annotated by `mitom run`
$P $M/mito_orf.py -i $NEG -o $W/eval/neg --threads 32
$P $M/mito_orf.py -i $W/eval/spionidae.gb -o $W/eval/spio --threads 16
$P $M/mito_orf.py -i $W/eval/ann335_in -o $W/eval/ann --threads 32
# 3. MITOS2 and MitoM on three Polydorini; NC_061377 alone and with a relative (--ref = mitom run --orf_ref)
$P $M/benchmark/mito_bench_run.py --set $W/mitos2_spio/set.tsv --root $W/mitos2_spio --split all --tools mitos2,mitom_native --jobs 6 --threads 2
G=$W/mitos2_spio/runs/mitom_native/NC_061377.1/raw/annotation/NC_061377.1.gb
$P $M/mito_orf.py -i $G -o $W/mitos2_spio/nc061377_ref --ref $W/mitos2_spio/gb/OQ078742.1.gb
