#!/bin/bash
# End-to-end `mitom run` vs the intron gold standard (2026-10-29), main-tree code with de novo r9 + site transfer:
#   mitom_lo     = leave-one-out: the record is left out of the core's references and of the intron-site library
#                  (MITOM_EXCLUDE_REF), on all 270 records (193 gold + 76 intron-free + KT989324.1)
#   mitom_native = pipeline default (whole library), on the 193 gold records + KT989324.1
# then eval_introns.py against gold/gb (gold records) and the NCBI records (intron-free, KT989324.1).
R=$(cd "$(dirname "$0")" && pwd)
P=/media/deyuan/18T_B/yzz/MitoM/benchmark
export OMP_NUM_THREADS=1
POS=$(awk -F'\t' 'NR>1 && $2!="neg"{print $1}' $R/groups.tsv | paste -sd,)
python3 $P/mito_bench_run.py --set $R/set.tsv --root $R --split all --tools mitom_lo --jobs 24 --threads 2 > $R/run_lo.log 2>&1
python3 $P/mito_bench_run.py --set $R/set.tsv --root $R --split all --tools mitom_native --acc $POS --jobs 24 --threads 2 > $R/run_native.log 2>&1
python3 $R/score.py > $R/score.txt 2>&1
cat $R/score.txt
