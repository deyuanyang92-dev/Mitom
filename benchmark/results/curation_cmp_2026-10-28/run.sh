#!/bin/bash
# three-way comparison 2026-10-28: MitoM (default / --curate no), MitoZ and MITOS2 raw (bench_mini runs), MitoZ / MITOS2 + MitoM curation
R=$(cd "$(dirname "$0")" && pwd)
B=/media/deyuan/18T_B/yzz/MitoM/benchmark
export OMP_NUM_THREADS=1
python3 $B/mito_bench_run.py --set $B/sets/mito_bench_mini_2026-10.tsv --root $R --split all --tools mitoz_cur,mitos2_cur --jobs 12 --threads 2 > $R/run_cur.log 2>&1
python3 $B/mito_bench_run.py --set $B/sets/mito_bench_mini_2026-10.tsv --root $R --split all --tools mitom_native,mitom_nocur --jobs 12 --threads 4 > $R/run_mitom.log 2>&1
python3 $B/mito_bench_eval.py --set $B/sets/mito_bench_mini_2026-10.tsv --root $R --split all --tools mitom_native,mitom_nocur,mitoz,mitos2,mitoz_cur,mitos2_cur --table2asn > $R/eval.log 2>&1
echo done > $R/done.txt
