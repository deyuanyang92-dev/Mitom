#!/bin/bash
# NCBI animal mitogenomes with introns not in the intron gold standard (2026-10-29): set, `mitom run`, scoring.
#   mitom_native = `mitom run` defaults on the bare record (source feature only; clade and code of the record)
#   mitom_lo     = the same with the record left out of the core's references and of the intron-site library
# Machine cody, MitoM f04d839 + the benchmark layout fix of this commit, mitotools environment.
set -e
H=$(cd "$(dirname "$0")" && pwd)
M=$(cd "$H/../../.." && pwd)
W=${W:-/media/cody/18T_R1/yzz/_realtest/intron_ncbi_2026-11-09}
P=${P:-python3}
export OMP_NUM_THREADS=1
$P $H/make_set.py $W --fetch
for t in mitom_native mitom_lo; do
  $P $M/benchmark/mito_bench_run.py --set $W/root/set.tsv --root $W/root --split all --tools $t --jobs 24 --threads 2
  $P $M/benchmark/eval_introns.py --gb $W/root/gb --pred $W/root/runs/$t --label $t -o $H
done
# the spliced products of both placements against the refseq89m proteins of the gene (blastp)
$P $H/check_products.py $W/root mitom_native ~/mitos_refs/refseq89m/featureProt > $H/products_native.tsv
