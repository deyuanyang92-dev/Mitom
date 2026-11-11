#!/bin/bash
# mito_rnaseq_annot.py on every pair of pairs/datasets.tsv (reads as in run_rna.sh). Re-runnable: done pairs skipped.
set -u
M=$(cd "$(dirname "$0")/../../.." && pwd)
D=${D:-/media/cody/18T_R1/yzz/_realtest/annelid_rna_2026-11-10}
BIN=${BIN:-$HOME/miniforge3/envs/mapbench/bin}
cd $D; mkdir -p rnaseq_annot
awk -F'\t' 'NR>1{print $1"\t"$8}' pairs/datasets.tsv | while IFS=$'\t' read mito run; do
  o=rnaseq_annot/${mito%.*}; [ -s $o/rnaseq_evidence.tsv ] && continue
  r1=""; for d in rna_local rna_sra/fastq; do
    for x in $d/${run}_1.fastq.gz $d/${run}.fastq.gz; do [ -s $x ] && r1=$x && break 2; done; done
  [ -z "$r1" ] && continue
  r2=${r1/_1./_2.}
  echo "python3 $M/mito_rnaseq_annot.py -g mitom_out/annotation/$mito.gb -1 $r1 $([ "$r2" != "$r1" ] && [ -s $r2 ] && echo -2 $r2) -o $o --run $run --threads 16 --bin $BIN > $o.log 2>&1"
done > rnaseq_annot/jobs.sh
echo "[run_rnaseq_annot] $(wc -l < rnaseq_annot/jobs.sh) pairs"
xargs -P 4 -I{} bash -c "{}" < rnaseq_annot/jobs.sh
echo "[run_rnaseq_annot] finished"
