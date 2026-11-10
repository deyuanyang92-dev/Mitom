#!/bin/bash
# mitom rna on the annelid genome-project mitogenomes (long-read assemblies) with the same species' RNA-seq.
# pairs/datasets.tsv (select_pairs.py): mitogenome + one RNA-seq run of the same genome project (ncbi-genome-inventory
# links: same specimen > same project, other specimen > same project, specimen unknown); every accession recorded.
# Reads: LAN server copies (first 8 M pairs, rna_local/) or organise-sra downloads (rna_sra/fastq/). Re-runnable:
# a species with rna_out/<mito>/rna_features.tsv is skipped; species whose reads are not here yet are listed.
set -u
M=$(cd "$(dirname "$0")/../../.." && pwd)
D=${D:-/media/cody/18T_R1/yzz/_realtest/annelid_rna_2026-11-10}
BIN=${BIN:-$HOME/miniforge3/envs/mapbench/bin}
cd $D; mkdir -p rna_out; : > rna_out/waiting.tsv
awk -F'\t' 'NR>1{print $1"\t"$8}' pairs/datasets.tsv | while IFS=$'\t' read mito run; do
  acc=${mito%.*}; o=rna_out/$acc
  [ -s $o/rna_features.tsv ] && continue
  r1=""; for d in rna_local rna_sra/fastq; do
    for x in $d/${run}_1.fastq.gz $d/${run}_1.head8M.fastq.gz $d/${run}.fastq.gz; do [ -s $x ] && r1=$x && break 2; done; done
  if [ -z "$r1" ]; then echo -e "$acc\t$run\tno reads yet" >> rna_out/waiting.tsv; continue; fi
  [ -e ${r1%/*}/.done/$run ] || [ -e ${r1%/*}/.done/$run.json ] || { echo -e "$acc\t$run\tdownloading" >> rna_out/waiting.tsv; continue; }
  r2=${r1/_1./_2.}
  echo "python3 $M/mito_rna.py -g mitom_out/annotation/$mito.gb -1 $r1 $([ -s $r2 ] && echo -2 $r2) -o $o --threads 16 --bin $BIN > $o.log 2>&1"
done > rna_out/jobs.sh
echo "[run_rna] $(wc -l < rna_out/jobs.sh) species to run, $(wc -l < rna_out/waiting.tsv) waiting for reads"
xargs -P 4 -I{} bash -c "{}" < rna_out/jobs.sh
