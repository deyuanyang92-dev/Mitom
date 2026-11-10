#!/bin/bash
# After the downloads: map the remaining pairs, then --refine every pair (run accession in the notes) and summarise.
B=$(cd "$(dirname "$0")" && pwd); M=$(cd $B/../../.. && pwd)
D=${D:-/media/cody/18T_R1/yzz/_realtest/annelid_rna_2026-11-10}; cd $D
# (run after the downloads have finished; a pgrep wait here matched its own launcher command and never ended)
$B/run_rna.sh
awk -F'\t' 'NR>1{print $1"\t"$8}' pairs/datasets.tsv | while IFS=$'\t' read mito run; do
  o=rna_out/${mito%.*}; [ -s $o/rna_features.tsv ] || continue
  python3 $M/mito_rna.py -g mitom_out/annotation/$mito.gb -o $o --refine --run $run 2>> rna_out/refine.log
done
python3 $B/summarise.py $D -o $B/results.tsv
