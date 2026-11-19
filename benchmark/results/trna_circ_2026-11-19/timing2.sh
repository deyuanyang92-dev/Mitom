#!/bin/bash
# usage: timing2.sh REPO OUTDIR records.gb... : mito_core annotate per record, MiTFi always vs skip-when-unused;
# prints record, mode, tRNA-step seconds (from the log timestamps), whole-core seconds, skipped?, tRNA names
REPO=$1; OUT=$2; shift 2
export PATH=/home/cody/miniforge3/envs/mitotools/bin:/home/cody/mybiosofts/NCBI/blast/bin:$PATH
one() {
  f=$1; b=$(basename $f .gb); mkdir -p $OUT
  read code lin phy topo < <(python -c "
from Bio import SeqIO; r=SeqIO.read('$f','genbank'); t=r.annotations.get('taxonomy',[])
c=[f.qualifiers['transl_table'][0] for f in r.features if f.type=='CDS' and 'transl_table' in f.qualifiers] or ['5']
open('$OUT/$b.fa','w').write('>'+r.id+'\n'+str(r.seq)+'\n'); print(c[0], ';'.join(t).replace(' ',''), t[2] if t[1]=='Metazoa' and len(t)>2 else 'x', r.annotations.get('topology'))")
  for mode in always skip; do
    o=$OUT/$b.$mode; rm -rf $o; mkdir -p $o
    if [ $mode = always ]; then export MITOM_TRNA_ALWAYS_MITFI=1; else unset MITOM_TRNA_ALWAYS_MITFI; fi
    s=$(date +%s.%N)
    python -W ignore $REPO/mito_core.py annotate -i $OUT/$b.fa -o $o --genetic_code $code --threads 2 --lineage "$lin" --log $o/core.log > $o/stdout 2>&1
    e=$(date +%s.%N)
    printf "%s\t%s\t%.1f\t%s\t%s\t%s\n" $b $mode $(echo "$e - $s" | bc) "$(grep -c 'MiTFi skipped' $o/core.log $o/stdout | awk -F: '{s+=$NF} END{print s}')" "$(grep -hv '^#' $o/*.core.tsv | md5sum | cut -c1-8)" "$(python3 $OUT/../mitfi_t.py $o/core.log)"
  done
}
export -f one; export OUT REPO
printf "%s\n" "$@" | xargs -P 6 -I{} bash -c 'one {}'
