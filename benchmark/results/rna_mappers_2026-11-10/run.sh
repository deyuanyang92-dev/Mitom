#!/bin/bash
# RNA-seq mapper benchmark on mitochondrial genomes (machine cody). Truth from simulate.py; mappers in conda env
# mapbench (mamba create -n mapbench -c conda-forge -c bioconda minimap2 bwa bwa-mem2 bowtie2 hisat2 star bbmap
# subread segemehl samtools; with ALL_PROXY=socks://... set, libmamba's TLS fails - unset ALL_PROXY / all_proxy).
set -e
B=$(cd "$(dirname "$0")" && pwd); M=$(cd $B/../../.. && pwd)
W=${W:-/media/cody/18T_R1/yzz/_realtest/rna_mappers_2026-11-10}
A=/media/cody/18T_R1/yzz/_realtest/annelid_rna_2026-11-10/mitom_out/annotation
cd $W
cp -n /media/cody/18T_R1/yzz/_realtest/rnaseq_2026-11-10/map/OQ078742.mitom.gb OQ078742.gb
cp -n $A/OZ022248.1.gb OZ022248.gb; cp -n $A/OZ414225.1.gb OZ414225.gb
out=$B/simulated.tsv; first=1; : > $out.tmp
for g in OQ078742 OZ022248 OZ414225; do
  python3 -c "from Bio import SeqIO; r=SeqIO.read('$g.gb','genbank'); r.description=''; SeqIO.write(r,'$g.fa','fasta')"
  for div in 0 0.01 0.03; do for len in 61 150; do
    s=sim/${g}_d${div}_l${len}; m=maps/${g}_d${div}_l${len}
    [ -s ${s}.genes.tsv ] || python3 $B/simulate.py -g $g.gb -o $s -n 200000 -l $len --div $div --seed 7
    python3 $B/run_mappers.py -r $g.fa -1 ${s}_1.fq.gz -2 ${s}_2.fq.gz -o $m -t 8
    python3 $B/evaluate.py --mitom $M -g $g.gb --sim $s --maps $m --label "${g}_d${div}_l${len}" \
      $([ $first = 1 ] && echo --header) >> $out.tmp
    first=0
  done; done
done
mv $out.tmp $out
