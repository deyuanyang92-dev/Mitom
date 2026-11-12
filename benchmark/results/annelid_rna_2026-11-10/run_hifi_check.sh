#!/bin/bash
# RNA-vs-genome indel sites (mito_rnaseq_annot step 4) checked in same-specimen HiFi reads, 5 records (hifi_site_check.py).
set -u
B=$(cd "$(dirname "$0")" && pwd)
D=${D:-/media/cody/18T_R1/yzz/_realtest/annelid_rna_2026-11-10}
W=${W:-/media/cody/18T_R1/yzz/_realtest/hifi_check_2026-11-12}
BIN=${BIN:-$HOME/miniforge3/envs/mapbench/bin}
cd $W
for x in OZ026631:ERR12760787 OZ024941:ERR12319376 OZ483948:ERR14835533 OW387151:ERR8575393 OZ461320:ERR13725949; do
  a=${x%:*}; hifi=${x#*:}; o=$W/$a; mkdir -p $o
  run=$(awk -F'\t' -v m="$a.1" '$1==m{print $8}' $B/datasets.tsv)
  r1=$(ls $D/rna_local/${run}_1.fastq.gz $D/rna_sra/fastq/${run}_1.fastq.gz 2>/dev/null | head -1)
  python3 -c "
from Bio import SeqIO
r=SeqIO.read('$D/mitom_out/annotation/$a.1.gb','genbank'); r.description=''; SeqIO.write(r,'$o/ref.fa','fasta')"
  $BIN/bwa index $o/ref.fa 2>/dev/null
  [ -s $o/rna.bam ] || $BIN/bwa mem -t 32 $o/ref.fa <(zcat $r1 | head -n 16000000) <(zcat ${r1/_1./_2.} | head -n 16000000) \
     2>/dev/null | $BIN/samtools view -b -F 0x904 -q 20 - | $BIN/samtools sort -o $o/rna.bam - 2>/dev/null
  $BIN/samtools index $o/rna.bam; $BIN/samtools fasta $o/rna.bam 2>/dev/null > $o/rna.fa
  [ -s $o/hifi.fa ] || $BIN/minimap2 -ax map-hifi -t 32 --secondary=no $o/ref.fa $W/fastq/$hifi.fastq.gz 2>/dev/null \
     | $BIN/samtools fasta -F 0x904 - 2>/dev/null > $o/hifi.fa
  echo "== $a (RNA $run, HiFi $hifi): $(grep -c '>' $o/hifi.fa) HiFi reads mapped"
  # the HiFi reads' edit rate: ~0.1 % for true mitochondrial reads; ~10 % = nuclear copies (NUMTs), sites not judged
  $BIN/minimap2 -ax map-hifi -t 16 --secondary=no $o/ref.fa $o/hifi.fa 2>/dev/null | $BIN/samtools view -F 0x904 - \
    | python3 -c "
import sys,re
e=a=0
for l in sys.stdin:
    c=l.split('\t'); nm=[t for t in c[11:] if t.startswith('NM:i:')]
    if nm: e+=int(nm[0][5:]); a+=sum(int(k) for k,o in re.findall(r'(\d+)([MI=X])',c[5]))
print(f'   HiFi edit rate {e/max(a,1):.4f}')"
  python3 $B/hifi_site_check.py -g $D/mitom_out/annotation/$a.1.gb --evidence $D/rnaseq_annot/$a/rnaseq_evidence.tsv \
     --rna-bam $o/rna.bam --rna-fa $o/rna.fa --hifi $o/hifi.fa --samtools $BIN/samtools --ref $o/ref.fa -o $o/sites.tsv
done
