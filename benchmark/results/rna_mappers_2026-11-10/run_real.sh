#!/bin/bash
# Mappers on real RNA-seq, first 4 M pairs: Polydora cf. ciliata SRR26830859 (61 bp, random primers) on MitoM's
# OQ078742; Erpobdella octoculata SRR10997416 (cDNA) on MitoM's OZ022248 (the same species' long-read assembly).
set -e
B=$(cd "$(dirname "$0")" && pwd); M=$(cd $B/../../.. && pwd)
W=${W:-/media/cody/18T_R1/yzz/_realtest/rna_mappers_2026-11-10}
cd $W; out=$B/real.tsv; : > $out.tmp; first=1
for x in polydora:OQ078742 erpobdella:OZ022248; do
  s=${x%:*}; g=${x#*:}
  python3 $B/run_mappers.py -r $g.fa -1 real/${s}_1.fq.gz -2 real/${s}_2.fq.gz -o realmaps/$s -t 8 --drop_unmapped
  python3 $B/evaluate_real.py --mitom $M -g $g.gb --maps realmaps/$s --label $s $([ $first = 1 ] && echo --header) >> $out.tmp
  first=0
done
mv $out.tmp $out
