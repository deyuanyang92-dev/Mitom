#!/bin/bash
# every case through MitoM gb2tbl and through NCBI's gbf2tbl + gbf2fsa, each validated by table2asn (2026-10-30)
R=$(cd "$(dirname "$0")" && pwd)
M=/media/deyuan/18T_B/yzz/MitoM
B=/home/deyuan/miniforge3/envs/busco/bin
T2A=/home/deyuan/miniforge3/envs/mitotools/bin/table2asn
rm -rf $R/mitom $R/ncbi; mkdir -p $R/mitom $R/ncbi
for f in $R/in/*.gb; do
  c=$(basename $f .gb)
  # MitoM: one file -> one pair named after the file, validated
  $M/mitom gb2tbl -i $f -o $R/mitom/$c --validate > $R/mitom/$c.log 2>&1; echo "mitom $c exit=$?" >> $R/mitom/exits.txt
  # NCBI EDirect: gbf2tbl + gbf2fsa, then table2asn on the pair
  mkdir -p $R/ncbi/$c
  ( PATH=$B:$PATH gbf2tbl < $f > $R/ncbi/$c/$c.tbl 2> $R/ncbi/$c/gbf2tbl.err; echo "gbf2tbl exit=$?" > $R/ncbi/$c/exits.txt
    PATH=$B:$PATH gbf2fsa < $f > $R/ncbi/$c/$c.fsa 2> $R/ncbi/$c/gbf2fsa.err; echo "gbf2fsa exit=$?" >> $R/ncbi/$c/exits.txt )
  if [ -s $R/ncbi/$c/$c.fsa ]; then $T2A -i $R/ncbi/$c/$c.fsa -V vb -o $R/ncbi/$c/$c.sqn > $R/ncbi/$c/t2a.log 2>&1; fi
done
echo done
