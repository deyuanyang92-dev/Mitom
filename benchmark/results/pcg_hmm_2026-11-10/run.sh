#!/bin/bash
# Protein genes the annotation missed, by profile HMM of the nearest relatives (mito_pcg_hmm.py). Machine cody,
# mitotools environment. W = work data (not in git); the downloads are those of ../orf_2026-11-09/ (RefSeq) and
# ../intron_annelid_2026-11-09/ (NCBI annelids).
set -e
H=$(cd "$(dirname "$0")" && pwd)
M=$(cd "$H/../../.." && pwd)
W=${W:-/media/cody/18T_R1/yzz/_realtest/orf_2026-11-09}
P=${P:-python3}
$P $M/scripts/build_pcg_refs.py $W/refseq $W/annelid                       # -> models/pcg_refs.faa.gz
# leave-genus-out: each gene removed from 189 records (15 per phylum + 86 Sabellida with atp8) and searched again;
# atp8 searched in 200 flatworm / nematode records without an annotated atp8
$P $H/eval_pcg_hmm.py $W/refseq $W/annelid $H/leave_genus_out.tsv --jobs 40
# atp8 / atp6 by real phylum (up to 15 records per phylum with the gene), leave-genus-out
$P $H/eval_atp.py $W/refseq $W/annelid $H/atp_by_phylum.tsv --jobs 40
# the 335 NCBI annelid records without annotation, annotated by `mitom run` (../intron_annelid_2026-11-09/)
mkdir -p $W/pcg335 && cp $W/eval/ann335_in/*.gb $W/pcg335/
ls $W/pcg335/*.gb | xargs -P 40 -I{} $P $M/mito_pcg_hmm.py search -i {} --report {}.pcg.tsv
cat $W/pcg335/*.pcg.tsv | awk 'NR==1 || !/^record\t/' > $H/annelid_unannotated335_pcg_hmm.tsv
