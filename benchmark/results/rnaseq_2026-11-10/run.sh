#!/bin/bash
# Transcript evidence (mito_rna.py) on the Polydora cf. ciliata RNA-seq of Selifanova et al. 2023 (SRR26830859,
# paired 61 bp, random-primed cDNA). Machine cody. Downloaded with the user's organise-sra (kingfisher, md5-checked).
set -e
M=$(cd "$(dirname "$0")/../../.." && pwd)
D=${D:-/media/cody/18T_R1/yzz/_realtest/rnaseq_2026-11-10}
organise-sra download polydora_rnaseq -r SRR26830859 -o $D/polydora_rnaseq_sra -j 1
# the reference: MitoM's own annotation of OQ078742 (mitom run, ../orf_2026-11-09/), so coordinates agree
G=$D/map/OQ078742.mitom.gb
python3 $M/mito_rna.py -g $G -1 $D/polydora_rnaseq_sra/fastq/SRR26830859_1.fastq.gz \
    -2 $D/polydora_rnaseq_sra/fastq/SRR26830859_2.fastq.gz -o $D/rna_out --threads 32 --bin <minimap2 dir>
# (here: minimap2 2.28 -ax sr -Y -t 32 to mito.bam once, then `samtools view -h mito.bam | mito_rna.py --sam -`)
