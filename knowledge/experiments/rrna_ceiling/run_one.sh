#!/usr/bin/env bash
# run_one.sh RECORD.gb RUN_ACC OUTDIR [NREADS=2000000] [THREADS=8]
# Streams the first NREADS read pairs of an ENA run, trims adapters (no quality trimming), maps end-to-end to the
# circularised mitogenome with bowtie2, writes fragments (frags.tsv) and the rRNA end table (ends.tsv).
set -euo pipefail
gb=$1; run=$2; out=$3; nreads=${4:-2000000}; thr=${5:-8}
here="$(cd "$(dirname "$0")" && pwd)"
E=/home/deyuan/miniforge3/envs
mkdir -p "$out"; cd "$out"
if [ ! -s r1.fq.gz ] || [ ! -s r2.fq.gz ]; then
  urls=$(curl -s -m 120 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$run&result=read_run&fields=fastq_ftp&format=tsv" | tail -1 | cut -f2)
  u1=$(echo "$urls" | tr ';' '\n' | grep -E '_1\.fastq\.gz$' | head -1)
  u2=$(echo "$urls" | tr ';' '\n' | grep -E '_2\.fastq\.gz$' | head -1)
  [ -n "$u1" ] && [ -n "$u2" ] || { echo "no paired fastq for $run" >&2; exit 3; }
  lines=$((nreads * 4))
  ( curl --noproxy "*" -s -m 3600 "http://$u1" | gzip -dc 2>/dev/null | head -n $lines | gzip -1 > r1.fq.gz ) &
  ( curl --noproxy "*" -s -m 3600 "http://$u2" | gzip -dc 2>/dev/null | head -n $lines | gzip -1 > r2.fq.gz ) &
  wait
fi
python3 - "$gb" ref.fa <<'EOF'
import sys
from Bio import SeqIO
rec = SeqIO.read(sys.argv[1], "genbank")
s = str(rec.seq).upper()
open(sys.argv[2], "w").write(f">ref\n{s + s[:1500]}\n")
open("ref.len", "w").write(str(len(s)))
EOF
$E/fastp/bin/fastp -i r1.fq.gz -I r2.fq.gz -o t1.fq.gz -O t2.fq.gz --detect_adapter_for_pe -Q -L --length_required 25 \
  -w "$thr" -j fastp.json -h /dev/null 2> fastp.log
$E/mitgard/bin/bowtie2-build -q ref.fa ref > /dev/null
$E/mitgard/bin/bowtie2 --end-to-end --very-sensitive -X 2000 --no-mixed --no-discordant -p "$thr" -x ref \
  -1 t1.fq.gz -2 t2.fq.gz 2> bowtie2.log | $E/mitgard/bin/samtools view -b -o aln.bam -
# strandedness: fraction of read-1 alignments antisense to the annotated rRNAs is decided in python
$E/mitgard/bin/samtools view -f 2 -F 0x904 -q 10 aln.bam > aln.sam
python3 - "$gb" <<'EOF' > strand.txt
import sys, re
from Bio import SeqIO
rec = SeqIO.read(sys.argv[1], "genbank")
n = len(rec)
rr = [(int(f.location.start), int(f.location.end), f.location.strand) for f in rec.features if f.type == "rRNA"]
same = anti = 0
for line in open("aln.sam"):
    f = line.split("\t", 10)
    flag = int(f[1])
    if not flag & 0x40:
        continue
    pos = (int(f[3]) - 1) % n
    for s, e, st in rr:
        if s <= pos < e:
            rev = bool(flag & 0x10)
            if (st == 1) == (not rev):
                same += 1
            else:
                anti += 1
tot = same + anti
mode = "."
if tot >= 1000:
    if anti / tot > 0.9: mode = "R"
    elif same / tot > 0.9: mode = "F"
print(mode, same, anti)
EOF
mode=$(cut -d' ' -f1 strand.txt)
python3 "$here/frags_from_bam.py" "$mode" < aln.sam > frags.tsv
python3 "$here/ends_from_rnaseq.py" "$gb" frags.tsv ends.tsv > ends.log
rm -f aln.sam t1.fq.gz t2.fq.gz ref.*.bt2
echo "$run mode=$(cat strand.txt) fragments=$(wc -l < frags.tsv) $(grep -m1 'overall alignment rate' bowtie2.log)"
