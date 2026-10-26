"""Fragments (left, right, transcript strand) from a name-sorted or coordinate-sorted SAM stream of proper pairs.

Reads only properly paired, end-to-end (no 'S'/'H' in CIGAR) primary alignments with MAPQ >= 10, takes the pair
from read 1 (left = min start, right = max end).  Transcript strand: for dUTP/TruSeq-stranded libraries read 1 is
antisense ('R' mode), for ligation libraries read 1 is sense ('F'); '.' = unstranded.  usage:
samtools view -f 2 -F 0x904 -q 10 in.bam | frags_from_bam.py MODE > frags.tsv
"""
import re
import sys

mode = sys.argv[1] if len(sys.argv) > 1 else "."
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def ref_len(cigar):
    return sum(int(n) for n, op in CIG.findall(cigar) if op in "MDN=X")


for line in sys.stdin:
    f = line.split("\t", 11)
    flag = int(f[1])
    if not flag & 0x40:            # read 1 only
        continue
    cig = f[5]
    if "S" in cig or "H" in cig:
        continue
    tlen = int(f[8])
    pos = int(f[3]) - 1
    if tlen == 0:
        continue
    if tlen > 0:
        left, right = pos, pos + tlen
    else:
        right = pos + ref_len(cig)
        left = right + tlen
    rev = bool(flag & 0x10)
    if mode == "R":
        strand = "+" if rev else "-"
    elif mode == "F":
        strand = "-" if rev else "+"
    else:
        strand = "."
    sys.stdout.write(f"{left}\t{right}\t{strand}\n")
