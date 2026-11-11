#!/usr/bin/env python3
"""Travisia forbesii (OZ479805, DToL wuTraForb1): three RNA-seq observations checked in DNA reads of the same specimen.

1. COX3 T-run (genome 6 T after position 4227): 32 % of RNA reads have 7 T. DNA with 7 T in a similar share =
   heteroplasmy (two mtDNA types); DNA all 6 T = an RNA-level change (transcriptional slippage / editing).
2. ATP8 position 3706 (MitoM skipped it in a frameshift join): is the base in the DNA?
3. COX1 group II intron (274..1525): DNA should hold exon|intron and intron|exon, not exon|exon.
The tests match flanking sequences inside each read (both strands), independent of how the aligner placed indels.

Usage: travisia_dna_check.py -g OZ479805.gb --reads mito_reads.fa [--label HiFi]
"""
import argparse
import re
from collections import Counter

from Bio import SeqIO

RC = str.maketrans("ACGT", "TGCA")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-g", "--genbank", required=True)
    p.add_argument("--reads", required=True, help="FASTA / FASTQ of reads mapped to the mitogenome")
    p.add_argument("--label", default="")
    a = p.parse_args()
    g = str(SeqIO.read(a.genbank, "genbank").seq).upper()
    fmt = "fastq" if a.reads.endswith(("q", "q.gz")) else "fasta"
    seqs = []
    for r in SeqIO.parse(a.reads, fmt):
        s = str(r.seq).upper()
        seqs += [s, s[::-1].translate(RC)]
    print(f"[{a.label}] {len(seqs) // 2} mitochondrial reads")
    # 1. COX3 T-run: 12-nt flanks around the run that starts after 1-based position 4227
    i = 4227
    while g[i] == "T":
        i += 1
    left, right = g[4227 - 12:4227], g[i:i + 12]
    runs = Counter(len(m.group(1)) for s in seqs for m in re.finditer(left + "(T{1,15})" + right, s))
    print(f"[{a.label}] COX3 T-run (genome {i - 4227} T): reads by length {dict(sorted(runs.items()))}")
    # 2. ATP8 base 3706 (1-based): 15-nt flanks with / without it
    w = g[3706 - 16:3706 + 15]
    print(f"[{a.label}] ATP8 3706: with the base {sum(w in s for s in seqs)}, "
          f"without {sum((w[:15] + w[16:]) in s for s in seqs)}")
    # 3. COX1 intron 274..1525: 15 + 15 nt junctions
    ee, ei, ie = g[258:273] + g[1525:1540], g[258:288], g[1510:1540]
    print(f"[{a.label}] COX1 intron: exon|intron {sum(ei in s for s in seqs)}, intron|exon "
          f"{sum(ie in s for s in seqs)}, exon|exon (spliced) {sum(ee in s for s in seqs)}")


if __name__ == "__main__":
    main()
