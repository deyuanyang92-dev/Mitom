"""ORF finding by length only (as NCBI ORFfinder / Geneious 'Find ORFs': every ORF from a start codon of at least L
codons, here in the sequence no exon covers): recall on make_truth.py and ORFs it would write in the 362 intron-free
RefSeq records (../intron_algo_2026-11-09/negatives_set.tsv).
Usage: python3 length_baseline.py <work> <intron-free gb folder>
"""
import csv
import glob
import sys
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import mito_orf as m  # noqa: E402

W, NEG = Path(sys.argv[1]), sys.argv[2]
truth = [t for t in csv.DictReader(open(W / "recall" / "truth.tsv"), delimiter="\t")
         if t["record"] != "NC_137528.1" and "atpase subunit" not in t["product"].lower()]
cands = {}
for f in glob.glob(str(W / "recall" / "in" / "*.gb")):
    r = SeqIO.read(f, "genbank")
    cands[r.id] = m.candidates(r, m.code_of(r, 5), 25)
neg = []
for f in glob.glob(f"{NEG}/*.gb"):
    r = SeqIO.read(f, "genbank")
    neg.append(m.candidates(r, m.code_of(r, 5), 25))
print("min_codons\trecall\tof\tintron_free_orfs\tin_records")
for L in (25, 50, 75, 100, 150):
    hit = sum(any(o.strand == (1 if t["strand"] == "1" else -1) and o.codons >= L and
                  max(0, min(o.end, int(t["end"])) - max(o.start, int(t["start"]) - 1)) >=
                  0.5 * (int(t["end"]) - int(t["start"]) + 1) for o in cands[t["record"]]) for t in truth)
    print(f"{L}\t{hit}\t{len(truth)}\t{sum(sum(o.codons >= L for o in c) for c in neg)}\t"
          f"{sum(any(o.codons >= L for o in c) for c in neg)}")
