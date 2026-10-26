#!/usr/bin/env python3
"""nuc_cm_speed.py - speed and agreement of the rRNA end placement settings (mito_nuc.cm_ends).

The ends to place are the complete 18S / 5.8S / 28S ends of a mitom nuc run made without the model step
(--no_cm, i.e. the reference projection). Each setting places them; the reference is the covariance model on
+-600 nt windows (the first implementation). Reported per gene: ends, placed, equal to the reference, time.

Usage:
  python3 benchmark/nuc_cm_speed.py --seqs records.gb [more.gb ...] --pred nuc_features.tsv [more.tsv ...] -o dir
       [--threads 4]
"""
import argparse
import csv
import re
import sys
import time
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_nuc as nuc  # noqa: E402

SETTINGS = [("CM 600", dict(win=600)), ("CM 150", dict(win=150)),
            ("HMM 150 (<= 20 col)", dict(win=150, hmm_only=True, max_extend=20)),
            ("HMM 150 (<= 50 col)", dict(win=150, hmm_only=True, max_extend=50))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seqs", nargs="+", type=Path, required=True)
    ap.add_argument("--pred", nargs="+", type=Path, required=True)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    seqs = {}
    for f in a.seqs:
        fmt = "genbank" if f.suffix.lower() in (".gb", ".gbk") else "fasta"
        seqs.update({re.sub(r"[^A-Za-z0-9_.-]", "_", r.id): r.seq for r in SeqIO.parse(str(f), fmt)})
    jobs = {g: [] for g in nuc.RFAM}
    for f in a.pred:
        for r in csv.DictReader(open(f), delimiter="\t"):
            g = r["gene"]
            if g not in nuc.RFAM:
                continue
            st = 1 if r["strand"] == "+" else -1
            for side, p, part in ((5, r["start"] if st == 1 else r["end"], r["partial_5"]),
                                  (3, r["end"] if st == 1 else r["start"], r["partial_3"])):
                if part == "0":
                    jobs[g].append((r["contig"], st, side, int(p)))
    rows = []
    for g, (acc, _) in nuc.RFAM.items():
        cm = nuc.LIB / "rfam" / f"{acc}.cm"
        res = {}
        for name, kw in SETTINGS:
            t0 = time.time()
            res[name] = nuc.cm_ends(jobs[g], seqs, g, cm, a.out, a.threads, **kw)
            dt = time.time() - t0
            ref = res[SETTINGS[0][0]]
            same = sum(1 for i, p in res[name].items() if ref.get(i) == p)
            rows.append([g, name, len(jobs[g]), len(res[name]), same, len(ref), round(dt, 1)])
            print("\t".join(map(str, rows[-1])), flush=True)
    with open(a.out / "report.tsv", "w") as o:
        o.write("gene\tsetting\tends\tplaced\tequal_to_CM600\tCM600_placed\tseconds\n")
        o.writelines("\t".join(map(str, r)) + "\n" for r in rows)


if __name__ == "__main__":
    main()
