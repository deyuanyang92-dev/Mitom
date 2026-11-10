#!/usr/bin/env python3
"""Score each mapper's SAM against the simulation truth (simulate.py).

Per mapper: share of reads aligned; placed correctly (leftmost aligned base within TOL of the truth); aligned with
MAPQ >= 20 (what mito_rna uses) and of those placed wrong; reads with a spliced (N) alignment - false, animal
mitochondrial transcripts are not spliced; reads with a poly(A) tail >= 6 nt whose tail mito_rna sees (clipped A / T
run); poly(A) sites (mito_rna's own rule, >= 3 reads) at the true 3' ends (recall) and called sites at a true end
(precision), both within SITE_TOL; per-gene depth (MAPQ >= 20) against the true depth: median |log2 ratio|; reads
crossing the origin aligned; time (s) and peak memory (MB).

Usage: evaluate.py --mitom <Mitom dir> -g genome.gb --sim <prefix> --maps <dir> [--label x] >> table.tsv
"""
from __future__ import annotations

import argparse
import gzip
import math
import re
import statistics
import sys
from pathlib import Path

import numpy as np

TOL = 5
SITE_TOL = 3
COLS = ["label", "mapper", "reads", "aligned", "correct", "mapq20", "mapq20_wrong", "spliced", "tail_reads",
        "tail_seen", "polya_true_sites", "polya_recall", "polya_called", "polya_precision", "depth_log2_err",
        "origin_reads", "origin_aligned", "seconds", "max_mb"]


def time_of(p: Path):
    if not p.exists():
        return "", ""
    t = p.read_text()
    mw = re.search(r"Elapsed \(wall clock\) time.*: (\S+)", t)
    mm = re.search(r"Maximum resident set size \(kbytes\): (\d+)", t)
    if not (mw and mm):
        return "", ""
    sec = sum(float(x) * 60 ** i for i, x in enumerate(reversed(mw.group(1).split(":"))))
    mb = int(mm.group(1)) / 1024
    return f"{sec:.1f}", f"{mb:.0f}"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mitom", required=True)
    p.add_argument("-g", "--genbank", required=True)
    p.add_argument("--sim", required=True)
    p.add_argument("--maps", required=True)
    p.add_argument("--label", default="")
    p.add_argument("--header", action="store_true")
    a = p.parse_args()
    sys.path.insert(0, a.mitom)
    import mito_rna
    from Bio import SeqIO
    rec = SeqIO.read(a.genbank, "genbank")
    n = len(rec.seq)
    truth = {}
    with gzip.open(a.sim + ".truth.tsv.gz", "rt") as fh:
        next(fh)
        for ln in fh:
            r, m, left, right, tail, _kind = ln.rstrip("\n").split("\t")
            truth[(r, int(m))] = (int(left), int(right), int(tail))
    genes = [ln.rstrip("\n").split("\t") for ln in open(a.sim + ".genes.tsv")][1:]
    true_sites = [(int(e) - 1, 1) if g[1] == "1" else (int(s), -1) for g in genes for _, _, s, e in [g[:4]]]
    # true depth (both mates, genomic part)
    tdep = np.zeros(n, dtype=np.int64)
    for (left, right, _) in truth.values():
        if left < 0:
            continue
        if right <= n:
            tdep[left:right] += 1
        else:
            tdep[left:] += 1
            tdep[: right - n] += 1
    gene_iv = [(int(g[2]), int(g[3])) for g in genes]
    if a.header:
        print("\t".join(COLS))
    for sam in sorted(Path(a.maps).glob("*.sam")):
        name = sam.stem
        seen = set()
        aligned = correct = q20 = q20wrong = spliced = tailr = tails = orig = origa = 0
        lines = []
        for ln in open(sam):
            if ln.startswith("@"):
                lines.append(ln)
                continue
            c = ln.split("\t", 10)
            flag = int(c[1])
            if flag & 0x900:
                continue
            rid = re.sub(r"/[12]$", "", c[0])
            mate = 2 if flag & 0x80 else 1
            key = (rid, mate)
            if key in seen or key not in truth:
                continue
            seen.add(key)
            left, right, tail = truth[key]
            if flag & 0x4:
                continue
            lines.append(ln)
            aligned += 1
            if right > n:
                origa += 1
            pos = int(c[3]) - 1
            ok = left >= 0 and min(abs(pos - left), n - abs(pos - left)) <= TOL
            correct += ok
            if int(c[4]) >= mito_rna.MIN_MAPQ:
                q20 += 1
                q20wrong += (not ok and left >= 0)
            if "N" in c[5]:
                spliced += 1
            if tail >= mito_rna.POLYA_MIN:
                ops = mito_rna.CIGAR.findall(c[5])
                seq = c[9]
                if (ops and ops[-1][1] == "S" and seq[len(seq) - int(ops[-1][0]):].startswith("A" * mito_rna.POLYA_MIN)) \
                        or (ops and ops[0][1] == "S" and seq[: int(ops[0][0])].endswith("T" * mito_rna.POLYA_MIN)):
                    tails += 1
        res = mito_rna.parse_sam(lines, n)
        called = mito_rna.polya_sites(res["polya"])
        hit = lambda p_, s_, pool: any(abs(p_ - q) <= SITE_TOL and s_ == t for q, t in pool)
        rec_ = sum(hit(p_, s_, [(x, y) for x, y, _ in called]) for p_, s_ in true_sites)
        prec = sum(hit(p_, s_, true_sites) for p_, s_, _ in called)
        dep = res["depth"][1] + res["depth"][-1]
        errs = [abs(math.log2((dep[s:e].mean() + 1) / (tdep[s:e].mean() + 1))) for s, e in gene_iv]
        total = len(truth)
        tailr = sum(1 for _, _, t in truth.values() if t >= mito_rna.POLYA_MIN)
        orig = sum(1 for _, r_, _ in truth.values() if r_ > n)
        sec, mb = time_of(sam.with_suffix(".time"))
        pct = lambda x, d: f"{100 * x / d:.2f}" if d else ""
        row = [a.label, name, total, pct(aligned, total), pct(correct, total), pct(q20, total), pct(q20wrong, q20),
               spliced, tailr, pct(tails, tailr), len(true_sites), f"{rec_}/{len(true_sites)}", len(called),
               f"{prec}/{len(called)}", f"{statistics.median(errs):.3f}", orig, pct(origa, orig), sec, mb]
        print("\t".join(map(str, row)), flush=True)


if __name__ == "__main__":
    main()
