#!/usr/bin/env python3
"""Simulate paired RNA-seq reads of an annotated mitochondrial genome with known truth (mapper benchmark).

What mitochondrial RNA-seq looks like (Ojala et al. 1981; Selifanova et al. 2023):
  - mature transcripts of each protein gene and rRNA, with a poly(A) tail (20-60 A) right after the annotated 3' end
    (it completes incomplete stop codons T / TA); expression differs by gene (log-normal);
  - precursor fragments (polycistronic, unprocessed: 1-4 kb windows anywhere, also across the origin of a circular
    genome) - PRECURSOR share of the fragments, no tail.
The sequenced individual differs from the reference: SNPs at rate --div, 1-3 nt indels at --div / 10; sequencing
errors at --err. Fragments 220 +- 40 nt (whole transcript when shorter); read 1 antisense (dUTP library).

Writes <out>_1.fq.gz, <out>_2.fq.gz and <out>.truth.tsv.gz (read, mate, ref_left (0-based), ref_right (exclusive,
may exceed the genome length when a read crosses the origin), tail_nt in the read, kind), <out>.genes.tsv (gene,
strand, start, end, weight: the true 3' ends that carry a tail).
"""
from __future__ import annotations

import argparse
import gzip
import random

import numpy as np
from Bio import SeqIO

COMP = str.maketrans("ACGTN", "TGCAN")
PRECURSOR = 0.15


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def mutate(ref: str, div: float, rng: random.Random):
    """Individual genome and, for each of its positions, the reference position (-1 = inserted base)."""
    out, pos = [], []
    i, n = 0, len(ref)
    while i < n:
        r = rng.random()
        if r < div:                                           # SNP
            out.append(rng.choice([b for b in "ACGT" if b != ref[i]])); pos.append(i); i += 1
        elif r < div * 1.05:                                  # insertion (1-3)
            for _ in range(rng.randint(1, 3)):
                out.append(rng.choice("ACGT")); pos.append(-1)
            out.append(ref[i]); pos.append(i); i += 1
        elif r < div * 1.10:                                  # deletion (1-3)
            i += rng.randint(1, 3)
        else:
            out.append(ref[i]); pos.append(i); i += 1
    return "".join(out), np.array(pos)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-g", "--genbank", required=True)
    p.add_argument("-o", "--out", required=True)
    p.add_argument("-n", "--pairs", type=int, default=200000)
    p.add_argument("-l", "--read_len", type=int, default=100)
    p.add_argument("--div", type=float, default=0.0)
    p.add_argument("--err", type=float, default=0.003)
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args()
    rng = random.Random(a.seed)
    nrng = np.random.default_rng(a.seed)
    rec = SeqIO.read(a.genbank, "genbank")
    ref = str(rec.seq).upper()
    n = len(ref)
    ind, ipos = mutate(ref, a.div, rng)
    m = len(ind)
    # reference position -> individual position (deleted reference bases map to the next kept base)
    r2i = np.full(n + 1, m, dtype=np.int64)
    kept = ipos >= 0
    r2i[ipos[kept]] = np.nonzero(kept)[0]
    for k in range(n - 1, -1, -1):
        if r2i[k] == m and k < n:
            r2i[k] = r2i[k + 1]
    genes = []
    for f in rec.features:
        if f.type not in ("CDS", "rRNA") or len(f.location.parts) > 1:
            continue
        s, e, st = int(f.location.start), int(f.location.end), f.location.strand or 1
        name = f.qualifiers.get("gene", f.qualifiers.get("product", ["?"]))[0]
        genes.append((name, st, s, e, float(nrng.lognormal(0, 1))))
    with open(a.out + ".genes.tsv", "w") as fh:
        fh.write("gene\tstrand\tstart\tend\tweight\n")
        for g in genes:
            fh.write(f"{g[0]}\t{g[1]}\t{g[2]}\t{g[3]}\t{g[4]:.3f}\n")
    w = np.array([g[4] for g in genes])
    w = w / w.sum()
    ind2 = ind + ind                                          # circular
    f1 = gzip.open(a.out + "_1.fq.gz", "wt", compresslevel=1)
    f2 = gzip.open(a.out + "_2.fq.gz", "wt", compresslevel=1)
    tr = gzip.open(a.out + ".truth.tsv.gz", "wt", compresslevel=1)
    tr.write("read\tmate\tref_left\tref_right\ttail_nt\tkind\n")
    L = a.read_len
    q = "I" * L

    def err(s: str) -> str:
        b = list(s)
        for k in range(len(b)):
            if rng.random() < a.err:
                b[k] = rng.choice([x for x in "ACGT" if x != b[k]])
        return "".join(b)

    def span(il: int, ir: int):
        """Reference interval of individual interval [il, ir) (ir may pass m: circular)."""
        idx = [k % m for k in range(il, ir)]
        rp = [int(ipos[k]) for k in idx if ipos[k] >= 0]
        if not rp:
            return -1, -1
        left = rp[0]
        right = left + sum(1 for k in idx if ipos[k] >= 0)    # reference length covered (incl. deletions approx.)
        # deletions inside: use last ref pos (unwrapped)
        last = rp[-1] if rp[-1] >= left else rp[-1] + n
        return left, max(right, last + 1)

    for i in range(a.pairs):
        flen = max(L, min(600, int(rng.gauss(220, 40))))
        if rng.random() < PRECURSOR:
            kind = "precursor"
            st = rng.choice((1, -1))
            tlen = rng.randint(1000, 4000)
            t0 = rng.randrange(m)
            tseq_g = ind2[t0: t0 + tlen] if t0 + tlen <= 2 * m else ind2[t0:]
            tail = 0
            gl = t0                                           # individual coordinate of transcript plus-strand left
            tseq = tseq_g if st == 1 else rc(tseq_g)
        else:
            g = genes[nrng.choice(len(genes), p=w)]
            kind = "mature:" + g[0]
            _, st, s, e, _ = g
            il, ir = int(r2i[s]), int(r2i[e])
            tseq_g = ind[il:ir]
            tail = rng.randint(20, 60)
            gl = il
            tseq = (tseq_g if st == 1 else rc(tseq_g)) + "A" * tail
        tl = len(tseq)
        if tl < L:
            continue
        flen = min(flen, tl)
        fs = rng.randrange(tl - flen + 1)
        frag = tseq[fs: fs + flen]
        glen = tl - tail                                       # genomic part of the transcript
        # transcript coordinate -> individual plus-strand coordinate
        def to_ind(a0: int, a1: int):
            a1 = min(a1, glen)
            if a1 <= a0:
                return None
            return (gl + a0, gl + a1) if st == 1 else (gl + glen - a1, gl + glen - a0)
        r2 = frag[:L]                                          # read 2: sense, fragment start
        r1 = rc(frag[-L:])                                     # read 1: antisense, fragment end
        for mate, seq, (a0, a1) in ((1, r1, (fs + flen - L, fs + flen)), (2, r2, (fs, fs + L))):
            tail_in = max(0, a1 - max(a0, glen))
            iv = to_ind(a0, a1)
            left, right = span(*iv) if iv else (-1, -1)
            tr.write(f"r{i}\t{mate}\t{left}\t{right}\t{tail_in}\t{kind}\n")
            (f1 if mate == 1 else f2).write(f"@r{i}/{mate}\n{err(seq)}\n+\n{q[:len(seq)]}\n")
    for fh in (f1, f2, tr):
        fh.close()


if __name__ == "__main__":
    main()
