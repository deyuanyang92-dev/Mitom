#!/usr/bin/env python3
"""nuc_cm_gold.py - independent gold standard for rDNA gene boundaries (18S, 5.8S, 28S) from Rfam covariance models.

Why: NCBI submitters' boundaries disagree between studies (5.8S / 28S conventions, 28S run into the spacer), so
they cannot be the truth. The mature ends of the curated RefSeq rRNAs (human, rat, macaque, rabbit, cow, hamster,
Drosophila, C. elegans) sit at fixed places of the Rfam models: aligned with cmalign, 18S starts at column 1 and
ends at the last column (0 / 0 nt outside), 5.8S has 1 nt before column 1 and 3 after the last column, 28S 2 and 1
(OFFSET, measured 2026-10-25 on the 28 curated records, Rfam RF01960 / RF00002 / RF02543). The gold boundary of a
gene is the residue aligned to the first / last consensus column, shifted by OFFSET, found with cmsearch on a
window around each end (FLANK nt outside, INNER inside; a hit cut by the window edge is allowed, so the 5' window
hit must start at model column 1 and the 3' window hit end at the last column). It is taken only when the
alignment is confident there: the end column and the next 4 consensus columns inward hold residues of posterior
probability >= 0.95 ('*'), and the window has sequence beyond the end (else the gene is cut by the sequence end
and that boundary is not gold). Windows are placed on mitom's genes and on the NCBI rRNA features (--ncbi).

cmsearch (structure-aware alignment to a covariance model) is independent of mitom nuc (BLAST chains of reference
units, boundary votes), so the gold does not inherit mitom's method.

Usage:
  python3 benchmark/nuc_cm_gold.py --seqs records.gb --pred mitom/nuc_features.tsv --rfam DIR -o gold.tsv
       [--ncbi records.gb] [--cpu 16]
  --rfam: directory with RF01960.cm, RF00002.cm, RF02543.cm (https://rfam.org/family/<acc>/cm)
Writes gold.tsv (record, gene, strand, start, end, start_ok, end_ok) and prints mitom (and NCBI) agreement.
"""
import argparse
import csv
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_env  # noqa: E402
import mito_nuc as nuc  # noqa: E402

MODEL = {"18S": "RF01960", "5.8S": "RF00002", "28S": "RF02543"}
CLEN = {"18S": 1831, "5.8S": 154, "28S": 3401}
OFFSET = {"18S": (0, 0), "5.8S": (1, 3), "28S": (2, 1)}
FLANK, INNER = 300, 400   # end windows: FLANK nt outside the gene, INNER nt inside (5.8S: the whole gene)
CONF = 5                  # consensus columns at the end that must hold residues of posterior probability >= 0.95


def read_sto(path):
    """Infernal -A Stockholm: name -> (aligned seq, PP); #=GC RF"""
    seq, pp, rf = defaultdict(str), defaultdict(str), ""
    for line in open(path):
        if line.startswith("#=GC RF"):
            rf += line.split()[-1]
        elif line.startswith("#=GR") and line.split()[2] == "PP":
            _, k, _, v = line.split()
            pp[k] += v
        elif line.strip() and not line.startswith("#") and not line.startswith("//"):
            k, v = line.split()
            seq[k] += v
    return rf, seq, pp


def confident(rf, s, pp, side):
    cons = [i for i, c in enumerate(rf) if c not in ".-~"]
    cols = cons[:CONF] if side == 5 else cons[::-1][:CONF]
    return all(s[i] not in ".-~" and pp[i] == "*" for i in cols)


def windows(seqs, pred, ncbi):
    """(record, gene, strand, gene 5' end, gene 3' end) of every predicted gene and NCBI rRNA feature, one per
    place (the first of overlapping ones on a strand)"""
    sp = []
    for r in pred:
        if r["gene"] in MODEL:
            st = 1 if r["strand"] == "+" else -1
            a, b = int(r["start"]), int(r["end"])
            sp.append((r["contig"], r["gene"], st, a if st == 1 else b, b if st == 1 else a))
    for rec in ncbi or []:
        for f in rec.features:
            g = nuc.rrna_name(f)
            if g in MODEL:
                st = f.location.strand or 1
                a, b = int(f.location.start) + 1, int(f.location.end)
                sp.append((rec.id, g, st, a if st == 1 else b, b if st == 1 else a))
    out = []
    for w in sp:
        if not any(w[:3] == o[:3] and min(w[3], w[4]) <= max(o[3], o[4]) and min(o[3], o[4]) <= max(w[3], w[4])
                   for o in out):
            out.append(w)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seqs", type=Path, required=True)
    ap.add_argument("--pred", type=Path, required=True)
    ap.add_argument("--ncbi", type=Path, help="GenBank records whose own rRNA features are also scored")
    ap.add_argument("--rfam", type=Path, required=True)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--cpu", type=int, default=16)
    a = ap.parse_args()
    fmt = "genbank" if a.seqs.suffix.lower() in (".gb", ".gbk", ".genbank") else "fasta"
    seqs = {re.sub(r"[^A-Za-z0-9_.-]", "_", r.id): r.seq for r in SeqIO.parse(str(a.seqs), fmt)}
    pred = list(csv.DictReader(open(a.pred), delimiter="\t"))
    ncbi = list(SeqIO.parse(str(a.ncbi), "genbank")) if a.ncbi else None
    work = a.out.with_suffix(".work")
    work.mkdir(parents=True, exist_ok=True)
    cmsearch = mito_env.find("cmsearch")
    gold = []
    for g, acc in MODEL.items():
        ws = [w for w in windows(seqs, pred, ncbi) if w[1] == g]
        if not ws:
            continue
        fa, meta = work / f"{g}.fa", {}
        with open(fa, "w") as o:
            for i, (c, _, st, p5, p3) in enumerate(ws):
                n = len(seqs[c])
                inner = INNER if g != "5.8S" else abs(p3 - p5) + 1
                for side, p in ((5, p5), (3, p3)):
                    # window in transcript orientation: 5' window [p5 - FLANK, p5 + inner], 3' [p3 - inner, p3 + FLANK]
                    lo_t, hi_t = (p - FLANK, p + inner) if side == 5 else (p - inner, p + FLANK)
                    if st == 1:
                        lo, hi = max(1, lo_t), min(n, hi_t)
                        s = seqs[c][lo - 1:hi]
                    else:
                        lo, hi = max(1, 2 * p - hi_t), min(n, 2 * p - lo_t)
                        s = seqs[c][lo - 1:hi].reverse_complement()
                    k = f"w{i}_{side}"
                    meta[k] = (i, side, st, lo, hi)
                    o.write(f">{k}\n{s}\n")
        tbl, sto = work / f"{g}.tbl", work / f"{g}.sto"
        if not tbl.exists():
            subprocess.run([cmsearch, "--cpu", str(a.cpu), "--notextw", "-E", "1e-5", "--tblout", str(tbl),
                            "-A", str(sto), str(a.rfam / f"{acc}.cm"), str(fa)], check=True, capture_output=True)
        rf, al, pp = read_sto(sto) if sto.exists() else ("", {}, {})
        best = {}
        for line in open(tbl):
            if line.startswith("#"):
                continue
            x = line.split()
            k, mf, mt, sf, st_, strand, sc = x[0], int(x[5]), int(x[6]), int(x[7]), int(x[8]), x[9], float(x[14])
            if strand == "+" and (k not in best or sc > best[k][-1]):
                best[k] = (mf, mt, sf, st_, sc)
        res = defaultdict(dict)
        o5, o3 = OFFSET[g]
        for k, (mf, mt, sf, st_, sc) in best.items():
            i, side, st, lo, hi = meta[k]
            wl = hi - lo + 1
            name = f"{k}/{sf}-{st_}"
            if side == 5 and mf == 1 and sf - o5 > 1 and name in al and confident(rf, al[name], pp[name], 5):
                t = sf - o5                                  # window position (1-based, transcript orientation)
            elif side == 3 and mt == CLEN[g] and st_ + o3 < wl and name in al and confident(rf, al[name], pp[name], 3):
                t = st_ + o3
            else:
                continue
            res[i][side] = lo + t - 1 if st == 1 else hi - t + 1
        for i, (c, _, st, p5, p3) in enumerate(ws):
            r5, r3 = res[i].get(5), res[i].get(3)
            lo_, hi_ = (r5, r3) if st == 1 else (r3, r5)
            gold.append({"record": c, "gene": g, "strand": "+" if st == 1 else "-", "start": lo_ or "",
                         "end": hi_ or "", "start_ok": int(lo_ is not None), "end_ok": int(hi_ is not None)})
    with open(a.out, "w", newline="") as o:
        w = csv.DictWriter(o, list(gold[0]), delimiter="\t")
        w.writeheader()
        w.writerows(gold)
    print(f"{len(gold)} genes, gold boundaries: {sum(r['start_ok'] + r['end_ok'] for r in gold)}")
    score("mitom", gold, [(r["contig"], r["gene"], r["strand"], int(r["start"]), int(r["end"])) for r in pred])
    if ncbi:
        score("NCBI submitters", gold, [(rec.id, nuc.rrna_name(f), "+" if (f.location.strand or 1) == 1 else "-",
                                         int(f.location.start) + 1, int(f.location.end))
                                        for rec in ncbi for f in rec.features if nuc.rrna_name(f) in MODEL])


def score(name, gold, calls):
    """per gene and side: gold boundaries, exact, <= 1 nt, called at all (a call of that gene overlapping the gold
    gene on its strand)"""
    by = defaultdict(list)
    for c, g, st, s, e in calls:
        by[(c, g, st)].append((s, e))
    tab = defaultdict(lambda: [0, 0, 0, 0])
    for r in gold:
        for side, ok in (("start", r["start_ok"]), ("end", r["end_ok"])):
            if not ok:
                continue
            t = tab[(r["gene"], side)]
            t[0] += 1
            gpos = int(r[side])
            cand = by.get((r["record"], r["gene"], r["strand"]), [])
            if not cand:
                continue
            p = min((x[0] if side == "start" else x[1] for x in cand), key=lambda x: abs(x - gpos))
            t[1] += p == gpos
            t[2] += abs(p - gpos) <= 1
            t[3] += 1
    print(f"## {name} vs CM gold   gene side: gold n, exact %, <=1 nt %, called %")
    tot = [0, 0, 0, 0]
    for (g, side), t in sorted(tab.items()):
        print(f"{g:5s} {side:5s} {t[0]:5d}  {100 * t[1] / t[0]:5.1f}  {100 * t[2] / t[0]:5.1f}  {100 * t[3] / t[0]:5.1f}")
        tot = [x + y for x, y in zip(tot, t)]
    if tot[0]:
        print(f"all         {tot[0]:5d}  {100 * tot[1] / tot[0]:5.1f}  {100 * tot[2] / tot[0]:5.1f}  {100 * tot[3] / tot[0]:5.1f}")


if __name__ == "__main__":
    main()
