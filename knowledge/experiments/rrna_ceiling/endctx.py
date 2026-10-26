"""Per rRNA end: truth neighbour + gap, MitoM neighbour + gap, and the offset; tabulates the causes of end errors.

usage: endctx.py SET.tsv ROOT [TOOL_DIR=mitom_native] > ends.tsv  (summary on stderr)
"""
import csv, os, re, sys
from collections import Counter
sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
from mito_bench_eval import norm_features, match   # noqa: E402
from Bio import SeqIO                               # noqa: E402

set_tsv, root = sys.argv[1], sys.argv[2]
tool = sys.argv[3] if len(sys.argv) > 3 else "mitom_native"
meta = {r["accession"]: r for r in csv.DictReader(open(set_tsv), delimiter="\t")}


def neighbours(feats, f, n):
    """(name5, gap5, bnd5, name3, gap3, bnd3) in transcript direction; gap 0 = contiguous, < 0 overlap.
    bnd = the neighbour's boundary coordinate facing the rRNA (0-based genome position of its last/first base)."""
    best_l = best_r = None
    for g in feats:
        if g is f:
            continue
        for k in (-n, 0, n):
            a, b = g.start + k, g.end + k
            if b <= f.start + 60 and a < f.start:            # ends before (or slightly inside) our start
                d = f.start - b
                if best_l is None or abs(d) < abs(best_l[0]):
                    best_l = (d, g, b - 1)
            if a >= f.end - 60 and b > f.end:
                d = a - f.end
                if best_r is None or abs(d) < abs(best_r[0]):
                    best_r = (d, g, a)
    def lab(x):
        return ("none", None, None) if x is None else (f"{x[1].kind}:{x[1].name}", x[0], x[2])
    L, R = lab(best_l), lab(best_r)
    return (L, R) if f.strand > 0 else (R, L)


out = csv.writer(sys.stdout, delimiter="\t")
out.writerow(["accession", "phylum", "gene", "end", "truth_nb", "truth_gap", "tool_nb", "tool_gap", "delta",
              "nb_bnd_diff"])
summary = Counter()
for acc in sorted(meta):
    gb = os.path.join(root, "gb", acc + ".gb")
    pr = os.path.join(root, "runs", tool, acc, "final.gbf")
    if not (os.path.exists(gb) and os.path.exists(pr)):
        continue
    try:
        trec = SeqIO.read(gb, "genbank")
        precs = list(SeqIO.parse(pr, "genbank"))
    except Exception:
        continue
    if len(precs) != 1:
        continue
    prec = precs[0]
    code = int(meta[acc]["genetic_code"] or 5)
    T = norm_features(trec, code)
    P = norm_features(prec, code)
    n = len(trec)
    if len(prec) != n:
        continue
    for t in [x for x in T if x.kind == "rRNA"]:
        cands = [p for p in P if p.kind == "rRNA" and p.name == t.name and p.strand == t.strand and
                 min(p.end, t.end) - max(p.start, t.start) > 0.5 * (t.end - t.start)]
        if not cands:
            continue
        p = cands[0]
        tn = neighbours(T, t, n)
        pn = neighbours(P, p, n)
        if t.strand > 0:
            d5, d3 = p.start - t.start, p.end - t.end
        else:
            d5, d3 = t.end - p.end, t.start - p.start
        for end, d, ti, pi in (("5'", d5, tn[0], pn[0]), ("3'", d3, tn[1], pn[1])):
            nbd = "" if ti[2] is None or pi[2] is None or ti[0] != pi[0] else (pi[2] - ti[2]) * (1 if t.strand > 0 else -1)
            out.writerow([acc, meta[acc]["phylum"], t.name, end, ti[0], ti[1], pi[0], pi[1], d, nbd])
            tg, pg = ti[1], pi[1]
            cat = ("exact" if d == 0 else
                   "same_nb_both_contig_nbdiff" if tg == 0 and pg == 0 and ti[0] == pi[0] else
                   "truth_contig_tool_not" if tg == 0 else
                   "tool_contig_truth_not" if pg == 0 else "neither_contig")
            summary[(end, cat, "le3" if abs(d) <= 3 else "gt3" if d else "0")] += 1
for k in sorted(summary):
    print(k, summary[k], file=sys.stderr)
