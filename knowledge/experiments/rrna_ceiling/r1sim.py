"""Replay MitoM's rRNA end rule (R1) with variants on finished runs and score the rRNAs against NCBI.

Pre-R1 rRNA ends are rebuilt from reports/curation_report.tsv (first R1 'old' -> last R1 'new', as 5' / 3'
shifts in transcript direction, which survive the re-orientation lift); every other feature is MitoM's final one.
Policies:  current  = trim overlaps, fill gaps <= 300 nt to the nearest feature (any type, any strand);
           context  = as current, but fill only where the RefSeq contiguity rate of (phylum, gene, end, neighbour
                      type) is >= THETA (conv_ends.tsv, benchmark species excluded); otherwise keep the pre-R1 end.
usage: r1sim.py SET.tsv ROOT conv_ends.tsv [split]
"""
import csv, os, re, sys
from collections import Counter, defaultdict
sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
from mito_bench_eval import norm_features, match, Feat   # noqa: E402
from Bio import SeqIO                                      # noqa: E402

MAX_FILL = 300


def load_rates(path):
    c = defaultdict(Counter)
    for r in csv.DictReader(open(path), delimiter="\t"):
        k = (r["phylum"], r["gene"], r["end"], r["nb_type"])
        c[k]["n"] += 1
        c[k]["c0"] += r["gap"] == "0"
    return {k: (v["c0"] / v["n"], v["n"]) for k, v in c.items()}


def pre_r1(report, n_rrna_names):
    """name -> list of (d5, d3) shifts applied by R1 (transcript direction), in report order."""
    out = defaultdict(list)
    if not os.path.exists(report):
        return out
    rows = [r for r in csv.DictReader(open(report), delimiter="\t") if r["feature"] == "rRNA" and r["rule"] == "R1"]
    first, last = {}, {}
    for r in rows:
        k = (r["gene"], r["strand"])
        first.setdefault(k, r["old"])
        last[k] = r["new"]
    for (g, st), old in first.items():
        a0, b0 = map(int, old.split(".."))
        a1, b1 = map(int, last[(g, st)].split(".."))
        if st == "+":
            out[g].append((a1 - a0, b1 - b0))
        else:
            out[g].append((b0 - b1, a0 - a1))
    return out


def nb_type(kind):
    return {"PCG": "CDS", "tRNA": "tRNA", "rRNA": "rRNA"}.get(kind, "none")


def apply_r1(rr, others, n, policy, phylum, rates, theta):
    """rr = (start0, end_excl, strand, name) pre-R1; others = Feats. Returns new (start0, end_excl)."""
    s, e = rr[0] + 1, rr[1]                     # 1-based inclusive, e may exceed n
    neigh = []
    for f in others:
        if f.kind not in ("PCG", "tRNA", "rRNA"):
            continue
        for k in (-n, 0, n):
            neigh.append((f.start + 1 + k, f.end + k, f.kind))
    changed = True
    while changed:
        changed = False
        for ns, ne, _ in sorted(neigh):
            if ne < s or ns > e:
                continue
            if ns <= s and ne >= e:
                continue
            if ns <= s:
                s = ne + 1
            elif ne >= e:
                e = ns - 1
            else:
                keep_left = ns - s >= e - ne
                lost = e - ns + 1 if keep_left else ne - s + 1
                if lost * 10 > e - s + 1 and lost > 250:
                    continue
                if keep_left:
                    e = ns - 1
                else:
                    s = ne + 1
            changed = True
    left = max([(ne, k) for ns, ne, k in neigh if ne < s] or [(0, "none")])
    right = min([(ns, k) for ns, ne, k in neigh if ns > e] or [(10 ** 9, "none")])
    fill_l = fill_r = True
    if policy == "context":
        strand = rr[2]
        # left side is the 5' end on +, the 3' end on -
        for side, nb in (("L", left), ("R", right)):
            end = ("5'" if side == "L" else "3'") if strand > 0 else ("3'" if side == "L" else "5'")
            rate, cnt = rates.get((phylum, rr[3], end, nb_type(nb[1])), (None, 0))
            ok = rate is None or cnt < 10 or rate >= theta
            if side == "L":
                fill_l = ok
            else:
                fill_r = ok
    if fill_l and 0 < s - left[0] - 1 <= MAX_FILL:
        s = left[0] + 1
    if fill_r and 0 < right[0] - e - 1 <= MAX_FILL:
        e = right[0] - 1
    return s - 1, e


def main():
    set_tsv, root, conv = sys.argv[1:4]
    split = sys.argv[4] if len(sys.argv) > 4 else "all"
    meta = [r for r in csv.DictReader(open(set_tsv), delimiter="\t") if split == "all" or r["split"] == split]
    rates = load_rates(conv)
    variants = [("final", None, None), ("current", "current", 0)] + \
               [(f"context>={t}", "context", t) for t in (0.2, 0.3, 0.4, 0.5, 0.6)]
    score = {v[0]: Counter() for v in variants}
    sanity = Counter()
    for m in meta:
        acc = m["accession"]
        gb = os.path.join(root, "gb", acc + ".gb")
        fin = os.path.join(root, "runs", "mitom_native", acc, "final.gbf")
        rep = os.path.join(root, "runs", "mitom_native", acc, "raw", "reports", "curation_report.tsv")
        if not os.path.exists(fin):
            continue
        try:
            trec = SeqIO.read(gb, "genbank")
            precs = list(SeqIO.parse(fin, "genbank"))
        except Exception:
            continue
        if len(precs) != 1 or len(precs[0]) != len(trec):
            continue
        code = int(m["genetic_code"] or 5)
        T = [f for f in norm_features(trec, code) if f.kind == "rRNA"]
        P = norm_features(precs[0], code)
        n = len(trec)
        shifts = pre_r1(rep, None)
        used = Counter()
        rr_pre = []
        for f in P:
            if f.kind != "rRNA":
                continue
            lst = shifts.get(f.name, [])
            d5, d3 = lst[used[f.name]] if used[f.name] < len(lst) else (0, 0)
            used[f.name] += 1
            if f.strand > 0:
                rr_pre.append((f.start - d5, f.end - d3, f.strand, f.name, f))
            else:
                rr_pre.append((f.start + d3, f.end + d5, f.strand, f.name, f))
        for vname, pol, theta in variants:
            if pol is None:
                preds = [f for f in P if f.kind == "rRNA"]
            else:
                preds = []
                for s0, e0, st, name, f in rr_pre:
                    others = [g for g in P if g is not f and g.kind in ("PCG", "tRNA", "rRNA")]
                    ns, ne = apply_r1((s0, e0, st, name), others, n, pol, m["phylum"], rates, theta)
                    preds.append(Feat("rRNA", name, st, ns, ne, "", n))
                    if pol == "current":
                        sanity[(ns, ne) == (f.start, f.end)] += 1
            for row in match(T, preds):
                if row["truth"]:
                    score[vname][row["class"]] += 1
                    score[vname]["n"] += 1
                    if row["class"] in ("exact", "le3", "gt3") and row["delta"] != "" and int(row["delta"]) <= 10:
                        score[vname]["le10"] += 1
    print(f"records {len(meta)} split {split}; replay of the current rule reproduces MitoM's final rRNA: "
          f"{sanity[True]} of {sanity[True] + sanity[False]}")
    for v, c in score.items():
        n = c["n"]
        print(f"{v:14s} n={n} exact {c['exact']} ({100*c['exact']/n:.1f}%)  within3 {100*(c['exact']+c['le3'])/n:.1f}%"
              f"  within10 {100*c['le10']/n:.1f}%  missed {c['missed']}")


if __name__ == "__main__":
    main()
