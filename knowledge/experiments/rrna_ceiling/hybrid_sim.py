"""Offline simulation: MitoM rRNA ends, replaced per end by the relatives' consensus projection when they agree.

Rule(min_votes, max_shift, min_ident): among the K projected positions (other species), the most common position
p (by count, ties by summed bitscore) is adopted for an end when it has >= min_votes votes, the voters' rRNA
identity is >= min_ident, and |p - MitoM end| <= max_shift; otherwise MitoM's end is kept.  Scored against ALL
NCBI ends (the truth is never used to choose).  usage: hybrid_sim.py proj.tsv features.tsv
"""
import csv, re, sys
from collections import defaultdict, Counter

proj, feats = sys.argv[1:3]


def loc(t):
    m = re.match(r"(\d+)\.\.(\d+)([+-])", t or "")
    return (int(m.group(1)) - 1, int(m.group(2)), 1 if m.group(3) == "+" else -1) if m else None


mit = {}
for r in csv.DictReader(open(feats), delimiter="\t"):
    if r["tool"] != "mitom_native" or r["kind"] != "rRNA" or not r["truth"]:
        continue
    t = loc(r["truth"])
    p = loc(r["pred"]) if r["class"] in ("exact", "le3", "gt3") else None
    key = f"{r['accession']}|{r['name']}|{t[0]}"
    if p is None:
        mit[key] = None
    elif t[2] == 1:
        mit[key] = (p[0] - t[0], p[1] - t[1])
    else:
        mit[key] = (t[1] - p[1], t[0] - p[0])
refs = defaultdict(list)
for r in csv.DictReader(open(proj), delimiter="\t"):
    if r["rank"] == "0" or r["d5"] == "" or r["d3"] == "":
        continue
    refs[r["key"]].append((int(r["d5"]), int(r["d3"]), float(r["bits"]), float(r["ident"])))


def run(min_votes, max_shift, min_ident):
    n = both = le3 = changed = 0
    per = Counter()
    for key, m in mit.items():
        n += 1
        if m is None:
            continue
        d = list(m)
        for i in (0, 1):
            votes = Counter()
            wts = Counter()
            for r in refs.get(key, []):
                if r[3] >= min_ident:
                    votes[r[i]] += 1
                    wts[r[i]] += r[2]
            if votes:
                pos = max(votes, key=lambda x: (votes[x], wts[x]))
                if votes[pos] >= min_votes and abs(pos - m[i]) <= max_shift and pos != m[i]:
                    d[i] = pos
                    changed += 1
            per[i] += d[i] == 0
        both += d[0] == 0 and d[1] == 0
        le3 += abs(d[0]) <= 3 and abs(d[1]) <= 3
    return n, both, le3, changed, per


n, b0, l0, _, p0 = run(99, 0, 1.1)
print(f"baseline MitoM: n={n} both exact {b0} ({100*b0/n:.1f}%), within 3 {l0} ({100*l0/n:.1f}%), 5' {p0[0]} 3' {p0[1]}")
for mv in (2, 3):
    for ms in (2, 3, 5, 10):
        for mi in (0.80, 0.85, 0.90):
            n, b, l, c, p = run(mv, ms, mi)
            print(f"votes>={mv} shift<={ms:2d} ident>={mi:.2f}: both exact {b} ({100*b/n:.1f}%, {b-b0:+d}), "
                  f"within 3 {100*l/n:.1f}%, ends changed {c}, 5' {p[0]} 3' {p[1]}")
