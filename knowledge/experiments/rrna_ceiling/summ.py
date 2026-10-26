"""Summaries of proj_<mode>.tsv: top-1 / vote projections vs MitoM's current ends, per end and both ends."""
import csv, sys
from collections import defaultdict, Counter

path = sys.argv[1]
rows = list(csv.DictReader(open(path), delimiter="\t"))
by = defaultdict(list)
for r in rows:
    by[r["key"]].append(r)


def num(x):
    return None if x in ("", None) else int(x)


def ex(d5, d3, tol=0):
    return d5 is not None and d3 is not None and abs(d5) <= tol and abs(d3) <= tol


tot = len(by)
mitom_exact = sum(1 for k, v in by.items() if v[0]["mitom_class"] == "exact")
mitom_le3 = sum(1 for k, v in by.items() if v[0]["mitom_class"] in ("exact", "le3"))
top1 = Counter(); vote = Counter(); oracle = Counter(); end5 = Counter(); end3 = Counter()
rel_tab = defaultdict(Counter)
id_tab = defaultdict(Counter)
for k, v in by.items():
    c = [r for r in v if r["rank"] != "0" and r["d5"] != "" and r["d3"] != ""]
    if not c:
        top1["no_ref"] += 1
        continue
    r0 = c[0]
    d5, d3 = num(r0["d5"]), num(r0["d3"])
    for tol in (0, 1, 3, 10):
        if ex(d5, d3, tol):
            top1[tol] += 1
    end5[abs(d5) == 0] += 1; end3[abs(d3) == 0] += 1
    rel_tab[r0["relation"]]["n"] += 1
    rel_tab[r0["relation"]]["exact"] += ex(d5, d3)
    idb = float(r0["ident"])
    b = ">=0.95" if idb >= .95 else ">=0.90" if idb >= .9 else ">=0.85" if idb >= .85 else ">=0.80" if idb >= .8 else "<0.80"
    id_tab[b]["n"] += 1; id_tab[b]["exact"] += ex(d5, d3); id_tab[b]["mitom"] += v[0]["mitom_class"] == "exact"
    # bitscore-weighted vote per end (separately)
    w5, w3 = Counter(), Counter()
    for r in c:
        w5[num(r["d5"])] += float(r["bits"]); w3[num(r["d3"])] += float(r["bits"])
    v5, v3 = w5.most_common(1)[0][0], w3.most_common(1)[0][0]
    for tol in (0, 1, 3, 10):
        if ex(v5, v3, tol):
            vote[tol] += 1
    # oracle: any of the K references exact at both ends
    oracle["any"] += any(ex(num(r["d5"]), num(r["d3"])) for r in c)
print(f"truth rRNAs: {tot}")
print(f"MitoM now: exact {mitom_exact} ({100*mitom_exact/tot:.1f}%), within 3 {mitom_le3} ({100*mitom_le3/tot:.1f}%)")
print("top-1 projection: " + ", ".join(f"<= {t}: {top1[t]} ({100*top1[t]/tot:.1f}%)" for t in (0, 1, 3, 10)) + f"; no ref {top1['no_ref']}")
print(f"  per end exact: 5' {end5[True]}  3' {end3[True]}")
print("weighted vote:   " + ", ".join(f"<= {t}: {vote[t]} ({100*vote[t]/tot:.1f}%)" for t in (0, 1, 3, 10)))
print(f"oracle (any of K refs exact): {oracle['any']} ({100*oracle['any']/tot:.1f}%)")
print("top-1 by relation:", {k: (v['n'], f"{100*v['exact']/v['n']:.0f}%") for k, v in sorted(rel_tab.items())})
print("top-1 by identity:", {k: (v['n'], f"proj {100*v['exact']/v['n']:.0f}%", f"mitom {100*v['mitom']/v['n']:.0f}%") for k, v in sorted(id_tab.items())})
