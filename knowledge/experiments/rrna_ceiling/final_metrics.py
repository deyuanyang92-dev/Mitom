"""rRNA metrics per tool from mito_bench_eval.py's features.tsv, on the records every listed tool finished.

detected = gene found on the right strand with the right name (>= 50 % overlap); exact / within 3 / within 10 =
both ends; consensus = ends where >= 2 of the 5 closest annotated relatives (other species, project.py) put their
NCBI end at the same homologous position (selection never looks at a tool's output).
usage: final_metrics.py features.tsv proj_species.tsv TOOL [TOOL ...]
"""
import csv
import re
import sys
from collections import Counter, defaultdict

feats, proj, tools = sys.argv[1], sys.argv[2], sys.argv[3:]


def loc(t):
    m = re.match(r"(\d+)\.\.(\d+)([+-])", t or "")
    return (int(m.group(1)) - 1, int(m.group(2)), 1 if m.group(3) == "+" else -1) if m else None


rows = [r for r in csv.DictReader(open(feats), delimiter="\t") if r["tool"] in tools]
recs = {t: {r["accession"] for r in rows if r["tool"] == t and r["kind"] in ("PCG", "rRNA", "tRNA") and r["pred"]} for t in tools}
common = set.intersection(*recs.values())
sup = defaultdict(lambda: [0, 0])
for r in csv.DictReader(open(proj), delimiter="\t"):
    if r["rank"] != "0" and r["d5"] != "" and r["d3"] != "":
        sup[r["key"]][0] += r["d5"] == "0"
        sup[r["key"]][1] += r["d3"] == "0"
print(f"records finished by every tool: {len(common)}")
print("tool\ttruth rRNAs\tdetected\tboth exact\twithin 3\twithin 10\tconsensus ends exact\tconsensus ends within 3")
for t in tools:
    c = Counter()
    for r in rows:
        if r["tool"] != t or r["kind"] != "rRNA" or not r["truth"] or r["accession"] not in common:
            continue
        c["n"] += 1
        found = r["class"] in ("exact", "le3", "gt3")
        c["det"] += found
        if found:
            d = int(r["delta"])
            c["ex"] += d == 0
            c["le3"] += d <= 3
            c["le10"] += d <= 10
        tr = loc(r["truth"])
        key = f"{r['accession']}|{r['name']}|{tr[0]}"
        s = sup.get(key)
        if not s:
            continue
        p = loc(r["pred"]) if found else None
        if p is None:
            d5 = d3 = None
        elif tr[2] == 1:
            d5, d3 = p[0] - tr[0], p[1] - tr[1]
        else:
            d5, d3 = tr[1] - p[1], tr[0] - p[0]
        for ok, dv in ((s[0] >= 2, d5), (s[1] >= 2, d3)):
            if ok:
                c["cn"] += 1
                c["cex"] += dv == 0
                c["cle3"] += dv is not None and abs(dv) <= 3
    n, cn = c["n"], max(c["cn"], 1)
    print(f"{t}\t{n}\t{100 * c['det'] / n:.1f}%\t{100 * c['ex'] / n:.1f}%\t{100 * c['le3'] / n:.1f}%\t{100 * c['le10'] / n:.1f}%"
          f"\t{100 * c['cex'] / cn:.1f}% (n={c['cn']})\t{100 * c['cle3'] / cn:.1f}%")
