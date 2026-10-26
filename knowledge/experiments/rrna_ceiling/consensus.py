"""Per-end accuracy of a tool on 'consensus-supported' NCBI rRNA ends.

An NCBI end is consensus-supported when >= MIN of the K closest relatives (other species, proj_<mode>.tsv) put
their own NCBI end at exactly the same homologous position.  The selection uses only NCBI annotations,
never the tool output.  usage: consensus.py proj.tsv features.tsv TOOL [MIN]
"""
import csv, re, sys
from collections import defaultdict

proj, feats, tool = sys.argv[1:4]
MIN = int(sys.argv[4]) if len(sys.argv) > 4 else 2


def loc(t):
    m = re.match(r"(\d+)\.\.(\d+)([+-])", t or "")
    return (int(m.group(1)) - 1, int(m.group(2)), 1 if m.group(3) == "+" else -1) if m else None


tool_end = {}
for r in csv.DictReader(open(feats), delimiter="\t"):
    if r["tool"] != tool or r["kind"] != "rRNA" or not r["truth"]:
        continue
    t = loc(r["truth"])
    p = loc(r["pred"]) if r["class"] in ("exact", "le3", "gt3") else None
    key = f"{r['accession']}|{r['name']}|{t[0]}"
    if p is None:
        tool_end[key] = (None, None)
    elif t[2] == 1:
        tool_end[key] = (p[0] - t[0], p[1] - t[1])
    else:
        tool_end[key] = (t[1] - p[1], t[0] - p[0])
sup = defaultdict(lambda: [0, 0, 0])          # key -> [n refs, refs exact 5', refs exact 3']
for r in csv.DictReader(open(proj), delimiter="\t"):
    if r["rank"] == "0" or r["d5"] == "" or r["d3"] == "":
        continue
    s = sup[r["key"]]
    s[0] += 1
    s[1] += r["d5"] == "0"
    s[2] += r["d3"] == "0"
ends = {"5'": [0, 0, 0], "3'": [0, 0, 0]}
both = [0, 0, 0]
allk = [0, 0]
for key, (n, e5, e3) in sup.items():
    d = tool_end.get(key)
    if d is None:
        continue
    ok5 = e5 >= MIN
    ok3 = e3 >= MIN
    for name, ok, dv in (("5'", ok5, d[0]), ("3'", ok3, d[1])):
        if ok:
            ends[name][0] += 1
            ends[name][1] += dv == 0
            ends[name][2] += dv is not None and abs(dv) <= 3
    if ok5 and ok3:
        both[0] += 1
        both[1] += d[0] == 0 and d[1] == 0
        both[2] += d[0] is not None and d[1] is not None and abs(d[0]) <= 3 and abs(d[1]) <= 3
    allk[0] += 1
print(f"{tool}: rRNAs with relatives {allk[0]}; ends supported by >= {MIN} relatives:")
for name, (n, ex, le3) in ends.items():
    print(f"  {name} end: n={n} exact {ex} ({100*ex/max(n,1):.1f}%), within 3 {le3} ({100*le3/max(n,1):.1f}%)")
n, ex, le3 = both
print(f"  both ends supported: n={n} both exact {ex} ({100*ex/max(n,1):.1f}%), within 3 {le3} ({100*le3/max(n,1):.1f}%)")
