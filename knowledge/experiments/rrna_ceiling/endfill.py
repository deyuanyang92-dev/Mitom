"""Per rRNA end: nearest neighbour (kind, same/opposite strand), gap from MitoM's pre-R1 end, and where NCBI put the
end: at the neighbour (filled, <= 3 nt), at the pre-R1 end (<= 3 nt), or elsewhere.
usage: endfill.py SET.tsv ROOT TOOL SPLIT   (-> result_endfill_*.txt)"""
import csv, os, sys
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1cap import norm_features, match, Feat, pre_r1, model_cols, loc, SeqIO, KIND

set_tsv, root, tool, split = sys.argv[1:5]
meta = [r for r in csv.DictReader(open(set_tsv), delimiter="\t") if split == "all" or r["split"] == split]
stat = defaultdict(Counter)
for m in meta:
    acc = m["accession"]
    run = os.path.join(root, "runs", tool, acc)
    fin = os.path.join(run, "final.gbf")
    if not os.path.exists(fin):
        continue
    try:
        trec = SeqIO.read(os.path.join(root, "gb", acc + ".gb"), "genbank")
        precs = list(SeqIO.parse(fin, "genbank"))
    except Exception:
        continue
    if len(precs) != 1 or len(precs[0]) != len(trec):
        continue
    code = int(m["genetic_code"] or 5)
    T = [f for f in norm_features(trec, code) if f.kind == "rRNA"]
    P = norm_features(precs[0], code)
    n = len(trec)
    shifts = pre_r1(os.path.join(run, "raw", "reports", "curation_report.tsv"), None)
    used = Counter()
    preds, pre = [], {}
    for f in P:
        if f.kind != "rRNA":
            continue
        lst = shifts.get(f.name, [])
        i = used[f.name]
        d5, d3 = lst[i] if i < len(lst) else (0, 0)
        used[f.name] += 1
        s0, e0 = (f.start - d5, f.end - d3) if f.strand > 0 else (f.start + d3, f.end + d5)
        preds.append(f)
        pre[id(f)] = (s0 + 1, e0)       # 1-based inclusive
    for row in match(T, preds):
        if not row["truth"] or row["class"] not in ("exact", "le3", "gt3"):
            continue
        tr, p = loc(row["truth"]), loc(row["pred"])
        f = next(x for x in preds if (x.start, x.end, x.strand) == (p[0], p[1], p[2]))
        s, e = pre[id(f)]
        others = [g for g in P if g is not f and g.kind in ("PCG", "tRNA", "rRNA")]
        nk = [(g.start + 1 + k, g.end + k, KIND[g.kind], g.strand) for g in others for k in (-n, 0, n)]
        # trimming as R1 (only overlaps), to get the post-trim span before any fill
        neigh = [(a, b) for a, b, _, _ in nk]
        changed = True
        while changed:
            changed = False
            for ns, ne in sorted(neigh):
                if ne < s or ns > e or (ns <= s and ne >= e):
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
        left = max([b for a, b in neigh if b < s] or [0])
        right = min([a for a, b in neigh if a > e] or [10 ** 9])
        for side, gap, pos_fill, pos_pre, tpos, nb in (
                ("L", s - left - 1, left + 1, s, tr[0] + 1, [(k, st) for a, b, k, st in nk if b == left]),
                ("R", right - e - 1, right - 1, e, tr[1], [(k, st) for a, b, k, st in nk if a == right])):
            if gap <= 0 or not nb:
                continue
            k = "t" if any(x[0] == "t" for x in nb) else nb[0][0]
            same = any(x[1] == f.strand for x in nb if x[0] == k)
            gb = "1-10" if gap <= 10 else "11-30" if gap <= 30 else "31-100" if gap <= 100 else "101-300" if gap <= 300 else ">300"
            if abs(tpos - pos_fill) <= 3:
                cls = "ncbi_at_neighbour"
            elif abs(tpos - pos_pre) <= 3:
                cls = "ncbi_at_model_end"
            else:
                cls = "other"
            stat[(k + ("_same" if same else "_opp"), gb)][cls] += 1
print("neighbour\tgap\tn\tNCBI at neighbour (<=3)\tNCBI at model end (<=3)\tother")
order = ["1-10", "11-30", "31-100", "101-300", ">300"]
for key in sorted(stat, key=lambda k: (k[0], order.index(k[1]))):
    c = stat[key]
    t = sum(c.values())
    print(f"{key[0]}\t{key[1]}\t{t}\t{100*c['ncbi_at_neighbour']/t:.0f}%\t{100*c['ncbi_at_model_end']/t:.0f}%\t{100*c['other']/t:.0f}%")
