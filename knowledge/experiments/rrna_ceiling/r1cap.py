"""Replay MitoM's rRNA end rule (R1) with the fill capped by the model columns the core's span is missing.

R1 now trims an rRNA off its neighbours, then fills a gap of <= 300 nt to the next feature of any kind. Next to a
control region that fill swallows up to 300 nt of it (insect rrnS 5' ends: NC_066676.1, span model 21-768/776,
filled 300 nt to trnK, NCBI end = the span's end). Policies, per end, after the unchanged trimming:
  current   fill if gap <= 300
  capS      fill if gap <= missing + S (missing = model columns the span lacks at that end: hmm_from - 1 at the 5'
            end, LENG - hmm_to at the 3' end, from the core's note "model a-b/LENG"); otherwise keep the end
  capS+x    as capS, but a longer gap gets the missing columns only (extend by `missing`, at most to the neighbour)
  current+free_x  current, and an end without a neighbour within 300 nt is extended by `missing`
  fill:KINDS[:ss][+x]  fill (gap <= 300) only when the nearest neighbour at that end is of KINDS (t = tRNA, r = rRNA,
            p = protein gene), with :ss only when it is on the rRNA's strand (tRNA punctuation, Ojala et al. 1981: the
            rRNA is cut out of its own strand's transcript at the flanking tRNAs); +x: an end left unfilled gets the
            missing model columns (at most to the neighbour)
rRNAs without a model note (added by the rRNA HMM step) keep the current rule.
Pre-R1 ends come from reports/curation_report.tsv as in r1sim.py; every other feature is MitoM's final one.
Scores: all NCBI rRNAs (both ends exact / within 3 / within 10) and the consensus-supported ends of
final_metrics.py (proj_species.tsv: >= 2 of the 5 closest relatives put their NCBI end at the same position).
usage: r1cap.py SET.tsv ROOT TOOL proj_species.tsv SPLIT
"""
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mito_bench_eval import norm_features, match, Feat   # noqa: E402
from r1sim import pre_r1                                  # noqa: E402
from Bio import SeqIO                                      # noqa: E402

MAX_FILL = 300
GENE = {"s-rrna": "rrnS", "rrns": "rrnS", "l-rrna": "rrnL", "rrnl": "rrnL"}


def model_cols(raw):
    """gene -> list of (hmm_from, hmm_to, leng) in file order, from the core's rRNA notes (04.final_gb)."""
    out = defaultdict(list)
    for gb in glob.glob(os.path.join(raw, "work", "04.final_gb", "*.gbf")):
        for rec in SeqIO.parse(gb, "genbank"):
            for f in rec.features:
                if f.type != "rRNA":
                    continue
                g = GENE.get(f.qualifiers.get("gene", [""])[0].lower())
                m = re.search(r"model\s+(\d+)-(\d+)/(\d+)", " ".join(f.qualifiers.get("note", [])))
                if g:
                    out[g].append(tuple(map(int, m.groups())) if m else None)
    return out


KIND = {"tRNA": "t", "rRNA": "r", "PCG": "p"}


def r1(s, e, strand, others, n, policy, slack, cols):
    """s, e: 1-based inclusive pre-R1 span (e may exceed n). Returns the 0-based start, exclusive end."""
    neigh = [(f.start + 1 + k, f.end + k) for f in others for k in (-n, 0, n)]
    nk = [(f.start + 1 + k, f.end + k, KIND[f.kind], f.strand) for f in others for k in (-n, 0, n)]
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
    left = max([ne for ns, ne in neigh if ne < s] or [0])
    right = min([ns for ns, ne in neigh if ns > e] or [10 ** 9])
    gl, gr = s - left - 1, right - e - 1
    if policy == "current" or cols is None:
        if 0 < gl <= MAX_FILL:
            s = left + 1
        if 0 < gr <= MAX_FILL:
            e = right - 1
        return s - 1, e
    if policy.startswith("fill:"):
        kinds, ss, x = policy[5:].split(":")[0].rstrip("+x"), ":ss" in policy, policy.endswith("+x")
        ok_l = any(k in kinds and (not ss or st == strand) for ns, ne, k, st in nk if ne == left)
        ok_r = any(k in kinds and (not ss or st == strand) for ns, ne, k, st in nk if ns == right)
        ml = mr = 0
        if x and cols is not None:
            a, b, leng = cols
            m5, m3 = max(a - 1, 0), max(leng - b, 0)
            ml, mr = (m5, m3) if strand > 0 else (m3, m5)
        if 0 < gl <= MAX_FILL and ok_l:
            s = left + 1
        elif gl > 0:
            s -= min(ml, gl)
        if 0 < gr <= MAX_FILL and ok_r:
            e = right - 1
        elif gr > 0:
            e += min(mr, gr)
        return s - 1, e
    a, b, leng = cols
    m5, m3 = max(a - 1, 0), max(leng - b, 0)
    ml, mr = (m5, m3) if strand > 0 else (m3, m5)
    if policy == "freex":                     # current fill; a free end (no neighbour within MAX_FILL) gets the
        if 0 < gl <= MAX_FILL:                # missing model columns
            s = left + 1
        elif gl > MAX_FILL:
            s -= ml
        if 0 < gr <= MAX_FILL:
            e = right - 1
        elif gr > MAX_FILL:
            e += mr
        return s - 1, e
    if 0 < gl <= min(MAX_FILL, ml + slack):
        s = left + 1
    elif policy == "capx" and gl > 0 and ml > 0:
        s -= min(ml, gl)
    if 0 < gr <= min(MAX_FILL, mr + slack):
        e = right - 1
    elif policy == "capx" and gr > 0 and mr > 0:
        e += min(mr, gr)
    return s - 1, e


def loc(t):
    m = re.match(r"(\d+)\.\.(\d+)([+-])", t or "")
    return (int(m.group(1)) - 1, int(m.group(2)), 1 if m.group(3) == "+" else -1) if m else None


def main():
    set_tsv, root, tool, proj, split = sys.argv[1:6]
    meta = [r for r in csv.DictReader(open(set_tsv), delimiter="\t") if split == "all" or r["split"] == split]
    sup = defaultdict(lambda: [0, 0])
    for r in csv.DictReader(open(proj), delimiter="\t"):
        if r["rank"] != "0" and r["d5"] != "" and r["d3"] != "":
            sup[r["key"]][0] += r["d5"] == "0"
            sup[r["key"]][1] += r["d3"] == "0"
    variants = [("final", None, 0), ("current", "current", 0)] + \
               [(f"cap{s}", "cap", s) for s in (0, 10, 20, 30, 50, 100)] + \
               [(f"cap{s}+x", "capx", s) for s in (0, 10, 20, 30, 50, 100)] + [("current+free_x", "freex", 0)] + \
               [(p, p, 0) for p in ("fill:t", "fill:tr", "fill:t:ss", "fill:tr:ss", "fill:trp:ss", "fill:t:ss+x",
                                    "fill:tr:ss+x")]
    score = {v[0]: Counter() for v in variants}
    sanity, nocols = Counter(), Counter()
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
        cols = model_cols(os.path.join(run, "raw"))
        used = Counter()
        rr_pre = []
        for f in P:
            if f.kind != "rRNA":
                continue
            lst, cl = shifts.get(f.name, []), cols.get(f.name, [])
            i = used[f.name]
            d5, d3 = lst[i] if i < len(lst) else (0, 0)
            c = cl[i] if i < len(cl) else None
            nocols[c is None] += 1
            used[f.name] += 1
            if f.strand > 0:
                rr_pre.append((f.start - d5, f.end - d3, f.strand, f.name, f, c))
            else:
                rr_pre.append((f.start + d3, f.end + d5, f.strand, f.name, f, c))
        for vname, pol, slack in variants:
            if pol is None:
                preds = [f for f in P if f.kind == "rRNA"]
            else:
                preds = []
                for s0, e0, st, name, f, c in rr_pre:
                    others = [g for g in P if g is not f and g.kind in ("PCG", "tRNA", "rRNA")]
                    ns, ne = r1(s0 + 1, e0, st, others, n, pol, slack, c)
                    preds.append(Feat("rRNA", name, st, ns, ne, "", n))
                    if pol == "current":
                        sanity[(ns, ne) == (f.start, f.end)] += 1
            sc = score[vname]
            for row in match(T, preds):
                if not row["truth"]:
                    continue
                sc["n"] += 1
                found = row["class"] in ("exact", "le3", "gt3")
                sc["det"] += found
                if found:
                    d = int(row["delta"])
                    sc["ex"] += d == 0
                    sc["le3"] += d <= 3
                    sc["le10"] += d <= 10
                tr = loc(row["truth"])
                s = sup.get(f"{acc}|{row['name']}|{tr[0]}")
                if not s:
                    continue
                p = loc(row["pred"]) if found else None
                if p is None:
                    d5 = d3 = None
                elif tr[2] == 1:
                    d5, d3 = p[0] - tr[0], p[1] - tr[1]
                else:
                    d5, d3 = tr[1] - p[1], tr[0] - p[0]
                for ok, dv in ((s[0] >= 2, d5), (s[1] >= 2, d3)):
                    if ok:
                        sc["cn"] += 1
                        sc["cex"] += dv == 0
                        sc["cle3"] += dv is not None and abs(dv) <= 3
    print(f"{tool} records {len(meta)} split {split}; replay of the current rule reproduces MitoM's final rRNA: "
          f"{sanity[True]} of {sanity[True] + sanity[False]}; rRNAs without a model note: {nocols[True]} of {sum(nocols.values())}")
    print("variant\tn\tdetected\tboth exact\twithin 3\twithin 10\tconsensus ends\texact\twithin 3")
    for v, c in score.items():
        n, cn = c["n"], max(c["cn"], 1)
        print(f"{v}\t{n}\t{100 * c['det'] / n:.1f}%\t{100 * c['ex'] / n:.1f}%\t{100 * c['le3'] / n:.1f}%\t"
              f"{100 * c['le10'] / n:.1f}%\t{c['cn']}\t{100 * c['cex'] / cn:.1f}%\t{100 * c['cle3'] / cn:.1f}%")


if __name__ == "__main__":
    main()
