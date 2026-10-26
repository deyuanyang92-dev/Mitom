#!/usr/bin/env python3
"""Score tRNA predictions against the RefSeq annotation.

Predictions TSV columns: acc, aa (one letter), family (A..Y, L1/L2/S1/S2 or ''), start, end, strand.
Match: same strand, overlap >= 50 % of the shorter feature, same amino acid; for Leu/Ser the family must
agree when both sides know it. Each truth and each prediction is matched at most once (best overlap first).

Usage: evaluate.py <dataset_dir> <pred.tsv> [<pred2.tsv> ...] [--by phylum|class] [--accs file]
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import read_tsv  # noqa: E402


def ov(a, b, L):
    def span(x):
        return int(x["start"]), int(x["end"])
    (s1, e1), (s2, e2) = span(a), span(b)
    best = 0
    for sh in (0, L, -L):
        best = max(best, min(e1, e2 + sh) - max(s1, s2 + sh) + 1)
    return best


def match(truth, pred, L):
    cand = []
    for i, t in enumerate(truth):
        for j, p in enumerate(pred):
            if int(t["strand"]) != int(p["strand"]):
                continue
            o = ov(t, p, L)
            short = min(int(t["end"]) - int(t["start"]), int(p["end"]) - int(p["start"])) + 1
            if o >= 0.5 * short:
                cand.append((o, i, j))
    cand.sort(reverse=True)
    ut, up, pairs = set(), set(), []
    for o, i, j in cand:
        if i in ut or j in up:
            continue
        ut.add(i)
        up.add(j)
        pairs.append((i, j))
    return pairs


def score(ds, pred_file, accs=None, by="phylum"):
    recs = {r["acc"]: r for r in read_tsv(ds / "records.tsv")}
    truth = defaultdict(list)
    for t in read_tsv(ds / "trnas_labeled.tsv"):
        if t["acc"] in recs and t["aa"]:
            truth[t["acc"]].append(t)
    pred = defaultdict(list)
    for p in read_tsv(pred_file):
        pred[p["acc"]].append(p)
    accs = sorted(accs or pred.keys())
    agg = defaultdict(lambda: defaultdict(int))
    for a in accs:
        L = int(recs[a]["length"])
        T, P = truth[a], pred.get(a, [])
        g = recs[a].get(by, "") or "?"
        for key in ("ALL", g):
            s = agg[key]
            s["rec"] += 1
            s["truth"] += len(T)
            s["pred"] += len(P)
        pairs = match(T, P, L)
        for i, j in pairs:
            t, p = T[i], P[j]
            name_ok = t["aa"] == p["aa"] and not (t["family"] and p.get("family") and t["aa"] in "LS"
                                                   and t["family"] != p["family"])
            ds5 = abs(int(t["start"]) - int(p["start"]))
            de3 = abs(int(t["end"]) - int(p["end"]))
            for key in ("ALL", g):
                s = agg[key]
                if not name_ok:
                    s["name_err"] += 1
                    continue
                s["tp"] += 1
                s["exact"] += ds5 == 0 and de3 == 0
                s["le3"] += ds5 <= 3 and de3 <= 3
        nm = {i for i, _ in pairs}
        for key in ("ALL", g):
            agg[key]["missed"] += len(T) - len(nm)
            agg[key]["fp"] += len(P) - len(pairs)
    return agg


def fmt(name, s):
    t = max(s["truth"], 1)
    return (f"{name:28s} rec {s['rec']:5d} truth {s['truth']:6d} recall {100 * s['tp'] / t:5.1f}% "
            f"exact {100 * s['exact'] / t:5.1f}% <=3nt {100 * s['le3'] / t:5.1f}% name_err {s['name_err']:4d} "
            f"missed {s['missed']:5d} FP {s['fp']:5d} FP/rec {s['fp'] / max(s['rec'], 1):.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("preds", nargs="+")
    ap.add_argument("--by", default="phylum")
    ap.add_argument("--accs", default="")
    ap.add_argument("--groups", action="store_true")
    a = ap.parse_args()
    ds = Path(a.dataset)
    accs = set(Path(a.accs).read_text().split()) if a.accs else None
    for pf in a.preds:
        agg = score(ds, pf, accs, a.by)
        print(fmt(Path(pf).stem, agg["ALL"]))
        if a.groups:
            for k in sorted(agg, key=lambda k: -agg[k]["truth"]):
                if k != "ALL":
                    print("   " + fmt(k, agg[k]))


if __name__ == "__main__":
    main()
