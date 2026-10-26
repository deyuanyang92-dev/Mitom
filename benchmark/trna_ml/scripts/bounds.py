#!/usr/bin/env python3
"""Learned tRNA boundaries.

For each called tRNA, candidate (start, end) pairs are formed from the ends reported by every source at the
locus (pan CM, clade CMs, ARWEN, ARAGORN) and their +-1 shifts. Each pair is scored by a gradient-boosting
model trained on the classifier-training genomes (label: both ends equal RefSeq) from: how many sources and
which kinds report each end, their best bit score, whether the CM alignment reached the model's first / last
column, shifts, length, and acceptor-stem pairing (7 bp, discriminator 0-2 nt) computed on the sequence.

Usage:
  bounds.py pairs <features.tsv> <dataset_dir> <out_pairs.tsv> [--positives_only]
  bounds.py train <pairs.tsv> <model.pkl>
  bounds.py apply <model.pkl> <pairs.tsv> <pred_in.tsv> <pred_out.tsv>
"""
import argparse
import pickle   # models are produced and read locally by this pipeline only
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import rc, read_tsv, write_tsv  # noqa: E402

PAIR = {("G", "C"), ("C", "G"), ("A", "T"), ("T", "A"), ("G", "T"), ("T", "G")}
WC = {("G", "C"), ("C", "G"), ("A", "T"), ("T", "A")}
F = ["s_n", "s_bits", "s_pan", "s_clade", "s_arwen", "s_aragorn", "s_full", "s_shift", "s_off",
     "e_n", "e_bits", "e_pan", "e_clade", "e_arwen", "e_aragorn", "e_full", "e_shift", "e_off",
     "len", "acc0", "acc1", "acc2", "wc1", "first_nt_T", "last_nt_A"]


def stem(seq5, seq3):
    """Pairs and Watson-Crick pairs between 7 nt 5' and the reversed 7 nt 3' strand."""
    p = w = 0
    for a, b in zip(seq5, reversed(seq3)):
        p += (a, b) in PAIR
        w += (a, b) in WC
    return p, w


def oriented(gen, s, e, strand):
    L = len(gen)
    g2 = gen + gen
    sub = g2[s - 1:e] if e <= 2 * L else ""
    return sub if strand == 1 else rc(sub)


def pairs_for(row, gen):
    strand = int(row["strand"])
    opts = []
    for o in row["opts"].split("|"):
        m, clade, depth, s, e, bits, mf, mt, cl = o.split(":")
        opts.append({"m": m, "clade": clade, "s": int(s), "e": int(e), "bits": float(bits),
                     "full5": int(mf) == 1, "full3": int(mt) == int(cl)})
    # 5' end of the tRNA = start on + strand, end on - strand
    S, E = defaultdict(list), defaultdict(list)
    for o in opts:
        S[o["s"]].append(o)
        E[o["e"]].append(o)
    med_s = statistics.median(o["s"] for o in opts)
    med_e = statistics.median(o["e"] for o in opts)
    cs = {s + d: (s, d) for s in S for d in (-1, 0, 1)}
    ce = {e + d: (e, d) for e in E for d in (-1, 0, 1)}
    for s, (s0, ds) in list(cs.items()):
        if s in S:
            cs[s] = (s, 0)
    for e, (e0, de) in list(ce.items()):
        if e in E:
            ce[e] = (e, 0)
    out = []
    for s, (s0, ds) in cs.items():
        for e, (e0, de) in ce.items():
            n = e - s + 1
            if not 40 <= n <= 110:
                continue
            seq = oriented(gen, s, e, strand)
            if len(seq) != n:
                continue
            accs = [stem(seq[:7], seq[-7 - d:len(seq) - d])[0] for d in (0, 1, 2)]
            wc1 = stem(seq[:7], seq[-8:-1])[1]
            so, eo = S[s0], E[e0]
            r = {"key": f"{row['acc']}|{row['locus']}|{row['family']}", "s": s, "e": e,
                 "s_n": len(so), "s_bits": max(o["bits"] for o in so),
                 "s_pan": int(any(o["clade"] == "pan" for o in so)),
                 "s_clade": int(any(o["clade"] not in ("pan", "ext") for o in so)),
                 "s_arwen": int(any(o["m"] == "ARWEN" for o in so)), "s_aragorn": int(any(o["m"] == "ARAGORN" for o in so)),
                 "s_full": int(any(o["full5" if strand == 1 else "full3"] for o in so if o["clade"] != "ext")),
                 "s_shift": ds, "s_off": s - med_s,
                 "e_n": len(eo), "e_bits": max(o["bits"] for o in eo),
                 "e_pan": int(any(o["clade"] == "pan" for o in eo)),
                 "e_clade": int(any(o["clade"] not in ("pan", "ext") for o in eo)),
                 "e_arwen": int(any(o["m"] == "ARWEN" for o in eo)), "e_aragorn": int(any(o["m"] == "ARAGORN" for o in eo)),
                 "e_full": int(any(o["full3" if strand == 1 else "full5"] for o in eo if o["clade"] != "ext")),
                 "e_shift": de, "e_off": e - med_e, "len": n, "acc0": accs[0], "acc1": accs[1], "acc2": accs[2],
                 "wc1": wc1, "first_nt_T": int(seq[0] == "T"), "last_nt_A": int(seq[-1] == "A")}
            if row.get("truth_se"):
                ts, te = map(int, row["truth_se"].split(":"))
                r["label"] = int(s == ts and e == te)
            else:
                r["label"] = 0
            out.append(r)
    return out


def make_pairs(ft, ds, out, positives_only, keys_from=None):
    keys = {f"{r['acc']}|{r['locus']}|{r['family']}" for r in read_tsv(keys_from)} if keys_from else None
    gen = {}
    rows = []
    for row in read_tsv(ft):
        if positives_only and row["label"] != "1":
            continue
        if keys is not None and f"{row['acc']}|{row['locus']}|{row['family']}" not in keys:
            continue
        a = row["acc"]
        if a not in gen:
            gen[a] = "".join((Path(ds) / "fasta" / f"{a}.fa").read_text().split("\n")[1:]).strip()
        rows += pairs_for(row, gen[a])
    write_tsv(out, rows, ["key", "s", "e", "label"] + F)


def train(pf, model):
    from sklearn.ensemble import HistGradientBoostingClassifier
    rows = read_tsv(pf)
    clf = HistGradientBoostingClassifier(max_iter=500, learning_rate=0.05, max_leaf_nodes=31, random_state=0)
    clf.fit([[float(r[k]) for k in F] for r in rows], [int(r["label"]) for r in rows])
    pickle.dump(clf, open(model, "wb"))
    print("trained", len(rows), "pairs")


def apply(model, pf, pin, pout):
    clf = pickle.load(open(model, "rb"))
    rows = read_tsv(pf)
    best = {}
    if rows:
        p = clf.predict_proba([[float(r[k]) for k in F] for r in rows])[:, 1]
        for r, q in zip(rows, p):
            if r["key"] not in best or q > best[r["key"]][0]:
                best[r["key"]] = (q, int(r["s"]), int(r["e"]))
    out = []
    for r in read_tsv(pin):
        k = f"{r['acc']}|{r['locus']}|{r['family']}"
        if k in best:
            r["start"], r["end"] = best[k][1], best[k][2]
            r["bound_p"] = round(best[k][0], 4)
        out.append(r)
    write_tsv(pout, out)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a1 = sub.add_parser("pairs")
    a1.add_argument("features")
    a1.add_argument("dataset")
    a1.add_argument("out")
    a1.add_argument("--positives_only", action="store_true")
    a1.add_argument("--keys", default="", help="only candidates called in this prediction TSV (needs locus column)")
    a2 = sub.add_parser("train")
    a2.add_argument("pairs")
    a2.add_argument("model")
    a3 = sub.add_parser("apply")
    a3.add_argument("model")
    a3.add_argument("pairs")
    a3.add_argument("pred_in")
    a3.add_argument("pred_out")
    a = ap.parse_args()
    if a.cmd == "pairs":
        make_pairs(a.features, a.dataset, a.out, a.positives_only, a.keys or None)
    elif a.cmd == "train":
        train(a.pairs, a.model)
    else:
        apply(a.model, a.pairs, a.pred_in, a.pred_out)


if __name__ == "__main__":
    main()
