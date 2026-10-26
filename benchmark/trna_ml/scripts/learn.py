#!/usr/bin/env python3
"""Learned tRNA caller on top of the CM hits.

Candidates: for every locus (same-strand hits overlapping >= 50 %) and every tRNA family hitting it, one
candidate with features from the pan model and the clade models of that family, the competing families at
the locus, and the genome context (rank of the candidate among all candidates of its family, number of
families already called with high confidence). Label (training only): the candidate overlaps a RefSeq tRNA
of the same family by >= 50 %.

A gradient-boosting classifier is trained on cross-fitted hits of training genomes (scanned with models
built without the genome's own family) and applied to test hits. Decision: keep the best candidate per
locus with p >= --threshold; boundaries from the hit chosen by --boundary (best E-value by default).

Usage:
  learn.py features <hits_dir> <dataset_dir> <out.tsv>
  learn.py train <features.tsv> <model.pkl>
  learn.py apply <model.pkl> <features.tsv> <pred.tsv> [--threshold 0.5]
"""
import argparse
import math
import pickle
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import FAMS, read_tsv, write_tsv  # noqa: E402

FEATS = ["pan_bits", "pan_logE", "clade_bits", "clade_logE", "best_logE", "depth", "n_models", "n_clade_models",
         "len", "cov", "trunc", "delta_bits_next_fam", "fam_rank", "fam_n_cands", "fam_gap_to_best", "best_fam_at_locus",
         "rel_to_fam_best", "is_LS", "overlap_stronger", "n_strong_genome", "gc", "ac_known", "ac_match", "arwen", "aragorn",
         "only_ext", "overlap_stronger_opp", "opp_bits_ratio"]


def loc_overlap(a, b):
    return min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1


def loci(hits):
    hits = sorted(hits, key=lambda h: h["evalue"])
    out = []
    for h in hits:
        for L in out:
            r = L[0]
            if r["strand"] == h["strand"] and loc_overlap(r, h) >= 0.5 * min(r["end"] - r["start"] + 1,
                                                                             h["end"] - h["start"] + 1):
                L.append(h)
                break
        else:
            out.append([h])
    return out


def ac_aa(ac, code):
    """Amino acid (one letter) the anticodon decodes under the genetic code; '' if unreadable."""
    from Bio.Data.CodonTable import unambiguous_dna_by_id
    if len(ac) != 3 or set(ac) - set("acgt"):
        return ""
    codon = ac[::-1].translate(str.maketrans("acgt", "TGCA"))
    t = unambiguous_dna_by_id.get(int(code or 5), unambiguous_dna_by_id[5])
    return t.forward_table.get(codon, "*")


LS_FAM = {"tag": "L1", "taa": "L2", "tct": "S1", "gct": "S1", "tga": "S2", "cct": "S1", "tcc": "S1", "act": "S1"}


def lg(e):
    return math.log10(max(e, 1e-300))


EXT = ("ARWEN", "ARAGORN")


def ext_hits(ext_dir):
    """ARWEN / ARAGORN predictions as pseudo-hits (clade 'ext'), by accession."""
    by = defaultdict(list)
    if not ext_dir:
        return by
    for src in EXT:
        f = Path(ext_dir) / f"{src}.tsv"
        if not f.exists():
            continue
        for p in read_tsv(f):
            fam = p["family"] or ("" if p["aa"] in "LS" else p["aa"])
            if not fam:
                continue
            by[p["acc"]].append({"acc": p["acc"], "model": src, "clade": "ext", "depth": -1, "family": fam,
                                 "start": int(p["start"]), "end": int(p["end"]), "strand": int(p["strand"]),
                                 "evalue": 1.0, "bits": 0.0, "gc": 0.0, "mfrom": 1, "mto": 1, "clen": 1,
                                 "trunc": "no", "anticodon": "", "ac_known": 1, "ac_match": 1})
    return by


def features(hits_dir, ds, out, ext_dir=None):
    ext = ext_hits(ext_dir)
    truth = defaultdict(list)
    lab = Path(ds) / "trnas_labeled.tsv"
    for t in read_tsv(lab):
        truth[t["acc"]].append(t)
    code = {r["acc"]: r["code"] for r in read_tsv(Path(ds) / "records.tsv")}
    rows = []
    for f in sorted(Path(hits_dir).glob("*.hits.tsv")):
        hs = read_tsv(f)
        for h in hs:
            a = ac_aa(h.get("anticodon", ""), code.get(h["acc"], "5"))
            h["ac_known"] = int(bool(a))
            if a and h["family"][0] == a and a in "LS" and LS_FAM.get(h["anticodon"], h["family"]) != h["family"]:
                h["family"] = LS_FAM[h["anticodon"]]          # L1/L2, S1/S2 by anticodon, not by model
            h["ac_match"] = int(a == h["family"][0]) if a else 0
            for k in ("start", "end", "strand", "depth", "clen"):
                h[k] = int(h[k])
            for k in ("evalue", "bits", "gc"):
                h[k] = float(h[k])
            h["mfrom"], h["mto"] = int(h["mfrom"]), int(h["mto"])
        acc = f.name.replace(".hits.tsv", "")
        hs = hs + ext.get(acc, [])
        if not hs:
            continue
        cands = []
        for li, L in enumerate(loci(hs)):
            byfam = defaultdict(list)
            for h in L:
                byfam[h["family"]].append(h)
            fam_best = {fa: max(v, key=lambda h: h["bits"]) for fa, v in byfam.items()}
            order = sorted(fam_best, key=lambda fa: -fam_best[fa]["bits"])
            for fa, v0 in byfam.items():
                ex = {h["model"] for h in v0 if h["clade"] == "ext"}
                v = [h for h in v0 if h["clade"] != "ext"] or v0
                pan = [h for h in v if h["clade"] == "pan"]
                cl = [h for h in v if h["clade"] not in ("pan", "ext")]
                deep = max(cl, key=lambda h: (h["depth"], h["bits"])) if cl else None
                best = min(v, key=lambda h: (h["evalue"], h["model"] != "ARWEN"))
                others = [fam_best[o]["bits"] for o in order if o != fa]
                c = {"acc": acc, "locus": li, "family": fa, "strand": best["strand"], "start": best["start"],
                     "end": best["end"], "pan_start": pan[0]["start"] if pan else "", "pan_end": pan[0]["end"] if pan else "",
                     "deep_start": deep["start"] if deep else "", "deep_end": deep["end"] if deep else "",
                     "pan_bits": max((h["bits"] for h in pan), default=-20),
                     "pan_logE": min((lg(h["evalue"]) for h in pan), default=2),
                     "clade_bits": deep["bits"] if deep else -20, "clade_logE": lg(deep["evalue"]) if deep else 2,
                     "best_logE": lg(best["evalue"]), "depth": deep["depth"] if deep else 0, "n_models": len(v),
                     "n_clade_models": len(cl), "len": best["end"] - best["start"] + 1,
                     "cov": (best["mto"] - best["mfrom"] + 1) / max(best["clen"], 1), "trunc": int(best["trunc"] != "no"),
                     "delta_bits_next_fam": fam_best[fa]["bits"] - (max(others) if others else 0),
                     "best_fam_at_locus": int(order[0] == fa), "is_LS": int(fa[0] in "LS"), "gc": best["gc"],
                     "ac_known": best["ac_known"], "ac_match": best["ac_match"],
                     "arwen": int("ARWEN" in ex), "aragorn": int("ARAGORN" in ex), "only_ext": int(best["clade"] == "ext"),
                     "_bits": fam_best[fa]["bits"],
                     "opts": "|".join(f"{h['model']}:{h['clade']}:{h['depth']}:{h['start']}:{h['end']}:{h['bits']}:"
                                      f"{h['mfrom']}:{h['mto']}:{h['clen']}" for h in v0),
                     "truth_se": ""}
                cands.append(c)
        fam_c = defaultdict(list)
        for c in cands:
            fam_c[c["family"]].append(c)
        for fa, v in fam_c.items():
            v.sort(key=lambda c: -c["_bits"])
            for r, c in enumerate(v):
                c["fam_rank"] = r
                c["fam_n_cands"] = len(v)
                c["fam_gap_to_best"] = v[0]["_bits"] - c["_bits"]
                c["rel_to_fam_best"] = c["_bits"] / max(v[0]["_bits"], 1)
        strong = [c for c in cands if c["best_logE"] < -5 and c["best_fam_at_locus"]]
        n_strong = len({c["family"] for c in strong})
        for c in cands:
            c["n_strong_genome"] = n_strong
            c["overlap_stronger"] = max([loc_overlap(c, s) for s in strong if s is not c and s["strand"] == c["strand"]
                                         and s["_bits"] > c["_bits"] and s["locus"] != c["locus"]] or [0])
            opp = [s for s in strong if s["strand"] != c["strand"] and loc_overlap(c, s) > 0.5 * c["len"]]
            c["overlap_stronger_opp"] = max([loc_overlap(c, s) for s in opp] or [0])
            c["opp_bits_ratio"] = max([s["_bits"] for s in opp] or [0]) / max(c["_bits"], 1)
            T = truth.get(acc, [])
            for t in T:
                if int(t["strand"]) == c["strand"] and (t["family"] == c["family"] or (
                        not t["family"] and t["aa"] == c["family"][0])) and loc_overlap(
                        {"start": int(t["start"]), "end": int(t["end"])}, c) >= 0.5 * min(int(t["length"]), c["len"]):
                    c["truth_se"] = f"{t['start']}:{t['end']}"
            c["label"] = int(any(int(t["strand"]) == c["strand"] and (t["family"] == c["family"] or (
                not t["family"] and t["aa"] == c["family"][0])) and loc_overlap(
                {"start": int(t["start"]), "end": int(t["end"])}, c) >= 0.5 * min(int(t["length"]), c["len"]) for t in T))
        rows += cands
    write_tsv(out, rows, ["acc", "locus", "family", "strand", "start", "end", "pan_start", "pan_end", "deep_start",
                          "deep_end", "label", "opts", "truth_se"] + FEATS)


def X(rows):
    return [[float(r[k]) for k in FEATS] for r in rows]


def train(ft, model, neg="", neg_weight=0.1):
    """neg: candidates of shuffled genomes (all label 0), down-weighted: weight 1 removed every shuffled-genome
    call but also the weak real tRNAs of divergent genomes (Acanthocephala, Kinorhyncha); 0.1 kept 0 / 150
    calls on unseen shuffled genomes and 15 of 25 of those tRNAs (2026-10-26)."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    rows = read_tsv(ft)
    nrows = read_tsv(neg) if neg else []
    clf = HistGradientBoostingClassifier(max_iter=400, learning_rate=0.05, max_leaf_nodes=31, l2_regularization=1.0,
                                         random_state=0)
    clf.fit(X(rows + nrows), [int(r["label"]) for r in rows + nrows],
            sample_weight=[1.0] * len(rows) + [neg_weight] * len(nrows))
    pickle.dump(clf, open(model, "wb"))
    print("trained on", len(rows), "candidates,", sum(int(r["label"]) for r in rows), "positive")


def apply(model, ft, out, thr, boundary):
    clf = pickle.load(open(model, "rb"))
    rows = read_tsv(ft)
    p = clf.predict_proba(X(rows))[:, 1] if rows else []
    for r, q in zip(rows, p):
        r["p"] = q
    best = {}
    for r in rows:
        k = (r["acc"], r["locus"])
        if r["p"] >= thr and (k not in best or r["p"] > best[k]["p"]):
            best[k] = r
    keep = sorted(best.values(), key=lambda r: -r["p"])
    final = []
    for r in keep:                                    # different loci, heavy overlap: keep more probable
        s, e = (int(r[f"{boundary}_start"]), int(r[f"{boundary}_end"])) if boundary != "best" and r.get(
            f"{boundary}_start") else (int(r["start"]), int(r["end"]))
        q = {"acc": r["acc"], "locus": r["locus"], "aa": r["family"][0], "family": r["family"], "start": s, "end": e,
             "strand": int(r["strand"]), "p": round(float(r["p"]), 4)}
        if any(x["acc"] == q["acc"] and x["strand"] == q["strand"] and
               loc_overlap(x, q) > 0.5 * min(x["end"] - x["start"], q["end"] - q["start"]) for x in final):
            continue
        final.append(q)
    write_tsv(out, final, ["acc", "locus", "aa", "family", "start", "end", "strand", "p"])


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a1 = sub.add_parser("features")
    a1.add_argument("hits")
    a1.add_argument("dataset")
    a1.add_argument("out")
    a1.add_argument("--ext", default="", help="dir with ARWEN.tsv / ARAGORN.tsv predictions")
    a2 = sub.add_parser("train")
    a2.add_argument("features")
    a2.add_argument("model")
    a2.add_argument("--neg", default="", help="features of shuffled genomes (negatives)")
    a2.add_argument("--neg_weight", type=float, default=0.1)
    a3 = sub.add_parser("apply")
    a3.add_argument("model")
    a3.add_argument("features")
    a3.add_argument("out")
    a3.add_argument("--threshold", type=float, default=0.5)
    a3.add_argument("--boundary", default="best", choices=["best", "pan", "deep"])
    a = ap.parse_args()
    if a.cmd == "features":
        features(a.hits, a.dataset, a.out, a.ext)
    elif a.cmd == "train":
        train(a.features, a.model, a.neg, a.neg_weight)
    else:
        apply(a.model, a.features, a.out, a.threshold, a.boundary)


if __name__ == "__main__":
    assert FAMS
    main()
