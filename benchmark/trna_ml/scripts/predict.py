#!/usr/bin/env python3
"""Turn CM hits into tRNA predictions.

v0 (rule): loci = same-strand hits overlapping >= 50 % of the shorter; the locus takes the family of its
best hit (lowest E-value over all models, pan or clade); a locus is kept when that E-value <= --evalue.
Loci of different families overlapping > 50 % of the shorter: the weaker is dropped.
--use pan|clade|all restricts which models count (pan = MiTFi-like baseline with the same machinery).

Usage: predict.py <hits_dir> <out.tsv> [--evalue 1e-3] [--use all]
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import read_tsv, write_tsv  # noqa: E402

AA = {f: f[0] for f in ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L1", "L2", "M", "N", "P", "Q", "R", "S1", "S2",
                        "T", "V", "W", "Y"]}


def overlap(a, b):
    return min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1


def loci(hits):
    hits = sorted(hits, key=lambda h: h["evalue"])
    out = []
    for h in hits:
        for L in out:
            r = L[0]
            if r["strand"] == h["strand"]:
                o = overlap(r, h)
                if o >= 0.5 * min(r["end"] - r["start"] + 1, h["end"] - h["start"] + 1):
                    L.append(h)
                    break
        else:
            out.append([h])
    return out


def load(f, use, code=None):
    from learn import LS_FAM, ac_aa
    hs = []
    for h in read_tsv(f):
        a = ac_aa(h.get("anticodon", ""), (code or {}).get(h["acc"], "5"))
        if a and a in "LS" and h["family"][0] == a and LS_FAM.get(h["anticodon"], h["family"]) != h["family"]:
            h["family"] = LS_FAM[h["anticodon"]]       # L1/L2, S1/S2 from the anticodon, not the model
        if use == "pan" and h["clade"] != "pan":
            continue
        if use == "clade" and h["clade"] == "pan":
            continue
        for k in ("start", "end", "strand", "depth"):
            h[k] = int(h[k])
        h["evalue"] = float(h["evalue"])
        h["bits"] = float(h["bits"])
        hs.append(h)
    return hs


def predict(hits, evalue):
    preds = []
    for L in loci(hits):
        best = L[0]
        if best["evalue"] > evalue:
            continue
        preds.append({"acc": best["acc"], "aa": AA[best["family"]], "family": best["family"], "start": best["start"],
                      "end": best["end"], "strand": best["strand"], "evalue": best["evalue"], "bits": best["bits"],
                      "model": best["model"], "n_models": len({h["model"] for h in L})})
    preds.sort(key=lambda p: p["evalue"])
    keep = []
    for p in preds:
        if any(q["strand"] == p["strand"] and overlap(p, q) > 0.5 * min(p["end"] - p["start"], q["end"] - q["start"])
               for q in keep):
            continue
        keep.append(p)
    return sorted(keep, key=lambda p: p["start"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("hits")
    ap.add_argument("out")
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--use", default="all", choices=["all", "pan", "clade"])
    ap.add_argument("--dataset", default="", help="dataset dir (genetic code per record, for anticodon naming)")
    a = ap.parse_args()
    code = {r["acc"]: r["code"] for r in read_tsv(a.dataset + "/records.tsv")} if a.dataset else {}
    rows = []
    for f in sorted(Path(a.hits).glob("*.hits.tsv")):
        rows += predict(load(f, a.use, code), a.evalue)
    write_tsv(a.out, rows, ["acc", "aa", "family", "start", "end", "strand", "evalue", "bits", "model", "n_models"])


if __name__ == "__main__":
    main()
