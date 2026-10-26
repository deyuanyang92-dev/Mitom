#!/usr/bin/env python3
"""Cross-fitted boundary models for MitoM-ML tRNA calls (from trna_dump.py dumps), on top of the MiTFi rescue of
rescue_sim.py: the records are split in two halves by family (splits/test_tune_dev.txt, test_tune_test.txt); a model
trained on one half chooses the ends in the other, so every record is judged by a model that never saw its family.

  bnd_ship   the shipped model (models/trna_ml/bounds_train.tsv.gz), = MitoM-ML as it is
  bnd_x      cross-fit on the other half + the shipped table, the shipped features (BFEATS)
  bnd_xn     same, plus neighbour features: gap from each candidate end to the nearest other tRNA call (transcript
             direction, either strand; NaN in the shipped rows, which have no genome context)
Each writes runs/<variant>/<acc>/final.gbf (+ status.json): the variant's ends for the MitoM-ML calls, the MiTFi
rescue calls as in res_mitfi_0p1. Labels: a pair is positive when both ends equal the matched NCBI tRNA.

usage: bounds_xfit.py --set SET.tsv --root ROOT [--variants bnd_ship,bnd_x,bnd_xn]
"""
import argparse
import csv
import gzip
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_core  # noqa: E402
import trna_ml  # noqa: E402
from rescue_sim import mitfi_calls, rescue_mitfi  # noqa: E402
from Bio import SeqIO  # noqa: E402

NB = ["nb5", "nb3", "nb5_same", "nb3_same"]
CAP = 30


def ncbi_trnas(rec):
    n = len(rec)
    out = []
    for f in rec.features:
        if f.type != "tRNA":
            continue
        parts = f.location.parts
        if len(parts) == 2 and {int(p.start) for p in parts} & {0} and {int(p.end) for p in parts} & {n}:
            tail = next(p for p in parts if int(p.end) == n)
            head = next(p for p in parts if int(p.start) == 0)
            s, e = int(tail.start) + 1, n + int(head.end)
        else:
            s, e = int(f.location.start) + 1, int(f.location.end)
        txt = " ".join(f.qualifiers.get("product", []) + f.qualifiers.get("gene", []))
        aa = next((v for k, v in trna_ml_aa().items() if f"tRNA-{k}" in txt), "")
        out.append({"s": s, "e": e, "strand": f.location.strand, "aa": aa})
    return out


_AA = None


def trna_ml_aa():
    global _AA
    if _AA is None:
        _AA = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
               "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
               "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
    return _AA


def ov(a0, a1, b0, b1, n):
    return max(min(a1, b1 + k) - max(a0, b0 + k) + 1 for k in (-n, 0, n))


def truth_for(c, truth, n):
    best, bo = None, 0
    for t in truth:
        if t["strand"] != c["strand"] or t["aa"] != c["family"][0]:
            continue
        o = ov(c["start"], c["end"], t["s"], t["e"], n)
        if o >= 0.5 * (t["e"] - t["s"] + 1) and o > bo:
            best, bo = t, o
    return best


def neighbour_feats(s, e, strand, others, n):
    """Gap (nt; negative = overlap) from the pair's 5' / 3' end to the nearest other call in transcript direction."""
    down, up = [], []                    # (gap, same strand) of calls after e / before s in genome order
    for x in others:
        for k in (-n, 0, n):
            xs, xe = x["start"] + k, x["end"] + k
            if xs > s and xe > e:
                down.append((xs - e - 1, x["strand"] == strand))
            if xe < e and xs < s:
                up.append((s - xe - 1, x["strand"] == strand))
    d = min(down) if down else (CAP, False)
    u = min(up) if up else (CAP, False)
    g5, g3 = (u, d) if strand == 1 else (d, u)
    return {"nb5": min(g5[0], CAP), "nb3": min(g3[0], CAP), "nb5_same": int(g5[1]), "nb3_same": int(g3[1])}


def pairs_of(rec_d, genome):
    """[(call, [pair dict with BFEATS + NB])] for the MitoM-ML calls of one record."""
    n, out = rec_d["n"], []
    cands = {(c["locus"], c["family"]): c for c in rec_d["cands"]}
    for c in rec_d["calls"]:
        cand = cands[(c["locus"], c["family"])]
        pr = trna_ml.end_pairs(dict(cand, _opts=cand["opts"]), genome)
        others = [x for x in rec_d["calls"] if x is not c]
        for p in pr:
            p.update(neighbour_feats(p["s"], p["e"], c["strand"], others, n))
        out.append((c, pr))
    return out


def fit(rows, feats):
    from sklearn.ensemble import HistGradientBoostingClassifier
    clf = HistGradientBoostingClassifier(random_state=0, max_iter=500, learning_rate=0.05, max_leaf_nodes=31)
    clf.fit([[r.get(k, math.nan) for k in feats] for r in rows], [r["label"] for r in rows])
    return clf


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", required=True)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--variants", default="bnd_ship,bnd_x,bnd_xn")
    a = p.parse_args(argv)
    sp = Path(__file__).resolve().parents[1] / "splits"
    half = {x.strip(): "dev" for x in (sp / "test_tune_dev.txt").read_text().split()}
    half.update({x.strip(): "test" for x in (sp / "test_tune_test.txt").read_text().split()})
    rows = [r for r in csv.DictReader(open(a.set), delimiter="\t") if (a.root / "cands" / f"{r['accession']}.json.gz").exists()]
    recs = {}
    for r in rows:
        acc = r["accession"]
        rec = SeqIO.read(str(a.root / "gb" / f"{acc}.gb"), "genbank")
        d = json.load(gzip.open(a.root / "cands" / f"{acc}.json.gz", "rt"))
        genome = str(rec.seq).upper()
        recs[acc] = (r, rec, d, genome, pairs_of(d, genome), ncbi_trnas(rec))
    train = {"dev": [], "test": []}
    for acc, (r, rec, d, genome, prs, truth) in recs.items():
        for c, pr in prs:
            t = truth_for(c, truth, d["n"])
            if t is None:
                continue
            for q in pr:
                ok = any((q["s"] + k, q["e"] + k) == (t["s"], t["e"]) for k in (-d["n"], 0, d["n"]))
                train[half[acc]].append(dict(q, label=int(ok)))
    with gzip.open(trna_ml.MODELS / "bounds_train.tsv.gz", "rt") as fh:
        ship = [{k: (int(v) if k == "label" else float(v)) for k, v in x.items()} for x in csv.DictReader(fh, delimiter="\t")]
    print(f"pairs: dev {len(train['dev'])} ({sum(x['label'] for x in train['dev'])} +), test {len(train['test'])} "
          f"({sum(x['label'] for x in train['test'])} +), shipped {len(ship)}", flush=True)
    _, bship = trna_ml.classifiers()
    models = {}
    for v in a.variants.split(","):
        if v == "bnd_ship":
            models[v] = {"dev": (bship, trna_ml.BFEATS), "test": (bship, trna_ml.BFEATS)}
        else:
            feats = trna_ml.BFEATS + (NB if v == "bnd_xn" else [])
            models[v] = {h: (fit(ship + train[o], feats), feats) for h, o in (("dev", "test"), ("test", "dev"))}
        print(f"{v}: trained", flush=True)
    for v, mh in models.items():
        for acc, (r, rec, d, genome, prs, truth) in recs.items():
            clf, feats = mh[half[acc]]
            calls = []
            for c, pr in prs:
                c = dict(c)
                if pr:
                    pp = clf.predict_proba([[q.get(k, math.nan) for k in feats] for q in pr])[:, 1]
                    i = max(range(len(pr)), key=lambda k: pp[k])
                    c["start"], c["end"] = pr[i]["s"], pr[i]["e"]
                calls.append(c)
            mg = a.root / "runs" / "trna_mitfi" / acc / "final.gbf"
            n = d["n"]
            calls = rescue_mitfi(calls, mitfi_calls(mg, n), 0.1, 10, n) if mg.exists() else calls
            feats_out = [mito_core.Feature("tRNA", mito_core.trna_name(c["family"], c.get("anticodon", ""), ""), c["start"] - 1,
                                           c["end"], c["strand"], c.get("src", "MitoM-ML"), anticodon=c.get("anticodon", ""))
                         for c in sorted(calls, key=lambda c: c["start"]) if d["circular"] or c["end"] <= n]
            out = a.root / "runs" / v / acc
            out.mkdir(parents=True, exist_ok=True)
            SeqIO.write(mito_core.genbank_record(genome, feats_out, int(r["genetic_code"] or 5), acc.split(".")[0][:16],
                                                 "circular" if d["circular"] else "linear"), str(out / "final.gbf"), "genbank")
            (out / "status.json").write_text(json.dumps({"tool": v, "accession": acc, "status": "ok", "seconds": 0}))
        print(f"{v}: {len(recs)} records written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
