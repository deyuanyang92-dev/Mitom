#!/usr/bin/env python3
"""Cross-fitted MitoM-ML caller: the shipped caller table plus the candidates of the other half of the unseen-family test
set (splits/test_tune_dev.txt / test_tune_test.txt; label = same amino acid and strand as an NCBI tRNA, >= 50 % of it),
then the shipped ends and the MiTFi rescue; every record is called by a model that never saw its family.
Writes runs/cal_x/<acc>/final.gbf + status.json.

usage: caller_xfit.py --set SET.tsv --root ROOT [--name cal_x] [--weight 1.0]
"""
import argparse
import csv
import gzip
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_core  # noqa: E402
import trna_ml  # noqa: E402
from bounds_xfit import ncbi_trnas, ov  # noqa: E402
from rescue_sim import mitfi_calls, rescue_mitfi  # noqa: E402
from Bio import SeqIO  # noqa: E402


def hit_ends(opts):
    """The hit whose ends trna_ml.candidates gave the candidate (CM hits before ARWEN / ARAGORN, then lowest E)."""
    v = [h for h in opts if h["clade"] != "ext"] or opts
    return min(v, key=lambda h: (h["evalue"], h["model"] != "ARWEN"))


def label(c, truth, n):
    return int(any(t["strand"] == c["strand"] and t["aa"] == c["family"][0] and
                   ov(c["start"], c["end"], t["s"], t["e"], n) >= 0.5 * (t["e"] - t["s"] + 1) for t in truth))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", required=True)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--name", default="cal_x")
    p.add_argument("--weight", type=float, default=1.0, help="weight of the new rows")
    a = p.parse_args(argv)
    from sklearn.ensemble import HistGradientBoostingClassifier
    sp = Path(__file__).resolve().parents[1] / "splits"
    half = {x.strip(): "dev" for x in (sp / "test_tune_dev.txt").read_text().split()}
    half.update({x.strip(): "test" for x in (sp / "test_tune_test.txt").read_text().split()})
    recs, rows = {}, {"dev": [], "test": []}
    for r in csv.DictReader(open(a.set), delimiter="\t"):
        acc = r["accession"]
        dp = a.root / "cands" / f"{acc}.json.gz"
        if not dp.exists():
            continue
        rec = SeqIO.read(str(a.root / "gb" / f"{acc}.gb"), "genbank")
        d = json.load(gzip.open(dp, "rt"))
        truth = ncbi_trnas(rec)
        for c in d["cands"]:
            # the dumped start / end of a called candidate are its refined ends; the hit ends are what the caller saw
            best = hit_ends(c["opts"])
            rows[half[acc]].append(([float(c[k]) for k in trna_ml.FEATS],
                                    label(dict(c, start=best["start"], end=best["end"]), truth, d["n"])))
        recs[acc] = (r, rec, d)
    with gzip.open(trna_ml.MODELS / "caller_train.tsv.gz", "rt") as fh:
        ship = [([float(x[k]) for k in trna_ml.FEATS], int(x["label"]), float(x.get("weight") or 1))
                for x in csv.DictReader(fh, delimiter="\t")]
    print(f"candidates: dev {len(rows['dev'])} ({sum(y for _, y in rows['dev'])} +), test {len(rows['test'])} "
          f"({sum(y for _, y in rows['test'])} +), shipped {len(ship)}", flush=True)
    models = {}
    for h, o in (("dev", "test"), ("test", "dev")):
        X = [x for x, _, _ in ship] + [x for x, _ in rows[o]]
        y = [y for _, y, _ in ship] + [y for _, y in rows[o]]
        w = [w for _, _, w in ship] + [a.weight] * len(rows[o])
        clf = HistGradientBoostingClassifier(random_state=0, max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
                                             l2_regularization=1.0)
        models[h] = clf.fit(X, y, sample_weight=w)
    _, bclf = trna_ml.classifiers()
    for acc, (r, rec, d) in recs.items():
        genome, n = str(rec.seq).upper(), d["n"]
        cands = [dict(c, _opts=c["opts"]) for c in d["cands"]]
        for c in cands:                                   # back to the hit ends the caller saw
            best = hit_ends(c["_opts"])
            c["start"], c["end"] = best["start"], best["end"]
        calls = trna_ml.call(cands, models[half[acc]], trna_ml.THRESHOLD)
        trna_ml.refine_ends(calls, genome, bclf)
        mg = a.root / "runs" / "trna_mitfi" / acc / "final.gbf"
        if mg.exists():
            calls = rescue_mitfi(calls, mitfi_calls(mg, n), 0.1, 10, n)
        feats = [mito_core.Feature("tRNA", mito_core.trna_name(c["family"], c.get("anticodon", ""), ""), c["start"] - 1,
                                   c["end"], c["strand"], c.get("src", "MitoM-ML"), anticodon=c.get("anticodon", ""))
                 for c in sorted(calls, key=lambda c: c["start"]) if d["circular"] or c["end"] <= n]
        out = a.root / "runs" / a.name / acc
        out.mkdir(parents=True, exist_ok=True)
        SeqIO.write(mito_core.genbank_record(genome, feats, int(r["genetic_code"] or 5), acc.split(".")[0][:16],
                                             "circular" if d["circular"] else "linear"), str(out / "final.gbf"), "genbank")
        (out / "status.json").write_text(json.dumps({"tool": a.name, "accession": acc, "status": "ok", "seconds": 0}))
    print(f"{a.name}: {len(recs)} records", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
