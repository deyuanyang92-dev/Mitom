#!/usr/bin/env python3
"""tRNA ensembles replayed offline from trna_dump.py dumps (MitoM-ML candidates, ARWEN / ARAGORN predictions) and the
MiTFi calls of trna_only.py; writes runs/<variant>/<acc>/final.gbf + status.json for mito_bench_eval.py.

  thr:T     MitoM-ML with caller threshold T (0.5 shipped), then the MiTFi rescue (mito_core.rescue_trnas)
  vote:K    every call of MitoM-ML (p >= 0.5), MiTFi, ARWEN, ARAGORN; a (family, strand, locus) is kept when >= K of the
            four agree (same family and strand, >= 50 % of the shorter); ends from MitoM-ML, else MiTFi, ARWEN, ARAGORN;
            same-strand calls overlapping > 50 % keep the better supported
  aa        thr:0.5 (= MitoM-ML + MiTFi rescue), then families still missing taken from loci where ARWEN and ARAGORN
            agree, when they overlap no call by > 10 nt (high-sensitivity callers corrected by agreement)
  corr:P    MitoM-ML calls that no other caller supports (MiTFi, ARWEN, ARAGORN: same family and strand, >= 50 %) are
            dropped when p < P, then the MiTFi rescue (weak calls need a second opinion)

usage: ensemble_sim.py --set SET.tsv --root ROOT --variants thr:0.4,vote:2,aa
"""
import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_core  # noqa: E402
import trna_ml  # noqa: E402
from rescue_sim import mitfi_calls, rescue_mitfi, overlap  # noqa: E402
from Bio import SeqIO  # noqa: E402

ORDER = ["ML", "MiTFi", "ARWEN", "ARAGORN"]


def ext_calls(d, src):
    out = []
    for e in d["ext"]:
        if e["source"] != src:
            continue
        fam = trna_ml.EXT_LS.get(e["anticodon"], "") if e["aa"] in "LS" else e["aa"]
        if fam:
            out.append({"family": fam, "start": e["start"], "end": e["end"], "strand": e["strand"],
                        "anticodon": e["anticodon"], "src": src})
    return out


def agree(a, b, n):
    return a["family"] == b["family"] and a["strand"] == b["strand"] and \
        overlap(a, b, n) >= 0.5 * min(a["end"] - a["start"] + 1, b["end"] - b["start"] + 1)


def ml_calls(d, genome, thr, clf, bclf):
    if abs(thr - trna_ml.THRESHOLD) < 1e-9:
        return [dict(c, src="ML") for c in d["calls"]]
    cands = [dict(c, _opts=c["opts"]) for c in d["cands"]]
    calls = trna_ml.call(cands, clf, thr)
    trna_ml.refine_ends(calls, genome, bclf)
    return [dict(c, src="ML") for c in calls]


def vote(sources, k, n):
    """sources: {name: [calls]}; clusters of agreeing calls, kept when >= k sources."""
    clusters = []
    for name in ORDER:
        for c in sources.get(name, []):
            for cl in clusters:
                if agree(cl[0], c, n) and name not in {x["src"] for x in cl}:
                    cl.append(dict(c, src=name))
                    break
            else:
                clusters.append([dict(c, src=name)])
    kept = sorted((cl for cl in clusters if len({x["src"] for x in cl}) >= k), key=lambda cl: -len(cl))
    out = []
    for cl in kept:
        rep = min(cl, key=lambda x: ORDER.index(x["src"]))
        if any(x["strand"] == rep["strand"] and overlap(x, rep, n) > 0.5 * min(x["end"] - x["start"], rep["end"] - rep["start"])
               for x in out):
            continue
        out.append(rep)
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", required=True)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--variants", required=True)
    p.add_argument("--mitfi_runs", default="trna_mitfi")
    a = p.parse_args(argv)
    clf, bclf = trna_ml.classifiers()
    rows = list(csv.DictReader(open(a.set), delimiter="\t"))
    for v in a.variants.split(","):
        kind, *par = v.split(":")
        name = "ens_" + v.replace(":", "_").replace(".", "p")
        done = 0
        for r in rows:
            acc = r["accession"]
            dp, mg = a.root / "cands" / f"{acc}.json.gz", a.root / "runs" / a.mitfi_runs / acc / "final.gbf"
            if not dp.exists() or not mg.exists():
                continue
            d = json.load(gzip.open(dp, "rt"))
            rec = SeqIO.read(str(a.root / "gb" / f"{acc}.gb"), "genbank")
            genome, n, code = str(rec.seq).upper(), len(rec), int(r["genetic_code"] or 5)
            mitfi = [dict(c, src="MiTFi") for c in mitfi_calls(mg, n)]
            if kind == "thr":
                calls = rescue_mitfi(ml_calls(d, genome, float(par[0]), clf, bclf), mitfi, 0.1, 10, n)
            elif kind == "vote":
                calls = vote({"ML": ml_calls(d, genome, 0.5, clf, bclf), "MiTFi": mitfi,
                              "ARWEN": ext_calls(d, "ARWEN"), "ARAGORN": ext_calls(d, "ARAGORN")}, int(par[0]), n)
            elif kind == "corr":
                aw, ag = ext_calls(d, "ARWEN"), ext_calls(d, "ARAGORN")
                ml = [c for c in ml_calls(d, genome, 0.5, clf, bclf)
                      if c["p"] >= float(par[0]) or any(agree(c, x, n) for src in (mitfi, aw, ag) for x in src)]
                calls = rescue_mitfi(ml, mitfi, 0.1, 10, n)
            elif kind == "aa":
                calls = rescue_mitfi(ml_calls(d, genome, 0.5, clf, bclf), mitfi, 0.1, 10, n)
                have = {c["family"] for c in calls}
                ag = ext_calls(d, "ARAGORN")
                for c in ext_calls(d, "ARWEN"):
                    if c["family"] not in have and any(agree(c, x, n) for x in ag) and \
                            all(overlap(c, x, n) <= 10 for x in calls):
                        calls.append(c)
                        have.add(c["family"])
            else:
                raise SystemExit(f"unknown variant {v}")
            feats = [mito_core.Feature("tRNA", mito_core.trna_name(c["family"], c.get("anticodon", ""), ""), c["start"] - 1,
                                       c["end"], c["strand"], c.get("src", ""), anticodon=c.get("anticodon", ""))
                     for c in sorted(calls, key=lambda c: c["start"]) if d["circular"] or c["end"] <= n]
            out = a.root / "runs" / name / acc
            out.mkdir(parents=True, exist_ok=True)
            SeqIO.write(mito_core.genbank_record(genome, feats, code, acc.split(".")[0][:16],
                                                 "circular" if d["circular"] else "linear"), str(out / "final.gbf"), "genbank")
            (out / "status.json").write_text(json.dumps({"tool": name, "accession": acc, "status": "ok", "seconds": 0}))
            done += 1
        print(f"{name}: {done} records", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
