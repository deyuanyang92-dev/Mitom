#!/usr/bin/env python3
"""Replay MitoM-ML on dumped hits (through trna_ml.drop_wrap_duplicates), optionally without some clade models (e.g.
the invertebrate-built Chordata models on vertebrate records), then the MiTFi rescue as in mito_core; writes
runs/<name>/<acc>/final.gbf + status.json (records already written are skipped).

usage: replay_noclade.py --set SET.tsv --root ROOT --drop Chordata[,Tunicata] --name ml_nochordata [--lineage Vertebrata]
       (--lineage: only records whose lineage has this token are replayed without the models; the others as dumped)
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
import mito_trna  # noqa: E402
import trna_ml  # noqa: E402
from rescue_sim import mitfi_calls, rescue_mitfi  # noqa: E402
from Bio import SeqIO  # noqa: E402


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", required=True)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--drop", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--lineage", default="")
    p.add_argument("--all", action="store_true", help="replay every record (hits as trna_ml.scan returns them now)")
    p.add_argument("--thr", type=float, default=None, help="caller threshold (default trna_ml.THRESHOLD)")
    p.add_argument("--rescue_e", type=float, default=0.1, help="E cut-off of the MiTFi calls that may be added")
    p.add_argument("--mitfi_runs", default="trna_mitfi", help="runs/<this>/<acc>/final.gbf = MiTFi calls for the rescue")
    a = p.parse_args(argv)
    drop = set(a.drop.split(","))
    clf, bclf = trna_ml.classifiers()
    done = 0
    for r in csv.DictReader(open(a.set), delimiter="\t"):
        acc = r["accession"]
        dp = a.root / "cands" / f"{acc}.json.gz"
        if not dp.exists() or (a.root / "runs" / a.name / acc / "status.json").exists():
            continue
        d = json.load(gzip.open(dp, "rt"))
        rec = SeqIO.read(str(a.root / "gb" / f"{acc}.gb"), "genbank")
        genome, n, code = str(rec.seq).upper(), len(rec), int(r["genetic_code"] or 5)
        if not a.lineage or a.lineage in d["lineage"].split(";") or a.all:
            hits = trna_ml.drop_wrap_duplicates([h for h in d["hits"] if h["clade"] not in drop], n)
            ext = [mito_trna.Pred(**e) for e in d["ext"]]
            cands = trna_ml.candidates(hits, trna_ml.ext_hits(ext), code)
            calls = trna_ml.call(cands, clf, trna_ml.THRESHOLD if a.thr is None else a.thr)
            trna_ml.refine_ends(calls, genome, bclf)
        else:
            calls = d["calls"]
        mg = a.root / "runs" / a.mitfi_runs / acc / "final.gbf"
        if mg.exists():
            calls = rescue_mitfi(calls, mitfi_calls(mg, n), a.rescue_e, 10, n)
        feats = [mito_core.Feature("tRNA", mito_core.trna_name(c["family"], c.get("anticodon", ""), ""), c["start"] - 1,
                                   c["end"], c["strand"], c.get("src", "MitoM-ML"), anticodon=c.get("anticodon", ""))
                 for c in sorted(calls, key=lambda c: c["start"]) if d["circular"] or c["end"] <= n]
        out = a.root / "runs" / a.name / acc
        out.mkdir(parents=True, exist_ok=True)
        SeqIO.write(mito_core.genbank_record(genome, feats, code, acc.split(".")[0][:16],
                                             "circular" if d["circular"] else "linear"), str(out / "final.gbf"), "genbank")
        (out / "status.json").write_text(json.dumps({"tool": a.name, "accession": acc, "status": "ok", "seconds": 0}))
        done += 1
    print(f"{a.name}: {done} records", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
