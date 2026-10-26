#!/usr/bin/env python3
"""Single-tool tRNA baselines as benchmark runs: ARWEN and ARAGORN from the trna_dump.py dumps (the predictions MitoM-ML
used as evidence), tRNAscan-SE 2 run here (-O organellar mode, mito_trna.trnascan). Writes runs/<tool>/<acc>/final.gbf.

usage: ext_tools.py --set SET.tsv --root ROOT --tools arwen,aragorn,trnascan [--jobs 12]
"""
import argparse
import csv
import gzip
import json
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(HERE))
import mito_core  # noqa: E402
import mito_env  # noqa: E402
import mito_trna  # noqa: E402
import trna_ml  # noqa: E402
from Bio import SeqIO  # noqa: E402

SRC = {"arwen": "ARWEN", "aragorn": "ARAGORN"}


def one(job):
    tool, row, root = job
    acc = row["accession"]
    out = Path(root) / "runs" / tool / acc
    if (out / "status.json").exists():
        return f"{tool} {acc} skipped"
    t0 = time.time()
    rec = SeqIO.read(str(Path(root) / "gb" / f"{acc}.gb"), "genbank")
    genome, code, n = str(rec.seq).upper(), int(row["genetic_code"] or 5), len(rec)
    circular = rec.annotations.get("topology", "circular") != "linear"
    try:
        if tool in SRC:
            d = json.load(gzip.open(Path(root) / "cands" / f"{acc}.json.gz", "rt"))
            preds = [mito_trna.Pred(**e) for e in d["ext"] if e["source"] == SRC[tool]]
        else:
            with tempfile.TemporaryDirectory(prefix="ext_tools.") as td:
                fa = Path(td) / "g.fa"
                fa.write_text(f">g\n{genome}\n")
                preds = mito_trna.unwrap(mito_trna.trnascan(fa, code, {"tRNAscan-SE": mito_env.find("tRNAscan-SE", None)},
                                                            Path(td), None), n)
        feats = []
        for p in preds:
            fam = trna_ml.EXT_LS.get(p.anticodon, p.aa) if p.aa in "LS" else p.aa
            if not fam or (p.end > n and not circular):
                continue
            feats.append(mito_core.Feature("tRNA", mito_core.trna_name(fam, p.anticodon or "", ""), p.start - 1, p.end, p.strand,
                                           p.source, anticodon=p.anticodon or ""))
        out.mkdir(parents=True, exist_ok=True)
        SeqIO.write(mito_core.genbank_record(genome, sorted(feats, key=lambda f: f.start), code, acc.split(".")[0][:16],
                                             "circular" if circular else "linear"), str(out / "final.gbf"), "genbank")
        st = {"tool": tool, "accession": acc, "status": "ok", "seconds": round(time.time() - t0, 1), "n_trna": len(feats)}
    except Exception as e:                       # one failing record must not stop the others
        out.mkdir(parents=True, exist_ok=True)
        st = {"tool": tool, "accession": acc, "status": "failed", "seconds": round(time.time() - t0, 1),
              "error": f"{type(e).__name__}: {e}"[:300]}
    (out / "status.json").write_text(json.dumps(st, indent=1))
    return f"{tool} {acc} {st['status']}"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", required=True)
    p.add_argument("--root", required=True)
    p.add_argument("--tools", default="arwen,aragorn,trnascan")
    p.add_argument("--jobs", type=int, default=12)
    a = p.parse_args(argv)
    rows = list(csv.DictReader(open(a.set), delimiter="\t"))
    jobs = [(t, r, a.root) for t in a.tools.split(",") for r in rows]
    with ProcessPoolExecutor(a.jobs) as ex:
        for f in as_completed([ex.submit(one, j) for j in jobs]):
            print(f.result(), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
