#!/usr/bin/env python3
"""Run only MitoM's tRNA step on benchmark records, as mito_core runs it, and write runs/<tool>/<acc>/final.gbf (tRNAs
only, NCBI orientation) + status.json, so mito_bench_eval.py scores them like any tool.

  trna_ml     mito_core.annotate_trnas_ml: MitoM-ML (trna_ml.py) with ARWEN / ARAGORN evidence, lineage from the record
  trna_mitfi  mito_core.annotate_trnas: cmsearch with the MiTFi models + MiTFi (the programs MITOS2 runs)

usage: trna_only.py --set SET.tsv --root ROOT --tools trna_ml,trna_mitfi [--jobs 24] [--threads 1] [--acc A,B]
"""
import argparse
import csv
import json
import logging
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(HERE))
import mito_core  # noqa: E402
import mito_env  # noqa: E402
from Bio import SeqIO  # noqa: E402

log = logging.getLogger("trna_only")


def tools_for(method: str) -> dict:
    t = {"cmsearch": mito_core.program("cmsearch")}
    if method == "trna_ml":
        for x in ("cmalign", "arwen", "aragorn"):
            t[x] = mito_env.find(x, None) or ""
    else:
        t["java"] = mito_core.program("java")
        t["mitfi"] = mito_core.mitfi_files(None)
    return t


def one(job) -> str:
    method, row, root, threads = job
    acc = row["accession"]
    out = Path(root) / "runs" / method / acc
    if (out / "status.json").exists():
        return f"{method} {acc} skipped"
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    try:
        rec = SeqIO.read(str(Path(root) / "gb" / f"{acc}.gb"), "genbank")
        genome, code = str(rec.seq).upper(), int(row["genetic_code"] or 5)
        n, circular = len(genome), rec.annotations.get("topology", "circular") != "linear"
        refdir = Path(mito_env.mitos_refdir(None, "refseq89m")) / "refseq89m"
        tools = tools_for(method)
        if method == "trna_ml":
            tools["refdir"] = str(refdir)
            lineage = ";".join(rec.annotations.get("taxonomy", []))
            feats = mito_core.annotate_trnas_ml(genome, n, code, threads, tools, lineage, circular, log)
        else:
            searched = genome + genome[:min(mito_core.ORIGIN_EXT, n)] if circular else genome
            with tempfile.TemporaryDirectory(prefix="trna_only.") as td:
                gfa = Path(td) / "g.fa"
                gfa.write_text(f">g\n{searched}\n")
                feats = mito_core.annotate_trnas(searched, gfa, refdir, code, threads, tools, Path(td), log, n,
                                                 row.get("phylum") or None)
        gb = mito_core.genbank_record(genome, feats, code, acc.split(".")[0][:16], "circular" if circular else "linear")
        SeqIO.write(gb, str(out / "final.gbf"), "genbank")
        st = {"tool": method, "accession": acc, "status": "ok", "seconds": round(time.time() - t0, 1), "n_trna": len(feats)}
    except Exception as e:                      # one failing record must not stop the others
        st = {"tool": method, "accession": acc, "status": "failed", "seconds": round(time.time() - t0, 1),
              "error": f"{type(e).__name__}: {e}"[:300]}
    (out / "status.json").write_text(json.dumps(st, indent=1))
    return f"{method} {acc} {st['status']} {st['seconds']}s"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", required=True)
    p.add_argument("--root", required=True)
    p.add_argument("--tools", default="trna_ml,trna_mitfi")
    p.add_argument("--jobs", type=int, default=24)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--acc", default="")
    a = p.parse_args(argv)
    rows = list(csv.DictReader(open(a.set), delimiter="\t"))
    if a.acc:
        keep = set(a.acc.split(","))
        rows = [r for r in rows if r["accession"] in keep]
    jobs = [(m, r, a.root, a.threads) for m in a.tools.split(",") for r in rows]
    with ProcessPoolExecutor(a.jobs) as ex:
        for f in as_completed([ex.submit(one, j) for j in jobs]):
            print(f.result(), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
