#!/usr/bin/env python3
"""MitoM-ML tRNA step of mito_core (annotate_trnas_ml: ARWEN / ARAGORN evidence, trna_ml.predict) on benchmark records,
keeping everything the caller saw, so caller / end rules can be replayed offline without cmsearch:
  runs/trna_ml/<acc>/final.gbf + status.json   the production calls (as trna_only.py writes them)
  cands/<acc>.json.gz                          hits, ARWEN / ARAGORN predictions, candidates with p, final calls

usage: trna_dump.py --set SET.tsv --root ROOT [--jobs 16] [--threads 2] [--acc A,B]
"""
import argparse
import csv
import gzip
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
import mito_trna  # noqa: E402
import trna_ml  # noqa: E402
from Bio import SeqIO  # noqa: E402

log = logging.getLogger("trna_dump")
TOOL = "trna_ml"


def ml_with_dump(genome, n, code, threads, tools, lineage, circular):
    """mito_core.annotate_trnas_ml, with trna_ml.predict unrolled to keep its intermediates."""
    from threadpoolctl import threadpool_limits
    ext = []
    with tempfile.TemporaryDirectory(prefix="trna_dump.") as td:
        fa = Path(td) / "g.fa"
        fa.write_text(f">g\n{genome}\n")
        for t in ("arwen", "aragorn"):
            if tools.get(t):
                ext += mito_trna.unwrap(mito_trna.TOOLS[t](fa, code, {t: tools[t]}, Path(td), None), n)
    pan = Path(tools["refdir"]) / "ncRNA"
    with threadpool_limits(limits=max(1, threads)):
        clf, bclf = trna_ml.classifiers()
        with tempfile.TemporaryDirectory() as td:
            hits = trna_ml.scan(genome, lineage, pan, trna_ml.clade_dir(), tools["cmsearch"], tools["cmalign"],
                                Path(td), threads)
        cands = trna_ml.candidates(hits, trna_ml.ext_hits(ext), code)
        calls = trna_ml.call(cands, clf, trna_ml.THRESHOLD)
        trna_ml.refine_ends(calls, genome.upper(), bclf)
    feats = []
    for c in sorted(calls, key=lambda c: c["start"]):
        if c["end"] > n and not circular:
            continue
        feats.append(mito_core.Feature("tRNA", mito_core.trna_name(c["family"], c["anticodon"], ""), c["start"] - 1,
                                       c["end"], c["strand"], f"MitoM-ML p {c['p']:.3f} ends p {c.get('bound_p', 0):.3f}",
                                       anticodon=c["anticodon"], score=c["p"], model="MitoM-ML"))
    keep = lambda d: {k: v for k, v in d.items() if k != "_opts"}       # noqa: E731
    dump = {"lineage": lineage, "n": n, "circular": circular,
            "ext": [p._asdict() if hasattr(p, "_asdict") else vars(p) for p in ext],
            "hits": hits,
            "cands": [dict(keep(c), opts=[h for h in c["_opts"]]) for c in cands],
            "calls": [keep(c) for c in calls]}
    return feats, dump


def one(job) -> str:
    row, root, threads = job
    acc = row["accession"]
    out = Path(root) / "runs" / TOOL / acc
    dump_p = Path(root) / "cands" / f"{acc}.json.gz"
    if (out / "status.json").exists() and dump_p.exists():
        return f"{acc} skipped"
    out.mkdir(parents=True, exist_ok=True)
    dump_p.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    try:
        rec = SeqIO.read(str(Path(root) / "gb" / f"{acc}.gb"), "genbank")
        genome, code = str(rec.seq).upper(), int(row["genetic_code"] or 5)
        n, circular = len(genome), rec.annotations.get("topology", "circular") != "linear"
        tools = {"cmsearch": mito_core.program("cmsearch"),
                 "refdir": str(Path(mito_env.mitos_refdir(None, "refseq89m")) / "refseq89m")}
        for x in ("cmalign", "arwen", "aragorn"):
            tools[x] = mito_env.find(x, None) or ""
        feats, dump = ml_with_dump(genome, n, code, threads, tools, ";".join(rec.annotations.get("taxonomy", [])), circular)
        gb = mito_core.genbank_record(genome, feats, code, acc.split(".")[0][:16], "circular" if circular else "linear")
        SeqIO.write(gb, str(out / "final.gbf"), "genbank")
        with gzip.open(dump_p, "wt") as fh:
            json.dump(dump, fh)
        st = {"tool": TOOL, "accession": acc, "status": "ok", "seconds": round(time.time() - t0, 1), "n_trna": len(feats)}
    except Exception as e:                      # one failing record must not stop the others
        st = {"tool": TOOL, "accession": acc, "status": "failed", "seconds": round(time.time() - t0, 1),
              "error": f"{type(e).__name__}: {e}"[:300]}
    (out / "status.json").write_text(json.dumps(st, indent=1))
    return f"{acc} {st['status']} {st['seconds']}s"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", required=True)
    p.add_argument("--root", required=True)
    p.add_argument("--jobs", type=int, default=16)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--acc", default="")
    a = p.parse_args(argv)
    rows = list(csv.DictReader(open(a.set), delimiter="\t"))
    if a.acc:
        keep = set(a.acc.split(","))
        rows = [r for r in rows if r["accession"] in keep]
    with ProcessPoolExecutor(a.jobs) as ex:
        for f in as_completed([ex.submit(one, (r, a.root, a.threads)) for r in rows]):
            print(f.result(), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
