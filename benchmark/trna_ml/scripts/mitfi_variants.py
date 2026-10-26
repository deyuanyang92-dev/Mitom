#!/usr/bin/env python3
"""MiTFi (Juhling et al. 2012, the tRNA step of MITOS2) with other parameters or models, as mito_core.annotate_trnas runs it
otherwise: cmsearch fast pass per model, MiTFi at E <= FAST; models whose tRNA name is missing searched again with
cmsearch --max and MiTFi re-run at E <= SENS. Writes runs/<variant>/<acc>/final.gbf + status.json.

variant = comma-free key list joined by '+', e.g. ovl50+glocal:
  ovlN       MiTFi -overlap N (default 10)
  fastE      E threshold of the fast pass (default 0.01), e.g. fast0.001
  sensE      E threshold after the --max pass (default 0.1), e.g. sens1
  allmax     --max pass for every model, not only the missing ones
  glocal     cmsearch -g and MiTFi -method g (glocal; MITOS2 --locandgloc runs local and glocal)
  onlycut    MiTFi -onlycutoff (no top-22 round; MitoZ's setting)
  clade      the MitoM-ML clade CMs on the record's lineage are searched too and given to MiTFi with the 22 pan CMs
  default    as mito_core

usage: mitfi_variants.py --set SET.tsv --root ROOT --variants default,ovl50,clade [--jobs 10] [--threads 2]
"""
import argparse
import csv
import json
import logging
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(HERE))
import mito_core  # noqa: E402
import mito_env  # noqa: E402
import trna_ml  # noqa: E402
from Bio import SeqIO  # noqa: E402

log = logging.getLogger("mitfi_variants")


def parse_variant(v):
    o = {"ovl": 10, "fast": mito_core.CM_EVALUE, "sens": mito_core.CM_EVALUE_SENSITIVE, "allmax": False,
         "glocal": False, "onlycut": False, "clade": False}
    for k in v.split("+"):
        if k == "default":
            continue
        if k.startswith("ovl"):
            o["ovl"] = int(k[3:])
        elif k.startswith("fast"):
            o["fast"] = float(k[4:])
        elif k.startswith("sens"):
            o["sens"] = float(k[4:])
        elif k in ("allmax", "glocal", "onlycut", "clade"):
            o[k] = True
        else:
            raise SystemExit(f"unknown key {k}")
    return o


def cmsearch(cm, fa, out, sensitive, o):
    """mito_core.cmsearch with the variant's E thresholds and optional glocal mode."""
    cmd = [mito_core.program("cmsearch"), "--cpu", "1", "--notextw", "--smxsize", "80000"] + (["-g"] if o["glocal"] else []) + \
          (["--max", "-E", str(o["sens"])] if sensitive else ["--noF4b", "-E", str(o["fast"])])
    out.write_text(mito_core.run(cmd + [str(cm), str(fa)], log).stdout)
    return out


def mitfi_run(o, outs, fa, n, code, ev):
    jar, codes = mito_core.mitfi_files(None)
    cmd = [mito_core.program("java"), "-jar", str(jar), "-code", str(code), "-codefile", str(codes), "-structure", "-genome",
           str(fa), "-glength", str(n), "-evalue", str(ev), "-overlap", str(o["ovl"])] + \
        (["-method", "g"] if o["glocal"] else []) + (["-onlycutoff"] if o["onlycut"] else [])
    for x in outs:
        cmd += ["-in", str(x)]
    return mito_core.parse_mitfi(mito_core.run(cmd, log).stdout, log, n)


def one(job):
    variant, row, root, threads = job
    acc = row["accession"]
    out = Path(root) / "runs" / ("mf_" + variant.replace("+", "_").replace(".", "p")) / acc
    if (out / "status.json").exists():
        return f"{variant} {acc} skipped"
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    try:
        o = parse_variant(variant)
        rec = SeqIO.read(str(Path(root) / "gb" / f"{acc}.gb"), "genbank")
        genome, code = str(rec.seq).upper(), int(row["genetic_code"] or 5)
        n, circular = len(genome), rec.annotations.get("topology", "circular") != "linear"
        searched = genome + genome[:min(mito_core.ORIGIN_EXT, n)] if circular else genome
        pan = Path(mito_env.mitos_refdir(None, "refseq89m")) / "refseq89m" / "ncRNA"
        models = sorted(pan.glob("trn*.cm"))
        if o["clade"]:
            models += [p for p, clade, _ in trna_ml.model_paths(pan, trna_ml.clade_dir(), ";".join(rec.annotations.get("taxonomy", [])))
                       if clade != "pan"]
        with tempfile.TemporaryDirectory(prefix="mitfi_var.") as td:
            fa = Path(td) / "g.fa"
            fa.write_text(f">g\n{searched}\n")
            with ThreadPoolExecutor(max(1, threads)) as ex:
                outs = list(ex.map(lambda m: cmsearch(m, fa, Path(td) / f"{m.stem}.out", o["allmax"], o), models))
            feats = mitfi_run(o, outs, fa, n, code, o["sens"] if o["allmax"] else o["fast"])
            if not o["allmax"]:
                have = {f.name for f in feats}
                missing = [m for m in models if m.stem.split("__")[-1].replace("trn", "") not in {h[3:] for h in have}]
                if row.get("phylum", "").lower() in mito_core.NO_SENSITIVE_TRNA:
                    missing = []
                if missing:
                    with ThreadPoolExecutor(max(1, threads)) as ex:
                        list(ex.map(lambda m: cmsearch(m, fa, Path(td) / f"{m.stem}.out", True, o), missing))
                    feats = mitfi_run(o, outs, fa, n, code, o["sens"])
        gb = mito_core.genbank_record(genome, feats, code, acc.split(".")[0][:16], "circular" if circular else "linear")
        SeqIO.write(gb, str(out / "final.gbf"), "genbank")
        st = {"tool": out.parent.name, "accession": acc, "status": "ok", "seconds": round(time.time() - t0, 1), "n_trna": len(feats)}
    except Exception as e:                       # one failing record must not stop the others
        st = {"tool": out.parent.name, "accession": acc, "status": "failed", "seconds": round(time.time() - t0, 1),
              "error": f"{type(e).__name__}: {e}"[:300]}
    (out / "status.json").write_text(json.dumps(st, indent=1))
    return f"{variant} {acc} {st['status']} {st['seconds']}s"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", required=True)
    p.add_argument("--root", required=True)
    p.add_argument("--variants", required=True)
    p.add_argument("--jobs", type=int, default=10)
    p.add_argument("--threads", type=int, default=2)
    a = p.parse_args(argv)
    rows = list(csv.DictReader(open(a.set), delimiter="\t"))
    jobs = [(v, r, a.root, a.threads) for v in a.variants.split(",") for r in rows]
    with ProcessPoolExecutor(a.jobs) as ex:
        for f in as_completed([ex.submit(one, j) for j in jobs]):
            print(f.result(), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
