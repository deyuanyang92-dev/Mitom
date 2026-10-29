#!/usr/bin/env python3
"""Benchmark mito_intron_ref.py (intron-site transfer) on core outputs of NCBI records.

Inputs: a harness table (columns acc, set, group, code, topology, phylum, lineage, src = the NCBI GenBank record),
a folder of core outputs <core>/<acc>/<acc>.gbf (the annotation core's record, rotated as the pipeline rotates it,
without introns), and the site library (mito_intron_ref.py build).
Modes (the query's NCBI coordinates never reach the library in modes loo / lso / lgo):
  own  library = the query's own NCBI record only (`mitom run` on an annotated GenBank file): must be 100 % exact
  all  the whole library, nothing excluded (pipeline default; the query's own record is one of the references)
  loo  full library without the query record (and its RefSeq / GenBank copy)
  lso  without any record of the query's species
  lgo  without any record of the query's genus
Outputs: <out>/<mode>/<acc>.gbf (changed or copied records), <out>/<mode>/report.tsv (every reference junction
tried), <out>/<mode>/times.tsv, then benchmark/eval_introns.py on all records -> <out>/eval/<mode>_*.tsv.

Usage:
  python3 benchmark/bench_intron_ref.py --sets sets.tsv --core core_rot --lib models/intron_sites.tsv.gz \
      -o out --modes own,loo,lso,lgo --jobs 24 --bin ~/miniforge3/envs/mitotools/bin
"""
import argparse
import csv
import shutil
import subprocess
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_intron_ref as M  # noqa: E402
from Bio import SeqIO  # noqa: E402

_JOB = {}


def truth_of(r) -> Path:
    """The record that is the truth: --truth_gb/<acc>.gb (gold standard) when there, else the NCBI record."""
    gold = _JOB.get("truth_gb")
    if gold and (Path(gold) / f"{r['acc']}.gb").exists():
        return Path(gold) / f"{r['acc']}.gb"
    return Path(r.get("src") or r["truth"])


def split(tsv: Path, rows) -> str:
    """exact / le6 / gt6 / missed / fp per harness group (and group:phylum) of an eval_introns.py table."""
    import collections
    group = {r["acc"]: r.get("group", "all") for r in rows}
    phylum = {r["acc"]: r.get("phylum", "") for r in rows}
    tab = collections.defaultdict(collections.Counter)
    for x in csv.DictReader(open(tsv), delimiter="\t"):
        g = group.get(x["accession"], "?")
        for key in (g, f"{g}:{phylum.get(x['accession'], '')}", "all_pos" if g in ("dev", "new", "frag") else None,
                    f"class:{x['truth_class']}" if x.get("truth_class") else None):
            if key:
                tab[key][x["class"]] += 1
    cols = ["exact", "le6", "gt6", "missed", "fp"]
    lines = ["group\tn_truth\t" + "\t".join(cols)]
    for k in sorted(tab, key=lambda k: (k.split(":")[0], k)):
        c = tab[k]
        lines.append(f"{k}\t{c['exact'] + c['le6'] + c['gt6'] + c['missed']}\t" + "\t".join(str(c[x]) for x in cols))
    return "\n".join(lines) + "\n"


def one(args):
    mode, r = args
    lib_all, core, out, tools = _JOB["lib"], _JOB["core"], _JOB["out"], _JOB["tools"]
    acc = r["acc"]
    src = core / acc / f"{acc}.gbf"
    dest = out / mode / f"{acc}.gbf"
    if not src.exists():
        return mode, acc, "no core output", 0.0, []
    ncbi = SeqIO.read(str(truth_of(r)), "genbank")              # own sites / exclusions: the truth record
    exclude, own = None, None
    if mode == "own":
        lib, own = [], M.ref_genes(ncbi)
    elif mode == "all":                                         # pipeline default: whole library, nothing excluded
        lib = lib_all
    else:
        lib = lib_all
        org = ncbi.annotations.get("organism", "")
        exclude = M.Exclude({acc}, M.species_of(org) if mode in ("lso", "lgo") else "",
                            M.genus_of(org) if mode == "lgo" else "").expand(lib)
    rec = SeqIO.read(str(src), "genbank")
    t0 = time.time()
    try:
        lineage = None if _JOB["no_lineage"] else r["lineage"].split(";")
        models, rows = M.annotate_record(rec, lib, int(r["code"]), lineage=lineage, exclude=exclude, tools=tools,
                                         own=own)
        written = M.apply_models(rec, models, int(r["code"]))
        status = "ok"
    except Exception as exc:                       # report and keep the core record
        import traceback
        traceback.print_exc()
        models, rows, written, status = [], [], [], f"error {exc!r}"
    dt = time.time() - t0
    if written:
        SeqIO.write(rec, str(dest), "genbank")
    else:
        shutil.copy(src, dest)
    for m in models:
        rows.append({"record": acc, "gene": m.gene, "state": "model", "why": f"{m.confidence}: {m.evidence}",
                     "ref": ",".join(m.refs)})
    # evaluation copy: introns across the origin also written as intron features (the valid output leaves them to
    # the CDS join)
    if _JOB["no_calls"]:
        return mode, acc, status, dt, rows
    rec2 = SeqIO.read(str(src), "genbank")
    if models and M.apply_models(rec2, models, int(r["code"]), origin_introns=True):
        SeqIO.write(rec2, str(out / f"{mode}_calls" / f"{acc}.gbf"), "genbank")
    else:
        shutil.copy(src, out / f"{mode}_calls" / f"{acc}.gbf")
    return mode, acc, status, dt, rows


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sets", required=True)
    p.add_argument("--core", required=True, help="Folder <acc>/<acc>.gbf of core outputs (pipeline input of the step)")
    p.add_argument("--lib", default=str(M.DEFAULT_LIB))
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--modes", default="own,loo,lso,lgo")
    p.add_argument("--only", default="", help="comma-separated accessions")
    p.add_argument("--jobs", type=int, default=24)
    p.add_argument("--bin", default=None)
    p.add_argument("--eval", default=str(HERE / "eval_introns.py"), help="eval_introns.py (intron evaluator)")
    p.add_argument("--no_lineage", action="store_true", help="do not give the query's lineage (ranking by alignment "
                   "score, phylum of the closest library relative)")
    p.add_argument("--truth_gb", default="", help="folder <acc>.gb of the truth (gold standard): scored against and the "
                   "own-annotation sites taken from it; records not in it: the NCBI record (column src / truth)")
    p.add_argument("--no_calls", action="store_true", help="skip the evaluation copy (introns across the origin as "
                   "intron features)")
    a = p.parse_args(argv)
    rows = list(csv.DictReader(open(a.sets), delimiter="\t"))
    if a.only:
        keep = set(a.only.split(","))
        rows = [r for r in rows if r["acc"] in keep]
    out = Path(a.out)
    modes = a.modes.split(",")
    for m in modes:
        (out / m).mkdir(parents=True, exist_ok=True)
        (out / f"{m}_calls").mkdir(parents=True, exist_ok=True)
    _JOB.update(lib=M.load_library(Path(a.lib)), core=Path(a.core), out=out, tools=M.find_tools(a.bin),
                no_lineage=a.no_lineage, truth_gb=a.truth_gb, no_calls=a.no_calls)
    jobs = [(m, r) for m in modes for r in rows]
    res = []
    t0 = time.time()
    with ProcessPoolExecutor(a.jobs) as ex:
        for mode, acc, status, dt, rep in ex.map(one, jobs):
            res.append((mode, acc, status, dt, rep))
    print(f"[bench_intron_ref] {len(jobs)} runs in {time.time() - t0:.0f} s", file=sys.stderr)
    cols = M.REPORT_COLS
    for m in modes:
        with open(out / m / "report.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, cols, delimiter="\t", restval="", extrasaction="ignore")
            w.writeheader()
            for mode, acc, status, dt, rep in res:
                if mode == m:
                    w.writerows(rep)
        with open(out / m / "times.tsv", "w") as fh:
            for mode, acc, status, dt, rep in res:
                if mode == m:
                    fh.write(f"{acc}\t{status}\t{dt:.2f}\n")
    # evaluation with the shared evaluator, all records at once (truth folder of links to the NCBI records)
    truth = out / "truth_gb"
    truth.mkdir(exist_ok=True)
    for r in rows:
        link = truth / f"{r['acc']}.gb"
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(truth_of(r).resolve())
    ev = a.eval
    accs = ",".join(r["acc"] for r in rows)
    for m in modes:
        for label in (m,) if a.no_calls else (m, f"{m}_calls"):
            subprocess.run([sys.executable, ev, "--gb", str(truth), "--pred", str(out / label), "--label", label,
                            "--acc", accs, "-o", str(out / "eval")], check=False, capture_output=True)
            txt = split(out / "eval" / f"{label}_introns.tsv", rows)
            (out / "eval" / f"{label}_groups.tsv").write_text(txt)
            print(f"== {label}\n{txt}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
