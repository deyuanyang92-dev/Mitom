"""Rerun only the intron step of `mitom run` (mito_intron.py, then mito_intron_ref.py, as mitom_steps.run_intron_fix)
on the pre-intron records of mitom_prein runs (runs/mitom_prein/<acc>/raw/work/04.final_gb/), leave-one-out
(MITOM_EXCLUDE_REF=<acc>: the record out of the site library), for working on the intron code.

Usage: python3 intron_only.py <root> <label> [--jobs 32] [--denovo_only] [--env KEY=VALUE ...] [--denovo_arg=--x ...]
Writes <root>/intron_only/<label>/<acc>.gb (scored with benchmark/eval_introns.py --pred <root>/intron_only/<label>)
and <root>/intron_only/<label>/<acc>.report.tsv (the de novo intron_report.tsv).
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

M = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(M))
from mitom_steps import intron_denovo_cmd, intron_ref_cmd  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument("root")
p.add_argument("label")
p.add_argument("--jobs", type=int, default=32)
p.add_argument("--denovo_only", action="store_true")
p.add_argument("--env", action="append", default=[])
p.add_argument("--denovo_arg", action="append", default=[], help="extra mito_intron.py argument (repeat)")
a = p.parse_args()
root = Path(a.root)
out = root / "intron_only" / a.label
out.mkdir(parents=True, exist_ok=True)
rows = list(csv.DictReader(open(root / "set.tsv"), delimiter="\t"))


def one(r):
    acc = r["accession"]
    src = root / "runs" / "mitom_prein" / acc / "raw" / "work" / "04.final_gb"
    if not src.is_dir():
        return acc, "no input"
    env = {**os.environ, "MITOM_EXCLUDE_REF": acc, "OMP_NUM_THREADS": "1",
           **dict(x.split("=", 1) for x in a.env)}
    with tempfile.TemporaryDirectory() as td:
        d = Path(td) / "gb"
        shutil.copytree(src, d)
        code = int(r["genetic_code"])
        steps = [intron_denovo_cmd(d, d, code, 1, 1) + a.denovo_arg] + ([] if a.denovo_only else [intron_ref_cmd(d, d, code, 1)])
        log = []
        for cmd in steps:
            pr = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=3600)
            log.append(pr.stderr[-2000:])
            if pr.returncode:
                (out / f"{acc}.log").write_text("\n".join(log))
                return acc, f"exit {pr.returncode}"
        gbs = [g for g in d.glob("*.gbf") if g.stat().st_size]
        if gbs:
            shutil.copy(max(gbs, key=lambda g: g.stat().st_size), out / f"{acc}.gb")
        if (d / "intron_report.tsv").exists():
            shutil.copy(d / "intron_report.tsv", out / f"{acc}.report.tsv")
    return acc, "ok"


with ThreadPoolExecutor(a.jobs) as ex:
    res = list(ex.map(one, rows))
bad = [x for x in res if x[1] != "ok"]
print(f"{len(res) - len(bad)} ok, {len(bad)} failed: {bad[:10]}", file=sys.stderr)
