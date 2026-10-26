#!/usr/bin/env python3
"""Run the existing tools (via MitoM's mito_trna.py) on test records, one process per record, and convert the
review GenBank into one prediction TSV per tool (+ mito_trna 'selected').

Usage: run_baselines.py <dataset_dir> <out_dir> --tools arwen,aragorn,trnascan [--accs file] [--jobs 16]
"""
import argparse
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import L_S, read_tsv, sh, write_tsv  # noqa: E402

MT = "/mnt/n/yzz-分类工作流/temp_scripts/Mitoz-annotate/mito_trna.py"
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}


def run_one(job):
    acc, fa, code, odir, tools = job
    od = Path(odir) / acc
    if (od / "trna_report.tsv").exists():
        return acc
    sh(["python3", MT, "-i", fa, "-o", str(od), "--genetic_code", code or "5", "--tools", tools,
        "--mitos_mode", "full", "--mitos_bin", "/root/miniconda3/envs/mitos/bin", "--mitos_refdir", "/opt/mitos_refs",
        "--bin", "/opt/micromamba/envs/invertrnax/bin"], timeout=7200)
    return acc


def feats(gbf, acc, L):
    out = []
    for r in SeqIO.parse(str(gbf), "genbank"):
        for f in r.features:
            if f.type != "tRNA":
                continue
            note = " ".join(f.qualifiers.get("note", []))
            src = re.search(r"source: ([^;]+);", note)
            m = re.search(r"tRNA-([A-Z][a-z]{2})", f.qualifiers.get("product", [""])[0])
            ac = re.search(r"anticodon ([acgt]{3})", note)
            sc = re.search(r"score ([-0-9.e]+);", note)
            if not m or m.group(1) not in AA3:
                continue
            aa = AA3[m.group(1)]
            ac = ac.group(1) if ac else ""
            parts = sorted(f.location.parts, key=lambda p: p.start)
            s, e = int(f.location.start) + 1, int(f.location.end)
            if len(parts) > 1 and int(parts[0].start) == 0 and int(parts[-1].end) == L:   # origin-spanning
                s, e = int(parts[-1].start) + 1, L + int(parts[0].end)
            out.append({"acc": acc, "aa": aa, "family": L_S.get(ac, "") if aa in "LS" else aa, "start": s, "end": e,
                        "strand": f.location.strand or 1, "source": src.group(1) if src else "selected",
                        "score": sc.group(1) if sc else "", "anticodon": ac})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("out")
    ap.add_argument("--tools", required=True)
    ap.add_argument("--accs", default="")
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--convert_only", action="store_true")
    a = ap.parse_args()
    ds, out = Path(a.dataset), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    want = set(Path(a.accs).read_text().split()) if a.accs else None
    recs = [r for r in read_tsv(ds / "records.tsv") if (r["acc"] in want if want else r["split"] == "test")]
    if not a.convert_only:
        jobs = [(r["acc"], str(ds / "fasta" / f"{r['acc']}.fa"), r["code"], str(out / "runs"), a.tools) for r in recs]
        with ProcessPoolExecutor(a.jobs) as ex:
            for i, _ in enumerate(ex.map(run_one, jobs), 1):
                if i % 50 == 0:
                    print(i, flush=True)
    by = {}
    done = []
    for r in recs:
        od = out / "runs" / r["acc"]
        if not (od / "trna_report.tsv").exists():
            continue
        done.append(r["acc"])
        L = int(r["length"])
        for p in feats(od / "review" / f"{r['acc']}.gbf", r["acc"], L):
            by.setdefault(p["source"], []).append(p)
        for p in feats(od / "selected" / f"{r['acc']}.gbf", r["acc"], L):
            p["source"] = "selected"
            by.setdefault("mito_trna_selected", []).append(p)
    for src, rows in by.items():
        write_tsv(out / f"{src}.tsv", rows, ["acc", "aa", "family", "start", "end", "strand", "source", "score",
                                             "anticodon"])
    (out / "done_accs.txt").write_text("\n".join(done) + "\n")
    print("done", len(done), {k: len(v) for k, v in by.items()})


if __name__ == "__main__":
    main()
