#!/usr/bin/env python3
"""Real MiTFi baseline: the tRNA step of MitoM's mito_core.py (cmsearch fast pass + --max pass for missing models,
MiTFi jar for selection and anticodons - the procedure of MITOS2 runmitos), tRNAs only, one process per record.

Usage: run_mitfi.py <dataset_dir> <out_tsv> --accs file [--jobs 10]
"""
import argparse
import logging
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import read_tsv, write_tsv  # noqa: E402

MT = Path("/mnt/n/yzz-分类工作流/temp_scripts/Mitoz-annotate")
sys.path.insert(0, str(MT))
REF = Path("/opt/mitos_refs/refseq89m")
ORIGIN_EXT = 300


def one(job):
    import mito_core
    acc, fa, code, phylum, cache = job
    out = Path(cache) / f"{acc}.tsv"
    if out.exists():
        return read_tsv(out)
    seq = "".join(Path(fa).read_text().split("\n")[1:]).strip().upper()
    n = len(seq)
    tools = {"cmsearch": "/opt/micromamba/envs/invertrnax/bin/cmsearch",
             "java": mito_core.program("java"), "mitfi": mito_core.mitfi_files(None)}
    log = logging.getLogger("mitfi")
    with tempfile.TemporaryDirectory() as td:
        g = Path(td) / "g.fa"
        searched = seq + seq[:ORIGIN_EXT]
        g.write_text(f">g\n{searched}\n")
        feats = mito_core.annotate_trnas(searched, g, REF, int(code or 5), 1, tools, Path(td), log, n, phylum)
    rows = []
    for f in feats:
        name = f.name
        aa = name[3]
        fam = name[3:] if name[3:] in ("L1", "L2", "S1", "S2") else ("" if aa in "LS" else aa)
        rows.append({"acc": acc, "aa": aa, "family": fam, "start": f.start + 1, "end": f.end, "strand": f.strand,
                     "source": "MiTFi", "score": f.score, "anticodon": f.anticodon})
    write_tsv(out, rows, ["acc", "aa", "family", "start", "end", "strand", "source", "score", "anticodon"])
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("out")
    ap.add_argument("--accs", required=True)
    ap.add_argument("--jobs", type=int, default=10)
    a = ap.parse_args()
    want = set(Path(a.accs).read_text().split())
    recs = [r for r in read_tsv(Path(a.dataset) / "records.tsv") if r["acc"] in want]
    cache = Path(a.out + ".d")
    cache.mkdir(parents=True, exist_ok=True)
    jobs = [(r["acc"], str(Path(a.dataset) / "fasta" / f"{r['acc']}.fa"), r["code"], r["phylum"], str(cache)) for r in recs]
    rows = []
    with ProcessPoolExecutor(a.jobs) as ex:
        for i, rr in enumerate(ex.map(one, jobs), 1):
            rows += rr
            if i % 100 == 0:
                print(i, flush=True)
    write_tsv(a.out, rows, ["acc", "aa", "family", "start", "end", "strand", "source", "score", "anticodon"])
    print("done", len(recs), len(rows))


if __name__ == "__main__":
    main()
