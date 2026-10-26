#!/usr/bin/env python3
"""Third scan pass, as MITOS2's sensitive pass: cmsearch --max with the pan CM of every family that still has no hit
at E <= 1e-3 (any model) after scan.py --rescue. Adds the new hits (mode 'max_pan', anticodon read as in scan.py) to a
copy of the hits files.

Usage: rescue_max.py <dataset_dir> <hits_in_dir> <hits_out_dir> [--jobs 12]
"""
import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import FAMS, INF, PAN, read_tsv, sh, write_tsv  # noqa: E402
from scan import FIELDS, RESCUE_E, WRAP, add_anticodons, clen_of, parse_tbl  # noqa: E402


def one(job):
    ds, f, odir = job
    acc = f.name.replace(".hits.tsv", "")
    out = Path(odir) / f.name
    if out.exists():
        return
    rows = read_tsv(f)
    for r in rows:
        for k in ("start", "end", "strand"):
            r[k] = int(r[k])
    good = {r["family"] for r in rows if float(r["evalue"]) <= RESCUE_E}
    miss = [x for x in FAMS if x not in good]
    if miss:
        seq = "".join((Path(ds) / "fasta" / f"{acc}.fa").read_text().split("\n")[1:]).strip()
        tmp = Path(odir) / f"{acc}.tmp"
        tmp.mkdir(exist_ok=True)
        (tmp / "q.fa").write_text(f">{acc}\n{seq + seq[:WRAP]}\n")
        info = {}
        with open(tmp / "m.cm", "w") as fh:
            for x in miss:
                p = PAN / f"trn{x}.cm"
                fh.write(p.read_text())
                nm = next(line.split()[1] for line in open(p) if line.startswith("NAME"))
                info[nm] = ("pan", 0, clen_of(p), x)
        sh([str(INF / "cmsearch"), "--cpu", "1", "--tblout", str(tmp / "h.tbl"), "-E", "10", "--notextw", "--max",
            str(tmp / "m.cm"), str(tmp / "q.fa")])
        seen = {(r["model"], r["start"], r["end"]) for r in rows}
        new = [r for r in parse_tbl(tmp / "h.tbl", info, acc, len(seq), "max_pan") if (r["model"], r["start"], r["end"]) not in seen]
        if new:
            add_anticodons(new, seq + seq[:WRAP], tmp)
            rows += new
        for x in tmp.iterdir():
            x.unlink()
        tmp.rmdir()
    write_tsv(out, rows, FIELDS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("hits_in")
    ap.add_argument("hits_out")
    ap.add_argument("--jobs", type=int, default=12)
    a = ap.parse_args()
    Path(a.hits_out).mkdir(parents=True, exist_ok=True)
    jobs = [(a.dataset, f, a.hits_out) for f in sorted(Path(a.hits_in).glob("*.hits.tsv"))]
    with ProcessPoolExecutor(a.jobs) as ex:
        for i, _ in enumerate(ex.map(one, jobs), 1):
            if i % 100 == 0:
                print(i, flush=True)
    print("done", len(jobs))


if __name__ == "__main__":
    main()
