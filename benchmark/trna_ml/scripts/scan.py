#!/usr/bin/env python3
"""Scan genomes with the clade CMs on the genome's lineage path + the MITOS2 pan-metazoan CMs.

The sequence is scanned as circular (first WRAP nt appended); hits are reported in 1-based coordinates
on the original sequence (end may exceed L for origin-spanning hits).

Usage: scan.py <dataset_dir> <models_dir> <out_dir> [--split test] [--accs a,b] [--jobs 16] [--max]
Output: <out_dir>/<acc>.hits.tsv
"""
import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import FAMS, INF, PAN, anticodon_from_alignment, parse_sto, rc, read_tsv, sh, write_tsv  # noqa: E402

WRAP = 300
FIELDS = ["acc", "model", "clade", "depth", "family", "start", "end", "strand", "mfrom", "mto", "clen", "trunc",
          "bits", "evalue", "gc", "anticodon", "mode"]


def model_paths(models_dir, lineage, meta):
    """[(path, clade, depth)] for every model on the lineage path, deepest last; pan models depth 0."""
    out = [(PAN / f"trn{f}.cm", "pan", 0) for f in FAMS]
    if models_dir is None:
        return out
    lin = lineage.split(";")
    for name, m in meta.items():
        if m["clade"] in lin and m.get("status") in ("ok", "cached"):
            out.append((Path(models_dir) / f"{name}.cm", m["clade"], lin.index(m["clade"])))
    return out


def clen_of(cm):
    for line in open(cm):
        if line.startswith("CLEN"):
            return int(line.split()[1])
    return 0


def add_anticodons(rows, seq, tmp, emax=1.0):
    """Anticodon of every hit (E <= emax), read from its alignment to the pan CM of the hit's family."""
    by = {}
    for i, r in enumerate(rows):
        r["anticodon"] = ""
        if float(r["evalue"]) <= emax:
            by.setdefault(r["family"], []).append(i)
    for fam, idx in by.items():
        fa = tmp / f"ac_{fam}.fa"
        with open(fa, "w") as fh:
            for i in idx:
                r = rows[i]
                sub = seq[r["start"] - 1:r["end"]]
                fh.write(f">h{i}\n{sub if r['strand'] == 1 else rc(sub)}\n")
        sto = tmp / f"ac_{fam}.sto"
        sh([str(INF / "cmalign"), "--cpu", "1", "--outformat", "Stockholm", "-o", str(sto), str(PAN / f"trn{fam}.cm"), str(fa)])
        if sto.exists():
            al, ss, _ = parse_sto(sto.read_text())
            for k, v in al.items():
                rows[int(k[1:])]["anticodon"] = anticodon_from_alignment(v, ss)


RESCUE_E = 1e-3
RESCUE_MODE = "--mid"   # --max is 10x slower (MITOS2 sensitive pass); --mid drops the HMM filter strictness


def info_fam(p):
    return Path(p).stem.split("__")[-1].replace("trn", "")


def parse_tbl(tbl, info, acc, L, mode):
    rows = []
    for line in tbl.read_text().splitlines() if tbl.exists() else []:
        if line.startswith("#"):
            continue
        c = line.split()
        clade, depth, clen, fam = info.get(c[2], ("?", -1, 0, "?"))
        s, e, strand = int(c[7]), int(c[8]), c[9]
        lo, hi = min(s, e), max(s, e)
        if lo > L:                                   # duplicate in the appended part
            continue
        rows.append({"acc": acc, "model": c[2], "clade": clade, "depth": depth, "family": fam, "start": lo, "end": hi,
                     "strand": 1 if strand == "+" else -1, "mfrom": c[5], "mto": c[6], "clen": clen, "trunc": c[10],
                     "bits": c[14], "evalue": c[15], "gc": c[12], "mode": mode})
    return rows


def scan_one(job):
    acc, fa, lineage, models_dir, meta, odir, sens, rescue = job
    out = Path(odir) / f"{acc}.hits.tsv"
    if out.exists():
        return acc
    seq = "".join(Path(fa).read_text().split("\n")[1:]).strip()
    L = len(seq)
    tmp = Path(odir) / f"{acc}.tmp"
    tmp.mkdir(exist_ok=True)
    q = tmp / "q.fa"
    q.write_text(f">{acc}\n{seq + seq[:WRAP]}\n")
    rows = []
    paths = model_paths(models_dir, lineage, meta)
    cat = tmp / "all.cm"
    with open(cat, "w") as fh:
        for p, _, _ in paths:
            fh.write(Path(p).read_text())
    info = {}
    for p, clade, depth in paths:
        nm = next(line.split()[1] for line in open(p) if line.startswith("NAME"))
        info[nm if clade != "pan" else nm] = (clade, depth, clen_of(p), p.stem.split("__")[-1].replace("trn", ""))
    tbl = tmp / "hits.tbl"
    cmd = [str(INF / "cmsearch"), "--cpu", "1", "--tblout", str(tbl), "-E", "10", "--notextw"]
    if sens:
        cmd += ["--max"]
    sh(cmd + [str(cat), str(q)])
    rows = parse_tbl(tbl, info, acc, L, "")
    if rescue:                                        # second pass: --max for families without a good hit
        good = {r["family"] for r in rows if float(r["evalue"]) <= RESCUE_E}
        miss = [f for f in FAMS if f not in good]
        if miss:
            cat2 = tmp / "miss.cm"
            with open(cat2, "w") as fh:
                for p, clade, depth in paths:
                    if info_fam(p) in miss:
                        fh.write(Path(p).read_text())
            sh([str(INF / "cmsearch"), "--cpu", "1", "--tblout", str(tbl), "-E", "10", "--notextw", RESCUE_MODE,
                str(cat2), str(q)])
            seen = {(r["model"], r["start"], r["end"]) for r in rows}
            rows += [r for r in parse_tbl(tbl, info, acc, L, "max") if (r["model"], r["start"], r["end"]) not in seen]
    add_anticodons(rows, seq + seq[:WRAP], tmp)
    write_tsv(out, rows, FIELDS)
    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()
    return acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("models", help="models dir with models.json, or 'pan' for MITOS2 CMs only")
    ap.add_argument("out")
    ap.add_argument("--split", default="test")
    ap.add_argument("--accs", default="")
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--max", action="store_true")
    ap.add_argument("--rescue", action="store_true", help="--max second pass for families without a hit E <= 1e-3")
    a = ap.parse_args()
    ds, out = Path(a.dataset), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    md = None if a.models == "pan" else Path(a.models)
    meta = json.load(open(md / "models.json"))["models"] if md else {}
    want = set(a.accs.split(",")) if a.accs else None
    recs = [r for r in read_tsv(ds / "records.tsv") if (r["acc"] in want if want else r["split"] == a.split)]
    jobs = [(r["acc"], str(ds / "fasta" / f"{r['acc']}.fa"), r["lineage"], str(md) if md else None, meta, str(out), a.max,
             a.rescue)
            for r in recs]
    n = 0
    with ProcessPoolExecutor(a.jobs) as ex:
        for _ in ex.map(scan_one, jobs):
            n += 1
            if n % 100 == 0:
                print(n, flush=True)
    print("scanned", n)


if __name__ == "__main__":
    main()
