"""Homolog table for the intron verification: per record (organism, lineage) and per standard gene the protein
(CDS /translation, else translated) or rRNA sequence, and whether that gene carries an intron there (an intron
feature inside it or a same-strand join gap >= 142 nt).

Sources: refpool_metazoa/Metazoa.gb (15,803 RefSeq records) and the 193 NCBI records of the gold set.
usage: homologs.py <out.pkl> [--jobs 16]
"""
import argparse
import os
import pickle
import sys
from io import StringIO
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
from mito_names import gene_name  # noqa: E402
from Bio import SeqIO  # noqa: E402

POOL = Path("/media/deyuan/18T_B/yzz/mitom_work/refpool_metazoa/Metazoa.gb")
GOLD = Path(__file__).resolve().parents[2]
MIN_INTRON = 142


def gaps(loc):
    """Same-strand gaps between consecutive parts (biological order); None for a part pair on two strands."""
    out = []
    for p, q in zip(loc.parts, loc.parts[1:]):
        if p.strand != q.strand:
            out.append(None)
            continue
        plus = (p.strand or 1) == 1
        s, e = (int(p.end), int(q.start)) if plus else (int(q.end), int(p.start))
        out.append(e - s)
    return out


def summarise(rec, source: str) -> dict:
    intr = [(int(f.location.start), int(f.location.end)) for f in rec.features if f.type == "intron"]
    out = {"acc": rec.id, "organism": rec.annotations.get("organism", ""), "lineage": rec.annotations.get("taxonomy", []),
           "source": source, "prot": {}, "rrna": {}}
    for f in rec.features:
        if f.type not in ("CDS", "rRNA"):
            continue
        g = gene_name(f.qualifiers, f.type)
        if not g:
            continue
        s, e = int(f.location.start), int(f.location.end)
        has_intron = any(x is not None and x >= MIN_INTRON for x in gaps(f.location)) or \
            any(s <= a and b <= e for a, b in intr)
        if f.type == "CDS":
            p = f.qualifiers.get("translation", [""])[0]
            if not p:
                try:
                    table = int(f.qualifiers.get("transl_table", ["1"])[0])
                    nt = f.location.extract(rec.seq)
                    p = str(nt[:len(nt) - len(nt) % 3].translate(table=table)).rstrip("*")
                except Exception:
                    p = ""
            if p and "*" not in p and len(p) >= 30 and g not in out["prot"]:
                out["prot"][g] = (p, has_intron)
        elif g not in out["rrna"]:
            try:
                out["rrna"][g] = (str(f.location.extract(rec.seq)).upper(), has_intron)
            except Exception:                  # record without sequence (CONTIG)
                pass
    return out


def chunk(args):
    path, start, end = args
    with open(path, "rb") as fh:
        fh.seek(start)
        text = fh.read(end - start).decode("utf-8", errors="ignore")
    return [summarise(r, "refpool") for r in SeqIO.parse(StringIO(text), "genbank")]


def offsets(path: Path, n: int):
    """n byte ranges that start at a LOCUS line."""
    size = path.stat().st_size
    cuts = [0]
    with open(path, "rb") as fh:
        for k in range(1, n):
            fh.seek(size * k // n)
            fh.readline()
            while True:
                pos = fh.tell()
                line = fh.readline()
                if not line or line.startswith(b"LOCUS"):
                    cuts.append(pos)
                    break
    cuts.append(size)
    return [(str(path), a, b) for a, b in zip(cuts, cuts[1:]) if b > a]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("out")
    p.add_argument("--jobs", type=int, default=16)
    a = p.parse_args()
    with Pool(a.jobs) as pool:
        data = [x for part in pool.map(chunk, offsets(POOL, a.jobs * 4)) for x in part]
    gold = [summarise(SeqIO.read(str(f), "genbank"), "gold_set") for f in sorted((GOLD / "ncbi_gb").glob("*.gb"))]
    seen = {x["acc"] for x in gold}
    data = gold + [x for x in data if x["acc"] not in seen]
    pickle.dump(data, open(a.out, "wb"))
    print(f"{len(data)} records ({len(gold)} gold-set records) -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
