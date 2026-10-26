#!/usr/bin/env python3
"""Take MiTFi's ends for the tRNAs MiTFi also calls (same name and strand, >= 50 % of the shorter one), for the records
whose lineage has a given token (Vertebrata: the learned ends never saw a vertebrate), from runs/<from>/ (MitoM-ML calls)
and runs/<mitfi>/ (MiTFi calls); writes runs/<name>/<acc>/final.gbf + status.json.

usage: mitfi_ends.py --set SET.tsv --root ROOT --from res_mitfi_0p1 --mitfi trna_mitfi --name ml_vertends [--lineage Vertebrata]
"""
import argparse
import csv
import json
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation


def iv(f, n):
    parts = f.location.parts
    if len(parts) == 2 and any(int(p.start) == 0 for p in parts) and any(int(p.end) == n for p in parts):
        t = next(p for p in parts if int(p.end) == n)
        h = next(p for p in parts if int(p.start) == 0)
        return int(t.start), n + int(h.end)
    return int(f.location.start), int(f.location.end)


def ov(a, b, n):
    return max(min(a[1], b[1] + k) - max(a[0], b[0] + k) for k in (-n, 0, n))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", required=True)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--from", dest="src", required=True)
    p.add_argument("--mitfi", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--lineage", default="Vertebrata")
    a = p.parse_args(argv)
    done = changed = 0
    for r in csv.DictReader(open(a.set), delimiter="\t"):
        acc = r["accession"]
        fp, mp = a.root / "runs" / a.src / acc / "final.gbf", a.root / "runs" / a.mitfi / acc / "final.gbf"
        if not fp.exists():
            continue
        rec = SeqIO.read(str(fp), "genbank")
        lineage = SeqIO.read(str(a.root / "gb" / f"{acc}.gb"), "genbank").annotations.get("taxonomy", [])
        n = len(rec)
        if mp.exists() and a.lineage in lineage:
            mf = [(f.qualifiers.get("gene", [""])[0], f.location.strand, iv(f, n)) for f in SeqIO.read(str(mp), "genbank").features
                  if f.type == "tRNA"]
            for f in rec.features:
                if f.type not in ("tRNA", "gene"):
                    continue
                name, st, x = f.qualifiers.get("gene", [""])[0], f.location.strand, iv(f, n)
                m = [y for nm, s2, y in mf if nm == name and s2 == st and ov(x, y, n) >= 0.5 * min(x[1] - x[0], y[1] - y[0])]
                if m and m[0] != x:
                    s, e = m[0]
                    f.location = FeatureLocation(s, e, st) if e <= n else CompoundLocation(
                        [FeatureLocation(s, n, st), FeatureLocation(0, e - n, st)][::1 if st == 1 else -1])
                    changed += f.type == "tRNA"
        out = a.root / "runs" / a.name / acc
        out.mkdir(parents=True, exist_ok=True)
        SeqIO.write(rec, str(out / "final.gbf"), "genbank")
        (out / "status.json").write_text(json.dumps({"tool": a.name, "accession": acc, "status": "ok", "seconds": 0}))
        done += 1
    print(f"{a.name}: {done} records, {changed} tRNA ends taken from MiTFi", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
