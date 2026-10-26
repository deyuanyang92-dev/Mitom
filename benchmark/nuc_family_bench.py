#!/usr/bin/env python3
"""nuc_family_bench.py - reproduce the Polynoidae / Nereididae rDNA comparison (mitom nuc vs the NCBI records).

The record list is fixed (accession.version, benchmark/sets/nuc_polynoidae_nereididae_2026-10-25.tsv), so a
re-run downloads the same sequences. References of each record's own genus are left out. For every record the
NCBI rRNA / ITS features are compared with mitom's: d=<5' shift>/<3' shift> in nt (mitom minus NCBI), MISSED,
mitom-only parts; records annotated only as one "contains ..." misc_RNA are compared by span. The 8 nt on each
side of mitom's 18S 3' end and 28S 3' end are listed so the conserved ends (...GATCATTA, ...GATTTGT) can be checked.

Usage:
  python3 benchmark/nuc_family_bench.py -o nuc_family [--email you@example.org] [--threads 16]
Writes in -o: records.gb, leave_out.tsv, mitom/ (mitom nuc output, incl. submission/), comparison.tsv.
"""
import argparse
import csv
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_fetch  # noqa: E402

SET = HERE / "sets" / "nuc_polynoidae_nereididae_2026-10-25.tsv"
PAT = (("18S", r"18S|small subunit"), ("5.8S", r"5\.8S"), ("28S", r"28S|large subunit"),
       ("ITS1", r"internal transcribed spacer 1|ITS1"), ("ITS2", r"internal transcribed spacer 2|ITS2"))


def ncbi_parts(rec):
    out = []
    for f in rec.features:
        if f.type not in ("rRNA", "misc_RNA"):
            continue
        lab = " ".join((f.qualifiers.get(k) or [""])[0] for k in ("product", "note", "gene"))
        g = "unit" if f.type == "misc_RNA" and lab.startswith("contains") else \
            next((g for g, p in PAT if re.search(p, lab)), None)
        if g:
            out.append((g, int(f.location.start) + 1, int(f.location.end), "<" in str(f.location), ">" in str(f.location)))
    return out


def compare(recs, pred):
    rows = []
    for rec in recs:
        P = {r["gene"]: r for r in pred.get(rec.id, [])}
        s = rec.seq
        cmp = []
        for g, a, b, p5, p3 in ncbi_parts(rec):
            tag = f"{'<' if p5 else ''}{a}-{b}{'>' if p3 else ''}"
            if g == "unit":
                lo = min((int(r["start"]) for r in P.values()), default=None)
                hi = max((int(r["end"]) for r in P.values()), default=None)
                cmp.append(f"unit {tag} vs mitom {lo}-{hi}")
            elif g not in P:
                cmp.append(f"{g} {tag} MISSED")
            else:
                d5, d3 = int(P[g]["start"]) - a, int(P[g]["end"]) - b
                cmp.append(f"{g} {tag} d={d5:+d}/{d3:+d}")
        ends = {g: str(s[int(P[g]["end"]) - 8:int(P[g]["end"]) + 8]).lower() for g in ("18S", "28S") if g in P}
        rows.append({"accession": rec.id, "organism": rec.annotations.get("organism", ""), "length": len(rec),
                     "ncbi_vs_mitom": "; ".join(cmp),
                     "mitom_only": ",".join(f"{g}:{r['start']}-{r['end']}" for g, r in P.items()
                                            if g not in {x[0] for x in ncbi_parts(rec)}),
                     "end_18S": ends.get("18S", ""), "end_28S": ends.get("28S", "")})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--email")
    ap.add_argument("--api_key")
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    gb = a.out / "records.gb"
    accs = [r["accession"] for r in csv.DictReader(open(SET), delimiter="\t") if r["has_sequence"] == "1"]
    if not gb.exists():
        mito_fetch.setup_entrez(a.email, a.api_key)
        mito_fetch.download(accs, gb)
    recs = list(SeqIO.parse(str(gb), "genbank"))
    lo = a.out / "leave_out.tsv"
    lo.write_text("".join(f"{r.id}\t{r.annotations.get('organism', '').split(' ')[0]}\n" for r in recs))
    od = a.out / "mitom"
    subprocess.run([sys.executable, str(HERE.parent / "mito_nuc.py"), "annotate", "-i", str(gb), "-o", str(od),
                    "--genes", "rdna", "--threads", str(a.threads), "--leave_out", str(lo)], check=True)
    pred = defaultdict(list)
    for r in csv.DictReader(open(od / "nuc_features.tsv"), delimiter="\t"):
        pred[r["contig"]].append(r)
    rows = compare(recs, pred)
    with open(a.out / "comparison.tsv", "w", newline="") as o:
        w = csv.DictWriter(o, list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with_rdna = [r for r in rows if r["ncbi_vs_mitom"]]
    found = [r for r in with_rdna if r["accession"] in pred]
    print(f"{len(recs)} records, {len(with_rdna)} with NCBI rDNA features, {len(found)} annotated by mitom; "
          f"{sum('MISSED' in r['ncbi_vs_mitom'] for r in rows)} with a missed part  ({a.out / 'comparison.tsv'})")
    print((od / "submission" / "validation.txt").read_text())


if __name__ == "__main__":
    main()
