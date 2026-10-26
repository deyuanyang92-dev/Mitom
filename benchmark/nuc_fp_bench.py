#!/usr/bin/env python3
"""nuc_fp_bench.py - precision of mitom nuc: features it reports that the GenBank record does not annotate.

nuc_bench.py scores the boundaries of the longest feature per record and gene, so duplicate copies and false
features are never counted there. Here every reported 18S / 5.8S / 28S / histone feature is classified:
  annotated   overlaps a feature of the same gene in the record (rRNA, misc_RNA "contains 18S ...", gene or CDS,
              partial ones included)
  duplicate   overlaps an earlier reported feature of the same gene on the same strand
  unannotated none of these: a false feature, or a real one the submitter did not annotate (check by hand; in
              negative sets such as complete mitogenomes every one is false)
Missed annotated genes are counted too (recall), for complete features only.
References of the record's own genus are left out (mito_nuc annotate --leave_out).

Usage:
  python3 benchmark/nuc_fp_bench.py --gb rdna_heldout.gb histone_heldout.gb mito.gb -o nuc_fp [--lib models/nuc]
Writes in -o: <set>/ (mitom nuc output), <set>.features.tsv (each reported feature with its class), report.txt.
"""
import argparse
import csv
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_nuc as nuc  # noqa: E402

GENES = nuc.RRNA + nuc.HISTONES


def annotated(rec):
    """gene -> [(lo, hi, complete)] of the record's own features, any type that names the gene"""
    out = defaultdict(list)
    for f in rec.features:
        if f.type not in ("rRNA", "misc_RNA", "gene", "CDS", "misc_feature"):
            continue
        lab = " ".join(f.qualifiers.get(k, [""])[0] for k in ("gene", "product", "note", "label"))
        names = {nuc.rrna_name(f), nuc.histone_name(f)}
        for g, pat in (("18S", r"\b18S\b"), ("5.8S", r"\b5\.8S\b"), ("28S", r"\b28S\b")):
            if re.search(pat, lab, re.I):
                names.add(g)
        for g in names - {None}:
            out[g].append((int(f.location.start) + 1, int(f.location.end), not any(nuc._partial(f))))
    return out


def classify(rows, recs):
    seen = defaultdict(list)
    for r in sorted(rows, key=lambda r: (r["contig"], int(r["start"]))):
        g, a, b = r["gene"], int(r["start"]), int(r["end"])
        if g not in GENES:
            continue
        ann = annotated(recs[r["contig"]]).get(g, [])
        if any(a <= y and x <= b for x, y, _ in ann):
            r["class"] = "annotated"
        elif any(a <= y and x <= b for x, y in seen[(r["contig"], g, r["strand"])]):
            r["class"] = "duplicate"
        else:
            r["class"] = "unannotated"
        if r["class"] != "duplicate" and any(a <= y and x <= b for x, y in seen[(r["contig"], g, r["strand"])]):
            r["class"] += "+duplicate"
        seen[(r["contig"], g, r["strand"])].append((a, b))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gb", nargs="+", type=Path, required=True)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--lib", type=Path, default=nuc.LIB)
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    report = []
    for gb in a.gb:
        name = gb.stem
        recs = {re.sub(r"[^A-Za-z0-9_.-]", "_", r.id): r for r in SeqIO.parse(str(gb), "genbank")}
        lo = a.out / f"{name}.leave_out.tsv"
        lo.write_text("".join(f"{k}\t{r.annotations.get('organism', '').split(' ')[0]}\n" for k, r in recs.items()))
        od = a.out / name
        if not (od / "nuc_features.tsv").exists():
            subprocess.run([sys.executable, str(HERE.parent / "mito_nuc.py"), "annotate", "-i", str(gb), "-o", str(od),
                            "--lib", str(a.lib), "--threads", str(a.threads), "--leave_out", str(lo)], check=True)
        rows = classify(list(csv.DictReader(open(od / "nuc_features.tsv"), delimiter="\t")), recs)
        with open(a.out / f"{name}.features.tsv", "w") as o:
            w = csv.DictWriter(o, list(rows[0]) if rows else ["contig"], delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        c = Counter((r["gene"], r["class"].split("+")[0]) for r in rows if r.get("class"))
        dup = Counter(r["gene"] for r in rows if "duplicate" in r.get("class", ""))
        miss = Counter()
        tot = Counter()
        got = defaultdict(list)
        for r in rows:
            got[(r["contig"], r["gene"])].append((int(r["start"]), int(r["end"])))
        for k, rec in recs.items():
            for g, ann in annotated(rec).items():
                for x, y, full in ann:
                    if g in GENES and full:
                        tot[g] += 1
                        miss[g] += not any(p <= y and x <= q for p, q in got[(k, g)])
        report.append(f"## {name}: {len(recs)} records, {sum(c.values())} reported features")
        report.append("gene\treported\tannotated\tunannotated\tduplicate\tcomplete annotated\tmissed")
        for g in GENES:
            n = sum(v for (gg, _), v in c.items() if gg == g)
            if n or tot[g]:
                report.append(f"{g}\t{n}\t{c[(g, 'annotated')]}\t{c[(g, 'unannotated')]}\t{dup[g]}\t{tot[g]}\t{miss[g]}")
        report.append("")
    (a.out / "report.txt").write_text("\n".join(report))
    print("\n".join(report))


if __name__ == "__main__":
    main()
