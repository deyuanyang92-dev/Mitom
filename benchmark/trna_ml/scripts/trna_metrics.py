#!/usr/bin/env python3
"""tRNA precision / recall / boundary accuracy per tool from mito_bench_eval.py's features.tsv, on the records every
listed tool had finished when the evaluation ran (status ok in its timing.tsv).

found     right amino acid (Leu / Ser class by anticodon where both sides have one), right strand, >= 50 % of the
          NCBI tRNA covered (mito_bench_eval.match); 'name unverifiable' (Leu / Ser without anticodon) counts as found
precision found / every predicted tRNA (found + wrong strand + wrong name + not matching any NCBI tRNA)
exact     both ends as NCBI; within 1 / within 3 = both ends within that many nt (share of the NCBI tRNAs)
usage: trna_metrics.py ROOT TOOL [TOOL ...] [--by phylum --set SET.tsv] [--accs FILE]
"""
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

FOUND = ("exact", "le3", "gt3", "name_unverifiable")


def score(rows):
    c = Counter()
    for r in rows:
        if r["truth"]:
            c["truth"] += 1
            if r["class"] in FOUND:
                c["found"] += 1
                d = int(r["delta"] or 0)
                c["exact"] += d == 0
                c["le1"] += d <= 1
                c["le3"] += d <= 3
            elif r["class"] == "strand_error":
                c["strand"] += 1
            elif r["class"] == "name_error":
                c["name"] += 1
        elif r["class"] == "fp":
            c["fp"] += 1
    c["pred"] = c["found"] + c["strand"] + c["name"] + c["fp"]
    return c


def line(name, c, nrec):
    t, p = max(c["truth"], 1), max(c["pred"], 1)
    rec, prec = c["found"] / t, c["found"] / p
    f1 = 2 * rec * prec / max(rec + prec, 1e-9)
    return (f"{name}\t{nrec}\t{c['truth']}\t{c['pred']}\t{100 * rec:.2f}%\t{100 * prec:.2f}%\t{100 * f1:.2f}%\t"
            f"{100 * c['exact'] / t:.2f}%\t{100 * c['le1'] / t:.2f}%\t{100 * c['le3'] / t:.2f}%\t"
            f"{c['fp'] / max(nrec, 1):.3f}\t{c['strand']}\t{c['name']}")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("root", type=Path)
    p.add_argument("tools", nargs="+")
    p.add_argument("--features", type=Path, help="default ROOT/eval/all/features.tsv")
    p.add_argument("--timing", type=Path, help="status table of the same evaluation (default: timing.tsv next to it)")
    p.add_argument("--set", type=Path)
    p.add_argument("--by", choices=["phylum"], default=None)
    p.add_argument("--accs", type=Path, help="only these records")
    a = p.parse_args(argv)
    feats = a.features or a.root / "eval" / "all" / "features.tsv"
    rows = [r for r in csv.DictReader(open(feats), delimiter="\t") if r["kind"] == "tRNA" and r["tool"] in a.tools]
    accs = {r["accession"] for r in rows}
    if a.accs:
        keep = {x.strip() for x in a.accs.read_text().split() if x.strip()}
        accs = {x for x in accs if x in keep or x.split(".")[0] in keep}
    ok = {(r["tool"], r["accession"]) for r in csv.DictReader(open(a.timing or feats.parent / "timing.tsv"), delimiter="\t")
          if r["status"] == "ok"}
    common = {x for x in accs if all((t, x) in ok for t in a.tools)}
    group = {}
    if a.by:
        group = {r["accession"]: r[a.by] for r in csv.DictReader(open(a.set), delimiter="\t")}
    print(f"records every tool finished: {len(common)}")
    print("tool\trecords\tNCBI tRNAs\tpredicted\trecall\tprecision\tF1\tboth ends exact\twithin 1\twithin 3\t"
          "FP per record\twrong strand\twrong name")
    for t in a.tools:
        sel = [r for r in rows if r["tool"] == t and r["accession"] in common]
        print(line(t, score(sel), len(common)))
    if a.by:
        by = defaultdict(list)
        for r in rows:
            if r["accession"] in common:
                by[(group.get(r["accession"], "?"), r["tool"])].append(r)
        nrec = Counter(group.get(x, "?") for x in common)
        for g in sorted(nrec, key=lambda g: -nrec[g]):
            for t in a.tools:
                print(line(f"{g}:{t}", score(by[(g, t)]), nrec[g]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
