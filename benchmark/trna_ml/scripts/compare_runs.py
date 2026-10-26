#!/usr/bin/env python3
"""Do two runs give the same tRNAs? Per record (both runs finished): tRNAs as (name, start, end, strand) after
mito_bench_eval.norm_features; reports identical records and the tRNAs only one run has.

usage: compare_runs.py ROOT TOOL_A TOOL_B [--set SET.tsv]
"""
import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mito_bench_eval import norm_features  # noqa: E402
from Bio import SeqIO  # noqa: E402


def trnas(path: Path, code: int):
    rec = SeqIO.read(str(path), "genbank")
    n = len(rec)
    return {(f.name, f.start % n, f.end % n or n, f.strand) for f in norm_features(rec, code) if f.kind == "tRNA"}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("root", type=Path)
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--set", type=Path)
    a = p.parse_args(argv)
    code = {r["accession"]: int(r["genetic_code"] or 5) for r in csv.DictReader(open(a.set), delimiter="\t")} if a.set else {}
    same = diff = 0
    only_a, only_b = [], []
    for st in sorted((a.root / "runs" / a.a).glob("*/status.json")):
        acc = st.parent.name
        sb = a.root / "runs" / a.b / acc / "status.json"
        if not sb.exists() or json.loads(st.read_text()).get("status") != "ok" or json.loads(sb.read_text()).get("status") != "ok":
            continue
        ta = trnas(st.parent / "final.gbf", code.get(acc, 5))
        tb = trnas(sb.parent / "final.gbf", code.get(acc, 5))
        if ta == tb:
            same += 1
        else:
            diff += 1
            only_a += [(acc, *x) for x in sorted(ta - tb)]
            only_b += [(acc, *x) for x in sorted(tb - ta)]
    print(f"records with the same tRNAs: {same} of {same + diff}; tRNAs only in {a.a}: {len(only_a)}, only in {a.b}: {len(only_b)}")
    for x in only_a[:20]:
        print(f"  only {a.a}:", *x)
    for x in only_b[:20]:
        print(f"  only {a.b}:", *x)
    return 0


if __name__ == "__main__":
    sys.exit(main())
