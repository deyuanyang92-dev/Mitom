#!/usr/bin/env python3
"""frozen_stats.py - the pre-registered metrics of the frozen evaluation (PLAN.md section 5) with confidence intervals.

Reads <eval>/features.tsv and <eval>/errors.tsv of mito_bench_eval.py and the set TSV. Only records every tool
finished (status ok) are compared. Per tool and metric: the value over all genes, and a 95 % interval from a
bootstrap that resamples families (records of one family are not independent); per pair "<tool> - <ref>": the paired
difference with its interval (same resampled families for both tools).

Metrics (fixed before the run; do not add one after seeing the results):
  pcg_exact     protein genes with both ends identical to NCBI / NCBI protein genes      (primary)
  pcg_le3       ... both ends within 3 nt
  pcg_found     protein genes found (exact, le3, gt3, name_unverifiable)
  rrna_le3      rRNAs with both ends within 3 nt
  trna_recall   tRNAs found / NCBI tRNAs
  trna_precision tRNAs found / tRNAs predicted (found + fp)
  err_free      records without an obvious error (mito_bench_eval errors.tsv, incl. table2asn ERROR)

Usage: python3 frozen_stats.py --eval <root>/eval/test --set frozen_set.tsv --tools mitom_native,mitos2,mitos2_cur
                               --ref mitos2 [--boot 2000 --seed 1] -o stats.tsv
"""
import argparse
import csv
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List

FOUND = {"exact", "le3", "gt3", "name_unverifiable"}
METRICS = ("pcg_exact", "pcg_le3", "pcg_found", "rrna_le3", "trna_recall", "trna_precision", "err_free")


def read(path: Path) -> List[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def counts(feats: List[dict], errors: List[dict], tools: List[str]) -> Dict[str, Dict[str, Counter]]:
    """{tool: {accession: Counter(numerators and denominators)}}"""
    c: Dict[str, Dict[str, Counter]] = {t: defaultdict(Counter) for t in tools}
    for f in feats:
        if f["tool"] not in c:
            continue
        x, k, cl = c[f["tool"]][f["accession"]], f["kind"], f["class"]
        if k == "PCG" and cl != "fp":
            x["pcg_n"] += 1
            x["pcg_exact"] += cl == "exact"
            x["pcg_le3"] += cl in ("exact", "le3")
            x["pcg_found"] += cl in FOUND
        elif k == "rRNA" and cl != "fp":
            x["rrna_n"] += 1
            x["rrna_le3"] += cl in ("exact", "le3")
        elif k == "tRNA":
            if cl == "fp":
                x["trna_fp"] += 1
            else:
                x["trna_n"] += 1
                x["trna_recall"] += cl in FOUND
    bad = {(e["tool"], e["accession"]) for e in errors}
    for t in tools:
        for acc, x in c[t].items():
            x["rec_n"] = 1
            x["err_free"] = int((t, acc) not in bad)
    return c


def value(xs: List[Counter], m: str) -> float:
    s = Counter()
    for x in xs:
        s.update(x)
    num, den = {"pcg_exact": ("pcg_exact", "pcg_n"), "pcg_le3": ("pcg_le3", "pcg_n"), "pcg_found": ("pcg_found", "pcg_n"),
                "rrna_le3": ("rrna_le3", "rrna_n"), "trna_recall": ("trna_recall", "trna_n"),
                "err_free": ("err_free", "rec_n")}.get(m, (None, None))
    if m == "trna_precision":
        d = s["trna_recall"] + s["trna_fp"]
        return s["trna_recall"] / d if d else float("nan")
    return s[num] / s[den] if s[den] else float("nan")


def ci(vals: List[float]) -> str:
    v = sorted(x for x in vals if x == x)
    if not v:
        return ""
    return f"{100 * v[int(0.025 * (len(v) - 1))]:.1f}..{100 * v[int(0.975 * (len(v) - 1))]:.1f}"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--eval", required=True, type=Path)
    p.add_argument("--set", required=True, type=Path)
    p.add_argument("--tools", required=True)
    p.add_argument("--ref", default="mitos2")
    p.add_argument("--boot", type=int, default=2000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("-o", "--out", type=Path)
    a = p.parse_args(argv)
    tools = a.tools.split(",")
    fam = {r["accession"]: r.get("family") or r["accession"] for r in read(a.set)}
    c = counts(read(a.eval / "features.tsv"), read(a.eval / "errors.tsv") if (a.eval / "errors.tsv").stat().st_size else [],
               tools)
    common = sorted(set(fam).intersection(*[set(c[t]) for t in tools]))
    groups: Dict[str, List[str]] = defaultdict(list)
    for acc in common:
        groups[fam[acc]].append(acc)
    names = sorted(groups)
    rng = random.Random(a.seed)
    boots = [[acc for g in (rng.choice(names) for _ in names) for acc in groups[g]] for _ in range(a.boot)]
    rows = []
    for m in METRICS:
        for t in tools:
            v = value([c[t][x] for x in common], m)
            bt = [value([c[t][x] for x in b], m) for b in boots]
            rows.append({"metric": m, "tool": t, "records": len(common), "families": len(names),
                         "value_pct": f"{100 * v:.1f}", "ci95_pct": ci(bt)})
            if t != a.ref and a.ref in tools:
                d = v - value([c[a.ref][x] for x in common], m)
                bd = [value([c[t][x] for x in b], m) - value([c[a.ref][x] for x in b], m) for b in boots]
                rows.append({"metric": m, "tool": f"{t} - {a.ref}", "records": len(common), "families": len(names),
                             "value_pct": f"{100 * d:+.1f}", "ci95_pct": ci(bd)})
    out = open(a.out, "w", newline="") if a.out else sys.stdout
    w = csv.DictWriter(out, list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
    if a.out:
        out.close()
        print(f"[frozen_stats] {len(common)} records of {len(names)} families, {len(tools)} tools -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
