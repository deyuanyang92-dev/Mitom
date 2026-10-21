"""Markdown summary of one benchmark evaluation (mito_bench_eval.py output) for several tools.

Tables: records (ok / failed / with an obvious error), per feature kind (recall, both ends exact, within 3 nt,
name / strand errors, missed, extra), obvious-error counts by type, median run time, and the records where
each tool has errors (so the clusters can be looked at).

Usage: python3 bench_report.py --eval <root>/eval/all --tools mitom,mitoz,mitos2,mitofinder [--set sets/x.tsv] [--out report.md]
"""
import argparse
import csv
import gzip
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


def read(path: Path):
    """Rows of a TSV; a gzip-compressed copy (features.tsv.gz in benchmark/results/) is read when the plain file is absent."""
    gz = path.with_name(path.name + ".gz")
    if not path.exists() and gz.exists():
        with gzip.open(gz, "rt", encoding="utf-8") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pct(a, b):
    return f"{100 * a / b:.1f}%" if b else "-"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--eval", required=True)
    p.add_argument("--tools", default="mitom,mitoz,mitos2,mitofinder")
    p.add_argument("--stratum", default="all", help="stratum name in summary.tsv (default all)")
    p.add_argument("--set", help="benchmark set TSV: adds a per-phylum table (protein genes, both ends exact)")
    p.add_argument("--out")
    a = p.parse_args(argv)
    E = Path(a.eval)
    tools = a.tools.split(",")
    summ = [r for r in read(E / "summary.tsv") if r["stratum"] == a.stratum and r["tool"] in tools]
    errors = [r for r in read(E / "errors.tsv") if r["tool"] in tools]
    timing = [r for r in read(E / "timing.tsv") if r["tool"] in tools]
    lines = [f"# Benchmark: {', '.join(tools)} ({E})", ""]

    lines += ["## Records", "", "| tool | records | ok | failed | timeout | unsupported | with obvious error | "
              "median s/record | median peak RAM (MB) | max peak RAM (MB) |", "|---|---|---|---|---|---|---|---|---|---|"]
    for t in tools:
        r = next((x for x in summ if x["tool"] == t and x["kind"] == "record"), None)
        secs = [float(x["seconds"]) for x in timing if x["tool"] == t and x["status"] == "ok"]
        mem = [float(x["peak_rss_mb"]) for x in timing if x["tool"] == t and x["status"] == "ok" and x.get("peak_rss_mb")]
        if r:
            lines.append(f"| {t} | {r['records']} | {r['ok']} | {r['failed']} | {r['timeout']} | {r['unsupported']} | "
                         f"{r['with_obvious_error']} | {f'{statistics.median(secs):.0f}' if secs else '-'} | "
                         f"{f'{statistics.median(mem):.0f}' if mem else '-'} | {f'{max(mem):.0f}' if mem else '-'} |")
    lines.append("")

    for kind in ("PCG", "rRNA", "tRNA"):
        lines += [f"## {kind} (records where every tool finished: stratum `common_ok` below; this table: {a.stratum})", "",
                  "| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |",
                  "|---|---|---|---|---|---|---|---|---|---|---|"]
        for t in tools:
            r = next((x for x in summ if x["tool"] == t and x["kind"] == kind), None)
            if not r:
                continue
            n = int(r["n_truth"])
            lines.append(f"| {t} | {n} | {r['found']} | {pct(int(r['found']), n)} | {r['exact']} ({pct(int(r['exact']), n)}) | "
                         f"{r['le3']} ({pct(int(r['exact']) + int(r['le3']), n)} cum.) | {r['gt3']} | {r['strand_error']} | "
                         f"{r['name_error']} | {r['missed']} | {r['fp']} |")
        lines.append("")
    common = [r for r in read(E / "summary.tsv") if r["stratum"] == "common_ok" and r["tool"] in tools and r["kind"] != "record"]
    if common:
        lines += ["## Same records for every tool (common_ok)", "", "| tool | kind | truth | recall | both ends exact | within 3 nt |",
                  "|---|---|---|---|---|---|"]
        for kind in ("PCG", "rRNA", "tRNA"):
            for t in tools:
                r = next((x for x in common if x["tool"] == t and x["kind"] == kind), None)
                if r:
                    n = int(r["n_truth"])
                    lines.append(f"| {t} | {kind} | {n} | {pct(int(r['found']), n)} | {pct(int(r['exact']), n)} | "
                                 f"{pct(int(r['exact']) + int(r['le3']), n)} |")
        lines.append("")

    if a.set:
        phylum = {r["accession"]: r["phylum"] for r in read(Path(a.set))}
        n_rec, n, ok = Counter(), Counter(), Counter()
        for r in read(E / "features.tsv"):
            ph = phylum.get(r["accession"])            # evaluation may hold records outside this set
            if ph and r["tool"] in tools and r["kind"] == "PCG" and r["truth"]:
                n[(ph, r["tool"])] += 1
                ok[(ph, r["tool"])] += r["class"] == "exact"
        for acc, ph in phylum.items():
            n_rec[ph] += 1
        lines += ["## Protein genes, both ends exact, per phylum", "",
                  "| phylum | records | " + " | ".join(tools) + " |", "|---|---|" + "---|" * len(tools)]
        for ph in sorted(n_rec, key=lambda x: (-n_rec[x], x)):
            lines.append(f"| {ph} | {n_rec[ph]} | " + " | ".join(
                pct(ok[(ph, t)], n[(ph, t)]) if n[(ph, t)] else "-" for t in tools) + " |")
        lines.append("")

    ok_by_acc = defaultdict(set)
    for r in timing:
        if r["status"] == "ok":
            ok_by_acc[r["accession"]].add(r["tool"])
    common_acc = {a for a, ts in ok_by_acc.items() if set(tools) <= ts}
    for title, rows in (("Obvious errors by type", errors),
                        (f"Obvious errors by type, records every tool finished ({len(common_acc)} records)",
                         [r for r in errors if r["accession"] in common_acc])):
        lines += [f"## {title}", "", "| error | " + " | ".join(tools) + " |", "|---|" + "---|" * len(tools)]
        by = defaultdict(Counter)
        for r in rows:
            by[r["error"]][r["tool"]] += 1
        for err in sorted(by, key=lambda e: -sum(by[e].values())):
            lines.append(f"| {err} | " + " | ".join(str(by[err][t]) for t in tools) + " |")
        recs = Counter(t for t, a in {(r["tool"], r["accession"]) for r in rows})
        lines.append("| records with any error | " + " | ".join(str(recs[t]) for t in tools) + " |")
        lines.append("")

    lines += ["## Error records per tool", ""]
    for t in tools:
        recs = defaultdict(list)
        for r in errors:
            if r["tool"] == t:
                recs[r["accession"]].append(f"{r['error']}:{r['gene']}" if r["gene"] else r["error"])
        lines.append(f"**{t}** ({len(recs)} records)")
        lines.append("")
        for acc in sorted(recs):
            lines.append(f"- {acc}: {', '.join(recs[acc])}")
        lines.append("")
    text = "\n".join(lines)
    if a.out:
        Path(a.out).write_text(text, encoding="utf-8")
        print(f"[report] {a.out}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
