"""Score the runs of run_all.sh with MitoM's benchmark/eval_introns.py against the NCBI records, and compare each
rotated run with the unrotated one of the same code (new MitoM).
    score.py WORKDIR MITOM_BENCHMARK_DIR > score.txt"""
import collections
import csv
import json
import statistics
import subprocess
import sys
from pathlib import Path

W, B = Path(sys.argv[1]), Path(sys.argv[2])
groups = {r["accession"]: r for r in csv.DictReader(open(B / "results/intron_gold_2026-10-29/e2e/groups.tsv"),
                                                     delimiter="\t")}
SETS = [s for s in ["orig", "exon", "old_exon", "intron", "old_intron", "near3", "old_near3", "near5", "old_near5"]
        if (W / s / "runs" / "mitom_lo").is_dir()]
rows = {}
for s in SETS:
    runs = W / s / "runs" / "mitom_lo"
    st = {p.parent.name: json.loads(p.read_text()) for p in runs.glob("*/status.json")}
    accs = sorted(a for a, x in st.items() if x["status"] == "ok")
    bad = {a: x["error"][:80] for a, x in st.items() if x["status"] != "ok"}
    out = W / "eval"
    out.mkdir(exist_ok=True)
    subprocess.run([sys.executable, str(B / "eval_introns.py"), "--gb", str(W / s / "gb"), "--pred", str(runs),
                    "--label", s, "--acc", ",".join(accs), "-o", str(out)], check=True, capture_output=True)
    rows[s] = list(csv.DictReader(open(out / f"{s}_introns.tsv"), delimiter="\t"))
    secs = [x["seconds"] for x in st.values() if x["status"] == "ok"]
    print(f"{s}: {len(st)} runs, ok {len(accs)}, not ok {len(bad)} {bad if bad else ''}; "
          f"median {statistics.median(secs):.0f} s, max {max(secs):.0f} s per record")


def table(sets, accs, title):
    print(f"\n## {title} ({len(accs)} records)")
    print("set\tgroup\tNCBI_introns\tfound\texact\twithin6nt\tmore6nt\tmissed\tfalse_introns")
    for s in sets:
        by = collections.defaultdict(collections.Counter)
        for r in rows[s]:
            if r["accession"] not in accs:
                continue
            g = groups[r["accession"]]
            for k in ("all", f"{g['group']}", f"{g['group']}:{g['phylum']}"):
                by[k][r["class"]] += 1
        for k in sorted(by, key=lambda k: (k != "all", k)):
            c = by[k]
            f = c["exact"] + c["le6"] + c["gt6"]
            print(f"{s}\t{k}\t{f + c['missed']}\t{f}\t{c['exact']}\t{c['le6']}\t{c['gt6']}\t{c['missed']}\t{c['fp']}")


def done(s):
    return {r["accession"] for r in rows[s]}


for rot in [s for s in ("exon", "intron", "near3", "near5") if s in SETS]:
    accs = done("orig") & done(rot) & done(f"old_{rot}")
    table(["orig", rot, f"old_{rot}"], accs, f"origin inside the gene: {rot}")
    # consistency: the same NCBI intron found / placed identically with and without rotation (new code)
    key = lambda r: (r["accession"], r["truth"] or r["pred"], r["class"] == "fp")
    a = {key(r): r for r in rows["orig"] if r["accession"] in accs}
    b = {key(r): r for r in rows[rot] if r["accession"] in accs}
    same = sum(1 for k in a if k in b and a[k]["pred"] == b[k]["pred"] and a[k]["class"] == b[k]["class"])
    diff = [(k, a[k]["class"], a[k]["pred"], b.get(k, {}).get("class"), b.get(k, {}).get("pred"))
            for k in a if not (k in b and a[k]["pred"] == b[k]["pred"] and a[k]["class"] == b[k]["class"])]
    diff += [(k, None, None, b[k]["class"], b[k]["pred"]) for k in b if k not in a]
    print(f"\nnew MitoM, unrotated vs {rot}: {same} identical rows, {len(diff)} differ")
    for d in diff:
        print("  ", *d)

print("\n## rows not found / false in the rotated runs of the new MitoM")
for s in [s for s in ("exon", "intron", "near3", "near5") if s in SETS]:
    for r in rows[s]:
        if r["class"] in ("missed", "gt6", "fp"):
            print(s, r["accession"], groups[r["accession"]]["group"], groups[r["accession"]]["phylum"], r["organism"],
                  r["class"], r["gene"] or r["pred_gene"], r["truth"], r["pred"], r["report"][:120], sep="\t")

print("\n## topology verdicts in the runs (reports/topology_decision.tsv)")
for s in SETS:
    c = collections.Counter()
    for f in (W / s / "runs" / "mitom_lo").glob("*/raw/**/topology_decision.tsv"):
        for r in csv.DictReader(open(f), delimiter="\t"):
            c[r["topology"] + (" via intron" if "intron" in r["reason"] else "")] += 1
    print(s, dict(c))
