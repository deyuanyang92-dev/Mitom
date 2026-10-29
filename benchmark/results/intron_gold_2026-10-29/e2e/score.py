"""Score the end-to-end runs (run.sh) with benchmark/eval_introns.py: gold records against gold/gb, the intron-free
records and KT989324.1 against their NCBI record; per tool, split by group (dev / new / frag / neg / unann), phylum and
gold class; plus run status and NCBI (table2asn) ERRORs."""
import collections
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

R = Path(__file__).resolve().parent
GOLD = Path("/media/deyuan/18T_B/yzz/mitom_work/intron_gold_2026-10-29/gold/gb")
EV = "/media/deyuan/18T_B/yzz/MitoM/benchmark/eval_introns.py"
groups = {r["accession"]: r for r in csv.DictReader(open(R / "groups.tsv"), delimiter="\t")}
truth = R / "truth"
truth.mkdir(exist_ok=True)
for acc, r in groups.items():
    t = truth / f"{acc}.gb"
    if t.is_symlink() or t.exists():
        t.unlink()
    t.symlink_to(GOLD / f"{acc}.gb" if r["group"] in ("dev", "new", "frag") else (R / "gb" / f"{acc}.gb").resolve())
cols = ["exact", "le6", "gt6", "missed", "fp"]
for tool in ("mitom_lo", "mitom_native"):
    runs = R / "runs" / tool
    if not runs.is_dir():
        continue
    st = {p.parent.name: json.loads(p.read_text()) for p in runs.glob("*/status.json")}
    accs = sorted(a for a in st if a in groups)
    status = collections.Counter(s["status"] for s in st.values())
    errs = collections.Counter()
    for a in accs:
        log = runs / a / "raw" / "run.log"
        m = re.search(r"NCBI ERROR in (\d+)", log.read_text(errors="ignore")) if log.exists() else None
        errs["records with NCBI ERROR" if m and int(m.group(1)) else "0 NCBI ERROR"] += 1
    out = R / "eval"
    out.mkdir(exist_ok=True)
    subprocess.run([sys.executable, EV, "--gb", str(truth), "--pred", str(runs), "--label", tool, "--acc", ",".join(accs),
                    "-o", str(out)], check=True, capture_output=True)
    tab = collections.defaultdict(collections.Counter)
    for x in csv.DictReader(open(out / f"{tool}_introns.tsv"), delimiter="\t"):
        g = groups[x["accession"]]
        keys = [g["group"], f"{g['group']}:{g['phylum']}", "all gold" if g["group"] in ("dev", "new", "frag") else None,
                f"class {x.get('truth_class') or '-'}" if g["group"] in ("dev", "new", "frag") else None]
        for k in keys:
            if k:
                tab[k][x["class"]] += 1
    print(f"== {tool}: {len(accs)} records, status {dict(status)}, {dict(errs)}")
    print("group\tn_truth\t" + "\t".join(cols))
    for k in sorted(tab, key=lambda k: (not k.startswith("all"), k.startswith("class"), k)):
        c = tab[k]
        print(f"{k}\t{c['exact'] + c['le6'] + c['gt6'] + c['missed']}\t" + "\t".join(str(c[x]) for x in cols))
    print()
