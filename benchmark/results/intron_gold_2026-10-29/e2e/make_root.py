"""End-to-end intron check against the gold standard (2026-10-29): the 193 gold records + 76 intron-free records +
KT989324.1 of intron-algo's harness (sets.tsv) as a mito_bench_run root: gb/<acc>.gb (NCBI record), bare/<acc>/<acc>.gbf
(MitoM input), set.tsv (mito_bench_run columns) and groups.tsv (dev / new / frag / neg / unann). Symlinks only."""
import csv
import sys
from pathlib import Path

sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
from mito_sample import COLUMNS, MITOZ_CLADE  # noqa: E402

R = Path(__file__).resolve().parent
H = Path("/tmp/claude-1000/-media-deyuan-18T-B-yzz/6bcce157-6435-4fde-a47f-55b33ff6135c/scratchpad/intron-gold/sets.tsv")
rows = list(csv.DictReader(open(H), delimiter="\t"))
for d in ("gb", "bare"):
    (R / d).mkdir(exist_ok=True)
out = []
for r in rows:
    acc = r["acc"]
    g = R / "gb" / f"{acc}.gb"
    if not g.exists():
        g.symlink_to(r["truth"])
    b = R / "bare" / acc
    b.mkdir(exist_ok=True)
    if not (b / f"{acc}.gbf").exists():
        (b / f"{acc}.gbf").symlink_to(r["bare"])
    out.append({**{c: "" for c in COLUMNS}, "accession": acc, "phylum": r["phylum"], "genetic_code": r["code"],
                "topology": r["topology"], "length": r["length"], "mitoz_clade": MITOZ_CLADE.get(r["phylum"]) or "",
                "split": "test", "retrieved": "local"})
with open(R / "set.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, COLUMNS, delimiter="\t")
    w.writeheader()
    w.writerows(out)
with open(R / "groups.tsv", "w") as fh:
    fh.write("accession\tgroup\tphylum\n" + "".join(f"{r['acc']}\t{r['group']}\t{r['phylum']}\n" for r in rows))
print(len(out), "records", file=sys.stderr)
