#!/usr/bin/env python3
"""mito_bench_run.py root for the RNA-seq pairs (datasets.tsv): <work>/tools_root/gb, bare, set.tsv - to annotate the
same mitogenomes with MITOS2 (and others) for the comparison in compare_tools.py.
Usage: python3 make_tool_set.py <work>
"""
import csv
import sys
from pathlib import Path

from Bio import SeqIO

H = Path(__file__).resolve()
sys.path.insert(0, str(H.parents[3]))
sys.path.insert(0, str(H.parents[2]))
from mito_sample import COLUMNS, MITOZ_CLADE, write_bare  # noqa: E402

W = Path(sys.argv[1])
pairs = list(csv.DictReader(open(H.parent / "datasets.tsv"), delimiter="\t"))
recs = {r.id: r for f in (W / "mito/ncbi.gb", W / "mito/OW387151.gb") for r in SeqIO.parse(str(f), "genbank")}
fam = {r["species"]: (r["family"], r["genus"]) for r in csv.DictReader(
    open(W / "Annelida_genome_inventory/Annelida_species.tsv"), delimiter="\t")}
root = W / "tools_root"
(root / "gb").mkdir(parents=True, exist_ok=True)
rows = []
for p in pairs:
    r = recs[p["mito"]]
    SeqIO.write(r, str(root / "gb" / f"{r.id}.gb"), "genbank")
    write_bare(r, root / "bare")
    f, g = fam.get(p["organism"], ("", ""))
    rows.append({**{c: "" for c in COLUMNS}, "accession": r.id, "phylum": "Annelida", "family": f, "genus": g,
                 "genetic_code": 5, "topology": r.annotations.get("topology", "circular"), "length": len(r),
                 "mitoz_clade": MITOZ_CLADE.get("Annelida") or "", "split": "test", "retrieved": "2026-11-11"})
with open(root / "set.tsv", "w", newline="") as fh:
    cw = csv.DictWriter(fh, COLUMNS, delimiter="\t", lineterminator="\n")
    cw.writeheader()
    cw.writerows(rows)
print(f"{len(rows)} records -> {root}", file=sys.stderr)
