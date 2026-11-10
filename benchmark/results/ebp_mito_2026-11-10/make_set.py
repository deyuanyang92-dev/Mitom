"""Test set of genome-project mitogenomes (NCBI title 'genome assembly, organelle: mitochondrion': DToL, ERGA, VGP ...,
mostly PacBio HiFi + MitoHiFi, so complete, circular sequences).

From <work>/gb/chunk_*.gb (fetch.py): every record's phylum, genus, length, topology and the protein genes the project
annotated (13 standard names; most records have no CDS at all). Test set: per phylum up to --per_phylum circular
records of >= 12 kb (complete genomes: every gene of the lineage should be found; genera spread, one record per
genus first; the project's annotation is not needed), written as a
mito_bench_run.py root: <work>/root/gb/<acc>.gb, bare/<acc>/<acc>.gbf (source only), set.tsv.
Usage: python3 make_set.py <work> [--per_phylum 12]  > records.tsv
"""
import argparse
import csv
import random
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

H = Path(__file__).resolve()
sys.path.insert(0, str(H.parents[3]))
sys.path.insert(0, str(H.parents[2]))
from mito_names import PCG, gene_name  # noqa: E402
from mito_sample import COLUMNS, DEFAULT_CODE, MITOZ_CLADE, write_bare  # noqa: E402

PHYLA = ["Porifera", "Cnidaria", "Ctenophora", "Placozoa", "Xenacoelomorpha", "Chordata", "Echinodermata",
         "Hemichordata", "Arthropoda", "Nematoda", "Tardigrada", "Onychophora", "Priapulida", "Nematomorpha",
         "Kinorhyncha", "Annelida", "Mollusca", "Brachiopoda", "Bryozoa", "Nemertea", "Platyhelminthes", "Rotifera",
         "Phoronida", "Entoprocta", "Chaetognatha", "Gastrotricha", "Cycliophora"]

p = argparse.ArgumentParser()
p.add_argument("work")
p.add_argument("--per_phylum", type=int, default=25)
a = p.parse_args()
W = Path(a.work)
info, recs = [], {}
for f in sorted((W / "gb").glob("chunk_*.gb")):
    for r in SeqIO.parse(str(f), "genbank"):
        tax = r.annotations.get("taxonomy", [])
        genes = {gene_name(x.qualifiers, "CDS") for x in r.features if x.type == "CDS"} - {""}
        codes = [x.qualifiers.get("transl_table", [""])[0] for x in r.features if x.type == "CDS"]
        codes = [int(c) for c in codes if str(c).isdigit()]
        ph = next((x for x in PHYLA if x in tax), "")
        row = {"accession": r.id, "phylum": ph, "genus": r.annotations.get("organism", "").split(" ")[0],
               "organism": r.annotations.get("organism", ""), "length": len(r),
               "topology": r.annotations.get("topology", ""), "n_cds": sum(x.type == "CDS" for x in r.features),
               "genes": len(genes), "atp8": "atp8" in genes, "atp6": "atp6" in genes,
               "code": max(set(codes), key=codes.count) if codes else DEFAULT_CODE[ph],
               "family": next((t for t in tax if t.endswith("idae")), "")}
        info.append(row)
        recs[r.id] = r
w = csv.DictWriter(sys.stdout, list(info[0]), delimiter="\t", lineterminator="\n")
w.writeheader()
w.writerows(info)
rnd = random.Random(1)
pick = []
by = defaultdict(list)
for x in info:
    if x["phylum"] and x["topology"] == "circular" and x["length"] >= 12000:
        by[x["phylum"]].append(x)
for ph, xs in sorted(by.items()):
    rnd.shuffle(xs)
    seen, first, rest = set(), [], []
    for x in xs:
        (rest if x["genus"] in seen else first).append(x)
        seen.add(x["genus"])
    pick += (first + rest)[: a.per_phylum]
root = W / "root"
(root / "gb").mkdir(parents=True, exist_ok=True)
rows = []
for x in pick:
    r = recs[x["accession"]]
    SeqIO.write(r, str(root / "gb" / f"{r.id}.gb"), "genbank")
    write_bare(r, root / "bare")
    rows.append({**{c: "" for c in COLUMNS}, "accession": r.id, "phylum": x["phylum"], "genus": x["genus"],
                 "family": x["family"], "genetic_code": x["code"], "topology": "circular", "length": x["length"],
                 "mitoz_clade": MITOZ_CLADE.get(x["phylum"]) or "", "split": "test", "retrieved": "2026-11-10"})
with open(root / "set.tsv", "w", newline="") as fh:
    cw = csv.DictWriter(fh, COLUMNS, delimiter="\t", lineterminator="\n")
    cw.writeheader()
    cw.writerows(rows)
print(f"{len(info)} records; test set {len(rows)} in {len(by)} phyla", file=sys.stderr)
