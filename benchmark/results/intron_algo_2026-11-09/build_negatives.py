"""Intron-free animal mitogenomes for the false-intron check: NCBI RefSeq complete mitochondrial genomes without an
intron feature, sampled per phylum (every phylum, the large ones capped), records of the other test sets left out.

NCBI esearch per phylum: <phylum>[Organism] AND mitochondrion[filter] AND refseq[filter] AND 10000:100000[SLEN]
NOT intron[FKEY]; ids spread evenly over the result list (retstart steps), then efetch gbwithparts (tool=mitom, no e-mail).
Records whose CDS joins leave a gap >= 142 nt (an intron NCBI gives only as a join) are dropped.

Usage: python3 build_negatives.py <work> [--per_phylum 40]   -> <work>/gb, <work>/bare, <work>/set.tsv
"""
import csv
import re
import subprocess
import sys
import time
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from mito_sample import COLUMNS, MITOZ_CLADE, write_bare  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PHYLA = {"Arthropoda": 60, "Chordata": 50, "Mollusca": 40, "Annelida": 40, "Nematoda": 25, "Platyhelminthes": 25,
         "Cnidaria": 30, "Porifera": 20, "Echinodermata": 20, "Bryozoa": 15, "Brachiopoda": 10, "Nemertea": 10,
         "Rotifera": 10, "Tardigrada": 6, "Onychophora": 6, "Chaetognatha": 6, "Hemichordata": 6,
         "Xenacoelomorpha": 6, "Ctenophora": 6, "Priapulida": 4, "Phoronida": 4, "Entoprocta": 4,
         "Gastrotricha": 4, "Kinorhyncha": 4, "Nematomorpha": 4, "Cycliophora": 2, "Placozoa": 4}
work = Path(sys.argv[1])
(work / "gb").mkdir(parents=True, exist_ok=True)
skip = set()
for f in (HERE.parent / "intron_gold_2026-10-29" / "e2e" / "groups.tsv", HERE.parent / "intron_ncbi_2026-11-09" / "ncbi_records.tsv"):
    skip |= {r[0].split(".")[0] for r in csv.reader(open(f), delimiter="\t")}


def esearch(term, retstart=0, retmax=1):
    x = subprocess.run(["curl", "-s", "--max-time", "60", "-G", f"{EUTILS}/esearch.fcgi", "--data-urlencode", "db=nuccore",
                        "--data-urlencode", "tool=mitom", "--data-urlencode", f"retstart={retstart}",
                        "--data-urlencode", f"retmax={retmax}", "--data-urlencode", f"term={term}"],
                       capture_output=True, text=True).stdout
    n = re.search(r"<Count>(\d+)</Count>", x)
    return int(n.group(1)) if n else 0, re.findall(r"<Id>(\d+)</Id>", x)


ids = []
for phylum, cap in PHYLA.items():
    term = (f"{phylum}[Organism] AND mitochondrion[filter] AND refseq[filter] AND 10000:100000[SLEN] "
            f"NOT intron[FKEY]")
    n, _ = esearch(term, 0, 0)
    time.sleep(0.4)
    if not n:
        continue
    k = min(cap, n)
    _, all_ids = esearch(term, 0, min(n, 5000))
    step = len(all_ids) / k
    ids += [all_ids[int(i * step)] for i in range(k)]
    print(phylum, n, k, file=sys.stderr)
    time.sleep(0.4)
for i in range(0, len(ids), 100):
    txt = subprocess.run(["curl", "-s", "--max-time", "300", f"{EUTILS}/efetch.fcgi?db=nuccore&tool=mitom"
                          f"&rettype=gbwithparts&retmode=text&id={','.join(ids[i:i + 100])}"],
                         capture_output=True, text=True).stdout
    (work / f"chunk_{i // 100:03d}.gb").write_text(txt)
    time.sleep(1)

out = []
for f in sorted(work.glob("chunk_*.gb")):
    for rec in SeqIO.parse(str(f), "genbank"):
        if rec.id.split(".")[0] in skip:
            continue
        gaps = [int(b.start) - int(a.end) for c in rec.features if c.type == "CDS" and len(c.location.parts) > 1
                for a, b in zip(c.location.parts, c.location.parts[1:]) if a.strand == b.strand == 1]
        if any(g >= 142 for g in gaps):
            continue
        tax = rec.annotations.get("taxonomy", [])
        phylum = next((t for t in tax if t in PHYLA), "")
        codes = sorted({c.qualifiers["transl_table"][0] for c in rec.features
                        if c.type == "CDS" and "transl_table" in c.qualifiers})
        if not codes:
            continue
        SeqIO.write(rec, str(work / "gb" / f"{rec.id}.gb"), "genbank")
        write_bare(rec, work / "bare")
        out.append({**{c: "" for c in COLUMNS}, "accession": rec.id, "phylum": phylum, "genetic_code": codes[0],
                    "topology": rec.annotations.get("topology", "circular"), "length": len(rec),
                    "mitoz_clade": MITOZ_CLADE.get(phylum) or "", "split": "test", "retrieved": "2026-11-09"})
with open(work / "set.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, COLUMNS, delimiter="\t")
    w.writeheader()
    w.writerows(out)
print(f"{len(out)} intron-free records", file=sys.stderr)
