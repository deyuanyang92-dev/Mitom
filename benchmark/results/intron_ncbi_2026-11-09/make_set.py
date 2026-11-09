"""NCBI animal mitochondrial records with introns that the intron gold standard (2026-10-29) does not hold.

1. NCBI esearch (tool=mitom, no e-mail): txid33208[Organism:exp] AND mitochondrion[filter] AND intron[FKEY] AND
   10000:100000[SLEN] (animal mitochondrial records of genome size with an intron feature); efetch gbwithparts.
2. Each record: in the gold standard (../intron_gold_2026-10-29/records.tsv), an INSDC copy of a gold NC_ record (the
   NC_ COMMENT "identical to ..."), an INSDC copy of a new NC_ record (the NC_ is tested), or tested.
3. The tested records as a mito_bench_run.py root: <root>/gb/<acc>.gb, <root>/bare/<acc>/<acc>.gbf (source only),
   <root>/set.tsv (genetic code from the CDS transl_table).

Usage: python3 make_set.py <work> [--fetch]     (work/gb/chunk_*.gb are fetched with --fetch, else reused)
Writes ncbi_records.tsv (here; every record and its status) and <work>/root/.
"""
import csv
import re
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from mito_sample import COLUMNS, MITOZ_CLADE, write_bare  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
QUERY = "txid33208[Organism:exp] AND mitochondrion[filter] AND intron[FKEY] AND 10000:100000[SLEN]"
PHYLA = ("Cnidaria", "Porifera", "Placozoa", "Annelida", "Bryozoa", "Chordata", "Mollusca", "Arthropoda")
work = Path(sys.argv[1])
(work / "gb").mkdir(parents=True, exist_ok=True)

if "--fetch" in sys.argv:
    xml = subprocess.run(["curl", "-s", "-G", f"{EUTILS}/esearch.fcgi", "--data-urlencode", "db=nuccore",
                          "--data-urlencode", "tool=mitom", "--data-urlencode", "retmax=5000",
                          "--data-urlencode", f"term={QUERY}"], capture_output=True, text=True).stdout
    ids = re.findall(r"<Id>(\d+)</Id>", xml)
    for i in range(0, len(ids), 100):
        gb = subprocess.run(["curl", "-s", f"{EUTILS}/efetch.fcgi?db=nuccore&tool=mitom&rettype=gbwithparts"
                             f"&retmode=text&id={','.join(ids[i:i + 100])}"], capture_output=True, text=True).stdout
        (work / "gb" / f"chunk_{i // 100:03d}.gb").write_text(gb)
        time.sleep(1)

gold = {r["accession"].split(".")[0] for r in csv.DictReader(
    open(HERE.parent / "intron_gold_2026-10-29" / "records.tsv"), delimiter="\t")}
recs = {r.id: r for f in sorted((work / "gb").glob("chunk_*.gb")) for r in SeqIO.parse(str(f), "genbank")}
ref_of = {}                                                       # INSDC accession -> its NC_ copy
for a, r in recs.items():
    m = re.search(r"identical to\s+([A-Z]{1,2}_?\d+)", r.annotations.get("comment", "").replace("\n", " "))
    if a.startswith("NC_") and m:
        ref_of[m.group(1)] = a.split(".")[0]
nc_in_set = {a.split(".")[0] for a in recs if a.startswith("NC_")}

root = work / "root"
(root / "gb").mkdir(parents=True, exist_ok=True)
out, rows = [], []
for a, r in sorted(recs.items()):
    b = a.split(".")[0]
    tax = r.annotations.get("taxonomy", [])
    phy = next((t for t in tax if t in PHYLA), tax[1] if len(tax) > 1 else "")
    gold_nc = {k for k, v in ref_of.items() if v == b}            # this NC_ is the copy of these INSDC records
    if b in gold or gold_nc & gold:
        status = "in_gold"
    elif b in ref_of and ref_of[b] in gold:
        status = "insdc_copy_of_gold"
    elif b in ref_of and ref_of[b] in nc_in_set:
        status = "insdc_copy_of_tested_nc"
    else:
        status = "tested"
    introns = [f for f in r.features if f.type == "intron"]
    rows.append([a, phy, tax[tax.index(phy) + 1] if phy in tax and tax.index(phy) + 1 < len(tax) else "",
                 r.annotations.get("organism", ""), len(r), r.annotations.get("topology", ""), len(introns),
                 ",".join(sorted({f.qualifiers.get("gene", ["?"])[0] for f in introns})), status,
                 r.annotations.get("date", "")])
    if status != "tested":
        continue
    codes = sorted({f.qualifiers["transl_table"][0] for f in r.features
                    if f.type == "CDS" and "transl_table" in f.qualifiers})
    code = codes[0] if codes else ("4" if phy in ("Cnidaria", "Porifera", "Placozoa") else "5")
    SeqIO.write(r, str(root / "gb" / f"{a}.gb"), "genbank")
    write_bare(r, root / "bare")
    out.append({**{c: "" for c in COLUMNS}, "accession": a, "phylum": phy, "genetic_code": code,
                "topology": r.annotations.get("topology", "circular"), "length": len(r),
                "mitoz_clade": MITOZ_CLADE.get(phy) or "", "split": "test", "retrieved": "2026-11-09"})
with open(root / "set.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, COLUMNS, delimiter="\t")
    w.writeheader()
    w.writerows(out)
with open(HERE / "ncbi_records.tsv", "w") as fh:
    fh.write("accession\tphylum\tclass\torganism\tlength\ttopology\tn_intron_features\tintron_genes\tstatus\tdate\n")
    fh.writelines("\t".join(map(str, x)) + "\n" for x in rows)
print(Counter(x[8] for x in rows), Counter(x["phylum"] for x in out), file=sys.stderr)
