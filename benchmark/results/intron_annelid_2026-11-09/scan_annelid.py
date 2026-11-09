"""Every NCBI annelid mitochondrial record of 8-100 kb: which have an intron (intron feature or a CDS join gap
>= 142 nt, benchmark/eval_introns.introns_of), and are they in the intron gold standard (2026-10-29)?

NCBI esearch Annelida[Organism] AND mitochondrion[filter] AND 8000:100000[SLEN] (tool=mitom, no e-mail), efetch
gbwithparts, 150 records per request, into <work>/chunk_*.gb.
Usage: python3 scan_annelid.py <work> [--fetch]  > annelid_records.tsv
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
from eval_introns import introns_of  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
work = Path(sys.argv[1])
work.mkdir(parents=True, exist_ok=True)
if "--fetch" in sys.argv:
    x = subprocess.run(["curl", "-s", "-G", f"{EUTILS}/esearch.fcgi", "--data-urlencode", "db=nuccore",
                        "--data-urlencode", "tool=mitom", "--data-urlencode", "retmax=5000", "--data-urlencode",
                        "term=Annelida[Organism] AND mitochondrion[filter] AND 8000:100000[SLEN]"],
                       capture_output=True, text=True).stdout
    ids = re.findall(r"<Id>(\d+)</Id>", x)
    for i in range(0, len(ids), 150):
        (work / f"chunk_{i // 150:03d}.gb").write_text(subprocess.run(
            ["curl", "-s", "--max-time", "600", f"{EUTILS}/efetch.fcgi?db=nuccore&tool=mitom&rettype=gbwithparts"
             f"&retmode=text&id={','.join(ids[i:i + 150])}"], capture_output=True, text=True).stdout)
        time.sleep(1)
gold = {r["accession"].split(".")[0] for r in csv.DictReader(
    open(HERE.parent / "intron_gold_2026-10-29" / "records.tsv"), delimiter="\t")}
w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
w.writerow(["accession", "organism", "length", "in_gold", "identical_to", "introns"])
n = 0
for f in sorted(work.glob("chunk_*.gb")):
    for r in SeqIO.parse(str(f), "genbank"):
        n += 1
        it = introns_of(r)
        if not it:
            continue
        m = re.search(r"identical to\s+([A-Z]{1,2}_?\d+)", r.annotations.get("comment", "").replace("\n", " "))
        w.writerow([r.id, r.annotations.get("organism", ""), len(r), "yes" if r.id.split(".")[0] in gold else "no",
                    m.group(1) if m else "", "; ".join(f"{i.gene} {i.end - i.start + 1} nt ({i.source})" for i in it)])
print(f"{n} records scanned", file=sys.stderr)
