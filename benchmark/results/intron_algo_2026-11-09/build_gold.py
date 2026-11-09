"""Rebuild the intron gold standard of 2026-10-29 from NCBI (the gold GenBank records are on the other machine).

Records: ../intron_gold_2026-10-29/records.tsv (193 with introns) and e2e/groups.tsv (76 intron-free 'neg' + KT989324.1),
fetched from NCBI (efetch gbwithparts, tool=mitom, no e-mail). Per gold record every NCBI intron feature is removed and
one intron feature is written per gold row of introns.tsv with gold_handling 'keep' (gold coordinates, /note
"gold_class=<class>"); the CDS joins stay as NCBI has them (eval_introns.py counts a join gap only when no intron
feature covers >= 50 % of it, so a moved gold intron replaces NCBI's gap). Check: eval_introns.py of the gold against
itself gives every gold intron exact and no other truth intron.

Usage: python3 build_gold.py <work>     -> <work>/gb/<acc>.gb, <work>/bare/<acc>/<acc>.gbf, <work>/set.tsv
"""
import csv
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

HERE = Path(__file__).resolve().parent
GOLD = HERE.parent / "intron_gold_2026-10-29"
sys.path.insert(0, str(HERE.parents[1]))
from mito_sample import COLUMNS, MITOZ_CLADE, write_bare  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
work = Path(sys.argv[1])
raw, gb = work / "ncbi", work / "gb"
raw.mkdir(parents=True, exist_ok=True)
gb.mkdir(exist_ok=True)

groups = {r["accession"]: r for r in csv.DictReader(open(GOLD / "e2e" / "groups.tsv"), delimiter="\t")}
rows = defaultdict(list)
for r in csv.DictReader(open(GOLD / "introns.tsv"), delimiter="\t"):
    rows[r["accession"]].append(r)
meta = {r["accession"]: r for r in csv.DictReader(open(GOLD / "records.tsv"), delimiter="\t")}
accs = sorted(groups)
todo = [a for a in accs if not (raw / f"{a}.gb").exists()]
for i in range(0, len(todo), 50):
    txt = subprocess.run(["curl", "-s", "--max-time", "300", f"{EUTILS}?db=nuccore&tool=mitom&rettype=gbwithparts"
                          f"&retmode=text&id={','.join(todo[i:i + 50])}"], capture_output=True, text=True).stdout
    for chunk in txt.split("//\n"):
        if "LOCUS" not in chunk:
            continue
        chunk = chunk[chunk.index("LOCUS"):]
        acc = next(ln.split()[1] for ln in chunk.splitlines() if ln.startswith("VERSION"))
        (raw / f"{acc}.gb").write_text(chunk + "//\n")
    time.sleep(1)


def loc(coord: str, strand: int, n: int):
    s, e = (int(x) for x in coord.split()[0].split(".."))
    if s <= e:
        return FeatureLocation(s - 1, e, strand=strand)
    a, b = FeatureLocation(s - 1, n, strand=strand), FeatureLocation(0, e, strand=strand)
    return CompoundLocation([a, b] if strand == 1 else [b, a])


out = []
for acc in accs:
    p = raw / f"{acc}.gb"
    if not p.exists():
        print(f"missing {acc}", file=sys.stderr)
        continue
    rec = SeqIO.read(str(p), "genbank")
    n = len(rec)
    if acc in rows:
        rec.features = [f for f in rec.features if f.type != "intron"]
        for g in rows[acc]:
            if g["gold_handling"] != "keep":
                continue
            strand = 1 if g["strand"] == "+" else -1
            rec.features.append(SeqFeature(loc(g["gold"], strand, n), type="intron", qualifiers={
                "gene": [g["gene"]], "note": [f"gold_class={g['class']}; NCBI {g['ncbi']}; {g['reason'][:200]}"]}))
        rec.features.sort(key=lambda f: (int(f.location.start), f.type != "source", f.type != "gene"))
    SeqIO.write(rec, str(gb / f"{acc}.gb"), "genbank")
    write_bare(rec, work / "bare")
    tax = rec.annotations.get("taxonomy", [])
    phylum = groups[acc]["phylum"]
    codes = sorted({f.qualifiers["transl_table"][0] for f in rec.features
                    if f.type == "CDS" and "transl_table" in f.qualifiers})
    code = (meta.get(acc, {}).get("genetic_code") or (codes[0] if codes else "")) or \
        ("4" if phylum in ("Cnidaria", "Porifera", "Placozoa") else "5")
    out.append({**{c: "" for c in COLUMNS}, "accession": acc, "phylum": phylum, "genetic_code": code,
                "topology": rec.annotations.get("topology", "circular"), "length": n,
                "mitoz_clade": MITOZ_CLADE.get(phylum) or "", "split": "test", "retrieved": "2026-11-09"})
with open(work / "set.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, COLUMNS, delimiter="\t")
    w.writeheader()
    w.writerows(out)
with open(work / "groups.tsv", "w") as fh:
    fh.write("accession\tgroup\tphylum\n" + "".join(f"{a}\t{groups[a]['group']}\t{groups[a]['phylum']}\n" for a in accs))
print(f"{len(out)} records", file=sys.stderr)
