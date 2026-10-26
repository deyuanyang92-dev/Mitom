#!/usr/bin/env python3
"""Parse downloaded RefSeq mitogenomes -> records.tsv + trnas.tsv + per-record FASTA, and a train/test split.

Split: test = records NOT in MITOS2 refseq89m (so MiTFi never saw them) AND whose family is held out
(about 25 % of families with such records, hashed); train = every record of the remaining families.
No family is in both sets, so per-clade models are judged on unseen families.

Usage: build_dataset.py <fetch_dir> <refseq89m_accs.txt> <out_dir>
"""
import csv
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
PHYLA = {"Annelida", "Arthropoda", "Mollusca", "Nematoda", "Platyhelminthes", "Cnidaria", "Porifera", "Echinodermata",
         "Chordata", "Hemichordata", "Bryozoa", "Brachiopoda", "Nemertea", "Rotifera", "Acanthocephala", "Tardigrada",
         "Onychophora", "Chaetognatha", "Ctenophora", "Placozoa", "Xenacoelomorpha", "Priapulida", "Kinorhyncha",
         "Nematomorpha", "Gastrotricha", "Entoprocta", "Phoronida", "Cycliophora", "Gnathostomulida", "Loricifera",
         "Orthonectida", "Dicyemida", "Micrognathozoa", "Sipuncula", "Echiura"}


def rc(s):
    return s.translate(COMP)[::-1]


def aa_of(f):
    txt = " ".join(f.qualifiers.get("product", []) + f.qualifiers.get("gene", []))
    m = re.search(r"tRNA-([A-Z][a-z]{2})", txt)
    if m and m.group(1) in AA3:
        return AA3[m.group(1)]
    m = re.match(r"trn([A-Z])", f.qualifiers.get("gene", [""])[0], re.I)
    return m.group(1).upper() if m and m.group(1).upper() in AA3.values() else ""


def anticodon_of(f):
    for q in f.qualifiers.get("anticodon", []):
        m = re.search(r"seq:([acgtu]{3})", q, re.I)
        if m:
            return m.group(1).lower().replace("u", "t")
    txt = " ".join(f.qualifiers.get("note", []) + f.qualifiers.get("gene", []) + f.qualifiers.get("product", []))
    m = re.search(r"\(([acgtuACGTU]{3})\)", txt) or re.search(r"anticodon[: ]+([acgtuACGTU]{3})", txt)
    return m.group(1).lower().replace("u", "t") if m else ""


def main(fetch, accs, out):
    out = Path(out)
    (out / "fasta").mkdir(parents=True, exist_ok=True)
    in89 = {a.strip() for a in open(accs) if a.strip()}
    recs, trnas = [], []
    for gb in sorted(Path(fetch).glob("batch_*.gb")):
        for r in SeqIO.parse(str(gb), "genbank"):
            acc = r.id
            tax = r.annotations.get("taxonomy", [])
            phylum = next((t for t in tax if t in PHYLA), "")
            ranks = tax[tax.index(phylum) + 1:] if phylum else []
            fam = next((t for t in ranks if t.endswith("idae")), "")
            order = next((t for t in ranks if t.endswith("ida") or t.endswith("ptera") or t.endswith("formes")), "")
            code = ""
            for f in r.features:
                if f.type == "CDS" and "transl_table" in f.qualifiers:
                    code = f.qualifiers["transl_table"][0]
                    break
            seq = str(r.seq).upper()
            if set(seq) - set("ACGTN") and sum(c not in "ACGT" for c in seq) > 0.01 * len(seq):
                continue
            n = 0
            for f in r.features:
                if f.type != "tRNA":
                    continue
                aa = aa_of(f)
                parts = f.location.parts
                s = str(f.location.extract(r.seq))
                trnas.append({"acc": acc, "aa": aa, "anticodon": anticodon_of(f),
                              "gene": f.qualifiers.get("gene", [""])[0],
                              "start": int(f.location.start) + 1, "end": int(f.location.end),
                              "strand": f.location.strand or 1, "parts": len(parts), "length": len(s),
                              "seq": s.upper()})
                n += 1
            sp = r.annotations.get("organism", "")
            recs.append({"acc": acc, "phylum": phylum, "class": ranks[0] if ranks else "", "order": order,
                         "family": fam, "genus": sp.split()[0] if sp else "", "species": sp,
                         "lineage": ";".join(tax), "code": code, "length": len(seq),
                         "topology": r.annotations.get("topology", ""), "n_trna": n,
                         "in_refseq89": int(acc.split(".")[0] in in89)})
            (out / "fasta" / f"{acc}.fa").write_text(f">{acc}\n{seq}\n")
    fam_new = defaultdict(int)
    for x in recs:
        if not x["in_refseq89"]:
            fam_new[x["family"] or x["genus"]] += 1
    for x in recs:
        key = x["family"] or x["genus"]
        held = key in fam_new and int(hashlib.md5(key.encode()).hexdigest(), 16) % 4 == 0
        x["split"] = ("test" if not x["in_refseq89"] else "excluded") if held else "train"
    for name, rows in (("records.tsv", recs), ("trnas.tsv", trnas)):
        with open(out / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    from collections import Counter
    print(len(recs), "records", len(trnas), "tRNAs", Counter(x["split"] for x in recs))


if __name__ == "__main__":
    main(*sys.argv[1:4])
