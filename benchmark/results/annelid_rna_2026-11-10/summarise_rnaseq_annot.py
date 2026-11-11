#!/usr/bin/env python3
"""Summary of mito_rnaseq_annot.py over the annelid pairs: verdicts per step for the learning set (clean_pairs.txt,
the 23 pairs the rules were learned on) and the test set (all other pairs); each corrected protein's length against
the median of annotated relatives of its family (models/pcg_refs.faa.gz, NCBI family from the inventory).

Usage: summarise_rnaseq_annot.py <work dir> -o rnaseq_annot_summary.tsv
"""
import argparse
import csv
import gzip
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
GENE = {"ATP6": "atp6", "ATP8": "atp8", "COX1": "cox1", "COX2": "cox2", "COX3": "cox3", "CYTB": "cob", "ND1": "nad1",
        "ND2": "nad2", "ND3": "nad3", "ND4": "nad4", "ND4L": "nad4l", "ND5": "nad5", "ND6": "nad6"}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("work")
    p.add_argument("-o", "--out", required=True)
    a = p.parse_args()
    w = Path(a.work)
    clean = {x.strip() for x in open(HERE / "clean_pairs.txt") if x.strip()}
    pairs = list(csv.DictReader(open(HERE / "datasets.tsv"), delimiter="\t"))
    fam = {}
    for r in csv.DictReader(open(w / "Annelida_genome_inventory/Annelida_species.tsv"), delimiter="\t"):
        fam[r.get("species", r.get("ScientificName", ""))] = r.get("family", "")
    refs = defaultdict(list)
    for r in SeqIO.parse(gzip.open(HERE.parents[2] / "models/pcg_refs.faa.gz", "rt"), "fasta"):
        d = r.description.split("|")
        if len(d) >= 4:
            refs[d[1].lower()].append((d[3], len(r.seq)))
    counts = {"learning": Counter(), "test": Counter()}
    rows = []
    for pr in pairs:
        acc = pr["mito"].split(".")[0]
        ev = w / "rnaseq_annot" / acc / "rnaseq_evidence.tsv"
        if not ev.exists():
            continue
        group = "learning" if pr["mito"] in clean else "test"
        gb = SeqIO.read(w / "rnaseq_annot" / acc / "rnaseq_annotated.gb", "genbank")
        prot = {f.qualifiers.get("gene", [""])[0]: len(f.qualifiers.get("translation", [""])[0])
                for f in gb.features if f.type == "CDS"}
        f_ = fam.get(pr["organism"], "")
        for e in csv.DictReader(open(ev), delimiter="\t"):
            counts[group][(e["step"], e["verdict"])] += 1
            if e["verdict"] in ("corrected", "flagged"):
                g = e["feature"].split()[0].split(",")[0]
                xs = [n for lin, n in refs.get(GENE.get(g, "?"), []) if f_ and f_ in lin]
                rows.append({"group": group, "mito": pr["mito"], "organism": pr["organism"], "family": f_,
                             "step": e["step"], "feature": e["feature"], "verdict": e["verdict"],
                             "protein_aa_now": prot.get(g, ""),
                             "family_median_aa": statistics.median(xs) if xs else "", "n_relatives": len(xs),
                             "detail": e["detail"]})
    with open(a.out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, list(rows[0]) if rows else ["group"], delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows)
    steps = sorted({k for c in counts.values() for k in c})
    print("step\tverdict\tlearning\ttest")
    for k in steps:
        print(f"{k[0]}\t{k[1]}\t{counts['learning'][k]}\t{counts['test'][k]}")


if __name__ == "__main__":
    main()
