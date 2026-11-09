"""Truth for the recall test: every CDS of a RefSeq animal mitogenome (>= 10 CDS) that is not one of the 13 genes and
>= 41 codons (mutS, atp9, F-ORF, homing endonucleases, hypothetical proteins ...), from the library download
(scripts/build_orf_library.py <work>). The records are written without these CDS (and their gene / mRNA features)
to <work>/recall/in/<acc>.gb, the CDS to <work>/recall/truth.tsv.
Usage: python3 make_truth.py <work>
"""
import csv
import glob
import sys
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from mito_orf import gene_key  # noqa: E402

W = Path(sys.argv[1])
(W / "recall" / "in").mkdir(parents=True, exist_ok=True)
rows = []
for f in sorted(glob.glob(str(W / "refseq" / "chunk_*.gb"))):
    for r in SeqIO.parse(f, "genbank"):
        if sum(x.type == "CDS" for x in r.features) < 10:
            continue
        odd = [x for x in r.features if x.type == "CDS" and not gene_key(x) and len(x) >= 123]
        if not odd:
            continue
        spans = {(int(x.location.start), int(x.location.end)) for x in odd}
        r.features = [x for x in r.features if not (x.type in ("CDS", "gene", "mRNA") and
                                                    (int(x.location.start), int(x.location.end)) in spans)]
        tax = r.annotations.get("taxonomy", [])
        for x in odd:
            rows.append([r.id, tax[tax.index("Metazoa") + 1] if "Metazoa" in tax else "", r.annotations["organism"],
                         x.qualifiers.get("gene", [""])[0], x.qualifiers.get("product", [""])[0],
                         int(x.location.start) + 1, int(x.location.end), x.location.strand, len(x) // 3,
                         len(x.location.parts)])
        SeqIO.write(r, str(W / "recall" / "in" / f"{r.id}.gb"), "genbank")
with open(W / "recall" / "truth.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["record", "phylum", "organism", "gene", "product", "start", "end", "strand", "codons", "parts"])
    w.writerows(rows)
