#!/usr/bin/env python3
"""build_ref_taxa.py - models/ref_taxa.tsv.gz for --clade auto without a taxonomy (mito_classify.py): the NCBI lineage
and mitochondrial genetic code (CDS /transl_table) of every genome whose proteins are in the MITOS2 reference data
refseq89m (featureProt/*.fas headers: >ACCESSION:gene-... Species name).

The lineages come from a local GenBank file of RefSeq mitogenomes (no download here), e.g. the output of
  mitom fetch --taxon Metazoa --refseq_only -o refpool
Genomes of refseq89m that are not in that file are listed on the screen and left out (their hits are not counted).

Usage: python scripts/build_ref_taxa.py --refdir ~/mitos_refs/refseq89m --gb refpool/Metazoa.gb
"""
from __future__ import annotations

import argparse
import gzip
from collections import Counter
from pathlib import Path

from Bio import SeqIO

OUT = Path(__file__).resolve().parent.parent / "models" / "ref_taxa.tsv.gz"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refdir", required=True, help="refseq89m folder (with featureProt/)")
    ap.add_argument("--gb", required=True, help="GenBank file of RefSeq mitogenomes with their taxonomy")
    ap.add_argument("-o", "--out", default=str(OUT))
    a = ap.parse_args()
    accs = set()
    for fas in sorted((Path(a.refdir).expanduser() / "featureProt").glob("*.fas")):
        with fas.open() as fh:
            accs.update(line[1:].split(":")[0] for line in fh if line.startswith(">"))
    rows = {}
    for rec in SeqIO.parse(a.gb, "genbank"):
        acc = rec.id.split(".")[0]
        if acc not in accs or acc in rows:
            continue
        codes = Counter(f.qualifiers["transl_table"][0] for f in rec.features
                        if f.type == "CDS" and "transl_table" in f.qualifiers)
        tax = [t for t in rec.annotations.get("taxonomy", []) if t]
        tax = tax[tax.index("Metazoa"):] if "Metazoa" in tax else tax
        rows[acc] = (codes.most_common(1)[0][0] if codes else "", ";".join(tax))
    with gzip.open(a.out, "wt") as fh:
        fh.write("accession\tgenetic_code\tlineage\n")
        for acc in sorted(rows):
            fh.write(f"{acc}\t{rows[acc][0]}\t{rows[acc][1]}\n")
    missing = sorted(accs - set(rows))
    print(f"{len(rows)} of {len(accs)} reference genomes written to {a.out}; not in {a.gb}: {len(missing)}"
          + (f" (e.g. {', '.join(missing[:5])})" if missing else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
