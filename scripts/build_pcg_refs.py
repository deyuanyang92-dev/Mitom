"""Reference proteins of the 13 mitochondrial protein genes with their lineages, for mito_pcg_hmm.py
(models/pcg_refs.faa.gz): the profile HMM of a missing gene is built from the references nearest the record in the
taxonomy.

Sources: GenBank files of annotated animal mitogenomes (the RefSeq download of scripts/build_orf_library.py, plus
any INSDC sets, e.g. all NCBI annelids). Per gene, a CDS is used when its translation (its own /transl_table, else
the record's) has no internal stop and is 0.5-2 x the gene's usual length; at most --per_genus records per genus
and gene (the first in file order: RefSeq first); identical proteins of one genus once.

Header: >accession|gene|genetic code|lineage (taxonomy joined by ';', spaces as '_').
Usage: python3 scripts/build_pcg_refs.py <GenBank files or folders ...> [--out models/pcg_refs.faa.gz] [--per_genus 2]
"""
import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mito_common import iter_gb_files  # noqa: E402
from mito_names import PCG, gene_name  # noqa: E402

USUAL = {"cox1": 513, "cox2": 228, "cox3": 259, "atp6": 225, "atp8": 52, "cob": 379, "nad1": 317, "nad2": 335,
         "nad3": 115, "nad4": 449, "nad4l": 95, "nad5": 580, "nad6": 160}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("inputs", nargs="+")
    p.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "models" / "pcg_refs.faa.gz"))
    p.add_argument("--per_genus", type=int, default=2)
    a = p.parse_args()
    files = [f for x in a.inputs for f in ([Path(x)] if Path(x).is_file() else iter_gb_files(Path(x)))]
    kept = defaultdict(int)                        # (genus, gene) -> records kept
    seen = set()                                   # (genus, gene, protein)
    accs = set()
    n = 0
    with gzip.open(a.out, "wt") as out:
        for f in files:
            for rec in SeqIO.parse(str(f), "genbank"):
                if rec.id in accs:
                    continue
                accs.add(rec.id)
                tax = rec.annotations.get("taxonomy", [])
                if "Metazoa" not in tax:
                    continue
                genus = rec.annotations.get("organism", "").split(" ")[0]
                codes = [c.qualifiers.get("transl_table", [""])[0] for c in rec.features if c.type == "CDS"]
                codes = [int(c) for c in codes if str(c).isdigit()]
                if not codes:
                    continue
                default = max(set(codes), key=codes.count)
                lineage = ";".join(t.replace(" ", "_") for t in tax + [rec.annotations.get("organism", "")])
                done = set()
                for c in rec.features:
                    g = gene_name(c.qualifiers, "CDS") if c.type == "CDS" else ""
                    if g not in PCG or g in done or kept[(genus, g)] >= a.per_genus:
                        continue
                    code = int(c.qualifiers.get("transl_table", [default])[0])
                    try:
                        nt = c.extract(rec.seq)
                        prot = str(nt[: len(nt) // 3 * 3].translate(table=code)).rstrip("*")
                    except Exception:
                        continue
                    if "*" in prot or not 0.5 * USUAL[g] <= len(prot) <= 2 * USUAL[g] or (genus, g, prot) in seen:
                        continue
                    seen.add((genus, g, prot))
                    kept[(genus, g)] += 1
                    done.add(g)
                    out.write(f">{rec.id}|{g}|{code}|{lineage}\n{prot}\n")
                    n += 1
    print(f"{len(accs)} records, {n} proteins -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
