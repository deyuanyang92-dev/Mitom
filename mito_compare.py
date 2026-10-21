#!/usr/bin/env python3
"""mito_compare.py - compare the original annotation of each record (e.g. as downloaded from NCBI) with the
MitoM annotation of the same sequence, gene by gene.

MitoM rotates / reverse-complements a record so that cox1 comes first; the MitoM genes are mapped back onto the
original sequence first, so both annotations are compared in the original coordinates.

One row per gene (protein genes, rRNAs, tRNAs, introns):
  result  same            same name, strand and both ends
          ends_differ     same gene, at least one end differs (difference_nt = larger end shift)
          name_differs    same place, different gene name
          strand_differs  same place, other strand
          only_original   annotated only in the original record
          only_mitom      annotated only by MitoM
          same_place      tRNA-Leu / tRNA-Ser at the same place, which of the two genes cannot be told (no anticodon)

Usage:  mitom compare -o original_vs_mitom.tsv --original original_annotation/ --mitom final_annotation/ [--genetic_code 5]
        (files are paired by name: original_annotation/<acc>.gb <-> final_annotation/<acc>.gb)
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "benchmark"))
import mito_bench_eval as me  # noqa: E402  (gene normalisation and matching)
from mito_bench_run import LiftError, lift_to_ncbi  # noqa: E402

GB_EXT = {".gb", ".gbk", ".gbf", ".genbank", ".gbff"}
RESULT = {"exact": "same", "le3": "ends_differ", "gt3": "ends_differ", "name_error": "name_differs",
          "strand_error": "strand_differs", "missed": "only_original", "fp": "only_mitom",
          "name_unverifiable": "same_place"}
COLUMNS = ["record", "organism", "type", "gene", "original_location", "mitom_location", "difference_nt", "result"]
KIND = {"PCG": "protein", "rRNA": "rRNA", "tRNA": "tRNA", "intron": "intron"}


def _gene_label(kind: str, name: str) -> str:
    return f"trn{name}" if kind == "tRNA" else name


def compare_record(original, mitom, code: int = 5) -> List[dict]:
    """Rows for one record; `mitom` may be rotated / reverse-complemented relative to `original`."""
    organism = str(original.annotations.get("organism", ""))
    lifted = lift_to_ncbi(mitom, Seq(str(original.seq)))
    rows = []
    for m in me.match(me.norm_features(original, code), me.norm_features(lifted, code)):
        rows.append({"record": original.id, "organism": organism, "type": KIND.get(m["kind"], m["kind"]),
                     "gene": _gene_label(m["kind"], m["name"]), "original_location": m["truth"],
                     "mitom_location": m["pred"], "difference_nt": m["delta"], "result": RESULT[m["class"]]})
    return rows


def gb_by_stem(path: Path) -> Dict[str, Path]:
    files = [path] if path.is_file() else [p for p in path.rglob("*") if p.suffix.lower() in GB_EXT]
    return {p.stem: p for p in sorted(files)}


def compare_dirs(original_dir: Path, mitom_dir: Path, out_tsv: Path, code: int = 5) -> Counter:
    """Write out_tsv for every record present in both directories; returns counts per result
    (+ 'records', 'records_without_original_annotation', 'errors')."""
    orig, mit = gb_by_stem(original_dir), gb_by_stem(mitom_dir)
    rows: List[dict] = []
    n = Counter()
    for stem in sorted(set(orig) & set(mit)):
        try:
            o = SeqIO.read(str(orig[stem]), "genbank")
            m = SeqIO.read(str(mit[stem]), "genbank")
            if not any(f.type in ("CDS", "rRNA", "tRNA") for f in o.features):
                n["records"] += 1                      # nothing to compare: every gene would be 'only_mitom'
                n["records_without_original_annotation"] += 1
                continue
            got = compare_record(o, m, code)
        except (LiftError, ValueError) as exc:
            n["errors"] += 1
            rows.append({"record": stem, "organism": "", "type": "", "gene": "", "original_location": "",
                         "mitom_location": "", "difference_nt": "", "result": f"not_compared: {exc}"})
            continue
        n["records"] += 1
        n.update(r["result"] for r in got)
        rows += got
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return n


def summary_text(n: Counter) -> str:
    genes = sum(n[k] for k in set(RESULT.values()))
    parts = [f"{k} {n[k]}" for k in ("same", "ends_differ", "only_original", "only_mitom", "name_differs",
                                       "strand_differs", "same_place") if n[k]]
    return (f"{n['records']} records, {genes} genes: " + ", ".join(parts)
            + (f"; {n['records_without_original_annotation']} records had no original annotation"
               if n["records_without_original_annotation"] else "")
            + (f"; {n['errors']} records not compared" if n["errors"] else ""))


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--original", required=True, help="Original GenBank file or directory.")
    ap.add_argument("--mitom", required=True, help="MitoM GenBank file or directory.")
    ap.add_argument("-o", "--out", default="original_vs_mitom.tsv")
    ap.add_argument("--genetic_code", type=int, default=5)
    a = ap.parse_args(argv)
    n = compare_dirs(Path(a.original), Path(a.mitom), Path(a.out), a.genetic_code)
    print(summary_text(n))
    print(f"table: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
