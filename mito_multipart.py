#!/usr/bin/env python3
"""mito_multipart.py - records that are parts of one mitochondrial genome (multipartite mtDNA, or fragments).

Some animal mitochondrial genomes are several circles: the sponge Lycopodina hypogea has three chromosomes, nad4l alone
on a 291-bp one (PP657140-PP657142); one HiFi assembly of the nematode Allodiplogaster pararmata gave two circles with
complementary genes (OZ057410 / OZ057411). Read one by one, each record "misses" genes. Records of one genome are
recognised by what the submitter states:
  1. the same assembly (GenBank COMMENT "The assembly <name> is based ..." of genome projects), else
  2. the same BioSample (DBLINK), else
  3. (only with group_by_organism) the same organism in one run.
A group of >= 2 distinct records (identical sequences - an INSDC record and its RefSeq copy - counted once) where
each record lacks protein genes but together they hold at least MIN_UNION of the 13 is reported as "parts of one
genome": a multipartite mitochondrial genome, or one genome in several pieces. A gene complete in one record and
present as a fragment (< FRAGMENT of the usual length) in another is listed too.

Report multipart_report.tsv (one row per group of >= 2 records); `mitom run`: step 9 (NCBI check and summary), a
WARNING per group found. Measured: benchmark/results/ebp_mito_2026-11-10/README.md.

Usage: python3 mito_multipart.py -i <GenBank file or folder> [-o multipart_report.tsv] [--by_organism]
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from mito_common import cli_guard, iter_gb_files, require_path, write_tsv
from mito_names import PCG, gene_name

MIN_UNION = 12              # protein genes the parts hold together (atp8 is lost in several lineages)
FRAGMENT = 0.5              # shorter than this share of the usual length: a fragment of the gene
USUAL = {"cox1": 1540, "cox2": 685, "cox3": 780, "atp6": 678, "atp8": 160, "cob": 1140, "nad1": 954, "nad2": 1008,
         "nad3": 348, "nad4": 1350, "nad4l": 288, "nad5": 1740, "nad6": 483}
COLUMNS = ["group", "basis", "records", "distinct", "genes_each", "together", "missing_together", "verdict",
           "fragments"]


def group_key(rec: SeqRecord, by_organism: bool = False) -> Optional[tuple]:
    comment = rec.annotations.get("comment", "").replace("\n", " ")
    m = re.search(r"The assembly (\S+?) is based", comment)
    if m:
        return ("assembly", m.group(1))
    bs = next((x for x in rec.dbxrefs if x.startswith("BioSample:")), None)
    if bs:
        return ("BioSample", bs.split(":", 1)[1])
    if by_organism and rec.annotations.get("organism"):
        return ("organism", rec.annotations["organism"])
    return None


def genes_of(rec: SeqRecord) -> Dict[str, int]:
    """Standard protein gene -> longest CDS length (nt) in the record."""
    out: Dict[str, int] = {}
    for f in rec.features:
        g = gene_name(f.qualifiers, "CDS") if f.type == "CDS" else ""
        if g in PCG:
            out[g] = max(out.get(g, 0), len(f))
    return out


def kmers(rec: SeqRecord, k: int = 21) -> set:
    """Canonical k-mers of a record (strand and start position do not matter)."""
    s = str(rec.seq).upper()
    s2 = s + s[: k - 1]                                  # (circular)
    out = set()
    for i in range(len(s)):
        x = s2[i: i + k]
        y = x[::-1].translate(str.maketrans("ACGT", "TGCA"))
        out.add(min(x, y))
    return out


def same_sequence(a: set, b: set) -> bool:
    """Copies of one sequence (an INSDC record and its RefSeq copy, another rotation or strand): >= 99 % shared."""
    return len(a & b) >= 0.99 * max(len(a), len(b))


def analyse(records: Iterable[SeqRecord], by_organism: bool = False) -> List[dict]:
    groups: Dict[tuple, List[SeqRecord]] = defaultdict(list)
    for r in records:
        k = group_key(r, by_organism)
        if k:
            groups[k].append(r)
    rows = []
    for (basis, name), recs in sorted(groups.items()):
        if len(recs) < 2:
            continue
        distinct, kms = [], []
        for r in recs:
            km = kmers(r)
            if not any(same_sequence(km, o) for o in kms):
                kms.append(km)
                distinct.append(r)
        each = {r.id: genes_of(r) for r in distinct}
        full = {r.id: {g for g, n in each[r.id].items() if n >= FRAGMENT * USUAL[g]} for r in distinct}
        together = set().union(*full.values()) if full else set()
        frags = sorted({f"{g} ({rid})" for rid, gs in each.items() for g, n in gs.items()
                        if n < FRAGMENT * USUAL[g] and any(g in full[o] for o in full if o != rid)})
        if len(distinct) < 2:
            verdict = "copies of one sequence"
        elif len(together) >= MIN_UNION and all(len(full[r.id]) < len(together) for r in distinct):
            verdict = ("parts of one genome: the genes are complete only together (multipartite mitochondrial "
                       "genome, or one genome in pieces)")
        elif len(together) >= MIN_UNION:
            verdict = "one record holds all genes; the others are extra copies or pieces"
        else:
            verdict = "incomplete even together"
        rows.append({"group": name, "basis": basis, "records": ",".join(r.id for r in recs),
                     "distinct": len(distinct),
                     "genes_each": "; ".join(f"{r.id}:{len(full[r.id])}" for r in distinct),
                     "together": len(together), "missing_together": ",".join(sorted(set(PCG) - together)),
                     "verdict": verdict, "fragments": ", ".join(frags)})
    return rows


def main(argv: Sequence[str] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="annotated GenBank file or folder")
    p.add_argument("-o", "--out", default="multipart_report.tsv")
    p.add_argument("--by_organism", action="store_true",
                   help="also group records of the same organism that name no assembly / BioSample")
    a = p.parse_args(argv)
    src = require_path(a.input)
    recs = [r for f in ([src] if src.is_file() else iter_gb_files(src)) for r in SeqIO.parse(str(f), "genbank")]
    rows = analyse(recs, a.by_organism)
    write_tsv(a.out, COLUMNS, rows)
    for r in rows:
        if r["verdict"].startswith("parts of one genome"):
            print(f"[mito_multipart] WARNING {r['basis']} {r['group']} ({r['records']}): {r['verdict']}; together "
                  f"{r['together']} of 13 protein genes", file=sys.stderr)
    print(f"[mito_multipart] {len(rows)} groups of >= 2 records -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(cli_guard(main))
