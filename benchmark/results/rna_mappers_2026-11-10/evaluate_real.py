#!/usr/bin/env python3
"""Mappers on real RNA-seq (no truth): what mito_rna gets from each mapper's SAM.

Per mapper: reads aligned with MAPQ >= 20 (per 1,000 pairs), spliced (N) alignments, poly(A) reads and sites,
sites at an annotated protein-gene / rRNA 3' end (within 3 nt), protein-gene median depth, the relative depth of each
feature against the median over the mappers (median |log2 ratio| over genes / ORFs / non-coding regions; how far a
mapper's mito_rna table is from the consensus), time and peak memory.

Usage: evaluate_real.py --mitom <Mitom dir> -g genome.gb --maps <dir> --label x [--pairs N] [--header]
"""
from __future__ import annotations

import argparse
import math
import statistics
import sys
from pathlib import Path

COLS = ["label", "mapper", "mapq20_per_1000_pairs", "spliced", "polya_reads", "polya_sites", "sites_at_gene_ends",
        "pcg_median_depth", "library", "relative_depth_vs_consensus", "seconds", "max_mb"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mitom", required=True)
    p.add_argument("-g", "--genbank", required=True)
    p.add_argument("--maps", required=True)
    p.add_argument("--label", default="")
    p.add_argument("--pairs", type=int, default=4000000)
    p.add_argument("--header", action="store_true")
    a = p.parse_args()
    sys.path.insert(0, a.mitom)
    sys.path.insert(0, str(Path(__file__).parent))
    import mito_rna
    from evaluate import time_of
    from Bio import SeqIO
    rec = SeqIO.read(a.genbank, "genbank")
    res = {}
    for sam in sorted(Path(a.maps).glob("*.sam")):
        with open(sam) as fh:
            lines = [ln for ln in fh if ln.startswith("@") or not int(ln.split("\t", 2)[1]) & 0x900]
        spliced = sum(1 for ln in lines if not ln.startswith("@") and "N" in ln.split("\t", 6)[5])
        rows, sites, info = mito_rna.analyse(rec, lines)
        ends = {(r["end"] - 1 if r["strand"] == "+" else r["start"] - 1, r["strand"]) for r in rows
                if r["type"] in ("CDS", "rRNA")}
        at_end = sum(1 for s in sites if any(abs(s["position"] - 1 - e) <= 3 and s["strand"] == st for e, st in ends))
        res[sam.stem] = (rows, sites, info, spliced, at_end, time_of(sam.with_suffix(".time")))
    feats = {}
    for name, (rows, *_r) in res.items():
        for r in rows:
            if r["relative_depth"]:
                feats.setdefault((r["feature"], r["start"]), {})[name] = float(r["relative_depth"])
    cons = {k: statistics.median(v.values()) for k, v in feats.items()}
    if a.header:
        print("\t".join(COLS))
    for name, (rows, sites, info, spliced, at_end, (sec, mb)) in res.items():
        errs = [abs(math.log2((feats[k][name] + 0.01) / (cons[k] + 0.01))) for k in feats if name in feats[k]]
        lib = ("reverse" if info["reverse"] else "stranded") if info["stranded"] else "unstranded"
        print("\t".join(map(str, [a.label, name, f"{1000 * info['reads_used'] / a.pairs:.1f}", spliced,
                                  info["polya_reads"], len(sites), at_end, f"{info['pcg_median_depth']:.0f}", lib,
                                  f"{statistics.median(errs):.3f}" if errs else "", sec, mb])), flush=True)


if __name__ == "__main__":
    main()
