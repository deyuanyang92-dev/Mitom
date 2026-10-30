#!/usr/bin/env python3
"""mito_ncr.py - non-coding regions (NCR) of an annotated mitochondrial genome: the stretches no feature covers.

A mitogenome often holds sizeable stretches no gene, rRNA, tRNA or control region is annotated on - the control
region when it is not named, or sequence that may hold unidentified genes. They do not show in a GenBank record or a
map unless they are features. `mitom run` (step after the curation, --ncr yes by default, --ncr_min_len N) writes
every such stretch of >= N nt as

    misc_feature    <start>..<end>   (join(<start>..<length>,1..<end>) across the origin of a circular record)
                    /note="non-coding region NCR<k>, <length> bp: no gene annotated (MitoM)"

numbered in record order; a joined gene covers its introns, so an intron is no NCR. The maps (mito_plot.py) shade
and label the same stretches. NCRs of an earlier run are replaced (the records can go through it again).

Usage (also `python3 mito_ncr.py -i <GenBank file or folder> -o <out> [--min_len 100]`): add_ncr(rec, min_len).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

MIN_LEN = 100                                 # default shortest NCR (--ncr_min_len, the maps' --min_noncoding)
NOTE = "non-coding region NCR"                # start of the /note that marks MitoM's NCR features


def is_ncr(f: SeqFeature) -> bool:
    """A non-coding region feature written by add_ncr."""
    return f.type == "misc_feature" and any(x.startswith(NOTE) for x in f.qualifiers.get("note", []))


def extent(loc, length: int, circular: bool) -> Tuple[int, int]:
    """(start, length in nt) of a feature from its first to its last part in transcript order (GenBank order of a
    join): a gene across the origin, and an intron given by its two ends only (NCBI hexacoral nad5 introns hold
    other genes: 10 kb of a 17 kb record) - span() would take the shorter arc the other way round. Parts on both
    strands (trans-splicing) or out of order on a linear record: span()."""
    parts = list(loc.parts)
    if len(parts) == 1:
        return int(parts[0].start), int(parts[0].end) - int(parts[0].start)
    strands = {p.strand for p in parts}
    if len(strands) == 1:
        s, e = (int(parts[0].start), int(parts[-1].end)) if strands != {-1} else \
            (int(parts[-1].start), int(parts[0].end))
        if e > s:
            return s, e - s
        if circular:
            return s, e - s + length
    return span([(int(p.start), int(p.end)) for p in parts], length, circular)


def span(parts: Sequence[Tuple[int, int]], length: int, circular: bool) -> Tuple[int, int]:
    """(start, length in nt) of the shortest arc covering all `parts`; on a circular record it may run across the
    origin (start + length > record length)."""
    ps = sorted(parts)
    start, end = ps[0][0], max(e for _, e in ps)
    if not circular:
        return start, end - start
    best, reach = (start + length - end, 0), ps[0][1]          # widest empty arc: across the origin, or a gap
    for i in range(1, len(ps)):
        best = max(best, (ps[i][0] - reach, i))
        reach = max(reach, ps[i][1])
    gap, i = best
    return (ps[i][0], length - gap) if i else (start, end - start)


def noncoding(rec: SeqRecord, min_len: int = MIN_LEN) -> List[Tuple[int, int]]:
    """(start, end) of the stretches of >= min_len nt no feature covers (a joined gene covers its introns too);
    end > record length for a stretch across the origin of a circular record. min_len 0: none."""
    n = len(rec.seq)
    if min_len <= 0 or not n:
        return []
    circular = rec.annotations.get("topology") == "circular"
    covered = np.zeros(n, dtype=bool)
    for f in rec.features:
        if f.type == "source" or f.location is None or is_ncr(f):     # (NCRs of an earlier run: found again)
            continue
        s, length = extent(f.location, n, circular)
        idx = np.arange(s, s + length) % n
        covered[idx] = True
    if covered.all():
        return []
    runs, start = [], None
    for i in range(n):
        if not covered[i] and start is None:
            start = i
        elif covered[i] and start is not None:
            runs.append([start, i])
            start = None
    if start is not None:
        runs.append([start, n])
    if circular and len(runs) > 1 and runs[0][0] == 0 and runs[-1][1] == n:   # one stretch across the origin
        runs[-1][1] = n + runs.pop(0)[1]
    return [(s, e) for s, e in runs if e - s >= min_len]


def add_ncr(rec: SeqRecord, min_len: int = MIN_LEN) -> int:
    """Replace the record's NCR features with one misc_feature per stretch of >= min_len nt no feature covers (see
    the module notes); returns the number written."""
    n = len(rec.seq)
    rec.features = [f for f in rec.features if not is_ncr(f)]
    regions = noncoding(rec, min_len)
    for k, (s, e) in enumerate(regions, 1):
        loc = FeatureLocation(s, e) if e <= n else CompoundLocation([FeatureLocation(s, n), FeatureLocation(0, e - n)])
        rec.features.append(SeqFeature(loc, type="misc_feature", qualifiers={
            "note": [f"{NOTE}{k}, {e - s:,} bp: no gene annotated (MitoM)"]}))
    # in record order; an NCR across the origin where it begins (its first part), not at position 1
    rec.features.sort(key=lambda f: (f.type != "source", 0 if f.location is None else
                                     int(f.location.parts[0].start) if is_ncr(f) else int(f.location.start)))
    return len(regions)


def main(argv: Sequence[str] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="GenBank file or folder (.gb .gbk .gbf .gbff .genbank)")
    p.add_argument("-o", "--out", required=True, help="output GenBank file (one file in) or folder (a folder in)")
    p.add_argument("--min_len", type=int, default=MIN_LEN, help=f"shortest NCR in nt (default {MIN_LEN})")
    a = p.parse_args(argv)
    src, out = Path(a.input), Path(a.out)
    files = [src] if src.is_file() else sorted(x for x in src.rglob("*") if x.suffix in
                                               (".gb", ".gbk", ".gbf", ".gbff", ".genbank"))
    total = 0
    for f in files:
        recs = list(SeqIO.parse(str(f), "genbank"))
        total += sum(add_ncr(r, a.min_len) for r in recs)
        dest = out if src.is_file() else out / f.relative_to(src)
        dest.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(recs, str(dest), "genbank")
    print(f"[mito_ncr] {total} non-coding regions of >= {a.min_len} nt in {len(files)} file(s) -> {out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
