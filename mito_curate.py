#!/usr/bin/env python3
"""mito_curate.py - rule-based correction of mitochondrial genome annotation (PCGs, 12S/16S rRNA).

Applies the annotation criteria of Fourdrilis et al. 2018 (Sci. Rep. 8:17964) as summarised by
Ghiselli et al. 2021 (Phil. Trans. R. Soc. B 376:20200159, doi:10.1098/rstb.2020.0159, section 3)
to automatic annotation (e.g. MitoZ / MITOS output):

  P1  PCG start: the first eligible in-frame start codon at the 5' end, i.e. the start codon
      nearest the preceding gene WITHOUT overlapping it (criterion i). Only applied when the
      annotated start overlaps the preceding same-strand feature (move downstream), or when an
      in-frame start lies upstream before the preceding feature and the open frame reaches it
      (move upstream, at most --max_start_shift codons; only when the annotated start codon is
      not a valid start, or for --aggressive_start_genes; annelids: nad5 and nad4L). A start that
      overlaps the preceding gene is moved downstream only when a start codon lies within 15 nt.
  P2  Two PCGs on the same strand and in the same frame cannot overlap (criterion ii) - reported.
  D1  A protein gene annotated twice at different places (e.g. MitoZ on a duplicated region): only one
      copy is kept - the complete one (start + stop, no internal stop), the longer, then the first in
      the record; the other CDS and its gene feature are removed (NCBI records annotate one copy).
      An overlapping copy on the other strand that is incomplete or has internal stops is removed as well.
  P3  PCG end: the first in-frame complete stop codon, or an abbreviated stop (T / TA) ending
      immediately before the downstream same-strand tRNA, completed by polyadenylation
      (criterion iii; tRNA punctuation model, Ojala et al. 1981). /transl_except + note added.
      When no tRNA follows and the first complete stop lies > 30 nt inside the next same-strand
      protein gene (not the bicistronic pairs atp8-atp6, nad4L-nad4; rRNA ends are trimmed by R1
      instead; > 10 nt for any other feature), the gene is not
      extended into it: the annotated end is kept when it is an abbreviated T / TA stop (the end of
      the protein homology, e.g. cox1 before cox2), else the T / TA right before the next gene
      (e.g. nad6 before cob) is used. An end that stops short of the first complete stop with
      nothing annotated in between (core homology lost at the C-terminus, 100-500 nt in divergent
      phyla) is extended to that stop (at most 600 nt).
  P7  start and end in different frames (internal stops): replaced by the longer ORF - annotated start
      to its first stop, or annotated stop with the in-frame start nearest the annotated start - when
      it keeps >= half the length (metazoa benchmark 2026-10: 14 CDS with internal stops, NCBI 0).
  O1  circular record: a CDS that ends <= 2 nt before the transcript end without a stop is continued
      across the origin to its first in-frame stop (join; genes cut by the cox1-first rotation); at most 60 nt
      where another feature starts at the origin (genuine overlaps are <= ~50 nt), else 300 nt.
  P5/P6  illegal start / no stop codon: CDS marked 5'- / 3'-partial with 'start / stop codon not
      determined' (NCBI practice; table2asn SEQ_FEAT.StartCodon / NoStop otherwise).
  R1  rRNA (12S rrnS, 16S rrnL) must not overlap flanking genes; boundaries are trimmed to the
      neighbours (criterion vi) and, by default, extended to the flanking genes (gaps <= --max_fill
      nt, the GenBank convention; --rrna_bounds fill). --rrna_bounds trim only removes overlaps.
  R3  rRNA names standardised (rrnS / rrnL; 12S / 16S ribosomal RNA); rRNA with a large CDS
      overlap is reported.
tRNA boundaries are kept as annotated (criterion v). Every change is listed in
curation_report.tsv with the rule, old and new coordinates and the reason.

Each run validates every record with NCBI table2asn before and after curation (when table2asn is
available) and writes validation_compare.tsv: every validator message type per feature marked
fixed / remaining / new, so the effect of the curation is measured, not assumed.

Usage: python3 mito_curate.py -i gbf_file_or_dir -o out_dir --genetic_code 5 [--rrna_bounds trim|fill]
"""
from __future__ import annotations

import argparse
import copy
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, ExactPosition, FeatureLocation, SeqFeature

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402  (finds external programs: config, PATH, conda environments)
import mito_trna  # noqa: E402  (tRNA names: trnL1/L2, trnS1/S2 by anticodon)

GB_EXT = {".gb", ".gbk", ".gbf", ".genbank", ".gbff"}
RRNA_NAMES = {"rrns": "rrnS", "srrna": "rrnS", "12s": "rrnS", "12srrna": "rrnS", "rns": "rrnS", "rrn12": "rrnS",
              "rrnl": "rrnL", "lrrna": "rrnL", "16s": "rrnL", "16srrna": "rrnL", "rnl": "rrnL", "rrn16": "rrnL"}
RRNA_PRODUCT = {"rrnS": "12S ribosomal RNA", "rrnL": "16S ribosomal RNA"}
# N1: NCBI-style protein gene symbols and products (as in RefSeq, e.g. NC_012920.1), with common synonyms
PCG_NAMES = {
    "COX1": ("cytochrome c oxidase subunit I", ("cox1", "coi", "co1", "coxi", "cytochromecoxidasesubuniti")),
    "COX2": ("cytochrome c oxidase subunit II", ("cox2", "coii", "co2", "coxii", "cytochromecoxidasesubunitii")),
    "COX3": ("cytochrome c oxidase subunit III", ("cox3", "coiii", "co3", "coxiii", "cytochromecoxidasesubunitiii")),
    "ATP6": ("ATP synthase F0 subunit 6", ("atp6", "atpase6", "atp", "atpsynthasef0subunit6")),
    "ATP8": ("ATP synthase F0 subunit 8", ("atp8", "atpase8", "atpsynthasef0subunit8")),
    "CYTB": ("cytochrome b", ("cytb", "cob", "ctb", "cyb", "cytochromeb")),
    "ND1": ("NADH dehydrogenase subunit 1", ("nd1", "nad1", "nadh1", "nadhdehydrogenasesubunit1")),
    "ND2": ("NADH dehydrogenase subunit 2", ("nd2", "nad2", "nadh2", "nadhdehydrogenasesubunit2")),
    "ND3": ("NADH dehydrogenase subunit 3", ("nd3", "nad3", "nadh3", "nadhdehydrogenasesubunit3")),
    "ND4": ("NADH dehydrogenase subunit 4", ("nd4", "nad4", "nadh4", "nadhdehydrogenasesubunit4")),
    "ND4L": ("NADH dehydrogenase subunit 4L", ("nd4l", "nad4l", "nadh4l", "nadhdehydrogenasesubunit4l")),
    "ND5": ("NADH dehydrogenase subunit 5", ("nd5", "nad5", "nadh5", "nadhdehydrogenasesubunit5")),
    "ND6": ("NADH dehydrogenase subunit 6", ("nd6", "nad6", "nadh6", "nadhdehydrogenasesubunit6")),
}
PCG_SYN = {syn: name for name, (_, syns) in PCG_NAMES.items() for syn in syns}
SPACERS = (0, 1, 2, 3, -1)  # nt from an abbreviated stop to the downstream tRNA; -1 last: a T on the tRNA's first base is often chance
NO_START_NOTE = "start codon not determined"
NO_STOP_NOTE = "stop codon not determined"
MAX_ORIGIN_EXT = 300        # O1: nt searched past the origin for the stop of a CDS cut at the record end
MAX_INSIDE_TRIM = 250       # R1: an rRNA end may be cut back this far for a gene lying inside it (NC_023929.1: trnV 199 nt
                            # inside MitoZ's l-rRNA overhang); larger losses are logged (R3), not trimmed
MAX_ORIGIN_OVERLAP = 60     # O1: ... when another feature starts there (cox1-first rotation): genuine overlaps are <= ~50 nt
MAX_START_MOVE = 15         # P1: a start overlapping the preceding gene is moved downstream only when a start codon lies
                            # within this distance (173-record metazoa benchmark 2026-10: moves <= 15 nt right 31 / wrong 12,
                            # longer moves right 14 / wrong 28 - NCBI keeps the overlap when the next start is far)
MAX_END_EXTENSION = 600     # P3: longest extension of an end that stopped short of the first complete stop (metazoa
                            # benchmark 2026-10: the 50 confirmed cases were 100-512 nt)
MAX_FEATURE_OVERLAP = 10    # P3: a complete stop may reach this far into the next feature that is not a same-strand
                            # protein gene (173-record truth: into an antisense tRNA 1-10 nt 41 times, deeper 0; into a
                            # same-strand rRNA at most 2 nt)
MAX_PCG_OVERLAP = 30        # P3: a complete stop may reach this far into the next same-strand protein gene (173-record
                            # truth: 11-30 nt inside, complete stop 17 vs T/TA 6; 31-50 nt 6 vs 6; deeper T/TA) (NCBI annelid
                            # records: cox1-cox2, nad6-cob end with T/TA instead; atp8-atp6 / nad4L-nad4 overlap by design)
BICISTRONIC = {("ATP8", "ATP6"), ("ND4L", "ND4")}
POLYA_NOTE = "TAA stop codon is completed by the addition of 3' A residues to the mRNA"


def key(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (name or "").lower())


def span(f: SeqFeature) -> Tuple[int, int]:
    return int(f.location.start) + 1, int(f.location.end)


def spans(f: SeqFeature) -> List[Tuple[int, int]]:
    """1-based (start, end) of every part. span() is the hull of a join: an origin-spanning feature
    join(L-29..L,1..35) is not a feature ending at L (review C finding 6)."""
    return [(int(p.start) + 1, int(p.end)) for p in f.location.parts]


def simple(f: SeqFeature) -> bool:
    """Single-interval feature with exact ends (multi-exon / partial / origin-spanning genes are left alone)."""
    loc = f.location
    return not isinstance(loc, CompoundLocation) and isinstance(loc.start, ExactPosition) \
        and isinstance(loc.end, ExactPosition)


def set_span(f: SeqFeature, s: int, e: int) -> None:
    f.location = FeatureLocation(s - 1, e, strand=f.location.strand)


def gene_feature(rec, f: SeqFeature) -> Optional[SeqFeature]:
    """The 'gene' feature that belongs to f (same name, same strand, overlapping)."""
    name = f.qualifiers.get("gene", [""])[0]
    s, e = span(f)
    for g in rec.features:
        if g.type == "gene" and g.qualifiers.get("gene", [""])[0] == name and g.location.strand == f.location.strand:
            gs, ge = span(g)
            if gs <= e and ge >= s:
                return g
    return None


class Curator:
    def __init__(self, rec, code: int, a):
        self.rec, self.a = rec, a
        self.seq = str(rec.seq).upper()
        self.L = len(self.seq)
        self.code = code
        t = CodonTable.unambiguous_dna_by_id[code]
        self.starts, self.stops = set(t.start_codons), set(t.stop_codons)
        self.rows: List[Dict[str, str]] = []
        self.noted: set = set()
        genes = getattr(a, "aggressive_start_genes", None)
        self.aggressive_genes = {PCG_SYN.get(key(x), x.upper()) for x in genes.split(",") if x} if genes else set()

    # ── helpers in transcript orientation ──
    def log(self, f: SeqFeature, rule: str, old: Tuple[int, int], new: Tuple[int, int], why: str) -> None:
        self.rows.append({"record": self.rec.id, "feature": f.type, "gene": f.qualifiers.get("gene", ["?"])[0],
                          "strand": "+" if f.location.strand == 1 else "-", "rule": rule,
                          "old": f"{old[0]}..{old[1]}", "new": f"{new[0]}..{new[1]}", "reason": why})

    def others(self, f: SeqFeature, types=("tRNA", "rRNA", "CDS"), same_strand: bool = False) -> List[SeqFeature]:
        return [g for g in self.rec.features if g is not f and g.type in types
                and (not same_strand or g.location.strand == f.location.strand)
                and gene_feature(self.rec, g) is not gene_feature(self.rec, f)]

    def frame_differs(self, f: SeqFeature, other: Tuple[int, int]) -> bool:
        """Reading frames of two same-strand CDS: anchored at the genome start on +, at the genome END on -
        (NC_015829.1 nad4L/nad4 on -: 7-nt overlap, starts 297 nt apart = 'same frame' by starts, but different)."""
        s, e = span(f)
        gs, ge = other
        return ((gs - s) if f.location.strand == 1 else (ge - e)) % 3 != 0

    def upstream_end(self, f: SeqFeature) -> int:
        """Transcript-order end of the nearest preceding same-strand feature that the start may not
        overlap: tRNA, rRNA, or a PCG in the same frame. PCGs in a different frame may overlap
        (criterion ii; e.g. atp8-atp6, nad4L-nad4 are translated from bicistronic mRNAs)."""
        s, e = span(f)
        best = 0
        for g in self.others(f, same_strand=True):
            for gs, ge in spans(g):                        # per part: a join is several features here
                if g.type == "CDS" and len(g.location.parts) == 1 and self.frame_differs(f, (gs, ge)):
                    continue                               # (the frame of a join's part is not given by its own start)
                if f.location.strand == 1 and gs < s:
                    best = max(best, ge)
                if f.location.strand == -1 and ge > e:
                    best = max(best, self.L + 1 - gs)      # in reverse coordinates
        return best

    def downstream_gene_start(self, f: SeqFeature) -> Optional[int]:
        """Transcript position where the next same-strand protein gene begins, ignoring the bicistronic
        partner (atp8-atp6, nad4L-nad4) whose overlap is genuine."""
        s, e = span(f)
        me = PCG_SYN.get(key(f.qualifiers.get("gene", [""])[0]), "")
        cand = []
        for g in self.others(f, types=("CDS",), same_strand=True):
            other = PCG_SYN.get(key(g.qualifiers.get("gene", [""])[0]), "")
            if (me, other) in BICISTRONIC:
                continue
            for gs, ge in spans(g):
                if f.location.strand == 1 and gs > s:
                    cand.append(gs)
                if f.location.strand == -1 and ge < e:
                    cand.append(self.L + 1 - ge)
        return min(cand) if cand else None

    def downstream_trna_start(self, f: SeqFeature) -> Optional[int]:
        """Transcript position (reverse coordinates on -) where the next same-strand tRNA begins."""
        s, e = span(f)
        cand = []
        for g in self.others(f, types=("tRNA",), same_strand=True):
            for gs, ge in spans(g):
                if f.location.strand == 1 and gs > s:
                    cand.append(gs)
                if f.location.strand == -1 and ge < e:
                    cand.append(self.L + 1 - ge)
        return min(cand) if cand else None

    def next_feature(self, f: SeqFeature) -> Tuple[Optional[int], Optional[str]]:
        """Transcript position where the nearest downstream feature (tRNA, rRNA, CDS, either strand) begins, and
        its type."""
        s, e = span(f)
        cand = []
        for g in self.others(f):
            for gs, ge in spans(g):
                if f.location.strand == 1 and gs > s:
                    cand.append((gs, g.type))
                if f.location.strand == -1 and ge < e:
                    cand.append((self.L + 1 - ge, g.type))
        return min(cand) if cand else (None, None)

    def next_feature_start(self, f: SeqFeature) -> Optional[int]:
        return self.next_feature(f)[0]

    def tx(self, f: SeqFeature) -> Tuple[int, int]:
        s, e = span(f)
        return (s, e) if f.location.strand == 1 else (self.L + 1 - e, self.L + 1 - s)

    def from_tx(self, f: SeqFeature, a: int, b: int) -> Tuple[int, int]:
        return (a, b) if f.location.strand == 1 else (self.L + 1 - b, self.L + 1 - a)

    def tseq(self, f: SeqFeature) -> str:
        s = self.seq if f.location.strand == 1 else str(self.rec.seq.reverse_complement()).upper()
        return s

    # ── PCG rules ──
    def pcg(self, f: SeqFeature) -> None:
        if not simple(f):
            return
        s = self.tseq(f)
        a, b = self.tx(f)
        old = span(f)
        # P1: start nearest the preceding gene without overlapping it
        up = self.upstream_end(f)
        new_a = a
        if a <= up:                                    # start overlaps preceding feature -> move downstream (if near)
            for c in range(a, min(b - 30, a + MAX_START_MOVE) + 1, 3):
                if c > up and s[c - 1:c + 2] in self.starts:
                    new_a = c
                    break
        elif self.a.aggressive_start or s[a - 1:a + 2] not in self.starts or \
                PCG_SYN.get(key(f.qualifiers.get("gene", [""])[0]), "") in self.aggressive_genes:
            # look upstream, in frame, until stop / preceding gene (criterion i). Without homolog
            # evidence only done when the annotated first codon is not a valid start codon, or for the
            # genes of --aggressive_start_genes (nad5, nad4L in annelids: NCBI annelid records start them at the
            # most upstream in-frame start, 7 / 7 in the annelid benchmark 2026-10; other phyla do not: wrong 52,
            # right 17 on 173 metazoa records, so `mitom run` sets it for the Annelida clade only)
            c = a - 3
            limit = max(up + 1, a - 3 * self.a.max_start_shift)
            while c >= limit and s[c - 1:c + 2] not in self.stops:
                if s[c - 1:c + 2] in self.starts:
                    new_a = c                          # keep the most upstream (nearest the preceding gene)
                c -= 3
        # P3: first complete stop, or T / TA ending right before the downstream tRNA
        first_stop = None
        for c in range(new_a, len(s) - 1, 3):
            if s[c - 1:c + 2] in self.stops:
                first_stop = c + 2
                break
        new_b = b
        trna = self.downstream_trna_start(f)
        down = self.downstream_gene_start(f)
        nfs, nfs_type = self.next_feature(f)           # nearest feature downstream, either strand
        reach = None                                   # furthest position a complete stop may reach (None: no limit)
        if nfs is not None:
            reach = nfs - 1 + (MAX_PCG_OVERLAP if down is not None and nfs == down else MAX_FEATURE_OVERLAP)
        tail_b = (b - new_a + 1) % 3
        has_tail = tail_b in (1, 2) and s[b - tail_b:b] == ("T", "TA")[tail_b - 1]   # annotated end is an abbreviated stop
        why_end = ""
        if down and (trna is None or down < trna) and (nfs is None or nfs >= down) and b < down \
                and first_stop is not None and new_a < down and first_stop - 2 > down - 1 + MAX_PCG_OVERLAP:
            # no feature before the next gene, the annotated end lies before it and the first complete stop
            # lies well inside it: do not extend (an annotated complete stop that already overlaps the next
            # gene is kept: NCBI records show both patterns, 24 complete-stop overlaps > 10 nt and 28 T/TA
            # ends whose complete stop lies > 10 nt inside the next gene, metazoa benchmark 2026-10)
            tail_b = (b - new_a + 1) % 3
            if b <= down - 1 and tail_b in (1, 2) and s[b - tail_b:b] == ("T", "TA")[tail_b - 1]:
                new_b, why_end = b, "abbreviated stop T/TA kept: the first complete stop lies inside the downstream gene"
            else:
                for gap in SPACERS:
                    end = down - 1 - gap
                    tail = (end - new_a + 1) % 3
                    if tail in (1, 2) and s[end - tail:end] == ("T", "TA")[tail - 1]:
                        new_b, why_end = end, "abbreviated stop T/TA immediately before downstream gene" + \
                            (f" ({gap:+d} nt spacer)" if gap else "")
                        break
                    if tail == 0 and s[end - 3:end] in self.stops:
                        new_b, why_end = end, "complete stop codon immediately before downstream gene" + \
                            (f" ({gap:+d} nt spacer)" if gap else "")
                        break
                else:
                    self.log(f, "P4", old, old, "complete stop inside the downstream gene, no abbreviated stop before it - check manually")
        elif trna and new_a <= trna - 1 and (first_stop is None or first_stop - 2 >= trna):
            # the first complete stop lies inside/after the downstream tRNA: the mRNA is cut at the
            # tRNA, so the gene ends with T / TA right before it (criterion iii)
            # NCBI places such ends 0 (143x), -1 (9x) and 1-3 nt (10x) before the tRNA (2026-10-20 dev set,
            # 93 RefSeq records). Order (0, 1, 2, 3, -1) measured best on the dev set: PCG both ends exact
            # 466 -> 468 / 580 (-1 tried second: 467, it matches chance T's on the tRNA's first base)
            for gap in SPACERS:
                end = trna - 1 - gap
                tail = (end - new_a + 1) % 3
                if tail in (1, 2) and s[end - tail:end] == ("T", "TA")[tail - 1]:
                    new_b, why_end = end, "abbreviated stop T/TA immediately before downstream tRNA" + \
                        (f" ({gap:+d} nt spacer)" if gap else "")
                    break
                if tail == 0 and s[end - 3:end] in self.stops:
                    new_b, why_end = end, "complete stop codon immediately before downstream tRNA" + \
                        (f" ({gap:+d} nt spacer)" if gap else "")
                    break
            else:
                self.log(f, "P4", old, old, "no complete or abbreviated stop before the downstream tRNA - check manually")
        elif first_stop is not None and first_stop != b and abs(first_stop - b) <= 3 * self.a.max_start_shift:
            if reach is None or first_stop <= reach:
                # complete stop overlapping the tRNA by 1-2 nt is within tRNA boundary uncertainty: kept
                new_b, why_end = first_stop, "first in-frame complete stop codon"
            elif not has_tail and nfs_type == "rRNA":
                # no stop signal at all and only an rRNA start in the way: that start is soft (R1 fills rRNAs to the
                # flanking genes and trims them back), so the gene is completed and R1 moves the rRNA (review C case)
                new_b, why_end = first_stop, "first in-frame complete stop codon"
            elif not has_tail:
                # the first complete stop lies deeper inside the next feature than NCBI records go (NC_067900.1 nad1:
                # 50 nt into the antisense tRNA): the gene ends with T / TA (or a stop) right before that feature
                for gap in SPACERS:
                    end = nfs - 1 - gap
                    tail = (end - new_a + 1) % 3
                    if end > new_a + 30 and tail in (1, 2) and s[end - tail:end] == ("T", "TA")[tail - 1]:
                        new_b, why_end = end, "abbreviated stop T/TA immediately before the next feature" + \
                            (f" ({gap:+d} nt spacer)" if gap else "")
                        break
                    if end > new_a + 30 and tail == 0 and s[end - 3:end] in self.stops:
                        new_b, why_end = end, "complete stop codon immediately before the next feature" + \
                            (f" ({gap:+d} nt spacer)" if gap else "")
                        break
                else:
                    self.log(f, "P4", old, old, "first complete stop lies inside the next feature, no abbreviated stop "
                             "before it - check manually")
            # else: the annotated end already is an abbreviated stop before that feature: kept
        elif first_stop is not None and first_stop > b and first_stop - b <= MAX_END_EXTENSION and \
                (reach is None or first_stop <= reach):
            # the core's end stops short (C-terminal homology lost: nad5 / nad4 / nad2 / cox1 in divergent phyla,
            # 100-500 nt) and nothing is annotated up to the first complete stop: extend to it
            # (metazoa benchmark 2026-10, 173 RefSeq records: 50 of 52 ends > 100 nt short match NCBI exactly)
            new_b, why_end = first_stop, "first in-frame complete stop codon (annotated end stopped short, no feature in between)"
        if (new_a, new_b) != (a, b) and new_b - new_a + 1 >= 60:
            gs, ge = self.from_tx(f, new_a, new_b)
            set_span(f, gs, ge)
            g = gene_feature(self.rec, f)
            if g is not None and simple(g):
                set_span(g, gs, ge)
            if new_a != a:
                self.log(f, "P1", old, (gs, ge), "start codon nearest the preceding gene without overlap")
            if new_b != b:
                self.log(f, "P3", old, (gs, ge), why_end)
        else:
            new_a, new_b = a, b                        # a rejected move must not drive the stop qualifiers
        # stop qualifiers
        nt = s[new_a - 1:new_b]
        tail = len(nt) % 3
        for k in ("transl_except",):
            f.qualifiers.pop(k, None)
        if tail and nt[-tail:] in ("T", "TA"):
            gs, ge = span(f)
            pos = (f"{ge - tail + 1}..{ge}" if tail > 1 else f"{ge}") if f.location.strand == 1 else \
                (f"complement({gs}..{gs + tail - 1})" if tail > 1 else f"complement({gs})")
            f.qualifiers["transl_except"] = [f"(pos:{pos},aa:TERM)"]
            notes = [n for n in f.qualifiers.get("note", []) if "completed by the addition" not in n]
            f.qualifiers["note"] = notes + [POLYA_NOTE]
        f.qualifiers.pop("translation", None)          # stale after coordinate changes

    def across_origin(self, f: SeqFeature) -> None:
        """O1: CDS cut at the end of a circular record (dev benchmark: ATP8 overlapping COX1 in NC_006902.1 /
        NC_022829.1, ND4L in NC_025335.1 stopped 2 nt before the end) -> joined across the origin to its stop."""
        if self.rec.annotations.get("topology") != "circular" or not simple(f):
            return
        s = self.tseq(f)
        a, b = self.tx(f)
        if b < self.L - 2 or (b - a + 1) % 3 and s[a - 1:b][-((b - a + 1) % 3):] in ("T", "TA"):
            return
        if (b - a + 1) % 3 == 0 and s[b - 3:b] in self.stops:
            return
        first = 1 if f.location.strand == 1 else self.L          # genome position where the continuation starts
        covered = any(gs <= first <= ge for g in self.others(f) for gs, ge in spans(g))
        ss = s + s[:MAX_ORIGIN_OVERLAP if covered else MAX_ORIGIN_EXT]
        end = next((c + 3 for c in range(b - (b - a + 1) % 3, len(ss) - 2, 3) if ss[c:c + 3] in self.stops), None)
        if end is None or end <= self.L:
            return
        old, k = span(f), end - self.L
        if f.location.strand == 1:
            f.location = CompoundLocation([FeatureLocation(a - 1, self.L, 1), FeatureLocation(0, k, 1)])
        else:
            f.location = CompoundLocation([FeatureLocation(0, self.L + 1 - a, -1),
                                           FeatureLocation(self.L - k, self.L, -1)])
        g = gene_feature(self.rec, f)
        if g is not None:
            g.location = f.location
        f.qualifiers.pop("translation", None)
        self.log(f, "O1", old, old, f"continued {k} nt across the origin to the first in-frame stop codon")

    def frame_repair(self, f: SeqFeature) -> None:
        """P7: a CDS whose start and end are in different frames (internal stops) is replaced by the longer of two
        open reading frames: the annotated start to its first in-frame stop, or the annotated end (a complete stop)
        with the in-frame start codon nearest the annotated start. Only when that ORF keeps >= half the length."""
        if not simple(f):
            return
        s = self.tseq(f)
        a, b = self.tx(f)
        nt = s[a - 1:b]
        if len(nt) % 3 == 0:                           # start and end in frame: an internal stop is reported by P4
            return
        if not any(nt[i:i + 3] in self.stops for i in range(0, len(nt) - len(nt) % 3, 3)):
            return
        c = a
        while c + 2 <= len(s) and s[c - 1:c + 2] not in self.stops:
            c += 3
        orfs = [(a, c + 2)] if c + 2 <= len(s) else []
        if s[b - 3:b] in self.stops:
            best, c = None, b - 5
            limit = max(self.upstream_end(f) + 1, a - 3 * self.a.max_start_shift)   # as P1: not into the preceding gene
            while c >= limit and s[c - 1:c + 2] not in self.stops:
                if s[c - 1:c + 2] in self.starts and (best is None or abs(c - a) < abs(best - a)):
                    best = c
                c -= 3
            if best is not None:
                orfs.append((best, b))
        if not orfs:
            return
        na, nb = max(orfs, key=lambda x: x[1] - x[0])
        if nb - na + 1 < (b - a + 1) / 2:
            return
        old = span(f)
        gs, ge = self.from_tx(f, na, nb)
        set_span(f, gs, ge)
        g = gene_feature(self.rec, f)
        if g is not None and simple(g):
            set_span(g, gs, ge)
        self.log(f, "P7", old, (gs, ge), "start and end in different frames (internal stops): "
                 + ("start kept, first in-frame stop" if na == a else "stop kept, in-frame start nearest the old start"))

    def pcg_check(self, f: SeqFeature) -> None:
        """P4: report problems that rules cannot fix (no change is made)."""
        if not simple(f):
            self.log(f, "P4", span(f), span(f), "multi-interval or partial CDS (intron / origin / sequence end) - not curated")
            return
        s = self.tseq(f)
        a, b = self.tx(f)
        nt = s[a - 1:b]
        prot = str(__import__("Bio.Seq", fromlist=["Seq"]).Seq(nt[: len(nt) - len(nt) % 3]).translate(
            table=self.a.genetic_code))
        n_stop = prot.rstrip("*").count("*")
        if n_stop:
            self.log(f, "P4", span(f), span(f), f"{n_stop} internal stop codon(s): wrong frame, frameshift "
                     "(assembly indel?) or unannotated intron - check manually")
        if nt[:3] not in self.starts:
            # P5: NCBI rejects an illegal start (SEQ_FEAT.StartCodon); RefSeq writes such CDS 5'-partial with
            # "start codon not determined" (e.g. echinoderm nad4l ATT/ATC under code 9: NC_001627.1, NC_061414.1)
            loc = f.location
            f.location = FeatureLocation(BeforePosition(loc.start), loc.end, loc.strand) if loc.strand == 1 else \
                FeatureLocation(loc.start, AfterPosition(loc.end), loc.strand)
            f.qualifiers["note"] = [n for n in f.qualifiers.get("note", []) if n != NO_START_NOTE] + [NO_START_NOTE]
            self.log(f, "P5", span(f), span(f), f"first codon {nt[:3]} is not a start codon of table "
                     f"{self.a.genetic_code}: CDS marked 5'-partial, 'start codon not determined' (NCBI practice)")
        tail = len(nt) % 3
        has_stop = nt[-3:] in self.stops if not tail else nt[-tail:] in ("T", "TA")   # last 3 nt are out of frame if tail
        if tail and not has_stop:
            self.log(f, "P4", span(f), span(f), f"length {len(nt)} is not a multiple of 3 and there is no T/TA tail - "
                     "frameshift, wrong end or unannotated intron - check manually")
        if not n_stop and not has_stop:
            # P6: no complete or abbreviated stop -> table2asn SEQ_FEAT.NoStop; written 3'-partial like P5
            loc = f.location
            f.location = FeatureLocation(loc.start, AfterPosition(loc.end), 1) if loc.strand == 1 else \
                FeatureLocation(BeforePosition(loc.start), loc.end, -1)
            f.qualifiers["note"] = [n for n in f.qualifiers.get("note", []) if n != NO_STOP_NOTE] + [NO_STOP_NOTE]
            self.log(f, "P6", span(f), span(f), "no complete or abbreviated stop codon: CDS marked 3'-partial, "
                     "'stop codon not determined' (NCBI practice) - check manually")

    def translation(self, f: SeqFeature) -> None:
        """Review C finding 9: /translation of every CDS (simple or join) recomputed from its final location, never
        dropped (a record is always rewritten). A start codon of the table gives M unless the CDS is 5'-partial
        (NCBI convention). Not written when an internal stop remains (P4 reports it)."""
        f.qualifiers.pop("translation", None)
        nt = str(f.extract(self.rec.seq)).upper()[int(f.qualifiers.get("codon_start", ["1"])[0]) - 1:]
        nt = nt[: len(nt) - len(nt) % 3]
        if not nt:
            return
        prot = str(Seq(nt).translate(table=self.code))
        prot = prot[:-1] if prot.endswith("*") else prot
        p5 = f.location.parts[0]                       # 5'-most part (transcript order)
        partial5 = isinstance(p5.end, AfterPosition) if f.location.strand == -1 else isinstance(p5.start, BeforePosition)
        if nt[:3] in self.starts and not partial5:
            prot = "M" + prot[1:]
        if "*" not in prot:
            f.qualifiers["translation"] = [prot]

    def pcg_overlaps(self) -> None:
        cds = [f for f in self.rec.features if f.type == "CDS" and simple(f)]
        for i, f in enumerate(cds):
            for g in cds[i + 1:]:
                if f.location.strand != g.location.strand:
                    continue
                fs, fe = span(f)
                gs, ge = span(g)
                if fs <= ge and gs <= fe and not self.frame_differs(f, (gs, ge)):
                    self.log(f, "P2", (fs, fe), (gs, ge),
                             f"same-strand same-frame overlap with {g.qualifiers.get('gene', ['?'])[0]} - check manually")

    def duplicates(self) -> None:
        """D1: keep one copy of a protein gene annotated at several non-overlapping places."""
        by: Dict[str, List[SeqFeature]] = {}
        for f in self.rec.features:
            if f.type == "CDS" and simple(f):
                name = PCG_SYN.get(key(f.qualifiers.get("gene", [""])[0]), "")
                if name:
                    by.setdefault(name, []).append(f)
        for name, fs in by.items():
            if len(fs) < 2:
                continue

            def quality(f: SeqFeature) -> tuple:
                nt = self.tseq(f)[self.tx(f)[0] - 1:self.tx(f)[1]]
                tail = len(nt) % 3
                codons = [nt[i:i + 3] for i in range(0, len(nt) - tail, 3)]
                complete = codons and codons[0] in self.starts and \
                    (nt[-tail:] in ("T", "TA") if tail else codons[-1] in self.stops)
                internal = any(c in self.stops for c in (codons if tail else codons[:-1]))
                return (bool(complete) and not internal, len(nt), -int(f.location.start))

            keep = max(fs, key=quality)
            for f in fs:
                if f is keep or any(a <= d and c <= b for a, b in spans(f) for c, d in spans(keep)):
                    continue                            # overlapping copies are the P2 case, not duplicates
                self.log(f, "D1", span(f), span(keep), f"duplicate {name}: second copy removed, {span(keep)[0]}..{span(keep)[1]} kept")
                g = gene_feature(self.rec, f)
                self.rec.features = [x for x in self.rec.features if x is not f and x is not g]
            # an overlapping copy on the other strand, incomplete or with internal stops, next to a complete copy: spurious
            # (MITOS2 / MitoZ hit on the antisense of the real gene: NC_009082.1 nad3, NC_014848.1 nad1)
            if quality(keep)[0]:
                for f in fs:
                    if f is keep or f not in self.rec.features or f.location.strand == keep.location.strand \
                            or quality(f)[0]:
                        continue
                    self.log(f, "D1", span(f), span(keep), f"antisense copy of {name} (incomplete / internal stops) removed, "
                             f"{span(keep)[0]}..{span(keep)[1]} kept")
                    g = gene_feature(self.rec, f)
                    self.rec.features = [x for x in self.rec.features if x is not f and x is not g]

    # ── rRNA rules ──
    def rrna(self, f: SeqFeature) -> None:
        if not simple(f):
            return
        g = gene_feature(self.rec, f)
        s, e = old = span(f)
        nf = [(sp, x) for x in self.others(f, types=("tRNA", "CDS", "rRNA")) for sp in spans(x)]
        neigh = [sp for sp, _ in nf]
        who = {sp: x.qualifiers.get("gene", ["?"])[0] for sp, x in nf}
        noted = self.noted            # R3 'inside' rows once per rRNA pass pair (rrna() runs twice in run())
        changed = True
        while changed:                                # R1: no overlap with any other gene (repeat until stable)
            changed = False
            for ns, ne in sorted(neigh):
                if ne < s or ns > e:
                    continue
                if ns <= s and ne >= e:               # rRNA inside another gene: left to R3 / manual check
                    continue
                if ns <= s:                           # overlaps the 5' (left) end
                    s = ne + 1
                elif ne >= e:                         # overlaps the 3' (right) end
                    e = ns - 1
                else:                                 # gene inside the rRNA: keep the longer side, but only if what is
                    keep_left = ns - s >= e - ne      # lost is <= 10 % or <= MAX_INSIDE_TRIM nt (a tRNA-sized overhang
                    lost = e - ns + 1 if keep_left else ne - s + 1   # of the core; a spurious tRNA call deep inside
                    if lost * 10 > e - s + 1 and lost > MAX_INSIDE_TRIM:   # must not cut a core rRNA)
                        if (ns, ne) not in noted:
                            noted.add((ns, ne))
                            self.log(f, "R3", old, old, f"gene {who[(ns, ne)]} inside rRNA, not trimmed (would lose "
                                     f"{lost} nt) - check manually (spurious tRNA call?)")
                        continue
                    if keep_left:
                        e = ns - 1
                    else:
                        s = ne + 1
                changed = True
        if self.a.rrna_bounds == "fill":
            left = max([ne for ns, ne in neigh if ne < s] or [0])
            right = min([ns for ns, ne in neigh if ns > e] or [self.L + 1])
            if 0 < s - left - 1 <= self.a.max_fill:
                s = left + 1
            if 0 < right - e - 1 <= self.a.max_fill:
                e = right - 1
        cds_ov = max([min(e, ne) - max(s, ns) + 1 for x in self.others(f, types=("CDS",)) for ns, ne in spans(x)] or [0])
        if cds_ov > 50:
            self.log(f, "R3", old, (s, e), f"rRNA overlaps a CDS by {cds_ov} nt - boundary or assembly problem?")
        if (s, e) != old and e - s + 1 >= 200:
            set_span(f, s, e)
            if g is not None and simple(g):
                set_span(g, s, e)
            self.log(f, "R1", old, (s, e), "trimmed to flanking genes" + (" and filled to them" if self.a.rrna_bounds == "fill" else ""))

    # ── N1: NCBI-style names (every feature, also joins / origin-spanning ones) ──
    def names(self, f: SeqFeature) -> None:
        g = gene_feature(self.rec, f)
        cur = f.qualifiers.get("gene", [""])[0]
        q: Dict[str, List[str]] = {}
        if f.type == "CDS":
            name = PCG_SYN.get(key(cur)) or PCG_SYN.get(key(f.qualifiers.get("product", [""])[0]))
            if name:
                q = {"gene": [name], "product": [PCG_NAMES[name][0]]}
        elif f.type == "rRNA":
            name = RRNA_NAMES.get(key(cur)) or \
                {"12sribosomalrna": "rrnS", "16sribosomalrna": "rrnL"}.get(key(f.qualifiers.get("product", [""])[0]))
            if name:
                q = {"gene": [name], "product": [RRNA_PRODUCT[name]]}
        elif f.type == "tRNA":
            txt = " ".join(f.qualifiers.get("product", []) + [cur])
            m = re.search(r"tRNA-([A-Z][a-z]{2})", txt)
            aa = mito_trna.AA3.get(m.group(1), "") if m else ""
            if not aa and re.match(r"trn([A-Z])", cur):
                aa = cur[3]
            ac = ""
            m = re.search(r"\(([acgtuACGTU]{3})\)", cur) or \
                re.search(r"anticodon:\s*([acgtuACGTU]{3})", " ".join(f.qualifiers.get("note", []))) or \
                re.search(r"seq:([acgtuACGTU]{3})", " ".join(f.qualifiers.get("anticodon", [])))
            if m:
                ac = m.group(1).lower().replace("u", "t")
            if aa in mito_trna.ONE2THREE:
                name = mito_trna.trna_name(aa, ac)
                notes = [n for n in f.qualifiers.get("note", []) if not n.startswith("anticodon:")]
                q = {"gene": [name], "product": [f"tRNA-{mito_trna.ONE2THREE[aa]}"],
                     "note": ([f"anticodon: {ac.upper().replace('T', 'U')}"] if ac else []) + notes}
        if not q:
            if f.type in ("CDS", "rRNA", "tRNA"):
                self.log(f, "N1", span(f), span(f), f"gene name '{cur}' not recognised - kept")
            return
        for k, v in q.items():
            if v:
                f.qualifiers[k] = v
            else:
                f.qualifiers.pop(k, None)
        if g is not None:
            g.qualifiers["gene"] = q["gene"]

    def run(self) -> List[Dict[str, str]]:
        for f in list(self.rec.features):
            if f.type in ("CDS", "rRNA", "tRNA"):
                self.names(f)
        for f in list(self.rec.features):              # R1 first: a raw rRNA overhang must not push a CDS start (P1)
            if f.type == "rRNA":
                self.rrna(f)
        self.duplicates()                              # D1 before the boundary rules: a copy must not act as neighbour
        for f in list(self.rec.features):
            if f.type == "CDS":
                self.frame_repair(f)
                self.pcg(f)
                self.across_origin(f)
                self.pcg_check(f)
        self.pcg_overlaps()
        for f in list(self.rec.features):              # R1 again: P3 / O1 moved CDS ends, trim the rRNAs to them
            if f.type == "rRNA":
                self.rrna(f)
        for f in self.rec.features:
            if f.type == "CDS":
                self.translation(f)
        for f in self.rec.features:                    # gene feature = exactly its CDS / RNA (SEQ_FEAT.CDSgeneRange)
            if f.type in ("CDS", "rRNA", "tRNA"):
                g = gene_feature(self.rec, f)
                if g is not None and str(g.location) != str(f.location):
                    g.location = f.location
        self.rec.features.sort(key=lambda x: (self.sort_start(x), x.type != "source", x.type != "gene"))
        return self.rows

    def sort_start(self, f: SeqFeature) -> int:
        """0-based sort position; an origin-spanning feature (a part from 0, a part to L) sorts by the start of the
        part that ends at L, not at 0 (review C finding 20)."""
        parts = f.location.parts
        if len(parts) > 1 and any(int(p.start) == 0 for p in parts):
            tails = [int(p.start) for p in parts if int(p.end) == self.L and int(p.start) > 0]
            if tails:
                return min(tails)
        return min(int(p.start) for p in parts)


FEATURE_RE = re.compile(r"FEATURE:\s*([\w-]+):\s*([^<\[]+)")


def msg_key(r: Dict[str, str]) -> Tuple[str, str, str, str]:
    """Match validator messages before/after curation by record, code and feature (type + name);
    coordinates are left out because curation changes them."""
    m = FEATURE_RE.search(r["message"])
    feat = f"{m.group(1)}: {m.group(2).strip()}" if m else ""
    return r["record"], r["severity"], r["code"], feat


def compare_validation(before: List[Dict[str, str]], after: List[Dict[str, str]]) -> List[Dict[str, str]]:
    from collections import Counter
    b, a = Counter(msg_key(r) for r in before), Counter(msg_key(r) for r in after)
    rows = []
    for k in sorted(set(b) | set(a)):
        nb, na = b[k], a[k]
        status = "fixed" if na < nb and na == 0 else "remaining" if na and nb else "new" if na > nb else "fixed"
        rows.append({"record": k[0], "severity": k[1], "code": k[2], "feature": k[3],
                     "before": str(nb), "after": str(na), "status": status})
    return rows


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--genetic_code", type=int, default=5)
    p.add_argument("--rrna_bounds", choices=["trim", "fill"], default="fill",
                   help="fill (default): remove overlaps with flanking genes and extend rRNA ends to them (gaps <= "
                        "--max_fill), the GenBank convention; trim: only remove overlaps (Fourdrilis criterion vi). "
                        "2026-10-21 dev benchmark (90 rRNAs vs NCBI): both ends exact trim 22, fill(300) 33")
    p.add_argument("--max_fill", type=int, default=300,
                   help="Largest gap filled with --rrna_bounds fill (nt; default 300: exact rRNAs 30 / 33 / 33 / 31 "
                        "at 100 / 300 / 500 / 1000 on the dev set)")
    p.add_argument("--max_start_shift", type=int, default=30, help="Largest start/stop move (codons)")
    p.add_argument("--aggressive_start", action="store_true",
                   help="Always move PCG starts to the in-frame start codon nearest the preceding gene "
                        "(criterion i literally); default only when the annotated start codon is invalid")
    p.add_argument("--aggressive_start_genes", default="",
                   help="Genes whose start is always moved to the in-frame start codon nearest the preceding gene "
                        "(default none; `mitom run` passes ND5,ND4L for the Annelida clade). Benchmark 2026-10: right "
                        "for nad5 / nad4L in every differing annelid record (40), wrong in most other phyla (173 "
                        "metazoa records: 17 right, 52 wrong).")
    p.add_argument("--table2asn", default=None,
                   help="NCBI table2asn for validation before/after curation (default: found by mito_env.py)")
    p.add_argument("--no_validate", action="store_true", help="Skip the NCBI validation before/after curation")
    a = p.parse_args(argv)
    src, out = Path(a.input), Path(a.out)
    t2a = None if a.no_validate else mito_env.find("table2asn", a.table2asn)
    if t2a:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import ncbi_validate as nv
    val_before: List[Dict[str, str]] = []
    val_after: List[Dict[str, str]] = []
    files = sorted(x for x in src.rglob("*") if x.suffix.lower() in GB_EXT) if src.is_dir() else [src]
    out.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, str]] = []
    for f in files:
        recs = list(SeqIO.parse(str(f), "genbank"))
        for i, rec in enumerate(recs):
            if t2a:
                val_before += nv.validate_record(copy.deepcopy(rec), out / "validation" / "before", t2a)
            backup = copy.deepcopy(rec)                  # run() edits in place: a failure must not leave a half-curated record
            try:
                rows += Curator(rec, a.genetic_code, a).run()
            except Exception as exc:
                recs[i] = backup
                rows.append({"record": rec.id, "rule": "error", "reason": f"left unchanged: {exc!r}"})
        if t2a:
            for rec in recs:
                val_after += nv.validate_record(rec, out / "validation" / "after", t2a)
        dest = out / (f.relative_to(src) if src.is_dir() else f.name)
        dest.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(recs, str(dest), "genbank")      # always: N1 names / gene ranges change without a report row
    cols = ["record", "feature", "gene", "strand", "rule", "old", "new", "reason"]
    with open(out / "curation_report.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    c = Counter(r["rule"] for r in rows)
    print(f"[mito_curate] {len(files)} file(s); changes: " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())))
    if t2a:
        cmp_rows = compare_validation(val_before, val_after)
        vcols = ["record", "severity", "code", "feature", "before", "after", "status"]
        with open(out / "validation_compare.tsv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=vcols, delimiter="\t")
            w.writeheader()
            w.writerows(cmp_rows)
        nb = Counter(r["severity"] for r in val_before)
        na = Counter(r["severity"] for r in val_after)
        st = Counter((r["severity"], r["status"]) for r in cmp_rows)
        print(f"[mito_curate] NCBI table2asn before -> after: ERROR {nb['ERROR']} -> {na['ERROR']}, "
              f"WARNING {nb['WARNING']} -> {na['WARNING']}; ERROR types fixed {st[('ERROR', 'fixed')]}, "
              f"remaining {st[('ERROR', 'remaining')]}, new {st[('ERROR', 'new')]} -> {out / 'validation_compare.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
