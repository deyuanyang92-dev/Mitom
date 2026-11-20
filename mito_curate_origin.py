"""Rules of mito_curate.Curator for genes across the origin of a circular record: O1 / O2 and the checks of joined CDS.

Code moved verbatim from mito_curate.py; the methods use the Curator helpers (self.log, self.tx, ...)."""
from __future__ import annotations

from typing import Optional, Tuple

from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    FeatureLocation,
    SeqFeature,
)

from mito_common import POLYA_NOTE, codon_translate, ends_with_stop, transl_except_term
from mito_curate_base import (
    MAX_FEATURE_OVERLAP,
    MAX_ORIGIN_EXT,
    MAX_ORIGIN_OVERLAP,
    NO_START_NOTE,
    NO_STOP_NOTE,
    PSEUDO_NOTE,
    apply_excepts,
    gene_feature,
    internal_excepts,
    set_span,
    simple,
    span,
    spans,
)


class OriginRules:
    """Mixin of mito_curate.Curator."""

    def across_origin(self, f: SeqFeature) -> None:
        """O1: CDS cut at the end of a circular record (dev benchmark: ATP8 overlapping COX1 in NC_006902.1 /
        NC_022829.1, ND4L in NC_025335.1 stopped 2 nt before the end) -> joined across the origin to its stop."""
        if self.rec.annotations.get("topology") != "circular" or not simple(f):
            return
        s = self.tseq(f)
        a, b = self.tx(f)
        tail = (b - a + 1) % 3
        if tail and s[b - tail:b] in ("T", "TA") or not tail and s[b - 3:b] in self.stops:
            return
        if b < self.L - 2:
            # the native core ends at the last aligned codon, which may lie further from the sequence end (review record
            # F: 12 nt, stop 6 nt past the origin): continued when the open frame reaches the end with nothing annotated
            # in between (an in-frame stop before the end is P3's business; a feature there bounds the gene)
            nfs = self.next_feature_start(f)
            if nfs is not None and nfs > b or any(s[c:c + 3] in self.stops for c in range(b - tail, self.L - 2, 3)):
                return
        first = 1 if f.location.strand == 1 else self.L          # genome position where the continuation starts
        covered = any(gs <= first <= ge for g in self.others(f) for gs, ge in spans(g))
        ss = s + s[:MAX_ORIGIN_OVERLAP if covered else MAX_ORIGIN_EXT]
        end = next((c + 3 for c in range(b - (b - a + 1) % 3, len(ss) - 2, 3) if ss[c:c + 3] in self.stops), None)
        if end is not None and end <= self.L or end is None and not covered:
            return
        if end is None or covered and end - self.L > MAX_FEATURE_OVERLAP or end - self.L <= 2:
            # (no stop within MAX_ORIGIN_OVERLAP nt past the origin, a gene starting there: annelid nad2 before cox1 in
            # the cox1-first frame, KT429007.1 / NC_012618.1 / EF494507.1 - NCBI ends it with T at the sequence end)
            # (a stop completed by only 1-2 nt past the origin: NCBI ends such genes with T / TA at the sequence end,
            # 11 of 12 cases in 255 RefSeq records; NC_023772.1 cob)
            # the stop lies well inside the feature at the origin: an in-frame T / TA ending at the sequence end is the
            # abbreviated stop NCBI annotates (NC_065213.1, NC_065650.1 nad2 before cox1: 23 / 29 nt continuations)
            for tail in (2, 1):
                if (self.L - tail - a + 1) % 3 == 0 and s[self.L - tail:self.L] == ("T", "TA")[tail - 1]:
                    old = span(f)
                    gs, ge = self.from_tx(f, a, self.L)
                    set_span(f, gs, ge)
                    g = gene_feature(self.rec, f)
                    if g is not None and simple(g):
                        set_span(g, gs, ge)
                    f.qualifiers.pop("translation", None)
                    f.qualifiers["transl_except"] = internal_excepts(f) + [transl_except_term(
                        gs, ge, 1 if f.location.strand == 1 else -1, tail)]
                    f.qualifiers["note"] = [n for n in f.qualifiers.get("note", []) if "completed by the addition" not in n] \
                        + [POLYA_NOTE]
                    self.log(f, "O1", old, (gs, ge), f"abbreviated stop {('T', 'TA')[tail - 1]} at the sequence end "
                             + (f"(the first complete stop lies {end - self.L} nt inside the gene at the origin)" if end
                                else f"(no stop within {MAX_ORIGIN_OVERLAP} nt past the origin, a gene starts there)"))
                    return
            if end is None:
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

    def origin_cut(self, f: SeqFeature) -> None:
        """O2: a CDS the core wrote as a join across the origin of a circular record, with a short part past the origin
        (the gene before cox1 in MitoM's cox1-first frame: annelid nad2 running 8-23 nt into cox1 without a stop codon,
        table2asn NoStop in 14 records of the 2026-10 runs), is cut back to its part before the origin; P1 / P3 / O1 then
        place its end like that of any gene cut at the sequence end (T / TA at the end, or the first stop past it)."""
        origin = self.origin_span(f)
        if origin is None:
            return
        head, tail = f.location.parts                  # transcript order: the part before the origin first
        if len(tail) > MAX_ORIGIN_EXT or len(head) < 60:
            return
        nt = str(f.location.extract(self.rec.seq)).upper()
        if ends_with_stop(nt, self.stops) and "*" not in apply_excepts(f, codon_translate(nt[:len(nt) - len(nt) % 3], self.code))[:-1]:
            # a complete gene across the origin (NCBI / MITOS2 records, a rotated input): not the core's run-on into
            # cox1 (2026-10-31 review: ~25 real genes were cut in half, some given an invented T / TA stop)
            return
        f.location = FeatureLocation(head.start, head.end, head.strand)
        g = gene_feature(self.rec, f)
        if g is not None:
            g.location = f.location
        f.qualifiers.pop("translation", None)
        self.log(f, "O2", origin, span(f), f"the core's join across the origin cut back to the sequence end ({len(tail)} nt "
                 "past it); the end is placed again by P3 / O1")

    def origin_join_check(self, f: SeqFeature) -> None:
        """P4 / P5 / P6 / P10 on a CDS joined across the origin that O1 / O2 did not end (the core's join whose part
        before the origin is short, so O2 cannot cut it back: leave-out runs, NC_006386.1 nad4L 38 nt before the
        origin on -, NC_023783.1 cox3; table2asn NoStop): the joined sequence is checked like a single interval -
        internal stops /pseudo, no start codon 5'-partial, no stop 3'-partial with the NCBI notes."""
        nt = str(f.extract(self.rec.seq)).upper()
        prot = apply_excepts(f, codon_translate(nt, self.code))
        parts = list(f.location.parts)
        minus = f.location.strand == -1
        if prot.rstrip("*").count("*"):
            for x in (f, gene_feature(self.rec, f)):
                if x is not None:
                    x.qualifiers["pseudo"] = [None]
            f.qualifiers["note"] = [n for n in f.qualifiers.get("note", []) if n != PSEUDO_NOTE] + [PSEUDO_NOTE]
            self.log(f, "P10", span(f), span(f), "CDS across the origin with internal stop codons written /pseudo")
            return
        has_stop = ends_with_stop(nt, self.stops)
        if nt[:3] not in self.starts:
            p = parts[0]
            parts[0] = FeatureLocation(p.start, AfterPosition(p.end), -1) if minus else \
                FeatureLocation(BeforePosition(p.start), p.end, 1)
            f.qualifiers["note"] = [n for n in f.qualifiers.get("note", []) if n != NO_START_NOTE] + [NO_START_NOTE]
            self.log(f, "P5", span(f), span(f), f"CDS across the origin: first codon {nt[:3]} is not a start codon: "
                     "5'-partial, 'start codon not determined'")
        if not has_stop:
            p = parts[-1]
            parts[-1] = FeatureLocation(BeforePosition(p.start), p.end, -1) if minus else \
                FeatureLocation(p.start, AfterPosition(p.end), 1)
            f.qualifiers["note"] = [n for n in f.qualifiers.get("note", []) if n != NO_STOP_NOTE] + [NO_STOP_NOTE]
            self.log(f, "P6", span(f), span(f), "CDS across the origin without a stop codon: 3'-partial, 'stop codon "
                     "not determined' - check manually")
        f.location = CompoundLocation(parts)

    # ── rRNA rules ──
    def origin_span(self, f: SeqFeature) -> Optional[Tuple[int, int]]:
        """(start, end) of a two-part join across the origin of a circular record as one interval, end > L."""
        parts = f.location.parts
        if self.rec.annotations.get("topology") != "circular" or len(parts) != 2:
            return None
        tail = next((p for p in parts if int(p.end) == self.L), None)
        head = next((p for p in parts if int(p.start) == 0), None)
        if tail is None or head is None or tail is head:
            return None
        return int(tail.start) + 1, self.L + int(head.end)

    def set_circular_span(self, f: SeqFeature, s: int, e: int) -> None:
        """set_span on the circle: s..e (1-based, e may exceed L) as one interval or as a join across the origin."""
        if s < 1:                                       # filled back across the origin
            s, e = s + self.L, e + self.L
        if e > self.L and s > self.L:
            s, e = s - self.L, e - self.L
        if e <= self.L:
            set_span(f, s, e)
            return
        strand = f.location.strand
        tail, head = FeatureLocation(s - 1, self.L, strand), FeatureLocation(0, e - self.L, strand)
        f.location = CompoundLocation([tail, head] if strand != -1 else [head, tail])
