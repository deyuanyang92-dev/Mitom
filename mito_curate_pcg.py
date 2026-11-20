"""Protein-gene rules of mito_curate.Curator: P1 / P3 / P11 starts and ends, P7 / P8 frame repair, P9 frameshift, P4 / P10 checks, D1 duplicates, /translation.

Code moved verbatim from mito_curate.py; the methods use the Curator helpers (self.log, self.tx, ...)."""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from Bio.Seq import Seq
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    ExactPosition,
    FeatureLocation,
    SeqFeature,
)

import mito_start
from mito_common import (
    POLYA_NOTE,
    abbreviated_stop,
    codon_translate,
    ends_with_stop,
    transl_except_term,
)
from mito_curate_base import (
    FRAMESHIFT_MAX_TRIM,
    FRAMESHIFT_NOTE,
    FRAMESHIFT_WINDOW,
    MAX_END_EXTENSION,
    MAX_FEATURE_OVERLAP,
    MAX_PCG_OVERLAP,
    MAX_START_MOVE,
    ND3_EXTRA,
    ND3_NOTE,
    NO_START_NOTE,
    NO_STOP_NOTE,
    OVERLAP_FALLBACK,
    PSEUDO_NOTE,
    SPACERS,
    STOP_TRIM_MIN,
    apply_excepts,
    gene_feature,
    internal_excepts,
    set_span,
    simple,
    span,
    spans,
)
from mito_names import PCG_SYN
from mito_names import name_key as key


class PcgRules:
    """Mixin of mito_curate.Curator."""

    def learned_start(self, f: SeqFeature, s: str, b: int, a: int) -> int:
        """P11: of the in-frame start codons near P1's start a (reading frame open in between), the one the start
        model scores best (mito_start.py); with --dump_starts every candidate is also kept for training."""
        gene = PCG_SYN.get(key(f.qualifiers.get("gene", [""])[0]), "")
        if not gene or id(f) not in self.aln:          # only genes of MitoM's core (its alignment note): the model knows
            return a                                   # no other annotator's starts (MitoZ / MITOS2 / mitom import)
        up = None                                      # the preceding feature, either strand: (end, type, same strand)
        for g in self.others(f):
            for gs, ge in spans(g):
                ts, te = (gs, ge) if f.location.strand == 1 else (self.L + 1 - ge, self.L + 1 - gs)
                if ts < a and te <= a + 60 and (up is None or te > up[0]):
                    up = (te, g.type, g.location.strand == f.location.strand)
        rows = mito_start.candidates(s, a, b, self.starts, self.stops, self.core_start.get(id(f)), self.aln.get(id(f)),
                                     up, gene, self.phylum, self.code)
        if getattr(self.a, "dump_starts", None):
            for r in rows:
                self.start_rows.append({"record": self.rec.id, "strand": "+" if f.location.strand == 1 else "-",
                                        "genome_pos": r["pos"] if f.location.strand == 1 else self.L + 1 - r["pos"], **r})
        if self.start_model is None or len(rows) < 2:
            return a
        return rows[self.start_model.best(rows)[0]]["pos"]

    def pcg(self, f: SeqFeature) -> None:
        if not simple(f):
            return
        s = self.tseq(f)
        a, b = self.tx(f)
        old = span(f)
        # P1: start nearest the preceding gene without overlapping it
        up = self.upstream_end(f)
        new_a = a
        why_start = "start codon nearest the preceding gene without overlap"
        enc = [te for ts, te in self.enclosed_trnas(f, a, b) if ts - a < (b - a + 1) / 3]
        if enc:
            # a whole tRNA in the 5' third of the gene: the start was placed upstream across it (unionid nad1 across trnL2
            # into trnV, NC_085464.1; table2asn SEQ_FEAT.BadCDScomponentOverlapTRNA). The first in-frame start codon after
            # it within 20 codons (NCBI: 0-2 nt after the tRNA in 13 of 14 such records), else the first in-frame codon
            # after it (P5 then marks the gene 'start codon not determined')
            c0 = a + -(-(max(enc) + 1 - a) // 3) * 3
            new_a = next((c for c in range(c0, min(b - 30, c0 + 60) + 1, 3) if s[c - 1:c + 2] in self.starts), c0)
            why_start = "start moved past the whole tRNA it contained"
        elif a <= up:                                  # start overlaps preceding feature -> move downstream (if near)
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
        elif s[a - 1:a + 2] != "ATG" and not getattr(self.a, "no_upstream_atg", False):
            # a start codon other than ATG (homology alone cannot place a divergent N-terminus): the most upstream
            # in-frame ATG before the preceding gene, open frame in between (leave-out core, 255 records: +37 exact
            # starts; not for nematodes / arthropods, whose genes often start with TTG / ATT / ATA: 6 right, 14 wrong)
            c, limit = a - 3, max(up + 1, a - 3 * self.a.max_start_shift)
            while c >= limit and s[c - 1:c + 2] not in self.stops:
                if s[c - 1:c + 2] == "ATG":
                    new_a = c
                c -= 3
        start_rule = "P1"
        if self.start_model is not None or getattr(self.a, "dump_starts", None):
            chosen = self.learned_start(f, s, b, new_a)
            if chosen != new_a:
                new_a, start_rule = chosen, "P11"
                why_start = "start codon where NCBI records put it (start model, mito_start.py)"
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
        pcg_overlap = MAX_FEATURE_OVERLAP if getattr(self.a, "abbreviated_pcg_overlap", False) else MAX_PCG_OVERLAP
        if nfs is not None:
            reach = nfs - 1 + (pcg_overlap if down is not None and nfs == down else MAX_FEATURE_OVERLAP)
        tail_b = (b - new_a + 1) % 3
        has_tail = tail_b in (1, 2) and s[b - tail_b:b] == ("T", "TA")[tail_b - 1]   # annotated end is an abbreviated stop
        reach_hard = reach                             # the same bound ignoring an rRNA start (branch 4, as in branch 3): R1
        if nfs_type == "rRNA" and not has_tail:        # fills the rRNA to the short CDS end first (review record B: end
            nfs_h = self.next_feature(f, types=("tRNA", "CDS"))[0]   # 150 nt short, reach = rRNA start + 10 blocked the
            reach_hard = None if nfs_h is None else \
                nfs_h - 1 + (pcg_overlap if nfs_h == down else MAX_FEATURE_OVERLAP)   # extension) and trims it back after
        why_end = ""
        if down and (trna is None or down < trna) and (nfs is None or nfs >= down) and b < down \
                and first_stop is not None and new_a < down and first_stop - 2 > down - 1 + pcg_overlap:
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
                    # (annelids, --abbreviated_pcg_overlap: no T / TA before the gene - the complete stop as elsewhere)
                    if not self.overlap_fallback(first_stop, down) and first_stop - 2 > down - 1 + MAX_PCG_OVERLAP:
                        self.log(f, "P4", old, old, "complete stop inside the downstream gene, no abbreviated stop before it - check manually")
                    else:
                        new_b, why_end = first_stop, f"first complete stop, {first_stop - down + 1} nt into the downstream gene (no abbreviated stop before it)"
        elif trna and new_a <= trna - 1 and (first_stop is None or first_stop - 2 >= trna) \
                and not (first_stop == b and first_stop - trna + 1 <= OVERLAP_FALLBACK):
            # (a complete stop the core's alignment reached, running <= OVERLAP_FALLBACK nt into the predicted tRNA, is
            # kept: NC_023772.1 nad4 TAG 10 nt into tRNA-His; leave-out, 255 records: +10 exact ends)
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
                if not self.overlap_fallback(first_stop, trna) and any(ts == trna and b - te < (b - new_a + 1) / 3
                                                                       for ts, te in self.enclosed_trnas(f, new_a, b)):
                    # the gene would contain the whole tRNA (table2asn SEQ_FEAT.BadCDScomponentOverlapTRNA): the mRNA is cut
                    # at the tRNA, the gene ends right before it without a stop codon (P6: 3'-partial)
                    new_b, why_end = trna - 1, "ends before the whole tRNA it contained (no stop codon before it)"
                elif not self.overlap_fallback(first_stop, trna):
                    self.log(f, "P4", old, old, "no complete or abbreviated stop before the downstream tRNA - check manually")
                else:
                    new_b, why_end = first_stop, f"first complete stop, {first_stop - trna + 1} nt into the downstream tRNA (no abbreviated stop before it)"
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
                    if not self.overlap_fallback(first_stop, nfs):
                        self.log(f, "P4", old, old, "first complete stop lies inside the next feature, no abbreviated stop "
                                 "before it - check manually")
                    else:
                        new_b, why_end = first_stop, f"first complete stop, {first_stop - nfs + 1} nt into the next feature (no abbreviated stop before it)"
            # else: the annotated end already is an abbreviated stop before that feature: kept
        elif first_stop is not None and first_stop > b and first_stop - b <= MAX_END_EXTENSION and \
                (reach_hard is None or first_stop <= reach_hard):
            # the core's end stops short (C-terminal homology lost: nad5 / nad4 / nad2 / cox1 in divergent phyla,
            # 100-500 nt) and nothing is annotated up to the first complete stop: extend to it
            # (metazoa benchmark 2026-10, 173 RefSeq records: 50 of 52 ends > 100 nt short match NCBI exactly)
            new_b, why_end = first_stop, "first in-frame complete stop codon (annotated end stopped short, no feature in between)"
        nxt = trna                                     # (a tRNA only: before a protein gene NCBI keeps the complete stop)
        if getattr(self.a, "abbreviated_overlap", False) and nxt is not None and 1 <= new_b - nxt + 1 <= 2 \
                and (new_b - new_a + 1) % 3 == 0 and s[new_b - 3:new_b] in self.stops:
            # vertebrate RefSeq convention: a complete stop running 1-2 nt into the next same-strand tRNA is
            # annotated as the abbreviated stop T / TA ending right before it (leave-out core, 255 records: Chordata
            # +14 exact ends, 1 lost; other phyla keep the complete stop: 82 records, 63 lost). Not before a protein
            # gene: atp6's TAA sharing its A with cox3's ATG is kept (1,719 leave-out records 2026-10-30: 294 more
            # exact ends, 70 fewer, mostly ray-finned fishes, which write it both ways: 104 / 60)
            tail = (nxt - 1 - new_a + 1) % 3
            if tail in (1, 2) and s[nxt - 1 - tail:nxt - 1] == ("T", "TA")[tail - 1]:
                new_b, why_end = nxt - 1, "abbreviated stop T/TA immediately before the next tRNA " \
                    "(the complete stop would overlap it)"
        if (new_a, new_b) != (a, b) and new_b - new_a + 1 >= 60:
            gs, ge = self.from_tx(f, new_a, new_b)
            set_span(f, gs, ge)
            g = gene_feature(self.rec, f)
            if g is not None and simple(g):
                set_span(g, gs, ge)
            if new_a != a:
                self.log(f, start_rule, old, (gs, ge), why_start)
            if new_b != b:
                self.log(f, "P3", old, (gs, ge), why_end)
        else:
            new_a, new_b = a, b                        # a rejected move must not drive the stop qualifiers
        # stop qualifiers
        nt = s[new_a - 1:new_b]
        tail = len(nt) % 3
        internal_excepts(f)                             # the old TERM exception goes; codon exceptions inside stay
        loc = f.location                               # an end left 3'-partial (cut by a linear sequence end) may continue
        cut = isinstance(loc.end, AfterPosition) if loc.strand == 1 else isinstance(loc.start, BeforePosition)
        if abbreviated_stop(nt) and not cut:           # (table2asn PartialProblemHasStop with a T / TA stop)
            gs, ge = span(f)
            f.qualifiers["transl_except"] = internal_excepts(f) + [transl_except_term(gs, ge, 1 if f.location.strand == 1
                                                                                      else -1, tail)]
            notes = [n for n in f.qualifiers.get("note", []) if "completed by the addition" not in n]
            f.qualifiers["note"] = notes + [POLYA_NOTE]
        f.qualifiers.pop("translation", None)          # stale after coordinate changes

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

    def internal_stops(self, s: str, a: int, b: int) -> List[int]:
        """Transcript positions of the in-frame stop codons of a..b before its last codon (frame of a)."""
        last = b - 2 if (b - a + 1) % 3 == 0 else b + 1    # the last complete codon may be the stop
        return [c for c in range(a, last - 2, 3) if s[c - 1:c + 2] in self.stops]

    def single_frame_repairs(self, s: str, a: int, b: int) -> List[Tuple[int, int]]:
        """Open reading frames P7 / P8 would keep: the annotated start to its first in-frame stop; with start and end in
        one frame the first start codon after the last internal stop to the annotated end; with different frames
        (annotated end a complete stop) the in-frame start nearest the annotated start in the end's frame (P7)."""
        st = self.internal_stops(s, a, b)
        out = [(a, st[0] + 2)] if st else []
        if (b - a + 1) % 3 == 0 and st:
            c = next((c for c in range(st[-1] + 3, b - 2, 3) if s[c - 1:c + 2] in self.starts), None)
            if c is not None:
                out.append((c, b))
        elif (b - a + 1) % 3 and s[b - 3:b] in self.stops:
            best, c = None, b - 5
            while c >= max(1, a - 3 * self.a.max_start_shift) and s[c - 1:c + 2] not in self.stops:
                if s[c - 1:c + 2] in self.starts and (best is None or abs(c - a) < abs(best - a)):
                    best = c
                c -= 3
            if best is not None:
                out.append((best, b))
        return out

    def other_frame_orf(self, f: SeqFeature, s: str, a: int, b: int) -> Optional[Tuple[int, int]]:
        """(start, end) in transcript coordinates of one uninterrupted ORF, in a frame other than the annotated start's,
        that covers >= FRAMESHIFT_MAX_TRIM of the annotated CDS: start codon from max_start_shift codons upstream of a
        (not into the preceding gene) to FRAMESHIFT_WINDOW nt after it, complete stop codon within FRAMESHIFT_WINDOW nt
        of b; the start nearest a. A real extra nucleotide (bird / turtle nad3) leaves no such ORF - each frame holds
        only part of the gene."""
        n, L = b - a + 1, len(s)
        lo = max(1, self.upstream_end(f) + 1, a - 3 * self.a.max_start_shift)
        best = None
        for c in range(lo, min(L - 2, a + FRAMESHIFT_WINDOW) + 1):
            if (c - a) % 3 == 0 or s[c - 1:c + 2] not in self.starts:
                continue
            q = next((x for x in range(c + 3, L - 1, 3) if s[x - 1:x + 2] in self.stops), None)
            if q is None:
                continue
            e = q + 2                                  # last base of the stop codon
            if abs(e - b) > FRAMESHIFT_WINDOW or min(e, b) - max(c, a) + 1 < FRAMESHIFT_MAX_TRIM * n:
                continue
            if best is None or abs(c - a) < abs(best[0] - a):
                best = (c, e)
        return best

    def frameshift(self, f: SeqFeature) -> None:
        """P9: a gene with one untranslated extra nucleotide - bird / turtle nad3 at position 174 (Mindell et al. 1998),
        sometimes other genes (NC_010095.1 grebe nad1, NC_017607.1 turtle nad4, NC_010769.1 sponge cox3 / nad6) - is
        annotated by NCBI as join(start..p-1,p+1..end): a single-frame CDS keeps internal stops (table2asn InternalStop,
        109 nad3 in the 2026-10 runs) or loses half the gene. Applied when no single-frame repair (P7 / P8) keeps
        FRAMESHIFT_MAX_TRIM of the CDS and a start codon and a stop codon within FRAMESHIFT_WINDOW nt of the annotated
        ends make an open frame over one skipped nucleotide (both parts >= 30 nt). The ends nearest the annotated ones are
        taken; the skip at nad3 position ND3_EXTRA, else the middle of the possible positions. 2026-10 runs: 97 of 107
        nad3 both ends as NCBI, 9 more 1 nt off before --abbreviated_overlap (now applied here as in P3)."""
        if not simple(f):
            return
        s = self.tseq(f)
        a, b = self.tx(f)
        if not self.internal_stops(s, a, b):
            return
        n = b - a + 1
        if any(nb - na + 1 >= FRAMESHIFT_MAX_TRIM * n for na, nb in self.single_frame_repairs(s, a, b)):
            return
        alt = self.other_frame_orf(f, s, a, b)
        if alt:                                        # the core's frame was wrong: no frameshift (Travisia ATP8)
            S, E = alt
            old = span(f)
            gs, ge = self.from_tx(f, S, E)
            f.location = FeatureLocation(gs - 1, ge, f.location.strand)
            g = gene_feature(self.rec, f)
            if g is not None:
                g.location = FeatureLocation(gs - 1, ge, f.location.strand)
            self.log(f, "P9", old, span(f), f"no frameshift: one uninterrupted ORF in another frame covers "
                     f"{(min(E, b) - max(S, a) + 1) / n:.0%} of the CDS")
            return
        name = PCG_SYN.get(key(f.qualifiers.get("gene", [""])[0]), "")
        L = len(s)
        starts = [c for c in range(max(1, a - FRAMESHIFT_WINDOW), min(L - 2, a + FRAMESHIFT_WINDOW) + 1)
                  if s[c - 1:c + 2] in self.starts]
        ends = [e for e in range(max(3, b - FRAMESHIFT_WINDOW), min(L, b + FRAMESHIFT_WINDOW) + 1) if s[e - 3:e] in self.stops]
        best = None
        for S in starts:
            q5 = next((c for c in range(S + 3, L - 1, 3) if s[c - 1:c + 2] in self.stops), None)   # first stop, S's frame
            if q5 is None:
                continue
            for E in ends:
                if (E - S) % 3 or E - S < 120:
                    continue
                # codons after the skipped nucleotide p are read in E's frame; the 5' part in S's frame up to p
                l3 = next((c for c in range(E - 5, S, -3) if s[c - 1:c + 2] in self.stops), S)   # last internal stop, E's frame
                ps = []
                for p in range(max(l3, S + 30), min(q5 + 2, E - 30) + 1):
                    k = (p - S) % 3                    # codon broken by the skip: (p-1, p+1, p+2) or (p-2, p-1, p+1)
                    mixed = "" if k == 0 else (s[p - 2] + s[p:p + 2] if k == 1 else s[p - 3:p - 1] + s[p])
                    if mixed not in self.stops:
                        ps.append(p)
                if not ps:
                    continue
                target = S + ND3_EXTRA - 1 if name == "ND3" else (ps[0] + ps[-1]) / 2
                p = min(ps, key=lambda x: (abs(x - target), x))
                cand = (abs(S - a) + abs(E - b), abs(p - target), S, p, E)
                if best is None or cand < best:
                    best = cand
        if best is None:
            return
        _, _, S, p, E = best
        nt = s[S - 1:p - 1] + s[p:E]
        prot = str(Seq(nt).translate(table=self.code))
        if len(nt) % 3 or "*" in prot[:-1] or not prot.endswith("*"):    # cannot happen: the search checked every codon
            return
        why = f"one nucleotide (CDS position {p - S + 1}) not translated: join over a frameshift"
        nxt = self.downstream_trna_start(f)
        tail = 0
        if getattr(self.a, "abbreviated_overlap", False) and nxt is not None and 1 <= E - nxt + 1 <= 2:
            t = (nxt - 1 - S) % 3                      # vertebrate convention as in P3: T / TA before the next tRNA
                                                       # (the join is nxt-1-S nt long: one nucleotide skipped)
            if t in (1, 2) and s[nxt - 1 - t:nxt - 1] == ("T", "TA")[t - 1]:
                E, tail = nxt - 1, t
        old = span(f)
        g = gene_feature(self.rec, f)
        if f.location.strand == 1:
            f.location = CompoundLocation([FeatureLocation(S - 1, p - 1, 1), FeatureLocation(p, E, 1)])
        else:
            f.location = CompoundLocation([FeatureLocation(self.L + 1 - p, self.L + 1 - S, -1),
                                           FeatureLocation(self.L - E, self.L - p, -1)])
        if g is not None:                              # the gene spans the whole join (NCBI)
            gs, ge = span(f)
            g.location = FeatureLocation(gs - 1, ge, strand=f.location.strand)
        internal_excepts(f)
        notes = [x for x in f.qualifiers.get("note", []) if x not in (ND3_NOTE, FRAMESHIFT_NOTE, POLYA_NOTE)]
        notes.append(ND3_NOTE if name == "ND3" and self.code == 2 and abs(p - S + 1 - ND3_EXTRA) <= 10 else FRAMESHIFT_NOTE)
        if tail:
            gs, ge = span(f)
            f.qualifiers["transl_except"] = internal_excepts(f) + [transl_except_term(gs, ge, 1 if f.location.strand == 1
                                                                                      else -1, tail)]
            notes.append(POLYA_NOTE)
        f.qualifiers["note"] = notes
        self.joined.add(id(f))
        self.log(f, "P9", old, span(f), why)

    def stop_trim(self, f: SeqFeature) -> None:
        """P8: start and end in one frame with internal stop codons - a start placed upstream across stops (NC_027667.1
        nad1 72 nt, insect nad4L ~30 nt), an end beyond the gene's stop (NC_024927.1 nad4 348 nt, NC_073527.1 atp8 into
        atp6): the longer of the annotated start to its first stop and the first start codon after the last internal
        stop to the annotated end, when it keeps STOP_TRIM_MIN of the length (2026-10 runs: 30 of 33 such genes then
        start or end as NCBI; a gene still with internal stops is written /pseudo by P10)."""
        if not simple(f):
            return
        s = self.tseq(f)
        a, b = self.tx(f)
        if (b - a + 1) % 3 or not self.internal_stops(s, a, b):
            return
        na, nb = max(self.single_frame_repairs(s, a, b), key=lambda x: x[1] - x[0])
        if nb - na + 1 < max(60, STOP_TRIM_MIN * (b - a + 1)):
            return
        old = span(f)
        gs, ge = self.from_tx(f, na, nb)
        set_span(f, gs, ge)
        g = gene_feature(self.rec, f)
        if g is not None and simple(g):
            set_span(g, gs, ge)
        self.log(f, "P8", old, (gs, ge), "internal stop codons: " + ("start kept, first in-frame stop" if na == a else
                                                                     "end kept, first start codon after the last stop"))

    def pcg_check(self, f: SeqFeature) -> None:
        """P4: report problems that rules cannot fix (no change is made). P10: a CDS that still has internal stops is
        written /pseudo (gene and CDS; no translation) - NCBI rejects a CDS that cannot be translated (table2asn
        InternalStop), and its frame needs a manual check."""
        if id(f) in self.joined:
            return
        if not simple(f):
            if self.origin_span(f) is not None:
                self.origin_join_check(f)
            else:
                self.log(f, "P4", span(f), span(f), "multi-interval CDS (intron / origin) - not curated")
            return
        s = self.tseq(f)
        a, b = self.tx(f)
        nt = s[a - 1:b]
        prot = apply_excepts(f, codon_translate(nt, self.code))      # a codon read as an amino acid is no stop
        n_stop = prot.rstrip("*").count("*")
        if n_stop:
            self.log(f, "P4", span(f), span(f), f"{n_stop} internal stop codon(s): wrong frame, frameshift "
                     "(assembly indel?) or unannotated intron - check manually")
            for x in (f, gene_feature(self.rec, f)):
                if x is not None:
                    x.qualifiers["pseudo"] = [None]    # value-less /pseudo
            f.qualifiers["note"] = [n for n in f.qualifiers.get("note", []) if n != PSEUDO_NOTE] + [PSEUDO_NOTE]
            self.log(f, "P10", span(f), span(f), "CDS with internal stop codons written /pseudo (not translated)")
            return
        if nt[:3] not in self.starts:
            # P5: NCBI rejects an illegal start (SEQ_FEAT.StartCodon); RefSeq writes such CDS 5'-partial with
            # "start codon not determined" (e.g. echinoderm nad4l ATT/ATC under code 9: NC_001627.1, NC_061414.1)
            loc = f.location
            f.location = FeatureLocation(BeforePosition(loc.start), loc.end, loc.strand) if loc.strand == 1 else \
                FeatureLocation(loc.start, AfterPosition(loc.end), loc.strand)
            f.qualifiers["note"] = [n for n in f.qualifiers.get("note", []) if n != NO_START_NOTE] + [NO_START_NOTE]
            self.log(f, "P5", span(f), span(f), f"first codon {nt[:3]} is not a start codon of table "
                     f"{self.code}: CDS marked 5'-partial, 'start codon not determined' (NCBI practice)")
        tail = len(nt) % 3
        has_stop = ends_with_stop(nt, self.stops)              # (the last 3 nt are out of frame when there is a tail)
        if tail and not has_stop:
            self.log(f, "P4", span(f), span(f), f"length {len(nt)} is not a multiple of 3 and there is no T/TA tail - "
                     "frameshift, wrong end or unannotated intron - check manually")
        if not n_stop and not has_stop:
            # P6: no complete or abbreviated stop -> table2asn SEQ_FEAT.NoStop; written 3'-partial like P5
            loc, old, why = f.location, span(f), ""
            if self.L - 2 <= b < self.L:
                # cut by the sequence end (linear record; circular when O1 found no stop): the trailing 1-2 nt belong to
                # the cut codon, and a partial end short of the sequence end is SEQ_FEAT.PartialProblemNotSpliceConsensus3Prime
                # (stress case 2026-10: nad2 15558..>16559 on a 16560-nt linear record)
                loc = FeatureLocation(loc.start, ExactPosition(self.L), 1) if loc.strand == 1 else \
                    FeatureLocation(ExactPosition(0), loc.end, -1)
                why = f"extended {self.L - b} nt to the sequence end; "
            f.location = FeatureLocation(loc.start, AfterPosition(loc.end), 1) if loc.strand == 1 else \
                FeatureLocation(BeforePosition(loc.start), loc.end, -1)
            f.qualifiers["note"] = [n for n in f.qualifiers.get("note", []) if n != NO_STOP_NOTE] + [NO_STOP_NOTE]
            self.log(f, "P6", old, span(f), why + "no complete or abbreviated stop codon: CDS marked 3'-partial, "
                     "'stop codon not determined' (NCBI practice) - check manually")

    def translation(self, f: SeqFeature) -> None:
        """Review C finding 9: /translation of every CDS (simple or join) recomputed from its final location, never
        dropped (a record is always rewritten). A start codon of the table gives M unless the CDS is 5'-partial
        (NCBI convention). Not written when an internal stop remains (P4 reports it)."""
        f.qualifiers.pop("translation", None)
        if "pseudo" in f.qualifiers or "pseudogene" in f.qualifiers:
            return                                      # a /pseudo CDS has no /translation
        nt = str(f.extract(self.rec.seq)).upper()[int(f.qualifiers.get("codon_start", ["1"])[0]) - 1:]
        nt = nt[: len(nt) - len(nt) % 3]
        if not nt:
            return
        prot = str(Seq(nt).translate(table=self.code))
        prot = apply_excepts(f, prot, int(f.qualifiers.get("codon_start", ["1"])[0]) - 1)
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
        """D1: keep one copy of a protein gene annotated at several non-overlapping places, unless the other copy is
        complete and >= 95 % identical to it (length within 3 nt): a real gene copy, which NCBI annotates too (1,050
        records 2026-10-29: 30 of 38 identical protein-gene copies; heron / cormorant / turaco nad6, Hexamermis atp6).
        An antisense copy next to a complete gene that is incomplete or has internal stops is a spurious hit."""
        by: Dict[str, List[SeqFeature]] = {}
        for f in self.rec.features:
            if f.type == "CDS" and simple(f):
                name = PCG_SYN.get(key(f.qualifiers.get("gene", [""])[0]), "")
                if name:
                    by.setdefault(name, []).append(f)
        def quality(f: SeqFeature) -> tuple:
            nt = self.tseq(f)[self.tx(f)[0] - 1:self.tx(f)[1]]
            tail = len(nt) % 3
            codons = [nt[i:i + 3] for i in range(0, len(nt) - tail, 3)]
            complete = codons and codons[0] in self.starts and ends_with_stop(nt, self.stops)
            internal = any(c in self.stops for c in (codons if tail else codons[:-1]))
            return (bool(complete) and not internal, len(nt), -int(f.location.start))

        def sound(f: SeqFeature) -> bool:
            """A start codon and no stop inside: the stop may still be missing, as the core ends a gene at its last
            aligned codon and the stop rules (P3) run later (NC_044480.1 cob, NC_010636.1 cox2 copies)."""
            nt = self.tseq(f)[self.tx(f)[0] - 1:self.tx(f)[1]]
            codons = [nt[i:i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]
            return bool(codons) and codons[0] in self.starts and not any(c in self.stops for c in codons[:-1])

        def true_copy(f: SeqFeature, keep: SeqFeature) -> bool:
            fn, kn = (self.tseq(x)[self.tx(x)[0] - 1:self.tx(x)[1]] for x in (f, keep))
            return sound(keep) and sound(f) and abs(len(fn) - len(kn)) <= 3 and \
                sum(x == y for x, y in zip(fn, kn)) >= 0.95 * max(len(fn), len(kn))

        keeps = {name: max(fs, key=quality) for name, fs in by.items() if len(fs) >= 2}
        for name, fs in by.items():
            if len(fs) < 2:
                continue
            keep = keeps[name]
            for f in fs:
                if f is keep or any(a <= d and c <= b for a, b in spans(f) for c, d in spans(keep)):
                    continue                            # overlapping copies are the P2 case, not duplicates
                if true_copy(f, keep):
                    continue                            # a complete, near-identical second copy: a real gene copy
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
