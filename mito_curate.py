#!/usr/bin/env python3
"""mito_curate.py - rule-based correction of mitochondrial genome annotation (PCGs, 12S/16S rRNA).

Applies the annotation criteria of Fourdrilis et al. 2018 (Sci. Rep. 8:17964) as summarised by
Ghiselli et al. 2021 (Phil. Trans. R. Soc. B 376:20200159, doi:10.1098/rstb.2020.0159, section 3)
to automatic annotation (MitoM's core, or another annotator's GenBank such as MITOS output):

  P1  PCG start: the first eligible in-frame start codon at the 5' end, i.e. the start codon
      nearest the preceding gene WITHOUT overlapping it (criterion i). Only applied when the
      annotated start overlaps the preceding same-strand feature (move downstream), or when an
      in-frame start lies upstream before the preceding feature and the open frame reaches it
      (move upstream, at most --max_start_shift codons; only when the annotated start codon is
      not a valid start, or for --aggressive_start_genes; annelids: nad5 and nad4L; or, unless
      --no_upstream_atg, to the most upstream in-frame ATG when the start codon is not ATG). A start that
      overlaps the preceding gene is moved downstream only when a start codon lies within 15 nt. A start placed
      upstream of a whole same-strand tRNA (in the gene's first third) moves to the first in-frame start codon after
      it (table2asn BadCDScomponentOverlapTRNA: unionid nad1 across trnL2).
  P11 PCG start, after P1: of the in-frame start codons within 30 codons of P1's start (reading frame open in
      between), the one NCBI records would choose, as scored by gradient-boosted trees learned from 1,683 NCBI
      records (mito_start.py, models/start_model/start_model.json; --start_model no turns it off): codon, distance to
      the preceding feature, the core's alignment (residues missing at the N-terminus), gene, phylum. Cross-validated
      by family: exact starts 89.1 -> 90.7 % (arthropods 77.8 -> 82.3 %, cnidarians 79.8 -> 83.2 %).
  P2  Two PCGs on the same strand and in the same frame cannot overlap (criterion ii) - reported.
  D1  A protein gene annotated twice at different places (e.g. MitoZ on a duplicated region): only one
      copy is kept - the complete one (start + stop, no internal stop), the longer, then the first in
      the record; the other CDS and its gene feature are removed (NCBI records annotate one copy), unless
      the other copy is complete too and its sequence is >= 95 % identical with the same length +-3 nt, for
      >= 2 genes (a duplicated block, annotated twice by NCBI: NC_009093.1 mite, 9 genes; one repeated
      atp8 is annotated once: NC_049135.1).
      An overlapping copy on the other strand that is incomplete or has internal stops is removed as well.
  P3  PCG end: the first in-frame complete stop codon, or an abbreviated stop (T / TA) ending
      immediately before the downstream same-strand tRNA, completed by polyadenylation
      (criterion iii; tRNA punctuation model, Ojala et al. 1981). /transl_except + note added.
      When no tRNA follows and the first complete stop lies > 30 nt inside the next same-strand
      protein gene (not the bicistronic pairs atp8-atp6, nad4L-nad4; rRNA ends are trimmed by R1
      instead; > 10 nt for any other feature; an rRNA start does not bound an end without any stop signal, R1 moves
      the rRNA), the gene is not extended into it: the annotated end is kept when it is an abbreviated T / TA stop (the end of
      the protein homology, e.g. cox1 before cox2), else the T / TA right before the next gene
      (e.g. nad6 before cob) is used. An end that stops short of the first complete stop with
      nothing annotated in between (core homology lost at the C-terminus, 100-500 nt in divergent
      phyla) is extended to that stop (at most 600 nt). A gene that would contain a whole same-strand tRNA and has
      no T / TA before it ends right before the tRNA (3'-partial).
  P7  start and end in different frames (internal stops): replaced by the longer ORF - annotated start
      to its first stop, or annotated stop with the in-frame start nearest the annotated start - when
      it keeps >= half the length (metazoa benchmark 2026-10: 14 CDS with internal stops, NCBI 0).
  P8  start and end in one frame with internal stops (a start placed upstream across a stop, an end past the gene's
      stop): the longer of start .. first stop and first start codon after the last internal stop .. end (>= half).
  P9  one untranslated extra nucleotide (bird / turtle nad3 at position 174, Mindell et al. 1998; rarely other genes):
      join(start..p-1,p+1..end) as NCBI writes it, with its note, when no single-frame repair keeps 85 % of the gene
      and a start / stop codon within 30 nt of the annotated ends make an open frame over one skipped nucleotide.
  P10 a CDS that still has internal stops is written /pseudo (gene and CDS, not translated; table2asn InternalStop).
  O2  a join across the origin of a circular record written by the core with <= 300 nt past the origin is cut back to
      the sequence end; P1 / P3 / O1 place its end again (annelid nad2 before cox1 ended without a stop: NoStop).
  T1  a weak tRNA call (MiTFi E > 1e-3) lying wholly inside a same-strand protein gene is removed when it is a second
      copy or lies in the gene's middle third; in its first / last third P1 / P3 move the gene end instead.
  C1  a tRNA still lying wholly inside a same-strand protein gene: the gene is removed when it has no start or stop
      codon or is /pseudo (a false call), else the conflict is reported.
  O1  circular record: a CDS that ends without a stop <= 2 nt before the transcript end, or further when its
      open frame runs to the end with nothing annotated in between (the core ends at the last aligned codon), is
      continued across the origin to its first in-frame stop (join; genes cut by the cox1-first rotation); at most
      60 nt where another feature starts at the origin (genuine overlaps are <= ~50 nt), else 300 nt.
  P5/P6  illegal start / no stop codon: CDS marked 5'- / 3'-partial with 'start / stop codon not
      determined' (NCBI practice; table2asn SEQ_FEAT.StartCodon / NoStop otherwise). A stop-less CDS ending 1-2 nt
      before the sequence end is first extended to it (the cut codon). A single-interval CDS the core already marked
      partial (< / >) goes through P1 / P3 / O1 like any other; only joins are left alone (P4).
  R1  rRNA (12S rrnS, 16S rrnL) must not overlap flanking genes; boundaries are trimmed to the
      neighbours (criterion vi) and, by default, extended to the flanking genes (gaps <= --max_fill
      nt, the GenBank convention; --rrna_bounds fill). --rrna_bounds trim only removes overlaps.
  R3  rRNA names standardised (rrnS / rrnL; 12S / 16S ribosomal RNA); rRNA with a large CDS
      overlap is reported.
  N1  NCBI gene names and products (COX1 ... ND4L, rrnS / rrnL, trnX with L1 / L2, S1 / S2 by anticodon); the
      notes the annotation core writes for itself (reference protein and score, model columns, MitoM-ML / MiTFi
      scores) are removed, so the GenBank and the NCBI submission files carry only biological notes.
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
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, ExactPosition, FeatureLocation, SeqFeature

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402  (finds external programs: config, PATH, conda environments)
# N1: NCBI gene symbols, products and synonyms; tRNA names trnL1/L2, trnS1/S2 by anticodon (mito_names.py)
from mito_names import AA1_TO_3, AA3, PCG_NAMES, PCG_SYN, RRNA_PRODUCT, name_key as key, trna_symbol  # noqa: E402
from mito_names import RRNA_ALIAS as RRNA_NAMES  # noqa: E402
from mito_common import (POLYA_NOTE, abbreviated_stop, codon_translate, ends_with_stop, iter_gb_files,  # noqa: E402
                         transl_except_term, write_tsv)
import mito_start  # noqa: E402  (P11: the start codon NCBI records choose, learned)
from mito_taxa import METAZOAN_PHYLA, PHYLUM_SYNONYMS  # noqa: E402
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
OVERLAP_FALLBACK = 20      # P3: no T / TA right before the next feature -> the first complete stop is taken when it runs at
                           # most this far into it (predicted tRNA ends are approximate; leave-out, 255 records: +14 exact
                           # ends, was marked 3'-partial; 40 nt: +3 more, kept small for the NCBI overlap checks)
MAX_END_EXTENSION = 600     # P3: longest extension of an end that stopped short of the first complete stop (metazoa
                            # benchmark 2026-10: the 50 confirmed cases were 100-512 nt)
MAX_FEATURE_OVERLAP = 10    # P3: a complete stop may reach this far into the next feature that is not a same-strand
                            # protein gene (173-record truth: into an antisense tRNA 1-10 nt 41 times, deeper 0; into a
                            # same-strand rRNA at most 2 nt)
MAX_PCG_OVERLAP = 30        # P3: a complete stop may reach this far into the next same-strand protein gene (173-record
                            # truth: 11-30 nt inside, complete stop 17 vs T/TA 6; 31-50 nt 6 vs 6; deeper T/TA) (NCBI annelid
                            # records: cox1-cox2, nad6-cob end with T/TA instead; atp8-atp6 / nad4L-nad4 overlap by design)
                            # - except annelids (--abbreviated_pcg_overlap: MAX_FEATURE_OVERLAP, as O1 at the origin):
                            # all local NCBI records 2026-10-30, a complete stop 11-30 nt inside the next same-strand
                            # protein gene: Annelida T/TA 142 vs complete 43 (nad2 before cox1 115 vs 32), Cnidaria 5 vs
                            # 109, Mollusca 12 vs 66, Arthropoda 40 vs 65; curation A/B on 1,718 leave-out records:
                            # annelid protein genes exact 7357 -> 7413 (nad2 +41, nad6 +21), no other phylum changed
BICISTRONIC = {("ATP8", "ATP6"), ("ND4L", "ND4")}
FRAMESHIFT_WINDOW = 30      # P9: nt around the annotated start / end searched for the start / stop codon of the join
FRAMESHIFT_MAX_TRIM = 0.85  # P9 only when no single-frame repair (P7 / P8) keeps this fraction of the CDS (2026-10 runs, 162 CDS
                            # with internal stops: the 107 bird / turtle nad3 and 5 other NCBI joins all lie below, one at
                            # 0.81 (grebe nad1); start / end errors NCBI annotates as one interval lie at 0.82-0.99, and
                            # those at 0.82-0.84 have no join)
ND3_EXTRA = 174             # P9: CDS position of the untranslated nucleotide of bird / turtle nad3 (NCBI: 174 in 66, 175 in
                            # 36 of 107 records; Mindell et al. 1998)
ND3_NOTE = "frameshift mechanism unknown (Mindell et al., 1998, Mol. Biol. Evol., 15:1568-1571)"   # NCBI's nad3 note
FRAMESHIFT_NOTE = "frameshift: 1 nt not translated (extra nucleotide or sequencing error)"
STOP_TRIM_MIN = 0.5         # P8: a CDS with internal stops is trimmed when the open part keeps this fraction (as P7)
PSEUDO_NOTE = "internal stop codon(s): not translated"
WEAK_TRNA_E = 1e-3          # T1: a MiTFi call with a larger E (or a second copy) lying inside a same-strand protein gene is
                            # dropped when it is a second copy or lies in the gene's middle third (2026-10 runs: 7 of 9 such calls absent from
                            # NCBI; the real ones lie at a gene end, where the gene end is moved instead)
# N1: notes the annotation core writes for itself (reference protein and score, model columns, MitoM-ML / MiTFi
# scores) - kept in the core's files, not in the GenBank / NCBI submission
INTERNAL_NOTE = re.compile(r"^(ref \S+ score \d|MitoM-ML p [\d.]|MiTFi \S+ bits [\d.]|nhmmer \S+ bits \d|cmsearch bits \d)")


def span(f: SeqFeature) -> Tuple[int, int]:
    return int(f.location.start) + 1, int(f.location.end)


def spans(f: SeqFeature) -> List[Tuple[int, int]]:
    """1-based (start, end) of every part. span() is the hull of a join: an origin-spanning feature
    join(L-29..L,1..35) is not a feature ending at L (review C finding 6)."""
    return [(int(p.start) + 1, int(p.end)) for p in f.location.parts]


def simple(f: SeqFeature) -> bool:
    """Single-interval feature (multi-exon / origin-spanning joins are left alone). Fuzzy ends (< / >) do not exclude
    it: the native core marks a CDS cut at a sequence end partial, and P1 / P3 / O1 must still run on it (review record
    A: `<100..700` + TAA got only P4). set_span writes exact positions; P5 / P6 re-mark what stays without start / stop."""
    return not isinstance(f.location, CompoundLocation)


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


def record_phylum(rec) -> str:
    """Phylum named in a record's taxonomy (a class / subphylum name counts for its phylum), '' if none."""
    tax = rec.annotations.get("taxonomy", []) or []
    return next((t for t in tax if t in METAZOAN_PHYLA), "") or next((PHYLUM_SYNONYMS[t] for t in tax
                                                                       if t in PHYLUM_SYNONYMS), "")


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
        self.weak: set = set()          # id() of tRNAs whose core note marks a weak call (T1), read before N1 drops the notes
        self.copies: set = set()        # ... and of those that are a second copy of a family
        self.joined: set = set()        # id() of CDS joined over an untranslated nucleotide (P9)
        genes = getattr(a, "aggressive_start_genes", None)
        self.aggressive_genes = {PCG_SYN.get(key(x), x.upper()) for x in genes.split(",") if x} if genes else set()
        self.aln: Dict[int, Tuple[int, int, int]] = {}   # id() of a CDS -> the core's alignment (first, last residue,
        self.core_start: Dict[int, int] = {}             # reference length) and its start before curation: P11
        self.start_model = getattr(a, "start_model_obj", None)
        self.start_rows: List[Dict] = []                 # --dump_starts
        self.phylum = record_phylum(rec) or getattr(a, "phylum", None) or ""

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

    def enclosed_trnas(self, f: SeqFeature, a: int, b: int) -> List[Tuple[int, int]]:
        """Transcript spans of the same-strand tRNAs lying wholly inside a..b (table2asn SEQ_FEAT.BadCDScomponentOverlapTRNA;
        T1 has dropped the weak ones in the middle of a gene)."""
        out = []
        for g in self.others(f, types=("tRNA",), same_strand=True):
            for gs, ge in spans(g):
                ts, te = (gs, ge) if f.location.strand == 1 else (self.L + 1 - ge, self.L + 1 - gs)
                if a <= ts and te <= b:
                    out.append((ts, te))
        return sorted(out)

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

    def next_feature(self, f: SeqFeature, types=("tRNA", "rRNA", "CDS")) -> Tuple[Optional[int], Optional[str]]:
        """Transcript position where the nearest downstream feature (tRNA, rRNA, CDS, either strand) begins, and
        its type."""
        s, e = span(f)
        cand = []
        for g in self.others(f, types=types):
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
    def overlap_fallback(self, first_stop, nxt) -> bool:
        return first_stop is not None and nxt is not None and first_stop - (nxt - 1) <= OVERLAP_FALLBACK

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
        for k in ("transl_except",):
            f.qualifiers.pop(k, None)
        loc = f.location                               # an end left 3'-partial (cut by a linear sequence end) may continue
        cut = isinstance(loc.end, AfterPosition) if loc.strand == 1 else isinstance(loc.start, BeforePosition)
        if abbreviated_stop(nt) and not cut:           # (table2asn PartialProblemHasStop with a T / TA stop)
            gs, ge = span(f)
            f.qualifiers["transl_except"] = [transl_except_term(gs, ge, 1 if f.location.strand == 1 else -1, tail)]
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
                    f.qualifiers["transl_except"] = [transl_except_term(gs, ge, 1 if f.location.strand == 1 else -1,
                                                                        tail)]
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
        f.location = FeatureLocation(head.start, head.end, head.strand)
        g = gene_feature(self.rec, f)
        if g is not None:
            g.location = f.location
        f.qualifiers.pop("translation", None)
        self.log(f, "O2", origin, span(f), f"the core's join across the origin cut back to the sequence end ({len(tail)} nt "
                 "past it); the end is placed again by P3 / O1")

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
        f.qualifiers.pop("transl_except", None)
        notes = [x for x in f.qualifiers.get("note", []) if x not in (ND3_NOTE, FRAMESHIFT_NOTE, POLYA_NOTE)]
        notes.append(ND3_NOTE if name == "ND3" and self.code == 2 and abs(p - S + 1 - ND3_EXTRA) <= 10 else FRAMESHIFT_NOTE)
        if tail:
            gs, ge = span(f)
            f.qualifiers["transl_except"] = [transl_except_term(gs, ge, 1 if f.location.strand == 1 else -1, tail)]
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

    def origin_join_check(self, f: SeqFeature) -> None:
        """P4 / P5 / P6 / P10 on a CDS joined across the origin that O1 / O2 did not end (the core's join whose part
        before the origin is short, so O2 cannot cut it back: leave-out runs, NC_006386.1 nad4L 38 nt before the
        origin on -, NC_023783.1 cox3; table2asn NoStop): the joined sequence is checked like a single interval -
        internal stops /pseudo, no start codon 5'-partial, no stop 3'-partial with the NCBI notes."""
        nt = str(f.extract(self.rec.seq)).upper()
        prot = codon_translate(nt, self.code)
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
        prot = codon_translate(nt, self.a.genetic_code)
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
                     f"{self.a.genetic_code}: CDS marked 5'-partial, 'start codon not determined' (NCBI practice)")
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

    def rrna(self, f: SeqFeature) -> None:
        # an rRNA across the origin is one interval s..e with e > L (MitoM's frame starts at cox1, so the echinoderm
        # rrnL before cox1 crosses it: NC_013881.1 ran 19 nt into cox1 and left 46 nt to nad2 unfilled)
        origin = self.origin_span(f)
        if not simple(f) and origin is None:
            return
        g = gene_feature(self.rec, f)
        s, e = old = origin or span(f)
        nf = [(sp, x) for x in self.others(f, types=("tRNA", "CDS", "rRNA")) for sp in spans(x)]
        if self.rec.annotations.get("topology") == "circular":        # neighbours across the origin
            nf += [((a + k, b + k), x) for (a, b), x in list(nf) for k in (-self.L, self.L)]
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
            # an end at the end of a linear sequence that R1 leaves in place keeps its < / >: the rRNA runs on beyond
            # the sequenced part (DQ517436.1, a partial genome starting inside rrnL)
            linear = self.rec.annotations.get("topology") != "circular"
            keep = linear and (isinstance(f.location.start, BeforePosition) and s == old[0] == 1,
                               isinstance(f.location.end, AfterPosition) and e == old[1] == self.L)
            self.set_circular_span(f, s, e)
            if keep and any(keep):
                f.location = FeatureLocation(BeforePosition(s - 1) if keep[0] else f.location.start,
                                             AfterPosition(e) if keep[1] else f.location.end, f.location.strand)
            if g is not None and (simple(g) or origin):
                g.location = f.location
            self.log(f, "R1", old, (s, e), "trimmed to flanking genes" + (" and filled to them" if self.a.rrna_bounds == "fill" else ""))

    def weak_trnas(self) -> None:
        """T1: a weak tRNA call (MiTFi E > WEAK_TRNA_E, or a second copy: the core's note) lying wholly inside a
        same-strand protein gene is dropped - tRNAs punctuate the transcripts, a whole one inside a gene on its strand is
        a false call or a wrong gene end (table2asn SEQ_FEAT.BadCDScomponentOverlapTRNA; 2026-10 runs: 7 of 9 such calls
        absent from NCBI). Only a second copy of a family, or a call in the middle third of the gene: one in its first or
        last third is more often a gene end placed across it (NC_060482.1 / NC_060708.1: nad6 extended to the stop
        inside trnD, NCBI ends nad6 right before trnD) and P1 / P3 move that end instead. Strong calls are never dropped.
        Run before the protein gene rules, else P3 ends a gene with the T / TA it finds right before such a call."""
        for t in [x for x in self.rec.features if x.type == "tRNA" and id(x) in self.weak]:
            host = next((c for c in self.rec.features if c.type == "CDS" and c.location.strand == t.location.strand
                         and all(any(cs <= ts and te <= ce for cs, ce in spans(c)) for ts, te in spans(t))), None)
            if host is None:
                continue
            a, b = self.tx(host) if simple(host) else (0, 0)
            ts, te = self.tx(t) if simple(t) else (0, 0)
            n = b - a + 1
            if simple(host) and simple(t) and id(t) not in self.copies and (ts - a < n / 3 or b - te < n / 3):
                continue                               # at a gene end: P1 / P3 move the end
            g = gene_feature(self.rec, t)
            self.rec.features = [x for x in self.rec.features if x is not t and x is not g]
            self.log(t, "T1", span(t), span(t), f"weak tRNA call inside {host.qualifiers.get('gene', ['?'])[0]} "
                     f"{span(host)[0]}..{span(host)[1]} (same strand) removed")

    def trna_conflicts(self) -> None:
        """C1: a tRNA still lying wholly inside a same-strand protein gene after P1 / P3 (table2asn
        BadCDScomponentOverlapTRNA): a gene without start or stop codon (P5 / P6) or not translatable (P10) is a false
        call and is removed (NC_061193.1 flatworm 'atp8' around tRNA-Glu, NC_013814.1 chaetognath 'atp6' around
        tRNA-Phe: neither phylum has the gene); a complete gene is kept and the conflict reported (P4)."""
        for c in [x for x in self.rec.features if x.type == "CDS"]:
            inside = [t for t in self.rec.features if t.type == "tRNA" and t.location.strand == c.location.strand
                      and span(c)[0] <= span(t)[0] and span(t)[1] <= span(c)[1] and simple(t)]
            if not inside:
                continue
            notes = c.qualifiers.get("note", [])
            if "pseudo" in c.qualifiers or NO_START_NOTE in notes or NO_STOP_NOTE in notes:
                g = gene_feature(self.rec, c)
                self.rec.features = [x for x in self.rec.features if x is not c and x is not g]
                self.log(c, "C1", span(c), span(c), f"incomplete or untranslatable gene around "
                         f"{inside[0].qualifiers.get('gene', ['?'])[0]} {span(inside[0])[0]}..{span(inside[0])[1]} removed")
            else:
                self.log(c, "P4", span(c), span(c), f"contains {inside[0].qualifiers.get('gene', ['?'])[0]} "
                         f"{span(inside[0])[0]}..{span(inside[0])[1]} on its strand - check manually")

    # ── N1: NCBI-style names (every feature, also joins / origin-spanning ones) ──
    def names(self, f: SeqFeature) -> None:
        g = gene_feature(self.rec, f)
        if f.type == "CDS" and simple(f):                  # P11 reads the core's alignment (its note is dropped below)
            m = re.search(r"residues (\d+)-(\d+)/(\d+)", " ".join(f.qualifiers.get("note", [])))
            if m:
                self.aln[id(f)] = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
            self.core_start[id(f)] = self.tx(f)[0]
        for n in f.qualifiers.get("note", []) if f.type == "tRNA" else []:
            m = re.match(r"MiTFi \S+ bits [\d.]+ E ([0-9.eE+-]+)( copy)?", n)
            if m and (float(m.group(1)) > WEAK_TRNA_E or m.group(2)):
                self.weak.add(id(f))                   # T1 (the note itself is dropped below)
                if m.group(2):
                    self.copies.add(id(f))
        notes = [n for n in f.qualifiers.get("note", []) if not INTERNAL_NOTE.match(n)]
        if notes:
            f.qualifiers["note"] = notes
        else:
            f.qualifiers.pop("note", None)
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
            # one of three readers of a tRNA's amino acid / anticodon, kept apart on purpose (see mito_names.trna_name)
            txt = " ".join(f.qualifiers.get("product", []) + [cur])
            m = re.search(r"tRNA-([A-Z][a-z]{2})", txt)
            aa = AA3.get(m.group(1), "") if m else ""
            if not aa and re.match(r"trn([A-Z])", cur):
                aa = cur[3]
            ac = ""
            m = re.search(r"\(([acgtuACGTU]{3})\)", cur) or \
                re.search(r"anticodon:\s*([acgtuACGTU]{3})", " ".join(f.qualifiers.get("note", []))) or \
                re.search(r"seq:([acgtuACGTU]{3})", " ".join(f.qualifiers.get("anticodon", [])))
            if m:
                ac = m.group(1).lower().replace("u", "t")
            if aa in AA1_TO_3:
                name = trna_symbol(aa, ac)
                if not ac and re.fullmatch(r"trn[LS][12]", cur):
                    name = cur                         # anticodon unknown: the core numbered the gene from the MiTFi model
                notes = [n for n in f.qualifiers.get("note", []) if not n.startswith("anticodon:")]
                q = {"gene": [name], "product": [f"tRNA-{AA1_TO_3[aa]}"],
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
        self.weak_trnas()
        for f in list(self.rec.features):
            if f.type == "CDS":
                self.origin_cut(f)
                self.frameshift(f)
                self.frame_repair(f)
                self.stop_trim(f)
                self.pcg(f)
                self.across_origin(f)
                self.pcg_check(f)
        self.trna_conflicts()
        self.pcg_overlaps()
        for f in list(self.rec.features):              # R1 again: P3 / O1 moved CDS ends, trim the rRNAs to them
            if f.type == "rRNA":
                self.rrna(f)
        for f in self.rec.features:
            if f.type == "CDS":
                self.translation(f)
        for f in self.rec.features:                    # gene feature = exactly its CDS / RNA (SEQ_FEAT.CDSgeneRange)
            if f.type in ("CDS", "rRNA", "tRNA") and id(f) not in self.joined:   # (a P9 gene spans its join)
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
    from Bio import SeqIO
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
    p.add_argument("--no_upstream_atg", action="store_true",
                   help="P1: keep a non-ATG start codon instead of moving to the most upstream in-frame ATG (mitom run "
                        "sets it for nematodes and arthropods)")
    p.add_argument("--abbreviated_overlap", action="store_true",
                   help="P3: a complete stop overlapping the next same-strand tRNA by 1-2 nt becomes the T / TA "
                        "before it (vertebrate RefSeq convention; before a protein gene the complete stop is kept: "
                        "atp6 TAA / cox3 ATG; mitom run sets it for the Chordata clade)")
    p.add_argument("--abbreviated_pcg_overlap", action="store_true",
                   help="P3: a complete stop more than 10 nt inside the next same-strand protein gene becomes the T / TA "
                        "before it, as O1 at the origin (NCBI annelid records; others allow 30 nt; mitom run sets it "
                        "for the Annelida clade)")
    p.add_argument("--aggressive_start_genes", default="",
                   help="Genes whose start is always moved to the in-frame start codon nearest the preceding gene "
                        "(default none; `mitom run` passes ND5,ND4L for the Annelida clade). Benchmark 2026-10: right "
                        "for nad5 / nad4L in every differing annelid record (40), wrong in most other phyla (173 "
                        "metazoa records: 17 right, 52 wrong).")
    p.add_argument("--start_model", default="auto",
                   help="P11: auto (default) = models/start_model/start_model.json when present; no = off; or a model "
                        "file (benchmark/train_start_model.py)")
    p.add_argument("--dump_starts", default=None,
                   help="write every P11 start candidate with its features to this TSV (training data)")
    p.add_argument("--phylum", default=None, help="phylum of records without a taxonomy (a feature of the start model)")
    p.add_argument("--table2asn", default=None,
                   help="NCBI table2asn for validation before/after curation (default: found by mito_env.py)")
    p.add_argument("--no_validate", action="store_true", help="Skip the NCBI validation before/after curation")
    a = p.parse_args(argv)
    src, out = Path(a.input), Path(a.out)
    a.start_model_obj = None if a.start_model == "no" else mito_start.load(None if a.start_model == "auto" else a.start_model)
    if a.start_model not in ("auto", "no") and a.start_model_obj is None:
        p.error(f"--start_model {a.start_model}: no such model")
    start_rows: List[Dict] = []
    t2a = None if a.no_validate else mito_env.find("table2asn", a.table2asn)
    if t2a:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import ncbi_validate as nv
    val_before: List[Dict[str, str]] = []
    val_after: List[Dict[str, str]] = []
    files = iter_gb_files(src) if src.is_dir() else [src]
    out.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, str]] = []
    for f in files:
        recs = list(SeqIO.parse(str(f), "genbank"))
        for i, rec in enumerate(recs):
            if t2a:
                val_before += nv.validate_record(copy.deepcopy(rec), out / "validation" / "before", t2a)
            backup = copy.deepcopy(rec)                  # run() edits in place: a failure must not leave a half-curated record
            try:
                cur = Curator(rec, a.genetic_code, a)
                rows += cur.run()
                start_rows += cur.start_rows
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
    write_tsv(out / "curation_report.tsv", cols, rows)
    if a.dump_starts:
        write_tsv(Path(a.dump_starts), list(start_rows[0]) if start_rows else ["record"], start_rows)
    from collections import Counter
    c = Counter(r["rule"] for r in rows)
    print(f"[mito_curate] {len(files)} file(s); changes: " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())),
          file=sys.stderr)
    if t2a:
        cmp_rows = compare_validation(val_before, val_after)
        vcols = ["record", "severity", "code", "feature", "before", "after", "status"]
        write_tsv(out / "validation_compare.tsv", vcols, cmp_rows)
        nb = Counter(r["severity"] for r in val_before)
        na = Counter(r["severity"] for r in val_after)
        st = Counter((r["severity"], r["status"]) for r in cmp_rows)
        print(f"[mito_curate] NCBI table2asn before -> after: ERROR {nb['ERROR']} -> {na['ERROR']}, "
              f"WARNING {nb['WARNING']} -> {na['WARNING']}; ERROR types fixed {st[('ERROR', 'fixed')]}, "
              f"remaining {st[('ERROR', 'remaining')]}, new {st[('ERROR', 'new')]} -> {out / 'validation_compare.tsv'}",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
