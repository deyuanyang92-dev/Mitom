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

from Bio.Data import CodonTable
from Bio.SeqFeature import SeqFeature

import mito_env  # noqa: E402  (finds external programs: config, PATH, conda environments)
# N1: NCBI gene symbols, products and synonyms; tRNA names trnL1/L2, trnS1/S2 by anticodon (mito_names.py)
from mito_names import PCG_SYN, name_key as key  # noqa: E402
from mito_common import iter_gb_files, write_gb, write_tsv  # noqa: E402
import mito_common  # noqa: E402
from mito_common import InputError, cli_guard, prepare_input, require_records  # noqa: E402
import mito_start  # noqa: E402  (P11: the start codon NCBI records choose, learned)
import mito_tbl  # noqa: E402  (record_code: the genetic code of a record's taxonomy)
from mito_curate_base import (  # noqa: F401  (re-exported: the rule constants and helpers)
    _EXCEPT_RE,
    apply_excepts,
    BICISTRONIC,
    FRAMESHIFT_MAX_TRIM,
    FRAMESHIFT_NOTE,
    FRAMESHIFT_WINDOW,
    gene_feature,
    internal_excepts,
    INTERNAL_NOTE,
    MAX_END_EXTENSION,
    MAX_FEATURE_OVERLAP,
    MAX_INSIDE_TRIM,
    MAX_ORIGIN_EXT,
    MAX_ORIGIN_OVERLAP,
    MAX_PCG_OVERLAP,
    MAX_START_MOVE,
    ND3_EXTRA,
    ND3_NOTE,
    NO_START_NOTE,
    NO_STOP_NOTE,
    OVERLAP_FALLBACK,
    PSEUDO_NOTE,
    record_phylum,
    set_span,
    simple,
    SPACERS,
    span,
    spans,
    STOP_TRIM_MIN,
    WEAK_TRNA_E,
)
from mito_curate_origin import OriginRules
from mito_curate_pcg import PcgRules
from mito_curate_rna import RnaRules


class Curator(OriginRules, PcgRules, RnaRules):
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
            if f.type == "CDS" and str(f.qualifiers.get("codon_start", ["1"])[0]) != "1":
                # the boundary rules read codons from the first base: a 5'-partial CDS with /codon_start 2 / 3 was
                # shifted, joined over 'frameshifts' or written /pseudo (Grok 1.11, 2026-11-05); only translated
                self.log(f, "P4", span(f), span(f), f"/codon_start={f.qualifiers['codon_start'][0]}: not curated")
                continue
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


def record_code(rec) -> Optional[int]:
    """The genetic code of a record: the most common CDS /transl_table (None without one)."""
    codes = [int(f.qualifiers["transl_table"][0]) for f in rec.features
             if f.type == "CDS" and str(f.qualifiers.get("transl_table", [""])[0]).isdigit()]
    return max(set(codes), key=codes.count) if codes else None


def build_parser() -> argparse.ArgumentParser:
    """The command line (main); options() gives the same options without one."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--genetic_code", type=int, default=None,
                   help="Genetic code (default: each record's CDS /transl_table, else the usual code of its phylum, else 5; 2026-10-31 review: code 5 was "
                        "used for every record - vertebrate AGG stops read through)")
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
    p.add_argument("--allow_empty", action="store_true", help=argparse.SUPPRESS)   # mitom run: a folder of 0 records
    return p


def options(**given) -> argparse.Namespace:
    """Options for run(): input and out are required, e.g. options(input="04.final_gb", out="05.curated")."""
    return mito_common.options(build_parser(), **given)


def run(a: argparse.Namespace) -> int:
    """Curate the GenBank file / folder a.input into a.out (curation_report.tsv; table2asn before / after unless
    a.no_validate). a: build_parser() or options(input=..., out=...). Returns the exit code."""
    from Bio import SeqIO
    src, out = Path(prepare_input(a.input)), Path(a.out)
    a.start_model_obj = None if a.start_model == "no" else mito_start.load(None if a.start_model == "auto" else a.start_model)
    if a.start_model not in ("auto", "no") and a.start_model_obj is None:
        raise InputError(f"--start_model {a.start_model}: no such model")
    start_rows: List[Dict] = []
    t2a = None if a.no_validate else mito_env.find("table2asn", a.table2asn)
    if t2a:
        import ncbi_validate as nv
    val_before: List[Dict[str, str]] = []
    val_after: List[Dict[str, str]] = []
    files = iter_gb_files(src) if src.is_dir() else [src]
    out.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, str]] = []
    n_recs = 0
    for f in files:
        recs = list(SeqIO.parse(str(f), "genbank"))
        n_recs += len(recs)
        codes: List[int] = []
        for i, rec in enumerate(recs):
            # no --genetic_code and no CDS /transl_table: the code of the record's taxonomy as gb2tbl takes it (Grok
            # 1.3, 2026-11-05: always 5 - a vertebrate's AGR stops read through); written on the CDS that lack one
            # the fallback 5 is used, never written: a guess must not look like the record's code (code review, 2026-11-06)
            known = a.genetic_code or record_code(rec) or mito_tbl.lineage_code(rec)[0]
            code = known or 5
            if t2a:
                val_before += nv.validate_record(copy.deepcopy(rec), out / "validation" / "before", t2a, code=code)
            for ft in rec.features if known else ():     # before the backup: kept when curation fails (Bugbot, PR #3)
                if ft.type == "CDS" and not str(ft.qualifiers.get("transl_table", [""])[0]).isdigit():
                    ft.qualifiers["transl_table"] = [str(code)]
            codes.append(code)
            backup = copy.deepcopy(rec)                  # run() edits in place: a failure must not leave a half-curated record
            try:
                cur = Curator(rec, code, a)
                rows += cur.run()
                start_rows += cur.start_rows
            except Exception as exc:
                recs[i] = backup
                rows.append({"record": rec.id, "rule": "error", "reason": f"left unchanged: {exc!r}"})
        if t2a:
            for rec, code in zip(recs, codes):            # the same code as before: compare like with like
                val_after += nv.validate_record(rec, out / "validation" / "after", t2a, code=code)
        dest = out / (f.relative_to(src) if src.is_dir() else f.name)
        dest.parent.mkdir(parents=True, exist_ok=True)
        write_gb(recs, dest)      # always: N1 names / gene ranges change without a report row
    if not n_recs and not (src.is_dir() and a.allow_empty):
        require_records(0, "GenBank records", a.input)
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
        # ERROR = every severity that stops a submission (REJECT / CRITICAL, FATAL too; code review, 2026-11-06)
        sev = lambda r: "ERROR" if r["severity"] in nv.BLOCKING else r["severity"]   # noqa: E731
        nb = Counter(sev(r) for r in val_before)
        na = Counter(sev(r) for r in val_after)
        st = Counter((sev(r), r["status"]) for r in cmp_rows)
        print(f"[mito_curate] NCBI table2asn before -> after: ERROR {nb['ERROR']} -> {na['ERROR']}, "
              f"WARNING {nb['WARNING']} -> {na['WARNING']}; ERROR types fixed {st[('ERROR', 'fixed')]}, "
              f"remaining {st[('ERROR', 'remaining')]}, new {st[('ERROR', 'new')]} -> {out / 'validation_compare.tsv'}",
              file=sys.stderr)
    return 0


def main(argv=None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    sys.exit(cli_guard(main))
