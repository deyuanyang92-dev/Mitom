"""Constants and feature helpers shared by the rules of mito_curate.py (limits, notes, span / gene_feature, /transl_except)."""
from __future__ import annotations

import re
from typing import List, Optional, Tuple

from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

from mito_taxa import METAZOAN_PHYLA, PHYLUM_SYNONYMS

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



_EXCEPT_RE = re.compile(r"pos:(complement\()?(\d+)\.\.(\d+)\)?,aa:(\w+)")


def has_internal_except(f: SeqFeature) -> bool:
    """Whether the CDS carries a /transl_except other than the TERM of its stop (an amino acid read for a codon)."""
    for q in f.qualifiers.get("transl_except", []):
        m = _EXCEPT_RE.search(q.replace(" ", ""))
        if m and m.group(4).upper() != "TERM":
            return True
    return False


def internal_excepts(f: SeqFeature) -> List[str]:
    """Keep the /transl_except of codons inside the gene (an ambiguous codon read as Trp, NCBI NC_001636.1 ATP8) and
    drop the TERM one of the old stop (written again by the caller) and any outside the final location; returns the
    kept ones (2026-10-31 review: every /transl_except was deleted - X in the translations)."""
    pos = {i for p in f.location.parts for i in range(int(p.start), int(p.end))}
    keep = []
    for q in f.qualifiers.get("transl_except", []):
        m = _EXCEPT_RE.search(q.replace(" ", ""))
        if m and m.group(4).upper() != "TERM" and int(m.group(2)) - 1 in pos and int(m.group(3)) - 1 in pos:
            keep.append(q)
    if keep:
        f.qualifiers["transl_except"] = keep
    else:
        f.qualifiers.pop("transl_except", None)
    return list(keep)


def apply_excepts(f: SeqFeature, prot: str, offset: int = 0) -> str:
    """prot with the amino acids of the feature's internal /transl_except (Trp for an NCC codon ...)."""
    from Bio.SeqUtils import seq1
    order = [i for p in f.location.parts for i in (range(int(p.start), int(p.end)) if p.strand != -1
                                                    else range(int(p.end) - 1, int(p.start) - 1, -1))]
    index = {g: k for k, g in enumerate(order)}
    out = list(prot)
    for q in f.qualifiers.get("transl_except", []):
        m = _EXCEPT_RE.search(q.replace(" ", ""))
        if not m or m.group(4).upper() == "TERM":
            continue
        first = int(m.group(3)) - 1 if m.group(1) else int(m.group(2)) - 1      # the codon's first base (transcript)
        k = index.get(first)
        if k is None or (k - offset) % 3:
            continue
        aa = seq1(m.group(4)) if len(m.group(4)) == 3 else m.group(4)[:1]
        if 0 <= (k - offset) // 3 < len(out) and aa and aa != "X":
            out[(k - offset) // 3] = aa
    return "".join(out)
