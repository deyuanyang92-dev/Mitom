#!/usr/bin/env python3
"""mito_rnaseq_annot.py - annotate with RNA-seq evidence: what the annelid genome-project RNA-seq taught, as one
algorithm. Separate from the MitoM annotation (`mitom run`) and from `mitom rna`; it reads an annotated mitogenome and
the same specimen's RNA-seq and writes a revised annotation with every decision labelled.

Knowledge behind each step: knowledge/topics/rnaseq_mito_annotation.md; measurements:
benchmark/results/annelid_rna_2026-11-10/README.md, benchmark/results/rna_mappers_2026-11-10/README.md.

Steps (each writes rows to rnaseq_evidence.tsv: step, feature, verdict, detail; verdicts confirmed / corrected /
flagged / not_testable):
  0  QC. Reads mapped with bwa mem (mito_rna.map_reads). Edit rate > DIVERGED: the RNA-seq is not this species (or
     contaminated / mislabelled) - stop, nothing is changed. Protein-gene median depth < MIN_DEPTH: expression not
     tested. Poly(A) reads < POLYA_FEW of the reads: 3' ends not tested. Features whose sequence occurs twice in the
     record (25-mers): depth not measurable.
  1  Introns from split reads (bwa primary + supplementary on one strand, MAPQ >= 20, gap > 20 nt), tallied without
     any annotation; gaps within JUNCTION_MERGE nt are one junction. An annotated intron with >= MIN_JUNCTION reads at
     its exact ends: confirmed; with none while its exons are covered: flagged. A frameshift join (gap < 20 nt, MitoM
     P9) whose skipped base is in the RNA: flagged, and an uninterrupted ORF in another frame covering >= 85 % of the
     gene is reported (Travisia forbesii ATP8). A strong junction not annotated: flagged as an intron candidate.
  2  Stop codon overlapping the next tRNA / gene: the mRNA is either cut at that gene's start, leaving T / TA (27 of
     50 cases in 23 clean pairs) or keeps the full stop (23). The poly(A) site decides: >= MIN_CUT_READS reads at the
     cut (+-1) and more than at the full stop: the CDS ends at the cut with an incomplete stop (corrected).
  3  3' ends: mito_rna.refine (confirmed at the end or after a short 3' UTR; shortened / extended to an in-frame T / TA;
     rRNA 3' ends) - rows copied here.
  4  RNA differing from the genome: an insertion / deletion in >= RNA_INDEL_SHARE of >= MIN_INDEL_DEPTH reads inside a
     gene (Travisia COX3: one T more in 35 % of the RNA, DNA of the specimen all 6 T) - flagged, the genome sequence
     is kept (a +1 / -1 in a gene frameshifts it).
  5  Expression. Short genes get fewer reads (library length bias: ATP8 ~0.07 of the median), so an ORF is called
     transcribed when its depth reaches the least-expressed of the short genes ATP8 / ND4L / ND3 x 0.5.

Usage: python3 mito_rnaseq_annot.py -g annotated.gb -1 R1.fq.gz [-2 R2.fq.gz] -o outdir --run SRR... [--threads 16]
       python3 mito_rnaseq_annot.py -g annotated.gb --sam mapped.sam -o outdir --run SRR...
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Sequence

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord

import mito_rna
from mito_common import (POLYA_NOTE, abbreviated_stop, cli_guard, codon_translate, transl_except_term, write_gb,
                         write_tsv)
from mito_names import gene_name

DIVERGED = 0.02            # edit rate: another species / contamination (same specimen <= 1.08 %, 30 annelid pairs)
MIN_DEPTH = 10             # protein-gene median depth for expression calls
MIN_POLYA_TOTAL = 100      # poly(A) reads below this: no 3' evidence (random primers / total RNA); an absolute
                           # count - a share of reads failed deep libraries (Sabella: 3,878 of 40.5 M reads)
MIN_JUNCTION = 10          # split reads for a junction ...
JUNCTION_RATIO = 0.2       # ... and >= this share of the depth at its flanks (annelid pairs: the 4 real introns
                           # 0.66-1.38; 305 of 311 other split-read gaps < 0.01 - chimeras / template switching)
MIN_FLANK_DEPTH = 50       # an unannotated junction is judged only at this flank depth
JUNCTION_MERGE = 5         # junction ends within 5 nt are one junction (alignment wobble at boundaries)
ADJACENT = 5              # nt of read allowed between a primary and a supplementary alignment (microhomology)
MIN_GAP = 20               # a split-read gap shorter than this is not counted as an intron
MIN_CUT_READS = 10         # poly(A) reads at a tRNA / gene start that move an overlapping stop there
RNA_INDEL_SHARE = 0.2      # share of reads with an indel at a position ...
SEQ_ERROR_SHARE = 0.8      # ... nearly all reads: the transcribed sequence differs from the genome (assembly error,
                           # or a polymorphism when the RNA is from another individual)
MIN_INDEL_DEPTH = 100      # ... at this depth or more
SINGLE_FRAME = 0.85        # an ORF in another frame covering this share of a frameshift-joined gene
SHORT_GENES = ("ATP8", "ND4L", "ND3")
CIGAR = mito_rna.CIGAR
COLUMNS = ["step", "feature", "verdict", "detail"]


class Tap:
    """Passes SAM lines on to mito_rna.analyse and tallies split-read junctions and indels on the way."""

    def __init__(self, lines: Iterable[str], n: int, ref: Optional[str] = None):
        self.lines, self.n, self.ref, self.n_sq = lines, n, ref, 0
        self.junctions: Counter = Counter()
        self.ins: Counter = Counter()          # (0-based position after which bases are inserted, length)
        self.dels: Counter = Counter()         # (0-based first deleted position, length)

    def __iter__(self) -> Iterator[str]:
        for ln in self.lines:
            if ln.startswith("@"):
                self.n_sq += ln.startswith("@SQ")
            else:
                self._read(ln)
            yield ln

    def _read(self, ln: str) -> None:
        c = ln.split("\t")
        flag = int(c[1])
        if flag & 0x904 or int(c[4]) < mito_rna.MIN_MAPQ:
            return
        if self.ref is not None and self.n_sq != 1 and c[2] != self.ref:
            return                                   # another reference sequence in the BAM: not this genome
        pos, ops = int(c[3]) - 1, CIGAR.findall(c[5])
        ref = pos
        for k, op in ops:
            k = int(k)
            if op == "I":
                self.ins[((ref - 1) % self.n, k)] += 1
            elif op == "D":
                self.dels[(ref % self.n, k)] += 1
            elif op == "N" and k > MIN_GAP:             # ordinary spliced alignment
                self.junctions[((ref - 1) % self.n, (ref + k) % self.n)] += 1
            if op in "MD=XN":
                ref += k
        sa = next((t for t in c[11:] if t.startswith("SA:Z:")), None)
        if not sa:
            return
        strand = "-" if flag & 16 else "+"
        lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
        qend = lead + sum(int(k) for k, op in ops if op in "MI=X")      # query interval of this alignment
        for s in sa[5:].strip().split(";"):
            if not s:
                continue
            _, sp, sst, scig, mq, _ = s.split(",")[:6]
            if s.split(",")[0] != c[2] or sst != strand or int(mq) < mito_rna.MIN_MAPQ:
                continue                             # other contig / other strand / weak
            sp = int(sp) - 1
            sops = CIGAR.findall(scig)
            slead = int(sops[0][0]) if sops and sops[0][1] in "SH" else 0
            sqend = slead + sum(int(k) for k, op in sops if op in "MI=X")
            slen = sum(int(k) for k, op in sops if op in "MD=X")
            # the two alignments must be neighbours on the read (no unexplained read bases between them)
            if lead < slead:
                if abs(qend - slead) > ADJACENT:
                    continue
                a, b = ref - 1, sp
            else:
                if abs(sqend - lead) > ADJACENT:
                    continue
                a, b = sp + slen - 1, pos                # last base before / first after the gap
            if b - a - 1 > MIN_GAP:
                self.junctions[(a % self.n, b % self.n)] += 1      # ends past the end of a circular genome wrap


def merge_junctions(j: Counter) -> List[tuple]:
    """(last exon-1 base, first exon-2 base, reads) - 0-based; ends within JUNCTION_MERGE nt joined at the most used."""
    out: List[list] = []
    for (a, b), c in sorted(j.items(), key=lambda t: -t[1]):
        for o in out:
            if abs(o[0] - a) <= JUNCTION_MERGE and abs(o[1] - b) <= JUNCTION_MERGE:
                o[2] += c
                break
        else:
            out.append([a, b, c])
    return [tuple(o) for o in out]


def _name(f) -> str:
    return f.qualifiers.get("gene", f.qualifiers.get("product", ["?"]))[0]


def step_introns(rec: SeqRecord, tap: Tap, depth, rows: list) -> None:
    junc = merge_junctions(tap.junctions)
    used = set()
    seq = str(rec.seq).upper()
    for f in rec.features:
        if f.type != "CDS" or len(f.location.parts) < 2:
            continue
        name = _name(f)
        parts = list(f.location.parts)                                 # transcript order
        if f.location.strand == -1:
            parts = parts[::-1]
        for p1, p2 in zip(parts, parts[1:]):
            if int(p1.end) == len(rec.seq) and int(p2.start) == 0:       # crosses the origin: not an intron
                continue
            a, b = int(p1.end) - 1, int(p2.start)                      # last base of the upstream part, first after
            gap = b - a - 1
            if gap < MIN_GAP:                                          # a frameshift join (MitoM P9), not an intron
                skipped = [x for x in range(a + 1, b)]
                dl = sum(c for (p, k), c in tap.dels.items() if p in skipped)
                d = int(depth[skipped[0] % len(depth)]) if skipped else 0
                verdict, alt = "flagged", single_frame_orf(rec, f)
                rows.append({"step": 1, "feature": name, "verdict": verdict, "detail":
                             f"frameshift join skipping {gap} nt at {a + 2}: the base is in {d - dl} of {d} RNA reads "
                             f"(not removed from the RNA); " + (f"an uninterrupted ORF in another frame covers it: "
                             f"{alt[0] + 1}..{alt[1]} ({alt[2]} aa) - check the reading frame" if alt else
                             "no single-frame ORF covers it: a real frameshift or an error in the sequence")})
                continue
            hit = [j for j in junc if abs(j[0] - a) <= JUNCTION_MERGE and abs(j[1] - b) <= JUNCTION_MERGE]
            cov = min(int(depth[a % len(depth)]), int(depth[b % len(depth)]))
            if hit and hit[0][2] >= MIN_JUNCTION:
                used.add(hit[0][:2])
                exact = hit[0][0] == a and hit[0][1] == b
                rows.append({"step": 1, "feature": f"{name} intron {a + 2}..{b}", "verdict": "confirmed", "detail":
                             f"{hit[0][2]} split reads skip it" + ("" if exact else
                             f" (ends {hit[0][0] + 2}..{hit[0][1]}, check the boundaries)")})
            else:
                rows.append({"step": 1, "feature": f"{name} intron {a + 2}..{b}",
                             "verdict": "flagged" if cov >= 50 else "not_testable",
                             "detail": f"{hit[0][2] if hit else 0} split reads at {cov}x exon depth"})
    for a, b, c in junc:
        flank = min(int(depth[a % len(depth)]), int(depth[b % len(depth)]))
        if c >= MIN_JUNCTION and (a, b) not in used and b > a and flank >= MIN_FLANK_DEPTH \
                and c >= JUNCTION_RATIO * flank:
            inside = [_name(f) for f in rec.features if f.type in ("CDS", "rRNA")
                      and int(f.location.start) <= a < int(f.location.end)]
            rows.append({"step": 1, "feature": f"junction {a + 2}..{b}", "verdict": "flagged", "detail":
                         f"{c} split reads ({c / flank:.0%} of the flank depth) skip {b - a - 1} nt, not annotated "
                         "as an intron"
                         + (f" (starts inside {', '.join(inside)})" if inside else "")
                         + f"; 5' {seq[a + 1:a + 7]}... 3' ...{seq[b - 4:b]}"})


def single_frame_orf(rec: SeqRecord, f) -> Optional[tuple]:
    """(start, end, aa) of an uninterrupted ORF in any frame covering >= SINGLE_FRAME of the span of f, start codon
    within 90 nt upstream of it, stop within 30 nt of its end. Coordinates are 0-based half-open on the forward
    genome for both strands (a minus-strand gene is searched on the reverse complement and mapped back, so start < end
    and 'upstream' is on the gene's own strand)."""
    code = int(f.qualifiers.get("transl_table", ["5"])[0])
    table = CodonTable.unambiguous_dna_by_id[code]
    starts, stops = set(table.start_codons) | {"ATG"}, set(table.stop_codons)
    g = str(rec.seq).upper()
    s0, e0 = int(f.location.start), int(f.location.end)
    minus = f.location.strand == -1
    if minus:
        g = str(Seq(g).reverse_complement())
        s0, e0 = len(g) - e0, len(g) - s0
    best = None
    for st in range(max(0, s0 - 90), s0 + 31):
        if g[st:st + 3] not in starts:
            continue
        e = next((x + 3 for x in range(st, min(len(g) - 2, e0 + 30), 3) if g[x:x + 3] in stops), None)
        if e is None or abs(e - e0) > 30:
            continue
        cover = (min(e, e0) - max(st, s0)) / (e0 - s0)
        if cover >= SINGLE_FRAME and (best is None or abs(st - s0) < abs(best[0] - s0)):
            best = (st, e, (e - st) // 3 - 1)
    if best and minus:
        best = (len(g) - best[1], len(g) - best[0], best[2])
    return best


def step_overlap_cut(rec: SeqRecord, sites: List[dict], stranded: bool, run: str, rows: list) -> None:
    """Complete stop codon overlapping the next same-strand tRNA / gene: end at the cut when the poly(A) says so."""
    src = f" ({run})" if run else ""
    feats = [f for f in rec.features if f.type in ("CDS", "tRNA", "rRNA") and f.location is not None
             and len(f.location.parts) == 1]
    for f in [x for x in feats if x.type == "CDS"]:
        st = f.location.strand or 1
        s0, e0 = int(f.location.start), int(f.location.end)
        code = int(f.qualifiers.get("transl_table", ["5"])[0])
        nt = str(f.location.extract(rec.seq)).upper()
        if len(nt) % 3 or nt[-3:] not in set(CodonTable.unambiguous_dna_by_id[code].stop_codons):
            continue
        nxt = [x for x in feats if x is not f and (x.location.strand or 1) == st and
               (s0 < int(x.location.start) < e0 if st == 1 else s0 < int(x.location.end) < e0)]
        if not nxt:
            continue
        x = min(nxt, key=lambda y: int(y.location.start)) if st == 1 else max(nxt, key=lambda y: int(y.location.end))
        ov = e0 - int(x.location.start) if st == 1 else int(x.location.end) - s0
        k = len(nt) - ov
        if ov > 30 or ov < 3 or k < 3:     # (overlap < 3: the cut and the full stop cannot be told apart; step 3)
            continue
        cut1 = s0 + k if st == 1 else e0 - k + 1                           # 1-based last base of the cut mRNA
        end1 = e0 if st == 1 else s0 + 1
        sign = "+" if st == 1 else "-"
        ok = [s for s in sites if not s.get("internal_priming") and s["strand"] == sign]
        near = [s for s in ok if abs(int(s["position"]) - cut1) <= 1]
        at_cut = sum(int(s["reads"]) for s in near)
        if near:                                       # the observed site, not the computed cut, ends the mRNA
            obs = int(max(near, key=lambda s: int(s["reads"]))["position"])
            k += (obs - cut1) * st
        if not abbreviated_stop(nt[:k], code):
            continue
        at_end = sum(int(s["reads"]) for s in ok if 0 <= (int(s["position"]) - end1) * st <= 2)
        name = _name(f)
        if at_cut >= MIN_CUT_READS and at_cut > at_end:
            loc = FeatureLocation(s0, s0 + k, st) if st == 1 else FeatureLocation(e0 - k, e0, st)
            tail = k % 3
            f.location = loc
            for gf in rec.features:
                if gf.type == "gene" and _name(gf) == name and int(gf.location.start) in (s0, int(loc.start)):
                    gf.location = loc
            f.qualifiers["transl_except"] = [transl_except_term(int(loc.start) + 1, int(loc.end), st, tail)]
            f.qualifiers["translation"] = [codon_translate(nt[:k], code).rstrip("*")]
            f.qualifiers["note"] = [n for n in f.qualifiers.get("note", []) if n != POLYA_NOTE] + [
                POLYA_NOTE, f"3' end at the start of {_name(x)} (mRNA cut there): {at_cut} poly(A) RNA-seq reads{src} "
                f"against {at_end} at the full stop codon (was ..{e0 if st == 1 else s0 + 1})"]
            rows.append({"step": 2, "feature": name, "verdict": "corrected", "detail":
                         f"stop overlapping {_name(x)} by {ov} nt: {at_cut} poly(A) reads at the cut, {at_end} at the "
                         f"full stop; ends with {nt[:k][-tail:]}"})
        elif at_end >= MIN_CUT_READS and at_end >= at_cut:
            rows.append({"step": 2, "feature": name, "verdict": "confirmed", "detail":
                         f"stop overlapping {_name(x)} by {ov} nt kept: {at_end} poly(A) reads at the full stop, "
                         f"{at_cut} at the cut"})


def step_rna_indels(rec: SeqRecord, tap: Tap, depth, rows: list) -> None:
    seq = str(rec.seq).upper()
    genes = [(int(f.location.start), int(f.location.end), _name(f),
              [(int(q.start), int(q.end)) for q in f.location.parts] if f.type == "CDS" else None)
             for f in rec.features if f.type in ("CDS", "rRNA") and f.location is not None]
    for kind, tally in (("insertion", tap.ins), ("deletion", tap.dels)):
        per, alleles = defaultdict(int), defaultdict(dict)
        for (p, k), c in tally.items():
            per[p] += c
            alleles[p][k] = alleles[p].get(k, 0) + c
        for p, c in per.items():
            al = alleles[p]
            d = int(depth[p])
            if d < MIN_INDEL_DEPTH or c < RNA_INDEL_SHARE * d:
                continue
            inside = [(n, parts) for s, e, n, parts in genes if s <= p < e]
            if not inside:
                continue
            coding = [n for n, parts in inside if parts and any(a <= p < b for a, b in parts)]
            where = (", ".join(coding) if coding else
                     ", ".join(f"{n} intron" if parts else n for n, parts in inside))
            run = re.match(r"(.)\1*", seq[p + 1:]) if kind == "insertion" else re.match(r"(.)\1*", seq[p:])
            hp = f"{run.group(1)}{len(run.group(0))}" if run else ""
            shift = any(k % 3 for k in al)               # any allele that moves the reading frame
            top = ", ".join(f"{k} nt: {n}" for k, n in sorted(al.items(), key=lambda t: -t[1])[:3])
            effect = (("a frameshift in those transcripts" if shift else
                       "a multiple of 3 nt (alleles " + top + "): in frame, adds / removes whole codons") if coding else
                      "not in a coding part (intron / rRNA): no effect on a protein")
            if shift and coding and len(al) > 1:
                effect += f" (alleles {top})"
            if c >= SEQ_ERROR_SHARE * d:
                effect = ("nearly all reads: the transcribed sequence differs from the genome here - a base error in "
                          "the assembly (same specimen) or a polymorphism (another individual)"
                          + ("; inside a protein gene, so the genome frame is suspect" if coding and shift else ""))
            rows.append({"step": 4, "feature": where, "verdict": "flagged", "detail":
                         f"RNA {kind} at {p + 1} in {c} of {d} reads ({c / d:.0%})"
                         + (f", homopolymer {hp}" if hp else "")
                         + f": {effect}; the genome sequence is kept (check DNA reads of the specimen: "
                           + ("an error if the DNA reads agree with the RNA)" if c >= SEQ_ERROR_SHARE * d else
                              "same length in DNA = slippage / RNA-level change)")})


def step_expression(features: List[dict], med: float, rows: list) -> None:
    short = [float(r["relative_depth"]) for r in features if r["type"] == "CDS" and r["relative_depth"]
             and gene_name({"gene": [r["feature"]]}, "CDS") in {g.lower() for g in SHORT_GENES}]
    floor = 0.5 * min(short) if short else None
    for r in features:
        if not r["feature"].lower().startswith("orf"):
            continue
        if r.get("duplicated_share") and float(r["duplicated_share"]) > mito_rna.DUPLICATED:
            rows.append({"step": 5, "feature": r["feature"], "verdict": "not_testable",
                         "detail": "sequence occurs twice in the record"})
        elif med < MIN_DEPTH or floor is None or not r["relative_depth"]:
            rows.append({"step": 5, "feature": r["feature"], "verdict": "not_testable", "detail": "too few reads"})
        else:
            rd = float(r["relative_depth"])
            rows.append({"step": 5, "feature": r["feature"], "verdict": "confirmed" if rd >= floor else "flagged",
                         "detail": f"depth {rd:.2f} x median, short genes' floor {floor:.2f}; "
                                   f"{float(r['covered_share']):.0%} of its length transcribed"})


def annotate(rec: SeqRecord, sam: Iterable[str], run: str = "") -> tuple:
    """-> (revised record, evidence rows, info). The record is changed in place."""
    tap = Tap(sam, len(rec.seq), rec.id)
    features, sites, info = mito_rna.analyse(rec, tap, rec.id)
    depth = info["depth"][1] + info["depth"][-1]
    rows: List[dict] = []
    er, med = info["edit_rate"], info["pcg_median_depth"]
    rows.append({"step": 0, "feature": rec.id, "verdict": "flagged" if er > DIVERGED else "confirmed",
                 "detail": f"edit rate {100 * er:.2f} %, {info['reads_used']} reads, protein-gene median depth "
                           f"{med:.0f}, {info['polya_reads']} poly(A) reads"})
    if er > DIVERGED:
        rows[0]["detail"] += ": the reads are not this species (contamination / mislabelled sample) - nothing changed"
        return rec, rows, info
    for f in features:
        if f.get("duplicated_share") and float(f["duplicated_share"]) > mito_rna.DUPLICATED and f["type"] == "CDS":
            rows.append({"step": 0, "feature": f["feature"], "verdict": "not_testable",
                         "detail": f"{float(f['duplicated_share']):.0%} of its sequence occurs elsewhere in the "
                                   "record"})
    step_introns(rec, tap, depth, rows)
    polya_ok = info["polya_reads"] >= MIN_POLYA_TOTAL
    stranded = info["stranded"]
    if polya_ok:
        step_overlap_cut(rec, sites, stranded, run, rows)
        for c in mito_rna.refine(rec, features, sites, stranded, run):
            rows.append({"step": 3, "feature": c["gene"], "verdict": "corrected",
                         "detail": f"{c['kind']}: {c['old_start']}..{c['old_end']} -> {c['new_start']}..{c['new_end']} "
                                   f"({c['site_reads']} poly(A) reads)"})
        conf = sum(any("confirmed by" in n for n in f.qualifiers.get("note", [])) for f in rec.features
                   if f.type in ("CDS", "rRNA"))
        rows.append({"step": 3, "feature": "3' ends", "verdict": "confirmed",
                     "detail": f"{conf} protein-gene / rRNA 3' ends confirmed by poly(A) reads"})
    else:
        rows.append({"step": 3, "feature": "3' ends", "verdict": "not_testable",
                     "detail": f"{info['polya_reads']} poly(A) reads: no 3' end evidence in this library"})
    step_rna_indels(rec, tap, depth, rows)
    step_expression(features, med, rows)
    return rec, rows, info


def main(argv: Sequence[str] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-g", "--genbank", required=True, help="annotated mitogenome (one record)")
    p.add_argument("-1", "--reads1", default=None)
    p.add_argument("-2", "--reads2", default=None)
    p.add_argument("--sam", default=None, help="SAM of reads already mapped to this record (bwa mem, with SA tags)")
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--run", default="", help="run accession quoted in the notes")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--bin", default=None, help="folder with bwa / samtools")
    a = p.parse_args(argv)
    rec = SeqIO.read(a.genbank, "genbank")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.sam:
        rec, rows, info = annotate(rec, sys.stdin if a.sam == "-" else open(a.sam), a.run)
    else:
        if not a.reads1:
            raise SystemExit("give reads (-1 [-2]) or --sam")
        ref = out / "ref.fa"
        ref.write_text(f">{rec.id}\n{rec.seq}\n")
        proc = mito_rna.map_reads(ref, a.reads1, a.reads2, a.threads, a.bin, "bwa")
        rec, rows, info = annotate(rec, proc.stdout, a.run)
        proc.wait()
    write_tsv(out / "rnaseq_evidence.tsv", COLUMNS, rows)
    write_gb(rec, out / "rnaseq_annotated.gb")
    n = Counter(r["verdict"] for r in rows)
    print(f"[mito_rnaseq_annot] {rec.id}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items()))
          + f" -> {out}/rnaseq_evidence.tsv, rnaseq_annotated.gb", file=sys.stderr)
    for r in rows:
        if r["verdict"] in ("corrected", "flagged"):
            print(f"  [{r['step']}] {r['verdict']:9s} {r['feature']}: {r['detail']}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(cli_guard(main))
