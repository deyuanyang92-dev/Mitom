#!/usr/bin/env python3
"""mito_orf.py - open reading frames of unknown function in the sequence no gene exon covers.

Many animal mitogenomes hold protein genes beyond the 13: atp9 of sponges, mutS of octocorals, F-ORF / M-ORF of
freshwater mussels, homing endonucleases and maturases inside introns, a conserved ~355-codon ORF between trnQ and
trnA of the spionid polychaetes Polydora and Boccardiella (NCBI leaves it unannotated). An ORF is a gene when
selection has kept it open. Length alone cannot tell: AT-rich sequence has few stops, so ORFs of 75 codons and more
are common by chance (an ORF finder at that length writes 357 ORFs in 362 intron-free RefSeq mitogenomes), and short
real proteins are no longer. The evidence, from first principles:

1. Length beyond chance: under the base composition of the record's unannotated sequence, the expected number of
   ORFs at least this long from a start codon of the code, E_len = 2 F p_start (1 - p_stop)^L (F nt; L codons).
2. Kept open: the ORF has homologs (blastp E <= 1e-5, low-complexity masked) that are ORFs themselves - the other
   records of the run, the ORF library of 16,073 RefSeq animal mitogenomes (models/orf_library.fna.gz,
   scripts/build_orf_library.py), --ref relatives. Up to 20 are placed on the ORF; the substitutions since they split
   are at least the larger of (a) one per extra base kind at each position and (b) half the minimum spanning tree of
   their differences (any history with ancestors is at least that long). In non-coding sequence a substitution makes
   a stop with probability f (single-nucleotide changes of the ORF's own codons), so all stayed open by chance with
   p_open = (1 - f)^substitutions. A close relative proves nothing; many relatives add up.
3. Coding frame: in local nucleotide alignments with the homologs (which do not know the frame; the protein alignment
   would favour differences at position 3), a coding frame changes most at codon position 3 and least at 2 (every
   change there alters the amino acid): chi-square against 1/3 each, p_frame, counted only in that order (else 1).
   The frame on the other strand whose third positions coincide changes least at its position 1, so it fails.
4. Tandem repeats (copies >= 85 % identical, period 3-150 nt) are masked first: an array of a stop-free unit is an
   ORF of any length and aligns with any other array (control regions).
5. Overlapping ORFs (other frames, other strand) keep the one whose position 2 changes least, then the larger share
   at position 3, then the stronger evidence; the rest are its shadows.

Classes (orf_report.tsv):
    conserved          p_open <= 0.05 and p_frame <= 0.05 (each necessary; both by chance ~0.0025)
                                                                    -> CDS orf<codons>, "hypothetical protein"
    known_family       homolog in a known family: models/mfannot_proteins (not the 13; atp9 gets its name) or a
                       Pfam family above its gathering threshold (models/orf_families.hmm.gz: homing endonucleases,
                       group II intron reverse transcriptase / maturase, polB, MutS; the product names it) -> CDS
    long_unconfirmed   E_len <= 0.01 without such evidence           -> a WARNING, NCR /note (--orf all: written)
    standard_fragment  homolog of one of the 13 protein genes        -> a WARNING (second copy, pseudogene, or a
                                                                       piece the annotation missed)
ORFs within chance length without such homologs are not reported (chance ORFs are in every AT-rich stretch).
Introns are open to the search (non-coding to their host gene). A written ORF runs from the first start codon of the
unannotated stretch to its stop, or - when the frame runs on into the next gene - to the base before it with an
incomplete stop T / TA (/transl_except); a frame running into a tRNA otherwise ends at its last whole codon before
the tRNA with an open 3' end (tRNA punctuation: the transcript is cut there; Selifanova et al. 2023 defined the
Polydora ORFans so); a frame running into a protein gene or rRNA is a WARNING naming it ("in frame": perhaps its start
lies further upstream). Records of an earlier run: their mito_orf genes are replaced.
Measured: benchmark/results/orf_2026-11-09/README.md.

Usage: python3 mito_orf.py -i <GenBank file or folder> -o <out folder> [--genetic_code 5] [--orf yes|all|report]
       [--ref relatives.gb]  (`mitom run` / `annotate` / `import`: after curation, before the non-coding regions)
"""
from __future__ import annotations

import argparse
import functools
import gzip
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from mito_common import cli_guard, iter_gb_files, require_path, run_tool, stop_qualifiers, write_gb, write_tsv  # noqa: E402
from mito_names import PCG, PRODUCTS, gene_name  # noqa: E402
from mito_ncr import extent  # noqa: E402

HERE = Path(__file__).resolve().parent
LIBRARY = HERE / "models" / "orf_library.fna.gz"
FAMILIES = HERE / "models" / "mfannot_proteins"
PFAM = HERE / "models" / "orf_families.hmm.gz"     # Pfam (CC0) HMMs of the extra mitochondrial protein families
PFAM_PRODUCT = {                                  # Pfam family -> product of a CDS written for it
    "LAGLIDADG_1": "LAGLIDADG homing endonuclease", "LAGLIDADG_2": "LAGLIDADG homing endonuclease",
    "GIY-YIG": "GIY-YIG homing endonuclease", "HNH": "HNH endonuclease",
    "RVT_1": "group II intron reverse transcriptase/maturase", "Intron_maturas2": "group II intron reverse "
    "transcriptase/maturase", "GIIM": "group II intron reverse transcriptase/maturase",
    "DNA_pol_B": "DNA polymerase family B", "DNA_pol_B_2": "DNA polymerase family B",
    "DNA_pol_B_exo1": "DNA polymerase family B",
    "MutS_I": "MutS-like protein", "MutS_II": "MutS-like protein", "MutS_III": "MutS-like protein",
    "MutS_V": "MutS-like protein"}
MIN_AA = 40                 # shortest ORF considered (codons): below this a homology search has no power
E_LEN = 0.01                # length beyond chance: expected chance ORFs this long per record
P_OPEN = 0.01               # kept open by selection: chance that both copies stayed open (alone)
P_EACH = 0.05               # kept open and coding pattern: each must hold (both by chance in non-coding: ~0.0025)
MAX_HOMOLOGS = 20           # homologs (best blastp E, distinct) joined in the selection tests
E_BLAST = 1e-5              # homolog: blastp E-value
REPEAT_ID = 0.85            # tandem repeat: copies >= 85 % identical ...
REPEAT_COPIES = 2           # ... over >= 2 periods (period 3-150 nt)
FREE = 0.9                  # share of a written ORF's codons in unannotated sequence (overlaps with a tRNA end)
NOTE = "found by MitoM mito_orf"      # start of the /note of every gene this module writes
EXON_TYPES = {"CDS", "rRNA", "tRNA"}
STANDARD = set(PCG)
COLUMNS = ["record", "orf", "class", "written", "start", "end", "strand", "codons", "e_len", "homolog",
           "homolog_evalue", "homolog_identity", "substitutions", "p_open", "pos3_share", "p_pos3", "family", "evidence"]


def gene_key(f: SeqFeature) -> str:
    """The standard name (cox1, nad4l ...) of a CDS, '' when it is not one of the 13."""
    return gene_name(f.qualifiers, "CDS")


def is_orf(f: SeqFeature) -> bool:
    """A CDS or gene written by this module."""
    return any(x.startswith(NOTE) for x in f.qualifiers.get("note", [])) or \
        (f.type == "gene" and f.qualifiers.get("note", [""])[0].startswith(NOTE))


def code_of(rec: SeqRecord, default: int) -> int:
    codes = [c.qualifiers.get("transl_table", [""])[0] for c in rec.features if c.type == "CDS"]
    codes = [int(c) for c in codes if str(c).isdigit()]
    return max(set(codes), key=codes.count) if codes else default


def covered_mask(rec: SeqRecord) -> np.ndarray:
    """True where an exon of a CDS, rRNA or tRNA lies. Introns are free: to their host gene they are non-coding, and
    they hold ORFs of their own (homing endonucleases of group I introns, maturases of group II; the mutS of
    octocorals). A gene feature covers only when no CDS / RNA feature gives its exons. This module's own ORFs and
    non-coding / control-region features do not cover."""
    n = len(rec.seq)
    cov = np.zeros(n, dtype=bool)
    exons = [f for f in rec.features if f.type in EXON_TYPES and f.location is not None and not is_orf(f)]
    named = {f.qualifiers.get("gene", [""])[0] for f in exons}
    for f in exons + [f for f in rec.features if f.type == "gene" and f.location is not None and not is_orf(f)
                      and f.qualifiers.get("gene", [""])[0] not in named and len(f.location.parts) == 1
                      and len(f) < 6000]:
        for part in f.location.parts:
            s, e = int(part.start), int(part.end)
            cov[np.arange(s, e if e >= s else e + n) % n] = True
    return cov


def uncovered(rec: SeqRecord) -> Callable[[int, int], float]:
    """free(s, e): share of positions s..e-1 (0-based, e may pass the origin) no gene covers."""
    cov, n = covered_mask(rec), len(rec.seq)
    return lambda s, e: float(1 - cov[np.arange(s, e) % n].mean()) if e > s else 0.0


def orfs_stop_to_stop(seq: str, code: int, circular: bool, min_aa: int = MIN_AA) -> Iterable[Tuple[int, int, int, str]]:
    """(start, end, strand, protein) of every stretch between two stop codons of >= min_aa codons in all six frames;
    start / end 0-based on the plus strand, end after the closing stop and > len(seq) across the origin. The open
    end of a linear record is not an ORF."""
    n = len(seq)
    s2 = seq + seq[: n - 1] if circular else seq
    m = len(s2)
    for strand in (1, -1):
        t = s2 if strand == 1 else str(Seq(s2).reverse_complement())
        for fr in range(3):
            sub = t[fr: m - (m - fr) % 3]
            prot = str(Seq(sub).translate(table=code))
            stops = [i for i, a in enumerate(prot) if a == "*"]
            for a, b in zip(stops, stops[1:]):
                if b - a - 1 < min_aa:
                    continue
                p0, p1 = fr + 3 * (a + 1), fr + 3 * (b + 1)            # codons after the stop .. closing stop
                s, e = (p0, p1) if strand == 1 else (m - p1, m - p0)
                if s >= n or e - s >= n:                               # (copy across the origin, or > genome)
                    continue
                yield s, e, strand, prot[a + 1: b]


def tandem_mask(nt: str, max_period: int = 150) -> np.ndarray:
    """True at the bases of tandem repeats in nt: a stretch of >= REPEAT_COPIES + 1 copies of a 3-150 nt unit with
    >= REPEAT_ID identity between neighbouring copies. Replication slippage builds such arrays in control regions; an
    array of a stop-free unit is an ORF of any length and aligns with any other such array, so it is evidence of
    neither length beyond chance nor selection."""
    a = np.frombuffer(nt.encode(), dtype=np.uint8)
    n = len(a)
    mask = np.zeros(n, dtype=bool)
    for p in range(3, min(max_period, n // (REPEAT_COPIES + 1)) + 1):
        eq = (a[:-p] == a[p:]).astype(np.int32)
        w = REPEAT_COPIES * p
        if len(eq) < w:
            break
        c = np.concatenate([[0], np.cumsum(eq)])
        win = (c[w:] - c[:-w]) >= REPEAT_ID * w                  # window i..i+w: i matches i+p
        for i in np.flatnonzero(win):
            mask[i: i + w + p] = True
    return mask


def composition(rec: SeqRecord, cov: np.ndarray) -> Dict[str, float]:
    s = np.frombuffer(str(rec.seq).upper().encode(), dtype="S1")[~cov]
    tot = sum(int((s == x.encode()).sum()) for x in "ACGT") or 1
    return {x: int((s == x.encode()).sum()) / tot for x in "ACGT"}


def codon_prob(codons: Iterable[str], freq: Dict[str, float]) -> float:
    return sum(freq[c[0]] * freq[c[1]] * freq[c[2]] for c in codons)


def stop_share(nt: str, code: int) -> float:
    """f: the share of single-nucleotide changes of the codons of nt that give a stop codon."""
    stops = set(CodonTable.unambiguous_dna_by_id[code].stop_codons)
    k = hit = 0
    for i in range(0, len(nt) - 2, 3):
        c = nt[i: i + 3]
        for j in range(3):
            for b in "ACGT":
                if b != c[j]:
                    k += 1
                    hit += (c[:j] + b + c[j + 1:]) in stops
    return hit / k if k else 0.0


@dataclass
class Orf:
    record: str
    start: int          # 0-based plus-strand start of the written ORF (first start codon), end after its last base
    end: int
    strand: int
    piece: str          # protein of the longest unannotated stretch of the reading frame (homology search)
    piece_nt: str       # its nucleotides (reading order)
    prot: str           # protein from the first start codon (written CDS)
    nt: str             # nucleotides of the written ORF (reading order; with the stop, or an incomplete T / TA one)
    e_len: float
    writable: bool
    hits: List[dict] = field(default_factory=list)
    cls: str = ""
    best: Optional[dict] = None
    p_open: float = 1.0
    family: str = ""
    written: bool = False
    shadow_of: str = ""
    runs_into: str = ""     # not writable: the features the reading frame runs on into (in frame with a CDS: named)
    trna_end: str = ""      # ends without a stop at this tRNA (tRNA punctuation): written with an open 3' end

    @property
    def codons(self) -> int:
        return len(self.prot)

    @property
    def name(self) -> str:
        return f"orf{self.codons}"


def candidates(rec: SeqRecord, code: int, min_aa: int = MIN_AA) -> List[Orf]:
    """ORFs of >= min_aa codons in the record's unannotated sequence, with their length E-value. The written ORF runs
    from the first start codon of the unannotated stretch to its stop codon, or - when the reading frame runs on into
    the next gene - to the base before that gene with an incomplete stop T / TA (as mitochondrial genes end before a
    tRNA); one ending otherwise is reported, not written."""
    n = len(rec.seq)
    seq = str(rec.seq).upper()
    circular = rec.annotations.get("topology") == "circular"
    cov = covered_mask(rec)
    if cov.all():
        return []
    freq = composition(rec, cov)
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    rfreq = {comp[k]: v for k, v in freq.items()}
    table = CodonTable.unambiguous_dna_by_id[code]
    starts = set(table.start_codons)
    F = int((~cov).sum())
    trna_at: Dict[int, str] = {}                            # position -> tRNA name (where an ORF may end, no stop)
    for f in rec.features:
        if f.type == "tRNA" and f.location is not None:
            name = f.qualifiers.get("gene", f.qualifiers.get("product", ["tRNA"]))[0]
            for part in f.location.parts:
                for i in range(int(part.start), int(part.end)):
                    trna_at[i % n] = name
    out = []
    for s, e, strand, prot in orfs_stop_to_stop(seq, code, circular, min_aa):
        k = len(prot)                                       # codons before the closing stop
        free = ~cov[np.arange(s, e) % n]
        if strand == -1:
            free = free[::-1]                               # reading order
        nt_all = seq[s:] + seq[: e - n] if e > n else seq[s:e]
        if strand == -1:
            nt_all = str(Seq(nt_all).reverse_complement())
        cod = free[: 3 * k].reshape(k, 3).all(1)
        best, run = (0, 0), 0                               # longest run of unannotated codons
        for i, ok in enumerate(cod):
            run = run + 1 if ok else 0
            if run > best[1] - best[0]:
                best = (i - run + 1, i + 1)
        if best[1] - best[0] < min_aa:
            continue
        first = next((i for i in range(best[0], best[1]) if nt_all[3 * i: 3 * i + 3] in starts), None)
        if first is None:
            continue
        a = 3 * first
        trna_end = ""
        if best[1] == k:                                    # runs to its own stop codon (which may touch a gene)
            b, writable = 3 * k + 3, True
        else:                                               # runs into a gene: end before it, incomplete stop
            b = a + int(np.argmin(free[a:])) if not free[a:].all() else 3 * k
            tail = nt_all[a:b][(b - a) // 3 * 3:]
            writable = tail in ("T", "TA")
            nxt = trna_at.get((s + b) % n if strand == 1 else (e - 1 - b) % n)
            if not writable and nxt:                        # tRNA punctuation: the transcript is cut at the tRNA
                b, writable, trna_end = a + (b - a) // 3 * 3, True, nxt
        L = (b - a) // 3 - (1 if b == 3 * k + 3 else 0)
        rep = tandem_mask(nt_all[: 3 * k])
        rcod = rep.reshape(k, 3).any(1)                    # codons in a tandem repeat
        piece = "".join("X" if rcod[i] else prot[i] for i in range(best[0], best[1]))
        if len(piece) - piece.count("X") < min_aa:
            continue                                        # (an ORF made by a repeat array)
        L -= int(rcod[first: first + L].sum())              # length beyond chance: the codons outside repeats
        f = freq if strand == 1 else rfreq
        e_len = 2 * F * codon_prob(starts, f) * (1 - codon_prob(table.stop_codons, f)) ** L
        ws, we = (s + a, s + b) if strand == 1 else (e - b, e - a)
        o = Orf(rec.id, ws, we, strand, piece, nt_all[3 * best[0]: 3 * best[1]],
                "M" + prot[first + 1: (b - a) // 3 + first - (1 if b == 3 * k + 3 else 0)], nt_all[a:b], e_len,
                writable and L >= min_aa)
        o.trna_end = trna_end
        if not o.writable:
            rest = (s + b, e) if strand == 1 else (s, e - b)           # the frame beyond the ORF, to its stop
            o.runs_into = runs_into(rec, rest, strand, (s if strand == 1 else e) % 3, n)
        out.append(o)
    return out


def runs_into(rec: SeqRecord, rest: Tuple[int, int], strand: int, frame: int, n: int) -> str:
    """The features in `rest` (the reading frame beyond an ORF, up to its stop); a CDS on the same strand in the same
    frame is marked 'in frame' (the ORF may be a part of it: a start further upstream, or a missed exon)."""
    names = []
    for f in rec.features:
        if f.type not in ("CDS", "tRNA", "rRNA") or f.location is None:
            continue
        fs, fl = extent(f.location, n, rec.annotations.get("topology") == "circular")
        if not (((rest[0] - fs) % n) < fl or ((fs - rest[0]) % n) < rest[1] - rest[0]):     # (circular overlap)
            continue
        name = f.qualifiers.get("gene", f.qualifiers.get("product", [f.type]))[0]
        if f.type == "CDS" and f.location.strand == strand:
            edge = int(f.location.start) if strand == 1 else int(f.location.end)
            if edge % 3 == frame:
                name += " (in frame)"
        if name not in names:
            names.append(name)
    return ", ".join(names) if names else "no annotated feature"


@functools.lru_cache(maxsize=1)
def nt_aligner():
    from Bio import Align
    return Align.PairwiseAligner(mode="local", match_score=2, mismatch_score=-3, open_gap_score=-5,
                                 extend_gap_score=-2)


def nt_differences(o: "Orf", s_nt: str) -> List[int]:
    """Differences at codon positions 1, 2, 3 (of the ORF's frame) in a local nucleotide alignment of the ORF's
    unannotated stretch with its homolog (blastn-like scores). A nucleotide alignment does not know the frame, so -
    unlike the protein alignment that found the homolog, which favours stretches of identical amino acids and so
    differences at position 3 - it leaves non-coding differences spread evenly over the three positions. Codons of a
    tandem repeat (X in the piece) are left out."""
    q = o.piece_nt
    aln = nt_aligner().align(q, s_nt)[0]
    d = [0, 0, 0]
    for (qa, qb), (sa, sb) in zip(*aln.aligned):
        for k in range(qb - qa):
            i = qa + k
            if o.piece[i // 3] != "X" and q[i] != s_nt[sa + k]:
                d[i % 3] += 1
    return d


def combined_evidence(o: "Orf", homologs: List[dict]) -> dict:
    """The selection tests over the ORF and all its homologs together. Each homolog is placed on the ORF's codons
    (by its blastp HSP for the substitution count, by a local nucleotide alignment for the codon positions); at every
    position of the ORF, k different nucleotides among the copies mean >= k - 1 substitutions on any tree that joins
    them. The sum over positions is a lower bound of the substitutions since the copies split, so that more relatives,
    each of them close, add up; identical copies add nothing."""
    by_p: Dict[int, set] = {}
    by_n: Dict[int, set] = {}
    q = o.piece_nt
    placed: List[Dict[int, str]] = [dict(enumerate(q))]     # each copy's base at each ORF position (HSP)
    for h in homologs:
        qa, sa, qi, si = h["hsp"]
        s_nt = h["s_nt"]
        mine: Dict[int, str] = {}
        for x, y in zip(qa, sa):
            if x != "-" and y != "-" and x != "X":
                for j in range(3):
                    if 3 * si + j < len(s_nt) and 3 * qi + j < len(q):
                        by_p.setdefault(3 * qi + j, {q[3 * qi + j]}).add(s_nt[3 * si + j])
                        mine[3 * qi + j] = s_nt[3 * si + j]
            qi += x != "-"
            si += y != "-"
        placed.append(mine)
        aln = nt_aligner().align(q, s_nt)[0]
        for (a, b), (c, _) in zip(*aln.aligned):
            for k in range(b - a):
                i = a + k
                if o.piece[i // 3] != "X":
                    by_n.setdefault(i, {q[i]}).add(s_nt[c + k])
    d = [0, 0, 0]
    for i, v in by_n.items():
        d[i % 3] += len(v) - 1
    columns = sum(len(v) - 1 for v in by_p.values())
    best = homologs[0]
    return {"name": best["name"] + (f" and {len(homologs) - 1} more" if len(homologs) > 1 else ""),
            "src": best["src"], "evalue": best["evalue"], "identity": best["identity"], "n_homologs": len(homologs),
            "subs": max(columns, mst_bound(placed)), "nt_diffs": sum(d),
            "pos3": d[2] / sum(d) if sum(d) else 0.0, "p_pos3": third_position_p(d), "d": d}


def mst_bound(copies: Sequence[Dict[int, str]]) -> int:
    """Half the minimum spanning tree over the copies, distances = differences at the positions both hold: a lower
    bound of the substitutions of any history joining them (a tree with ancestors at least half its MST long; tree
    paths at least the differences). With many close relatives it is far tighter than one difference per base kind
    and position."""
    n = len(copies)
    if n < 2:
        return 0

    def dist(a: Dict[int, str], b: Dict[int, str]) -> int:
        return sum(1 for i, x in a.items() if i in b and b[i] != x)
    inside, best = {0}, {j: dist(copies[0], copies[j]) for j in range(1, n)}
    total = 0
    while best:                                            # Prim
        j = min(best, key=best.get)
        total += best.pop(j)
        inside.add(j)
        for k in best:
            best[k] = min(best[k], dist(copies[j], copies[k]))
    return total // 2


def codon_differences(q_aln: str, s_aln: str, q_nt: str, q_from: int, s_nt: str, s_from: int) -> int:
    """Nucleotide differences between the codons of two aligned proteins (the blastp HSP; gaps and masked repeat
    codons skipped): a lower bound of the substitutions since the two copies split - synonymous ones included, they
    too could have made a stop. (The HSP favours stretches of similar protein, so it undercounts: conservative.)"""
    qi, si, d = q_from, s_from, 0
    for x, y in zip(q_aln, s_aln):
        if x != "-" and y != "-" and x != "X":
            d += sum(c1 != c2 for c1, c2 in zip(q_nt[3 * qi: 3 * qi + 3], s_nt[3 * si: 3 * si + 3]))
        qi += x != "-"
        si += y != "-"
    return d


def third_position_p(d: Sequence[int]) -> float:
    """p of the coding pattern of the differences at codon positions 1, 2, 3 of a frame: a protein-coding frame changes
    most at its synonymous third positions and least at the second (every change there alters the amino acid);
    non-coding sequence read in any frame changes the three alike. Chi-square goodness of fit against 1/3 each (2 df:
    p = exp(-X/2)), counted only when position 3 is the most and position 2 the least changed (else 1). The frame on
    the other strand whose third positions coincide with the gene's changes least at its position 1, so it fails."""
    from math import exp
    n = sum(d)
    if not n or d[2] < max(d) or d[1] > min(d):
        return 1.0
    x = sum((k - n / 3) ** 2 / (n / 3) for k in d)
    return min(1.0, exp(-x / 2))


# ── homology ─────────────────────────────────────────────────────────────────────────────────────────────────────
def find_tool(name: str, bin_dir: Optional[str] = None) -> Optional[str]:
    """bin_dir, then the bin folder of this Python (the MitoM environment, also when it is not activated), then PATH."""
    for d in (bin_dir, str(Path(sys.executable).parent)):
        if d and (Path(d) / name).exists():
            return str(Path(d) / name)
    return shutil.which(name)


def read_library(library: Optional[Path], exclude: set) -> Dict[str, Tuple[str, str]]:
    """ORF library entries (models/orf_library.fna.gz: nucleotide ORFs in reading order, header
    accession|start|end|strand|kind|code|phylum|organism) -> {id: (protein, nucleotides)}; records in `exclude` left
    out (the record itself; MITOM_EXCLUDE_REF in benchmarks)."""
    out: Dict[str, Tuple[str, str]] = {}
    if not library or not library.exists():
        return out
    with gzip.open(library, "rt") as fh:
        for r in SeqIO.parse(fh, "fasta"):
            f = r.id.split("|")
            if f[0] in exclude or f[0].split(".")[0] in exclude:
                continue
            nt = str(r.seq)
            out[r.id] = (str(Seq(nt[: len(nt) // 3 * 3]).translate(table=int(f[5]))).replace("*", "X"), nt)
    return out


def search(batch: List[Orf], library: Optional[Path], exclude: set, threads: int = 1,
           bin_dir: Optional[str] = None, extra: Optional[Dict[str, Tuple[str, str]]] = None) -> None:
    """blastp of every ORF piece against the other records' pieces, the MFannot families and the ORF library; fills
    Orf.hits (with the nucleotide substitutions between the two copies where both are ORFs)."""
    blastp, makedb = find_tool("blastp", bin_dir), find_tool("makeblastdb", bin_dir)
    if not (blastp and makedb):
        raise RuntimeError("blastp / makeblastdb not found (conda install -c bioconda blast)")
    lib = {**read_library(library, exclude), **(extra or {})}
    with tempfile.TemporaryDirectory() as td:
        db, q = Path(td) / "db.faa", Path(td) / "q.faa"
        names: List[Tuple[str, str]] = []               # db sequence k -> (kind, key): BLAST rewrites '|' ids
        with open(db, "w") as fh:
            def put(kind: str, key: str, prot: str) -> None:
                fh.write(f">s{len(names)}\n{prot}\n")
                names.append((kind, key))
            for i, o in enumerate(batch):
                put("b", str(i), o.piece)
            for f in sorted(FAMILIES.glob("*.faa")):
                for r in SeqIO.parse(str(f), "fasta"):
                    put("f", f.stem, str(r.seq).replace("*", ""))
            for k, (prot, _) in lib.items():
                put("l", k, prot)
        q.write_text("".join(f">{i}\n{o.piece}\n" for i, o in enumerate(batch)))
        run_tool([makedb, "-in", str(db), "-dbtype", "prot", "-out", str(Path(td) / "db")], check=True)
        p = run_tool([blastp, "-query", str(q), "-db", str(Path(td) / "db"), "-evalue", str(E_BLAST),
                            "-seg", "yes", "-max_target_seqs", "200", "-num_threads", str(threads), "-outfmt",
                            "6 qseqid sseqid evalue pident length mismatch qstart sstart qseq sseq"],
                     check=True)
    seen = set()
    for ln in p.stdout.splitlines():
        q, s, ev, pid, length, mm, qs, ss, qa, sa = ln.split("\t")
        if (q, s) in seen:                                   # (best HSP only)
            continue
        seen.add((q, s))
        o = batch[int(q)]
        kind, rest = names[int(s[1:])]
        s_nt = ""
        if kind == "b":
            other = batch[int(rest)]
            if other.record == o.record:
                continue
            name, src, s_nt = f"{other.record}:{other.start + 1}..{other.end}", "run", other.piece_nt
        elif kind == "f":
            name, src = rest, "family"
        else:
            name, src, s_nt = rest, "library", lib[rest][1]
        h = {"name": name, "src": src, "evalue": float(ev), "identity": float(pid), "length": int(length),
             "mismatches": int(mm)}
        if s_nt:
            h["s_nt"], h["hsp"] = s_nt, (qa, sa, int(qs) - 1, int(ss) - 1)
        o.hits.append(h)
    for o in batch:                                          # the selection tests on all homologs together
        homologs, seen_names = [], set()
        for h in sorted((h for h in o.hits if "s_nt" in h), key=lambda h: h["evalue"]):
            if h["name"] not in seen_names and len(homologs) < MAX_HOMOLOGS:
                seen_names.add(h["name"])
                homologs.append(h)
        if homologs:
            o.hits.append(combined_evidence(o, homologs))
        for h in o.hits:
            h.pop("s_nt", None)
            h.pop("hsp", None)


def pfam_search(batch: List["Orf"], threads: int = 1, bin_dir: Optional[str] = None) -> None:
    """hmmsearch of the ORF pieces against the Pfam families of extra mitochondrial proteins (homing endonucleases,
    group II intron reverse transcriptase / maturase, family-B DNA polymerase, MutS) with Pfam's own curated
    gathering thresholds (--cut_ga); a hit is a homolog of src 'pfam'."""
    hmmsearch = find_tool("hmmsearch", bin_dir)
    if not batch:
        return
    if not hmmsearch or not PFAM.exists():
        print(f"[mito_orf] WARNING {'hmmsearch (HMMER) not found' if not hmmsearch else f'{PFAM} missing'}: ORFs "
              "of known Pfam families (homing endonucleases, intron reverse transcriptases, polB, MutS) are not "
              "named", file=sys.stderr)
        return
    with tempfile.TemporaryDirectory() as td:
        hmm, q, tbl = Path(td) / "f.hmm", Path(td) / "q.faa", Path(td) / "t.tbl"
        with gzip.open(PFAM, "rt") as src:
            hmm.write_text(src.read())
        q.write_text("".join(f">{i}\n{o.piece}\n" for i, o in enumerate(batch)))
        run_tool([hmmsearch, "--cut_ga", "--noali", "--cpu", str(threads), "--tblout", str(tbl), str(hmm), str(q)],
                 check=True)
        for ln in tbl.read_text().splitlines():
            if ln.startswith("#"):
                continue
            c = ln.split()
            o, fam, ev = batch[int(c[0])], c[2], float(c[4])
            o.hits.append({"name": fam, "src": "pfam", "evalue": ev, "identity": 0.0, "length": 0, "mismatches": 0})


STANDARD_FAMILIES = {g.lower() for g in STANDARD}


def fisher(p1: float, p2: float) -> float:
    """Fisher's combined p of two independent tests: P(chi2_4 >= -2 ln(p1 p2)) = p1 p2 (1 - ln(p1 p2))."""
    from math import log
    x = max(p1 * p2, 1e-300)
    return min(1.0, x * (1 - log(x)))


def frame_rank(o: "Orf") -> tuple:
    """Order of overlapping ORFs: the coding frame changes least at codon position 2 (every change there alters the
    amino acid) and most at 3; the frame on the other strand whose third positions coincide with them changes least
    at its position 1. So: position 2 <= position 1 first, then the larger share at position 3, then the evidence."""
    h = o.best or {}
    d = h.get("d") or [0, 0, 0]
    return (0 if d[1] <= d[0] else 1, -h.get("pos3", 0.0), strength(o))


def classify(o: Orf, code: int) -> None:
    """The class of an ORF from its hits (see the module notes)."""
    std = [h for h in o.hits if h["src"] == "family" and h["name"].lower() in STANDARD_FAMILIES]
    fam = [h for h in o.hits if (h["src"] == "family" and h["name"].lower() not in STANDARD_FAMILIES)
           or h["src"] == "pfam"]
    orf = [h for h in o.hits if h["src"] not in ("family", "pfam")]
    f = stop_share(o.piece_nt, code)
    orf = [h for h in orf if "p_pos3" in h]
    for h in orf:
        h["p_open"] = (1 - f) ** h["subs"]
        h["p_joint"] = fisher(h["p_open"], h["p_pos3"])
        h["selected"] = h["p_open"] <= P_EACH and h["p_pos3"] <= P_EACH
    if std:
        o.best = min(std, key=lambda h: h["evalue"])
        o.cls, o.family = "standard_fragment", o.best["name"]
    elif fam:
        o.best = min(fam, key=lambda h: (h["src"] != "pfam", h["evalue"]))
        o.cls, o.family = "known_family", o.best["name"]
    elif any(h["selected"] for h in orf):
        o.best = min((h for h in orf if h["selected"]), key=lambda h: h["p_joint"])
        o.cls, o.p_open = "conserved", o.best["p_open"]
    elif o.e_len <= E_LEN:
        o.cls = "long_unconfirmed"
        if orf:
            o.best = min(orf, key=lambda h: h["p_open"] * h["p_pos3"])
            o.p_open = o.best["p_open"]


def evidence(o: Orf) -> str:
    text = _evidence(o)
    if o.runs_into and o.cls in ("conserved", "known_family", "long_unconfirmed"):
        text += (f"; not written: no stop codon before the frame runs on into {o.runs_into}"
                 + (" - perhaps a part of that gene (start further upstream?)" if "(in frame)" in o.runs_into else
                    "") + " - check")
    return text


def _evidence(o: Orf) -> str:
    h = o.best
    if o.cls == "conserved":
        return (f"homolog {h['name']} ({h['src']}), blastp E={h['evalue']:.0e}, {h['identity']:.0f}% identity; "
                f">= {h['subs']} substitutions among the copies: all open by chance p={o.p_open:.0e}, "
                f"{100 * h['pos3']:.0f}% of the {h['nt_diffs']} differences of a nucleotide alignment at codon position 3 "
                f"(p={h['p_pos3']:.0e})")
    if o.cls == "known_family" and h.get("src") == "pfam":
        return f"Pfam {o.family} domain (hmmsearch, above the Pfam gathering threshold, E={h['evalue']:.0e})"
    if o.cls == "known_family":
        return f"similar to {o.family} (MFannot protein family), blastp E={h['evalue']:.0e}"
    if o.cls == "standard_fragment":
        return (f"similar to {o.family} (E={h['evalue']:.0e}): a second copy, a pseudogene or a piece of {o.family} "
                "the annotation missed - check")
    if o.cls == "long_unconfirmed":
        tail = (f"; homolog {h['name']} shows no selection ({h['subs']} differences, p_open={h['p_open']:.2f}, "
                f"position 3 p={h['p_pos3']:.2f})" if h else "; no homolog")
        return f"longer than chance ORFs (E={o.e_len:.0e}){tail} - not confirmed, check"
    return f"homolog {h['name']} too similar to show selection (p={o.p_open:.2f}); length within chance"


# ── writing ──────────────────────────────────────────────────────────────────────────────────────────────────────
def location(o: Orf, n: int):
    """The ORF's location; one that ends at a tRNA without a stop has an open 3' end (> / <)."""
    from Bio.SeqFeature import AfterPosition, BeforePosition, ExactPosition
    lo, hi = o.start, o.end
    lo_pos = BeforePosition(lo) if (o.trna_end and o.strand == -1) else ExactPosition(lo)
    if hi <= n:
        hi_pos = AfterPosition(hi) if (o.trna_end and o.strand == 1) else ExactPosition(hi)
        return FeatureLocation(lo_pos, hi_pos, strand=o.strand)
    hi_pos = AfterPosition(hi - n) if (o.trna_end and o.strand == 1) else ExactPosition(hi - n)
    parts = [FeatureLocation(lo_pos, n, strand=o.strand), FeatureLocation(0, hi_pos, strand=o.strand)]
    return CompoundLocation(parts if o.strand == 1 else parts[::-1])


def write_orf(rec: SeqRecord, o: Orf, code: int) -> None:
    """A gene and CDS feature for the ORF: /gene="orf<codons>" /product="hypothetical protein"; a known family with a
    mitochondrial product name (atp9 of sponges and cnidarians) gets that gene name and product."""
    loc = location(o, len(rec.seq))
    name, product = o.name, "hypothetical protein"
    what = f"ORF of unknown function ({o.cls})"
    if o.cls == "known_family" and o.family in PRODUCTS:
        name, product, what = o.family, PRODUCTS[o.family], f"{o.family} by similarity (known protein family)"
    elif o.cls == "known_family" and o.family in PFAM_PRODUCT:
        product, what = PFAM_PRODUCT[o.family], f"Pfam {o.family} (gathering threshold)"
    note = f"{NOTE}: {what}; {evidence(o)}"
    if o.trna_end:
        note += (f"; no stop codon: the reading frame runs on into {o.trna_end} (mitochondrial transcripts are cut at "
                 f"tRNAs); 3' end at the last whole codon before it")
    quals = {"gene": [name], "codon_start": ["1"], "transl_table": [str(code)], "product": [product],
             "note": [note]}
    for k, v in stop_qualifiers(rec, loc, code).items():          # incomplete stop: /transl_except, polyA note
        quals[k] = quals.get(k, []) + v
    quals["translation"] = [o.prot]
    rec.features.append(SeqFeature(loc, type="gene", qualifiers={"gene": [name], "note": [note]}))
    rec.features.append(SeqFeature(loc, type="CDS", qualifiers=quals))
    o.written = True


def reference_orfs(refs: Iterable[SeqRecord], default_code: int, exclude: set) -> Dict[str, Tuple[str, str]]:
    """Library entries from GenBank records of relatives (--orf_ref): their ORFs in unannotated sequence and their
    non-standard CDS, as scripts/build_orf_library.py makes them."""
    out: Dict[str, Tuple[str, str]] = {}
    for rec in refs:
        if rec.id in exclude or rec.id.split(".")[0] in exclude:
            continue
        code = code_of(rec, default_code)
        for f in rec.features:
            if f.type == "CDS" and not gene_key(f) and not is_orf(f):
                nt = str(f.extract(rec.seq)).upper()
                if len(nt) >= 3 * MIN_AA:
                    out[f"{rec.id}|{int(f.location.start) + 1}|{int(f.location.end)}|cds|ref"] = (
                        str(Seq(nt[: len(nt) // 3 * 3]).translate(table=code)).replace("*", "X"), nt)
        for o in candidates(rec, code):
            out[f"{rec.id}|{o.start + 1}|{o.end}|ncr|ref"] = (o.piece, o.piece_nt)
    return out


def annotate(recs: List[SeqRecord], default_code: int, mode: str = "yes", library: Optional[Path] = LIBRARY,
             threads: int = 1, bin_dir: Optional[str] = None, exclude: Iterable[str] = (),
             refs: Iterable[SeqRecord] = ()) -> List[Orf]:
    """Find, test and (mode 'yes' / 'all'; 'report': none) write the ORFs of all records together (each record's ORFs
    are homologs for the others). Returns every reported ORF."""
    excl = {x for x in exclude if x} | {r.id for r in recs} | {r.id.split(".")[0] for r in recs}
    batch: List[Orf] = []
    codes = {}
    for rec in recs:
        rec.features = [f for f in rec.features if not is_orf(f)]
        codes[rec.id] = code_of(rec, default_code)
        batch += candidates(rec, codes[rec.id])
    if not batch:
        return []
    search(batch, library, excl, threads, bin_dir, reference_orfs(refs, default_code, excl))
    pfam_search(batch, threads, bin_dir)
    for o in batch:
        classify(o, codes[o.record])
    rank = {"standard_fragment": 0, "known_family": 1, "conserved": 2, "long_unconfirmed": 3}
    reported = sorted((o for o in batch if o.cls),
                      key=lambda o: (rank[o.cls], *(frame_rank(o) if o.cls == "conserved" else (0, 0, strength(o)))))
    by_id = {r.id: r for r in recs}
    kept: List[Orf] = []
    for o in reported:                       # strongest first; an ORF overlapping a stronger one is its shadow
        top = next((x for x in kept if x.record == o.record and overlap(x, o, len(by_id[o.record].seq)) >
                    0.3 * min(x.end - x.start, o.end - o.start)), None)
        if top is not None:
            o.shadow_of = f"{top.name} {top.start + 1}..{top.end}"
            continue
        kept.append(o)
        if mode != "report" and o.writable and o.shadow_of == "" and (o.cls in ("conserved", "known_family") or
                                                 (mode == "all" and o.cls == "long_unconfirmed")):
            write_orf(by_id[o.record], o, codes[o.record])
    for rec in recs:                          # the new features in record order (a record without one: untouched)
        if any(o.written and o.record == rec.id for o in kept):
            rec.features.sort(key=lambda f: (f.type != "source", int(f.location.parts[0].start) if f.location else 0,
                                              f.type != "gene"))
    return sorted(kept, key=lambda o: (o.record, o.start))


def strength(o: Orf) -> float:
    """Smaller is stronger: the chance explanation of the ORF's evidence."""
    if o.cls == "conserved":
        return o.best["p_joint"]
    if o.best is not None and o.cls in ("standard_fragment", "known_family"):
        return o.best["evalue"]
    return o.e_len


def overlap(x: Orf, y: Orf, n: int) -> int:
    """Shared bases of two ORFs of one record (either strand)."""
    a = set(np.arange(x.start, x.end) % n)
    return len(a.intersection(np.arange(y.start, y.end) % n))


def report_rows(orfs: Sequence[Orf]) -> List[dict]:
    rows = []
    for o in orfs:
        h = o.best or {}
        rows.append({"record": o.record, "orf": o.name, "class": o.cls, "written": "yes" if o.written else "no",
                     "start": o.start + 1, "end": o.end, "strand": "+" if o.strand == 1 else "-",
                     "codons": o.codons, "e_len": f"{o.e_len:.1e}", "homolog": h.get("name", ""),
                     "homolog_evalue": f"{h['evalue']:.0e}" if h else "", "homolog_identity": h.get("identity", ""),
                     "substitutions": h.get("subs", ""), "p_open": f"{o.p_open:.1e}" if "p_pos3" in h else "",
                     "pos3_share": f"{h['pos3']:.2f}" if "p_pos3" in h else "",
                     "p_pos3": f"{h['p_pos3']:.1e}" if "p_pos3" in h else "",
                     "family": o.family, "evidence": evidence(o)})
    return rows


def warnings(orfs: Sequence[Orf]) -> List[str]:
    """One line per ORF the user should look at: long ORFs not written, fragments of the 13 genes."""
    return [f"{o.record} {o.start + 1}..{o.end} ({'+' if o.strand == 1 else '-'}, {o.codons} codons): "
            f"{o.cls}: {evidence(o)}" for o in orfs if not o.written]


def ncr_notes(orfs: Sequence[Orf]) -> Dict[str, List[Tuple[int, int, str]]]:
    """record -> (start, end, text) of the ORFs not written that an NCR /note should name."""
    out: Dict[str, List[Tuple[int, int, str]]] = {}
    for o in orfs:
        if not o.written:
            out.setdefault(o.record, []).append(
                (o.start, o.end, f"contains ORF {o.start + 1}..{o.end} ({'+' if o.strand == 1 else '-'}, "
                                 f"{o.codons} codons): {o.cls}, check"))
    return out


def notes_from_report(rows: Iterable[dict]) -> Dict[str, List[Tuple[int, int, str]]]:
    """ncr_notes from the rows of orf_report.tsv (the pipeline runs this module as a program)."""
    out: Dict[str, List[Tuple[int, int, str]]] = {}
    for r in rows:
        if r.get("written") == "no":
            s, e = int(r["start"]) - 1, int(r["end"])
            out.setdefault(r["record"], []).append(
                (s, e, f"contains ORF {s + 1}..{e} ({r['strand']}, {r['codons']} codons): {r['class']}, check"))
    return out


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="GenBank file or folder")
    p.add_argument("-o", "--out", required=True, help="output folder (GenBank files, orf_report.tsv)")
    p.add_argument("--genetic_code", type=int, default=5, help="for records without /transl_table (default 5)")
    p.add_argument("--orf", choices=["yes", "all", "report"], default="yes",
                   help="yes (default): write the ORFs a homolog shows to be under selection and known families; "
                        "all: also ORFs longer than chance without such a homolog; report: write none")
    p.add_argument("--library", default=str(LIBRARY), help=f"ORF library (default {LIBRARY.name}; '' = none)")
    p.add_argument("--ref", action="append", default=[],
                   help="GenBank file or folder of relatives (repeat): their ORFs are homologs too (a genome whose "
                        "relatives are not in RefSeq)")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--bin", default=None, help="folder with blastp, makeblastdb and hmmsearch")
    return p


def main(argv: Sequence[str] = None) -> int:
    a = build_parser().parse_args(argv)
    src, out = require_path(a.input), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    files = [src] if src.is_file() else iter_gb_files(src)
    recs_by_file = {f: list(SeqIO.parse(str(f), "genbank")) for f in files}
    recs = [r for rs in recs_by_file.values() for r in rs]
    refs = [r for x in a.ref for f in ([Path(x)] if Path(x).is_file() else iter_gb_files(Path(x)))
            for r in SeqIO.parse(str(f), "genbank")]
    orfs = annotate(recs, a.genetic_code, a.orf, Path(a.library) if a.library else None, a.threads, a.bin,
                    os.environ.get("MITOM_EXCLUDE_REF", "").split(","), refs)
    for f, rs in recs_by_file.items():
        write_gb(rs, out / (f.name if src.is_file() else f.relative_to(src)))
    write_tsv(out / "orf_report.tsv", COLUMNS, report_rows(orfs))
    for w in warnings(orfs):
        print(f"[mito_orf] WARNING {w}", file=sys.stderr)
    n_w = sum(o.written for o in orfs)
    print(f"[mito_orf] {len(orfs)} ORFs reported, {n_w} written, {len(warnings(orfs))} to check, "
          f"in {len(recs)} records -> {out}/orf_report.tsv", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(cli_guard(main))
