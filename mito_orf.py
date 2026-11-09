#!/usr/bin/env python3
"""mito_orf.py - open reading frames of unknown function in the sequence no gene covers.

Many animal mitogenomes hold a protein gene that is not one of the 13 (e.g. a conserved ~355-codon ORF between trnQ
and trnA in the spionid polychaetes Polydora and Boccardiella, which NCBI leaves unannotated). An ORF is a gene when
selection has kept it open; length alone cannot tell that: AT-rich sequence has few stop codons, so 80-110-codon ORFs
are common by chance, and a short real protein is no longer than a chance ORF. Two pieces of evidence, both from
first principles, decide it:

1. Length beyond chance (E_len): under the base composition of the record's unannotated sequence, the expected number
   of ORFs at least this long that start at a start codon of the genetic code,
       E_len = 2 * F * p_start * (1 - p_stop) ** L      (F nt unannotated, both strands; L codons),
   tells whether the ORF is longer than chance ORFs are.
2. Kept open (p_open): the ORF's protein has a homolog (blastp E <= 1e-5, low-complexity masked) that is itself an
   ORF: in the other records of the run, or in a library of ORFs and non-standard CDS of RefSeq animal mitogenomes
   (models/orf_library.fna.gz, scripts/build_orf_library.py). The nucleotide differences between the aligned codons
   of the two copies, d, are a lower bound of the substitutions since they split; in non-coding sequence a
   substitution makes a stop codon with probability f (the share of single-nucleotide changes of the ORF's own
   codons that give one), so both copies stayed open by chance with probability
       p_open = (1 - f) ** d.
   A close relative (few differences) proves nothing; p_open <= 0.01 shows selection. A homolog in a known protein
   family (models/mfannot_proteins, mitochondrial and plastid genes of all eukaryotes) is evidence by itself.

Classes (orf_report.tsv: every ORF of >= 40 codons longer than chance or with a homolog that shows selection):
    conserved          p_open <= 0.01                                  -> CDS orf<codons>, "hypothetical protein"
    known_family       homolog in a known family (not the 13)           -> CDS orf<codons>, the /note names the family
    long_unconfirmed   E_len <= 0.01, no homolog showing selection      -> not written (--orf all: written);
                                                                           a WARNING, and the NCR /note names it
    standard_fragment  homolog of one of the 13 protein genes           -> not written; a WARNING (second copy,
                                                                           pseudogene, or a piece the annotation missed)
ORFs within chance length without such a homolog are not reported (chance ORFs: in all frames of any AT-rich stretch).
A written ORF is a CDS from the first start codon of the unannotated stretch to its stop codon, or - when the frame
runs on into the next gene - to the base before it with an incomplete stop T / TA (/transl_except, as mitochondrial
genes end before a tRNA); /gene="orf<codons>", /product="hypothetical protein", /note with the evidence. Records of
an earlier run: their mito_orf genes are replaced.

Usage: python3 mito_orf.py -i <GenBank file or folder> -o <out folder> [--genetic_code 5] [--orf yes|all|report]
       (in `mitom run` / `annotate` / `import`: the step after curation, before the non-coding regions; --orf no: off)
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
MIN_AA = 40                 # shortest ORF considered (codons): below this a homology search has no power
E_LEN = 0.01                # length beyond chance: expected chance ORFs this long per record
P_OPEN = 0.01               # kept open by selection: chance that both copies stayed open
TESTED = 5                  # homologs (best blastp E) given the nucleotide selection tests
E_BLAST = 1e-5              # homolog: blastp E-value
REPEAT_ID = 0.85            # tandem repeat: copies >= 85 % identical ...
REPEAT_COPIES = 2           # ... over >= 2 periods (period 3-150 nt)
FREE = 0.9                  # share of a written ORF's codons in unannotated sequence (overlaps with a tRNA end)
NOTE = "found by MitoM mito_orf"      # start of the /note of every gene this module writes
COVERING = {"CDS", "rRNA", "tRNA", "gene", "mRNA", "intron", "exon"}
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
    """True where a gene, CDS, rRNA, tRNA or intron lies (a joined gene covers its introns: mito_ncr.extent);
    this module's own ORFs and non-coding / control-region features do not cover."""
    n = len(rec.seq)
    circular = rec.annotations.get("topology") == "circular"
    cov = np.zeros(n, dtype=bool)
    for f in rec.features:
        if f.type not in COVERING or f.location is None or is_orf(f):
            continue
        s, length = extent(f.location, n, circular)
        cov[np.arange(s, s + length) % n] = True
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
        if best[1] == k:                                    # runs to its own stop codon (which may touch a gene)
            b, writable = 3 * k + 3, True
        else:                                               # runs into a gene: end before it, incomplete stop
            b = a + int(np.argmin(free[a:])) if not free[a:].all() else 3 * k
            tail = nt_all[a:b][(b - a) // 3 * 3:]
            writable = tail in ("T", "TA")
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
    """One-sided binomial p that >= d[2] of the sum(d) differences fall on codon position 3 if each position were
    equally likely (1/3, as in non-coding sequence read in any frame). A protein-coding frame changes mostly at its
    synonymous third positions; the frames beside it (and the other strand) do not."""
    from math import comb
    n, k = sum(d), d[2]
    return float(sum(comb(n, i) * (1 / 3) ** i * (2 / 3) ** (n - i) for i in range(k, n + 1))) if n else 1.0


# ── homology ─────────────────────────────────────────────────────────────────────────────────────────────────────
def find_tool(name: str, bin_dir: Optional[str] = None) -> Optional[str]:
    if bin_dir and (Path(bin_dir) / name).exists():
        return str(Path(bin_dir) / name)
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
            h["s_nt"] = s_nt
            h["subs"] = codon_differences(qa, sa, o.piece_nt, int(qs) - 1, s_nt, int(ss) - 1)
        o.hits.append(h)
    for o in batch:                                          # the selection tests on the best few homologs
        for h in sorted((h for h in o.hits if "s_nt" in h), key=lambda h: h["evalue"])[:TESTED]:
            d = nt_differences(o, h.pop("s_nt"))
            h.update(nt_diffs=sum(d), pos3=d[2] / sum(d) if sum(d) else 0.0, p_pos3=third_position_p(d))
        for h in o.hits:
            h.pop("s_nt", None)


STANDARD_FAMILIES = {g.lower() for g in STANDARD}


def classify(o: Orf, code: int) -> None:
    """The class of an ORF from its hits (see the module notes)."""
    std = [h for h in o.hits if h["src"] == "family" and h["name"].lower() in STANDARD_FAMILIES]
    fam = [h for h in o.hits if h["src"] == "family" and h["name"].lower() not in STANDARD_FAMILIES]
    orf = [h for h in o.hits if h["src"] != "family"]
    f = stop_share(o.piece_nt, code)
    orf = [h for h in orf if "p_pos3" in h]
    for h in orf:
        h["p_open"] = (1 - f) ** h["subs"]
        h["selected"] = h["p_open"] <= P_OPEN and h["p_pos3"] <= P_OPEN
    if std:
        o.best = min(std, key=lambda h: h["evalue"])
        o.cls, o.family = "standard_fragment", o.best["name"]
    elif fam:
        o.best = min(fam, key=lambda h: h["evalue"])
        o.cls, o.family = "known_family", o.best["name"]
    elif any(h["selected"] for h in orf):
        o.best = min((h for h in orf if h["selected"]), key=lambda h: h["p_open"] * h["p_pos3"])
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
        return (f"homolog {h['name']} ({h['src']}), blastp E={h['evalue']:.0e}, {h['identity']:.0f}% identity, "
                f"{h['subs']} nucleotide differences: both open by chance p={o.p_open:.0e}, "
                f"{100 * h['pos3']:.0f}% of the {h['nt_diffs']} differences of a nucleotide alignment at codon position 3 "
                f"(p={h['p_pos3']:.0e})")
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
    if o.end <= n:
        return FeatureLocation(o.start, o.end, strand=o.strand)
    parts = [FeatureLocation(o.start, n, strand=o.strand), FeatureLocation(0, o.end - n, strand=o.strand)]
    return CompoundLocation(parts if o.strand == 1 else parts[::-1])


def write_orf(rec: SeqRecord, o: Orf, code: int) -> None:
    """A gene and CDS feature for the ORF: /gene="orf<codons>" /product="hypothetical protein"; a known family with a
    mitochondrial product name (atp9 of sponges and cnidarians) gets that gene name and product."""
    loc = location(o, len(rec.seq))
    name, product = o.name, "hypothetical protein"
    what = f"ORF of unknown function ({o.cls})"
    if o.cls == "known_family" and o.family in PRODUCTS:
        name, product, what = o.family, PRODUCTS[o.family], f"{o.family} by similarity (known protein family)"
    note = f"{NOTE}: {what}; {evidence(o)}"
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
    for o in batch:
        classify(o, codes[o.record])
    rank = {"standard_fragment": 0, "known_family": 1, "conserved": 2, "long_unconfirmed": 3}
    reported = sorted((o for o in batch if o.cls), key=lambda o: (rank[o.cls], strength(o)))
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
        return o.best["p_open"] * o.best["p_pos3"]
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


def main(argv: Sequence[str] = None) -> int:
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
    p.add_argument("--bin", default=None, help="folder with blastp and makeblastdb")
    a = p.parse_args(argv)
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
