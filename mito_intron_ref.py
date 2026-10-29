#!/usr/bin/env python3
"""mito_intron_ref.py - mitochondrial introns at known sites, transferred from annotated relatives.

Animal mitochondrial introns are homing introns: group I / II introns at conserved insertion sites of a few genes
(cnidarian cox1 / nad5, annelid cox1 / nad1 / nad4, placozoan cox1 / nad5 / rrnL, sponge cox1). This module maps every
intron site of annotated relatives onto a query genome and transfers the annotated boundaries; MFannot's method
(mito_intron.py) predicts introns de novo from protein alignments and splice rules. The pipeline runs the de novo step
first and this one on its output: sites of related records re-place de novo boundaries, add introns de novo misses
(mini-exons, trans-spliced pieces, rRNA introns), and de novo introns at sites no relative has stay.

Library (build, models/intron_sites.tsv.gz + .fa.gz, 76 KB): every joined CDS / rRNA of annotated GenBank records ->
one row per exon-exon junction: host gene, insertion site in the spliced reference (offset: residue and codon phase),
cis / trans splicing, intron length and group (annotated, else Infernal cmsearch with the Rfam group I / II models),
30 nt of exon and intron on each side of both intron ends (junction probes), offsets of an intron feature that differs
from the CDS gap, half-introns of trans-spliced junctions, the spliced protein / rRNA and the lineage. A CDS gap is
an intron when >= 142 nt (eval_introns / MFannot) or covered by an intron feature; overlapping parts (frameshift
edits), misordered origin joins and split rRNAs without intron feature are not sites. The shipped library is built
from the gold standard of the 193 local intron records (intron_gold_2026-10-29: NCBI placements verified, slid to
the canonical placement or corrected; erroneous ones removed).

Transfer (annotate_record), per gene of the library:
  1. Anchors: tblastn of the reference proteins (blastn of the rRNAs) on the doubled genome (circular records),
     HSPs split at query insertions >= 30 nt and scored per column (BLOSUM62 of the aligned codon); blastn-short of
     the 60-nt junction probes. References ranked by shared lineage, then alignment score: the closest 8, plus up to
     24 that bring a site not seen yet.
  2. Each reference junction: the upstream side (alignment of the upstream reference exon reaching within 8 units of
     the site, or the donor probe) and the downstream side; exon ends extrapolated from anchors <= 25 residues away.
     Present: both sides on one strand >= 100 nt apart (two strands only where the reference junction is
     trans-spliced); absent: an alignment running on across the site (at least half as strong as the pair) or two
     alignment pieces next to each other (at least as strong). A side from an alignment that runs on across the site
     must beat that alignment's own continuation (nad4 paralogue joined to nad5; an intron start copying exon 2).
  3. Exon ends: the three strongest pairs are refined - every (donor, acceptor) near them that keeps the reading frame
     between the anchors and adds no stop codon, scored by (1) the spliced query across the junction aligned to the
     reference protein, (2) position-weighted identity of the four junction arms to the reference probes (the
     reference's annotated boundary, transferred), (3) group I / II splice rules - and the best score wins. A call
     needs junction identity >= 0.5 (>= 0.8 from another phylum; the query's phylum from --lineage, else from its
     closest library relative).
  4. One call per query site: weighted vote over the references (identity x 2^(lineage depth - deepest depth));
     sites chained in transcript order into gene models (separately annotated parts such as Trichoplax COX1a / COX1b
     stay separate; overlapping chains of one gene: the better supported one); gene start / stop from an existing
     in-frame CDS end, else the nearest start / first stop where the reference protein ends (partial when absent).
Write-back (apply_models): an existing (de novo) model of the gene is the base - at a de novo intron that is a library
site the transferred ends replace the de novo ones when the site comes from a near-identical relative (junction
identity >= 0.9 or the record's own annotation) or the de novo placement breaks the splice rules; a transferred site
inside a de novo exon splits it; de novo introns at no library site stay. A model reaching outside the de novo model
(trans-spliced pieces) is the base instead, unless from distant relatives and dropping at least as many de novo
introns as it adds. Then every intron not from the record's own annotation is placed by the gold-standard rules: as
it is when it meets its group's rule (I: exon ends T, intron ends G; II: intron GTGYG...AY); else the nearest
canonical placement in the sliding window (identical spliced sequence); else the nearest shift <= 60 nt that keeps
the protein (synonymous) and is canonical; then the group II 5' motif alone; else as it is. The group: the intron's
covariance-model hit (de novo note) when E <= 1e-5, else the reference's group, else either. A site whose library
members mostly sit at a placement meeting neither group rule (the hexacoral nad5 intron, exon ...T | intron ...A) is
placed by that clade consensus: the transferred placement stays (gold-standard family rule).
Own annotation: sites of the query's own annotated record (an NCBI GenBank given to `mitom run`; own=, CLI --ref_gb)
outrank everything and keep their quirks (intron feature off its CDS gap, half-introns); other records' quirks, even
of the same species, are not transferred.

Deviations and limits: a site no library record has comes only from the de novo step (boundaries by its rules);
placements the gold standard corrected for homology (class D-corrected) are reproduced only through relatives that
carry the corrected site; an intron across the origin of a circular record is given by the CDS join (GenBank allows
no multi-interval intron, table2asn MultiIntervalIntron); a start / stop not found gives a 5' / 3'-partial CDS.
Benchmark (benchmark/bench_intron_ref.py; 193 intron records + 77 intron-free, pipeline input, de novo first, final
gold standard 2026-10-29): whole library (pipeline default) and own annotation 309 / 309 exact; leave-one-out 298
exact, 4 within 6 nt, 1 more, 6 missed (de novo alone 295 / 4 / 4 / 6); 0 false introns in the intron-free records.

API: build_library(paths) / save_library / load_library; ref_genes(record) (sites of one annotated record);
annotate_record(record, library, genetic_code, lineage=, exclude=, own=) -> (gene models, report rows);
apply_models(record, models, genetic_code) writes gene / CDS (or rRNA) / intron features, merged with a de novo model.
Usage:
  python3 mito_intron_ref.py build -i gold_gb_dir ... -o models/intron_sites.tsv.gz --jobs 8 --bin env/bin
  python3 mito_intron_ref.py annotate -i final_gb_dir -o out_dir --genetic_code 4 [--lineage "Eukaryota;Metazoa;..."]
      [--ref_gb annotated_input.gb [--own_only]] [--bin env/bin] [--jobs N]
Outputs: out_dir/<same relative paths> (changed records rewritten, others copied), out_dir/intron_ref_report.tsv.
Leave-out runs (mitom_lo): MITOM_EXCLUDE_REF=<acc> leaves that record (and its RefSeq / INSDC copy) out of the library,
MITOM_EXCLUDE_LEVEL=species|genus its species / genus too.
"""
from __future__ import annotations

import argparse
import dataclasses
import gzip
import io
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mito_names  # noqa: E402
from mito_common import codon_translate, iter_gb_files, run_tool, stop_qualifiers, write_gb  # noqa: E402

DEFAULT_LIB = HERE / "models" / "intron_sites.tsv.gz"
MIN_INTRON = 142        # a CDS gap at least this long is an intron without an intron feature (eval_introns, MFannot)
PROBE = 30              # nt of exon / intron sequence kept on each side of an intron end

# ── gene names: mito_names.py (atp9 included: sponge / cnidarian introns) ─────────────────────────────────────────────
def gene_name(q: dict, ftype: str) -> str:
    """Normalised PCG / rRNA name (cox1, nad4l, rrnL ...); '' when not a standard one. COX1a / COX1b -> cox1."""
    return mito_names.gene_name(q, ftype, atp9=True)


# ── library ──────────────────────────────────────────────────────────────────
@dataclass
class Junction:
    """One exon-exon junction of a reference gene, i.e. one intron site."""
    k: int                  # 1-based junction number in the gene
    offset: int             # spliced nt of the gene upstream of the junction (from the first coding nt)
    splice: str             # cis / trans
    length: int             # cis: intron length (exon gap); trans: -1
    itype: str              # I / II / '' (annotated type)
    note: str               # annotated intron note
    e5: str = ""            # transcript-strand probes around the exon gap: exon end | intron start ... intron end | exon
    i5: str = ""
    i3: str = ""
    e3: str = ""
    d5: int = 0             # annotated intron feature start / end minus the exon gap's, transcript orientation
    d3: int = 0
    half5: int = 0          # trans: length of an annotated half-intron after the donor exon / before the acceptor exon
    half3: int = 0
    gap5: int = 0           # trans: nt between the donor exon and its half-intron (0 = adjacent) / acceptor likewise
    gap3: int = 0
    iseq: str = field(default="", repr=False)    # cis intron sequence while the library is built (typing; not saved)
    family: bool = False    # site whose annotated placement is the clade's rule, not the group I / II one (mark_family_sites)


@dataclass
class RefGene:
    gid: str
    acc: str
    same_as: List[str]      # other accessions of the same sequence (RefSeq copy / source record)
    organism: str
    lineage: List[str]
    code: int
    gene: str
    kind: str               # CDS / rRNA
    exon_lens: List[int]
    spliced: str            # spliced nt from the first coding nt (CDS) / the whole rRNA
    protein: str            # CDS translation ('' for rRNA)
    junctions: List[Junction] = field(default_factory=list)
    own: bool = False       # from the query record's own annotation (always the closest reference; not saved)


def iter_records(paths: Iterable[str]) -> Iterator:
    """GenBank records from files / folders (multi-record files allowed); text pre-filter on joins / introns."""
    for p in map(Path, paths):
        files = iter_gb_files(p) if p.is_dir() else [p]
        for f in files:
            buf: List[str] = []
            with open(f, errors="ignore") as fh:
                for line in fh:
                    buf.append(line)
                    if line.startswith("//"):
                        txt = "".join(buf)
                        buf = []
                        if "join(" in txt:
                            yield txt


def _circ(seq: str, start: int, end: int) -> str:
    """seq[start:end] on a circle (start may be negative, end may exceed len)."""
    L = len(seq)
    if 0 <= start and end <= L:
        return seq[start:end]
    return "".join(seq[i % L] for i in range(start, end))


def _rc(s: str) -> str:
    return str(Seq(s).reverse_complement())


def tx_probe(seq: str, pos: int, strand: int, before: int, after: int) -> Tuple[str, str]:
    """(before nt ending at transcript-position pos, after nt starting at pos+1) on strand; pos is the genomic
    0-based index of the last nt of the 'before' block in transcript order."""
    if strand == 1:
        return _circ(seq, pos - before + 1, pos + 1), _circ(seq, pos + 1, pos + 1 + after)
    a = _rc(_circ(seq, pos, pos + before))            # transcript goes to lower coordinates
    b = _rc(_circ(seq, pos - after, pos))
    return a, b


def intron_kind(text: str) -> str:
    """Intron group from an annotation note: 'I', 'II' or ''."""
    t = text.lower()
    if re.search(r"group[ _-]?(ii|2)|gpii|domain v", t):
        return "II"
    if re.search(r"group[ _-]?(i|1)(?![iv0-9])|gpi\b", t):
        return "I"
    return ""


def ref_genes(rec) -> List[RefGene]:
    """Joined CDS / rRNA features of one annotated record -> RefGene with its junctions (intron sites)."""
    if not rec.seq.defined:
        return []                                         # CONTIG-only record without sequence
    seq = str(rec.seq).upper()
    L = len(seq)
    circular = rec.annotations.get("topology", "circular") != "linear"
    same = sorted(set(re.findall(r"(?:identical to|derived from)\s+([A-Z]{1,2}_?\d+(?:\.\d+)?)",
                                 rec.annotations.get("comment", ""))))
    organism = rec.annotations.get("organism", "")
    lineage = list(rec.annotations.get("taxonomy", []))
    introns = [f for f in rec.features if f.type == "intron"]
    out: List[RefGene] = []
    seen: Dict[str, int] = {}
    for f in rec.features:
        if f.type not in ("CDS", "rRNA") or len(f.location.parts) < 2:
            continue
        g = gene_name(f.qualifiers, f.type)
        if not g:
            continue
        parts = [(int(p.start), int(p.end), p.strand or 1) for p in f.location.parts]
        try:
            nt = str(f.location.extract(rec.seq)).upper()
        except Exception:
            continue
        code = int(f.qualifiers.get("transl_table", ["1" if f.type == "CDS" else "0"])[0] or 0)
        cs = int(f.qualifiers.get("codon_start", ["1"])[0]) - 1 if f.type == "CDS" else 0
        prot = ""
        if f.type == "CDS":
            cod = nt[cs:]
            own = codon_translate(cod, code or 1).rstrip("*")
            prot = (f.qualifiers.get("translation", [""])[0] or "").upper() or own
            if len(own) < 0.9 * len(prot) or sum(a == b for a, b in zip(own[1:], prot[1:])) < 0.9 * (len(prot) - 1):
                continue                                  # join does not encode its /translation (misordered parts)
        seen[g] = seen.get(g, 0) + 1
        rg = RefGene(f"{rec.id}|{g}|{seen[g]}", rec.id, same, organism, lineage, code, g, f.type,
                     [e - s for s, e, _ in parts], nt[cs:], prot)
        off = -cs
        k = 0
        for (ps, pe, pst), (qs, qe, qst) in zip(parts, parts[1:]):
            off += pe - ps
            if min(pe - ps, qe - qs) < 4:
                continue                                      # 1-3 nt "exon": an annotation artefact
            if pst == qst:
                if pst == 1:
                    inside = ps <= qs < pe                    # next part starts inside this one: overlap
                    gap = qs - pe if qs >= pe else (L - pe + qs if circular else -1)
                    g0 = pe                                   # genomic start of the gap
                else:
                    inside = qs < pe <= qe or qs < ps < qe
                    gap = ps - qe if ps >= qe else (ps + L - qe if circular else -1)
                    g0 = qe
                if inside:
                    continue                                  # overlapping parts: frameshift / editing join
                if gap < 0:
                    splice, feat = "trans", None
                else:
                    span = (g0, g0 + gap)
                    feat = _overlapping_intron(introns, span, pst, L, g)
                    if feat is None and (gap < MIN_INTRON or f.type != "CDS"):
                        continue                              # frameshift join / split rRNA, not an intron
                    splice = "cis"
            else:
                splice, feat = "trans", None
            k += 1
            j = Junction(k, off, splice, gap if splice == "cis" else -1, "", "")
            donor = pe - 1 if pst == 1 else ps                # genomic index of the last nt of the upstream exon
            acc_ = qs - 1 if qst == 1 else qe                 # genomic index of the nt before the downstream exon
            j.e5, j.i5 = tx_probe(seq, donor, pst, PROBE, PROBE)
            j.i3, j.e3 = tx_probe(seq, acc_, qst, PROBE, PROBE)
            if splice == "cis":
                j.iseq = _circ(seq, g0, g0 + gap) if pst == 1 else _rc(_circ(seq, g0, g0 + gap))
            if splice == "cis" and feat is not None:
                fs, fe = _span(feat, L)
                if pst == 1:
                    j.d5, j.d3 = _wrap(fs - g0, L), _wrap(fe - (g0 + gap), L)
                else:
                    j.d5, j.d3 = _wrap((g0 + gap) - fe, L), _wrap(g0 - fs, L)
                note = " ".join(feat.qualifiers.get("note", []) + feat.qualifiers.get("product", []))
                j.note, j.itype = note[:120], intron_kind(note)
            elif splice == "trans":
                h5 = _adjacent_intron(introns, donor, pst, L, after=True)
                h3 = _adjacent_intron(introns, acc_, qst, L, after=False)
                notes = []
                if h5:
                    j.half5, j.gap5 = h5[0], h5[1]
                    notes.append(h5[2])
                if h3:
                    j.half3, j.gap3 = h3[0], h3[1]
                    notes.append(h3[2])
                j.note = " | ".join(notes)[:120]
                j.itype = intron_kind(j.note)
            rg.junctions.append(j)
        if rg.junctions:
            out.append(rg)
    return out


def _span(f, L: int) -> Tuple[int, int]:
    """Genomic span (start, end) of a feature, end > L when it is written across the origin."""
    parts = list(f.location.parts)
    s, e = int(f.location.start), int(f.location.end)
    strand = f.location.strand or 1
    if len(parts) > 1 and any((int(q.start) < int(p.start)) == (strand == 1) for p, q in zip(parts, parts[1:])):
        s, e = (int(parts[0].start), int(parts[-1].end) + L) if strand == 1 else (int(parts[-1].start), int(parts[0].end) + L)
    return s, e


def _wrap(d: int, L: int) -> int:
    """A coordinate difference on a circle of length L, in (-L/2, L/2]."""
    d %= L
    return d - L if d > L // 2 else d


def _ov(a: Tuple[int, int], b: Tuple[int, int], L: int) -> int:
    return max(max(0, min(a[1], b[1] + k) - max(a[0], b[0] + k)) for k in (-L, 0, L))


def _overlapping_intron(introns, span: Tuple[int, int], strand: int, L: int, gene: str):
    """The intron feature that is this exon gap: same strand, each covering >= 80 % of the other (an intron of
    another gene can lie inside the gap: cox1 intron in the corallimorpharian nad5 intron), not of another gene."""
    best, bo = None, 0
    for f in introns:
        if (f.location.strand or 1) != strand:
            continue
        g = gene_name({"gene": f.qualifiers.get("gene", [""])}, "CDS") or gene_name({"gene": f.qualifiers.get("gene", [""])}, "rRNA")
        if g and g != gene:
            continue
        fs = _span(f, L)
        o = _ov(span, fs, L)
        if o > bo and o >= 0.8 * max(span[1] - span[0], fs[1] - fs[0]):
            best, bo = f, o
    return best


def _adjacent_intron(introns, pos: int, strand: int, L: int, after: bool, tol: int = 3):
    """Intron feature starting just after transcript position pos (after=True) or ending just before it:
    (length, gap nt, note)."""
    for f in introns:
        if (f.location.strand or 1) != strand:
            continue
        s, e = _span(f, L)
        if strand == 1:
            gap = (s - (pos + 1)) % L if after else (pos + 1 - e) % L
        else:
            gap = (pos - e) % L if after else (s - pos) % L
        if gap <= tol:
            note = " ".join(f.qualifiers.get("note", []) + f.qualifiers.get("product", []))
            return e - s, gap, note
    return None


# ── library file ─────────────────────────────────────────────────────────────
SITE_COLS = ["gid", "acc", "same_as", "organism", "lineage", "code", "gene", "kind", "exon_lens", "k", "offset",
             "splice", "length", "itype", "note", "e5", "i5", "i3", "e3", "d5", "d3", "half5", "gap5", "half3", "gap3"]


def _open(path: Path, mode: str):
    return gzip.open(path, mode + "t", encoding="utf-8") if str(path).endswith(".gz") else open(path, mode, encoding="utf-8")


def fasta_path(tsv: Path) -> Path:
    name = tsv.name
    for ext in (".tsv.gz", ".tsv"):
        if name.endswith(ext):
            return tsv.with_name(name[: -len(ext)] + (".fa.gz" if ext == ".tsv.gz" else ".fa"))
    return tsv.with_suffix(".fa")


def save_library(genes: List[RefGene], tsv: Path) -> None:
    """Sites table (one row per junction) + FASTA of the spliced sequences (>gid protein, >gid|nt spliced nt)."""
    tsv.parent.mkdir(parents=True, exist_ok=True)
    with _open(tsv, "w") as fh:
        fh.write("\t".join(SITE_COLS) + "\n")
        for g in genes:
            for j in g.junctions:
                row = [g.gid, g.acc, ",".join(g.same_as), g.organism, ";".join(g.lineage), str(g.code), g.gene, g.kind,
                       ",".join(map(str, g.exon_lens)), str(j.k), str(j.offset), j.splice, str(j.length), j.itype,
                       j.note.replace("\t", " "), j.e5, j.i5, j.i3, j.e3, str(j.d5), str(j.d3), str(j.half5),
                       str(j.gap5), str(j.half3), str(j.gap3)]
                fh.write("\t".join(row) + "\n")
    with _open(fasta_path(tsv), "w") as fh:
        for g in genes:
            if g.protein:
                fh.write(f">{g.gid}\n{g.protein}\n")
            fh.write(f">{g.gid}|nt\n{g.spliced}\n")


def load_library(tsv: Path) -> List[RefGene]:
    seqs: Dict[str, str] = {}
    with _open(fasta_path(tsv), "r") as fh:
        for r in SeqIO.parse(fh, "fasta"):
            seqs[r.id] = str(r.seq)
    genes: Dict[str, RefGene] = {}
    with _open(tsv, "r") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            v = dict(zip(head, line.rstrip("\n").split("\t")))
            g = genes.get(v["gid"])
            if g is None:
                g = genes[v["gid"]] = RefGene(v["gid"], v["acc"], [x for x in v["same_as"].split(",") if x],
                                              v["organism"], [x for x in v["lineage"].split(";") if x], int(v["code"]),
                                              v["gene"], v["kind"], [int(x) for x in v["exon_lens"].split(",")],
                                              seqs.get(v["gid"] + "|nt", ""), seqs.get(v["gid"], ""))
            g.junctions.append(Junction(int(v["k"]), int(v["offset"]), v["splice"], int(v["length"]), v["itype"],
                                        v["note"], v["e5"], v["i5"], v["i3"], v["e3"], d5=int(v["d5"]), d3=int(v["d3"]),
                                        half5=int(v["half5"]), half3=int(v["half3"]), gap5=int(v["gap5"]),
                                        gap3=int(v["gap3"])))
    lib = list(genes.values())
    mark_family_sites(lib)
    return lib


FAMILY_MIN = 3          # members a site needs before its majority placement can be the clade's rule


def mark_family_sites(lib: List[RefGene]) -> None:
    """Flag the junctions of sites placed by a clade rule instead of the group I / II rules (gold standard, criterion
    4): per gene and phylum, junctions within 1 % of the protein length form a site; when most of its >= 3 members are
    annotated at a placement that meets neither group rule (the hexacoral nad5 intron: exon ...T | intron ...A), the
    annotated (consensus) placement is the rule and a transferred placement is not moved."""
    from collections import defaultdict
    by = defaultdict(list)
    for g in lib:
        n = max(1, len(g.protein) if g.kind == "CDS" else len(g.spliced))
        unit = 3 if g.kind == "CDS" else 1
        for j in g.junctions:
            if j.splice == "cis":
                by[(g.gene, phylum_of(g.lineage))].append((j.offset / unit / n, j))
    for items in by.values():
        items.sort(key=lambda x: x[0])
        clusters: List[list] = []
        for x in items:
            if clusters and x[0] - clusters[-1][-1][0] <= 0.01:
                clusters[-1].append(x)
            else:
                clusters.append([x])
        for cl in clusters:
            canon = sum(1 for _, j in cl if {"I", "II"} & set(canonical_kinds(j.e5[-1:], j.i5[:6], j.i3[-6:])))
            if len(cl) >= FAMILY_MIN and canon < len(cl) / 2:
                for _, j in cl:
                    j.family = True


def build_library(paths: List[str], jobs: int = 1) -> List[RefGene]:
    """RefGene for every joined CDS / rRNA with an intron in the GenBank sources; a record met twice is kept once."""
    texts = iter_records(paths)
    if jobs > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(jobs) as pool:
            found = [g for gs in pool.imap(_genes_of_text, texts, chunksize=20) for g in gs]
    else:
        found = [g for t in texts for g in _genes_of_text(t)]
    out: Dict[str, RefGene] = {}
    for g in found:
        out.setdefault(g.gid, g)
    return sorted(out.values(), key=lambda g: (g.gene, g.lineage, g.acc))


DEFAULT_CM = HERE / "models" / "rfam_group_I_II_introns.cm"


def type_introns(genes: List[RefGene], cmsearch: Optional[str], cm: Path = DEFAULT_CM, cpu: int = 4) -> int:
    """Group I / II of every library intron without an annotated type, by Infernal cmsearch with the Rfam models
    (RF00028 group I; RF00029 and the group II D1-D4 models), as mito_intron.py types introns; the type drives the
    splice rules when a site is transferred. Returns the number typed."""
    import tempfile
    todo = {f"i{n}": j for n, (g, j) in enumerate((g, j) for g in genes for j in g.junctions) if j.iseq and not j.itype}
    if not todo or not cmsearch or not cm.exists():
        return 0
    best: Dict[str, Tuple[float, str]] = {}
    with tempfile.TemporaryDirectory() as td:
        fa, tbl = Path(td) / "introns.fa", Path(td) / "hits.tbl"
        fa.write_text("".join(f">{k}\n{j.iseq}\n" for k, j in todo.items()))
        run_tool([cmsearch, "--noali", "--toponly", "--rfam", "--cpu", cpu, "--tblout", tbl, cm, fa])
        for line in tbl.read_text().splitlines() if tbl.exists() else []:
            if line.startswith("#"):
                continue
            c = line.split()
            name, model, e = c[0], c[2], float(c[15])
            if e <= 0.01 and (name not in best or e < best[name][0]):
                best[name] = (e, model)
    for k, (e, model) in best.items():
        j = todo[k]
        j.itype = "I" if model == "Intron_gpI" else "II"
        j.note = (j.note + " " if j.note else "") + f"[type by Infernal/Rfam {model}, E={e:.2g}]"
    return len(best)


def _genes_of_text(txt: str) -> List[RefGene]:
    try:
        return ref_genes(SeqIO.read(io.StringIO(txt), "genbank"))
    except Exception as exc:                              # one malformed record must not stop a library build
        print(f"[mito_intron_ref] record skipped: {exc!r}", file=sys.stderr)
        return []


# ── query genome ─────────────────────────────────────────────────────────────
class Genome:
    """Both strands of a record. Strand coordinates: on +1 the genomic 0-based index, on -1 the index into the
    reverse complement (genomic L-1-c); on a circular record they run past L (taken modulo L)."""

    def __init__(self, seq: str, circular: bool):
        self.seq, self.L, self.circular = seq.upper(), len(seq), circular
        self.rc = _rc(self.seq)

    def sub(self, strand: int, a: int, b: int) -> str:
        s = self.seq if strand == 1 else self.rc
        if self.circular:
            return _circ(s, a, b) if b > a else ""
        return s[max(a, 0):max(b, 0)]

    def norm(self, strand: int, a: int, b: int) -> Optional[Tuple[int, int, int]]:
        """(strand, a, b) with 0 <= a < L and a < b (b past L across the origin); None if empty / impossible."""
        n = (b - a) % self.L if self.circular else b - a
        if n <= 0 or (not self.circular and (a < 0 or b > self.L)):
            return None
        a0 = a % self.L if self.circular else a
        return strand, a0, a0 + n

    def genomic(self, strand: int, a: int, b: int) -> List[Tuple[int, int]]:
        """Strand interval [a, b) -> genomic 0-based half-open parts in transcript order (split at the origin)."""
        L = self.L
        if strand == -1:
            a, b = L - b, L - a
        a0 = a % L
        b0 = a0 + (b - a)
        parts = [(a0, b0)] if b0 <= L else [(a0, L), (0, b0 - L)]
        return parts if strand == 1 else [(x, y) for x, y in reversed(parts)]


@dataclass
class Hsp:
    gid: str
    strand: int
    bits: float
    cols: Dict[int, int]        # reference unit (residue / nt) -> strand coordinate of its first nt
    sc: Dict[int, float] = field(default_factory=dict)     # reference unit -> alignment score of that column

    def seg(self, lo: int, hi: int) -> float:
        """Alignment score of the HSP's columns for reference units lo .. hi-1."""
        return sum(v for u, v in self.sc.items() if lo <= u < hi)

    @property
    def q0(self) -> int:
        return min(self.cols)

    @property
    def q1(self) -> int:
        return max(self.cols) + 1


def _blast(exe: str, query: Path, subject: Path, extra: List[str]) -> List[List[str]]:
    import subprocess
    cmd = [exe, "-query", str(query), "-subject", str(subject), "-outfmt",
           "6 qseqid qstart qend sstart send sframe bitscore evalue qseq sseq", *extra]
    try:
        out = run_tool(cmd, timeout=TIMEOUT).stdout
    except (subprocess.TimeoutExpired, OSError) as exc:
        print(f"[mito_intron_ref] {Path(exe).name} failed: {exc}", file=sys.stderr)
        return []
    return [line.split("\t") for line in out.splitlines() if line.count("\t") == 9]


TIMEOUT = 900
SPLIT_GAP = 30                  # a subject insertion of >= 30 nt inside one HSP is an intron: split the HSP there


def parse_hsps(rows: List[List[str]], L: int, unit: int, ids: Dict[str, List[str]]) -> Dict[str, List[Hsp]]:
    """BLAST rows on the doubled genome -> HSPs per reference gene id, in strand coordinates; HSPs split at long
    subject insertions; copies of one hit in both halves of the doubled genome kept once."""
    out: Dict[str, List[Hsp]] = {}
    for qid, qs, qe, ss, se, frame, bits, ev, qseq, sseq in rows:
        strand = 1 if int(frame) > 0 or (unit == 1 and int(se) >= int(ss)) else -1
        t = int(ss) - 1 if strand == 1 else 2 * L - int(ss)
        q = int(qs) - 1
        step = 3 if unit == 3 else 1
        pieces: List[Dict[int, int]] = [{}]
        last_t = None
        for a, b in zip(qseq, sseq):
            if a != "-" and b != "-":
                if last_t is not None and t - last_t > SPLIT_GAP + step:
                    pieces.append({})
                pieces[-1][q] = t
                last_t = t
                q += 1
                t += step
            elif a == "-":
                t += step
            else:
                q += 1
        for cols in pieces:
            if len(cols) < 3:
                continue
            shift = (min(cols.values()) // L) * L
            cols = {k: v - shift for k, v in cols.items()}
            for gid in ids[qid]:
                out.setdefault(gid, []).append(Hsp(gid, strand, float(bits) * len(cols) / max(1, int(qe) - int(qs) + 1),
                                                   cols))
    # one hit is found in both halves of the doubled genome, possibly cut differently at its ends: same strand and
    # diagonal (modulo L) with overlapping reference ranges -> keep the longest
    for gid, hs in out.items():
        keep: List[Hsp] = []
        for h in sorted(hs, key=lambda h: -len(h.cols)):
            u = min(h.cols)
            diag = (h.cols[u] - (3 if unit == 3 else 1) * u) % L
            if not any(k.strand == h.strand and (k.cols[min(k.cols)] - (3 if unit == 3 else 1) * min(k.cols)) % L == diag
                       and min(k.cols) <= h.q1 and h.q0 <= max(k.cols) for k in keep):
                keep.append(h)
        out[gid] = keep
    return out


_B62: Dict[Tuple[str, str], float] = {}


def score_hsps(G: Genome, hsps: Dict[str, List[Hsp]], lib: Dict[str, RefGene], code: int) -> None:
    """Per-column scores of every HSP: BLOSUM62 of reference residue vs translated query codon (CDS), +2 / -3 per nt
    (rRNA). Sides of a junction are compared on these local scores, not on whole-HSP bit scores."""
    if not _B62:
        from Bio.Align import substitution_matrices
        m = substitution_matrices.load("BLOSUM62")
        for a in m.alphabet:
            for b in m.alphabet:
                _B62[(a, b)] = float(m[a][b])
    from Bio.Data import CodonTable
    t = CodonTable.unambiguous_dna_by_id[code]
    aa_of = dict(t.forward_table, **{c: "*" for c in t.stop_codons})
    for gid, hs in hsps.items():
        g = lib[gid]
        for h in hs:
            if g.kind == "CDS":
                for u, x in h.cols.items():
                    aa = aa_of.get(G.sub(h.strand, x, x + 3), "X")
                    h.sc[u] = _B62.get((g.protein[u] if u < len(g.protein) else "X", aa), -4.0)
            else:
                for u, t in h.cols.items():
                    h.sc[u] = 2.0 if u < len(g.spliced) and G.sub(h.strand, t, t + 1) == g.spliced[u] else -3.0


def parse_probes(rows: List[List[str]], L: int, ids: Dict[str, List[str]], flank: int = 6
                 ) -> Dict[str, List["ProbeHit"]]:
    """blastn hits of the 60-nt junction probes (30 nt each side of an exon end) -> the query position aligned to the
    probe's junction nt, for hits that align across the junction (+-flank nt)."""
    out: Dict[str, List[ProbeHit]] = {}
    seen = set()
    for qid, qs, qe, ss, se, frame, bits, ev, qseq, sseq in rows:
        strand = 1 if int(se) >= int(ss) else -1
        q, t = int(qs) - 1, (int(ss) - 1 if strand == 1 else 2 * L - int(ss))
        for pid in ids[qid]:
            target = PROBE - 1 if pid.endswith("D") else PROBE      # donor: last exon nt; acceptor: first exon nt
            if not (q <= target - flank and int(qe) - 1 >= target + flank):
                continue
            qq, tt, pos, m, n = q, t, None, 0, 0
            for x, y in zip(qseq, sseq):
                if x != "-" and y != "-":
                    if qq == target:
                        pos = tt
                    m += x == y
                    n += 1
                    qq += 1
                    tt += 1
                elif x == "-":
                    tt += 1
                else:
                    if qq == target:
                        pos = tt - 1 if pid.endswith("D") else tt
                    qq += 1
            if pos is None:
                continue
            pos %= L
            key = (pid, strand, pos)
            if key in seen:
                continue
            seen.add(key)
            out.setdefault(pid, []).append(ProbeHit(strand, pos, float(bits), m / max(n, 1)))
    return out


# ── junction calling ─────────────────────────────────────────────────────────
@dataclass
class JCall:
    """One reference junction mapped onto the query."""
    gid: str
    k: int
    state: str                  # present / absent / unknown
    r: int = 0                  # reference unit of the junction
    splice: str = ""            # cis / trans (in the query)
    sN: int = 0                 # strand and last nt of the upstream exon (strand coordinates)
    d: int = 0
    sC: int = 0                 # strand and first nt of the downstream exon
    a: int = 0
    xN: int = 0                 # anchors: strand coordinate of reference unit iN / iC
    iN: int = 0
    xC: int = 0
    iC: int = 0
    score: float = 0.0
    ident: float = 0.0          # identity of the four junction probes at (d, a)
    ident_ex: float = 0.0       # ... of the two exon arms
    ident_in: float = 0.0       # ... of the two intron arms
    sprot: float = 0.0
    snuc: float = 0.0
    srule: float = 0.0
    why: str = ""


ANCH = {3: 6, 1: 15}            # anchors this many units away from the junction (HSP ends near it are unreliable)
REACH = {3: 25, 1: 60}          # anchors at most this far from the junction (extrapolation over divergent sequence)
NEAR = {3: 8, 1: 20}            # an HSP serves a side only when it aligns a unit this close to the junction
SPAN = {3: 25, 1: 60}           # an HSP aligned contiguously this far on both sides: no intron at this site
MIN_QUERY_INTRON = 100          # shortest intron called in a query (NCBI metazoan mt introns: >= 159 nt)
MIN_IDENT = 0.5                 # junction-probe identity needed to call a site present (benchmark, leave-one-out:
                                # 98 % of false calls below, all true introns but one have a call above)
MIN_IDENT_OTHER_PHYLUM = 0.8    # ... from a reference of another phylum (every false call at 0.5-0.8 in the benchmark
                                # was a placozoan site put on a bilaterian / medusozoan; a recent horizontal transfer
                                # of a homing intron would keep near-identical junctions)
PHYLA = {"Porifera", "Placozoa", "Cnidaria", "Ctenophora", "Annelida", "Mollusca", "Arthropoda", "Chordata",
         "Echinodermata", "Hemichordata", "Nematoda", "Platyhelminthes", "Brachiopoda", "Bryozoa", "Nemertea",
         "Phoronida", "Rotifera", "Acanthocephala", "Chaetognatha", "Tardigrada", "Onychophora", "Priapulida",
         "Kinorhyncha", "Nematomorpha", "Entoprocta", "Cycliophora", "Gastrotricha", "Xenacoelomorpha", "Dicyemida",
         "Orthonectida"}


def phylum_of(lineage: List[str]) -> str:
    return next((t for t in lineage if t in PHYLA), "")


@dataclass
class Side:
    """One candidate end of a query exon at a reference junction: the donor (N side, pos = last exon nt) or the
    acceptor (C side, pos = first exon nt), with an anchor (strand coordinate x of reference unit i) that fixes the
    reading frame, from a protein HSP or derived from a junction-probe hit."""
    strand: int
    pos: int
    x: Optional[int]
    i: Optional[int]
    w: float
    hsp: Optional[Hsp] = None
    probe: bool = False


@dataclass
class ProbeHit:
    strand: int
    pos: int                    # strand coordinate of the probe's junction nt (donor: last exon nt; acceptor: first)
    bits: float
    ident: float


def exon_bounds(ref: RefGene, j: Junction, unit: int) -> Tuple[int, int]:
    """Spliced-sequence offsets of the reference exons around junction j: (start of the upstream exon, end of the
    downstream exon)."""
    offs = sorted(x.offset for x in ref.junctions)
    i = offs.index(j.offset)
    total = len(ref.protein) * 3 if unit == 3 else len(ref.spliced)
    return (offs[i - 1] if i else 0), (offs[i + 1] if i + 1 < len(offs) else total)


def _virtual_anchor(side: str, pos: int, phase: int, r: int, exon_len: int, unit: int) -> Tuple[Optional[int], Optional[int]]:
    """Anchor inside the exon next to a junction found by a probe: the reference's codon phase gives the frame."""
    if unit == 1:
        m = min(ANCH[1], exon_len - 1)
        if m < 1:
            return None, None
        return (pos - m, r - m) if side == "N" else (pos + m, r + m)
    if side == "N":
        m = min(ANCH[3], (exon_len - phase) // 3)
        return (pos + 1 - phase - 3 * m, r - m) if m >= 1 else (None, None)
    first = (3 - phase) % 3                               # nt of the split codon in the downstream exon
    m = min(ANCH[3], (exon_len - first - 3) // 3)
    return (pos + first + 3 * m, r + (1 if phase else 0) + m) if m >= 0 else (None, None)


def _spans(h: Hsp, r: int, unit: int) -> bool:
    """Reference units r-SPAN .. r+SPAN aligned contiguously (no subject insertion) in this HSP."""
    lo, hi = r - SPAN[unit], r + SPAN[unit]
    us = sorted(u for u in h.cols if lo <= u <= hi)
    return bool(us) and us[0] <= lo + 2 and us[-1] >= hi - 2


def _signed_gap(G: Genome, d: int, a: int, unit: int) -> int:
    """nt between a donor's last exon nt d and an acceptor's first nt a on one strand; on a circle a small overlap is
    negative (not an intron around the whole genome): gaps beyond L - 2 REACH codons count as overlaps."""
    if not G.circular:
        return a - d - 1
    g = (a - d - 1) % G.L
    return g - G.L if g > G.L - 2 * unit * REACH[unit] else g


def _near(G: Genome, p: int, q: int, tol: int) -> bool:
    return min((p - q) % G.L, (q - p) % G.L) <= tol if G.circular else abs(p - q) <= tol


def call_junction(G: Genome, ref: RefGene, j: Junction, hsps: List[Hsp], code: int, unit: int,
                  probes: Optional[Dict[str, List[ProbeHit]]] = None) -> JCall:
    """Map one reference junction onto the query: present (with exon ends), absent or unknown."""
    r, phase = (j.offset // 3, j.offset % 3) if unit == 3 else (j.offset, 0)
    c = JCall(ref.gid, j.k, "unknown", r)
    lo_off, hi_off = exon_bounds(ref, j, unit)
    up_len, down_len = j.offset - lo_off, hi_off - j.offset          # reference exon lengths (nt)
    lo_u = (lo_off + unit - 1) // unit                               # first whole unit of the upstream exon
    hi_u = hi_off // unit                                            # one past the last unit of the downstream exon
    ns: List[Side] = []
    cs: List[Side] = []
    for h in hsps:
        # anchors only in the two reference exons at this junction (another exon's HSP would be extrapolated over an
        # intron), weighted by the HSP's part in that exon: an HSP often runs a few units past an exon end into the
        # intron, which must not make it serve the next junction too (Placozoa cox1: 11-nt exon after a 720-nt one)
        up = sorted(u for u in h.cols if max(lo_u, r - REACH[unit]) <= u <= r - 1)
        if up and up[-1] >= min(r - NEAR[unit], lo_u):
            far = [u for u in up if u <= r - ANCH[unit]]
            i = far[-1] if far else up[-1]
            x = h.cols[i]
            w = h.seg(max(lo_u, r - REACH[unit]), r)
            if w > 0:
                ns.append(Side(h.strand, x + unit * (r - i) + phase - 1, x, i, w, h))
        first = r + (1 if phase else 0)
        dn = sorted(u for u in h.cols if first <= u <= min(hi_u - 1, r + REACH[unit]))
        if dn and dn[0] <= max(first + NEAR[unit], hi_u - 1):
            far = [u for u in dn if u >= r + ANCH[unit]]
            i = far[0] if far else dn[0]
            x = h.cols[i]
            w = h.seg(r, min(hi_u, r + REACH[unit]))
            if w > 0:
                cs.append(Side(h.strand, x - unit * (i - r) + phase, x, i, w, h))
    for kind, lst, ln in (("D", ns, up_len), ("A", cs, down_len)):
        for hit in (probes or {}).get(f"{ref.gid}#{j.k}{kind}", []):
            same = [s_ for s_ in lst if s_.strand == hit.strand and _near(G, s_.pos, hit.pos, 12)]
            if same:
                for s_ in same:
                    s_.probe, s_.w = True, s_.w + hit.bits
                continue
            x, i = _virtual_anchor("N" if kind == "D" else "C", hit.pos, phase, r, ln, unit)
            lst.append(Side(hit.strand, hit.pos, x, i, hit.bits, None, True))
    best = None
    pairs = []
    adjacent = 0.0              # evidence of absence: both sides found, contiguous on one strand
    K = REACH[unit]
    spans = {id(h): _spans(h, r, unit) for h in hsps}
    for n_ in ns:
        for c_ in cs:
            if n_.hsp is not None and n_.hsp is c_.hsp:
                continue
            # a side taken from an HSP that runs on across the site must beat that HSP's own continuation (the
            # nad4 paralogue's start joined to nad5's end; an intron start that mimics the next exon, Ricordea cox1)
            if n_.hsp is not None and spans[id(n_.hsp)] and c_.w <= n_.hsp.seg(r, r + K) + 10:
                continue
            if c_.hsp is not None and spans[id(c_.hsp)] and n_.w <= c_.hsp.seg(r - K, r) + 10:
                continue
            if n_.strand == c_.strand:
                gap = _signed_gap(G, n_.pos, c_.pos, unit)
                if gap < MIN_QUERY_INTRON:
                    # two alignment pieces next to each other: the gene runs on across the site (a junction probe
                    # can also hit inside the intron, near a duplicated exon end: no evidence of absence)
                    if gap >= -3 * unit * NEAR[unit] and n_.hsp is not None and c_.hsp is not None:
                        adjacent = max(adjacent, n_.w + c_.w)
                    continue
                kind = "cis"
            elif j.splice == "trans":
                kind = "trans"          # a query junction is called trans-spliced only where the reference's is
            else:
                continue
            val = n_.w + c_.w
            pairs.append((val, n_, c_, kind))
            if best is None or val > best[0]:
                best = (val, n_, c_, kind)
    # absent: an HSP aligned contiguously across the site (not one of the pair), or both sides adjacent, as strong as
    # the pair - a paralogue (nad4 / nad5) spanning the site is much weaker
    spanning = [h for h in hsps if spans[id(h)] and (best is None or (h is not best[1].hsp and h is not best[2].hsp))]
    span_w = max([h.seg(r - K, r + K) for h in spanning], default=0.0)
    # adjacent pieces share a side with the intron pair, so they must be at least as strong as it (an actinarian
    # cox1 intron starts with a 98-residue copy of exon 2 in frame with exon 1: NC_022470.1, NC_047221.1 ...)
    if (span_w > 0 and (best is None or span_w >= 0.5 * best[0])) or \
            (adjacent > 0 and (best is None or adjacent >= best[0])):
        c.state, c.why = "absent", ("reference aligned contiguously across the site" if spanning
                                    else "both sides found next to each other")
        return c
    if best is None:
        c.why = "no alignment on " + ("either side" if not ns and not cs else "the upstream side" if not ns
                                      else "the downstream side") + " of the site"
        return c
    # the three strongest pairs are refined and the best full score wins: an exact junction-probe hit beats an
    # alignment extrapolated over an intron the reference lacks (a site next to a mini-exon)
    tried, done = [], set()
    for _, n_, c_, kind in sorted(pairs, key=lambda x: -x[0]):
        key = (n_.strand, n_.pos, c_.strand, c_.pos)
        if key in done:
            continue
        done.add(key)
        t = JCall(ref.gid, j.k, "present", r, kind, n_.strand, 0, c_.strand, 0, n_.x, n_.i, c_.x, c_.i)
        refine(G, t, ref, j, n_, c_, code, unit)
        if t.state == "present":
            tried.append(t)
        if len(done) == 3:
            break
    if not tried:
        c.why = "no stop-free, in-frame junction near the extrapolated site"
        return c
    c = max(tried, key=lambda t: t.score)
    if c.state == "present" and c.ident < MIN_IDENT:
        c.state, c.why = "weak", f"junction identity {c.ident:.2f} < {MIN_IDENT}: no recognisably homologous junction"
    return c


# position weights of the junction probes, nearest nt first
_PW = [1.0] * 10 + [0.5] * 10 + [0.25] * 10


def _probe_score(q: str, ref: str, outward_left: bool) -> Tuple[float, int, int]:
    """Weighted +1/-1 identity of two equal-length arms counted from the junction outwards."""
    if outward_left:
        q, ref = q[::-1], ref[::-1]
    sc, m, n = 0.0, 0, 0
    for i, (x, y) in enumerate(zip(q, ref)):
        if x == "N" or y == "N":
            continue
        n += 1
        m += x == y
        sc += _PW[i] if x == y else -_PW[i]
    return sc, m, n


def rule_score(itype: str, G: Genome, sN: int, d: int, sC: int, a: int) -> float:
    """Splice rules: group I - exon-terminal U, intron-terminal omega-G; group II - GUGYG ... AY."""
    if itype == "I":
        return 2.0 * (G.sub(sN, d, d + 1) == "T") + 3.0 * (G.sub(sC, a - 1, a) == "G")
    if itype == "II":
        five = G.sub(sN, d + 1, d + 6)
        sc = sum(w for w, x, ok in zip((2, 1, 2, 1, 2), five, ("G", "T", "G", "CT", "G")) if x in ok)
        three = G.sub(sC, a - 2, a)
        return sc + (2.0 * (three[:1] == "A") + 2.0 * (three[1:] in ("C", "T")) if len(three) == 2 else 0)
    return 0.0


_ALN: Dict[str, object] = {}


def _aligner(kind: str):
    if kind not in _ALN:
        from Bio.Align import PairwiseAligner, substitution_matrices
        if kind == "prot":
            al = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"))
            al.open_gap_score, al.extend_gap_score = -11, -1
        else:
            al = PairwiseAligner(mode="global", match_score=2, mismatch_score=-3)
            al.open_gap_score, al.extend_gap_score = -5, -2
        _ALN[kind] = al
    return _ALN[kind]


def refine(G: Genome, c: JCall, ref: RefGene, j: Junction, n_: Side, c_: Side, code: int, unit: int) -> None:
    """Choose the exon ends (d, a) near the two sides: reading frame kept between the anchors, no stop codon, then the
    best sum of (1) alignment of the spliced query across the junction to the reference, (2) identity of the query's
    exon / intron junctions to the reference's (its annotated boundary, transferred), (3) splice rules of the intron
    group. A side found only by a junction probe is searched in a narrow window; without an anchor on one side (exon
    shorter than a codon) there is no frame test here (the whole spliced CDS is tested later)."""
    def width(s_: Side, dist: int) -> int:
        return 6 if s_.hsp is None else min(45, 15 + unit * max(0, dist - ANCH[unit]))
    WN = width(n_, c.r - c.iN if c.iN is not None else 0)
    WC = width(c_, c.iC - c.r if c.iC is not None else 0)
    framed = c.xN is not None and c.xC is not None
    refseq = ref.protein if unit == 3 else ref.spliced
    target = refseq[c.iN:c.iC] if framed else ""
    al = _aligner("prot" if unit == 3 else "nt")
    cache: Dict[str, float] = {}
    best = None
    d0, a0 = n_.pos, c_.pos
    for d in range(d0 - WN, d0 + WN + 1):
        if c.xN is not None and d < c.xN + unit - 1:
            continue
        up = G.sub(c.sN, c.xN, d + 1) if framed else ""
        e5 = G.sub(c.sN, d - 29, d + 1)
        i5 = G.sub(c.sN, d + 1, d + 31)
        p5 = _probe_score(e5, j.e5[-len(e5):] if e5 else "", True)
        p5i = _probe_score(i5, j.i5[:len(i5)], False)
        for a in range(a0 - WC, a0 + WC + 1):
            if c.xC is not None and a > c.xC:
                continue
            if c.sN == c.sC and _signed_gap(G, d, a, unit) < MIN_QUERY_INTRON:
                continue
            sp = 0.0
            if framed:
                nt = up + G.sub(c.sC, a, c.xC)
                if unit == 3:
                    if len(nt) % 3:
                        continue
                    prot = str(Seq(nt).translate(table=code))
                    if "*" in prot:
                        continue
                else:
                    prot = nt
                cached = cache.get(prot)
                if cached is None:
                    cached = cache[prot] = al.score(prot, target) if prot and target else 0.0
                sp = cached
            i3 = G.sub(c.sC, a - 30, a)
            e3 = G.sub(c.sC, a, a + 30)
            p3i = _probe_score(i3, j.i3[-len(i3):] if i3 else "", True)
            p3 = _probe_score(e3, j.e3[:len(e3)], False)
            sn = p5[0] + p5i[0] + p3i[0] + p3[0]
            sr = rule_score(j.itype, G, c.sN, d, c.sC, a)
            tot = sp + sn + sr
            key = (tot, -(abs(d - d0) + abs(a - a0)))
            if best is None or key > best[0]:
                m = p5[1] + p5i[1] + p3i[1] + p3[1]
                n = p5[2] + p5i[2] + p3i[2] + p3[2]
                best = (key, d, a, sp, sn, sr, m / max(n, 1), (p5[1] + p3[1]) / max(1, p5[2] + p3[2]),
                        (p5i[1] + p3i[1]) / max(1, p5i[2] + p3i[2]))
    if best is None:
        c.state, c.why = "unknown", "no stop-free, in-frame junction near the extrapolated site"
        return
    (tot, _), c.d, c.a, c.sprot, c.snuc, c.srule, c.ident, c.ident_ex, c.ident_in = best
    c.score = tot


# ── gene models ──────────────────────────────────────────────────────────────
@dataclass
class IntronCall:
    strand: int
    start: int                  # strand coordinates [start, end)
    end: int
    number: int
    itype: str
    ref: str                    # reference gene id and junction number
    ident: float
    note: str
    half: str = ""              # '', '5' / '3' (half-intron of a trans-spliced junction)
    fixed: bool = False         # placement from the record's own annotation: never moved
    family: bool = False        # site placed by its clade's consensus (family rule): the transferred placement stays
    cm: str = ""                # covariance-model typing of this query intron (the de novo note), if any


@dataclass
class GeneModel:
    gene: str
    kind: str                                   # CDS / rRNA
    exons: List[Tuple[int, int, int]]           # (strand, start, end) strand coordinates, transcript order
    introns: List[IntronCall]
    calls: List[JCall]
    refs: List[str]
    partial5: bool = False
    partial3: bool = False
    confidence: str = "high"
    evidence: str = ""
    protein: str = ""

    def trans(self) -> bool:
        return len({s for s, _, _ in self.exons}) > 1 or any(c.splice == "trans" for c in self.calls)


def ref_depth(q: List[str], g: RefGene) -> int:
    """Closeness of a reference to the query: shared lineage depth; the query's own annotation above everything."""
    return 1000 if g.own else lineage_score(q, g.lineage)


def lineage_score(q: List[str], r: List[str]) -> int:
    n = 0
    for a, b in zip(q, r):
        if a != b:
            break
        n += 1
    return n


def species_of(organism: str) -> str:
    return " ".join(organism.split()[:2]).lower()


def genus_of(organism: str) -> str:
    return organism.split()[0].lower() if organism else ""


def _base(acc: str) -> str:
    return acc.split(".")[0]


@dataclass
class Exclude:
    """References left out (benchmarks): accessions (any version; RefSeq copies via same_as), species, genus."""
    accs: set = field(default_factory=set)
    species: str = ""
    genus: str = ""

    def expand(self, lib: List[RefGene]) -> "Exclude":
        """Add the copies of the excluded records: a RefSeq record names its INSDC source ('identical to ...'), the
        source does not name the RefSeq copy, so both directions are needed (NC_087806.1 / OR123442.1)."""
        base = {_base(a) for a in self.accs}
        for g in lib:
            names = {_base(g.acc), *map(_base, g.same_as)}
            if names & base:
                self.accs |= names
        return self

    def __call__(self, g: RefGene) -> bool:
        names = {_base(g.acc), *map(_base, g.same_as)}
        if names & {_base(a) for a in self.accs}:
            return True
        if self.species and species_of(g.organism) == self.species.lower():
            return True
        return bool(self.genus) and genus_of(g.organism) == self.genus.lower()


def _same_site(x: JCall, y: JCall, L: int) -> bool:
    """Two calls of one query intron: same strands, one end within 30 nt (references annotate the other end of a
    homologous intron up to tens of nt apart)."""
    def near(p, q):
        return min((p - q) % L, (q - p) % L) <= 30
    return x.sN == y.sN and x.sC == y.sC and (near(x.d, y.d) or near(x.a, y.a))


def _codons(code: int):
    from Bio.Data import CodonTable
    t = CodonTable.unambiguous_dna_by_id[code]
    return set(t.start_codons), set(t.stop_codons)


def core_cds(rec, gene: str, G: Genome) -> List[Tuple[int, int, int]]:
    """Existing CDS / rRNA pieces of this gene as (strand, start, end) strand-coordinate intervals."""
    out = []
    for f in rec.features:
        if f.type not in ("CDS", "rRNA") or gene_name(f.qualifiers, f.type) != gene:
            continue
        for p in f.location.parts:
            st = p.strand or 1
            a, b = int(p.start), int(p.end)
            out.append((st, a, b) if st == 1 else (st, G.L - b, G.L - a))
    return out


def gene_start(G: Genome, strand: int, x: int, i: int, unit: int, code: int, core: List[Tuple[int, int, int]]
               ) -> Tuple[int, bool]:
    """First nt of the gene: extrapolated from anchor x (reference unit i) to unit 0; for a CDS an existing in-frame
    start of the gene (annotation core) near it, else the in-frame start codon nearest to it with no stop before x."""
    s0 = x - unit * i
    if unit == 1:
        return s0, (not G.circular and s0 < 0)
    starts, stops = _codons(code)
    for st, a, b in core:
        if st == strand and abs(a - s0) <= 90 and (x - a) % 3 == 0 and G.sub(strand, a, a + 3) in starts \
                and not any(G.sub(strand, p, p + 3) in stops for p in range(a, x, 3)):
            return a, False
    best = None
    p = x - 3
    while p >= s0 - 150:
        if not G.circular and p < 0:
            return (best if best is not None else p + 3), best is None
        cod = G.sub(strand, p, p + 3)
        if cod in stops:
            break
        if cod in starts and (best is None or abs(p - s0) < abs(best - s0)):
            best = p
        p -= 3
    # no start codon where the reference starts (a piece of a trans-spliced gene, a wrong chain): 5'-partial CDS
    return (best, False) if best is not None else (s0, True)


def gene_end(G: Genome, strand: int, x: int, i: int, n: int, unit: int, code: int,
             core: List[Tuple[int, int, int]]) -> Tuple[int, bool]:
    """One past the last nt: for a CDS the first in-frame stop codon after anchor x (reference unit i of n), an
    existing gene end there, or an incomplete T / TA stop where the reference ends; rRNA: extrapolated."""
    e0 = x + unit * (n - i)
    if unit == 1:
        return e0, (not G.circular and e0 > G.L)
    _, stops = _codons(code)
    for st, a, b in core:
        if st == strand and abs(b - e0) <= 90 and (b - x) % 3 == 0 and G.sub(strand, b - 3, b) in stops \
                and not any(G.sub(strand, p, p + 3) in stops for p in range(x, b - 3, 3)):
            return b, False
    p = x
    while p < e0 + 300:
        if not G.circular and p + 3 > G.L:
            return p, True
        if G.sub(strand, p, p + 3) in stops:
            if p >= e0 - 90:
                return p + 3, False
            break
        p += 3
    # no complete stop where the reference ends: incomplete stop (T / TA completed by polyadenylation)
    q = x + 3 * ((e0 - x) // 3)
    tail = G.sub(strand, q, q + 2)
    if tail == "TA":
        return q + 2, False
    if tail[:1] == "T":
        return q + 1, False
    # an early stop (truncated protein) ends the CDS there; no stop at all: 3'-partial CDS (table2asn NoStop otherwise)
    return (p + 3, False) if p < e0 + 300 else (q, True)


MAX_EXON = 2500         # longest exon between two sites of one chain (mt protein / rRNA genes are < 2 kb)


def vote(group: List[JCall], lib: Dict[str, RefGene], q_lin: List[str], rank: Dict[str, int]) -> JCall:
    """Exon ends of one query site: the (donor, acceptor) pair with the largest support, each reference counting its
    junction identity x 2^(lineage depth - deepest depth in the group) - the boundary most close relatives annotate,
    not one relative's idiosyncratic one (Diadumene lineata cox1 / nad5 ends 22-24 nt from the other sea anemones;
    Myriopathes ulex cox1 3' end 3 nt from the other black corals). The query's own annotation outweighs all."""
    depth = {c.gid: ref_depth(q_lin, lib[c.gid]) for c in group}
    top = max(depth.values())
    support: Dict[Tuple[int, int, int, int], float] = {}
    for c in group:
        key = (c.sN, c.d, c.sC, c.a)
        support[key] = support.get(key, 0.0) + c.ident * 2.0 ** (depth[c.gid] - top)
    win = max(support, key=lambda k: support[k])
    return max((c for c in group if (c.sN, c.d, c.sC, c.a) == win),
               key=lambda c: (depth[c.gid], c.ident, -rank[c.gid]))


def assemble(G: Genome, gene: str, kind: str, ranked: List[RefGene], calls: Dict[str, List[JCall]],
             q_lin: List[str], code: int, rec) -> Tuple[List[GeneModel], List[str]]:
    """Present sites of all ranked references -> gene models. One call per query site (closest lineage, then most
    identical junction, then reference rank); sites chained in transcript order when the exon between two of them
    lies on one strand with a plausible length - separately annotated parts of a gene (Trichoplax COX1a / COX1b)
    become separate models."""
    unit = 3 if kind == "CDS" else 1
    rank = {g.gid: n for n, g in enumerate(ranked)}
    lib = {g.gid: g for g in ranked}
    groups: List[List[JCall]] = []
    for g in ranked:
        for c in calls.get(g.gid, []):
            if c.state != "present":
                continue
            grp = next((x for x in groups if _same_site(x[0], c, G.L)), None)
            if grp is None:
                groups.append([c])
            else:
                grp.append(c)
    chosen = [vote(grp, lib, q_lin, rank) for grp in groups]

    def rel(c: JCall) -> float:                         # junction position relative to the reference length
        g = lib[c.gid]
        return c.r / max(1, len(g.protein) if unit == 3 else len(g.spliced))
    chains: List[List[JCall]] = []
    for c in sorted(chosen, key=rel):
        best = None
        for ch in chains:
            last = ch[-1]
            if last.sC != c.sN:
                continue
            ln = (c.d - last.a) % G.L if G.circular else c.d - last.a
            if ln < 0 or ln > MAX_EXON:
                continue
            pref = (last.gid == c.gid, -ln)
            if best is None or pref > best[0]:
                best = (pref, ch)
        if best is None:
            chains.append([c])
        else:
            best[1].append(c)
    models: List[GeneModel] = []
    whys: List[str] = []
    for ch in chains:
        m, why = gene_model(G, gene, kind, ch, lib, code, rec, unit)
        if m is not None:
            models.append(m)
        if why:
            whys.append(why)
    # two chains of one gene over the same exons (a weak site joining a paralogue piece): the better supported one
    models.sort(key=lambda m: (m.confidence != "low", sum(c.ident for c in m.calls) / len(m.calls), len(m.calls)),
                reverse=True)
    kept: List[GeneModel] = []
    for m in models:
        span = [(int(p.start), int(p.end)) for p in _location(G, m.exons).parts]
        clash = next((k for k in kept if any(st == ks for st, _, _ in m.exons for ks, _, _ in k.exons)
                      and _overlaps(_location(G, k.exons), span) > 0), None)
        if clash is not None:
            whys.append(f"{gene}: model via {', '.join(m.refs)} overlaps the better supported one via "
                        f"{', '.join(clash.refs)} - not written")
            continue
        kept.append(m)
    return kept, whys


def _end_anchor(c: JCall, ref: RefGene, side: str, unit: int) -> Tuple[int, int]:
    """(x, i) anchor for extrapolating a gene end; a junction found only by its probe (no anchor) uses the
    reference phase at the junction itself."""
    j = next(x for x in ref.junctions if x.k == c.k)
    if side == "N":
        if c.xN is not None:
            return c.xN, c.iN
        return (c.d + 1 - j.offset % 3, j.offset // 3) if unit == 3 else (c.d + 1, j.offset)
    if c.xC is not None:
        return c.xC, c.iC
    if unit == 3:
        first = (3 - j.offset % 3) % 3
        return c.a + first, j.offset // 3 + (1 if j.offset % 3 else 0)
    return c.a, j.offset


GROUP_II_5, GROUP_II_3 = re.compile(r"^G[CT]G[CT]G"), re.compile(r"A[CT]$")
CM_SURE = 1e-5          # a covariance-model hit this good fixes the intron group
SYN_SHIFT = 60          # shifts outside the sliding window tried for a synonymous canonical placement


def canonical_kinds(last: str, i5: str, i3: str) -> List[str]:
    """Rules a placement meets (gold standard, intron_gold_2026-10-29/gold/scripts/verify.py): I - the 5' exon ends T
    and the intron ends G (omega-G); II - the intron starts GYGYG (GUGYG, GCGCG of the annelid cox1 introns) and ends
    AY; II5 - the group II 5' motif alone."""
    out = []
    if last == "T" and i3.endswith("G"):
        out.append("I")
    if GROUP_II_5.match(i5):
        out.append("II" if GROUP_II_3.search(i3) else "II5")
    return out


def meets(kinds: List[str], allowed: set) -> str:
    """'full' when a placement meets the rule of an allowed group, '5p' for the group II 5' motif alone, else ''."""
    if any(k in allowed for k in kinds if k in ("I", "II")):
        return "full"
    return "5p" if "II5" in kinds and "II" in allowed else ""


def allowed_groups(cm: str, ref_kind: str) -> set:
    """Groups whose rule may place an intron (gold classify.py rules()): the intron's own covariance-model group when
    E <= 1e-5; a weaker hit, the group of the reference intron (NCBI's note in the gold standard) or either."""
    kind = intron_kind(cm)
    m = re.search(r"E=([0-9.]+(?:e[-+]?\d+)?)", cm)
    ev = float(m.group(1)) if m else None
    if kind and ev is not None and ev <= CM_SURE:
        return {kind}
    if kind:
        return {kind, ref_kind} if ref_kind else {"I", "II"}
    return {ref_kind} if ref_kind else {"I", "II"}


def gold_shift(G: Genome, exons: List[Tuple[int, int, int]], k: int, allowed: set, code: int, coding: bool
               ) -> Tuple[int, str]:
    """Shift (transcript direction, both ends) that puts the cis intron between exons k and k+1 where the gold standard
    puts NCBI's introns: as it is when it meets a rule; else the canonical placement nearest to it inside the sliding
    window (identical spliced sequence; on a tie the lower genome coordinate); else the nearest shift <= 60 nt that
    keeps the protein identical (synonymous change) and meets a rule; then the same for the group II 5' motif alone;
    else (gold class C) as it is. Homing introns sit in short direct repeats and submitters place them anywhere in
    the repeat; the rules make the placement independent of the reference's choice."""
    st, d, a = exons[k][0], exons[k][2] - 1, exons[k + 1][1]
    ln = _signed_gap(G, d, a, 3)
    if ln < 50:
        return 0, ""

    def rule(x: int) -> str:
        return meets(canonical_kinds(G.sub(st, d + x, d + x + 1), G.sub(st, d + x + 1, d + x + 7),
                                     G.sub(st, a + x - 6, a + x)), allowed)
    here = rule(0)
    if here == "full":
        return 0, ""
    kmax = 0
    while kmax < ln - 1 and G.sub(st, d + 1 + kmax, d + 2 + kmax) == G.sub(st, a + kmax, a + kmax + 1):
        kmax += 1
    kmin = 0
    while -kmin < ln - 1 and G.sub(st, d + kmin, d + kmin + 1) == G.sub(st, a - 1 + kmin, a + kmin):
        kmin -= 1
    win = {x: rule(x) for x in range(kmin, kmax + 1)}
    prot0 = None

    def synonymous(x: int) -> bool:
        nonlocal prot0
        if not coding:
            return False
        ex = list(exons)
        ex[k], ex[k + 1] = (st, ex[k][1], ex[k][2] + x), (st, ex[k + 1][1] + x, ex[k + 1][2])
        if ex[k][2] <= ex[k][1] or ex[k + 1][2] <= ex[k + 1][1]:
            return False

        def prot(e):
            nt = "".join(G.sub(*y) for y in e)
            return codon_translate(nt, code)
        if prot0 is None:
            prot0 = prot(exons)
        p2 = prot(ex)
        return p2 == prot0 and "*" not in p2.rstrip("*")
    order = sorted(win, key=lambda x: (abs(x), x * st))                  # tie: lower genome coordinate
    syn_order = sorted((x for x in range(-SYN_SHIFT, SYN_SHIFT + 1) if not kmin <= x <= kmax),
                       key=lambda x: (abs(x), x))
    for want in ("full", "5p"):
        if want == "5p" and here == "5p":
            return 0, ""
        x = next((x for x in order if x and win[x] == want), None)
        if x is not None:
            return x, f"same spliced sequence, {'canonical' if want == 'full' else 'group II 5-prime motif'} {x:+d} nt"
        x = next((x for x in syn_order if ln + 0 >= 50 and rule(x) == want and synonymous(x)), None)
        if x is not None:
            return x, f"same protein (synonymous), {'canonical' if want == 'full' else 'group II 5-prime motif'} {x:+d} nt"
    return 0, ""


def place_introns(G: Genome, m: GeneModel, code: int) -> None:
    """Gold-standard placement of every cis intron of m not taken from the record's own annotation."""
    for i in m.introns:
        if i.fixed or i.half or i.family:
            continue
        k = next((n for n in range(len(m.exons) - 1) if m.exons[n][0] == m.exons[n + 1][0] == i.strand
                  and _near(G, m.exons[n][2], i.start, 0) and _near(G, m.exons[n + 1][1], i.end, 0)), None)
        if k is None:
            continue
        x, why = gold_shift(G, m.exons, k, allowed_groups(i.cm, i.itype), code, m.kind == "CDS")
        if x:
            st = i.strand
            m.exons[k] = G.norm(st, m.exons[k][1], m.exons[k][2] + x) or m.exons[k]
            m.exons[k + 1] = G.norm(st, m.exons[k + 1][1] + x, m.exons[k + 1][2]) or m.exons[k + 1]
            i.start, i.end = G.norm(st, i.start + x, i.end + x)[1:]
            i.note += f"; placed by the gold-standard rules: {why} from the {'de novo' if i.ref == 'de novo' else 'transferred'} placement"


def gene_model(G: Genome, gene: str, kind: str, chain: List[JCall], lib: Dict[str, RefGene], code: int, rec,
               unit: int) -> Tuple[Optional[GeneModel], str]:
    core = core_cds(rec, gene, G)
    first, last = chain[0], chain[-1]
    xs, is_ = _end_anchor(first, lib[first.gid], "N", unit)
    s, p5 = gene_start(G, first.sN, xs, is_, unit, code, core)
    gl = lib[last.gid]
    xe, ie = _end_anchor(last, gl, "C", unit)
    e, p3 = gene_end(G, last.sC, xe, ie, len(gl.protein) if unit == 3 else len(gl.spliced), unit, code, core)
    exons = [(first.sN, s, first.d + 1)]
    for x, y in zip(chain, chain[1:]):
        exons.append((x.sC, x.a, y.d + 1))
    exons.append((last.sC, last.a, e))
    exons = [G.norm(st, a, b) for st, a, b in exons]
    if any(x is None for x in exons):
        return None, f"{gene}: empty or impossible exon"
    m = GeneModel(gene, kind, exons, [], chain, sorted({lib[c.gid].acc for c in chain}), p5, p3)
    for n, c in enumerate(chain, 1):
        g = lib[c.gid]
        j = next(x for x in g.junctions if x.k == c.k)
        # the input record's own annotation (an annotated GenBank file given to the pipeline): keep its quirks too
        # (intron feature off its CDS gap, half-introns of trans-spliced junctions); never another record's, even of
        # the same species (Decemunciger KY742027.1 intron 1 nt short of its CDS gap, KY774370.1 not)
        same = c.ident >= 0.99 and (g.own or _base(rec.id) in {_base(g.acc), *map(_base, g.same_as)})
        itype = j.itype
        src = f"{g.acc} {g.gene} intron {j.k} ({g.organism})"
        note = (f"{'group ' + itype + ' intron' if itype else 'intron'}; site transferred from {src}, "
                f"junction identity {c.ident:.2f}")
        if c.splice == "cis":
            st, a, b = c.sN, c.d + 1, c.a
            if same and j.splice == "cis":
                a, b = a + j.d5, b + j.d3
            st, a, b = G.norm(st, a, b) or (st, a, b)
            # an identical junction (the record itself under another name, a conspecific) keeps the reference's
            # placement: the gold standard already decided it for that very sequence
            m.introns.append(IntronCall(st, a, b, n, itype, f"{c.gid}#{c.k}", c.ident, note,
                                        fixed=same or c.ident >= IDENTICAL, family=j.family))
        else:
            if same and j.half5:
                a = c.d + 1 + j.gap5
                st, a, b = G.norm(c.sN, a, a + j.half5) or (c.sN, a, a + j.half5)
                m.introns.append(IntronCall(st, a, b, n, itype, f"{c.gid}#{c.k}", c.ident,
                                            note + "; 5' part of a trans-spliced intron", "5", fixed=True))
            if same and j.half3:
                b = c.a - j.gap3
                st, a, b = G.norm(c.sC, b - j.half3, b) or (c.sC, b - j.half3, b)
                m.introns.append(IntronCall(st, a, b, n, itype, f"{c.gid}#{c.k}", c.ident,
                                            note + "; 3' part of a trans-spliced intron", "3", fixed=True))
    if unit == 3:
        nt = "".join(G.sub(st, a, b) for st, a, b in exons)
        prot = codon_translate(nt, code)
        m.protein = prot.rstrip("*")
        starts, _ = _codons(code)
        if "*" in m.protein:
            m.confidence, m.evidence = "low", f"internal stop codon in the spliced CDS (residue {m.protein.index('*') + 1})"
        elif not p5 and nt[:3] not in starts:
            m.confidence, m.evidence = "low", f"no start codon ({nt[:3]})"
    if m.confidence != "low":
        weak = [c for c in chain if c.ident < 0.6]
        m.confidence = "medium" if weak else "high"
        m.evidence = (f"{len(chain)} site(s) transferred; junction identity "
                      + ", ".join(f"{c.ident:.2f}" for c in chain))
    return m, ""


# ── driver ───────────────────────────────────────────────────────────────────
TOP_REFS = 8            # references per gene (closest lineage, then best alignment) ...
MAX_REFS = 24           # ... plus further ones that bring a site not yet seen, up to this many


def select_refs(cands: List[RefGene], hsps: Dict[str, List[Hsp]], q_lin: List[str]) -> List[RefGene]:
    def key(g: RefGene):
        return (-ref_depth(q_lin, g), -sum(h.bits for h in hsps.get(g.gid, [])))
    out: List[RefGene] = []
    seen: set = set()
    for g in sorted(cands, key=key):
        n = max(1, len(g.protein) if g.kind == "CDS" else len(g.spliced))
        unit = 3 if g.kind == "CDS" else 1
        sig = {(round(j.offset / unit / n * 50)) for j in g.junctions}
        if len(out) < TOP_REFS or (sig - seen and len(out) < MAX_REFS):
            out.append(g)
            seen |= sig
    return out


def find_tools(bin_dir: Optional[str]) -> Dict[str, Optional[str]]:
    return {t: find_tool(t, bin_dir) for t in ("tblastn", "blastn")}


def annotate_record(rec, lib: List[RefGene], code: int, lineage: Optional[List[str]] = None,
                    exclude: Optional[Exclude] = None, tools: Optional[Dict[str, Optional[str]]] = None,
                    genes: Optional[List[str]] = None, own: Optional[List[RefGene]] = None
                    ) -> Tuple[List[GeneModel], List[dict]]:
    """Map every library intron site onto one record (SeqRecord; circular unless annotated linear). code: genetic
    code; lineage: the query's NCBI lineage (ranks references, phylum rule); exclude: references to leave out
    (benchmarks); own: sites of the record's own annotation (ref_genes of the annotated GenBank record given to the
    pipeline) - closest reference of all, their annotation quirks kept. Returns gene models (genes with >= 1 intron)
    and report rows (one per reference junction tried). The record is not changed (see apply_models)."""
    import tempfile
    tools = tools or find_tools(None)
    G = Genome(str(rec.seq), rec.annotations.get("topology", "circular") != "linear")
    q_lin = list(lineage or rec.annotations.get("taxonomy") or [])
    q_phylum = phylum_of(q_lin)
    for g in own or []:
        g.own = True
    # the query is itself a library record (an NCBI input re-annotated under its accession): its entries are its own
    # annotation (copies, the shared library is not changed)
    me = {_base(rec.id)}
    refs = [dataclasses.replace(g, own=True) if not g.own and {_base(g.acc), *map(_base, g.same_as)} & me else g
            for g in list(own or []) + list(lib) if (g.own or not (exclude and exclude(g)))
            and (not genes or g.gene in genes)]
    rows: List[dict] = []
    if not refs or G.L < 100:
        return [], rows
    hsps: Dict[str, List[Hsp]] = {}
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        subj = tdp / "g.fna"
        subj.write_text(">g\n" + (G.seq + G.seq if G.circular else G.seq) + "\n")
        for kind, unit, exe in (("CDS", 3, tools.get("tblastn")), ("rRNA", 1, tools.get("blastn"))):
            seqs: Dict[str, str] = {}
            ids: Dict[str, List[str]] = {}
            for g in refs:
                if g.kind != kind:
                    continue
                s = g.protein if kind == "CDS" else g.spliced
                if not s:
                    continue
                qid = seqs.setdefault(s, f"q{len(seqs)}")
                ids.setdefault(qid, []).append(g.gid)
            if not seqs or not exe:
                continue
            qf = tdp / f"{kind}.fa"
            qf.write_text("".join(f">{q}\n{s}\n" for s, q in seqs.items()))
            extra = (["-db_gencode", str(code), "-seg", "no", "-comp_based_stats", "0", "-evalue", "1"]
                     if kind == "CDS" else ["-task", "blastn", "-dust", "no", "-evalue", "1e-5"])
            hsps.update(parse_hsps(_blast(exe, qf, subj, extra), G.L, unit, ids))
        # junction probes (60 nt across each annotated exon end): find exon ends at the nucleotide level, needed for
        # exons too short for a protein alignment (Placozoa cox1: 11-nt and 42-nt exons) and for close relatives
        pseqs: Dict[str, str] = {}
        pids: Dict[str, List[str]] = {}
        for g in refs:
            for j in g.junctions:
                for side, sq in (("D", j.e5 + j.i5), ("A", j.i3 + j.e3)):
                    if len(sq) == 2 * PROBE and "N" not in sq:
                        qid = pseqs.setdefault(sq, f"p{len(pseqs)}")
                        pids.setdefault(qid, []).append(f"{g.gid}#{j.k}{side}")
        probes: Dict[str, List[ProbeHit]] = {}
        if pseqs and tools.get("blastn"):
            pf = tdp / "probes.fa"
            pf.write_text("".join(f">{q}\n{sq}\n" for sq, q in pseqs.items()))
            prow = _blast(tools["blastn"], pf, subj, ["-task", "blastn-short", "-dust", "no", "-evalue", "1e-4",
                                                     "-max_hsps", "5"])
            probes = parse_probes(prow, G.L, pids)
    score_hsps(G, hsps, {g.gid: g for g in refs}, code)
    if not q_phylum and hsps:
        # no lineage given: the phylum of the closest library relative (best total alignment score) stands in for
        # the query's, so that sites of other phyla still need near-identical junctions (benchmark without lineage:
        # Trichoplax sites on a scyphozoan and a sponge otherwise)
        best = max(refs, key=lambda g: sum(h.bits for h in hsps.get(g.gid, [])))
        q_phylum = phylum_of(best.lineage)
    models: List[GeneModel] = []
    for gene in sorted({g.gene for g in refs}):
        cands = [g for g in refs if g.gene == gene and g.gid in hsps]
        if not cands:
            continue
        kind = cands[0].kind
        unit = 3 if kind == "CDS" else 1
        ranked = select_refs(cands, hsps, q_lin)
        calls: Dict[str, List[JCall]] = {}
        for g in ranked:
            calls[g.gid] = [call_junction(G, g, j, hsps[g.gid], code, unit, probes) for j in g.junctions]
            other = q_phylum and phylum_of(g.lineage) not in ("", q_phylum)
            for c in calls[g.gid]:
                if other and c.state == "present" and c.ident < MIN_IDENT_OTHER_PHYLUM:
                    c.state, c.why = "weak", (f"reference of another phylum ({phylum_of(g.lineage)}) with junction "
                                              f"identity {c.ident:.2f} < {MIN_IDENT_OTHER_PHYLUM}")
            for c in calls[g.gid]:
                rows.append({"record": rec.id, "gene": gene, "ref": g.gid, "ref_organism": g.organism, "junction": c.k,
                             "lineage_shared": ref_depth(q_lin, g), "state": c.state, "splice": c.splice,
                             "donor": _gpos(G, c.sN, c.d) if c.state in ("present", "weak") else "",
                             "acceptor": _gpos(G, c.sC, c.a) if c.state in ("present", "weak") else "",
                             "ident": f"{c.ident:.3f}", "ident_ex": f"{c.ident_ex:.3f}", "ident_in": f"{c.ident_in:.3f}",
                             "score": f"{c.score:.1f}", "why": c.why})
        ms, whys = assemble(G, gene, kind, ranked, calls, q_lin, code, rec)
        models += ms
        for why in whys:
            rows.append({"record": rec.id, "gene": gene, "state": "rejected", "why": why})
    return models, rows


def _gpos(G: Genome, strand: int, c: int) -> str:
    """1-based genomic position of strand coordinate c, with the strand sign."""
    g = c % G.L if strand == 1 else G.L - 1 - (c % G.L)
    return f"{'+' if strand == 1 else '-'}{g + 1}"


# ── write-back ───────────────────────────────────────────────────────────────
PRODUCTS = mito_names.PRODUCTS


def _location(G: Genome, pieces: List[Tuple[int, int, int]], partial5: bool = False, partial3: bool = False):
    """(strand, start, end) strand-coordinate pieces in transcript order -> Biopython location (origin split)."""
    from Bio.SeqFeature import (AfterPosition, BeforePosition, CompoundLocation, ExactPosition, FeatureLocation)
    locs = []
    for st, a, b in pieces:
        for x, y in G.genomic(st, a, b):
            locs.append(FeatureLocation(ExactPosition(x), ExactPosition(y), strand=st))
    if partial5:
        f = locs[0]
        locs[0] = (FeatureLocation(BeforePosition(f.start), f.end, strand=f.strand) if f.strand == 1
                   else FeatureLocation(f.start, AfterPosition(f.end), strand=f.strand))
    if partial3:
        f = locs[-1]
        locs[-1] = (FeatureLocation(f.start, AfterPosition(f.end), strand=f.strand) if f.strand == 1
                    else FeatureLocation(BeforePosition(f.start), f.end, strand=f.strand))
    return locs[0] if len(locs) == 1 else CompoundLocation(locs)


def _span_text(loc) -> str:
    return ",".join(f"{'complement(' if p.strand == -1 else ''}{int(p.start) + 1}..{int(p.end)}"
                    f"{')' if p.strand == -1 else ''}" for p in loc.parts)


def _overlaps(loc, parts: List[Tuple[int, int]]) -> int:
    return sum(max(0, min(int(p.end), b) - max(int(p.start), a)) for p in loc.parts for a, b in parts)


def _spans_of(G: Genome, exons: List[Tuple[int, int, int]]) -> List[Tuple[int, int, int]]:
    """One span per run of exons on one strand (the gene feature; NCBI style for trans-spliced genes)."""
    spans: List[Tuple[int, int, int]] = []
    for st, a, b in exons:
        if spans and spans[-1][0] == st and ((a - spans[-1][2]) % G.L if G.circular else a - spans[-1][2]) >= 0:
            a0 = spans[-1][1]
            spans[-1] = (st, a0, a0 + ((b - a0) % G.L if G.circular else b - a0))
        else:
            spans.append((st, a, b))
    return spans


def _existing_models(G: Genome, rec, m: GeneModel) -> list:
    """The joined CDS / rRNA of the same gene already in the record (the de novo step, mito_intron.py) that overlap
    m, largest overlap first, each as (feature, exons in strand coordinates and transcript order)."""
    span = [(int(p.start), int(p.end)) for st, a, b in m.exons for p in _location(G, [(st, a, b)]).parts]
    feats = [f for f in rec.features if f.type == m.kind and gene_name(f.qualifiers, m.kind) == m.gene
             and len(f.location.parts) > 1 and _overlaps(f.location, span) > 0]
    return [_tx_model(G, f) for f in sorted(feats, key=lambda f: -_overlaps(f.location, span))]


def _tx_model(G: Genome, f):
    """(feature, exons in strand coordinates and transcript order; parts split at the origin joined)."""
    ex: List[Tuple[int, int, int]] = []
    for p in f.location.parts:
        st, a, b = p.strand or 1, int(p.start), int(p.end)
        st, a, b = (st, a, b) if st == 1 else (st, G.L - b, G.L - a)
        if ex and ex[-1][0] == st and ex[-1][2] % G.L == a:              # a part continued across the origin
            ex[-1] = (st, ex[-1][1], ex[-1][2] + (b - a))
        else:
            ex.append(G.norm(st, a, b) or (st, a, b))
    return f, ex


def _intron_note(G: Genome, rec, st: int, a: int, b: int) -> str:
    gl = [(int(p.start), int(p.end)) for p in _location(G, [(st, a, b)]).parts]
    return next((n.qualifiers.get("note", [""])[0] for n in rec.features if n.type == "intron" and
                 (n.location.strand or 1) == st and _overlaps(n.location, gl) >= 0.5 * (b - a)), "")


def _same_intron(G: Genome, x: Tuple[int, int, int], y: Tuple[int, int, int]) -> bool:
    """Same strand and overlapping by >= 50 % of the shorter (strand coordinates, circular)."""
    if x[0] != y[0]:
        return False
    o = _ov((x[1], x[2]), (y[1], y[2]), G.L) if G.circular else max(0, min(x[2], y[2]) - max(x[1], y[1]))
    return o >= 0.5 * min(x[2] - x[1], y[2] - y[1])


def _inside(G: Genome, x: int, y: int, a: int, b: int) -> Optional[int]:
    """x..y (strand coordinates) strictly inside a..b: x shifted onto the circle turn of a..b, else None."""
    x0 = x if x >= a else x + G.L
    return x0 if a < x0 and x0 + (y - x) < b else None


CLOSE_IDENT = 0.9       # junction identity of a site that re-places a canonical de novo placement (a near relative)
IDENTICAL = 0.99        # junction identity of an identical junction: the reference's placement is kept as it is


class KeepDeNovo(Exception):
    """The de novo model of a gene is better left as it is."""


def merge_existing(G: Genome, rec, m: GeneModel, code: int) -> str:
    """Combine m with the model the de novo step (mito_intron.py) already wrote for the gene. The de novo exons are the
    base: at a de novo intron that is a library site the transferred exon ends replace the de novo ones; a transferred
    site inside a de novo exon splits it; de novo introns at no library site stay. Only when a transferred site lies
    outside the de novo model (pieces of a trans-spliced gene it missed) is m the base, its exons split at the de novo
    introns they contain. Returns a note on the merge ('' if there was no de novo model)."""
    olds = _existing_models(G, rec, m)
    if not olds:
        return ""
    old = olds[0]
    ex = old[1]
    mine = [(k, a, b) for k, (a, b) in enumerate(zip(m.exons, m.exons[1:]))]
    nums = {i.number for i in m.introns}

    def calls_of(k: int) -> List[IntronCall]:
        return [i for i in m.introns if i.number == k + 1] if nums else []
    base, introns, done = list(ex), [], []
    def gaps_of(exs):
        return [(i, G.norm(a[0], a[2], b[1]) or (a[0], a[2], b[1])) for i, (a, b) in enumerate(zip(exs, exs[1:]))
                if a[0] == b[0] and _signed_gap(G, a[2] - 1, b[1], 3) >= MIN_QUERY_INTRON]
    de_novo = gaps_of(ex)
    on_base, skipped = len(olds) == 1, []                    # several de novo pieces: m is the base
    footprint = _spans_of(G, ex)

    def in_footprint(st: int, pos: int) -> bool:
        return any(fs == st and _inside(G, pos, pos + 1, fa - 1, fb + 1) is not None for fs, fa, fb in footprint)
    for k, a, b in (mine if on_base else []):
        if a[0] != b[0]:
            on_base = False                                   # trans-spliced: the de novo step has no such junction
            break
        gap = G.norm(a[0], a[2], b[1]) or (a[0], a[2], b[1])
        # the de novo intron this site re-places: same strand, overlapping by >= 50 % of the shorter, and the
        # transferred intron swallows no de novo exon (else it would delete one, e.g. the 9-nt Plakinastrella exon)
        swallows = any(es == gap[0] and _inside(G, ea, eb, gap[1] - 1, gap[2] + 1) is not None for es, ea, eb in base)
        hit = None if swallows or a[0] != b[0] else \
            next((i for i, g in de_novo if _same_intron(G, gap, g)), None)
        if hit is not None:
            g = de_novo[[i for i, _ in de_novo].index(hit)][1]
            cm = _intron_note(G, rec, *g)                         # the de novo step's typing of this query intron
            mine_k = calls_of(k)
            close = any(c_.fixed or c_.family or c_.ident >= CLOSE_IDENT for c_ in mine_k)
            st_, d_, a_ = g[0], g[1] - 1, g[2]
            de_novo_ok = meets(canonical_kinds(G.sub(st_, d_, d_ + 1), G.sub(st_, d_ + 1, d_ + 7),
                                               G.sub(st_, a_ - 6, a_)),
                               allowed_groups(cm, next((c_.itype for c_ in mine_k), ""))) == "full"
            done.append(hit)
            if close or not de_novo_ok:                      # re-place the de novo boundaries
                base[hit] = G.norm(base[hit][0], base[hit][1], a[2]) or base[hit]
                base[hit + 1] = G.norm(base[hit + 1][0], b[1], base[hit + 1][2]) or base[hit + 1]
                for c_ in mine_k:
                    c_.cm = cm
                introns += mine_k
            else:                                            # a distant relative's site: the canonical de novo ends stay
                src = ", ".join(sorted({c_.ref.split("|")[0] for c_ in mine_k}))
                introns.append(IntronCall(g[0], g[1], g[2], 0, "", "de novo", 0.0,
                                          (cm or "intron predicted de novo") + f"; site also in {src}", cm=cm))
            continue
        slot = next(((n, x0) for n, (es, ea, eb) in enumerate(base) if a[0] == b[0] == es
                     for x0 in [_inside(G, a[2], b[1], ea, eb)] if x0 is not None), None)
        if slot is None:
            if in_footprint(a[0], a[2] - 1) or in_footprint(b[0], b[1]):
                skipped.append(f"{_gpos(G, a[0], a[2] - 1)}/{_gpos(G, b[0], b[1])}")   # crosses a de novo exon
                continue
            on_base = False                                   # a site outside the de novo model
            break
        n, x0 = slot
        es, ea, eb = base[n]
        base[n:n + 1] = [(es, ea, x0), (es, x0 + (b[1] - a[2]), eb)]
        de_novo = [(i + (1 if i >= n else 0), g) for i, g in de_novo]
        done = [i + (1 if i >= n else 0) for i in done]
        introns += calls_of(k)
    if on_base:
        kept = [(i, g) for i, g in de_novo if i not in done]
        for i, (st, x, y) in kept:
            note = _intron_note(G, rec, st, x, y)
            introns.append(IntronCall(st, x, y, 0, "", "de novo", 0.0, note or "intron predicted de novo", cm=note))
        exons, msg = base, (f"de novo model kept, {len(done)} intron(s) re-placed at library sites, "
                            f"{len(kept)} de novo intron(s) at no library site"
                            + (f", transferred site(s) {', '.join(skipped)} not used (cross a de novo exon)"
                               if skipped else ""))
    else:
        exons, introns = list(m.exons), list(m.introns)
        gaps = [G.norm(a[0], a[2], b[1]) or (a[0], a[2], b[1]) for a, b in zip(m.exons, m.exons[1:]) if a[0] == b[0]]
        kept, lost = [], []
        for i, g in [x for o in olds for x in gaps_of(o[1])]:
            if any(_same_intron(G, g, x) for x in gaps):
                continue
            for n, (es, ea, eb) in enumerate(exons):
                x0 = _inside(G, g[1], g[2], ea, eb) if es == g[0] else None
                if x0 is not None:
                    exons[n:n + 1] = [(es, ea, x0), (es, x0 + (g[2] - g[1]), eb)]
                    note = _intron_note(G, rec, *g)
                    introns.append(IntronCall(g[0], x0, x0 + (g[2] - g[1]), 0, "", "de novo", 0.0,
                                              note or "intron predicted de novo", cm=note))
                    kept.append(g)
                    break
            else:
                lost.append(g)
        added = [i for i in m.introns if i.ref != "de novo" and not i.half
                 and not any(_same_intron(G, (i.strand, i.start, i.end), g) for o in olds for _, g in gaps_of(o[1]))]
        if len(lost) >= max(1, len(added)) and \
                not all(i.fixed or i.ident >= CLOSE_IDENT for i in m.introns if i.ref != "de novo"):
            # a model from more distant relatives that would drop at least as many de novo introns as it adds: the
            # de novo model stays
            raise KeepDeNovo(f"transferred model would drop {len(lost)} de novo intron(s) and add {len(added)}; "
                             "de novo model kept")
        msg = f"transferred model (a site outside the de novo model), {len(kept)} de novo intron(s) kept inside it"
    if m.kind == "CDS":
        nt = "".join(G.sub(st, a, b) for st, a, b in exons)
        prot = codon_translate(nt, code).rstrip("*")
        if "*" in prot:
            return "not merged with the de novo model: internal stop codon"
        m.protein = prot

    def junction(i: IntronCall) -> float:          # transcript position: after the exon ending at its start
        for k, (es, a, b) in enumerate(exons):
            if es == i.strand and _near(G, b, i.start, 15):
                return k + 0.5
        for k, (es, a, b) in enumerate(exons):
            if es == i.strand and _near(G, a, i.end, 15):
                return k - 0.5
        return len(exons)
    for i in introns:
        i.number = int(junction(i) + 0.5)
    m.exons, m.introns = exons, sorted(introns, key=junction)
    if on_base:                                        # ends as the de novo model wrote them (5' / 3' partial)
        from Bio.SeqFeature import AfterPosition, BeforePosition
        first, last = old[0].location.parts[0], old[0].location.parts[-1]
        m.partial5 = isinstance(first.start if (first.strand or 1) == 1 else first.end, (BeforePosition, AfterPosition))
        m.partial3 = isinstance(last.end if (last.strand or 1) == 1 else last.start, (BeforePosition, AfterPosition))
    return msg


def apply_models(rec, models: List[GeneModel], code: int, origin_introns: bool = False) -> List[str]:
    """Write the gene models (confidence high / medium) into the record: gene + CDS / rRNA replaced (features of the
    same gene overlapping the new exons), intron features added. An intron across the origin of a circular record is
    left to the CDS join (table2asn rejects multi-interval introns) unless origin_introns (evaluation only: it is
    then written as a join() intron feature). Returns the genes written."""
    G = Genome(str(rec.seq), rec.annotations.get("topology", "circular") != "linear")
    written = []
    for m in models:
        if m.confidence == "low":
            continue
        try:
            name = _apply_model(rec, G, m, code, origin_introns)
        except (ValueError, KeepDeNovo) as exc:            # impossible location / de novo model kept: not written
            m.confidence, m.evidence = "low", f"{m.evidence}; not written: {exc}"
            continue
        written.append(name)
    if written:
        rec.features.sort(key=lambda f: (f.type != "source", int(f.location.start), f.type != "gene"))
    return written


def _apply_model(rec, G: Genome, m: GeneModel, code: int, origin_introns: bool) -> str:
    """Write one model: merged with a de novo model of the gene, introns placed by the gold-standard rules, new
    gene / CDS (rRNA) / intron features built first, then the replaced features removed (a ValueError leaves the
    record as it was)."""
    from Bio.SeqFeature import SeqFeature
    ftype = m.kind
    merged = merge_existing(G, rec, m, code)
    if merged:
        m.evidence += "; " + merged
    place_introns(G, m, code)
    if ftype == "CDS":
        nt = "".join(G.sub(st, a, b) for st, a, b in m.exons)
        m.protein = codon_translate(nt, code).rstrip("*")
    loc = _location(G, m.exons, m.partial5, m.partial3)
    gparts = [(int(p.start), int(p.end)) for p in loc.parts]
    spans = _spans_of(G, m.exons)
    span = [(int(p.start), int(p.end)) for p in _location(G, spans).parts]
    old = [f for f in rec.features if f.type in ("gene", ftype, "mRNA") and
           gene_name(f.qualifiers, ftype if f.type != "gene" else ("rRNA" if ftype == "rRNA" else "CDS")) == m.gene
           and _overlaps(f.location, gparts) > 0]
    old += [f for f in rec.features if f.type == "intron" and _overlaps(f.location, span) > 0 and
            (gene_name({"gene": f.qualifiers.get("gene", [""])}, ftype) or m.gene) == m.gene]
    name = next((f.qualifiers["gene"][0] for f in old if f.type == ftype and f.qualifiers.get("gene")), m.gene)
    product = next((f.qualifiers["product"][0] for f in old if f.type == ftype and f.qualifiers.get("product")),
                   PRODUCTS.get(m.gene, m.gene))
    gq: Dict[str, List[str]] = {"gene": [name]}
    note = f"exon-intron structure from known intron sites ({', '.join(m.refs)}; mito_intron_ref); {m.evidence}"
    q: Dict[str, List[str]] = {"gene": [name]}
    if ftype == "CDS":
        q.update({"codon_start": ["1"], "transl_table": [str(code)], "product": [product], "note": [note]})
        for k, v in ({} if m.partial3 else stop_qualifiers(rec, loc, code)).items():
            q[k] = q.get(k, []) + v
        q["translation"] = [m.protein]
    else:
        q.update({"product": [product], "note": [note]})
    if m.trans():
        gq["trans_splicing"] = [None]                      # value-less qualifier: written as /trans_splicing
        q["trans_splicing"] = [None]
    new = [SeqFeature(_location(G, spans, m.partial5, m.partial3), type="gene", qualifiers=gq),
           SeqFeature(loc, type=ftype, qualifiers=q)]
    for x in m.introns:
        iloc = _location(G, [(x.strand, x.start, x.end)])
        if len(iloc.parts) > 1 and not origin_introns:
            # across the origin: GenBank allows no multi-interval intron (table2asn MultiIntervalIntron); NCBI
            # records leave such an intron implicit in the CDS join (Trichoplax NC_008151.3 COX1a)
            q["note"].append(f"intron {x.number} ({_span_text(iloc)}) crosses the origin: given by the join only; "
                             + x.note)
            continue
        new.append(SeqFeature(iloc, type="intron",
                              qualifiers={"gene": [name], "number": [str(x.number)], "note": [x.note]}))
    rec.features = [f for f in rec.features if f not in old] + new
    return name


# ── command line ─────────────────────────────────────────────────────────────
REPORT_COLS = ["record", "gene", "ref", "ref_organism", "junction", "lineage_shared", "state", "splice", "donor",
               "acceptor", "ident", "ident_ex", "ident_in", "score", "why"]
_JOB: Dict[str, object] = {}


def _annotate_file(f: Path) -> List[dict]:
    """Annotate every record of one GenBank file; a failing record is reported and written back unchanged."""
    import shutil
    a, lib, own, tools, src_root = _JOB["a"], _JOB["lib"], _JOB["own"], _JOB["tools"], _JOB["root"]
    rows: List[dict] = []
    recs = list(SeqIO.parse(str(f), "genbank"))
    changed = False
    for rec in recs:
        try:
            mine = [g for g in own if _base(g.acc) == _base(rec.id)]
            # --lineage only for a record without its own taxonomy (FASTA input): an NCBI record keeps its lineage
            lineage = list(rec.annotations.get("taxonomy") or []) or [x for x in a.lineage.split(";") if x] or None
            models, r = annotate_record(rec, lib, a.genetic_code, lineage=lineage, tools=tools, own=mine,
                                        exclude=env_exclude(lib, rec))
            rows += r
            written = apply_models(rec, models, a.genetic_code)
            changed |= bool(written)
            for m in models:
                rows.append({"record": rec.id, "gene": m.gene, "state": "model",
                             "why": f"{m.confidence}: {m.evidence}" + ("" if m.confidence != "low" else " - not written"),
                             "ref": ",".join(m.refs)})
        except Exception as exc:                          # keep going: record is written back unchanged
            rows.append({"record": rec.id, "state": "error", "why": f"record left unchanged: {exc!r}"})
    dest = Path(a.out) / f.relative_to(src_root)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if changed:
        write_gb(recs, dest)
    elif dest.resolve() != f.resolve():
        shutil.copy2(f, dest)
    return rows


def env_exclude(lib: List[RefGene], rec) -> Optional[Exclude]:
    """Leave-out benchmarks (mitom_lo): MITOM_EXCLUDE_REF=<acc> leaves that record and its RefSeq / INSDC copy out of
    the library; MITOM_EXCLUDE_LEVEL=species|genus also every record of its species / genus (the organism of the
    excluded record in the library, else the query's). Unset: the whole library."""
    import os
    acc = os.environ.get("MITOM_EXCLUDE_REF", "").strip()
    if not acc:
        return None
    ex = Exclude({acc}).expand(lib)
    level = os.environ.get("MITOM_EXCLUDE_LEVEL", "").strip().lower()
    if level in ("species", "genus"):
        org = next((g.organism for g in lib if ex(g)), "") or rec.annotations.get("organism", "")
        if org:
            ex.species = species_of(org)
            ex.genus = genus_of(org) if level == "genus" else ""
    return ex


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="intron-site library from annotated GenBank records")
    b.add_argument("-i", "--input", nargs="+", required=True, help="GenBank files / folders (multi-record files allowed)")
    b.add_argument("-o", "--out", default=str(DEFAULT_LIB), help=f"site table (default {DEFAULT_LIB}); the FASTA of "
                   "spliced sequences is written next to it")
    b.add_argument("--jobs", type=int, default=1)
    b.add_argument("--bin", default=None, help="folder with cmsearch (types the library introns without annotated type)")
    b.add_argument("--cm", default=str(DEFAULT_CM), help="Rfam group I / II intron covariance models")
    n = sub.add_parser("annotate", help="add introns at known sites to GenBank records (e.g. the annotation core's output)")
    n.add_argument("-i", "--input", required=True, help="GenBank file or folder")
    n.add_argument("-o", "--out", required=True, help="output folder (same relative paths; unchanged records copied)")
    n.add_argument("--genetic_code", type=int, default=5)
    n.add_argument("--lib", default=str(DEFAULT_LIB), help="site library (build)")
    n.add_argument("--ref_gb", action="append", default=[],
                   help="annotated GenBank file(s) / folder(s): the sites of a record with the query's accession are its own "
                        "annotation (closest reference, quirks kept); other records join the library")
    n.add_argument("--own_only", action="store_true", help="only the records' own annotation (--ref_gb), no library")
    n.add_argument("--lineage", default="", help="NCBI lineage, ';'-separated, for records without their own taxonomy "
                   "(FASTA input; e.g. 'Eukaryota;Metazoa;Cnidaria'); a record's taxonomy is always used first")
    n.add_argument("--bin", default=None, help="folder with tblastn / blastn")
    n.add_argument("--jobs", type=int, default=1, help="GenBank files processed in parallel")
    a = p.parse_args(argv)
    if a.cmd == "build":
        genes = build_library(a.input, a.jobs)
        cms = find_tool("cmsearch", a.bin)
        typed = type_introns(genes, cms, Path(a.cm), cpu=max(1, a.jobs))
        save_library(genes, Path(a.out))
        print(f"[mito_intron_ref] {len(genes)} genes, {sum(len(g.junctions) for g in genes)} sites from "
              f"{len({g.acc for g in genes})} records ({typed} intron types by cmsearch) -> {a.out}", file=sys.stderr)
        return 0
    tools = find_tools(a.bin)
    if not tools.get("tblastn") or not tools.get("blastn"):
        print("[mito_intron_ref] tblastn / blastn not found (conda install -c bioconda blast; or --bin)", file=sys.stderr)
        return 2
    lib = [] if a.own_only else load_library(Path(a.lib))
    own: List[RefGene] = []
    if a.ref_gb:
        ids = set()
        src = Path(a.input)
        for f in iter_gb_files(src):
            ids |= {_base(r.id) for r in SeqIO.parse(str(f), "genbank")}
        for g in build_library(a.ref_gb):
            (own if _base(g.acc) in ids else lib).append(g)
    src = Path(a.input)
    files = iter_gb_files(src)
    _JOB.update(a=a, lib=lib, own=own, tools=tools, root=src.parent if src.is_file() else src)
    Path(a.out).mkdir(parents=True, exist_ok=True)
    rows: List[dict] = []
    if a.jobs > 1 and len(files) > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(a.jobs) as pool:
            for r in pool.imap(_annotate_file, files):
                rows += r
    else:
        for f in files:
            rows += _annotate_file(f)
    import csv
    with open(Path(a.out) / "intron_ref_report.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, REPORT_COLS, delimiter="\t", restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    models = [r for r in rows if r.get("state") == "model"]
    print(f"[mito_intron_ref] {len(files)} file(s), library {len(lib)} genes + {len(own)} own; "
          f"{sum(not r['why'].startswith('low') for r in models)} gene model(s) written, "
          f"{sum(r['why'].startswith('low') for r in models)} low-confidence not written -> "
          f"{Path(a.out) / 'intron_ref_report.tsv'}", file=sys.stderr)
    return 0


def find_tool(name: str, bin_dir: Optional[str]) -> Optional[str]:
    try:
        import mito_env
        return mito_env.find(name, bin_dir)
    except Exception:
        import shutil
        return str(Path(bin_dir) / name) if bin_dir and (Path(bin_dir) / name).exists() else shutil.which(name)


if __name__ == "__main__":
    sys.exit(main())
