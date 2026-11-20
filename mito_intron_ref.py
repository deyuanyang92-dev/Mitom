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
Benchmark (benchmark/bench_intron_ref.py, gold standard 2026-10-29): docs/METHODS.md, intron table.

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
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO

import mito_names  # noqa: E402
from mito_common import codon_translate, iter_gb_files, stop_qualifiers, write_gb, write_tsv  # noqa: E402
import mito_common  # noqa: E402
from mito_intron_sites import (  # noqa: E402, F401  (GROUP_II_5, GROUP_II_3 and some others re-exported)
    GROUP_II_3,
    GROUP_II_5,
    _circ,
    _ov,
    _rc,
    build_library,
    canonical_kinds,
    DEFAULT_CM,
    DEFAULT_LIB,
    fasta_path,
    gene_name,
    intron_kind,
    Junction,
    load_library,
    phylum_of,
    PROBE,
    ref_genes,
    RefGene,
    save_library,
    type_introns,
)
from mito_intron_junction import (  # noqa: E402
    _blast,
    _near,
    _signed_gap,
    call_junction,
    Genome,
    Hsp,
    JCall,
    MIN_IDENT_OTHER_PHYLUM,
    MIN_QUERY_INTRON,
    parse_hsps,
    parse_probes,
    ProbeHit,
    score_hsps,
)
from mito_intron_models import (  # noqa: E402, F401  (re-exported: gold_shift, vote)
    _base,
    _codons,
    _location,
    _overlaps,
    allowed_groups,
    assemble,
    CLOSE_IDENT,
    Exclude,
    GeneModel,
    genus_of,
    gold_shift,
    IntronCall,
    meets,
    place_introns,
    ref_depth,
    species_of,
    vote,
)
from mito_common import cli_guard, prepare_input, require_records  # noqa: E402


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


def _span_text(loc) -> str:
    return ",".join(f"{'complement(' if p.strand == -1 else ''}{int(p.start) + 1}..{int(p.end)}"
                    f"{')' if p.strand == -1 else ''}" for p in loc.parts)


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
            # the merged exons do not translate: the de novo introns outside the library were wrong, so the transferred
            # model alone is used - unless that does not translate either, then the de novo model is the better base
            # (2026-11-20 review: both were replaced by a model with an internal stop)
            alone = "".join(G.sub(st, a, b) for st, a, b in m.exons)
            if "*" in codon_translate(alone, code).rstrip("*"):
                raise KeepDeNovo("neither the merged nor the transferred model translates (internal stop codon); "
                                 "de novo model kept")
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
        # the de novo step's description of an intron across the origin (mito_intron.INTRON_NOTE: type, E-value,
        # splice score) stays with the gene: mitom_orient restores the intron feature from it after rotation
        q["note"] += [x for f in old if f.type == "CDS" for x in f.qualifiers.get("note", [])
                      if x.startswith("intron across the origin of the record")]
        for k, v in ({} if m.partial3 else stop_qualifiers(rec, loc, code)).items():
            q[k] = q.get(k, []) + v
        q["translation"] = [("M" + m.protein[1:]) if not m.partial5 and nt[:3] in _codons(code)[0] else m.protein]   # (Kimi 3.3)
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


class _ExcludeAll:
    """Excludes every library record (MITOM_EXCLUDE_LEVEL=all)."""

    def __call__(self, g) -> bool:
        return True


def env_exclude(lib: List[RefGene], rec) -> Optional[Exclude]:
    """Leave-out benchmarks (mitom_lo): MITOM_EXCLUDE_REF=<acc> leaves that record and its RefSeq / INSDC copy out of
    the library; MITOM_EXCLUDE_LEVEL=species|genus also every record of its species / genus (the organism of the
    excluded record in the library, else the query's); MITOM_EXCLUDE_LEVEL=all leaves the whole library out (equal-data
    benchmarks: de novo introns only). Unset: the whole library."""
    import os
    if os.environ.get("MITOM_EXCLUDE_LEVEL", "").strip().lower() == "all":
        return _ExcludeAll()
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


def build_parser() -> argparse.ArgumentParser:
    """The command line (main); options() gives the same options without one."""
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
    return p


def options(cmd: str, **given) -> argparse.Namespace:
    """Options for run(cmd): "build" (input: list of GenBank paths) or "annotate" (input, out)."""
    return mito_common.options(build_parser(), cmd, **given)


def run(a: argparse.Namespace) -> int:
    """a.cmd build: the intron-site library from annotated records; annotate: introns at known sites added to the
    GenBank file / folder a.input -> a.out (intron_ref_report.tsv). Returns the exit code."""
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
        return 1
    lib = [] if a.own_only else load_library(Path(a.lib))
    own: List[RefGene] = []
    if a.ref_gb:
        ids = set()
        src = Path(a.input)
        for f in iter_gb_files(src):
            ids |= {_base(r.id) for r in SeqIO.parse(str(f), "genbank")}
        for g in build_library(a.ref_gb):
            (own if _base(g.acc) in ids else lib).append(g)
    src = Path(prepare_input(a.input))
    files = iter_gb_files(src)
    require_records(len(files), "GenBank files", a.input)
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
    write_tsv(Path(a.out) / "intron_ref_report.tsv", REPORT_COLS, rows)
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


def main(argv=None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    sys.exit(cli_guard(main))
