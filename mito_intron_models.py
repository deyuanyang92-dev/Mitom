"""Gene models of mito_intron_ref.py: intron calls voted across references, assembled into exon models with start / end, and the gold-standard shift of site boundaries. Code moved verbatim from mito_intron_ref.py."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from mito_common import codon_translate
from mito_intron_junction import (  # noqa: E402
    Genome,
    JCall,
    _near,
    _signed_gap,
)
from mito_intron_sites import (
    RefGene,
    canonical_kinds,
    gene_name,
    intron_kind,
)


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


CM_SURE = 1e-5          # a covariance-model hit this good fixes the intron group


SYN_SHIFT = 60          # shifts outside the sliding window tried for a synonymous canonical placement


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


def _location(G: Genome, pieces: List[Tuple[int, int, int]], partial5: bool = False, partial3: bool = False):
    """(strand, start, end) strand-coordinate pieces in transcript order -> Biopython location (origin split)."""
    from Bio.SeqFeature import (
        AfterPosition,
        BeforePosition,
        CompoundLocation,
        ExactPosition,
        FeatureLocation,
    )
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


def _overlaps(loc, parts: List[Tuple[int, int]]) -> int:
    return sum(max(0, min(int(p.end), b) - max(int(p.start), a)) for p in loc.parts for a, b in parts)


CLOSE_IDENT = 0.9       # junction identity of a site that re-places a canonical de novo placement (a near relative)


IDENTICAL = 0.99        # junction identity of an identical junction: the reference's placement is kept as it is
