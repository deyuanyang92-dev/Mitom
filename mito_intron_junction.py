"""Query genome (circular coordinates, BLAST HSP parsing) and junction calling of mito_intron_ref.py: the intron junctions of a reference gene placed on the query by probe alignment, scored and refined. Code moved verbatim from mito_intron_ref.py."""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio.Seq import Seq

from mito_common import run_tool
from mito_intron_sites import (
    PROBE,
    Junction,
    RefGene,
    _circ,
    _rc,
)


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
        p = run_tool(cmd, timeout=TIMEOUT)
        out = p.stdout
    except (subprocess.TimeoutExpired, OSError) as exc:
        print(f"[mito_intron_ref] {Path(exe).name} failed: {exc}", file=sys.stderr)
        return []
    if p.returncode != 0:                              # (2026-10-31 review: a failing tblastn / blastn was silent)
        print(f"[mito_intron_ref] WARNING: {Path(exe).name} failed (exit {p.returncode}): "
              f"{(p.stderr or '').strip()[-200:]}", file=sys.stderr)
        return []
    return [line.split("\t") for line in out.splitlines() if line.count("\t") == 9]


TIMEOUT = 900


SPLIT_GAP = 30                  # a subject insertion of > 30 nt (blastn; tblastn > 10 codons) inside one HSP is an intron: split there


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
