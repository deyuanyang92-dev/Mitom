"""Short, divergent exons that Exonerate does not align (issue #15).

Exonerate reports an alignment only when it scores >= MIN_REPORT_SCORE, so an exon of 28-81 aa of a divergent protein
(bryozoan nad4, cox2, cob, nad4L) is missed and the chain left has a single exon or a hole. The search here is
homology-constrained and tests significance:

* where: only the nt interval where the missing piece must lie (between the chain and the next chain, or up to
  `maxintron` before / after the chain), never the whole genome; no hit may overlap a coding sequence in `blocked`;
* what: the piece of every homolog of the gene that corresponds to the missing residues of the query protein, local
  alignment (BLOSUM62, gap open -11, ext -1) to the three translations of the interval;
* significance: one divergent homolog scores like noise, but the homologs agree on the place. A locus is accepted when
  >= MIN_VOTES homologs hit it (same frame, starts within VOTE_NT nt) and more than any locus obtained the same way
  from shuffled homolog pieces (empirical null, NULL_RUNS shuffles of every piece);
* then the model goes through the usual junction refinement, splicing gain and continuity gates of mito_intron.
"""
from __future__ import annotations

import random
from typing import List, Optional, Tuple

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

MIN_MISSING_AA = 12        # shorter holes are left to mito_intron_rules.add_mini_exons
MIN_VOTES = 4              # homologs that must agree on the locus
MIN_FRACTION = 0.15        # ... and this fraction of the homologs used
VOTE_NT = 90               # their hits start within this many nt (same frame)
NULL_RUNS = 6
MAX_HOMOLOGS = 30
MAX_SHORT_INTRON = 1000   # intron beside a new short exon (bryozoan: 250-632 nt); a larger window made false end exons
MIN_HIT_AA = 10


def _aligner() -> PairwiseAligner:
    return PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                           open_gap_score=-11, extend_gap_score=-1)


def homolog_pieces(prot: str, homs: List[str], q_lo: int, q_hi: int, al: PairwiseAligner) -> List[str]:
    """The part of each homolog that corresponds to query residues [q_lo, q_hi) (from the local alignment of the query
    to the homolog; proportional when the alignment does not reach the piece)."""
    out = [prot[q_lo:q_hi]]
    for h in homs:
        lo = hi = None
        try:
            aln = al.align(prot, h)[0]
            qs, hs = aln.aligned
            for (qa, qb), (ha, _hb) in zip(qs, hs):
                if lo is None and qb > q_lo:
                    lo = ha + max(0, q_lo - qa)
                if qa < q_hi:
                    hi = ha + min(qb - qa, q_hi - qa)
        except (IndexError, ValueError):
            pass
        if lo is None or hi is None or hi - lo < 4:
            lo, hi = int(len(h) * q_lo / len(prot)), int(len(h) * q_hi / len(prot)) + 1
        out.append(h[max(0, lo - 3):hi + 3])
    return [p for p in out if len(p) >= MIN_HIT_AA]


def _frames(use: str, nt_lo: int, nt_hi: int, gencode: int) -> List[Tuple[int, str]]:
    w = use[nt_lo:nt_hi]
    return [(nt_lo + f, str(Seq(w[f:len(w) - (len(w) - f) % 3]).translate(table=str(gencode)))) for f in range(3)]


def _best_hit(frames, piece: str, al: PairwiseAligner):
    best = None
    for off, t in frames:
        if not t:
            continue
        try:
            a = al.align(t, piece)[0]
        except (IndexError, ValueError):
            continue                  # no positive-scoring alignment in this frame (e.g. a window of N)
        sc = float(al.score(t, piece))
        if best is None or sc > best[0]:
            ts, ps = a.aligned
            best = (sc, off, int(ts[0][0]), int(ts[-1][1]), int(ps[0][0]), int(ps[-1][1]))
    return best


def _vote(hits) -> Tuple[int, List[tuple]]:
    """Largest set of hits agreeing on the locus (same frame, starts within VOTE_NT nt)."""
    best: List[tuple] = []
    for h in hits:
        if h is None or h[3] - h[2] < MIN_HIT_AA:
            continue
        s = h[1] + 3 * h[2]
        grp = [x for x in hits if x is not None and x[1] % 3 == h[1] % 3 and abs(x[1] + 3 * x[2] - s) <= VOTE_NT
               and x[3] - x[2] >= MIN_HIT_AA]
        if len(grp) > len(best):
            best = grp
    return len(best), best


def locate(use: str, prot: str, homs: List[str], q_lo: int, q_hi: int, nt_lo: int, nt_hi: int, gencode: int,
           blocked: List[Tuple[int, int]], min_intron: int, al: Optional[PairwiseAligner] = None, seed: int = 1
           ) -> List[Tuple[int, int]]:
    """Exons (0-based half-open nt intervals of `use`, whole codons) for query residues [q_lo, q_hi) inside the nt
    interval [nt_lo, nt_hi); recursive: after the best locus the residues and the interval on either side are searched
    again. Empty when no locus is significant."""
    if q_hi - q_lo < MIN_MISSING_AA or nt_hi - nt_lo < 3 * MIN_HIT_AA:
        return []
    al = al or _aligner()
    pieces = homolog_pieces(prot, homs[:MAX_HOMOLOGS], q_lo, q_hi, al)
    if len(pieces) < MIN_VOTES:
        return []
    frames = _frames(use, nt_lo, nt_hi, gencode)

    def free(h) -> bool:
        s, e = h[1] + 3 * h[2], h[1] + 3 * h[3]
        return not any(s < y and x < e for x, y in blocked)
    hits = [_best_hit(frames, p, al) for p in pieces]
    hits = [h if h and free(h) else None for h in hits]
    votes, grp = _vote(hits)
    if votes < max(MIN_VOTES, MIN_FRACTION * len(pieces)):
        return []
    null = 0
    rng = random.Random(seed)
    for _ in range(NULL_RUNS):
        sh = []
        for p in pieces:
            l = list(p)
            rng.shuffle(l)
            sh.append("".join(l))
        nh = [_best_hit(frames, p, al) for p in sh]
        null = max(null, _vote([h if h and free(h) else None for h in nh])[0])
    if votes <= null:
        return []
    starts = sorted(h[1] + 3 * h[2] for h in grp)
    ends = sorted(h[1] + 3 * h[3] for h in grp)
    f = grp[0][1] % 3
    s = starts[len(starts) // 2]
    s -= (s - f) % 3
    e = ends[len(ends) // 2]
    e -= (e - f) % 3
    if e - s < 3 * MIN_HIT_AA:
        return []
    # the residues of the query this exon holds: proportional to where the median hit lies inside the pieces
    frac_lo = sorted(h[4] / max(1, len(pieces[0])) for h in grp)[len(grp) // 2]
    frac_hi = sorted(min(1.0, h[5] / max(1, len(pieces[0]))) for h in grp)[len(grp) // 2]
    mid_lo = q_lo + int(frac_lo * (q_hi - q_lo))
    mid_hi = q_lo + int(frac_hi * (q_hi - q_lo))
    left = locate(use, prot, homs, q_lo, mid_lo, nt_lo, s - min_intron, gencode, blocked, min_intron, al, seed + 1)
    right = locate(use, prot, homs, mid_hi, q_hi, e + min_intron, nt_hi, gencode, blocked, min_intron, al, seed + 2)
    return left + [(s, e)] + right


FIT_WINDOW = 60            # nt each end of a new exon may move when its frame and boundaries are fitted


def fit_new_exons(use: str, exons: List[Tuple[int, int]], new: List[Tuple[int, int]], gencode: int, to_ref,
                  min_intron: int, window: int = FIT_WINDOW) -> List[Tuple[int, int]]:
    """Boundaries and frame of the exons found by `locate` (1-based inclusive transcript coordinates; `exons` = the whole
    chain, `new` the members found here). The vote gives only the place, its ends are whole codons of the hit and
    cannot know the phase at the junctions. Each new exon's ends are searched (+-window nt, intron >= min_intron to its
    neighbours; the outer end of a first / last exon stays): the placement kept is the one whose spliced protein has no
    internal stop codon and aligns best to the homologs (median of to_ref), ties to the smaller shift."""
    from Bio.Seq import Seq as _Seq
    ex = sorted(exons)
    for n in sorted(new):
        i = ex.index(n)
        lo_s = ex[i - 1][1] + min_intron + 1 if i > 0 else n[0]
        hi_e = ex[i + 1][0] - min_intron - 1 if i + 1 < len(ex) else n[1]
        best, best_key = None, None
        for a in range(-window if i > 0 else 0, (window if i > 0 else 0) + 1):
            st = n[0] + a
            if st < lo_s or st < 1:
                continue
            for b in range(-window if i + 1 < len(ex) else 0, (window if i + 1 < len(ex) else 0) + 1):
                en = n[1] + b
                if en > hi_e or en - st < 14 or en > len(use):
                    continue
                cand = ex[:i] + [(st, en)] + ex[i + 1:]
                nt = "".join(use[x - 1:y] for x, y in cand)
                prot = str(_Seq(nt[:len(nt) - len(nt) % 3]).translate(table=str(gencode)))
                if "*" in prot.rstrip("*"):
                    continue
                sc = to_ref(prot)
                sc = sorted(sc if isinstance(sc, list) else [sc])
                key = (sc[(len(sc) - 1) // 2], -(abs(a) + abs(b)))
                if best_key is None or key > best_key:
                    best, best_key = cand, key
        if best is not None:
            ex = best
    return ex


def complete_chain(use: str, tx_exons: List[Tuple[int, int]], reports: List[Tuple[int, int, int, int]], prot: str,
                   homs: List[str], gencode: int, blocked: List[Tuple[int, int]], maxintron: int, min_intron: int
                   ) -> Tuple[List[Tuple[int, int]], str]:
    maxintron = min(maxintron, MAX_SHORT_INTRON)
    """Exons of a chain of Exonerate alignments plus the short exons found for its holes. reports: (qab, qae, nt_start,
    nt_end) per alignment in 1-based transcript coordinates of `use`, sorted by qab. Returns (exons 1-based inclusive,
    note), the exons unchanged and an empty note when nothing significant is found."""
    reps = sorted(reports)
    holes = []
    if reps[0][0] >= MIN_MISSING_AA:
        holes.append((0, reps[0][0], max(0, reps[0][2] - 1 - maxintron), reps[0][2] - 1 - min_intron))
    for x, y in zip(reps, reps[1:]):
        if y[0] - x[1] >= MIN_MISSING_AA:
            holes.append((x[1], y[0], x[3] + min_intron, y[2] - 1 - min_intron))
    if len(prot) - reps[-1][1] >= MIN_MISSING_AA:
        holes.append((reps[-1][1], len(prot), reps[-1][3] + min_intron, min(len(use), reps[-1][3] + maxintron)))
    al = _aligner()
    found: List[Tuple[int, int]] = []
    notes = []
    for q_lo, q_hi, a, b in holes:
        got = locate(use, prot, homs, q_lo, q_hi, a, b, gencode, blocked, min_intron, al)
        found += got
        if got:
            notes.append(f"residues {q_lo + 1}-{q_hi}: {len(got)} short exon(s) found by homolog vote")
    if not found:
        return tx_exons, ""
    merged = sorted(list(tx_exons) + [(s + 1, e) for s, e in found])
    return merged, "; ".join(notes)
