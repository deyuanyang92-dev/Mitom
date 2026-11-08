"""mito_intron_rules.py - the gene model of mito_intron.py and MFannot's boundary rules on it (port of mfannot
Adjust_all_intronic_junctions): splice sites moved within ADJUST_RANGE, group II sites scored on the homologs, mini
exons, the CDS ends extended to the start / stop codon (adjust_introns).

Part of the Python port of MFannot's intron workflow in mito_intron.py (Beck & Lang,
https://github.com/BFL-lab/Mfannot, GPL-2.0; data https://github.com/BFL-lab/MFannot_data, GPL-3.0; Lang et al.
2023, Front. Plant Sci., doi:10.3389/fpls.2023.1222186): this file is therefore distributed under GPL-3.0, as
mito_intron.py (whose docstring describes the ported steps and the MFannot line numbers)."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, List, Optional, Tuple
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data import CodonTable
from Bio.Seq import Seq


TYPE_I_MIN, TYPE_II_MIN = 19, 40     # mfannot:5735 / Adjust_TypeII
REFINE_5, REFINE_3, REFINE_MIN_GAIN = 6, 24, 10   # refine_junction: windows (nt) and median homolog gain needed
# MFannot defaults (PirModels/Option.pir, mfannot:284)
MIN_INTRON, MAX_INTRON, MIN_REPORT_SCORE = 142, 5000, 200
COPY_KMER, COPY_MIN_SHARED = 20, 60   # copied_exons: k-mer length, shared k-mers that make an intron hold an exon copy
                                     # (Hexamermis atp6 copies 340; one 33-46-nt repeat in hexacoral nad5 introns 14-27)


ADJUST_RANGE = range(-6, 7)          # mfannot:5640 %adjustement, same k at both ends


II_AA_TOLERANCE = 5                  # homolog score (median, BLOSUM62) a group II placement may lose when none is free


WIDE_II = 15                         # same-shift window for a group II intron with no canonical placement within 6 nt


MINI_EXON_MAX_AA = 10                # longest homolog gap at an exon junction searched for a mini-exon (codons)


STOP_OVERLAP_MAX = 30                # nt a 3' end may run into the next feature to reach a complete stop (mito_curate P3)


@dataclass
class Model:
    query: str
    gene: str
    strand: int                       # +1 / -1
    exons: List[Tuple[int, int]]      # 1-based inclusive, ascending genome coords
    score: int
    qcov: float
    qmiss_c: int = 0                  # homolog residues left unaligned at the C-terminus
    qmiss_n: int = 0                  # ... and at the N-terminus
    introns: List[dict] = field(default_factory=list)
    core_only: bool = False         # candidate only because the core's CDS covers part of its reference
    open5: bool = False             # a piece of a trans-spliced gene: the protein continues in another piece
    open3: bool = False             #   before (open5) / after (open3) it, so that end is left open (partial)
    piece: str = ""

    @property
    def span(self) -> Tuple[int, int]:
        return self.exons[0][0], self.exons[-1][1]


# ── step 3: boundary adjustment (port of mfannot Adjust_all_intronic_junctions) ──
def _ch(s: str, pos1: int) -> str:
    return s[pos1 - 1] if 1 <= pos1 <= len(s) else ""


def score_type_i(s: str, start: int, end: int, k: int) -> float:
    """ScoreNtLTypeI + ScoreNtRTypeI: exon ends ...GGT | intron ends with omega-G."""
    sc: float = 5 * (_ch(s, start + k - 1) == "T") + 3 * (_ch(s, start + k - 2) == "G") + 3 * (_ch(s, start + k - 3) == "G")
    sc += 10 * (_ch(s, end + k) == "G")
    sc += sum(0.5 for i in range(1, 9) if _ch(s, end + k - i) in ("A", "T"))
    return sc


def score_type_ii(s: str, start: int, end: int, k: int) -> float:
    """ScoreNtLTypeII + ScoreNtRTypeII - |k|: intron GTGYG ... AY."""
    sc: float = 10 * (_ch(s, start + k) == "G") + 7.5 * (_ch(s, start + k + 1) == "T") + 9 * (_ch(s, start + k + 2) == "G")
    sc += 10 * (_ch(s, start + k + 3) in ("C", "T")) + 10 * (_ch(s, start + k + 4) == "G")
    sc += 9 * (_ch(s, end + k) in ("C", "T")) + 9 * (_ch(s, end + k - 1) == "A")
    return sc - abs(k)


def _translate(s: str, exons: List[Tuple[int, int]], gencode: int) -> str:
    nt = "".join(s[a - 1:b] for a, b in exons)
    return str(Seq(nt[: len(nt) - len(nt) % 3]).translate(table=gencode))


def homolog_scorer(*refs: str):
    """Local alignment scores of a protein against homologs of the gene (the one that produced the model first)."""
    al = PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                         open_gap_score=-11, extend_gap_score=-1)
    return lambda prot: [al.score(prot, r) for r in refs]


def _aligner(open_gap: float) -> PairwiseAligner:
    return PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                           open_gap_score=open_gap, extend_gap_score=-1)


def homolog_gap(prot: str, refs: List[str], junction_aa: int, window: int = 4) -> Tuple[float, float]:
    """Residues the homologs have but the spliced protein lacks at an exon junction
    (AnnotateMiniExonsByHMM idea: a missing mini-exon shows as a gap in our protein).
    A single pairwise alignment rarely opens a gap for a few residues, so a mild gap
    penalty is used and the signal is taken over many homologs of the same gene.
    Returns (median gap, fraction of homologs with a gap >= 1)."""
    al = _aligner(-6)
    gaps = []
    for ref in refs:
        aln = al.align(prot, ref)[0]
        blocks = list(zip(aln.aligned[0], aln.aligned[1]))
        g = max([int((rb[0] - ra[1]) - (pb[0] - pa[1])) for (pa, ra), (pb, rb) in zip(blocks, blocks[1:])
                 if abs(pa[1] - junction_aa) <= window] or [0])
        gaps.append(g)
    if not gaps:
        return 0.0, 0.0
    gaps.sort()
    return float(gaps[len(gaps) // 2]), sum(g >= 1 for g in gaps) / len(gaps)


def find_mini_exon(s: str, exons: List[Tuple[int, int]], i: int, gencode: int, refs: List[str],
                   k_max: int) -> Optional[Tuple[int, int, float]]:
    """Search intron i (transcript coords) for a short in-frame exon (1..k_max codons) whose
    insertion improves the alignment to every homolog in refs without adding stop codons;
    both flanking introns must stay >= MIN_INTRON. Returns (start, end, mean gain)."""
    al = _aligner(-11)
    start, end = exons[i][1] + 1, exons[i + 1][0] - 1
    base_p = _translate(s, exons, gencode)
    base = [al.score(base_p, r) for r in refs]
    stops0 = base_p.rstrip("*").count("*")
    # stage 1 (fast): one homolog, local alignment of a window around the junction only
    j = sum(b - a + 1 for a, b in exons[:i + 1]) // 3
    w0 = max(0, j - 50)
    base_w = al.score(base_p[w0:j + 50], refs[0])
    shortlist: List[Tuple[float, int, int]] = []
    for ln in range(3, 3 * k_max + 1, 3):
        for x in range(start + MIN_INTRON, end - MIN_INTRON - ln + 2):
            new = exons[:i + 1] + [(x, x + ln - 1)] + exons[i + 1:]
            prot = _translate(s, new, gencode)
            if prot.rstrip("*").count("*") > stops0:
                continue
            g1 = al.score(prot[w0:j + 50 + ln // 3], refs[0]) - base_w
            if g1 > 0:
                shortlist.append((g1, x, x + ln - 1))
    shortlist = sorted(shortlist, reverse=True)[:5]
    # stage 2 (exact): full proteins against every homolog
    best: Optional[Tuple[float, int, int]] = None
    for _, x, y in shortlist:
        prot = _translate(s, exons[:i + 1] + [(x, y)] + exons[i + 1:], gencode)
        gains = [al.score(prot, r) - b for r, b in zip(refs, base)]
        if min(gains) <= 0:
            continue                          # must help every homolog
        g = sum(gains) / len(gains)
        if best is None or g > best[0]:
            best = (g, x, y)
    return (best[1], best[2], best[0]) if best and best[0] >= 2 * (best[2] - best[1] + 1) / 3 else None


def add_mini_exons(s: str, exons: List[Tuple[int, int]], kinds: List[Any], notes: List[str],
                   scores: List[float], gencode: int, refs: List[str]):
    """Insert homolog-supported mini-exons; returns updated exons / kinds / notes / scores."""
    i = 0
    while i < len(exons) - 1 and refs:
        junction = sum(b - a + 1 for a, b in exons[:i + 1]) // 3
        med, frac = homolog_gap(_translate(s, exons, gencode), refs, junction)
        # a gap longer than MINI_EXON_MAX_AA is a missing exon, not a mini-exon; introns longer than MAX_INTRON (that
        # host other genes) are not searched: both made the search run for hours (NC_015309.1, 2026-10-29)
        if 1 <= med <= MINI_EXON_MAX_AA and frac >= 0.5 and exons[i + 1][0] - exons[i][1] - 1 <= MAX_INTRON:
            mx = find_mini_exon(s, exons, i, gencode, refs[:3], int(med) + 4)
            if mx:
                exons = exons[:i + 1] + [(mx[0], mx[1])] + exons[i + 1:]
                tag = (f"mini-exon {mx[1] - mx[0] + 1} nt inserted (gap of ~{med:g} aa in {frac:.0%} of "
                       f"{len(refs)} homologs; MFannot AnnotateMiniExons idea)")
                kinds = kinds[:i] + [kinds[i], kinds[i]] + kinds[i + 1:]
                notes = notes[:i] + [f"{notes[i]}; {tag}"] * 2 + notes[i + 1:]
                scores = scores[:i] + [scores[i], scores[i]] + scores[i + 1:]
                i += 2
                continue
        i += 1
    return exons, kinds, notes, scores


def aa_score(s: str, exons: List[Tuple[int, int]], new: List[Tuple[int, int]], gencode: int,
             to_ref, tolerance: float = 0) -> Optional[float]:
    """CheckForCorrespondance, judged against the homolog protein.

    MFannot compares the shifted protein with the starting model and rejects any
    change scoring < 0 in PAM70. The starting model here comes from Exonerate's
    built-in splice model and can be wrong at the junction (Metridium senile COX1:
    model HMFTSG, NCBI and homolog HMFTVG), so the shifted protein is scored
    against the homolog instead: rejected when it aligns worse than the starting
    model, otherwise the score gain is returned. With several homologs the median gain decides (one distant or
    flawed homolog no longer vetoes a shift: Decemunciger nad1 against Homo sapiens nad1 refused the GUGCG
    placement that the annelid homolog prefers).
    """
    old_p, new_p = _translate(s, exons, gencode), _translate(s, new, gencode)
    if old_p == new_p:
        return 9999.0
    olds, news = to_ref(old_p), to_ref(new_p)
    if not isinstance(olds, list):
        olds, news = [olds], [news]
    gains = sorted(n - o for n, o in zip(news, olds))
    gain = gains[(len(gains) - 1) // 2]
    return None if gain < -tolerance else gain


def refine_junction(s: str, exons: List[Tuple[int, int]], i: int, gencode: int, to_ref,
                    w5: int = REFINE_5, w3: int = REFINE_3) -> Tuple[List[Tuple[int, int]], float]:
    """Intron i moved with its two ends shifted independently (5' by -w5..w5, 3' by -w3..w3 nt, frame kept), to the
    placement whose spliced protein aligns best to the homologs (median gain over the current one): exonerate's
    junction can be a few codons off at one end only (hexacoral nad5 NC_015642.1, NC_065774.1: 3' end 11-12 nt
    upstream), which the same-shift adjustment cannot correct. Returns (exons, gain); unchanged when no placement
    gains REFINE_MIN_GAIN or more."""
    best, best_gain = exons, 0.0
    base = _translate(s, exons, gencode)
    stops0 = base.rstrip("*").count("*")
    j = sum(b - a + 1 for a, b in exons[:i + 1]) // 3          # junction in the protein: score a window around it
    w0 = max(0, j - 60)
    def scores(prot: str) -> list:
        x = to_ref(prot)
        return x if isinstance(x, list) else [x]
    base_scores = scores(base[w0:j + 60])
    for a in range(-w5, w5 + 1):
        for b in range(-w3, w3 + 1):
            if (a - b) % 3 or (a == 0 and b == 0):
                continue
            new = exons[:i] + [(exons[i][0], exons[i][1] + a), (exons[i + 1][0] + b, exons[i + 1][1])] + exons[i + 2:]
            if not valid_exons(new) or new[i][1] - new[i][0] < 5 or new[i + 1][1] - new[i + 1][0] < 5:
                continue
            prot = _translate(s, new, gencode)
            if prot == base or prot.rstrip("*").count("*") > stops0:
                continue
            gains = sorted(x - y for x, y in zip(scores(prot[w0:j + 60 + (len(prot) - len(base))]), base_scores))
            gain = gains[(len(gains) - 1) // 2]
            if gain > best_gain:
                best, best_gain = new, gain
    return (best, best_gain) if best_gain >= REFINE_MIN_GAIN else (exons, 0.0)


def adjust_introns(s: str, exons: List[Tuple[int, int]], kinds: List[str], gencode: int,
                   to_ref, only: Optional[set] = None, w5: int = REFINE_5
                   ) -> Tuple[List[Tuple[int, int]], List[str], List[float]]:
    """s and exons in transcript orientation. Returns adjusted exons and a note per intron.
    An untyped intron is placed by the group II rule when its best allowed shift reaches TYPE_II_MIN (degenerate
    annelid group II introns without a domain V/VI hit: Chaetopterus cox1 GUGCG...AC, NC_028710.1), else scored as
    group I as before. only: the introns to adjust (the others get an empty note); w5: refine_junction's 5' window."""
    exons = list(exons)
    notes: List[str] = []
    scores: List[float] = []
    for i, kind in enumerate(kinds):
        if only is not None and i not in only:
            notes.append("")
            scores.append(0.0)
            continue
        exons, _ = refine_junction(s, exons, i, gencode, to_ref, w5=w5)
        start, end = exons[i][1] + 1, exons[i + 1][0] - 1
        prev_len, next_len = exons[i][1] - exons[i][0] + 1, exons[i + 1][1] - exons[i + 1][0] + 1

        def shifts(mode: str, tolerance: float = 0, window=ADJUST_RANGE) -> list:
            scored = []
            for k in window:
                if (prev_len + k <= 5 or next_len - k <= 5) and not (mode != "II" and k == 0):
                    continue
                nt = score_type_ii(s, start, end, k) if mode == "II" else score_type_i(s, start, end, k)
                new = exons[:i] + [(exons[i][0], exons[i][1] + k), (exons[i + 1][0] + k, exons[i + 1][1])] + exons[i + 2:]
                aa = aa_score(s, exons, new, gencode, to_ref, tolerance)
                if aa is not None:
                    scored.append((nt, aa, abs(k), k, new))
            return sorted(scored, key=lambda x: (-(x[0] - x[2]), -x[1], x[2], x[3]))
        scored = shifts(kind or "II")
        like_ii = not kind and bool(scored) and scored[0][0] >= TYPE_II_MIN
        if not kind and not like_ii:
            scored = shifts(kind)
        if kind == "II" and (not scored or scored[0][0] < TYPE_II_MIN):
            # no group II placement without a homolog loss: a GUGYG / AY placement is taken at a small loss
            # (Decemunciger nad1: the GUGCG start 6 nt upstream costs 4 BLOSUM62 points against the closest homolog),
            # then up to WIDE_II nt away (Decemunciger nad4: GUGYG...AY 11 nt downstream of exonerate's junction)
            near = [x for x in shifts("II", II_AA_TOLERANCE) if x[0] >= TYPE_II_MIN] or \
                [x for x in shifts("II", II_AA_TOLERANCE, range(-WIDE_II, WIDE_II + 1)) if x[0] >= TYPE_II_MIN]
            scored = near + scored
        if not scored:
            notes.append("splice boundaries uncertain")
            scores.append(0.0)
            continue
        nt, _, _, k, new = scored[0]
        scores.append(nt)
        threshold = TYPE_II_MIN if kind == "II" or like_ii else TYPE_I_MIN
        how = "group II splice sites (GUGYG...AY), " if like_ii else ""
        if nt >= threshold or kind == "I":      # mfannot: typed group I adjusted even below 19
            exons = new
            notes.append(how + (f"boundaries adjusted by {k:+d} nt (splice score {nt:g})" if k else f"splice score {nt:g}"))
        else:
            notes.append(f"splice boundaries uncertain (best score {nt:g})")
    return exons, notes, scores


def extend_to_stop(s: str, exons: List[Tuple[int, int]], gencode: int, limit: int, qmiss_c: int,
                   max_codons: int = 100, rule: str = "first") -> Tuple[List[Tuple[int, int]], bool]:
    """Set the 3' end of the CDS.

    Alignment reached the homolog's C-terminus (qmiss_c <= 1): the protein is complete,
    so end on the next codon if it is a stop, else on an incomplete stop (T / TA,
    completed by polyadenylation) as NCBI mitochondrial records do.
    Alignment stopped short: extend in frame to the first complete stop within max_codons
    (overlapping the next feature is allowed, as mitochondrial genes often overlap);
    without one, end on an incomplete stop just before the next feature (limit). A first complete
    stop more than STOP_OVERLAP_MAX nt inside the next feature gives way to an incomplete stop
    right before it (mito_curate P3's rule; Glycera KT989324.1 cox1 ran 717 nt into cox2 and the
    model was refused as overlapping COX2).
    """
    stops = CodonTable.unambiguous_dna_by_id[gencode].stop_codons
    a, b = exons[-1]

    def abbreviated() -> Optional[int]:
        """End of an in-frame T / TA right before limit, None if a complete stop comes first or there is none."""
        for c0 in range(b, min(limit, len(s)), 3):
            partial = s[c0:limit]
            if len(partial) < 3:
                return limit if partial in ("T", "TA") else None
            if s[c0:c0 + 3] in stops:
                return None
        return None
    if s[b - 3:b] in stops:
        return (exons), False
    if rule == "homolog" and qmiss_c <= 1:
        if s[b:b + 3] in stops and b + 3 <= limit:
            return (exons[:-1] + [(a, b + 3)]), False
        if s[b:b + 2] == "TA" and b + 2 <= limit:
            return (exons[:-1] + [(a, b + 2)]), False
        return (exons[:-1] + [(a, b + 1)] if s[b:b + 1] == "T" and b + 1 <= limit else exons), False
    # mitochondrial genes often overlap their neighbours, so a complete stop codon may lie
    # inside the next feature; only when none is found end on an incomplete stop before it
    # "first": up to the ORF end, one step past the last complete codon so that a frame running exactly to the end
    # of the sequence is 3'-partial too (LC670763.1 fragment: CDS ended on its last codon without '>', NoStop)
    n_codons = (len(s) - b) // 3 + 1 if rule == "first" else max(max_codons, qmiss_c + 5)
    for i in range(n_codons):
        c0 = b + 3 * i
        if c0 + 3 > len(s):
            if rule == "first":            # open frame runs off the 3' end: 3'-partial CDS
                return (exons[:-1] + [(a, c0)]), True
            break
        if s[c0:c0 + 3] in stops:
            t = abbreviated() if c0 + 3 > limit + STOP_OVERLAP_MAX else None
            return (exons[:-1] + [(a, t or c0 + 3)]), False
    t = abbreviated()
    return (exons[:-1] + [(a, t)] if t else exons), False


def core_ends(rec, m: Model, exons_tx: List[Tuple[int, int]], L: int, window: int = 90
               ) -> Tuple[Optional[int], Optional[int]]:
    """5' start / 3' end (transcript coords) of existing CDS features (the core's model) on the
    same strand that overlap the model and are in frame with its first / last exon."""
    a, b = exons_tx[0][0], exons_tx[-1][1]
    g_span = to_genome([(a, b)], m.strand, L)[0]
    start = end = None
    for ft in rec.features:
        if ft.type != "CDS" or ft.location.strand != m.strand:
            continue
        fs, fe = int(ft.location.start) + 1, int(ft.location.end)
        if not overlap((fs, fe), g_span):
            continue
        s_tx, e_tx = (fs, fe) if m.strand == 1 else (L + 1 - fe, L + 1 - fs)
        if abs(s_tx - a) <= window and (a - s_tx) % 3 == 0 and s_tx <= exons_tx[0][1]:
            start = s_tx if start is None or abs(s_tx - a) < abs(start - a) else start
        if abs(e_tx - b) <= window and e_tx >= exons_tx[-1][0]:
            end = e_tx if end is None or abs(e_tx - b) < abs(end - b) else end
    return start, end


def ends_with_stop(s: str, exons: List[Tuple[int, int]], gencode: int) -> bool:
    """The spliced CDS ends on a complete stop codon or a T / TA tail (abbreviated stop) in its own frame."""
    nt = "".join(s[a - 1:b] for a, b in exons)
    tail = len(nt) % 3
    return nt[-tail:] in ("T", "TA") if tail else nt[-3:] in CodonTable.unambiguous_dna_by_id[gencode].stop_codons


def copied_exons(s: str, exons: List[Tuple[int, int]], k: int = COPY_KMER) -> str:
    """Reason not to model a gene whose intron holds a copy of its own exons ('' if none): exonerate splices the
    N-terminal part of one copy of a duplicated gene to the C-terminal part of another (Hexamermis NC_008828.1, three
    identical atp6: CDS 6888..11167 with a 3.7-kb 'intron'). s and exons in transcript orientation; a copy = at
    least COPY_MIN_SHARED k-mers of the exons inside an intron (self-splicing introns pair only 4-7 nt with exons;
    a single 33-46-nt repeat of exon sequence occurs in true nad5 introns, NC_047221.1 / NC_082103.1)."""
    kmers = {s[a - 1 + i:a - 1 + i + k] for a, b in exons for i in range(b - a + 2 - k)}
    for i in range(len(exons) - 1):
        intr = s[exons[i][1]:exons[i + 1][0] - 1]
        shared = sum(intr[j:j + k] in kmers for j in range(len(intr) - k + 1))
        if shared >= COPY_MIN_SHARED:
            return (f"intron {i + 1} ({exons[i][1] + 1}..{exons[i + 1][0] - 1} on the transcript) contains a copy of the "
                    f"exons ({shared} shared {k}-mers): exons of duplicated gene copies, not one gene - not modelled")
    return ""


def valid_exons(exons: List[Tuple[int, int]]) -> bool:
    """Exons non-empty, ordered and separated by at least MIN_INTRON nt."""
    return all(a <= b for a, b in exons) and all(exons[i][1] + MIN_INTRON <= exons[i + 1][0]
                                                for i in range(len(exons) - 1))


def extend_to_start(s: str, exons: List[Tuple[int, int]], gencode: int, lo: int, qmiss_n: int
                    ) -> Tuple[List[Tuple[int, int]], bool]:
    """Find the start codon. Used when the alignment stopped short of the homolog's N-terminus
    or the first codon is not a start codon: upstream in frame (no stop codon, not before the
    previous feature, lo) the start codon closest to where the homolog's first residue would be;
    if none, the nearest start codon a few codons downstream. If the open frame runs off the
    5' end of the sequence the CDS is 5'-partial (returns partial=True)."""
    table = CodonTable.unambiguous_dna_by_id[gencode]
    a, b = exons[0]
    if qmiss_n <= 1 and s[a - 1:a + 2] in table.start_codons:
        return exons, False
    target = a - 3 * max(qmiss_n, 0)
    best, ran_off = None, False
    for i in range(1, max(qmiss_n + 6, 30)):
        c0 = a - 1 - 3 * i
        if c0 < 0:
            ran_off = lo <= 1
            break
        if c0 < lo - 1 or s[c0:c0 + 3] in table.stop_codons:
            break
        if s[c0:c0 + 3] in table.start_codons and (best is None or abs(c0 + 1 - target) < abs(best - target)):
            best = c0 + 1
    if best is None and s[a - 1:a + 2] in table.start_codons:
        return exons, False
    if best is None and ran_off:
        c0 = (a - 1) % 3                          # first in-frame base of the sequence
        return [(c0 + 1, b)] + exons[1:], True
    if best is None:
        for i in range(1, 11):                    # downstream, inside the first exon
            c0 = a - 1 + 3 * i
            if c0 + 3 > b or s[c0:c0 + 3] in table.stop_codons:
                break
            if s[c0:c0 + 3] in table.start_codons:
                best = c0 + 1
                break
    return ([(best, b)] + exons[1:] if best else exons), False


def prev_feature_end(rec, m: Model, exons_tx: List[Tuple[int, int]], L: int) -> int:
    """Transcript position just after the nearest upstream annotated feature (other genes)."""
    start_tx = exons_tx[0][0]
    g_span = to_genome([(start_tx, exons_tx[-1][1])], m.strand, L)[0]
    lo = 1
    for ft in rec.features:
        if ft.type not in ("gene", "CDS", "tRNA", "rRNA"):
            continue
        # each part: a join's start / end are its min / max, so a neighbour across the origin or with this gene in
        # its intron covered the whole span and was skipped (Grok 2.2, 2026-11-05)
        parts = [(int(p.start) + 1, int(p.end)) for p in ft.location.parts]
        if sum(overlap(x, g_span) for x in parts) >= 30:
            continue
        for fs, fe in parts:
            end_tx = fe if m.strand == 1 else L + 1 - fs
            if end_tx < start_tx:
                lo = max(lo, end_tx + 1)
    return lo


def next_feature_start(rec, m: Model, exons_tx: List[Tuple[int, int]], L: int) -> int:
    """Transcript position just before the nearest downstream annotated feature (other genes)."""
    end_tx = exons_tx[-1][1]
    g_span = to_genome([(exons_tx[0][0], end_tx)], m.strand, L)[0]
    limit = L
    for ft in rec.features:
        if ft.type not in ("gene", "CDS", "tRNA", "rRNA"):
            continue
        parts = [(int(p.start) + 1, int(p.end)) for p in ft.location.parts]     # (as prev_feature_end)
        if sum(overlap(x, g_span) for x in parts) >= 30:  # the gene being re-modelled (or a nested ORF)
            continue
        for fs, fe in parts:
            start_tx = fs if m.strand == 1 else L + 1 - fe
            if start_tx > end_tx:
                limit = min(limit, start_tx - 1)
    return limit


def to_genome(exons: List[Tuple[int, int]], strand: int, L: int) -> List[Tuple[int, int]]:
    return list(exons) if strand == 1 else sorted((L + 1 - b, L + 1 - a) for a, b in exons)


def overlap(a: Tuple[int, int], b: Tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)
