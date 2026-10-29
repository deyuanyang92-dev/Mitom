"""Evidence for every NCBI intron of the gold set (criteria: gold/README.md).

Introns: the truth of benchmark/eval_introns.py (intron features + same-strand CDS-join gaps >= 142 nt; this includes
the cis introns of trans-spliced genes and the rRNA intron). Per intron:
  1 product   the host CDS spliced as NCBI gives it: internal stops, start (or 5' partial), stop (complete / T / TA /
              3' partial), length, agreement with /translation;
  2 junction  the spliced protein aligned (semi-global, BLOSUM62, gap -11/-1) to <= 12 proteins of the same gene from
              other species (refpool_metazoa + gold set; without an intron in that gene first, then by lineage):
              homologs with an insertion / deletion within +-5 aa of the junction, identity before / after it;
  3 signals   intron type by the Rfam + subtype covariance models (mito_intron.intron_types, E <= 0.01, all
              introns of a record in one call as MitoM runs it; --cm_only) and NCBI's note; the window of placements giving the identical spliced sequence; canonical placements in it
              (group I: exon 1 ends T, intron ends G; group II: intron starts GYGYG - GUGYG or GCGCG -, ends AY);
              alternative placements (5' / 3' ends moved independently by <= 60 nt, frame kept) that meet a
              canonical rule (or the group II 5' motif alone) and give a CDS without internal stop, scored against the
              same homologs (a placement giving the same protein through a synonymous change is kept as such);
  4 clade     junction position in a reference protein (Metridium senile NC_000933.1) and phase.
usage: verify.py [--jobs 16] [--acc A,B]   -> gold/work/evidence.pkl
"""
import argparse
import csv
import pickle
import re
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path
from statistics import mean, median
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM")
import eval_introns as E  # noqa: E402
from mito_names import gene_name  # noqa: E402
from Bio import SeqIO  # noqa: E402
from Bio.Align import PairwiseAligner, substitution_matrices  # noqa: E402
from Bio.Data import CodonTable  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from Bio.SeqFeature import BeforePosition, AfterPosition  # noqa: E402

GOLD = Path(__file__).resolve().parents[1]
SET = GOLD.parent
CMS = [Path("/media/deyuan/18T_B/yzz/MitoM/models/rfam_group_I_II_introns.cm")] + \
    sorted((Path.home() / ".cache/mito_intron/selles_vidal_2025_subtype_CMs").glob("*.cm"))
CMSEARCH = "/home/deyuan/miniforge3/envs/mitotools/bin/cmsearch"
REF_ACC = "NC_000933.1"            # Metridium senile: reference protein coordinates for the clade comparison
N_HOM, SHIFT, WIN_AA = 12, 60, 5
GROUP_II_5, GROUP_II_3 = re.compile(r"^G[CT]G[CT]G"), re.compile(r"A[CT]$")   # GUGYG, and GCGCG (annelid cox1)
COMP = str.maketrans("ACGTN", "TGCAN")

AL = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                     open_gap_score=-11, extend_gap_score=-1)
AL.end_gap_score = 0                 # semi-global: overhangs free (partial genes, different termini)


# ── gene models ──────────────────────────────────────────────────────────────
class Model:
    """A CDS / rRNA as exons (start, length, strand) in biological order on the circular (or linear) record."""

    def __init__(self, seq: str, n: int, exons: List[Tuple[int, int, int]], code: int, codon_start: int = 1):
        self.seq2, self.n, self.exons, self.code, self.cs = seq + seq, n, exons, code, codon_start

    def ext(self, start: int, length: int, strand: int) -> str:
        s = self.seq2[start % self.n:start % self.n + length]
        return s if strand == 1 else s[::-1].translate(COMP)

    def spliced(self, exons=None) -> str:
        return "".join(self.ext(*x) for x in (exons or self.exons))

    def intron(self, k: int, exons=None) -> Tuple[int, int, int]:
        """Genome interval (start mod n, length) and strand of the gap between exons k and k+1 (same strand)."""
        ex = exons or self.exons
        (a, la, st), (b, lb, _) = ex[k], ex[k + 1]
        if st == 1:
            return (a + la) % self.n, (b - (a + la)) % self.n, 1
        return (b + lb) % self.n, (a - (b + lb)) % self.n, -1

    def moved(self, k: int, d5: int, d3: int) -> Optional[list]:
        """Exons with intron k's 5' end moved by d5 and 3' end by d3 nt (transcript direction)."""
        ex = list(self.exons)
        (a, la, st), (b, lb, _) = ex[k], ex[k + 1]
        if st == 1:
            ex[k], ex[k + 1] = (a, la + d5, st), ((b + d3) % self.n, lb - d3, st)
        else:
            ex[k], ex[k + 1] = ((a - d5) % self.n, la + d5, st), (b, lb - d3, st)
        return ex if ex[k][1] > 0 and ex[k + 1][1] > 0 else None

    def protein(self, exons=None) -> str:
        nt = self.spliced(exons)[self.cs - 1:]
        return str(Seq(nt[:len(nt) - len(nt) % 3]).translate(table=self.code))

    def junction_nt(self, k: int, exons=None) -> int:
        return sum(x[1] for x in (exons or self.exons)[:k + 1]) - (self.cs - 1)


def window(m: Model, k: int) -> Tuple[int, int]:
    """Genome shifts g of intron k that keep the spliced sequence identical: [gmin, gmax]."""
    s, ln, _ = m.intron(k)
    S, gmax, gmin = m.seq2, 0, 0
    while gmax < ln - 1 and S[(s + gmax) % m.n] == S[(s + ln + gmax) % m.n]:
        gmax += 1
    while -gmin < ln - 1 and S[(s + gmin - 1) % m.n] == S[(s + ln + gmin - 1) % m.n]:
        gmin -= 1
    return gmin, gmax


def signals(m: Model, s: int, ln: int, strand: int) -> Tuple[str, str, str]:
    """(last exon-1 base, first 6 and last 6 intron bases) in transcript orientation."""
    intr = m.ext(s, ln, strand)
    last = m.ext(s - 1, 1, 1) if strand == 1 else m.ext((s + ln) % m.n, 1, -1)
    return last, intr[:6], intr[-6:]


def canonical(last: str, i5: str, i3: str) -> List[str]:
    """Rules a placement meets: I (exon 1 ends T, intron ends G), II (intron starts GYGYG, ends AY),
    II5 (group II 5' motif GYGYG only)."""
    out = []
    if last == "T" and i3.endswith("G"):
        out.append("I")
    if GROUP_II_5.match(i5):
        out.append("II" if GROUP_II_3.search(i3) else "II5")
    return out


# ── homologs ─────────────────────────────────────────────────────────────────
HOM: list = []


def homologs(gene: str, organism: str, lineage: List[str], kind: str = "prot") -> List[Tuple[str, str, bool]]:
    """Closest proteins (or rRNAs) of other species: without an intron in that gene first, then by shared lineage."""
    pool = []
    for x in HOM:
        v = x[kind].get(gene)
        if not v or x["organism"] == organism:
            continue
        depth = 0
        for a, b in zip(x["lineage"], lineage):
            if a != b:
                break
            depth += 1
        pool.append((v[1], -depth, x["acc"], v[0]))
    out, seen = [], set()
    for has, _, acc, seq in sorted(pool):
        if seq not in seen:
            seen.add(seq)
            out.append((acc, seq, has))
        if len(out) >= N_HOM:
            break
    return out


def junction_check(prot: str, q: float, homs: List[Tuple[str, str, bool]]) -> dict:
    """Indels within +-WIN_AA of junction q (aa, float) per homolog, identity in 15 aa before / after it."""
    lo, hi = int(q) - WIN_AA, int(q + 0.999) + WIN_AA
    n_indel, nets, idb, ida = 0, [], [], []
    for _, h, _ in homs:
        aln = AL.align(prot, h)[0]
        blocks = list(zip(aln.aligned[0], aln.aligned[1]))
        net, hit = 0, False
        for (pa, ha), (pb, hb) in zip(blocks, blocks[1:]):
            ins, dele = pb[0] - pa[1], hb[0] - ha[1]          # residues only in prot / only in homolog
            if (ins or dele) and lo <= pa[1] <= hi:
                hit, net = True, net + ins - dele
        n_indel += hit
        nets.append(net)
        same = {"b": [0, 0], "a": [0, 0]}
        for (p0, p1), (h0, h1) in blocks:
            for i in range(p1 - p0):
                pi = p0 + i
                side = "b" if q - 15 <= pi < q else "a" if q <= pi < q + 15 else ""
                if side:
                    same[side][0] += prot[pi] == h[h0 + i]
                    same[side][1] += 1
        idb.append(same["b"][0] / 15)
        ida.append(same["a"][0] / 15)
    return {"n": len(homs), "indel": n_indel, "net": median(nets) if nets else 0,
            "id_before": round(mean(idb), 2) if idb else 0, "id_after": round(mean(ida), 2) if ida else 0}


def scores(prot: str, homs) -> List[float]:
    return [AL.score(prot.rstrip("*"), h) for _, h, _ in homs]


# ── product ──────────────────────────────────────────────────────────────────
def product(m: Model, feat, partial5: bool, partial3: bool) -> dict:
    table = CodonTable.unambiguous_dna_by_id[m.code]
    nt = m.spliced()[m.cs - 1:]
    tail = len(nt) % 3
    prot = str(Seq(nt[:len(nt) - tail]).translate(table=m.code))
    body = prot[:-1] if prot.endswith("*") and not tail else prot
    stops = body.count("*")
    start = "5' partial" if partial5 or m.cs != 1 else ("ok" if nt[:3] in table.start_codons else f"bad ({nt[:3]})")
    if partial3:
        stop = "3' partial"
    elif not tail and prot.endswith("*"):
        stop = "complete"
    elif tail and nt[-tail:] in ("T", "TA"):
        stop = nt[-tail:]
    else:
        stop = f"none ({nt[-3:]})"
    given = feat.qualifiers.get("translation", [""])[0] if feat is not None else ""
    ours = ("M" + body[1:] if start == "ok" else body).replace("*", "")
    agree = "no /translation" if not given else "same" if given == ours else \
        f"differs ({sum(a != b for a, b in zip(given, ours)) + abs(len(given) - len(ours))} aa)"
    return {"stops": stops, "start": start, "stop": stop, "aa": len(body.rstrip("*")), "translation": agree,
            "ok": stops == 0 and not start.startswith("bad") and not stop.startswith("none")}


# ── per intron ───────────────────────────────────────────────────────────────
def host_of(rec, t: E.Intron):
    """(feature, gap index) of the CDS / rRNA whose same-strand gap is this intron (>= 50 % reciprocal overlap)."""
    n = len(rec)
    best = None
    for f in rec.features:
        if f.type not in ("CDS", "rRNA") or len(f.location.parts) < 2:
            continue
        for k, (p, q) in enumerate(zip(f.location.parts, f.location.parts[1:])):
            if p.strand != q.strand:
                continue
            plus = (p.strand or 1) == 1
            s, e = (int(p.end), int(q.start)) if plus else (int(q.end), int(p.start))
            if (int(q.start) < int(p.start)) == plus:
                e += n
            ov = E._overlap((s, e), (t.start, t.end), n)
            if (p.strand or 1) == t.strand and ov >= 0.5 * max(e - s, t.end - t.start):
                if best is None or ov > best[2]:
                    best = (f, k, ov)
    return (best[0], best[1]) if best else (None, None)


def model_of(rec, f, code: int) -> Model:
    seq = str(rec.seq).upper()
    exons = [(int(p.start), int(p.end) - int(p.start), p.strand or 1) for p in f.location.parts]
    return Model(seq, len(seq), exons, int(f.qualifiers.get("transl_table", [code])[0]),
                 int(f.qualifiers.get("codon_start", ["1"])[0]))


def note_type(text: str) -> str:
    t = text.lower()
    if re.search(r"group\s*(ii|2)\b", t):
        return "II"
    if re.search(r"group\s*(i|1)\b", t):
        return "I"
    return ""


def analyse(job) -> List[dict]:
    acc, info = job
    rec = SeqIO.read(str(SET / "ncbi_gb" / f"{acc}.gb"), "genbank")
    code, n = int(info["genetic_code"]), len(rec)
    lineage, organism = rec.annotations.get("taxonomy", []), rec.annotations.get("organism", "")
    ref = {g: p for x in HOM if x["acc"] == REF_ACC for g, (p, _) in x["prot"].items()}
    out = []
    for t in sorted(E.introns_of(rec), key=lambda x: x.start):
        r = {"accession": acc, "phylum": info["phylum"], "order": info["order"], "organism": organism, "code": code,
             "n": n, "gene": t.gene, "strand": t.strand, "ncbi_s": t.start, "ncbi_e": t.end, "source": t.source,
             "ncbi_note": t.note, "note_type": note_type(t.note), "fragment": info["fragment"] == "yes"}
        f, k = host_of(rec, t)
        if f is not None and gene_name(f.qualifiers, f.type) and gene_name(f.qualifiers, f.type) != t.gene:
            r["ncbi_gene"], r["gene"] = t.gene, gene_name(f.qualifiers, f.type)   # e.g. NC_086801.1: ND5 intron, /gene="COX1"
        r["host"] = "" if f is None else f"{f.type} {f.location}"
        r["trans_spliced"] = f is not None and len({p.strand for p in f.location.parts}) > 1
        r["host_kind"] = "none" if f is None else f.type
        if f is None:
            # intron feature without a matching join gap: what contains it?
            r["intron_seq"] = Model(str(rec.seq).upper(), n, [], code).ext(t.start, t.end - t.start, t.strand)
            r["containing"] = ", ".join(sorted({f"{x.type} {x.qualifiers.get('gene', x.qualifiers.get('product', ['?']))[0]}"
                                               for x in rec.features if x.type in ("CDS", "rRNA", "tRNA", "gene")
                                               and E._overlap((int(x.location.start), int(x.location.end)),
                                                              (t.start, t.end), n)}))
            out.append(r)
            continue
        m = model_of(rec, f, code)
        r["k"], r["n_exons"] = k, len(m.exons)
        s, ln, strand = m.intron(k)
        r["gap_s"], r["gap_len"] = s, ln            # the join gap (may differ from an intron feature by a few nt)
        r["intron_seq"] = m.ext(s, ln, strand)
        last, i5, i3 = signals(m, s, ln, strand)
        r.update(sig_last=last, sig_5=i5, sig_3=i3, sig_kinds=canonical(last, i5, i3))
        gmin, gmax = window(m, k)
        r["window"] = (gmin, gmax)
        r["window_sig"] = [(g, *signals(m, (s + g) % n, ln, strand), canonical(*signals(m, (s + g) % n, ln, strand)))
                           for g in range(gmin, gmax + 1)]
        if f.type == "rRNA":
            rh = homologs(r["gene"], organism, lineage, "rrna")
            r["homologs"] = [(a_, has) for a_, _, has in rh]
            r["rrna_junction"] = rrna_junction(m.spliced(), m.junction_nt(k), rh)
            out.append(r)
            continue
        loc = f.location
        p5 = isinstance(loc.parts[0].start if (loc.parts[0].strand or 1) == 1 else loc.parts[0].end,
                        (BeforePosition, AfterPosition))
        p3 = isinstance(loc.parts[-1].end if (loc.parts[-1].strand or 1) == 1 else loc.parts[-1].start,
                        (BeforePosition, AfterPosition))
        r["product"] = product(m, f, p5, p3)
        prot = m.protein()
        homs = homologs(r["gene"], organism, lineage)
        r["homologs"] = [(a, has) for a, _, has in homs]
        q = m.junction_nt(k) / 3
        r["junction_aa"], r["phase"] = q, m.junction_nt(k) % 3
        r["junction"] = junction_check(prot.rstrip("*"), q, homs)
        base = scores(prot, homs)
        r["score"], r["scores"] = round(mean(base), 1), base
        r["ctx"] = prot[max(0, int(q) - 5):int(q)] + "|" + prot[int(q):int(q) + 5]
        refp = ref.get(r["gene"])
        if refp:
            r["ref_pos"] = ref_position(AL.align(prot.rstrip("*"), refp)[0], q)
        # alternatives: both ends moved independently (<= SHIFT nt), frame kept, meeting a canonical rule (or the
        # group II 5' motif), CDS without internal stop; the same protein from another placement is kept (synonymous)
        alts: Dict[str, dict] = {}
        w0 = set(range(gmin, gmax + 1))
        for d5 in range(-SHIFT, SHIFT + 1):
            for d3 in range(-SHIFT, SHIFT + 1):
                if (d5 - d3) % 3 or (d5 == d3 and (d5 if strand == 1 else -d5) in w0):
                    continue
                ex = m.moved(k, d5, d3)
                if ex is None:
                    continue
                s2, ln2, _ = m.intron(k, ex)
                if ln2 < 50:
                    continue
                kinds = canonical(*signals(m, s2, ln2, strand))
                if not kinds:
                    continue
                p2 = m.protein(ex)
                if "*" in p2.rstrip("*"):
                    continue
                alts.setdefault(p2, []).append({"d5": d5, "d3": d3, "kinds": kinds, "start": s2, "len": ln2,
                                                "q": m.junction_nt(k, ex) / 3})
        ranked = []
        for p2, pl in alts.items():
            pl.sort(key=lambda x: (abs(x["d5"]) + abs(x["d3"]), x["d5"]))
            sc = base if p2 == prot else scores(p2, homs)
            ranked.append({"placements": pl[:6], "n_placements": len(pl), "same_protein": p2 == prot,
                           "score": round(mean(sc), 1), "scores": sc,
                           "wins": sum(x > y for x, y in zip(sc, base)), "losses": sum(x < y for x, y in zip(sc, base)),
                           "aa_diff": aa_diff(prot, p2), "prot": p2})
        ranked.sort(key=lambda x: (-x["score"], not x["same_protein"]))
        keep = [a for a in ranked if a["same_protein"]] + [a for a in ranked if not a["same_protein"]][:8]
        for a in keep:
            q2 = a["placements"][0]["q"]
            a["junction"] = r["junction"] if a["same_protein"] else junction_check(a["prot"].rstrip("*"), q2, homs)
            if refp:
                a["ref_pos"] = ref_position(AL.align(a["prot"].rstrip("*"), refp)[0], q2)
            a.pop("prot")
        r["alts"] = keep
        r["exons"] = m.exons
        out.append(r)
    return out


def aa_diff(a: str, b: str) -> str:
    if len(a) != len(b):
        return f"{len(b) - len(a):+d} aa"
    return f"{sum(x != y for x, y in zip(a, b))} subst"


def ref_position(aln, q: float) -> int:
    """Reference residue index aligned to the protein residue at the junction (nearest aligned one)."""
    j, best = int(q), None
    for (p0, p1), (h0, h1) in zip(*aln.aligned):
        for i in range(p1 - p0):
            d = abs(p0 + i - j)
            if best is None or d < best[0]:
                best = (d, h0 + i)
    return best[1] if best else -1


def rrna_junction(sp: str, j: int, homs) -> dict:
    """Nucleotide alignment of the spliced rRNA to homolog rRNAs: indels within +-15 nt of the junction."""
    al = PairwiseAligner(mode="global", match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2)
    al.end_gap_score = 0
    hit, nets = 0, []
    for _, h, _ in homs[:6]:
        aln = al.align(sp, h)[0]
        blocks = list(zip(aln.aligned[0], aln.aligned[1]))
        net, any_ = 0, False
        for (pa, ha), (pb, hb) in zip(blocks, blocks[1:]):
            ins, dele = pb[0] - pa[1], hb[0] - ha[1]
            if (ins or dele) and j - 15 <= pa[1] <= j + 15:
                any_, net = True, net + ins - dele
        hit += any_
        nets.append(net)
    return {"n": min(6, len(homs)), "indel": hit, "net": median(nets) if nets else 0}


# ── covariance-model typing ──────────────────────────────────────────────────
def cm_chunk(seqs: Dict[str, str]) -> Dict[str, tuple]:
    import mito_intron as MI
    with tempfile.TemporaryDirectory() as td:
        return MI.intron_types(CMSEARCH, CMS, seqs, 0.01, Path(td))


def init(h):
    global HOM
    HOM = h


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--jobs", type=int, default=16)
    p.add_argument("--acc", default="")
    p.add_argument("--no_cm", action="store_true")
    p.add_argument("--cm_only", action="store_true",
                   help="only (re)type the introns of an existing evidence.pkl, one cmsearch call per record")
    a = p.parse_args()
    if a.cm_only:
        path = GOLD / "work" / "evidence.pkl"
        rows = pickle.load(open(path, "rb"))         # written by this script (this folder)
        by: Dict[str, Dict[str, str]] = {}
        for i, r in enumerate(rows):
            by.setdefault(r["accession"], {})[f"i{i}"] = r["intron_seq"]
        with Pool(a.jobs) as pool:
            typed = {}
            for d in pool.imap_unordered(cm_chunk, sorted(by.values(), key=lambda d: -sum(map(len, d.values())))):
                typed.update(d)
        for i, r in enumerate(rows):
            r["cm"] = typed.get(f"i{i}", ("", None, ""))
        pickle.dump(rows, open(path, "wb"))
        print(f"{len(rows)} introns retyped per record -> {path}", file=sys.stderr)
        return 0
    info = {r["accession"]: r for r in csv.DictReader(open(SET / "records.tsv"), delimiter="\t")}
    accs = a.acc.split(",") if a.acc else sorted(info)
    hom = pickle.load(open(GOLD / "work" / "homologs.pkl", "rb"))     # written by homologs.py (this folder)
    with Pool(a.jobs, initializer=init, initargs=(hom,)) as pool:
        rows = [r for part in pool.imap(analyse, [(x, info[x]) for x in accs]) for r in part]
        if not a.no_cm:
            seqs = {f"i{i}": r["intron_seq"] for i, r in enumerate(rows)}
            keys = sorted(seqs, key=lambda x: -len(seqs[x]))
            chunks = [{x: seqs[x] for x in keys[j::a.jobs]} for j in range(a.jobs)]
            typed = {}
            for d in pool.imap_unordered(cm_chunk, chunks):
                typed.update(d)
            for i, r in enumerate(rows):
                r["cm"] = typed.get(f"i{i}", ("", None, ""))
    out = GOLD / "work" / ("evidence.pkl" if not a.acc else "evidence_test.pkl")
    pickle.dump(rows, open(out, "wb"))
    print(f"{len(rows)} introns -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
