"""mito_intron_sites.py - the intron-site library of mito_intron_ref.py: every intron of annotated GenBank records as a site
(gene, kind, exon lengths, junction probes) -> models/intron_sites.tsv.gz + .fa.gz (build_library, save_library,
load_library), the intron groups typed by the Rfam covariance models (type_introns)."""
from __future__ import annotations
import gzip
import io
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
import mito_names
from mito_common import codon_translate, iter_gb_files, run_tool


HERE = Path(__file__).resolve().parent


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
    codes = [f.qualifiers["transl_table"][0] for f in rec.features if f.type == "CDS" and "transl_table" in f.qualifiers]
    rec_code = max(sorted(set(codes)), key=codes.count) if codes else "1"   # a CDS without one: the record's (Kimi 2.4)
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
        code = int(f.qualifiers.get("transl_table", [rec_code if f.type == "CDS" else "0"])[0] or 0)
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


                                # was a placozoan site put on a bilaterian / medusozoan; a recent horizontal transfer
                                # of a homing intron would keep near-identical junctions)
PHYLA = {"Porifera", "Placozoa", "Cnidaria", "Ctenophora", "Annelida", "Mollusca", "Arthropoda", "Chordata",
         "Echinodermata", "Hemichordata", "Nematoda", "Platyhelminthes", "Brachiopoda", "Bryozoa", "Nemertea",
         "Phoronida", "Rotifera", "Acanthocephala", "Chaetognatha", "Tardigrada", "Onychophora", "Priapulida",
         "Kinorhyncha", "Nematomorpha", "Entoprocta", "Cycliophora", "Gastrotricha", "Xenacoelomorpha", "Dicyemida",
         "Orthonectida"}


def phylum_of(lineage: List[str]) -> str:
    return next((t for t in lineage if t in PHYLA), "")


GROUP_II_5, GROUP_II_3 = re.compile(r"^G[CT]G[CT]G"), re.compile(r"A[CT]$")


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
