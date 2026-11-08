"""mito_nuc_lib.py - shared part of mito_nuc.py (`mitom nuc`): settings, the external programs, the NCBI reference records
(fetch) and the reference library (build, Library), coordinate helpers of the rDNA and histone searches."""
from __future__ import annotations
import csv
import gzip
import logging
import re
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import BeforePosition, AfterPosition
import mito_env
from mito_common import run_tool


HERE = Path(__file__).resolve().parent


LIB = HERE / "models" / "nuc"


RRNA = ("18S", "5.8S", "28S")


UNIT = ("18S", "ITS1", "5.8S", "ITS2", "28S")


HISTONES = ("H1", "H2A", "H2B", "H3", "H4")


PRODUCT = {"18S": "18S ribosomal RNA", "5.8S": "5.8S ribosomal RNA", "28S": "28S ribosomal RNA",
           "ITS1": "internal transcribed spacer 1", "ITS2": "internal transcribed spacer 2",
           **{h: f"histone {h}" for h in HISTONES}}


MAX_PER_SPECIES = 3  # rDNA reference records per species (NCBI has up to ~100 of one species)


TOPN = 5            # references that vote on an rDNA boundary


PREFILTER_NT = 5_000_000   # inputs larger than this go through the prefilter


CHUNK_NT = 300_000  # query nt per blastn chunk (memory bound of the rDNA search)


BLAST_TARGETS = 100   # reference units kept per region (blastn; 100 vs 500: same gold result, 2 of 992 genes differ)


EXTRAP = 30         # projection past an HSP end that still counts as a vote (nt)


MIN_CHAIN_BITS = 100


MIN_GENE_OVERLAP = 80   # nt of a reference rRNA an alignment chain must cover for the gene to be reported


PIECE_MIN_ID = 0.7      # identity of an unaligned, projected piece at a contig end to the reference gene's end


MAX_ITS = 3000      # spacer between an rRNA and the contig end that is still reported as a partial ITS


MAX_RESCUE = 3000   # stretch beside an rRNA (or a whole sequence) searched again for a missing 5.8S (blastn word 11)


CURATED_MIN_ID = 90 # a curated reference decides boundaries only from chains of >= this % identity (close relatives)


QC_TOL = 15         # library QC: a reference boundary this far (nt) from the one other studies map onto it ...


QC_MIN_STUDIES = 3  # ... supported by >= this many other studies (and > 2x the studies that agree with it) is not used


TOPK = 3            # histone proteins tried with exonerate per locus


INTRON_PEN = 0.9


MIN_PROT_COVER = 0.3    # histone model must cover this fraction of the protein unless it reaches a contig end


MIN_H34_PARTIAL_ID = 80 # % identity (exonerate) below which a partial H3 / H4 model is a degenerate copy, not reported


SEED_COV = 0.1      # consensus seed: fragments with >= this fraction of the best fragment's contig coverage


END_TOL = 20        # consensus: a copy variant this close (nt) to a contig end does not keep it from joining the next


PAD = 20000         # flank around prefilter hits taken from long contigs


DB_PIECE = 1_000_000    # min. FASTA bytes per makeblastdb piece of the prefilter database (one piece per thread)


# a bounded amount of work per assembly, whatever the number of candidates (real genome skims had <= 48 rDNA and
# <= 3,737 histone regions, 2026-10-29); only where the candidates are < 10 % of the input sequences (an assembly,
# not a batch of records)
MAX_RDNA_REGIONS = 1000 # rDNA regions annotated, most aligned reference bits first


MAX_HIST_REGIONS = 500  # histone regions per gene, closest to a reference protein first


MAX_CM_JOBS = 100       # rRNA ends per gene the covariance model places (+-150 nt), best supported genes first; the


                        # others and the +-600 nt windows (~45 s each for an 18S end) keep the HMM / projection
MAX_CONS_SEED = 30      # best-covered fragments per gene the consensus backbone is built from (real assemblies: <= 8)


STOPS = ("TAA", "TAG", "TGA")


BAD_HISTONE = re.compile(r"factor|transferase|acetyl|chaperone|binding|kinase|demethylase|ligase|DOT1|lysine|"
                         r"regulator|interact|domain|methyl", re.I)


QUERIES = {
    "rdna": 'txid33208[Organism:exp] AND 18S[Title] AND 5.8S[Title] AND 28S[Title] AND "complete sequence"[Title] '
            'NOT partial[Title]',
    # curated RefSeq rRNA whose ends agree between species (excluded: mouse NR_003278-80, 18S shifted 1 nt and 28S /
    # 5.8S 3' ends unlike human; hamster 18S +4 nt; pig 18S 2302 nt; horse 5'-truncated or fragment; human N4 5.8S,
    # rat 5.8S NR_046238 and macaque NR_146159 1 nt short / longer; Drosophila 28S hidden break;
    # Saccoglossus unusual ends)
    "curated": " OR ".join(f"{a}[Accession]" for a in (
        "NR_003286", "NR_003285", "NR_003287", "NR_145820", "NR_145821", "NR_145822", "NR_146146", "NR_146147",
        "NR_146148", "NR_146152", "NR_146153", "NR_146154", "NR_146119", "NR_146118",       # human N1-N5
        "NR_046237", "NR_046246",                                                            # rat
        "NR_146166", "NR_146160", "NR_033238", "NR_036642", "NR_036643", "NR_045133",       # macaque, rabbit, cow, hamster
        "NR_133559", "NR_133555", "NR_133550", "NR_000053", "NR_000054", "NR_000057")),     # Drosophila, C. elegans
    "histone": 'txid33208[Organism:exp] AND ("histone H3"[Title] OR "histone H4"[Title] OR "histone H2A"[Title] OR '
               '"histone H2B"[Title] OR "histone H1"[Title]) AND "complete cds"[Title] AND 300:20000[SLEN]',
}


def tool(name: str) -> str:
    p = mito_env.find(name)
    if not p:
        raise SystemExit(f"{name} not found (mitom doctor; it is in the mitotools environment)")
    return p


def tool_ok(name: str) -> bool:
    return bool(mito_env.find(name))


def program_output(cmd: List[str]) -> str:
    """stdout of a program; its failure raises RuntimeError with the end of its stderr."""
    p = run_tool(cmd)                             # [CMD] line; the program's own error text kept (Kimi 1.8, 2026-11-05)
    if p.returncode:
        raise RuntimeError(f"{Path(cmd[0]).name} failed ({p.returncode}): {p.stderr.strip()[-400:]}")
    return p.stdout


# ── reference library ───────────────────────────────────────────────────────
def _label(f: SeqFeature) -> str:
    q = f.qualifiers
    return " ".join(q.get("product", []) + q.get("note", []) + q.get("gene", []))


def rrna_name(f: SeqFeature) -> Optional[str]:
    lab = _label(f)
    if f.type == "rRNA":
        for key, name in (("5.8S", "5.8S"), ("18S", "18S"), ("28S", "28S"),
                          ("small subunit", "18S"), ("large subunit", "28S")):
            if key in lab:
                return name
    return None


def histone_name(f: SeqFeature) -> Optional[str]:
    if f.type != "CDS" or BAD_HISTONE.search(_label(f)):
        return None
    m = re.search(r"histone (H2A|H2B|H3|H4|H1)\b", _label(f), re.I)
    return m.group(1).upper() if m else None


def _partial(f: SeqFeature) -> Tuple[bool, bool]:
    """(5' partial, 3' partial) of a feature, on its own strand"""
    lo, hi = isinstance(f.location.start, BeforePosition), isinstance(f.location.end, AfterPosition)
    return (lo, hi) if f.location.strand != -1 else (hi, lo)


def _study(rec) -> str:
    """first publication title of a record (records of one submission share it)"""
    refs = rec.annotations.get("references", [])
    t = [r.title for r in refs if r.title and r.title != "Direct Submission"]
    return re.sub(r"\s+", " ", t[0] if t else refs[-1].journal if refs else rec.id)[:150]


def build(gbs: List[Path], out: Path, threads: int = 8, curated: List[Path] = ()) -> Tuple[int, int]:
    """rdna_ref.fa.gz + rdna_ref.tsv: each record with >= 1 rRNA (<= MAX_PER_SPECIES per species, first ones in
    the files), oriented to the rDNA strand; ends marked partial are kept but not used as boundaries, nor are
    boundaries that other studies contradict (qc_rdna). Records of the curated files are trusted: no species cap,
    no QC, and where one aligns over a boundary its position is used instead of the vote.
    hist_ref.faa + hist_ref.tsv: complete histone CDS (<= 400 aa)."""
    out.mkdir(parents=True, exist_ok=True)
    n_r = n_h = 0
    with gzip.open(out / "rdna_ref.fa.gz", "wt") as rfa, open(out / "rdna_ref.tsv", "w", newline="") as rtsv, \
            open(out / "hist_ref.faa", "w") as hfa, open(out / "hist_ref.tsv", "w", newline="") as htsv:
        rw, hw = csv.writer(rtsv, delimiter="\t"), csv.writer(htsv, delimiter="\t")
        rw.writerow(["ref", "organism", "genus", "taxonomy", "study", "gene", "start", "end", "start_ok", "end_ok",
                     "curated"])
        hw.writerow(["ref", "accession", "organism", "genus", "taxonomy", "gene", "length_aa"])
        seen, per_species = set(), Counter()
        for gb in [*curated, *gbs]:
            cur = gb in curated
            for rec in SeqIO.parse(str(gb), "genbank"):
                if rec.id in seen:
                    continue
                seen.add(rec.id)
                org = rec.annotations.get("organism", "unknown")
                genus, tax = org.split()[0], ";".join(rec.annotations.get("taxonomy", [])[:6])
                rr = defaultdict(list)
                for f in rec.features:
                    g = rrna_name(f)
                    if g:
                        rr[g].append(f)
                rr = {g: fs[0] for g, fs in rr.items() if len(fs) == 1}
                strands = {f.location.strand for f in rr.values()}
                if rr and len(strands) == 1 and (cur or per_species[org] < MAX_PER_SPECIES):
                    per_species[org] += not cur
                    minus, n = strands.pop() == -1, len(rec)
                    rfa.write(f">{rec.id}\n{rec.seq.reverse_complement() if minus else rec.seq}\n")
                    for g, f in rr.items():
                        a, b = int(f.location.start) + 1, int(f.location.end)
                        p5, p3 = _partial(f)
                        if minus:
                            a, b = n - b + 1, n - a + 1
                        rw.writerow([rec.id, org, genus, tax, _study(rec), g, a, b, int(not p5), int(not p3), int(cur)])
                    n_r += 1
                for f in rec.features:
                    h = histone_name(f)
                    if not h or "translation" not in f.qualifiers or any(_partial(f)) or \
                            int(f.qualifiers.get("codon_start", ["1"])[0]) != 1:
                        continue
                    prot = f.qualifiers["translation"][0]
                    if len(prot) > 400:
                        continue
                    rid = f"{rec.id}_{int(f.location.start) + 1}"
                    hfa.write(f">{rid}\n{prot}\n")
                    hw.writerow([rid, rec.id, org, genus, tax, h, len(prot)])
                    n_h += 1
    qc_rdna(out, threads)
    return n_r, n_h


def qc_rdna(lib: Path, threads: int) -> None:
    """Submitters of one study often share a wrong boundary (e.g. a 28S 3' end 300-600 nt before the conserved
    ...GATTCGT end); their records then outvote everyone for a related query. Every reference is aligned with
    the references of other studies, their boundaries are mapped onto it, and a boundary is set not ok when the
    position most studies map to is > QC_TOL nt away. Flagged boundaries are listed in rdna_qc.tsv."""
    rows = list(csv.DictReader(open(lib / "rdna_ref.tsv"), delimiter="\t"))
    if not rows:
        return
    ann: Dict[str, Dict[str, tuple]] = defaultdict(dict)
    study, curated = {}, {r["ref"] for r in rows if r["curated"] == "1"}
    for r in rows:
        ann[r["ref"]][r["gene"]] = (int(r["start"]), int(r["end"]), r["start_ok"] == "1", r["end_ok"] == "1")
        study[r["ref"]] = r["study"]
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "ref.fa"
        with gzip.open(lib / "rdna_ref.fa.gz", "rt") as i, open(fa, "w") as o:
            shutil.copyfileobj(i, o)
        program_output([tool("makeblastdb"), "-in", str(fa), "-dbtype", "nucl", "-out", str(Path(td) / "db")])
        tsv = program_output([tool("blastn"), "-task", "megablast", "-word_size", "16", "-query", str(fa),
                   "-db", str(Path(td) / "db"), "-evalue", "1e-10", "-dust", "no", "-max_target_seqs", "1000",
                   "-num_threads", str(threads), "-outfmt", "6 qseqid sseqid qstart qend sstart send bitscore qseq sseq"])
    hits: Dict[str, Dict[str, list]] = defaultdict(lambda: defaultdict(list))
    for line in tsv.splitlines():
        q, s, qs, qe, ss, se, bs, qa, sa = line.split("\t")
        if study[q] != study[s]:
            hits[q][s].append((int(qs), int(qe), int(ss), int(se), float(bs), qa, sa))
    bad = []
    for q in [q for q in ann if q not in curated]:
        votes: Dict[tuple, Dict[int, set]] = defaultdict(lambda: defaultdict(set))   # (gene, side) -> pos -> studies
        for s, hs in hits[q].items():
            c = chains(hs)
            if not c or _strand(c[0][0]) != 1:
                continue
            for g, (a, b, ok_a, ok_b) in ann[s].items():
                for side, pos, ok in ((0, a, ok_a), (1, b, ok_b)):
                    m = next((x for x in (map_pos(h, pos) for h in c[0]) if x is not None), None) if ok else None
                    if m is not None:
                        votes[(g, side)][int(m)].add(study[s])
        for g, v in ann[q].items():
            for side in (0, 1):
                vv = votes[(g, side)]
                if not v[2 + side] or not vv:
                    continue
                def near(p, vv=vv):
                    return set().union(set(), *(st for x, st in vv.items() if abs(x - p) <= 5))
                best = max(vv, key=lambda p: (len(vv[p]), len(near(p))))    # position most studies map to
                n_best, n_own = len(near(best)), len(near(v[side]))
                if abs(best - v[side]) > QC_TOL and n_best >= QC_MIN_STUDIES and n_best > 2 * n_own:
                    bad.append((q, g, side, v[side], best, n_own, n_best))
    flag = {(q, g, side) for q, g, side, *_ in bad}
    with open(lib / "rdna_ref.tsv", "w", newline="") as o:
        w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        for r in rows:
            for side, col in ((0, "start_ok"), (1, "end_ok")):
                if (r["ref"], r["gene"], side) in flag:
                    r[col] = "0"
            w.writerow(r)
    with open(lib / "rdna_qc.tsv", "w") as o:
        o.write("ref\tgene\tside\tannotated\tother_studies\tstudies_agreeing\tstudies_other\n")
        for q, g, side, own, best, n_own, n_best in bad:
            o.write(f"{q}\t{g}\t{('start', 'end')[side]}\t{own}\t{best}\t{n_own}\t{n_best}\n")
    logging.info(f"[NUC] library QC: {len(bad)} reference boundaries contradicted by other studies (rdna_qc.tsv)")


def fetch(out: Path, email: Optional[str], api_key: Optional[str]) -> None:
    import mito_fetch
    out.mkdir(parents=True, exist_ok=True)
    mito_fetch.setup_entrez(email, api_key)
    for name, term in QUERIES.items():
        uids = mito_fetch.search_uids(term)
        gb = out / f"{name}.gb"
        have = mito_fetch.existing_accessions(gb)
        todo = [u for u in uids if u not in have]
        logging.info(f"[NUC] {name}: {len(uids)} records match, {len(todo)} to download")
        if todo:
            mito_fetch.download(todo, gb)
        (out / f"{name}.query.txt").write_text(term + "\n")


class Library:
    def __init__(self, lib: Path, work: Path):
        self.dir = lib
        self.rdna: Dict[str, Dict[str, Tuple[int, int, bool, bool]]] = defaultdict(dict)
        self.org: Dict[str, Tuple[str, str]] = {}
        for r in csv.DictReader(open(lib / "rdna_ref.tsv"), delimiter="\t"):
            self.rdna[r["ref"]][r["gene"]] = (int(r["start"]), int(r["end"]), r["start_ok"] == "1", r["end_ok"] == "1")
            self.org[r["ref"]] = (r["organism"], r["genus"])
        self.curated = {r["ref"] for r in csv.DictReader(open(lib / "rdna_ref.tsv"), delimiter="\t") if r["curated"] == "1"}
        self.hist: Dict[str, dict] = {r["ref"]: r for r in csv.DictReader(open(lib / "hist_ref.tsv"), delimiter="\t")}
        self.prot = {r.id: str(r.seq) for r in SeqIO.parse(str(lib / "hist_ref.faa"), "fasta")}
        work.mkdir(parents=True, exist_ok=True)
        self.rdna_fa = work / "rdna_ref.fa"
        with gzip.open(lib / "rdna_ref.fa.gz", "rt") as i, open(self.rdna_fa, "w") as o:
            shutil.copyfileobj(i, o)
        program_output([tool("makeblastdb"), "-in", str(self.rdna_fa), "-dbtype", "nucl", "-out", str(work / "rdna_ref")])
        self.rdna_db = str(work / "rdna_ref")
        self.curated_db = None
        if self.curated:
            cfa = work / "curated_ref.fa"
            SeqIO.write([r for r in SeqIO.parse(str(self.rdna_fa), "fasta") if r.id in self.curated], str(cfa), "fasta")
            program_output([tool("makeblastdb"), "-in", str(cfa), "-dbtype", "nucl", "-out", str(work / "curated_ref")])
            self.curated_db = str(work / "curated_ref")
        self.hist_faa = lib / "hist_ref.faa"
        # representatives for the prefilter
        seen, reps = set(), set()
        for ref, (_, genus) in self.org.items():
            if genus not in seen:
                seen.add(genus)
                reps.add(ref)
        self.rdna_rep = work / "rdna_rep.fa"
        SeqIO.write([r for r in SeqIO.parse(str(self.rdna_fa), "fasta") if r.id in reps], str(self.rdna_rep), "fasta")
        seen, self.hist_rep = set(), work / "hist_rep.faa"
        with open(self.hist_rep, "w") as o:
            for ref, r in self.hist.items():
                key = (r["gene"], ";".join(r["taxonomy"].split(";")[:4]))
                if key not in seen:
                    seen.add(key)
                    o.write(f">{ref}\n{self.prot[ref]}\n")


# ── rDNA ────────────────────────────────────────────────────────────────────
def _strand(h) -> int:
    return 1 if h[3] >= h[2] else -1


def chains(hsps: list, min_bits: float = MIN_CHAIN_BITS) -> List[list]:
    """collinear HSP chains of one reference, best first; HSPs of a chain do not overlap the query span of a
    better chain (copies on one contig become separate chains)"""
    left, out = list(hsps), []
    while left:
        best, best_sc = [], 0.0
        for st in (1, -1):
            H = sorted([h for h in left if _strand(h) == st], key=lambda h: h[0])
            sc, prev = [h[4] for h in H], [-1] * len(H)
            for i in range(len(H)):
                for j in range(i):
                    # both sequences advance (overlap <= 50 nt): stacked HSPs of one repeat never chain
                    ok = H[j][0] < H[i][0] and H[j][1] < H[i][1] and H[j][1] <= H[i][0] + 50 and \
                        (H[j][2] < H[i][2] and H[j][3] < H[i][3] and H[j][3] <= H[i][2] + 50 if st == 1 else
                         H[j][2] > H[i][2] and H[j][3] > H[i][3] and H[j][3] >= H[i][2] - 50)
                    if ok and sc[j] + H[i][4] > sc[i]:
                        sc[i], prev[i] = sc[j] + H[i][4], j
            if H:
                i = max(range(len(H)), key=lambda k: sc[k])
                if sc[i] > best_sc:
                    c = []
                    while i >= 0:
                        c.append(H[i])
                        i = prev[i]
                    best, best_sc = c[::-1], sc[max(range(len(H)), key=lambda k: sc[k])]
        if best_sc < min_bits:
            break
        out.append(best)
        lo, hi = best[0][0], best[-1][1]
        left = [h for h in left if h not in best and (h[1] < lo or h[0] > hi)]
    return out


def map_pos(h, b: int) -> Optional[float]:
    """reference position b -> query position through the HSP alignment"""
    qs, qe, ss, se, _, qa, sa = h
    st = _strand(h)
    if not min(ss, se) <= b <= max(ss, se):
        return None
    qp, sp = qs - 1, ss - st
    for x, y in zip(qa, sa):
        if x != "-":
            qp += 1
        if y != "-":
            sp += st
            if sp == b:
                return qp if x != "-" else qp + 0.5
    return None


# ── input, regions, output ─────────────────────────────────────────────────
def _contig(region_id: str) -> str:
    return region_id.rsplit(":", 1)[0] if re.search(r":\d+-\d+$", region_id) else region_id


def _cov(name: str) -> str:
    m = re.search(r"cov_([\d.]+)", name)
    return m.group(1) if m else ""


def _ref_gene(lib: "Library", ref: str, g: str) -> Optional[str]:
    if g not in lib.rdna.get(ref, {}):
        return None
    if not hasattr(lib, "_seqs"):
        lib._seqs = SeqIO.index(str(lib.rdna_fa), "fasta")
    a, b = lib.rdna[ref][g][:2]
    return str(lib._seqs[ref].seq[a - 1:b])
