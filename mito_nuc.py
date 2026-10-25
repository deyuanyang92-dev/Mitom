#!/usr/bin/env python3
"""mito_nuc.py - nuclear marker genes in contigs: rDNA (18S, ITS1, 5.8S, ITS2, 28S) and histones (H1, H2A, H2B, H3, H4).

Reference-based, like MitoFinder for mitogenomes: the contigs are searched with complete NCBI annotations and
each gene boundary is projected from the reference alignments onto the contig.

rDNA: blastn (megablast, word 16) of the candidate regions against the reference rDNA units; per reference the best
collinear chain of HSPs is kept (several chains = several copies on one contig, e.g. a tandem array); every
annotated reference boundary (18S/5.8S/28S start and end) is mapped through the alignment column by column
(or projected <= EXTRAP nt past the end of an HSP) and the TOPN best references vote, weighted by chain bit
score; the position with the most weight wins. A boundary that no reference reaches is projected without
alignment from the chain end; when that lands outside the contig the gene is partial (< / >). ITS1 / ITS2 are
the spacers between the rRNAs.
Histones: tblastn of the reference proteins, HSPs grouped into loci, exonerate protein2genome and protein2dna
(intronless) with the TOPK best proteins of the locus, best score wins (spliced models x INTRON_PEN per
intron); the CDS is extended in frame to the next stop codon and, when the model does not start with ATG at
protein position 1, to the farthest upstream in-frame ATG before a stop. Nuclear code (transl_table 1).
Massive assemblies: both searches first run with one reference per genus (rDNA) / per gene and phylum
(histones) against a BLAST database of all contigs; only the regions hit (+-PAD nt) go to the steps above.

Per-gene consensus (multi-copy rDNA split over many contigs of one assembly): see consensus(); output
consensus/<gene>.fasta and consensus.tsv (which contigs were used, which were dropped and why).

Evidence (2026-10-23, benchmark/nuc_bench.py, knowledge/topics/nuclear_markers.md): NCBI metazoan records with
complete 18S-ITS1-5.8S-ITS2-28S annotation, one per species (378), references of the same genus left out:
boundaries exact 85.2 %, within 5 nt 92.6 % (barrnap 0.9 --kingdom euk: 1.8 % / 42.1 %); 677 complete histone
CDS: whole model (start, stop, exons, name) exact 97.5 %.

Usage:
  mitom nuc fetch -o nuc_refs                           # NCBI GenBank records: rdna.gb, histone.gb
  mitom nuc build --gb nuc_refs/rdna.gb nuc_refs/histone.gb [my_geneious.gb] -o models/nuc
  mitom nuc annotate -i contigs.fasta -o out [--genes rdna,histone] [--threads 16] [--lib models/nuc]
Outputs of annotate in -o: nuc_features.tsv (one row per feature), nuc_annotation.gb (contigs with features),
genes/<gene>.fasta (every copy; + histone .faa), consensus/<gene>.fasta + consensus.tsv (one sequence per
gene), nuc.log.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import logging
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqFeature import BeforePosition, AfterPosition
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402

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
EXTRAP = 30         # projection past an HSP end that still counts as a vote (nt)
MIN_CHAIN_BITS = 100
MIN_GENE_OVERLAP = 80   # nt of a reference rRNA an alignment chain must cover for the gene to be reported
MAX_ITS = 3000      # spacer between an rRNA and the contig end that is still reported as a partial ITS
MAX_RESCUE = 3000   # stretch beside an rRNA (or a whole sequence) searched again for a missing 5.8S (blastn word 11)
CURATED_MIN_ID = 90 # a curated reference decides boundaries only from chains of >= this % identity (close relatives)
QC_TOL = 15         # library QC: a reference boundary this far (nt) from the one other studies map onto it ...
QC_MIN_STUDIES = 3  # ... supported by >= this many other studies (and > 2x the studies that agree with it) is not used
TOPK = 3            # histone proteins tried with exonerate per locus
INTRON_PEN = 0.9
MIN_PROT_COVER = 0.3    # histone model must cover this fraction of the protein unless it reaches a contig end
SEED_COV = 0.1      # consensus seed: fragments with >= this fraction of the best fragment's contig coverage
PAD = 20000         # flank around prefilter hits taken from long contigs
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


def run(cmd: List[str]) -> str:
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


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
        run([tool("makeblastdb"), "-in", str(fa), "-dbtype", "nucl", "-out", str(Path(td) / "db")])
        tsv = run([tool("blastn"), "-task", "megablast", "-word_size", "16", "-query", str(fa),
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
                near = lambda p: set().union(set(), *(st for x, st in vv.items() if abs(x - p) <= 5))
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
        run([tool("makeblastdb"), "-in", str(self.rdna_fa), "-dbtype", "nucl", "-out", str(work / "rdna_ref")])
        self.rdna_db = str(work / "rdna_ref")
        self.curated_db = None
        if self.curated:
            cfa = work / "curated_ref.fa"
            SeqIO.write([r for r in SeqIO.parse(str(self.rdna_fa), "fasta") if r.id in self.curated], str(cfa), "fasta")
            run([tool("makeblastdb"), "-in", str(cfa), "-dbtype", "nucl", "-out", str(work / "curated_ref")])
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


# ── prefilter: candidate regions of the input ──────────────────────────────
def prefilter(contigs: Path, lib: Library, genes: List[str], work: Path,
              threads: int) -> Dict[str, List[Tuple[int, int]]]:
    """contig -> merged (start, end) windows hit by a representative reference (+- PAD, clipped)"""
    db = str(work / "contigs")
    run([tool("makeblastdb"), "-in", str(contigs), "-dbtype", "nucl", "-out", db])
    fmt = "6 sseqid sstart send slen"
    hits = []
    if "rdna" in genes:
        hits += run([tool("blastn"), "-task", "megablast", "-word_size", "20", "-query", str(lib.rdna_rep), "-db", db,
                     "-evalue", "1e-20", "-outfmt", fmt, "-max_target_seqs", "100000", "-num_threads", str(threads),
                     "-dust", "no"]).splitlines()
    if "histone" in genes:
        hits += run([tool("tblastn"), "-query", str(lib.hist_rep), "-db", db, "-evalue", "1e-10", "-outfmt", fmt,
                     "-max_target_seqs", "100000", "-num_threads", str(threads), "-seg", "no"]).splitlines()
    win = defaultdict(list)
    for line in hits:
        s, a, b, n = line.split("\t")
        a, b, n = sorted((int(a), int(b))) + [int(n)]
        win[s].append((max(1, a - PAD), min(n, b + PAD)))
    merged = {}
    for s, ws in win.items():
        ws.sort()
        m = [list(ws[0])]
        for a, b in ws[1:]:
            if a <= m[-1][1] + 1:
                m[-1][1] = max(m[-1][1], b)
            else:
                m.append([a, b])
        merged[s] = [tuple(x) for x in m]
    return merged


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


def project(h, b: int, maxd: Optional[int]) -> Optional[int]:
    """reference position b outside the HSP -> query position, ungapped from the nearer HSP end"""
    qs, qe, ss, se = h[:4]
    st = _strand(h)
    d1, d2 = (ss - b) * st, (b - se) * st
    if d1 > 0 and (maxd is None or d1 <= maxd):
        return qs - d1
    if d2 > 0 and (maxd is None or d2 <= maxd):
        return qe + d2
    return None


def rdna_loci(hits: Dict[str, list], lib: Library, qlen: int, min_overlap: int = MIN_GENE_OVERLAP,
              min_bits: float = MIN_CHAIN_BITS) -> List[dict]:
    """one locus per rDNA copy: chains of all references grouped by query overlap; boundaries by weighted vote"""
    allc = []
    for ref, hs in hits.items():
        for c in chains(hs, min_bits):
            allc.append((sum(h[4] for h in c), ref, c))
    allc.sort(key=lambda x: -x[0])
    loci: List[dict] = []
    for sc, ref, c in allc:
        lo, hi, st = c[0][0], c[-1][1], _strand(c[0])
        for L in loci:
            ov = min(hi, L["hi"]) - max(lo, L["lo"])
            if L["st"] == st and ov > 0.5 * min(hi - lo, L["hi"] - L["lo"]):
                if ref not in {r for _, r, _ in L["chains"]}:
                    L["chains"].append((sc, ref, c))
                break
        else:
            loci.append({"lo": lo, "hi": hi, "st": st, "chains": [(sc, ref, c)]})
    out = []
    for L in loci:
        top = L["chains"][:TOPN]
        st = L["st"]
        feats = {}
        for g in RRNA:
            votes = {"start": Counter(), "end": Counter()}
            proj = {"start": [], "end": []}
            covered = 0
            for sc, ref, c in top:
                if g in lib.rdna[ref]:
                    a, b = lib.rdna[ref][g][:2]
                    cov = set()
                    for h in c:
                        cov.update(range(max(a, min(h[2], h[3])), min(b, max(h[2], h[3])) + 1))
                    covered = max(covered, len(cov))
            trusted = {}
            for side, k in (("start", 0), ("end", 1)):
                usable = [x for x in L["chains"] if g in lib.rdna[x[1]] and lib.rdna[x[1]][g][2 + k]]
                # a curated reference aligned over the boundary (or <= EXTRAP nt past its alignment) decides it
                for sc, ref, c in (x for x in usable if x[1] in lib.curated and _identity(x[2]) >= CURATED_MIN_ID):
                    pos = lib.rdna[ref][g][k]
                    m = next((x for x in (map_pos(h, pos) for h in c) if x is not None), None)
                    if m is None:
                        m = next((x for x in (project(h, pos, EXTRAP) for h in c) if x is not None), None)
                    if m is not None:
                        trusted[side] = int(m)
                        break
                # a side no top reference may vote on (partial or failed QC) takes the next references that may
                for sc, ref, c in [x for x in usable if x in top] or usable[:TOPN]:
                    pos = lib.rdna[ref][g][k]
                    m = next((x for x in (map_pos(h, pos) for h in c) if x is not None), None)
                    if m is None:
                        m = next((x for x in (project(h, pos, EXTRAP) for h in c) if x is not None), None)
                    if m is not None:
                        votes[side][m] += sc
                    else:
                        near = min(c, key=lambda h: min(abs(h[2] - pos), abs(h[3] - pos)))
                        proj[side].append(project(near, pos, None))
            ends = {}
            for side in ("start", "end"):
                if side in trusted:
                    ends[side] = (trusted[side], "curated")
                elif votes[side]:
                    ends[side] = (int(votes[side].most_common(1)[0][0]), "aligned")
                elif [p for p in proj[side] if p is not None]:
                    p = sorted(x for x in proj[side] if x is not None)
                    ends[side] = (p[len(p) // 2], "projected")
                else:
                    ends[side] = None
            if ends["start"] is None or ends["end"] is None:
                continue
            g5, g3 = ends["start"], ends["end"]           # gene 5' / 3' in query coordinates
            lo_, hi_ = (g5, g3) if st == 1 else (g3, g5)
            p_lo, p_hi = lo_[0] < 1, hi_[0] > qlen
            lo_p, hi_p = max(1, lo_[0]), min(qlen, hi_[0])
            if hi_p < lo_p or (hi_p - lo_p < 20 and not (p_lo or p_hi)):      # short pieces only at a contig end
                continue
            feats[g] = {"lo": lo_p, "hi": hi_p, "p_lo": p_lo, "p_hi": p_hi, "weak": covered < min_overlap,
                        "how": f"{lo_[1]}/{hi_[1]}" if st == 1 else f"{hi_[1]}/{lo_[1]}"}
        # a gene needs MIN_GENE_OVERLAP aligned nt, except a piece cut by the contig end next to a supported gene
        if all(f["weak"] for f in feats.values()):
            continue
        feats = {g: f for g, f in feats.items() if not f["weak"] or f["p_lo"] or f["p_hi"]}
        _spacers(feats, st, qlen)
        best = top[0][1]
        ident = _identity(top[0][2])
        out.append({"st": st, "feats": feats, "ref": best, "ref_org": lib.org[best][0], "identity": ident})
    return out


def _spacers(feats: dict, st: int, qlen: int) -> None:
    """ITS1 / ITS2 between the left and right rRNA in forward coordinates; a missing rRNA is allowed when the
    contig ends within MAX_ITS nt (partial spacer), not when 18S and 28S are there and only 5.8S was not found"""
    for its in ("ITS1", "ITS2"):
        feats.pop(its, None)
    if "5.8S" not in feats and "18S" in feats and "28S" in feats:
        return
    for its, (g1, g2) in (("ITS1", ("18S", "5.8S")), ("ITS2", ("5.8S", "28S"))):
        L, R = (feats.get(g1), feats.get(g2)) if st == 1 else (feats.get(g2), feats.get(g1))
        if (L and L["p_hi"]) or (R and R["p_lo"]) or not (L or R):
            continue
        a = L["hi"] + 1 if L else (1 if R["lo"] - 1 <= MAX_ITS else None)
        b = R["lo"] - 1 if R else (qlen if qlen - L["hi"] <= MAX_ITS else None)
        if a is not None and b is not None and b - a >= 9:
            feats[its] = {"lo": a, "hi": b, "p_lo": L is None, "p_hi": R is None, "how": "between rRNAs"}


def _rescue_58s(loci: Dict[str, List[dict]], regions: Path, lib: Library, work: Path, leave_out: Dict[str, str],
                threads: int, qlens: Dict[str, int]) -> None:
    """5.8S (~155 nt) of divergent taxa is often below the megablast seed. Searched again with blastn word 11 and
    voted on like a whole unit: the spacer where 18S and 28S were found but no 5.8S, and the whole sequence
    (<= MAX_RESCUE nt) where no 5.8S was found at all (5.8S RNA records, ITS amplicons, short contigs)"""
    spacers = {}
    for q, ls in loci.items():
        if any("5.8S" in L["feats"] for L in ls):
            continue
        n = qlens[q]
        for i, L in enumerate(ls):
            f = L["feats"]
            if "18S" in f and "28S" in f:
                a, b = (f["18S"]["hi"] + 1, f["28S"]["lo"] - 1) if L["st"] == 1 else \
                    (f["28S"]["hi"] + 1, f["18S"]["lo"] - 1)
                if b - a >= 100:
                    spacers[f"{q}|{i}|{a}"] = (q, a, b)
                continue
            # one rRNA: the stretch on the spacer side of it up to the contig end (ITS amplicons, contigs ending in ITS)
            for g, side in (("18S", 1), ("28S", -1)):
                if g in f:
                    a, b = (f[g]["hi"] + 1, min(n, f[g]["hi"] + MAX_RESCUE)) if side * L["st"] == 1 else \
                        (max(1, f[g]["lo"] - MAX_RESCUE), f[g]["lo"] - 1)
                    if b - a >= 100:
                        spacers[f"{q}|{i}|{a}"] = (q, a, b)
        if not ls and n <= MAX_RESCUE:      # rDNA hits but no locus: 5.8S RNA records, short ITS contigs
            spacers[f"{q}|-1|1"] = (q, 1, n)
    if not spacers:
        return
    seqs = SeqIO.index(str(regions), "fasta")
    fa = work / "rdna_58s.fa"
    with open(fa, "w") as o:
        for k, (q, a, b) in spacers.items():
            o.write(f">{k}\n{seqs[q].seq[a - 1:b]}\n")
    hits: Dict[str, Dict[str, list]] = defaultdict(lambda: defaultdict(list))
    for line in run([tool("blastn"), "-task", "blastn", "-word_size", "11", "-query", str(fa), "-db", lib.rdna_db,
                     "-evalue", "1e-5", "-dust", "no", "-num_threads", str(max(1, threads)),
                     "-max_target_seqs", "500",
                     "-outfmt", "6 qseqid sseqid qstart qend sstart send bitscore qseq sseq"]).splitlines():
        k, s, qs, qe, ss, se, bs, qa, sa = line.split("\t")
        q = spacers[k][0]
        if leave_out and leave_out.get(_contig(q)) == lib.org[s][1] or _contig(q) == s:
            continue
        hits[k][s].append((int(qs), int(qe), int(ss), int(se), float(bs), qa, sa))
    for k, h in hits.items():
        q, a, b = spacers[k]
        i = int(k.split("|")[1])
        for x in rdna_loci(h, lib, b - a + 1, min_overlap=40, min_bits=50):
            f = x["feats"].get("5.8S")
            if not f or (i >= 0 and x["st"] != loci[q][i]["st"]):
                continue
            # a piece cut by the searched stretch is only kept at the contig end
            if (f["p_lo"] and a > 1) or (f["p_hi"] and b < qlens[q]):
                continue
            f = dict(f, lo=f["lo"] + a - 1, hi=f["hi"] + a - 1, p_lo=f["p_lo"] and a == 1, p_hi=f["p_hi"] and b == qlens[q])
            if i < 0:
                L = dict(x, feats={})
                loci[q].append(L)
            else:
                L = loci[q][i]
            L["feats"]["5.8S"] = f
            _spacers(L["feats"], L["st"], qlens[q])
            break


def _identity(chain) -> float:
    m = n = 0
    for h in chain:
        m += sum(1 for x, y in zip(h[5], h[6]) if x == y)
        n += len(h[5])
    return round(100 * m / n, 1) if n else 0.0


def annotate_rdna(regions: Path, lib: Library, work: Path, threads: int, leave_out: Dict[str, str],
                  qlens: Dict[str, int]) -> Dict[str, List[dict]]:
    fmt = "6 qseqid sseqid qstart qend sstart send bitscore qseq sseq"
    # BLAST threads scale poorly over many short queries: one single-threaded blastn per chunk of the regions
    chunks: List[list] = [[] for _ in range(max(1, threads))]
    size = [0] * len(chunks)
    for r in sorted(SeqIO.parse(str(regions), "fasta"), key=lambda r: -len(r)):
        i = size.index(min(size))
        chunks[i].append(r)
        size[i] += len(r)

    def blast(i: int) -> Path:
        fa, tsv = work / f"rdna_q{i}.fa", work / f"rdna_blast{i}.tsv"
        SeqIO.write(chunks[i], str(fa), "fasta")
        with open(tsv, "w") as o:
            subprocess.run([tool("blastn"), "-task", "megablast", "-word_size", "16", "-query", str(fa),
                            "-db", lib.rdna_db,
                            "-evalue", "1e-10", "-outfmt", fmt, "-num_threads", "1",
                            "-max_target_seqs", "500",   # short curated references and left-out genera fall out of a top 50
                            "-dust", "no"], stdout=o, check=True)
            # curated references again with blastn word 11: their conserved ends must align in distant taxa too
            if lib.curated_db:
                subprocess.run([tool("blastn"), "-task", "blastn", "-word_size", "11", "-query", str(fa),
                                "-db", lib.curated_db, "-evalue", "1e-5", "-outfmt", fmt, "-num_threads", "1",
                                "-dust", "no"], stdout=o, check=True)
        return tsv

    with ThreadPoolExecutor(len(chunks)) as pool:
        tsvs = list(pool.map(blast, [i for i in range(len(chunks)) if chunks[i]]))
    hits: Dict[str, Dict[str, list]] = defaultdict(lambda: defaultdict(list))
    for line in (x for t in tsvs for x in open(t)):
        q, s, qs, qe, ss, se, bs, qa, sa = line.rstrip("\n").split("\t")
        if leave_out and leave_out.get(_contig(q)) == lib.org[s][1] or _contig(q) == s:
            continue
        h = (int(qs), int(qe), int(ss), int(se), float(bs), qa, sa)
        if h not in hits[q][s]:      # a curated reference found by both searches
            hits[q][s].append(h)
    loci = {q: rdna_loci(h, lib, qlens[q]) for q, h in hits.items()}
    _rescue_58s(loci, regions, lib, work, leave_out, threads, qlens)
    return loci


# ── histones ────────────────────────────────────────────────────────────────
def _exonerate(target: str, prot: str, model: str, td: Path) -> Optional[tuple]:
    (td / "p.fa").write_text(f">p\n{prot}\n")
    out = subprocess.run([tool("exonerate"), "--model", model, "--showalignment", "no", "--showvulgar", "no",
                          "--showtargetgff", "yes", "--bestn", "1", "--geneticcode", "1", "-q", str(td / "p.fa"),
                          "-t", target, "--ryo", "SCORE %s %qab %qae %ql\n"], capture_output=True, text=True).stdout
    sc, ex = None, []
    for line in out.splitlines():
        if line.startswith("SCORE"):
            _, sc, qab, qae, ql = line.split()
            sc, qab, qae, ql = int(sc), int(qab), int(qae), int(ql)
        x = line.split("\t")
        # protein2genome lists exons; protein2dna (no introns) only one similarity line
        if len(x) > 8 and (x[2] == "exon" or x[2] == "similarity" and model == "protein2dna"):
            ex.append((int(x[3]), int(x[4]), x[6]))
    return (sc, ex, qab, qae, ql) if sc and ex else None


def _orf(seq: Seq, ex: List[Tuple[int, int]], st: int, at_start: bool) -> Tuple[List[Tuple[int, int]], bool, bool]:
    """extend the exon model (forward coordinates) to a stop codon and an ATG; returns exons, 5' and 3' partial"""
    n = len(seq)
    cs = seq if st == 1 else seq.reverse_complement()
    tr = (lambda a, b: (a, b)) if st == 1 else (lambda a, b: (n - b + 1, n - a + 1))
    E = sorted(tr(a, b) for a, b in ex)
    p, p3 = E[-1][1], True
    for _ in range(300):
        c = str(cs[p:p + 3]).upper()
        if len(c) < 3:
            break
        p += 3
        if c in STOPS:
            p3 = False
            break
    E[-1] = (E[-1][0], min(p, n))
    a0, p5 = E[0][0], False
    if not (at_start and str(cs[a0 - 1:a0 + 2]).upper() == "ATG"):
        q, cand = a0, []
        while q - 3 >= 1:
            c = str(cs[q - 4:q - 1]).upper()
            if c in STOPS:
                break
            q -= 3
            if c == "ATG":
                cand.append(q)
        if cand:
            E[0] = (cand[-1], E[0][1])
        else:
            p5 = str(cs[a0 - 1:a0 + 2]).upper() != "ATG"
            if p5 and q - 3 < 1:
                E[0] = (q, E[0][1])                   # open to the contig end
    return sorted(tr(a, b) for a, b in E), (p5 if st == 1 else p3), (p3 if st == 1 else p5)


def annotate_histones(regions: Path, lib: Library, work: Path, threads: int, leave_out: Dict[str, str],
                      seqs: Dict[str, Seq]) -> Dict[str, List[dict]]:
    db = str(work / "regions")
    run([tool("makeblastdb"), "-in", str(regions), "-dbtype", "nucl", "-out", db])
    out = run([tool("tblastn"), "-query", str(lib.hist_faa), "-db", db, "-evalue", "1e-5", "-outfmt",
               "6 qseqid sseqid sstart send bitscore", "-max_target_seqs", "100000", "-num_threads", str(threads),
               "-seg", "no"])
    hits = defaultdict(list)
    for line in out.splitlines():
        ref, s, a, b, bs = line.split("\t")
        r = lib.hist[ref]
        if leave_out and leave_out.get(_contig(s)) == r["genus"] or _contig(s) == r["accession"]:
            continue
        a, b = int(a), int(b)
        hits[s].append({"ref": ref, "lo": min(a, b), "hi": max(a, b), "st": 1 if b >= a else -1, "bs": float(bs)})

    def job(s: str) -> List[dict]:
        seq, res = seqs[s], []
        loci: List[dict] = []
        for h in sorted(hits[s], key=lambda h: (h["st"], h["lo"])):
            if loci and loci[-1]["st"] == h["st"] and h["lo"] <= loci[-1]["hi"]:
                loci[-1]["hi"] = max(loci[-1]["hi"], h["hi"])
                loci[-1]["h"].append(h)
            else:
                loci.append({"st": h["st"], "lo": h["lo"], "hi": h["hi"], "h": [h]})
        with tempfile.TemporaryDirectory(dir=work) as td:
            td = Path(td)
            for L in loci:
                a, b = max(0, L["lo"] - 1 - 3000), min(len(seq), L["hi"] + 3000)
                (td / "t.fa").write_text(f">t\n{seq[a:b]}\n")
                refs = []
                for h in sorted(L["h"], key=lambda h: -h["bs"]):
                    if h["ref"] not in refs:
                        refs.append(h["ref"])
                models = []
                for ref in refs[:TOPK]:
                    for model in ("protein2genome", "protein2dna"):
                        m = _exonerate(str(td / "t.fa"), lib.prot[ref], model, td)
                        if m:
                            sc = m[0] * (INTRON_PEN ** (len(m[1]) - 1) if model == "protein2genome" else 1)
                            models.append((sc, ref, m))
                if not models:
                    continue
                sc, ref, (_, ex, qab, qae, ql) = max(models, key=lambda x: x[0])
                st = 1 if ex[0][2] == "+" else -1
                ex = [(x + a, y + a) for x, y, _ in ex]
                exons, p5, p3 = _orf(seq, ex, st, qab == 0)
                at_edge = exons[0][0] <= 3 or exons[-1][1] >= len(seq) - 2
                if (qae - qab) / ql < MIN_PROT_COVER and not at_edge:
                    continue
                res.append({"gene": lib.hist[ref]["gene"], "st": st, "exons": exons, "p5": p5, "p3": p3,
                            "ref": lib.hist[ref]["accession"], "ref_org": lib.hist[ref]["organism"], "score": sc})
        # same CDS found from two overlapping loci: keep the better one
        res.sort(key=lambda r: -r["score"])
        kept = []
        for r in res:
            if not any(k["st"] == r["st"] and k["exons"][0][0] <= r["exons"][-1][1] and
                       r["exons"][0][0] <= k["exons"][-1][1] for k in kept):
                kept.append(r)
        return kept

    with ThreadPoolExecutor(max(1, threads)) as pool:
        return dict(zip(list(hits), pool.map(job, list(hits))))


# ── input, regions, output ─────────────────────────────────────────────────
def _contig(region_id: str) -> str:
    return region_id.rsplit(":", 1)[0] if re.search(r":\d+-\d+$", region_id) else region_id


def _offset(region_id: str) -> int:
    m = re.search(r":(\d+)-\d+$", region_id)
    return int(m.group(1)) - 1 if m else 0


def read_input(path: Path, work: Path) -> Tuple[Path, Dict[str, str], int, int]:
    """FASTA or GenBank -> FASTA with unique BLAST-safe ids; returns it, organism by id (GenBank input),
    number of sequences and of nucleotides. Records are read back on demand (SeqIO.index), so an assembly of
    millions of contigs is never held in memory."""
    fmt = "genbank" if path.suffix.lower() in (".gb", ".gbk", ".gbf", ".genbank") else "fasta"
    ids, org, nt = set(), {}, 0
    fa = work / "input.fa"
    with open(fa, "w") as o:
        for r in SeqIO.parse(str(path), fmt):
            rid = re.sub(r"[^A-Za-z0-9_.-]", "_", r.id)
            if rid in ids:
                raise SystemExit(f"duplicate sequence id {rid} in {path}")
            ids.add(rid)
            nt += len(r.seq)
            if "organism" in r.annotations:
                org[rid] = r.annotations["organism"]
            desc = r.description[len(r.id):].strip() if r.description.startswith(r.id) else r.description
            o.write(f">{rid} {desc}\n{r.seq}\n" if desc else f">{rid}\n{r.seq}\n")
    return fa, org, len(ids), nt


def write_regions(windows: Dict[str, List[Tuple[int, int]]], recs, fa: Path) -> Dict[str, Seq]:
    """whole contig when one window covers it, else contig:start-end pieces"""
    seqs = {}
    with open(fa, "w") as o:
        for c, ws in windows.items():
            s = recs[c].seq
            for a, b in ws:
                rid = c if (a == 1 and b == len(s)) else f"{c}:{a}-{b}"
                seqs[rid] = s[a - 1:b]
                o.write(f">{rid}\n{seqs[rid]}\n")
    return seqs


def _cov(name: str) -> str:
    m = re.search(r"cov_([\d.]+)", name)
    return m.group(1) if m else ""


def _loc(a: int, b: int, st: int, p_lo: bool, p_hi: bool) -> FeatureLocation:
    return FeatureLocation(BeforePosition(a - 1) if p_lo else a - 1, AfterPosition(b) if p_hi else b, strand=st)


def collect(recs, rdna, hist) -> Tuple[List[dict], Dict[str, List[SeqFeature]]]:
    rows, feats = [], defaultdict(list)
    for rid, loci in (rdna or {}).items():
        c, off = _contig(rid), _offset(rid)
        for L in loci:
            for g in UNIT:
                if g not in L["feats"]:
                    continue
                f = L["feats"][g]
                a, b = f["lo"] + off, f["hi"] + off
                n = len(recs[c].seq)
                p_lo, p_hi = f["p_lo"], f["p_hi"]
                p5, p3 = (p_lo, p_hi) if L["st"] == 1 else (p_hi, p_lo)
                ftype = "rRNA" if g in RRNA else "misc_RNA"
                note = f"boundaries {f['how']}; best reference {L['ref']} {L['ref_org']} ({L['identity']} % identity)"
                feats[c].append(SeqFeature(_loc(a, b, L["st"], p_lo, p_hi), type=ftype,
                                           qualifiers={"product": [PRODUCT[g]], "note": [note]}))
                rows.append({"contig": c, "contig_len": n, "contig_cov": _cov(c), "type": ftype, "gene": g,
                             "start": a, "end": b, "strand": "+" if L["st"] == 1 else "-",
                             "partial_5": int(p5), "partial_3": int(p3), "length": b - a + 1, "exons": f"{a}-{b}",
                             "best_ref": L["ref"], "best_ref_organism": L["ref_org"], "identity": L["identity"],
                             "evidence": f["how"]})
    for rid, models in (hist or {}).items():
        c, off = _contig(rid), _offset(rid)
        for m in models:
            ex = [(x + off, y + off) for x, y in m["exons"]]
            lo_p, hi_p = (m["p5"], m["p3"]) if m["st"] == 1 else (m["p3"], m["p5"])
            parts = [_loc(x, y, m["st"], lo_p and i == 0, hi_p and i == len(ex) - 1) for i, (x, y) in enumerate(ex)]
            if m["st"] == -1:
                parts = parts[::-1]
            loc = parts[0] if len(parts) == 1 else CompoundLocation(parts)
            seq = loc.extract(recs[c].seq)
            prot = str(seq[: len(seq) // 3 * 3].translate(table=1)).rstrip("*")
            q = {"gene": [m["gene"]], "product": [PRODUCT[m["gene"]]], "codon_start": ["1"], "transl_table": ["1"],
                 "note": [f"exonerate model from {m['ref']} {m['ref_org']}"]}
            if "*" not in prot:
                q["translation"] = [prot]
            feats[c].append(SeqFeature(loc, type="CDS", qualifiers=q))
            rows.append({"contig": c, "contig_len": len(recs[c].seq), "contig_cov": _cov(c), "type": "CDS",
                         "gene": m["gene"], "start": ex[0][0], "end": ex[-1][1], "strand": "+" if m["st"] == 1 else "-",
                         "partial_5": int(m["p5"]), "partial_3": int(m["p3"]), "length": sum(y - x + 1 for x, y in ex),
                         "exons": ";".join(f"{x}-{y}" for x, y in ex), "best_ref": m["ref"],
                         "best_ref_organism": m["ref_org"], "identity": "",
                         "evidence": "internal stop codon" if "*" in prot else "exonerate"})
    return rows, feats


def write_outputs(out: Path, recs, org: Dict[str, str], rows: List[dict], feats: Dict[str, List[SeqFeature]]) -> None:
    rows.sort(key=lambda r: (r["contig"], r["start"]))
    cols = ["contig", "contig_len", "contig_cov", "type", "gene", "start", "end", "strand", "partial_5", "partial_3",
            "length", "exons", "best_ref", "best_ref_organism", "identity", "evidence"]
    with open(out / "nuc_features.tsv", "w", newline="") as o:
        w = csv.DictWriter(o, cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    gbrecs = []
    for c, fs in feats.items():
        r = recs[c]
        g = SeqRecord(r.seq, id=c[:16] if len(c) > 16 else c, name=c[:16], description=r.description,
                      annotations={"molecule_type": "DNA", "organism": org.get(c, "unknown"),
                                   "topology": "linear"})
        g.features = sorted(fs, key=lambda f: int(f.location.start))
        gbrecs.append(g)
    SeqIO.write(gbrecs, str(out / "nuc_annotation.gb"), "genbank")
    gd = out / "genes"
    gd.mkdir(exist_ok=True)
    by_gene = defaultdict(list)
    for c, fs in feats.items():
        for f in fs:
            g = f.qualifiers.get("gene", [None])[0] or \
                next(k for k, v in PRODUCT.items() if v == f.qualifiers["product"][0])
            s = f.location.extract(recs[c].seq)
            tag = f"{c}:{int(f.location.start) + 1}-{int(f.location.end)}({'+' if f.location.strand == 1 else '-'})"
            by_gene[g].append(SeqRecord(s, id=tag, description=""))
            if "translation" in f.qualifiers:
                by_gene[g + ".faa"].append(SeqRecord(Seq(f.qualifiers["translation"][0]), id=tag, description=""))
    for g, rs in by_gene.items():
        SeqIO.write(rs, str(gd / (g if g.endswith(".faa") else f"{g}.fasta")), "fasta")


def _gene_seq(r: dict, recs) -> Seq:
    s = recs[r["contig"]].seq
    seq = sum((s[int(a) - 1:int(b)] for a, b in (x.split("-") for x in r["exons"].split(";"))), Seq(""))
    return seq if r["strand"] == "+" else seq.reverse_complement()


def _ref_gene(lib: "Library", ref: str, g: str) -> Optional[str]:
    if g not in lib.rdna.get(ref, {}):
        return None
    if not hasattr(lib, "_seqs"):
        lib._seqs = SeqIO.index(str(lib.rdna_fa), "fasta")
    a, b = lib.rdna[ref][g][:2]
    return str(lib._seqs[ref].seq[a - 1:b])


def _pileup(frags: List[str], wt: List[float], backbone: str, work: Path) -> Tuple[str, Dict[int, Tuple[float, int]]]:
    """consensus of the fragments placed on the backbone by their best blastn HSP (plus strand): per backbone
    position the weighted majority base or deletion; an insertion is kept when most of the spanning weight has
    it. Uppercase backbone letters are sample sequence, lowercase ones only a reference placeholder: positions no
    fragment spans keep the backbone letter in lowercase. Returns the consensus and, per fragment, (identity,
    aligned columns) over the uppercase backbone positions."""
    (work / "bb.fa").write_text(f">bb\n{backbone}\n")
    with open(work / "fr.fa", "w") as o:
        for i, f in enumerate(frags):
            o.write(f">{i}\n{f}\n")
    out = run([tool("blastn"), "-task", "blastn", "-query", str(work / "fr.fa"), "-subject", str(work / "bb.fa"),
               "-strand", "plus", "-dust", "no", "-evalue", "1e-5",
               "-outfmt", "6 qseqid qstart qend sstart send bitscore qseq sseq"])
    n = len(backbone)
    best: Dict[int, tuple] = {}
    for line in out.splitlines():
        q, qs, qe, ss, se, bs, qa, sa = line.split("\t")
        q, qs, qe, ss, se = int(q), int(qs), int(qe), int(ss), int(se)
        if q in best and float(bs) <= best[q][2]:
            continue
        # local alignments stop before terminal mismatches: extend ungapped to the end of the overlap (<= 20 nt)
        f = frags[q]
        e5, e3 = min(qs - 1, ss - 1, 20), min(len(f) - qe, n - se, 20)
        qa = f[qs - 1 - e5:qs - 1] + qa + f[qe:qe + e3]
        sa = backbone[ss - 1 - e5:ss - 1] + sa + backbone[se:se + e3]
        best[q] = (ss - e5, se + e3, float(bs), qa.upper(), sa)
    base = [Counter() for _ in range(n + 1)]
    ins = [Counter() for _ in range(n + 1)]       # insertion after backbone position k
    span = [0.0] * (n + 2)
    ident: Dict[int, Tuple[float, int]] = {}
    for q, (ss, se, _, qa, sa) in best.items():
        w, k, buf, m, c = wt[q], ss - 1, "", 0, 0
        for x, y in zip(qa, sa):
            if y == "-":
                buf += x
                continue
            if buf:
                ins[k][buf] += w
            buf = ""
            k += 1
            base[k][x] += w
            if y.isupper():
                c += 1
                m += x == y
        for k in range(ss, se + 1):
            span[k] += w
        ident[q] = (m / c if c else 0.0, c)
    cons = []
    for k in range(1, n + 1):
        if span[k] == 0:
            cons.append(backbone[k - 1].lower())
            continue
        x = base[k].most_common(1)[0][0]
        cons.append("" if x == "-" else x)
        if k < n and ins[k]:
            s_, w_ = ins[k].most_common(1)[0]
            if w_ > min(span[k], span[k + 1]) / 2:
                cons.append(s_)
    return "".join(cons), ident


def _assemble(frags: List[str], wt: List[float], ref: Optional[str], work: Path) -> str:
    """backbone from fragments of one genome: contained fragments are removed and dovetailing ones merged
    (blastn, >= 98 % identity, >= 40 nt overlap) until nothing changes; several pieces left are ordered on the
    reference gene; neighbours sharing an exact end overlap of >= 15 nt are joined on it, others with lowercase n
    for the reference distance between them (a piece placed over a heavier one is dropped). Without a reference only the heaviest piece is kept."""
    pieces = sorted(zip(frags, wt), key=lambda x: -x[1])
    changed = True
    while changed and len(pieces) > 1:
        changed = False
        with open(work / "pc.fa", "w") as o:
            for i, (f, _) in enumerate(pieces):
                o.write(f">{i}\n{f}\n")
        out = run([tool("blastn"), "-task", "blastn", "-query", str(work / "pc.fa"), "-subject", str(work / "pc.fa"),
                   "-strand", "plus", "-dust", "no", "-perc_identity", "98", "-evalue", "1e-10",
                   "-outfmt", "6 qseqid sseqid qstart qend sstart send length"])
        for line in out.splitlines():
            q, t, qs, qe, ss, se, ln = map(int, line.split("\t"))
            if q == t or ln < 40:
                continue
            lq, lt = len(pieces[q][0]), len(pieces[t][0])
            if qs <= 3 and qe >= lq - 2 and lq <= lt:                  # q inside t
                pieces[t] = (pieces[t][0], pieces[t][1] + pieces[q][1])
                del pieces[q]
                changed = True
                break
            if qe >= lq - 2 and ss <= 3 and se < lt - 2:                 # q's end overlaps t's start
                pieces[q] = (pieces[q][0][:qe] + pieces[t][0][se:], pieces[q][1] + pieces[t][1])
                del pieces[t]
                changed = True
                break
    if len(pieces) == 1 or not ref:
        return pieces[0][0]
    with open(work / "pc.fa", "w") as o:
        for i, (f, _) in enumerate(pieces):
            o.write(f">{i}\n{f}\n")
    (work / "rf.fa").write_text(f">r\n{ref}\n")
    out = run([tool("blastn"), "-task", "blastn", "-query", str(work / "pc.fa"), "-subject", str(work / "rf.fa"),
               "-strand", "plus", "-dust", "no", "-evalue", "1e-5",
               "-outfmt", "6 qseqid qstart qend sstart send bitscore"])
    place: Dict[int, Tuple[int, int, float]] = {}
    for line in out.splitlines():
        q, qs, qe, ss, se, bs = line.split("\t")
        q, qs, qe, ss, se = int(q), int(qs), int(qe), int(ss), int(se)
        if q not in place or float(bs) > place[q][2]:
            place[q] = (ss - (qs - 1), se + (len(pieces[q][0]) - qe), float(bs))   # piece ends in reference positions
    kept: List[int] = []
    for q in sorted(place, key=lambda q: -pieces[q][1]):
        a, b, _ = place[q]
        if all(min(b, place[k][1]) - max(a, place[k][0]) < 50 for k in kept):
            kept.append(q)
    kept.sort(key=lambda q: place[q][0])
    seq = pieces[kept[0]][0]
    for prev, q in zip(kept, kept[1:]):
        nxt = pieces[q][0]
        # assembler contigs often share an exact end overlap shorter than the blastn one (k-mer overlap)
        k = next((k for k in range(min(len(seq), len(nxt), 300), 14, -1) if seq[-k:] == nxt[:k]), 0)
        seq += nxt[k:] if k else "n" * max(1, place[q][0] - place[prev][1] - 1) + nxt
    return seq


def consensus(rows: List[dict], recs, lib: "Library", work: Path, min_ident: float) -> List[dict]:
    """One sequence per gene from all contigs. rDNA genes / spacers: the fragments with >= SEED_COV of the best
    contig coverage (contig cov_ of SPAdes names, else 1 = all fragments) are joined into a backbone by their
    overlaps and ordered on the reference gene of the best fragment (_assemble); every fragment is placed on it
    and kept when it has >= min_ident over >= 50 nt of backbone sequence; the others (other organisms,
    contamination) are listed as dropped. Each position takes the coverage-weighted majority of the kept
    fragments (_pileup); stretches no kept fragment covers are N (ends stripped). Histones: the CDS of the
    best-covered contig, complete and >= 80 % of the median reference length first (copies are paralogs; a
    consensus would be a chimera)."""
    out = []
    for g in UNIT:
        frags = [r for r in rows if r["gene"] == g]
        if not frags:
            continue
        seqs = [str(_gene_seq(r, recs)).upper() for r in frags]
        wt = [float(r["contig_cov"] or 1) for r in frags]
        top = max(range(len(frags)), key=lambda i: (wt[i], int(frags[i]["length"])))
        ref = _ref_gene(lib, frags[top]["best_ref"], g) if g in RRNA else None
        # seed: the well-covered fragments (same genome as the best one in an assembly of mixed samples)
        seed = [i for i in range(len(frags)) if wt[i] >= SEED_COV * wt[top]]
        backbone = _assemble([seqs[i] for i in seed], [wt[i] for i in seed], ref, work)
        _, idt = _pileup(seqs, wt, backbone, work)
        ident = {i: idt.get(i, (0.0, 0)) for i in range(len(frags))}
        keep = sorted({i for i in range(len(frags)) if ident[i][0] >= min_ident and ident[i][1] >= 50} | {top})
        cons, _ = _pileup([seqs[i] for i in keep], [wt[i] for i in keep], backbone, work)
        seq = re.sub("[a-z]", "N", cons).strip("N")
        out.append({"gene": g, "length": len(seq), "n_fragments": len(keep), "n_dropped": len(frags) - len(keep),
                    "internal_N": seq.count("N"),
                    "partial_5": int(all(frags[i]["partial_5"] == "1" for i in keep)),
                    "partial_3": int(all(frags[i]["partial_3"] == "1" for i in keep)),
                    "contigs": ";".join(frags[i]["contig"] for i in keep),
                    "dropped": ";".join(f"{frags[i]['contig']}({ident[i][0] * 100:.1f}% over {ident[i][1]} nt)"
                                        for i in range(len(frags))
                                        if i not in keep),
                    "seq": seq})
    for h in HISTONES:
        cds = [r for r in rows if r["gene"] == h]
        if not cds:
            continue
        lens = sorted(int(r["length_aa"]) for r in lib.hist.values() if r["gene"] == h)
        typical = 3 * lens[len(lens) // 2] if lens else 0
        r = max(cds, key=lambda r: (r["partial_5"] == r["partial_3"] == "0" and int(r["length"]) >= 0.8 * typical,
                                    float(r["contig_cov"] or 0), int(r["length"])))
        out.append({"gene": h, "length": int(r["length"]), "n_fragments": 1, "n_dropped": 0, "internal_N": 0,
                    "partial_5": int(r["partial_5"]), "partial_3": int(r["partial_3"]), "contigs": r["contig"],
                    "dropped": f"{len(cds) - 1} other copies", "seq": str(_gene_seq(r, recs))})
    return out


def write_consensus(out: Path, cons: List[dict]) -> None:
    d = out / "consensus"
    d.mkdir(exist_ok=True)
    cols = ["gene", "length", "n_fragments", "n_dropped", "internal_N", "partial_5", "partial_3", "contigs", "dropped"]
    with open(out / "consensus.tsv", "w", newline="") as o:
        w = csv.DictWriter(o, cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(cons)
    for c in cons:
        part = "".join(x for x, y in (("5'-partial ", c["partial_5"]), ("3'-partial ", c["partial_3"])) if y)
        (d / f"{c['gene']}.fasta").write_text(
            f">{c['gene']} {part}fragments={c['n_fragments']} length={c['length']}\n{c['seq']}\n")
        if c["gene"] in HISTONES:
            (d / f"{c['gene']}.faa").write_text(f">{c['gene']}\n{str(Seq(c['seq']).translate(table=1)).rstrip('*')}\n")


def annotate(inp: Path, out: Path, lib_dir: Path, genes: List[str], threads: int,
             leave_out: Optional[Path] = None, consensus_identity: float = 0.97) -> List[dict]:
    out.mkdir(parents=True, exist_ok=True)
    work = out / "tmp"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir()
    lo = {}
    if leave_out:
        for line in open(leave_out):
            if line.strip():
                k, v = line.rstrip("\n").split("\t")[:2]
                lo[k] = v
    lib = Library(lib_dir, work)
    fa, org, n_seq, n_nt = read_input(inp, work)
    recs = SeqIO.index(str(fa), "fasta")
    logging.info(f"[NUC] {n_seq} sequences, {n_nt:,} nt")
    windows = prefilter(fa, lib, genes, work, threads)
    regions = work / "regions.fa"
    seqs = write_regions(windows, recs, regions)
    logging.info(f"[NUC] prefilter: {len(windows)} contig(s), {len(seqs)} region(s), "
                 f"{sum(len(s) for s in seqs.values()):,} nt")
    rdna = hist = None
    if seqs and "rdna" in genes:
        rdna = annotate_rdna(regions, lib, work, threads, lo, {k: len(v) for k, v in seqs.items()})
    if seqs and "histone" in genes:
        hist = annotate_histones(regions, lib, work, threads, lo, seqs)
    rows, feats = collect(recs, rdna, hist)
    write_outputs(out, recs, org, rows, feats)
    if rows and consensus_identity > 0:
        rows_s = [{k: str(v) for k, v in r.items()} for r in rows]
        cons = consensus(rows_s, recs, lib, work, consensus_identity)
        write_consensus(out, cons)
        logging.info(f"[NUC] consensus: " + ", ".join(f"{c['gene']} {c['length']} nt from {c['n_fragments']}"
                                                     for c in cons) + f"  ({out / 'consensus.tsv'})")
    recs.close()
    shutil.rmtree(work)
    counts = Counter(r["gene"] for r in rows)
    found = ", ".join(f"{g} {counts[g]}" for g in (*UNIT, *HISTONES) if counts[g])
    logging.info(f"[NUC] features: {found or 'none'}")
    logging.info(f"[NUC] {out / 'nuc_features.tsv'}  {out / 'nuc_annotation.gb'}  {out / 'genes'}/")
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch", help="download the NCBI reference records (rdna.gb, histone.gb)")
    f.add_argument("-o", "--out", type=Path, default=Path("nuc_refs"))
    f.add_argument("--email")
    f.add_argument("--api_key")
    b = sub.add_parser("build", help="build the reference library from annotated GenBank files")
    b.add_argument("--gb", nargs="+", type=Path, required=True)
    b.add_argument("-o", "--out", type=Path, default=LIB)
    b.add_argument("--threads", type=int, default=8)
    b.add_argument("--curated", nargs="*", type=Path, default=[],
                   help="GenBank files of trusted references (curated.gb of mitom nuc fetch, own checked records)")
    a = sub.add_parser("annotate", help="annotate rDNA and histone genes in contigs (FASTA or GenBank)")
    a.add_argument("-i", "--input", type=Path, required=True)
    a.add_argument("-o", "--out", type=Path, required=True)
    a.add_argument("--lib", type=Path, default=LIB, help="reference library (mitom nuc build)")
    a.add_argument("--genes", default="rdna,histone", help="rdna, histone or both (default)")
    a.add_argument("--threads", type=int, default=8)
    a.add_argument("--consensus_identity", type=float, default=0.97,
                   help="fragments below this identity to a gene's consensus are left out of it (default 0.97; "
                        "0 = no consensus)")
    a.add_argument("--leave_out", type=Path, help=argparse.SUPPRESS)   # benchmark: TSV sequence id -> genus to skip
    args = ap.parse_args(argv)
    handlers = [logging.StreamHandler(sys.stdout)]
    if args.cmd == "annotate":
        args.out.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(args.out / "nuc.log", mode="w"))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", handlers=handlers)
    if args.cmd == "fetch":
        fetch(args.out, args.email, args.api_key)
    elif args.cmd == "build":
        n_r, n_h = build(args.gb, args.out, args.threads, args.curated)
        logging.info(f"[NUC] library {args.out}: {n_r} rDNA records, {n_h} histone proteins")
    else:
        genes = [g.strip() for g in args.genes.split(",") if g.strip()]
        bad = set(genes) - {"rdna", "histone"}
        if bad:
            ap.error(f"--genes: unknown {', '.join(sorted(bad))}")
        annotate(args.input, args.out, args.lib, genes, args.threads, args.leave_out, args.consensus_identity)
    return 0


if __name__ == "__main__":
    sys.exit(main())
