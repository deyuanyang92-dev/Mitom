#!/usr/bin/env python3
"""mito_core.py - MitoM's own annotation core: protein genes, rRNAs and tRNAs of one mitogenome without MitoZ.

Protein genes: the genome is searched with blastx against the per-gene protein sets of the MITOS2 reference data
(refseq89m/featureProt/<gene>.fas, MIT licence: RefSeq proteins of all Metazoa); the best locus per gene is aligned
with exonerate protein2dna against its best reference proteins (frameshift-aware, intron-less); the start codon is the
nearest in-frame start codon upstream of the alignment (at most <residues before the alignment> + 10 codons), else the
nearest downstream (10 codons), else the gene is 5'-partial; the end is the last aligned codon (plus the stop codon
when it follows directly) - the stop (first in-frame stop, or an abbreviated T / TA before the next tRNA / gene) is
placed by mito_curate.py P3, which knows the neighbours. rRNA ends and NCBI names are left to it as before.
rRNAs: the per-phylum profile HMM of MitoM (models/rrna_hmm, nhmmer; --phylum) when the library has the phylum, else
Infernal cmsearch (fast mode, local) with the MITOS2 covariance models (refseq89m/ncRNA/rrnS.cm, rrnL.cm), best
cluster of hits; an rRNA still missing is left to MitoM's later rRNA HMM step.
tRNAs: cmsearch with the MiTFi tRNA models (refseq89m/ncRNA/trn*.cm; a sensitive pass for models without a hit) and
MiTFi (mitfi_import.jar of the mitos package) for hit selection and anticodons - the same programs MITOS2 runs.

Output: <out>/<prefix>.gbf (GenBank with gene / CDS / rRNA / tRNA features, MitoZ-style names so that every later
MitoM step works unchanged) and <out>/<prefix>.core.tsv (one row per feature with the evidence).

Usage: mito_core.py annotate -i genome.fa -o workdir [--prefix ID] --genetic_code 5 [--threads 2]
                             [--refdir ~/mitos_refs/refseq89m] [--no_trna]
Programs: blastx, exonerate, cmsearch (Infernal 1.1), java (MiTFi); found by mito_env.py.
"""
from __future__ import annotations

import argparse
import logging
import re
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mito_env  # noqa: E402
import mito_rrna_hmm  # noqa: E402  (per-phylum rRNA profile HMMs: nhmmer search, best span)

PCG = ("atp6", "atp8", "cob", "cox1", "cox2", "cox3", "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6")
NAMES = {                                   # MitoZ-style /gene and /product, understood by every later MitoM step
    "atp6": ("ATP6", "ATP synthase F0 subunit 6"), "atp8": ("ATP8", "ATP synthase F0 subunit 8"),
    "cox1": ("COX1", "cytochrome c oxidase subunit I"), "cox2": ("COX2", "cytochrome c oxidase subunit II"),
    "cox3": ("COX3", "cytochrome c oxidase subunit III"), "cob": ("CYTB", "cytochrome b"),
    "nad1": ("ND1", "NADH dehydrogenase subunit 1"), "nad2": ("ND2", "NADH dehydrogenase subunit 2"),
    "nad3": ("ND3", "NADH dehydrogenase subunit 3"), "nad4": ("ND4", "NADH dehydrogenase subunit 4"),
    "nad4l": ("ND4L", "NADH dehydrogenase subunit 4L"), "nad5": ("ND5", "NADH dehydrogenase subunit 5"),
    "nad6": ("ND6", "NADH dehydrogenase subunit 6"),
    "rrnS": ("s-rRNA", "12S ribosomal RNA"), "rrnL": ("l-rRNA", "16S ribosomal RNA"),
}
AA3 = {"A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe", "G": "Gly", "H": "His", "I": "Ile", "K": "Lys",
       "L": "Leu", "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg", "S": "Ser", "T": "Thr", "V": "Val",
       "W": "Trp", "Y": "Tyr"}
BLAST_EVALUE = 10          # blastx is only the locator; short divergent genes (atp8, nad6, nad4l) score weakly
BLAST_MIN_BITS = 35        # weakest locus accepted (exonerate then has to align it)
EXONERATE_MIN_SCORE = 40   # exonerate --score (default 100 loses short divergent genes)
LOCUS_GAP = 150            # blastx hits of one gene closer than this (same strand) form one locus
MIN_REF_COVER = 0.25       # the best alignment must cover this fraction of its reference protein (divergent genes
                           # align only in their conserved core; the start / stop rules complete them)
WINDOW = 300               # nt of genome on each side of a locus given to exonerate
START_UP_EXTRA = 10        # start codon search: upstream codons beyond the residues missing from the alignment
START_DOWN = 10            # ... and downstream codons when none upstream
CM_EVALUE = 0.01           # cmsearch fast pass and MiTFi threshold of that pass (MITOS2 runmitos default --ncev)
CM_EVALUE_SENSITIVE = 0.1  # cmsearch --max pass for models without a hit (MITOS2 sensitive pass)
RRNA_CLUSTER = 500         # rRNA hits of one model on one strand closer than this are one gene
EXONERATE_MODEL = "protein2dna"   # protein2dna:bestfit = whole reference protein aligned (set by --exonerate_model)


@dataclass
class Hit:
    gene: str
    ref: str
    start: int             # genome, 0-based, half-open, plus-strand coordinates
    end: int
    strand: int
    bits: float
    sstart: int            # reference protein residues (1-based)
    send: int
    slen: int


@dataclass
class Feature:
    kind: str              # CDS / rRNA / tRNA
    name: str              # MITOS2-style: cox1, rrnL, trnL1
    start: int             # 0-based half-open
    end: int
    strand: int
    note: str = ""
    partial5: bool = False
    anticodon: str = ""
    acpos: int = -1            # tRNA: 0-based position of the anticodon within the tRNA (MiTFi)
    model: str = ""            # tRNA: the covariance model that found it
    score: float = 0.0


# ── programs and data ────────────────────────────────────────────────────────
def program(name: str, explicit: Optional[str] = None) -> str:
    p = mito_env.find(name, explicit)
    if not p:
        raise SystemExit(f"[core] {name} not found (mito_env.py doctor; install.sh)")
    return p


def mitfi_files(mitos_bin: Optional[str]) -> Tuple[Path, Path]:
    """mitfi_import.jar and mitfiGeneticCodes of the mitos package (searched next to the MITOS2 python)."""
    bins = [Path(mitos_bin)] if mitos_bin else []
    mb = mito_env.mitos_bin(None)
    if mb:
        bins.append(Path(mb))
    for b in bins:
        for p in sorted((b.parent / "lib").glob("python3*/site-packages/mitos/mitfi/parser/mitfi_import.jar")):
            return p, p.parent / "mitfiGeneticCodes"
    for p in HERE.glob("models/mitfi/mitfi_import.jar"):
        return p, p.parent / "mitfiGeneticCodes"
    raise SystemExit("[core] MiTFi (mitfi_import.jar of the mitos package) not found")


def run(cmd: List[str], log: logging.Logger, **kw) -> subprocess.CompletedProcess:
    log.debug("[CMD] " + " ".join(cmd))
    p = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if p.returncode != 0:
        raise RuntimeError(f"{Path(cmd[0]).name} failed ({p.returncode}): {p.stderr.strip()[-400:]}")
    return p


# ── protein genes ────────────────────────────────────────────────────────────
def blastx_gene(blastx: str, genome_fa: Path, db: Path, code: int, gene: str, threads: int, log,
                sensitive: bool = False) -> List[Hit]:
    """blastx of the genome against one gene's protein set. Default: standard settings (1-2 s per gene); sensitive:
    word size 2, BLOSUM45, no composition statistics (3-4 s per gene), run only for genes the default pass misses."""
    extra = ["-seg", "no", "-comp_based_stats", "0", "-word_size", "2", "-matrix", "BLOSUM45"] if sensitive else ["-seg", "no"]
    p = run([blastx, "-query", str(genome_fa), "-db", str(db), "-query_gencode", str(code), "-evalue", str(BLAST_EVALUE),
             *extra, "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore slen",
             "-max_target_seqs", "50", "-num_threads", str(max(1, threads))], log)
    hits = []
    for line in p.stdout.splitlines():
        c = line.split("\t")
        qs, qe = int(c[4]), int(c[5])
        strand = 1 if qs <= qe else -1
        if float(c[9]) >= BLAST_MIN_BITS:
            hits.append(Hit(gene, c[1], min(qs, qe) - 1, max(qs, qe), strand, float(c[9]), int(c[6]), int(c[7]), int(c[10])))
    return hits


def loci(hits: List[Hit]) -> List[List[Hit]]:
    """Same-strand hits closer than LOCUS_GAP form one locus; loci sorted by their best bit score."""
    out: List[List[Hit]] = []
    for h in sorted(hits, key=lambda x: (x.strand, x.start)):
        if out and out[-1][0].strand == h.strand and h.start <= max(x.end for x in out[-1]) + LOCUS_GAP:
            out[-1].append(h)
        else:
            out.append([h])
    return sorted(out, key=lambda L: -max(x.bits for x in L))


def clean_protein(seq: str) -> str:
    """exonerate stops on letters outside the 20 amino acids + X (RefSeq proteins carry J, B, Z, U, O, *)."""
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", seq.upper().rstrip("*"))


def reference_proteins(db_fasta: Path, ids: List[str]) -> Dict[str, str]:
    want = set(ids)
    out = {}
    for r in SeqIO.parse(str(db_fasta), "fasta"):
        if r.id in want:
            out[r.id] = clean_protein(str(r.seq))
            if len(out) == len(want):
                break
    return out


def exonerate_best(exonerate: str, window: str, prots: Dict[str, str], code: int, tmp: Path, log) \
        -> Optional[Tuple[str, int, int, int, int, float]]:
    """Best protein2dna alignment of the reference proteins on the window (plus strand):
    (ref, qab, qae, tab, tae, score); coordinates 0-based half-open."""
    (tmp / "w.fa").write_text(f">w\n{window}\n")
    best = None
    for rid, seq in prots.items():
        (tmp / "q.fa").write_text(f">{rid}\n{seq}\n")
        try:
            p = run([exonerate, "--model", EXONERATE_MODEL, "--geneticcode", str(code), "--bestn", "1", "--showalignment", "no",
                 "--showvulgar", "no", "--score", str(EXONERATE_MIN_SCORE),
                 "--ryo", "RYO\t%qi\t%ql\t%qab\t%qae\t%tab\t%tae\t%s\n",
                 *(["--exhaustive", "yes"] if EXONERATE_MODEL.endswith("bestfit") else []),   # non-local models need full DP
                 str(tmp / "q.fa"), str(tmp / "w.fa")], log)
        except RuntimeError as exc:                        # one bad reference must not lose the gene
            log.warning(f"[core] exonerate with {rid}: {exc}")
            continue
        for line in p.stdout.splitlines():
            if line.startswith("RYO\t"):
                c = line.split("\t")
                score = float(c[7])
                if best is None or score > best[5]:
                    best = (c[1], int(c[3]), int(c[4]), int(c[5]), int(c[6]), score)
    return best


def coding_frame(window: str, tab: int, tae: int, table: CodonTable.CodonTable) -> int:
    """Start of the aligned region moved to the reading frame with the fewest stop codons over [tab, tae): a
    frameshift early in the alignment (exonerate protein2dna) otherwise shifts every codon of the gene
    (NC_086642.1 nad5: 78 internal stops, both ends 5 nt off)."""
    stops = set(table.stop_codons)
    best = min(range(3), key=lambda d: sum(1 for c in range(tab + d, tae - 2, 3) if window[c:c + 3] in stops))
    return tab + best


def place_cds(window: str, tab: int, tae: int, missing_n: int, table: CodonTable.CodonTable) \
        -> Tuple[int, int, bool]:
    """Start and end (0-based half-open) of the CDS on the plus-strand window from the alignment [tab, tae):
    MitoZ's rule for the start (nearest start codon upstream within missing_n + 10 codons, else downstream within
    10, else 5'-partial); the end is the last complete aligned codon - the stop codon (complete, or T / TA before
    the next gene or tRNA) is placed by mito_curate.py P3, which knows the neighbouring genes."""
    starts, stops = set(table.start_codons), set(table.stop_codons)
    a = tab                                                # exonerate's alignment start is the first aligned codon
    start, partial = a, True
    if missing_n > START_UP_EXTRA:
        # a local alignment that misses many N-terminal residues of a full-length reference (divergent N-terminus):
        # the start is expected ~3 x missing_n nt upstream; take the in-frame start codon nearest to that position
        # (within +-10 codons, no stop in between), else fall through to MitoZ's nearest-upstream rule
        target = a - 3 * missing_n
        best = None
        c = a - 3
        while c >= 0 and c >= target - 3 * START_UP_EXTRA and window[c:c + 3] not in stops:
            if window[c:c + 3] in starts and abs(c - target) <= 3 * START_UP_EXTRA and \
                    (best is None or abs(c - target) < abs(best - target)):
                best = c
            c -= 3
        if best is not None:
            start, partial = best, False
    for j in range(0, missing_n + START_UP_EXTRA + 1):     # upstream, nearest first (j = 0: the aligned start itself)
        if not partial:
            break
        c = a - 3 * j
        if c < 0:
            break
        codon = window[c:c + 3]
        if j and codon in stops:
            break
        if codon in starts:
            start, partial = c, False
            break
    if partial:
        for j in range(1, START_DOWN + 1):
            c = a + 3 * j
            if c + 3 > tae or window[c:c + 3] in stops:
                break
            if window[c:c + 3] in starts:
                start, partial = c, False
                break
    end = a + ((tae - a) // 3) * 3                           # last complete codon of the alignment
    if window[end:end + 3] in stops:                         # the stop right after the alignment belongs to the gene
        end += 3
    return start, min(end, len(window)), partial


def window_to_genome(ws: int, we: int, strand: int, s: int, e: int) -> Tuple[int, int]:
    """Plus-strand genome coordinates (0-based half-open) of [s, e) on the window that covers genome [ws, we):
    the window reads forward on +, reverse-complemented on -, so on - the window's start is the genome's end."""
    return (ws + s, ws + e) if strand == 1 else (we - e, we - s)


def annotate_pcgs(genome: str, genome_fa: Path, refdir: Path, code: int, threads: int, tools: dict, tmp: Path, log) \
        -> List[Feature]:
    table = CodonTable.unambiguous_dna_by_id[code]
    n = len(genome)
    rc = str(Seq(genome).reverse_complement())
    feats: List[Feature] = []
    with ThreadPoolExecutor(max(1, threads)) as ex:
        all_hits = dict(zip(PCG, ex.map(lambda g: blastx_gene(tools["blastx"], genome_fa, refdir / "featureProt" / f"{g}.fas",
                                                                 code, g, 1, log), PCG)))
        missing = [g for g in PCG if not all_hits[g]]
        if missing:
            log.info(f"[core] sensitive blastx for {len(missing)} gene(s) without a hit: " + " ".join(missing))
            all_hits.update(zip(missing, ex.map(lambda g: blastx_gene(tools["blastx"], genome_fa, refdir / "featureProt" / f"{g}.fas",
                                                                        code, g, 1, log, sensitive=True), missing)))
    for gene in PCG:
        hits = all_hits[gene]
        if not hits:
            log.info(f"[core] {gene}: no blastx hit")
            continue
        locus = loci(hits)[0]
        strand = locus[0].strand
        ls, le = min(h.start for h in locus), max(h.end for h in locus)
        top = sorted(locus, key=lambda h: -h.bits)
        refs: List[str] = []
        for h in top:
            if h.ref not in refs:
                refs.append(h.ref)
            if len(refs) == 3:
                break
        prots = reference_proteins(refdir / "featureProt" / f"{gene}.fas", refs)
        # window on the gene's own strand: plus interval [ws, we) of the genome, reverse-complemented for minus genes;
        # wide enough for a whole reference protein on either side of the blastx locus (blastx may find only the
        # conserved half of a divergent gene: NC_083122.1 nad4 was placed 552 nt late with a fixed 300-nt margin)
        margin = 3 * max(h.slen for h in locus) + WINDOW
        ws, we = max(0, ls - margin), min(n, le + margin)
        window = genome[ws:we] if strand == 1 else rc[n - we:n - ws]
        best = exonerate_best(tools["exonerate"], window, prots, code, tmp, log)
        if best is None:
            log.info(f"[core] {gene}: exonerate found no alignment")
            continue
        ref, qab, qae, tab, tae, score = best
        reflen = len(prots[ref])
        if (qae - qab) < MIN_REF_COVER * reflen:
            log.info(f"[core] {gene}: alignment covers {qae - qab} of {reflen} residues: skipped")
            continue
        tab = coding_frame(window, tab, tae, table)        # a frameshift near the start can put the alignment off frame
        s, e, partial = place_cds(window, tab, tae, qab, table)
        gs, ge = window_to_genome(ws, we, strand, s, e)
        feats.append(Feature("CDS", gene, gs, ge, strand, f"ref {ref} score {score:.0f} cover {qae - qab}/{reflen}",
                             partial5=partial, score=score))
        log.info(f"[core] {gene}: {gs + 1}..{ge} ({'+' if strand == 1 else '-'}) ref {ref} score {score:.0f}"
                 + (" start not found" if partial else ""))
    return feats


# ── RNAs ─────────────────────────────────────────────────────────────────────
def cmsearch(cm: str, model: Path, genome_fa: Path, out: Path, sensitive: bool, log, rrna: bool = False) -> Path:
    """cmsearch as MITOS2 runs it (local mode, --cpu 1 --notextw --smxsize 80000); fast pass --noF4b -E 0.01,
    sensitive pass --max -E 0.1 (rRNA: --nohmm). Local mode finds the truncated tRNAs of nematodes."""
    cmd = [cm, "--cpu", "1", "--notextw", "--smxsize", "80000", "--tblout", str(out.with_suffix(".tbl"))]
    if sensitive:
        cmd += ["--nohmm" if rrna else "--max", "-E", str(CM_EVALUE_SENSITIVE)]
    else:
        cmd += ["--noF4b", "-E", str(CM_EVALUE)]
    p = run(cmd + [str(model), str(genome_fa)], log)
    out.write_text(p.stdout)
    return out


def tbl_hits(tbl: Path) -> List[Tuple[int, int, int, float, float]]:
    """(start, end, strand, bits, evalue) of a cmsearch --tblout file, 0-based half-open plus coordinates."""
    hits = []
    for line in tbl.read_text().splitlines():
        if line.startswith("#"):
            continue
        c = line.split()
        a, b = int(c[7]), int(c[8])
        strand = 1 if c[9] == "+" else -1
        hits.append((min(a, b) - 1, max(a, b), strand, float(c[14]), float(c[15])))
    return hits


def rrna_hmm(lib: Path, phylum: Optional[str], model: str) -> Optional[Path]:
    """The per-phylum profile HMM of model (rrnS / rrnL) in lib for `phylum`, None when there is none."""
    if not phylum or not lib.is_dir():
        return None
    name = {p.name.split(".")[0].lower(): p.name.split(".")[0] for p in lib.glob("*.hmm")}.get(phylum.lower())
    p = lib / f"{name}.{model}.hmm" if name else None
    return p if p and p.is_file() else None


def annotate_rrnas(genome_fa: Path, refdir: Path, tools: dict, tmp: Path, log, phylum: Optional[str] = None,
                   rrna_lib: Optional[Path] = None) -> List[Feature]:
    """12S / 16S: the per-phylum profile HMM (nhmmer, ~0.1 GB, seconds) when the library has the record's phylum,
    else the MITOS2 covariance model (cmsearch, 0.4-0.6 GB, 10 s). On the 24-record benchmark both find the same
    genes (leave-out library); R1 of mito_curate.py then fits the ends to the flanking genes."""
    feats = []
    lib = rrna_lib or mito_rrna_hmm.LIB
    for model in ("rrnS", "rrnL"):
        hmm = rrna_hmm(lib, phylum, model)
        if hmm is not None:
            span = mito_rrna_hmm.best_span(mito_rrna_hmm.nhmmer(hmm, genome_fa))
            if span is not None:
                feats.append(Feature("rRNA", model, span.start, span.end, span.strand,
                                     f"nhmmer {hmm.name} bits {span.score:.0f}", score=span.score))
                log.info(f"[core] {model}: {span.start + 1}..{span.end} ({'+' if span.strand == 1 else '-'}) by {hmm.name}")
                continue
            log.info(f"[core] {model}: no nhmmer hit with {hmm.name}; trying the covariance model")
        cm = refdir / "ncRNA" / f"{model}.cm"
        hits = tbl_hits(cmsearch(tools["cmsearch"], cm, genome_fa, tmp / f"{model}.out", False, log, rrna=True).with_suffix(".tbl"))
        if not hits:
            # no sensitive (--nohmm) pass: 2.4 GB and 2 min per model; a missing rRNA is found later by MitoM's
            # per-phylum rRNA HMM step (mito_rrna_hmm.py, nhmmer)
            log.info(f"[core] {model}: no cmsearch hit (left to the rRNA HMM step)")
            continue
        clusters: List[List[tuple]] = []
        for h in sorted(hits, key=lambda x: (x[2], x[0])):
            if clusters and clusters[-1][0][2] == h[2] and h[0] <= max(x[1] for x in clusters[-1]) + RRNA_CLUSTER:
                clusters[-1].append(h)
            else:
                clusters.append([h])
        best = max(clusters, key=lambda c: sum(x[3] for x in c))
        s, e, strand = min(x[0] for x in best), max(x[1] for x in best), best[0][2]
        feats.append(Feature("rRNA", model, s, e, strand, f"cmsearch {len(best)} hit(s), bits {max(x[3] for x in best):.0f}",
                             score=max(x[3] for x in best)))
        log.info(f"[core] {model}: {s + 1}..{e} ({'+' if strand == 1 else '-'})")
    return feats


def annotate_trnas(genome: str, genome_fa: Path, refdir: Path, code: int, threads: int, tools: dict, tmp: Path, log) \
        -> List[Feature]:
    """MITOS2's tRNA step: cmsearch per MiTFi model (fast pass, then --max for models without a hit), MiTFi for
    selection, structure and anticodon."""
    models = sorted((refdir / "ncRNA").glob("trn*.cm"))
    seqname = genome_fa.name

    def search(model: Path, sensitive: bool) -> Path:
        return cmsearch(tools["cmsearch"], model, genome_fa, tmp / f"{seqname}_{model.stem}.cmout", sensitive, log)

    jar, codes = tools["mitfi"]

    def mitfi(outs: List[Path], evalue: float) -> List[Feature]:
        cmd = [tools["java"], "-jar", str(jar), "-code", str(code), "-codefile", str(codes), "-structure", "-genome",
               str(genome_fa), "-glength", str(len(genome)), "-evalue", str(evalue), "-overlap", "10"]
        for o in outs:
            cmd += ["-in", str(o)]
        return parse_mitfi(run(cmd, log).stdout, log)

    # MITOS2's procedure (runmitos.py): fast pass judged at E <= 0.001 - a model whose only hits are weaker counts
    # as empty and is searched again with --max; the final MiTFi run accepts E <= 0.1
    with ThreadPoolExecutor(max(1, threads)) as ex:
        outs = list(ex.map(lambda m: search(m, False), models))
    feats = mitfi(outs, CM_EVALUE)
    accepted = {f.model for f in feats}
    missing = [m for m in models if m.stem not in accepted]
    if missing:
        log.info(f"[core] tRNA sensitive pass for {len(missing)} model(s): " + " ".join(m.stem for m in missing))
        with ThreadPoolExecutor(max(1, threads)) as ex:
            list(ex.map(lambda m: search(m, True), missing))
        feats = mitfi(outs, CM_EVALUE_SENSITIVE)
    log.info(f"[core] tRNAs: {len(feats)}")
    return feats


def trna_name(aa: str, anticodon: str, model: str) -> str:
    """MITOS2 naming: the amino acid MiTFi derives from the anticodon; the two Leu / Ser genes by codon family
    (L1 = CUN: anticodon UAG; L2 = UUR: UAA; S1 = AGN: GCU / UCU; S2 = UCN: UGA), else the model's number."""
    aa = (aa or "").strip().upper()[:1]
    ac = anticodon.replace("u", "t")
    if aa == "L":
        return "trnL1" if ac == "tag" else "trnL2" if ac == "taa" else ("trnL" + (re.search(r"trnL([12])", model) or [None, ""])[1])
    if aa == "S":
        return "trnS1" if ac in ("gct", "tct") else "trnS2" if ac == "tga" else ("trnS" + (re.search(r"trnS([12])", model) or [None, ""])[1])
    if aa and aa in AA3:
        return "trn" + aa
    m = re.search(r"(trn[A-Z](?:[12])?)", model)
    return m.group(1) if m else ""


def parse_mitfi(text: str, log) -> List[Feature]:
    """MiTFi output (what mitos.mitfi.parse reads): three lines per feature, the first
    '>acc|code|start|stop|strand|qstart|qstop|bits|evalue|pvalue|gc|acpos|anticodon|aa|model|tophit|copy|g/l',
    then the sequence and the structure."""
    feats: List[Feature] = []
    for line in text.splitlines():
        line = line.strip()
        if not line.startswith((">", "|")):
            continue
        c = line[1:].split("|")
        try:
            start, stop, strand = int(c[2]) - 1, int(c[3]), int(c[4])
            bits, evalue = float(c[7]), float(c[8])
            acpos = int(c[11]) - 1
            anticodon = c[12].strip().lower().replace("u", "t") if acpos >= 0 else ""
            aa, model = c[13], c[14]
            tophit, copy = int(c[15]), int(c[16])
        except (IndexError, ValueError) as exc:
            log.warning(f"[core] MiTFi line not understood ({exc}): {line[:80]}")
            continue
        if not (tophit or (copy and evalue <= CM_EVALUE)):    # MITOS2: a copy (second hit of a model) only when E <= 0.001
            continue
        name = trna_name(aa, anticodon, model)
        if not name or start >= stop:
            continue
        mm = re.search(r"(trn[A-Z](?:[12])?)\.cm", model)
        feats.append(Feature("tRNA", name, start, stop, strand, f"MiTFi {model} bits {bits:.1f} E {evalue:.1e}"
                             + (" copy" if copy and not tophit else ""), anticodon=anticodon, acpos=acpos, score=bits,
                             model=mm.group(1) if mm else model))
    # overlapping tRNAs of different models: keep the better score (MITOS2 keeps both as copies of a locus)
    feats.sort(key=lambda f: -f.score)
    kept: List[Feature] = []
    for f in feats:
        if all(min(f.end, k.end) - max(f.start, k.start) <= 10 for k in kept):
            kept.append(f)
    return sorted(kept, key=lambda f: f.start)


# ── output ───────────────────────────────────────────────────────────────────
def write_genbank(genome: str, feats: List[Feature], code: int, prefix: str, topology: str, out: Path) -> Path:
    rec = SeqRecord(Seq(genome), id=prefix, name=prefix[:16], description="")
    rec.annotations.update(molecule_type="DNA", topology=topology, data_file_division="INV")
    rec.features.append(SeqFeature(FeatureLocation(0, len(genome), 1), type="source",
                                   qualifiers={"organism": [prefix], "organelle": ["mitochondrion"], "mol_type": ["genomic DNA"]}))
    for f in sorted(feats, key=lambda x: x.start):
        if f.kind == "CDS":
            gene, product = NAMES[f.name]
            start = BeforePosition(f.start) if (f.partial5 and f.strand == 1) else f.start
            loc = FeatureLocation(start, f.end, f.strand)
            q = {"gene": [gene], "product": [product], "codon_start": ["1"], "transl_table": [str(code)],
                 "note": [f.note]}
            if f.partial5:
                q["note"].append("start codon not determined")
        elif f.kind == "rRNA":
            gene, product = NAMES[f.name]
            loc = FeatureLocation(f.start, f.end, f.strand)
            q = {"gene": [gene], "product": [product], "note": [f.note]}
        else:
            gene = f.name
            loc = FeatureLocation(f.start, f.end, f.strand)
            aa3 = AA3.get(f.name[3], "Xaa")
            q = {"gene": [gene], "product": [f"tRNA-{aa3}"], "note": [f.note]}
            if f.anticodon:                                # as MITOS2 output is written: mito_curate reads this note
                q["note"] = [f"anticodon: {f.anticodon.upper().replace('T', 'U')}", f.note]
        rec.features.append(SeqFeature(loc, type="gene", qualifiers={"gene": [gene]}))
        rec.features.append(SeqFeature(loc, type=f.kind, qualifiers=q))
    SeqIO.write(rec, str(out), "genbank")
    return out


def annotate(args) -> int:
    log = logging.getLogger("core")
    t0 = time.time()
    fa = Path(args.input).resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    rec = next(SeqIO.parse(str(fa), "fasta"))
    genome = str(rec.seq).upper()
    topology = "linear" if re.search(r"topology\s*=\s*linear", rec.description, re.I) else "circular"
    prefix = args.prefix or re.sub(r"[^A-Za-z0-9_.-]", "_", rec.id)[:16]
    refdir = Path(args.refdir or mito_env.mitos_refdir(None) or "").expanduser()
    if (refdir / "refseq89m").is_dir():                     # mito_env gives the directory that holds refseq89m
        refdir = refdir / "refseq89m"
    if not (refdir / "featureProt").is_dir() or not (refdir / "ncRNA").is_dir():
        raise SystemExit(f"[core] reference data not found: {refdir} (refseq89m with featureProt/ and ncRNA/)")
    tools = {"blastx": program("blastx", args.blastx), "exonerate": program("exonerate", args.exonerate),
             "cmsearch": program("cmsearch", args.cmsearch)}
    if not args.no_trna:
        tools["java"] = program("java", args.java)
        tools["mitfi"] = mitfi_files(args.mitos_bin)
    code = int(args.genetic_code)
    with tempfile.TemporaryDirectory(prefix="mito_core.") as td:
        tmp = Path(td)
        gfa = tmp / fa.name
        gfa.write_text(f">{prefix}\n{genome}\n")
        feats = annotate_pcgs(genome, gfa, refdir, code, args.threads, tools, tmp, log)
        feats += annotate_rrnas(gfa, refdir, tools, tmp, log, args.phylum, Path(args.rrna_lib) if args.rrna_lib else None)
        if not args.no_trna:
            feats += annotate_trnas(genome, gfa, refdir, code, args.threads, tools, tmp, log)
        if args.keep_tmp:
            shutil.copytree(tmp, out / "core_tmp", dirs_exist_ok=True)
    gb = write_genbank(genome, feats, code, prefix, topology, out / f"{prefix}.gbf")
    with (out / f"{prefix}.core.tsv").open("w") as fh:
        fh.write("kind\tname\tstart\tend\tstrand\tscore\tanticodon\tnote\n")
        for f in sorted(feats, key=lambda x: x.start):
            fh.write(f"{f.kind}\t{f.name}\t{f.start + 1}\t{f.end}\t{'+' if f.strand == 1 else '-'}\t{f.score:.1f}\t"
                     f"{f.anticodon}\t{f.note}\n")
    n = {k: sum(1 for f in feats if f.kind == k) for k in ("CDS", "rRNA", "tRNA")}
    log.info(f"[core] done in {time.time() - t0:.0f} s: {n['CDS']} protein genes, {n['rRNA']} rRNAs, {n['tRNA']} tRNAs -> {gb}")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    global EXONERATE_MODEL
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("annotate", help="annotate one mitogenome FASTA")
    a.add_argument("-i", "--input", required=True, help="FASTA with one sequence")
    a.add_argument("-o", "--out", required=True, help="output directory")
    a.add_argument("--prefix", default="", help="record / file name (default: the FASTA id)")
    a.add_argument("--genetic_code", type=int, default=5)
    a.add_argument("--threads", type=int, default=2)
    a.add_argument("--refdir", default=None, help="MITOS2 reference data (refseq89m); default from mito_env")
    a.add_argument("--no_trna", action="store_true", help="skip the tRNA step")
    a.add_argument("--phylum", default=None, help="phylum of the record: selects the per-phylum rRNA profile HMM")
    a.add_argument("--rrna_lib", default=None, help="rRNA profile HMM library (default models/rrna_hmm)")
    a.add_argument("--keep_tmp", action="store_true", help="keep blastx / cmsearch / MiTFi files in <out>/core_tmp")
    a.add_argument("--exonerate_model", default=EXONERATE_MODEL, help="exonerate model (protein2dna, protein2dna:bestfit)")
    a.add_argument("--blastx"), a.add_argument("--exonerate"), a.add_argument("--cmsearch"), a.add_argument("--java")
    a.add_argument("--mitos_bin", default=None, help="MITOS2 environment bin (to find MiTFi); default from mito_env")
    a.add_argument("--log", default=None)
    args = p.parse_args(argv)
    EXONERATE_MODEL = args.exonerate_model
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stderr)] + ([logging.FileHandler(args.log)] if args.log else []))
    return annotate(args)


if __name__ == "__main__":
    sys.exit(main())
