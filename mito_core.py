#!/usr/bin/env python3
"""mito_core.py - MitoM's own annotation core: protein genes, rRNAs and tRNAs of one mitogenome.

Protein genes: the genome is searched with blastx against the per-gene protein sets of the MITOS2 reference data
(refseq89m/featureProt/<gene>.fas, MIT licence: RefSeq proteins of all Metazoa); the best locus per gene is aligned
with exonerate protein2dna against its best reference proteins (frameshift-aware, intron-less); the start codon is the
one nearest the expected position when the alignment misses N-terminal residues, else the nearest in-frame start codon
upstream (at most <residues before the alignment> + 10 codons), else the nearest downstream (10 codons), else the gene
is 5'-partial; a CDS lying mostly on a better-scoring CDS of another gene (paralog hit) is dropped; the end is the last aligned codon (plus the stop codon
when it follows directly) - the stop (first in-frame stop, or an abbreviated T / TA before the next tRNA / gene) is
placed by mito_curate.py P3, which knows the neighbours. rRNA ends and NCBI names are left to it as before.
rRNAs: the per-phylum profile HMM of MitoM (models/rrna_hmm, nhmmer; --phylum) when the library has the phylum, else
Infernal cmsearch (fast mode, local) with the MITOS2 covariance models (refseq89m/ncRNA/rrnS.cm, rrnL.cm), best
cluster of hits; an rRNA still missing is left to MitoM's later rRNA HMM step.
tRNAs (--trna auto, default): MitoM-ML (trna_ml.py) when its models are installed - MITOS2 tRNA CMs plus the clade CMs
on --lineage (default Metazoa;<--phylum>), ARWEN / ARAGORN as extra evidence when found, a learned caller and learned
ends (unseen-family benchmark: recall 98.5 % vs 94.9 % for the MiTFi CMs, projects/invert_trna/REPORT.md); MiTFi (as
below) gives a second opinion: calls with p < 0.9 need its call of the same tRNA, a weaker second call of a family
needs p >= 0.9 (trna_second_opinion), the families MitoM-ML does not call are taken from it (rescue_trnas), and
vertebrates take its ends (mitfi_ends); else, or
with --trna mitfi: cmsearch with the MiTFi tRNA models (refseq89m/ncRNA/trn*.cm; a sensitive pass for models without a
hit, between the protein genes / rRNAs) and MiTFi (mitfi_import.jar of the mitos package) for hit selection and
anticodons - the programs MITOS2 runs. The tRNA step is mito_trna.core_trnas (mito_trna.py holds all tRNA code).

Circular records are searched as sequence + its first ORIGIN_EXT nt, so genes across the origin are found whole and
written as joins (MITOS2 searches two rotations for the same reason). Several FASTA sequences = contigs of one sample.

Output: <out>/<prefix>.gbf (GenBank with gene / CDS / rRNA / tRNA features, the core's names such as nad1 and
trnN(guu), which every later MitoM step understands; mito_curate N1 writes NCBI's names) and <out>/<prefix>.core.tsv
(one row per feature with the evidence).

Usage: mito_core.py annotate -i genome.fa -o workdir [--prefix ID] --genetic_code 5 [--threads 2]
                             [--refdir ~/mitos_refs/refseq89m] [--no_trna]
Programs: blastx, exonerate, cmsearch (Infernal 1.1), java (MiTFi); found by mito_env.py.
"""
from __future__ import annotations

import argparse
import logging
import os
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

from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mito_env  # noqa: E402
import mito_rrna_hmm  # noqa: E402,F401  (the rRNA step: per-phylum profile HMMs, MITOS2's covariance models)
# the rRNA step, moved to mito_rrna_hmm.py, stays importable from here (tests)
from mito_rrna_hmm import annotate_rrnas, cm_len, rrna_hmm, tbl_hits  # noqa: E402,F401
from mito_common import contig_record_id, CoreFeature as Feature, debug_mode, positions, read_fasta, run_or_raise  # noqa: E402,F401
from mito_common import topology_from_header  # noqa: E402
import mito_trna  # noqa: E402  (the tRNA step: MitoM-ML checked by MiTFi, mito_trna.core_trnas)
# the tRNA step's names, moved to mito_trna.py, stay importable from here (benchmark scripts, tests)
from mito_trna import (annotate_trnas, annotate_trnas_ml, CM_ALN, CM_EVALUE, CM_EVALUE_SENSITIVE,  # noqa: E402,F401
                       CM_HIT, CM_ROW, CM_SENSITIVE, cmsearch, COPY_EVALUE, gene_free, mitfi_ends,  # noqa: F401
                       mitfi_files, MITFI_TRNA_OVL, mitfi_trna_name as trna_name, NO_SENSITIVE_TRNA,  # noqa: F401
                       parse_mitfi, rescue_trnas, same_trna, shift_cmout, TRNA_NAMES, TRNA_RESCUE_OVL,  # noqa: F401
                       TRNA_SENSITIVE_FLANK, TRNA_SURE_P, trna_second_opinion, vertebrate)  # noqa: F401
# the core's /gene and /product (understood by every later MitoM step) and the amino acids, one letter -> three
from mito_names import AA1_TO_3 as AA3, CORE_NAMES as NAMES  # noqa: E402

PCG = ("atp6", "atp8", "cob", "cox1", "cox2", "cox3", "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6")
BLAST_EVALUE = 10          # blastx is only the locator; short divergent genes (atp8, nad6, nad4l) score weakly
BLAST_MIN_BITS = 35        # weakest locus accepted (exonerate then has to align it)
EXONERATE_MIN_SCORE = 40   # exonerate --score (default 100 loses short divergent genes)
LOCUS_GAP = 150            # blastx hits of one gene closer than this (same strand) form one locus
MIN_REF_COVER = 0.25       # the best alignment must cover this fraction of its reference protein (divergent genes
                           # align only in their conserved core; the start / stop rules complete them)
WINDOW = 300               # nt of genome on each side of a locus given to exonerate
START_UP_EXTRA = 10        # start codon search: upstream codons beyond the residues missing from the alignment
START_DOWN = 10            # ... and downstream codons when none upstream
BLASTX_PLACE_COVER = 0.5   # a blastx hit covering this fraction of its reference places a gene exonerate cannot align
BLASTX_PLACE_REFS = 10     # ... or when this many references support the locus (covering 80 % of the best one)
BLASTX_PLACE_BITS = 50     # ... when it scores at least this (a 40-bit atp8 hit in a flatworm, which has no atp8, does not)
MAX_LOCI = 4               # loci tried per gene when the best one is a paralog hit on another gene
# a short gene no search places, ending where its neighbour starts in >= 90 % of the RefSeq records of these phyla (1,050
# records, 2026-10-23): the ORF whose stop codon lies 25 nt before to 10 nt after the neighbour's start, start codon
# chosen for the typical length (in RefSeq: 66 % of atp8, 81 % of nad4l exact, nearly all others overlapping)
ADJACENT = {"atp8": ("atp6", {"arthropoda", "echinodermata"}, 162, 90, 240),
            "nad4l": ("nad4", {"arthropoda", "annelida", "chordata", "nemertea"}, 297, 240, 330)}
COPY_FRAC = 0.9            # a second locus of a gene scoring this fraction of the placed copy (blastx and exonerate) and
                           # aligned complete is a second copy - kept when >= 2 genes have one (a duplicated block, both
                           # copies in NCBI: NC_009093.1 mite, 9 genes) or when it is a COPY_IDENTITY copy of the gene
COPY_IDENTITY = 0.95       # a single second copy this identical (length within 3 nt) is a real gene copy: 1,050 NCBI
                           # records 2026-10-29: 30 of 38 identical protein-gene copies annotated twice by NCBI (heron /
                           # cormorant / turaco nad6, Hexamermis atp6 ...), near-identical ones too (Ardeola, Botaurus
                           # nad6, Todarodes cox1 / cox2); dropping it lost 12 of MitoM's 55 missing genes
N_REFS = 3                 # reference proteins aligned per gene (the best blastx hits of the locus)
NO_START_DOWN = {"echinodermata": {"nad4l"}}   # phylum -> genes left 5'-partial at the alignment start instead of
                           # moved to a downstream start codon: echinoderm nad4l starts right after the preceding tRNA with
                           # no start codon (NCBI marks it 5'-partial; 5 of 5 echinoderm records, 2 of them were moved
                           # 27 nt downstream; for all genes of all phyla the fallback is right more often)
START_EXPECTED_MIN = 2     # residues missing from the alignment from which the start is sought at the expected position
EXONERATE_MODEL = "protein2dna"   # protein2dna:bestfit = whole reference protein aligned (set by --exonerate_model)
EXONERATE_CODES = {1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23}   # built into exonerate 2.4.0; others
                           # (24, 33: Hemichordata) are passed as the 64-letter table ("No built in genetic code ... [24]")
ORIGIN_EXT = 3000          # circular records are searched as seq + seq[:ORIGIN_EXT]: a gene across the origin is
                           # contiguous there (MITOS2 searches two rotations; review 2026-10-22: cox1 cut 60/40 by the
                           # origin lost 614 nt, trnF / rrnL across it were lost or truncated)
PARALOG_MIN_OVERLAP = 60   # a CDS overlapping a better-scoring CDS of another gene by more than this and by >= half of
                           # the shorter one is a paralog hit (atp8 inside atp6, nad2 inside nad5 when the real gene is
                           # absent; MitoZ genomic_cluster.pl / MITOS2 _nonoverlapping_greedy drop these too)
BAD_REFS = {"NC_008151:cob_0-a--1-18769-16822"}   # refseq89m cob.fas entry that is a cox1 protein (87.8 % identity
                           # to cox1, no cob hit): gave a false CYTB inside COX1 (stress test 2026-10-22)
EXCLUDE_REFS: set = set()  # accessions whose own proteins are not used (--exclude_ref / $MITOM_EXCLUDE_REF; benchmark
                           # leave-out: refseq89m holds the proteins of 8,012 RefSeq records)


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


# ── programs and data ────────────────────────────────────────────────────────
def program(name: str, explicit: Optional[str] = None) -> str:
    p = mito_env.find(name, explicit)
    if not p:
        raise SystemExit(f"[core] {name} not found (mito_env.py doctor; install.sh)")
    return p


def run(cmd: List[str], log: logging.Logger, **kw) -> subprocess.CompletedProcess:
    """A program of the core (mito_common.run_or_raise: logged at DEBUG; a failure raises RuntimeError with the end
    of its stderr); log is not used (kept for the scripts that call it)."""
    return run_or_raise(cmd, **kw)


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
        if float(c[9]) >= BLAST_MIN_BITS and c[1] not in BAD_REFS and c[1].split(":")[0] not in EXCLUDE_REFS:
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
    for r in read_fasta(db_fasta):
        if r.id in want:
            out[r.id] = clean_protein(str(r.seq))
            if len(out) == len(want):
                break
    return out


def exonerate_code(code: int) -> str:
    """--geneticcode value: the id when exonerate has the table built in, else the 64 amino acids in TCAG order."""
    if code in EXONERATE_CODES:
        return str(code)
    t = CodonTable.unambiguous_dna_by_id[code]
    return "".join("*" if c in t.stop_codons else t.forward_table[c]
                   for c in (a + b + d for a in "TCAG" for b in "TCAG" for d in "TCAG"))


def exonerate_best(exonerate: str, window: str, prots: Dict[str, str], code: int, tmp: Path, log,
                   model: Optional[str] = None) -> Optional[Tuple[str, int, int, int, int, float]]:
    """Best alignment of the reference proteins on the window (plus strand) with the EXONERATE_MODEL (or `model`,
    e.g. protein2genome for a gene split by an intron): (ref, qab, qae, tab, tae, score); 0-based half-open."""
    model = model or EXONERATE_MODEL
    (tmp / "w.fa").write_text(f">w\n{window}\n")
    best = None
    for rid, seq in prots.items():
        (tmp / "q.fa").write_text(f">{rid}\n{seq}\n")
        try:
            p = run([exonerate, "--model", model, "--geneticcode", exonerate_code(code), "--bestn", "1", "--showalignment", "no",
                 "--showvulgar", "no", "--score", str(EXONERATE_MIN_SCORE),
                 "--ryo", "RYO\t%qi\t%ql\t%qab\t%qae\t%tab\t%tae\t%s\t%tS\n",
                 *(["--exhaustive", "yes"] if model.endswith("bestfit") else []),   # non-local models need full DP
                 str(tmp / "q.fa"), str(tmp / "w.fa")], log)
        except RuntimeError as exc:                        # one bad reference must not lose the gene
            log.warning(f"[core] exonerate with {rid}: {exc}")
            continue
        for line in p.stdout.splitlines():
            if line.startswith("RYO\t"):
                c = line.split("\t")
                if c[8].strip() == "-":                    # reverse-strand alignment (another copy in the window): the
                    continue                               # window is already on the gene's strand (review: start > end crash)
                score = float(c[7])
                if best is None or score > best[5]:
                    best = (c[1], int(c[3]), int(c[4]), int(c[5]), int(c[6]), score)
    return best


def coding_frame(window: str, tab: int, tae: int, table: CodonTable.CodonTable) -> int:
    """Start of the aligned region moved to the reading frame with the fewest stop codons over [tab, tae): a
    frameshift early in the alignment (exonerate protein2dna) otherwise shifts every codon of the gene
    (NC_086642.1 nad5: 78 internal stops, both ends 5 nt off)."""
    stops = set(table.stop_codons)
    n = [sum(1 for c in range(tab + d, tae - 2, 3) if window[c:c + 3] in stops) for d in range(3)]
    if n[0] < 2:               # exonerate aligns through a single stop (premature stop / sequencing error): keep its frame,
        return tab             # mito_curate reports the internal stop (a short gene's other frames may have none)
    return tab + min(range(3), key=lambda d: n[d])


def place_cds(window: str, tab: int, tae: int, missing_n: int, table: CodonTable.CodonTable,
              downstream: bool = True) -> Tuple[int, int, bool]:
    """Start and end (0-based half-open) of the CDS on the plus-strand window from the alignment [tab, tae):
    the start codon nearest the expected position when >= START_EXPECTED_MIN residues are missing, else MitoZ's rule
    (nearest start codon upstream within missing_n + 10 codons, else downstream within 10, else 5'-partial); the end is the last complete aligned codon - the stop codon (complete, or T / TA before
    the next gene or tRNA) is placed by mito_curate.py P3, which knows the neighbouring genes."""
    starts, stops = set(table.start_codons), set(table.stop_codons)
    a = tab                                                # exonerate's alignment start is the first aligned codon
    start, partial = a, True
    if missing_n >= START_EXPECTED_MIN:
        # an alignment that misses N-terminal residues of a full-length reference (divergent N-terminus): the start is
        # expected ~3 x missing_n nt upstream; take the in-frame start codon nearest to that position (within +-10
        # codons, no stop in between), else fall through to MitoZ's nearest-upstream rule (whose j = 0 takes the
        # aligned codon itself when it is a start codon: 18 more exact starts of 1045 genes, 82 records, leave-out)
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
    if partial and downstream:
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


def feature_nt(genome: str, f: Feature) -> str:
    """Nucleotide sequence of a feature on the searched genome (plus-strand coordinates; minus strand reverse-complemented)."""
    s = genome[f.start:f.end]
    return s if f.strand == 1 else str(Seq(s).reverse_complement())


def copy_identity(genome: str, a: Feature, b: Feature) -> float:
    """Fraction of identical positions between two features of one gene, 0 when their lengths differ by more than 3 nt."""
    sa, sb = feature_nt(genome, a), feature_nt(genome, b)
    if not sa or not sb or abs(len(sa) - len(sb)) > 3:
        return 0.0
    return sum(x == y for x, y in zip(sa, sb)) / max(len(sa), len(sb))


def drop_paralogs(feats: List[Feature], n: int, log) -> List[Feature]:
    """Best score first; a CDS lying mostly (>= half of the shorter, > PARALOG_MIN_OVERLAP nt) on a better CDS of
    another gene is a paralog hit placed where the real gene is absent or undetectable."""
    kept: List[Feature] = []
    for f in sorted(feats, key=lambda x: -x.score):
        pf = positions(f, n)
        k = next((k for k in kept if len(pf & positions(k, n)) > PARALOG_MIN_OVERLAP
                  and len(pf & positions(k, n)) >= 0.5 * min(len(pf), k.end - k.start)), None)
        if k is None:
            kept.append(f)
        else:
            log.info(f"[core] {f.name}: {f.start + 1}..{f.end} (score {f.score:.0f}) lies on {k.name} "
                     f"(score {k.score:.0f}): paralog hit, dropped")
    return sorted(kept, key=lambda f: f.start)


def annotate_pcgs(genome: str, genome_fa: Path, refdir: Path, code: int, threads: int, tools: dict, tmp: Path, log,
                  n: Optional[int] = None, phylum: Optional[str] = None) -> List[Feature]:
    """genome = the searched sequence (a circular record with ORIGIN_EXT nt appended), n = the record's length:
    loci starting in the appended part repeat the first ORIGIN_EXT nt and are skipped; a gene across the origin ends
    beyond n (written as a join)."""
    table = CodonTable.unambiguous_dna_by_id[code]
    no_down = NO_START_DOWN.get((phylum or "").lower(), set())
    N = len(genome)
    n = n or N
    rc = str(Seq(genome).reverse_complement())
    # one blastx process at a time with all threads: as fast as 2 x 2 threads (13 genes ~15 s), one ~80-MB process
    # instead of two (peak memory)
    bt = max(1, threads)
    with ThreadPoolExecutor(1) as ex:
        all_hits = dict(zip(PCG, ex.map(lambda g: blastx_gene(tools["blastx"], genome_fa, refdir / "featureProt" / f"{g}.fas",
                                                                 code, g, bt, log), PCG)))
        missing = [g for g in PCG if not all_hits[g]]
        if missing:
            log.info(f"[core] sensitive blastx for {len(missing)} gene(s) without a hit: " + " ".join(missing))
            all_hits.update(zip(missing, ex.map(lambda g: blastx_gene(tools["blastx"], genome_fa, refdir / "featureProt" / f"{g}.fas",
                                                                        code, g, bt, log, sensitive=True), missing)))
    def place(gene: str, locus: List[Hit], mask: frozenset = frozenset()) -> Optional[Feature]:
        """The gene at one blastx locus; mask = genome positions (of genes already placed) hidden from exonerate."""
        strand = locus[0].strand
        ls, le = min(h.start for h in locus), max(h.end for h in locus)
        top = sorted(locus, key=lambda h: -h.bits)
        refs: List[str] = []
        for h in top:
            if h.ref not in refs:
                refs.append(h.ref)
            if len(refs) == 2 * N_REFS:
                break
        # window on the gene's own strand: plus interval [ws, we) of the genome, reverse-complemented for minus genes;
        # wide enough for a whole reference protein on either side of the blastx locus (blastx may find only the
        # conserved half of a divergent gene: NC_083122.1 nad4 was placed 552 nt late with a fixed 300-nt margin)
        margin = 3 * max(h.slen for h in locus) + WINDOW
        ws, we = max(0, ls - margin), min(N, le + margin)
        window = "".join("N" if i % n in mask else genome[i] for i in range(ws, we)) if mask else genome[ws:we]
        if strand == -1:
            window = str(Seq(window).reverse_complement()) if mask else rc[N - we:N - ws]
        placed, why = None, "exonerate found no alignment"
        for batch in (refs[:N_REFS], refs[N_REFS:]):
            # the next N_REFS references only when the best ones gave no complete gene (82 records, leave-out:
            # 3 more genes found and 4 more exact starts with 6 references; always 6 doubled the time)
            if not batch or placed is not None and not placed[-1]:
                break
            prots = reference_proteins(refdir / "featureProt" / f"{gene}.fas", batch)
            best = exonerate_best(tools["exonerate"], window, prots, code, tmp, log)
            if best is None:
                continue
            ref, qab, qae, tab, tae, score = best
            reflen = len(prots[ref])
            if (qae - qab) < MIN_REF_COVER * reflen:
                why = f"alignment covers {qae - qab} of {reflen} residues"
                continue
            tab = coding_frame(window, tab, tae, table)    # a frameshift near the start can put the alignment off frame
            s, e, partial = place_cds(window, tab, tae, qab, table, gene not in no_down)
            if placed is None or placed[-1] and (not partial or score > placed[5]):
                placed = (ref, qab, qae, reflen, s, e, score, partial)
        if placed is None:
            # exonerate aligns only a fragment of a divergent gene that blastx covers end to end (flatworm atp6 of
            # NC_062047.1: residues 4-188 of 189 by blastx, 17 by exonerate): the best blastx hit places it, as in MITOS2
            h = max(locus, key=lambda x: x.bits)
            cover = (h.send - h.sstart + 1) / h.slen
            # ... or a weaker locus that many references support (NC_045937.1 atp8: 47 bits, full length, 50 hits; a
            # false 40-bit atp8 in a flatworm, which has no atp8, had 1 reference)
            supported = cover >= 0.8 and len({x.ref for x in locus}) >= BLASTX_PLACE_REFS
            if not (cover >= BLASTX_PLACE_COVER and h.bits >= BLASTX_PLACE_BITS or supported):
                log.info(f"[core] {gene}: {why}: skipped")
                return None
            tab, tae = (h.start - ws, h.end - ws) if strand == 1 else (we - h.end, we - h.start)
            s, e, partial = place_cds(window, tab, tae, h.sstart - 1, table, gene not in no_down)
            placed = (f"{h.ref} (blastx)", h.sstart - 1, h.send, h.slen, s, e, h.bits, partial)
        ref, qab, qae, reflen, s, e, score, partial = placed
        gs, ge = window_to_genome(ws, we, strand, s, e)
        if gs >= n:                                        # placed in the appended copy of the first ORIGIN_EXT nt
            gs, ge = gs - n, ge - n
        f = Feature("CDS", gene, gs, ge, strand, f"ref {ref} score {score:.0f} residues {qab + 1}-{qae}/{reflen}",
                    partial5=partial, score=score)
        log.info(f"[core] {gene}: {gs + 1}..{ge} ({'+' if strand == 1 else '-'}) ref {ref} score {score:.0f}"
                 + (" start not found" if partial else ""))
        return f

    feats: List[Feature] = []
    all_loci: Dict[str, List[List[Hit]]] = {}
    for gene in PCG:
        hits = [h for h in all_hits[gene] if h.start < n]
        if not hits:
            log.info(f"[core] {gene}: no blastx hit")
            continue
        all_loci[gene] = loci(hits)
        f = place(gene, all_loci[gene][0])
        if f is None and gene not in missing and not any(h.bits > BLASTX_PLACE_BITS for h in hits):
            # weak standard hits that exonerate cannot align (NC_045937.1 atp8: one 35-bit hit, one reference): the
            # sensitive search brings the gene's other references (50 hits)
            log.info(f"[core] {gene}: sensitive blastx (weak hits, not placed)")
            hits = [h for h in blastx_gene(tools["blastx"], genome_fa, refdir / "featureProt" / f"{gene}.fas", code, gene,
                                           2 if threads >= 2 else 1, log, sensitive=True) if h.start < n]
            if hits:
                all_loci[gene] = loci(hits)
                f = place(gene, all_loci[gene][0])
        if f is not None:
            feats.append(f)
    kept = drop_paralogs(feats, n, log)
    # a gene whose best locus lay on another gene (atp8 hit inside atp6 of NC_086642.1, the real atp8 at a weaker locus):
    # its next loci that are not on an annotated gene
    for gene in sorted({f.name for f in feats} - {f.name for f in kept}):
        taken = set().union(*(positions(k, n) for k in kept))
        # hits lying mostly on the placed genes are dropped (one locus may join the real gene and a hit inside its
        # neighbour); the placed genes are masked for exonerate
        off = [h for h in all_hits[gene] if h.start < n
               and len({i % n for i in range(h.start, h.end)} & taken) <= 0.5 * (h.end - h.start)]
        if not off:                                        # the sensitive search, as for a gene without any hit
            log.info(f"[core] {gene}: sensitive blastx for loci off the other genes")
            off = [h for h in blastx_gene(tools["blastx"], genome_fa, refdir / "featureProt" / f"{gene}.fas", code, gene,
                                          1, log, sensitive=True) if h.start < n
                   and len({i % n for i in range(h.start, h.end)} & taken) <= 0.5 * (h.end - h.start)]
        for locus in loci(off)[:MAX_LOCI]:
            f = place(gene, locus, frozenset(taken))
            if f is not None and f in drop_paralogs(kept + [f], n, log):
                log.info(f"[core] {gene}: next locus {f.start + 1}..{f.end} taken")
                kept.append(f)
                break
    single, copies = kept, 0
    for f in [k for k in kept]:
        if f not in kept:                                  # dropped by a second copy placed before it
            continue
        pf = positions(f, n)
        on = [h.bits for loc in all_loci.get(f.name, []) for h in loc if {i % n for i in range(h.start, h.end)} & pf]
        for locus in all_loci.get(f.name, []):
            span = {i % n for h in locus for i in range(h.start, h.end)}
            if not on or span & pf or max(h.bits for h in locus) < COPY_FRAC * max(on):
                continue
            g = place(f.name, locus)
            if g is None or g.partial5 or g.score < COPY_FRAC * f.score or positions(g, n) & pf:
                continue
            new = drop_paralogs(kept + [g], n, log)
            if g in new:
                log.info(f"[core] {f.name}: second copy {g.start + 1}..{g.end} (score {g.score:.0f})")
                kept, copies = new, copies + 1
    if copies == 1:
        new = next(f for f in kept if f not in single)
        orig = next(f for f in single if f.name == new.name)
        ident = copy_identity(genome, orig, new)
        if ident >= COPY_IDENTITY:
            log.info(f"[core] {new.name}: second copy kept, {100 * ident:.1f} % identical to the gene")
        else:
            log.info(f"[core] {new.name}: one gene with a second copy ({100 * ident:.1f} % identical): not a duplicated "
                     "block, copy dropped")
            kept = single
    for gene, (nb, phyla, typical, shortest, longest) in ADJACENT.items():
        nbf = next((f for f in kept if f.name == nb), None)
        if (phylum or "").lower() not in phyla or nbf is None or any(f.name == gene for f in kept):
            continue
        f = orf_before(genome, rc, nbf, gene, table, typical, shortest, longest)
        if f is not None and f in drop_paralogs(kept + [f], n, log):
            log.info(f"[core] {gene}: {f.start + 1}..{f.end} ({'+' if f.strand == 1 else '-'}) {f.note}")
            kept.append(f)
    return sorted(kept, key=lambda f: f.start)


def orf_before(genome: str, rc: str, nb: Feature, gene: str, table: CodonTable.CodonTable, typical: int,
               shortest: int, longest: int) -> Optional[Feature]:
    """ADJACENT: the open reading frame ending with a stop codon 25 nt before to 10 nt after the start of `nb`, same
    strand, start codon giving the length nearest `typical` (shortest..longest nt)."""
    N = len(genome)
    seq, bs = (genome, nb.start) if nb.strand == 1 else (rc, N - nb.end)
    starts, stops = set(table.start_codons), set(table.stop_codons)
    best = None
    for e in range(max(3, bs - 25), min(N, bs + 11)):
        if seq[e - 3:e] not in stops:
            continue
        i = e - 6
        while i >= 0 and e - i <= longest and seq[i:i + 3] not in stops:
            if seq[i:i + 3] in starts and e - i >= shortest and (best is None or abs(e - i - typical) < best[0]):
                best = (abs(e - i - typical), i, e)
            i -= 3
    if best is None:
        return None
    s, e = (best[1], best[2]) if nb.strand == 1 else (N - best[2], N - best[1])
    return Feature("CDS", gene, s, e, nb.strand, f"open reading frame before {nb.name} (gene order)", score=1.0)


# ── sequence ends (mitom --topology auto) ────────────────────────────────────
END_TOUCH = 15             # nt between a CDS and a sequence end, beyond the room its missing residues need: checked
END_FLANK = 30             # nt of the reference alignment on both sides of the junction: it runs across (END_MEET /
                           # END_LONG decide whether the ends are contiguous)
END_CUT = 30               # residues of the reference missing where the alignment stops at the junction, the frame open
                           # up to it, not found across it: the gene is cut there (fewer can be a divergent end)
END_REACH = 60             # ... and the alignment stops at most this far from the junction (exonerate trims the last
                           # codons of a gene cut by a sequence end: 20-21 nt in the 2026-10-30 synthetic tests)
END_MEET = 10              # joined: the reference aligned on each side alone, the two pieces meet within this many residues
                           # (a chance extension across a false junction has none: 5 of 20 synthetic partial genomes ran
                           # 35-338 nt across with one alignment) ...
END_LONG = 90              # ... or >= this many nt aligned on both sides, the frame open across and no piece of one side
                           # re-aligning residues of the other (a divergent half aligns only with the rest: synthetic
                           # NC_045069.1 nad2, 510 / 372 nt). An alignment across that is neither: no evidence either way
ENDS_HEADER = "record\tgene\tstrand\tend\treference\tresidues\tbefore_nt\tafter_nt\tverdict\n"


def side_pieces(window: str, j: int, prots: Dict[str, str], code: int, tools: dict, tmp: Path, log) -> list:
    """The reference aligned on each side of the junction alone (exonerate score >= EXONERATE_MIN_SCORE), as
    (qab, qae, tab, tae) on the whole window; a side without an alignment is left out."""
    out = []
    for lo, hi in ((0, j), (j, len(window))):
        b = exonerate_best(tools["exonerate"], window[lo:hi], prots, code, tmp, log)
        if b is not None:
            out.append((b[1], b[2], b[3] + lo, b[4] + lo))
    return out


FAR_COVER, FAR_SCORE = 0.5, 0.5   # the part a spliced alignment adds across the junction, aligned alone: >= half
#                                   its residues, at >= half the score per residue of the gene's part at the junction


def far_part_aligns(exonerate: str, wide: str, wj: int, prot: str, spliced: Tuple[int, int], near: Tuple[int, int],
                    near_score: float, code: int, tmp: Path, log) -> bool:
    """The residues a spliced alignment (spliced: reference residues sab, sae) adds to the gene's part at the junction
    (near: qab, qae, protein2dna score near_score) align by themselves on the far side of the junction (the window's
    part before wj for missing 5' residues, after it for 3' ones): >= FAR_COVER of them, at >= FAR_SCORE of the near
    part's score per residue. protein2genome over the whole rest of a circle otherwise joins a cut gene to a chance
    match anywhere (contigs, 2026-10-31: atp6 with an 800-nt 'intron', annelid nad4 with an 11-kb one)."""
    per_res = near_score / max(1, near[1] - near[0])
    for a, b, region in ((spliced[0], near[0], wide[:wj]), (near[1], spliced[1], wide[wj:])):
        if b - a < END_CUT:                                   # nothing (or little) added on this side
            continue
        hit = exonerate_best(exonerate, region, {"far": prot[a:b]}, code, tmp, log)
        if hit is None:
            return False
        aligned = hit[2] - hit[1]
        if aligned < FAR_COVER * (b - a) or hit[5] / max(1, aligned) < FAR_SCORE * per_res:
            log.info(f"[core] spliced alignment across the junction rejected: residues {a + 1}-{b} align alone over "
                     f"{aligned} residues at score {hit[5]:.0f} (gene part {per_res:.1f} per residue)")
            return False
    return True


def end_evidence(genome: str, feats: List[Feature], refdir: Path, code: int, tools: dict, tmp: Path, log) -> List[str]:
    """CDS at the sequence ends - across the origin, or nearer an end than its missing residues need (+ END_TOUCH):
    the sequence joined end to start and the CDS's reference protein (the one its note names) aligned again across the
    junction (exonerate, as in annotate_pcgs). One row per CDS (ENDS_HEADER without the record): gene, strand, end
    (across / start / end), reference, residues aligned, nt aligned before / after the junction on the gene's strand,
    and the verdict: joined (>= END_FLANK nt on both sides, and the reference aligned on each side alone reaches the
    junction with pieces that meet, or END_LONG nt on both sides in an open frame: the ends are contiguous), cut (the
    alignment stops within END_REACH nt of the junction with >= END_CUT residues of the reference missing there and
    the reading frame open up to it: they are not) or complete (no evidence). mitom --topology auto decides the
    topology of the record with it."""
    n = len(genome)
    stops = set(CodonTable.unambiguous_dna_by_id[code].stop_codons)
    rows = []
    for f in feats:
        m = re.match(r"ref (\S+)(?: \(blastx\))? score \S+ residues (\d+)-(\d+)/(\d+)", f.note)
        if f.kind != "CDS" or not m:                       # (an ORF placed by gene order has no reference)
            continue
        a, b, reflen = int(m.group(2)), int(m.group(3)), int(m.group(4))
        room5, room3 = 3 * (a - 1) + END_TOUCH, 3 * (reflen - b) + END_TOUCH   # nt the missing residues need
        left, right = (room5, room3) if f.strand == 1 else (room3, room5)      # ... before / after it on the plus strand
        where = "across" if f.end > n else "end" if n - f.end < right else "start" if f.start < left else ""
        if not where:
            continue
        prots = reference_proteins(refdir / "featureProt" / f"{f.name}.fas", [m.group(1)])
        if not prots:
            continue
        s, e = (f.start + n, f.end + n) if where == "start" else (f.start, f.end)
        def window_of(lo: int, hi: int) -> Tuple[str, int]:
            """the joined sequence from lo nt before the gene to hi nt after it, and the junction on it (gene's
            strand); lo + hi <= the rest of the circle (n - gene length), so no position occurs twice"""
            ws, we = s + n - lo, e + n + hi                # on three copies: the junction at 2n
            win = (genome * 3)[ws:we]
            return (win, 2 * n - ws) if f.strand == 1 else (str(Seq(win).reverse_complement()), we - 2 * n)
        rest = max(0, n - (e - s))                         # the circle outside the gene
        flank = min(3 * reflen + WINDOW, rest // 2)
        window, j = window_of(flank, flank)
        best = exonerate_best(tools["exonerate"], window, prots, code, tmp, log)
        if best is None:
            continue
        ref, qab, qae, tab, tae, _ = best
        before, after = j - tab, tae - j
        if before >= END_FLANK and after >= END_FLANK:     # an alignment across the junction
            pieces = side_pieces(window, j, prots, code, tools, tmp, log)
            meet = len(pieces) == 2 and pieces[0][3] >= j - END_REACH and pieces[1][2] <= j + END_REACH \
                and abs(pieces[1][0] - pieces[0][1]) <= END_MEET
            apart = len(pieces) == 2 and abs(pieces[1][0] - pieces[0][1]) > END_MEET
            frame_open = not any(window[c:c + 3] in stops for c in range(tab, tae - 3, 3))
            verdict = "joined" if meet or before >= END_LONG and after >= END_LONG and frame_open and not apart \
                else "complete"
        elif j - END_REACH <= tae < j + END_FLANK and reflen - qae >= END_CUT \
                and not any(window[c:c + 3] in stops for c in range(tae, j - 2, 3)) \
                or j - END_FLANK < tab <= j + END_REACH and qab >= END_CUT \
                and not any(window[c:c + 3] in stops for c in range(tab - 3, j - 1, -3)):
            # the gene's 3' side runs into the junction, or its 5' side starts there, the frame open up to it: cut -
            # unless the missing residues lie across the junction behind an intron (protein2dna has no introns; a
            # record started inside a gene with an intron: Travisia cox1, every MJ20260820041 CHY record). The intron
            # may hold other genes and be many kb long (hexacoral nad5: 11-20 kb) and lies wholly on one side of the
            # gene, so the spliced alignment gets the whole rest of the circle on one side, then on the other: no
            # intron length limit, no assumption about the side.
            verdict = "cut"
            for lo, hi in ((rest, 0), (0, rest)):
                wide, wj = window_of(lo, hi)
                spl = exonerate_best(tools["exonerate"], wide, {ref: prots[ref]}, code, tmp, log, "protein2genome")
                if spl is None:
                    continue
                _, sab, sae, stab, stae, _ = spl
                gained = (qab - sab) + (sae - qae)          # residues aligned beyond the protein2dna alignment
                if stab < wj - END_FLANK and stae > wj + END_FLANK and gained >= END_CUT \
                        and far_part_aligns(tools["exonerate"], wide, wj, prots[ref], (sab, sae), (qab, qae), best[5],
                                            code, tmp, log):
                    qab, qae = sab, sae
                    before, after = wj - stab, stae - wj
                    verdict = "joined_intron"
                    break
        else:
            verdict = "complete"
        rows.append(f"{f.name}\t{'+' if f.strand == 1 else '-'}\t{where}\t{ref}\t{qab + 1}-{qae}/{reflen}\t"
                    f"{before}\t{after}\t{verdict}\n")
        log.info(f"[core] {f.name} at the sequence {where}: reference aligned {before} nt before / {after} nt after "
                 f"the junction, residues {qab + 1}-{qae}/{reflen}: {verdict}")
    return rows


# ── output ───────────────────────────────────────────────────────────────────
def location(f: Feature, n: int):
    """f's location; a feature ending beyond n (circular record) is a join across the origin in transcript order;
    a 5'-partial CDS gets < (plus strand) or > (minus strand) at its 5' end."""
    s = BeforePosition(f.start) if f.partial5 and f.strand == 1 else f.start
    if f.end <= n:
        return FeatureLocation(s, AfterPosition(f.end) if f.partial5 and f.strand == -1 else f.end, f.strand)
    tail = f.end - n
    first = FeatureLocation(s, n, f.strand)
    last = FeatureLocation(0, AfterPosition(tail) if f.partial5 and f.strand == -1 else tail, f.strand)
    return CompoundLocation([first, last] if f.strand == 1 else [last, first])


def write_genbank(genome: str, feats: List[Feature], code: int, prefix: str, topology: str, out: Path) -> Path:
    from Bio import SeqIO
    SeqIO.write(genbank_record(genome, feats, code, prefix, topology), str(out), "genbank")
    return out


def genbank_record(genome: str, feats: List[Feature], code: int, prefix: str, topology: str) -> SeqRecord:
    rec = SeqRecord(Seq(genome), id=prefix, name=prefix[:16], description="")
    rec.annotations.update(molecule_type="DNA", topology=topology, data_file_division="INV")
    rec.features.append(SeqFeature(FeatureLocation(0, len(genome), 1), type="source",
                                   qualifiers={"organism": [prefix], "organelle": ["mitochondrion"], "mol_type": ["genomic DNA"]}))
    for f in sorted(feats, key=lambda x: x.start):
        if f.kind == "CDS":
            gene, product = NAMES[f.name]
            loc = location(f, len(genome))
            q = {"gene": [gene], "product": [product], "codon_start": ["1"], "transl_table": [str(code)],
                 "note": [f.note]}
            if f.partial5:
                q["note"].append("start codon not determined")
        elif f.kind == "rRNA":
            gene, product = NAMES[f.name]
            loc = location(f, len(genome))
            q = {"gene": [gene], "product": [product], "note": [f.note]}
        else:
            gene = f.name
            loc = location(f, len(genome))
            aa3 = AA3.get(f.name[3], "Xaa")
            q = {"gene": [gene], "product": [f"tRNA-{aa3}"], "note": [f.note]}
            if f.anticodon:                                # as MITOS2 output is written: mito_curate reads this note
                q["note"] = [f"anticodon: {f.anticodon.upper().replace('T', 'U')}", f.note]
        rec.features.append(SeqFeature(loc, type="gene", qualifiers={"gene": [gene]}))
        rec.features.append(SeqFeature(loc, type=f.kind, qualifiers=q))
    return rec


def annotate(args) -> int:
    """Every sequence of the FASTA (contigs of one sample: one GenBank file with one record each)."""
    global EXCLUDE_REFS
    log = logging.getLogger("core")
    t0 = time.time()
    fa = Path(args.input).resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    seqs = list(read_fasta(fa))
    if not seqs:
        raise SystemExit(f"[core] no sequence in {fa}")
    prefix = args.prefix or re.sub(r"[^A-Za-z0-9_.-]", "_", seqs[0].id)[:16]
    refdir = Path(args.refdir or mito_env.mitos_refdir(None) or "").expanduser()
    if (refdir / "refseq89m").is_dir():                     # mito_env gives the directory that holds refseq89m
        refdir = refdir / "refseq89m"
    if not (refdir / "featureProt").is_dir() or not (refdir / "ncRNA").is_dir():
        raise SystemExit(f"[core] reference data not found: {refdir} (refseq89m with featureProt/ and ncRNA/)")
    tools = {"blastx": program("blastx", args.blastx), "exonerate": program("exonerate", args.exonerate),
             "cmsearch": program("cmsearch", args.cmsearch)}
    lineage = args.lineage or (f"Metazoa;{args.phylum}" if args.phylum else "")
    use_ml = False
    if not args.no_trna and args.trna in ("auto", "ml"):
        for t in ("cmalign", "arwen", "aragorn"):
            tools[t] = mito_env.find(t, getattr(args, t, None)) or ""
        tools["refdir"] = str(refdir)
        use_ml = mito_trna.mitom_ready(tools, str(refdir))
        if args.trna == "ml" and not use_ml:
            raise SystemExit("[core] --trna ml: MitoM-ML needs models/trna_ml, cmsearch and cmalign")
        if not use_ml:
            log.warning("[core] MitoM-ML not installed (models/trna_ml, cmalign, scikit-learn): tRNAs by MiTFi")
    if not args.no_trna and not use_ml:
        tools["java"] = program("java", args.java)
        tools["mitfi"] = mitfi_files(args.mitos_bin)
    elif not args.no_trna:                                   # MiTFi for the families MitoM-ML does not call
        try:
            tools["java"], tools["mitfi"] = program("java", args.java), mitfi_files(args.mitos_bin)
        except (SystemExit, OSError) as e:
            log.warning(f"[core] MiTFi not found ({e}): tRNA families MitoM-ML does not call are not searched again")
    code = int(args.genetic_code)
    EXCLUDE_REFS = {x.split(".")[0] for x in re.split(r"[,\s]+", args.exclude_ref or os.environ.get("MITOM_EXCLUDE_REF", "")) if x}
    if EXCLUDE_REFS:
        log.info(f"[core] reference proteins of {', '.join(sorted(EXCLUDE_REFS))} not used")
    recs, rows, ends, n_all = [], [], [], {"CDS": 0, "rRNA": 0, "tRNA": 0}
    for k, rec in enumerate(seqs):
        genome = str(rec.seq).upper()
        n = len(genome)
        topology = topology_from_header(rec.description)
        rid = contig_record_id(prefix, rec.id, k, len(seqs))
        searched = genome + genome[:min(ORIGIN_EXT, n)] if topology == "circular" else genome
        with tempfile.TemporaryDirectory(prefix="mito_core.") as td:
            tmp = Path(td)
            gfa = tmp / f"{rid}.fa"
            gfa.write_text(f">{rid}\n{searched}\n")
            feats = annotate_pcgs(searched, gfa, refdir, code, args.threads, tools, tmp, log, n, args.phylum)
            feats += annotate_rrnas(gfa, refdir, tools, tmp, log, args.phylum,
                                    Path(args.rrna_lib) if args.rrna_lib else None, n if topology == "circular" else None,
                                    args.rrna_method)
            if not args.no_trna:
                feats += mito_trna.core_trnas(genome, searched, gfa, refdir, code, args.threads, tools, tmp, log, n,
                                              args.phylum, lineage, topology == "circular", use_ml, feats)
            if args.check_ends:
                ends += [f"{rid}\t{row}" for row in end_evidence(genome, feats, refdir, code, tools, tmp, log)]
            if args.keep_tmp:
                shutil.copytree(tmp, out / "core_tmp" / rid, dirs_exist_ok=True)
        recs.append(genbank_record(genome, feats, code, rid, topology))
        for f in sorted(feats, key=lambda x: x.start):
            rows.append(f"{rid}\t{f.kind}\t{f.name}\t{f.start + 1}\t{f.end if f.end <= n else f.end - n}\t"
                        f"{'+' if f.strand == 1 else '-'}\t{f.score:.1f}\t{f.anticodon}\t{f.note}\n")
        for kind in n_all:
            n_all[kind] += sum(1 for f in feats if f.kind == kind)
    gb = out / f"{prefix}.gbf"
    from Bio import SeqIO
    SeqIO.write(recs, str(gb), "genbank")
    with (out / f"{prefix}.core.tsv").open("w") as fh:
        fh.write("record\tkind\tname\tstart\tend\tstrand\tscore\tanticodon\tnote\n")
        fh.writelines(rows)
    if args.check_ends:
        (out / f"{prefix}.ends.tsv").write_text(ENDS_HEADER + "".join(ends))
    log.info(f"[core] done in {time.time() - t0:.0f} s: {n_all['CDS']} protein genes, {n_all['rRNA']} rRNAs, "
             f"{n_all['tRNA']} tRNAs in {len(recs)} sequence(s) -> {gb}")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    global EXONERATE_MODEL
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("annotate", help="annotate one mitogenome FASTA")
    a.add_argument("-i", "--input", required=True, help="FASTA (one mitogenome; several sequences = contigs)")
    a.add_argument("-o", "--out", required=True, help="output directory")
    a.add_argument("--prefix", default="", help="record / file name (default: the FASTA id)")
    a.add_argument("--genetic_code", type=int, default=5)
    a.add_argument("--threads", type=int, default=2)
    a.add_argument("--refdir", default=None, help="MITOS2 reference data (refseq89m); default from mito_env")
    a.add_argument("--no_trna", action="store_true", help="skip the tRNA step")
    a.add_argument("--trna", choices=["auto", "ml", "mitfi"], default="auto",
                   help="tRNA caller: auto (default) = MitoM-ML when installed, else MiTFi; ml; mitfi (MITOS2's)")
    a.add_argument("--lineage", default="", help="NCBI lineage of the record (';'-separated) for the MitoM-ML clade "
                   "models; default Metazoa;<--phylum>")
    a.add_argument("--phylum", default=None, help="phylum of the record: selects the per-phylum rRNA profile HMM")
    a.add_argument("--rrna_lib", default=None, help="rRNA profile HMM library (default models/rrna_hmm)")
    a.add_argument("--rrna_method", choices=["hmm", "cm"], default="hmm",
                   help="12S / 16S search: hmm (default) = MitoM's per-phylum profile HMMs (faster), MITOS2's covariance "
                        "models where no model or hit; cm = MITOS2's covariance models only (same accuracy)")
    a.add_argument("--check_ends", action="store_true",
                   help="align the CDS at the sequence ends again across them: <out>/<prefix>.ends.tsv tells whether "
                        "the ends are contiguous (mitom --topology auto)")
    a.add_argument("--keep_tmp", action="store_true", help="keep blastx / cmsearch / MiTFi files in <out>/core_tmp")
    a.add_argument("--exonerate_model", default=EXONERATE_MODEL, help="exonerate model (protein2dna, protein2dna:bestfit)")
    a.add_argument("--blastx"), a.add_argument("--exonerate"), a.add_argument("--cmsearch"), a.add_argument("--java")
    a.add_argument("--cmalign"), a.add_argument("--arwen"), a.add_argument("--aragorn")
    a.add_argument("--mitos_bin", default=None, help="MITOS2 environment bin (to find MiTFi); default from mito_env")
    a.add_argument("--exclude_ref", default=None, help="accessions whose own reference proteins are not used "
                   "(benchmark leave-out; default $MITOM_EXCLUDE_REF)")
    a.add_argument("--log", default=None)
    args = p.parse_args(argv)
    EXONERATE_MODEL = args.exonerate_model
    logging.basicConfig(level=logging.DEBUG if debug_mode() else logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stderr)] + ([logging.FileHandler(args.log)] if args.log else []))
    return annotate(args)


if __name__ == "__main__":
    sys.exit(main())
