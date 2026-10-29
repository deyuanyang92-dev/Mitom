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
hit) and MiTFi (mitfi_import.jar of the mitos package) for hit selection and anticodons - the programs MITOS2 runs.

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
from typing import Dict, Iterator, List, Optional, Tuple

from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mito_env  # noqa: E402
import mito_rrna_hmm  # noqa: E402  (per-phylum rRNA profile HMMs: nhmmer search, best span)

PCG = ("atp6", "atp8", "cob", "cox1", "cox2", "cox3", "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6")
NAMES = {                                   # the core's /gene and /product, understood by every later MitoM step
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
CM_EVALUE = 0.01           # cmsearch fast pass and MiTFi threshold of that pass (MITOS2 runmitos default --ncev)
NO_SENSITIVE_TRNA = {"cnidaria"}   # phyla whose missing tRNA models are not searched again with cmsearch --max:
                           # cnidarian mitogenomes encode 1-6 tRNAs, the ~20-model --max pass found junk (18 records,
                           # leave-out: 13 of 17 extra tRNAs, 2 real ones; 344 -> 50 s median; a false trnE inside nad5
                           # was an NCBI error in the stress test)
CM_EVALUE_SENSITIVE = 0.1  # cmsearch --max pass for models without a hit (MITOS2 sensitive pass)
TRNA_NAMES = ["trn" + f for f in "A C D E F G H I K L1 L2 M N P Q R S1 S2 T V W Y".split()]
TRNA_RESCUE_OVL = 10       # a MiTFi tRNA of a family MitoM-ML did not call is added when it overlaps no call by more
                           # than this many nt (rescue_trnas)
TRNA_SURE_P = 0.9          # MitoM-ML calls below this p need MiTFi's call of the same tRNA, and a weaker second call
                           # of a family needs at least this p (trna_second_opinion; tuning halves: 37 % of the calls
                           # below 0.9 without MiTFi were wrong, 2 % above, 30 % of the weaker second copies below)
EXONERATE_MODEL = "protein2dna"   # protein2dna:bestfit = whole reference protein aligned (set by --exonerate_model)
EXONERATE_CODES = {1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23}   # built into exonerate 2.4.0; others
                           # (24, 33: Hemichordata) are passed as the 64-letter table ("No built in genetic code ... [24]")
ORIGIN_EXT = 3000          # circular records are searched as seq + seq[:ORIGIN_EXT]: a gene across the origin is
                           # contiguous there (MITOS2 searches two rotations; review 2026-10-22: cox1 cut 60/40 by the
                           # origin lost 614 nt, trnF / rrnL across it were lost or truncated)
PARALOG_MIN_OVERLAP = 60   # a CDS overlapping a better-scoring CDS of another gene by more than this and by >= half of
                           # the shorter one is a paralog hit (atp8 inside atp6, nad2 inside nad5 when the real gene is
                           # absent; MitoZ genomic_cluster.pl / MITOS2 _nonoverlapping_greedy drop these too)
COPY_EVALUE = 1e-3         # MiTFi copies (second hits of a model) only this good (MITOS2 runmitos.py)
MITFI_TRNA_OVL = 50        # tRNAs may overlap this much in MiTFi: runmitos --maxtrnaovl default (its function default 10 is
                           # not what the command runs); with 10, NC_060484.1 lost trnL2 overlapping trnP by 12 nt (MITOS2
                           # keeps both); MiTFi -overlap 50 -method l reproduces MITOS2's own tRNAout.nc line for line
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


def read_fasta(path: Path) -> Iterator[SeqRecord]:
    """FASTA records (id = first word, description = whole header), streamed, without Bio.SeqIO: it imports numpy
    (~27 MB), resident in this process while nhmmer / blastx / cmsearch run (peak memory); SeqIO only writes GenBank"""
    head, parts = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if head is not None:
                    yield SeqRecord(Seq("".join(parts)), id=(head.split() or [""])[0], description=head)
                head, parts = line[1:].strip(), []
            elif head is not None:
                parts.append("".join(line.split()))
    if head is not None:
        yield SeqRecord(Seq("".join(parts)), id=(head.split() or [""])[0], description=head)


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


def exonerate_best(exonerate: str, window: str, prots: Dict[str, str], code: int, tmp: Path, log) \
        -> Optional[Tuple[str, int, int, int, int, float]]:
    """Best protein2dna alignment of the reference proteins on the window (plus strand):
    (ref, qab, qae, tab, tae, score); coordinates 0-based half-open."""
    (tmp / "w.fa").write_text(f">w\n{window}\n")
    best = None
    for rid, seq in prots.items():
        (tmp / "q.fa").write_text(f">{rid}\n{seq}\n")
        try:
            p = run([exonerate, "--model", EXONERATE_MODEL, "--geneticcode", exonerate_code(code), "--bestn", "1", "--showalignment", "no",
                 "--showvulgar", "no", "--score", str(EXONERATE_MIN_SCORE),
                 "--ryo", "RYO\t%qi\t%ql\t%qab\t%qae\t%tab\t%tae\t%s\t%tS\n",
                 *(["--exhaustive", "yes"] if EXONERATE_MODEL.endswith("bestfit") else []),   # non-local models need full DP
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


def positions(f: Feature, n: int) -> set:
    return {i % n for i in range(f.start, f.end)}


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


# ── RNAs ─────────────────────────────────────────────────────────────────────
def cmsearch(cm: str, model: Path, genome_fa: Path, out: Path, sensitive: bool, log, rrna: bool = False,
             z: Optional[float] = None) -> Path:
    """cmsearch as MITOS2 runs it (local mode, --cpu 1 --notextw --smxsize 80000); fast pass --noF4b -E 0.01,
    sensitive pass --max -E 0.1 (rRNA: --nohmm). Local mode finds the truncated tRNAs of nematodes. z = search space
    (Mb) for the E-values: a circular record is searched with ORIGIN_EXT nt appended, MITOS2 searches it at its length."""
    cmd = [cm, "--cpu", "1", "--notextw", "--smxsize", "80000", "--tblout", str(out.with_suffix(".tbl"))]
    if z:
        cmd += ["-Z", f"{z:.6f}"]
    if rrna:                                   # same rrnL hit, 536 -> 463 MB, 3.8 -> 2.7 s (alignment matrix cap)
        cmd += ["--mxsize", "128"]
    if sensitive:
        cmd += ["--nohmm" if rrna else "--max", "-E", str(CM_EVALUE_SENSITIVE)]
    else:
        cmd += ["--noF4b", "-E", str(CM_EVALUE)]
    p = run(cmd + [str(model), str(genome_fa)], log)
    out.write_text(p.stdout)
    return out


def tbl_hits(tbl: Path) -> List["mito_rrna_hmm.Hit"]:
    """Hits of a cmsearch --tblout file as mito_rrna_hmm.Hit (0-based half-open plus coordinates, model columns,
    bits), so that the rRNA HMM rules join the pieces of one gene (consecutive, disjoint model columns)."""
    hits = []
    for line in tbl.read_text().splitlines():
        if line.startswith("#"):
            continue
        c = line.split()
        a, b = int(c[7]), int(c[8])
        hits.append(mito_rrna_hmm.Hit(min(a, b) - 1, max(a, b), 1 if c[9] == "+" else -1, int(c[5]), int(c[6]),
                                      float(c[14])))
    return hits


def cm_len(cm: Path) -> int:
    """CLEN (model columns) of an Infernal .cm file."""
    with open(cm) as fh:
        for line in fh:
            if line.startswith("CLEN"):
                return int(line.split()[1])
    return 0


def rrna_hmm(lib: Path, phylum: Optional[str], model: str) -> Optional[Path]:
    """The per-phylum profile HMM of model (rrnS / rrnL) in lib for `phylum`, None when there is none."""
    if not phylum or not lib.is_dir():
        return None
    name = {p.name.split(".")[0].lower(): p.name.split(".")[0] for p in lib.glob("*.hmm")}.get(phylum.lower())
    p = lib / f"{name}.{model}.hmm" if name else None
    return p if p and p.is_file() else None


def annotate_rrnas(genome_fa: Path, refdir: Path, tools: dict, tmp: Path, log, phylum: Optional[str] = None,
                   rrna_lib: Optional[Path] = None, n: Optional[int] = None, method: str = "hmm") -> List[Feature]:
    """12S / 16S: method hmm (default) = the per-phylum profile HMM (nhmmer, ~0.1 GB, seconds) when the library has
    the record's phylum, else the MITOS2 covariance model (cmsearch, 0.4-0.6 GB, 10 s); method cm = the MITOS2
    covariance model only. Both give the same accuracy (2026-10-28, full pipeline on 300 + 150 leave-out records:
    found 98.5 / 98.5 % and 99.3 / 99.3 %, both ends exact 47.4 / 47.6 % and 32.3 / 32.3 %); R1 of mito_curate.py
    then fits the ends to the flanking genes. Hits are joined by
    mito_rrna_hmm's rules (n = record length of a circular record searched with ORIGIN_EXT nt appended); a span over
    less than MIN_COVER of the model is not taken (MitoM's rRNA HMM step looks again). A second copy (best_spans'
    runner-up over the same model columns elsewhere) scoring >= COPY_FRAC of the gene is added too: duplicated rRNAs,
    both copies in NCBI (NC_044480.1 rrnS + rrnL, mermithid nematodes NC_008046.1, NC_008047.1; nhmmer with the
    leave-out libraries on the 1,050 + 669 benchmark records found 7 such copies, all of them NCBI's second copy)."""
    feats = []
    lib = rrna_lib or mito_rrna_hmm.LIB
    for model in ("rrnS", "rrnL"):
        hmm = rrna_hmm(lib, phylum, model) if method == "hmm" else None
        cm = refdir / "ncRNA" / f"{model}.cm"
        span, copy, how = None, None, ""
        if hmm is not None:
            span, copy = mito_rrna_hmm.best_spans(mito_rrna_hmm.nhmmer(hmm, genome_fa), n)
            leng, how = mito_rrna_hmm.model_len(hmm), f"nhmmer {hmm.name}"
            if span is None:
                log.info(f"[core] {model}: no nhmmer hit with {hmm.name}; trying the covariance model")
            elif (span.hmm_to - span.hmm_from + 1) / leng < mito_rrna_hmm.MIN_COVER:
                # a weak partial hit of a divergent genome (tunicates against the vertebrate-dominated Chordata
                # model: NC_009833.1, NC_013561.1, whose rRNAs MITOS2's covariance models find)
                log.info(f"[core] {model}: nhmmer span covers {(span.hmm_to - span.hmm_from + 1) / leng:.2f} of "
                         f"{hmm.name}; trying the covariance model")
                span = None
        if span is None:
            hits = tbl_hits(cmsearch(tools["cmsearch"], cm, genome_fa, tmp / f"{model}.out", False, log, rrna=True).with_suffix(".tbl"))
            span, copy = mito_rrna_hmm.best_spans(hits, n)
            leng, how = cm_len(cm), "cmsearch"
        if span is None:
            # no sensitive (--nohmm) pass: 2.4 GB and 2 min per model; a missing rRNA is found later by MitoM's
            # per-phylum rRNA HMM step (mito_rrna_hmm.py, nhmmer)
            log.info(f"[core] {model}: no hit (left to the rRNA HMM step)")
            continue
        cover = (span.hmm_to - span.hmm_from + 1) / leng if leng else 1.0
        if cover < mito_rrna_hmm.MIN_COVER:
            log.info(f"[core] {model}: best span {span.start + 1}..{span.end} covers {cover:.2f} of the model: "
                     "not taken (left to the rRNA HMM step)")
            continue
        feats.append(Feature("rRNA", model, span.start, span.end, span.strand,
                             f"{how} bits {span.score:.0f} model {span.hmm_from}-{span.hmm_to}/{leng}", score=span.score))
        log.info(f"[core] {model}: {span.start + 1}..{span.end} ({'+' if span.strand == 1 else '-'}) by {how}")
        if copy is not None and copy.score >= COPY_FRAC * span.score and \
                (copy.hmm_to - copy.hmm_from + 1) / leng >= mito_rrna_hmm.MIN_COVER:
            feats.append(Feature("rRNA", model, copy.start, copy.end, copy.strand,
                                 f"{how} bits {copy.score:.0f} model {copy.hmm_from}-{copy.hmm_to}/{leng}", score=copy.score))
            log.info(f"[core] {model}: second copy {copy.start + 1}..{copy.end} ({'+' if copy.strand == 1 else '-'}) by {how}")
    return feats


def annotate_trnas(genome: str, genome_fa: Path, refdir: Path, code: int, threads: int, tools: dict, tmp: Path, log,
                   n: Optional[int] = None, phylum: Optional[str] = None) -> List[Feature]:
    """MITOS2's tRNA step: cmsearch per MiTFi model (fast pass, then --max for the tRNA names without a hit), MiTFi
    for selection, structure and anticodon. genome / genome_fa = the searched sequence, n = record length (circular
    records carry ORIGIN_EXT appended nt; positions are mapped back modulo n)."""
    n = n or len(genome)
    models = sorted((refdir / "ncRNA").glob("trn*.cm"))
    seqname = genome_fa.name

    def search(model: Path, sensitive: bool) -> Path:
        # E-values for the record's length on both strands, as MITOS2 (which searches rotations of the record, not the
        # record + ORIGIN_EXT): the 3,000 appended nt made every E ~1.2x larger on a 15-kb record, so hits MITOS2 keeps
        # at E 0.08-0.1 were lost (NC_060482.1 / NC_060708.1 acanthocephalan trnD, E 0.096 -> 0.12)
        return cmsearch(tools["cmsearch"], model, genome_fa, tmp / f"{seqname}_{model.stem}.cmout", sensitive, log,
                        z=2 * n / 1e6)

    jar, codes = tools["mitfi"]

    def mitfi(outs: List[Path], evalue: float) -> List[Feature]:
        cmd = [tools["java"], "-jar", str(jar), "-code", str(code), "-codefile", str(codes), "-structure", "-genome",
               str(genome_fa), "-glength", str(n), "-evalue", str(evalue), "-method", "l", "-overlap",
               str(MITFI_TRNA_OVL)]
        for o in outs:
            cmd += ["-in", str(o)]
        return parse_mitfi(run(cmd, log).stdout, log, n)

    # MITOS2's procedure (runmitos.py): fast pass judged at E <= 0.001 - a model whose only hits are weaker counts
    # as empty and is searched again with --max; the final MiTFi run accepts E <= 0.1
    with ThreadPoolExecutor(max(1, threads)) as ex:
        outs = list(ex.map(lambda m: search(m, False), models))
    feats = mitfi(outs, CM_EVALUE)
    accepted = {f.name for f in feats}          # MITOS2 re-runs the models whose NAME is missing (MiTFi names by anticodon:
    missing = [m for m in models if m.stem not in accepted]   # an L1 gene is often found by trnL2.cm)
    if (phylum or "").lower() in NO_SENSITIVE_TRNA:
        log.info(f"[core] tRNA sensitive pass skipped for {phylum}")
        missing = []
    if missing:
        log.info(f"[core] tRNA sensitive pass for {len(missing)} model(s): " + " ".join(m.stem for m in missing))
        with ThreadPoolExecutor(max(1, threads)) as ex:
            list(ex.map(lambda m: search(m, True), missing))
        feats = mitfi(outs, CM_EVALUE_SENSITIVE)
    log.info(f"[core] tRNAs: {len(feats)}")
    return feats


def annotate_trnas_ml(genome: str, n: int, code: int, threads: int, tools: dict, lineage: str, circular: bool,
                      log) -> List[Feature]:
    """MitoM-ML tRNAs of one record (genome = the record without the appended ORIGIN_EXT nt)."""
    import mito_trna
    import trna_ml
    ext: list = []
    with tempfile.TemporaryDirectory(prefix="mito_core_trna.") as td:
        fa = Path(td) / "g.fa"
        fa.write_text(f">g\n{genome}\n")
        for t in ("arwen", "aragorn"):
            if tools.get(t):
                ext += mito_trna.unwrap(mito_trna.TOOLS[t](fa, code, {t: tools[t]}, Path(td), None), n)
    calls = trna_ml.predict(genome, lineage, code, Path(tools["refdir"]) / "ncRNA", tools["cmsearch"], tools["cmalign"], ext,
                            cpu=max(1, threads))
    feats = []
    for c in calls:
        if c["end"] > n and not circular:
            continue
        feats.append(Feature("tRNA", trna_name(c["family"], c["anticodon"], ""), c["start"] - 1, c["end"], c["strand"],
                             f"MitoM-ML p {c['p']:.3f} ends p {c['bound_p']:.3f}", anticodon=c["anticodon"],
                             score=c["p"], model="MitoM-ML"))
    log.info(f"[core] tRNAs (MitoM-ML, lineage {lineage or 'none'}, {'ARWEN/ARAGORN evidence' if ext else 'no ARWEN/ARAGORN'}):"
             f" {len(feats)}")
    return feats


def rescue_trnas(ml: List[Feature], mitfi: List[Feature], n: int) -> List[Feature]:
    """MitoM-ML calls plus MiTFi's calls of the families MitoM-ML did not call (lowest E first, one per family), each
    when it overlaps no call by more than TRNA_RESCUE_OVL nt (either strand, across the origin). The learned caller
    rejects weak hits in genomes with few strong tRNAs, where MiTFi still finds some (unseen-family benchmark,
    benchmark/trna_ml/scripts/rescue_sim.py: F1 98.41 -> 98.66 % on the tuning half, 98.77 -> 98.83 % on the other
    half; Acanthocephala recall 17 -> 43 %; no cnidarian record changed)."""
    def evalue(f: Feature) -> float:
        m = re.search(r" E ([0-9.eE+-]+)", f.note)
        return float(m.group(1)) if m else 1.0

    def overlap(a: Feature, b: Feature) -> int:
        return max(min(a.end, b.end + k) - max(a.start, b.start + k) for k in (-n, 0, n))
    out, have = list(ml), {f.name for f in ml}
    for f in sorted(mitfi, key=evalue):
        if f.name not in have and all(overlap(f, g) <= TRNA_RESCUE_OVL for g in out):
            out.append(f)
            have.add(f.name)
    return out


def same_trna(a: Feature, b: Feature, n: int) -> bool:
    """Same tRNA name and strand, overlapping >= 50 % of the shorter (across the origin too)."""
    o = max(min(a.end, b.end + k) - max(a.start, b.start + k) for k in (-n, 0, n))
    return a.name == b.name and a.strand == b.strand and o >= 0.5 * min(a.end - a.start, b.end - b.start)


def trna_second_opinion(ml: List[Feature], mitfi: Optional[List[Feature]], n: int) -> List[Feature]:
    """MitoM-ML calls kept when sure enough: a weaker second call of a family needs p >= TRNA_SURE_P, and any call with
    p < TRNA_SURE_P needs MiTFi's call of the same tRNA (not applied without MiTFi). Chosen on the tuning halves of
    the unseen-family set and the GenBank-only annelids (precision >= MITOS2 + 0.1 points, most recall), then
    checked on the held-out halves (benchmark/trna_ml/README.md)."""
    best: Dict[str, float] = {}
    for f in ml:
        best[f.name] = max(best.get(f.name, 0.0), f.score)
    out = [f for f in ml if f.score >= best[f.name] or f.score >= TRNA_SURE_P]
    if mitfi is not None:
        out = [f for f in out if f.score >= TRNA_SURE_P or any(same_trna(f, g, n) for g in mitfi)]
    return out


def mitfi_ends(ml: List[Feature], mitfi: List[Feature], n: int) -> List[Feature]:
    """MitoM-ML calls with MiTFi's ends wherever MiTFi calls the same tRNA (name and strand, >= 50 % of the shorter).
    For vertebrates: the learned ends never saw one and are 1 nt off NCBI more often than MiTFi's (434 large-set
    chordates: both ends exact 90.2 -> 91.4 %, MITOS2 91.1 %; 188 held-out vertebrates 91.3 -> 92.7 %)."""
    out = []
    for f in ml:
        m = next((g for g in mitfi if same_trna(f, g, n)), None)
        out.append(f if m is None else Feature(f.kind, f.name, m.start, m.end, f.strand, f.note + ", ends by MiTFi",
                                               anticodon=f.anticodon, score=f.score, model=f.model))
    return out


def vertebrate(lineage: str) -> bool:
    """Lineage of a vertebrate; a bare 'Metazoa;Chordata' (FASTA input, --phylum only) counts as one."""
    tokens = {t.strip() for t in lineage.split(";")}
    return bool(tokens & {"Vertebrata", "Craniata"}) or lineage.strip() == "Metazoa;Chordata"


def trna_name(aa: str, anticodon: str, model: str) -> str:
    """MITOS2 naming: the amino acid MiTFi derives from the anticodon; the two Leu / Ser genes by codon family
    (L1 = CUN: anticodon UAG; L2 = UUR: UAA; S1 = AGN: GCU / UCU; S2 = UCN: UGA), else the model's number."""
    aa = (aa or "").strip().upper()
    if re.fullmatch(r"[LS][12]", aa):                  # MiTFi names the two Leu / Ser genes itself (codon family)
        return "trn" + aa
    aa = aa[:1]
    ac = anticodon.replace("u", "t")
    if aa == "L":
        return "trnL1" if ac == "tag" else "trnL2" if ac == "taa" else ("trnL" + (re.search(r"trnL([12])", model) or [None, ""])[1])
    if aa == "S":
        return "trnS1" if ac in ("gct", "tct") else "trnS2" if ac == "tga" else ("trnS" + (re.search(r"trnS([12])", model) or [None, ""])[1])
    if aa and aa in AA3:
        return "trn" + aa
    m = re.search(r"(trn[A-Z](?:[12])?)", model)
    return m.group(1) if m else ""


def parse_mitfi(text: str, log, n: Optional[int] = None) -> List[Feature]:
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
            m_aa = re.match(r"trn(\w\d*)", aa)       # no anticodon found: MiTFi writes the model name here ("trnL1";
            aa = m_aa.group(1) if m_aa else aa        # MITOS2 mitfi parser: "workaround"); [:1] made it Thr (NC_058285.1)
            tophit, copy = int(c[15]), int(c[16])
        except (IndexError, ValueError) as exc:
            log.warning(f"[core] MiTFi line not understood ({exc}): {line[:80]}")
            continue
        if n and start >= n:                     # found again in the appended first ORIGIN_EXT nt of a circular record
            start, stop = start - n, stop - n
        if not (tophit or (copy and evalue <= COPY_EVALUE)):  # MITOS2: a copy (second hit of a model) only when E <= 0.001
            continue
        name = trna_name(aa, anticodon, model)
        if not name or start >= stop:
            continue
        mm = re.search(r"(trn[A-Z](?:[12])?)\.cm", model)
        feats.append(Feature("tRNA", name, start, stop, strand, f"MiTFi {model} bits {bits:.1f} E {evalue:.1e}"
                             + (" copy" if copy and not tophit else ""), anticodon=anticodon, acpos=acpos, score=bits,
                             model=mm.group(1) if mm else model))
    # tRNAs overlapping by more than MiTFi allows (the same tRNA found at i and at i + n in a circular record: complete
    # overlap): the better score is kept
    feats.sort(key=lambda f: -f.score)
    kept: List[Feature] = []
    L = n or max((f.end for f in feats), default=1)
    for f in feats:
        if all(len(positions(f, L) & positions(k, L)) <= MITFI_TRNA_OVL for k in kept):
            kept.append(f)
    return sorted(kept, key=lambda f: f.start)


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
        import mito_trna
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
    recs, rows, n_all = [], [], {"CDS": 0, "rRNA": 0, "tRNA": 0}
    for k, rec in enumerate(seqs):
        genome = str(rec.seq).upper()
        n = len(genome)
        topology = "linear" if re.search(r"topology\s*=\s*linear", rec.description, re.I) else "circular"
        rid = prefix if len(seqs) == 1 else (re.sub(r"[^A-Za-z0-9_.-]", "_", rec.id)[:16] or f"{prefix}_{k + 1}")
        searched = genome + genome[:min(ORIGIN_EXT, n)] if topology == "circular" else genome
        with tempfile.TemporaryDirectory(prefix="mito_core.") as td:
            tmp = Path(td)
            gfa = tmp / f"{rid}.fa"
            gfa.write_text(f">{rid}\n{searched}\n")
            feats = annotate_pcgs(searched, gfa, refdir, code, args.threads, tools, tmp, log, n, args.phylum)
            feats += annotate_rrnas(gfa, refdir, tools, tmp, log, args.phylum,
                                    Path(args.rrna_lib) if args.rrna_lib else None, n if topology == "circular" else None,
                                    args.rrna_method)
            if not args.no_trna and use_ml:
                trnas = annotate_trnas_ml(genome, n, code, args.threads, tools, lineage, topology == "circular", log)
                mitfi = annotate_trnas(searched, gfa, refdir, code, args.threads, tools, tmp, log, n, args.phylum) \
                    if tools.get("mitfi") else None
                kept = trna_second_opinion(trnas, mitfi, n)
                if len(kept) < len(trnas):
                    log.info(f"[core] MitoM-ML calls not sure enough without MiTFi: "
                             f"{' '.join(f.name for f in trnas if f not in kept)}")
                trnas = kept
                if mitfi is not None:
                    if vertebrate(lineage):
                        trnas = mitfi_ends(trnas, mitfi, n)
                        log.info(f"[core] vertebrate: tRNA ends of MiTFi for {sum('ends by MiTFi' in f.note for f in trnas)}")
                    missing = [x for x in TRNA_NAMES if x not in {f.name for f in trnas}]
                    if missing:
                        found = rescue_trnas(trnas, mitfi, n)
                        log.info(f"[core] MitoM-ML did not call {' '.join(missing)}; added from MiTFi: "
                                 f"{' '.join(f.name for f in found[len(trnas):]) or 'none'}")
                        trnas = found
                feats += trnas
            elif not args.no_trna:
                feats += annotate_trnas(searched, gfa, refdir, code, args.threads, tools, tmp, log, n, args.phylum)
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
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stderr)] + ([logging.FileHandler(args.log)] if args.log else []))
    return annotate(args)


if __name__ == "__main__":
    sys.exit(main())
