#!/usr/bin/env python3
"""mito_intron.py - intron-aware re-annotation of mitochondrial protein-coding genes.

The annotation core does not annotate introns. This script adds them to its GenBank output with
a Python port of the intron workflow of MFannot (Beck & Lang,
https://github.com/BFL-lab/Mfannot, GPL-2.0; data https://github.com/BFL-lab/MFannot_data,
GPL-3.0; Lang et al. 2023, Front. Plant Sci., doi:10.3389/fpls.2023.1222186).
This file is therefore distributed under GPL-3.0.

Ported steps (MFannot line numbers refer to BFL-lab/Mfannot/mfannot):
  0. Library: MFannot's curated protein collections (models/mfannot_proteins, from
     MFannot_data/protein_collections/identified); gene name from the header
     "<species> mt <gene>;" exactly as mfannot:1180 parses it.
  1. Pre-screen (RunFlip / BlastFlipVSGene / FillFLIPBLASTPROTArrayFromBlastResult /
     SelectProteinForExo, mfannot:962-1480): ORFs >= 40 aa starting at the first
     start codon (flip -m -l 40 -g <code>), blastp vs the library (BLOSUM62); hits
     kept at E <= 1e-10 (and E <= 1e-11 when the top HSP is < 50 % identical); an
     ORF is assigned to every gene whose hit covers a non-overlapping part of it;
     a gene with >= 2 ORFs on one strand separated by >= 142 nt (or ORFs on both
     strands) is a candidate intron-containing gene.
  2. Gene model (FindExonsInHypProtArray, SortAndGroupC4Report, DefineBestSol,
     mfannot:2455-2960): Exonerate protein2genome with the best-hit protein and
     MFannot's parameters (--minintron 142 --maxintron 5000 --intronpenalty -35,
     neutral splice PSSM, genetic code); all reports with raw score >= 200 are
     chained (query overlap < 5 % of protein length, mutually closest in the
     genome) and the group covering most of the protein is used - this joins
     exons around introns longer than --maxintron (up to MAX_CHAIN_INTRON).
  3. Intron type: Infernal cmsearch with Rfam RF00028 (group I), RF00029 + the seven
     group II D1-D4 models (RF01998, RF01999, RF02001, RF02003-05, RF02012), and the 14
     group I subtype models of Selles Vidal et al. 2025 (doi:10.1093/bioadv/vbaf020;
     downloaded once into ~/.cache), replacing MFannot's RNAweasel/ERPIN step.
  4. Boundary adjustment (Adjust_all_intronic_junctions, mfannot:5637-6100): both
     intron ends shifted by the same k in -6..6; group I scores exon-terminal U and
     intron-terminal omega-G (adjust when >= 19), group II GUGYG...AY (>= 40).

Documented deviations: step 3 uses Rfam instead of RNAweasel; in step 4 each
shift is judged against the homolog protein (MFannot compares with the starting
model, which Exonerate can get wrong at the junction - Metridium senile COX1:
model HMFTSG, NCBI and homolog HMFTVG - and would then reject the true boundary);
3' ends use the incomplete-stop (T / TA) convention of animal mtDNA (an existing CDS end is
reused only when it is a stop / T-TA tail); a trans-spliced gene is not joined: its cis-spliced
pieces are written as separate partial CDS (below); genes cut by both ends of a linear record are
not handled; chained alignments more than MAX_CHAIN_INTRON nt apart are not joined into an intron.
Since 2026-10-29 (benchmark: 37 NCBI records with introns + 72 without, core output rotated as the
MitoM pipeline gives it to this step):
  - step 1 also takes every standard gene whose CDS from MitoM's annotation core leaves
    >= STUB_MIN_AA residues of its reference protein unaligned at one end (the core is intron-less
    and annotates one exon; its CDS note names the protein, MITOS2 refseq89m), with that protein as the
    Exonerate query; such a gene is written only when Rfam types all its introns (Decemunciger nad4,
    exon 2 of 52 aa, is not found by the ORF pre-screen);
  - step 2 chains alignments only in protein order along the transcript (MFannot takes the closest
    report in either direction, which joined the permuted exons of corallimorpharian nad5 into a
    false 300-nt intron); on a circular record the chain may cross the origin and the gene is then
    modelled on a rotated copy of the record and written back across the origin (MitoM rotates
    records to start at the core's cox1 piece; when that piece is exon 2, exon 1 and the intron lie
    at the record's end - 6 of 12 annelid records with cox1 introns were not modelled);
  - step 4: an untyped intron (no Rfam hit) is placed by the group II rule when GUGYG...AY reaches
    TYPE_II_MIN within the shift range (degenerate Chaetopterus group II introns, GUGCG...AC); a model
    with an untyped intron is written only when it aligns >= UNTYPED_MIN_COVERAGE of its homolog;
  - mini-exons are searched only for homolog gaps <= MINI_EXON_MAX_AA in introns <= MAX_INTRON;
  - a model whose intron holds a copy of its own exons (copied_exons) is refused: exonerate splices the
    copies of a duplicated gene (Hexamermis NC_008828.1 atp6 x3: a 3.7-kb false intron);
  - a chain with exons on both strands, or two pieces more than MAX_INTRON apart that leave >= CHAIN_HOLE_AA
    homolog residues unaligned (an exon elsewhere: trans-spliced placozoan cox1), is split into runs of
    collinear same-strand alignments; each run with an intron is written as a partial CDS with its introns
    (gene_modules; MFannot reports such genes without a model);
  - 3' end: a first complete stop more than STOP_OVERLAP_MAX nt inside the next feature gives way to a T / TA
    right before it (Glycera KT989324.1 cox1 ran 717 nt into cox2 and was refused as overlapping COX2).
Since 2026-10-29 evening (gold standard of the 193 local NCBI intron records, intron_gold_2026-10-29):
  - step 4 judges a shift by the median homolog-score gain over up to five homologs (the model's, the core's
    reference of the gene and three of the library) instead of the model's homolog alone;
  - step 4 first moves the two intron ends independently (5' by up to REFINE_5, 3' by up to REFINE_3 nt, frame kept)
    when the spliced protein then aligns at least REFINE_MIN_GAIN better (refine_junction; exonerate put the 3' end of
    hexacoral nad5 introns a few codons off), then applies the same-shift splice rules as before;
  - a typed group II intron without a GUGYG / AY placement within 6 nt takes one that costs at most II_AA_TOLERANCE of
    homolog score, else one up to WIDE_II nt away; a group I subtype hit weaker than WEAK_SUBTYPE_E gives no type
    when the intron has group II splice sites scoring TYPE_II_STRONG (Endomyzostoma cox1);
  - chains across the origin may share WRAP_OVERLAP % of the protein (5 % elsewhere);
  - the introns around an inserted mini-exon get the boundary step again (5' window MINI_REFINE): exonerate had
    spread the mini-exon's residues over the neighbouring junction (Plakinastrella NC_010217.1 cox1, +3 aa).

Usage:
  python3 mito_intron.py -i final_gb_dir -o out_dir --genetic_code 5 [--bin env/bin]
Outputs: out_dir/<same relative paths> (changed records rewritten, others copied byte
for byte) and out_dir/intron_report.tsv.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from Bio import SeqIO
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402  (finds external programs: config, PATH, conda environments)
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import (AfterPosition, BeforePosition, CompoundLocation, ExactPosition, FeatureLocation,
                            SeqFeature)

HERE = Path(__file__).resolve().parent
DEFAULT_LIB = HERE / "models" / "mfannot_proteins"
DEFAULT_CM = HERE / "models" / "rfam_group_I_II_introns.cm"
SPLICE_PSSM = HERE / "models" / "neutral.pssm"      # MFannot_data/models/Splice/neutral.pssm
from mito_common import (GB_SUFFIXES as GB_EXT, move_position_qualifiers, stop_qualifiers, write_gb,  # noqa: E402
                         write_tsv)                                                                       # (CC0)
ADJUST_RANGE = range(-6, 7)          # mfannot:5640 %adjustement, same k at both ends
TYPE_I_MIN, TYPE_II_MIN = 19, 40     # mfannot:5735 / Adjust_TypeII
II_AA_TOLERANCE = 5                  # homolog score (median, BLOSUM62) a group II placement may lose when none is free
REFINE_5, REFINE_3, REFINE_MIN_GAIN = 6, 24, 10   # refine_junction: windows (nt) and median homolog gain needed
WIDE_II = 15                         # same-shift window for a group II intron with no canonical placement within 6 nt
MINI_REFINE = 12                     # refine_junction's 5' window for the introns around an inserted mini-exon
WEAK_SUBTYPE_E, TYPE_II_STRONG = 1e-3, 55   # a group I subtype hit weaker than this yields to group II sites this strong
# MFannot defaults (PirModels/Option.pir, mfannot:284)
MIN_ORF_AA, BLAST_EVALUE, SHORT_HSP_EVALUE = 40, 1e-10, 1e-11
MIN_INTRON, MAX_INTRON, MIN_REPORT_SCORE = 142, 5000, 200
MINI_EXON_MAX_AA = 10                # longest homolog gap at an exon junction searched for a mini-exon (codons)
STOP_OVERLAP_MAX = 30                # nt a 3' end may run into the next feature to reach a complete stop (mito_curate P3)
COPY_KMER, COPY_MIN_SHARED = 20, 60   # copied_exons: k-mer length, shared k-mers that make an intron hold an exon copy
                                     # (Hexamermis atp6 copies 340; one 33-46-nt repeat in hexacoral nad5 introns 14-27)
MAX_CHAIN_INTRON = 25000             # longest intron made by chaining separate alignments (exonerate itself stops at
                                     # --maxintron): metazoan mt group I / II introns are < 3 kb except the hexacoral
                                     # nad5 group I intron that hosts other genes, 11,485 nt in Montipora NC_006902.1 and
                                     # 18,013-19,887 nt in corallimorpharians (NC_027103.1, NC_027106.1 around the
                                     # origin); the two pieces of a gene cut by the ends of a linear / rotated record
                                     # (stress case 2026-10: 14,901 nt) are out of protein order and refused by
                                     # build_model, not by this length
CHAIN_OVERLAP = 30                   # nt two chained alignments may overlap on the genome and still count as collinear
WRAP_OVERLAP = 15                    # % of the protein two alignments chained across the origin may share (5 % else)
CHAIN_HOLE_AA = 20                   # homolog residues two alignments chained across more than MAX_INTRON may leave
                                     # unaligned: placozoan cox1 is trans-spliced, its 31-35-aa exon on the other strand
                                     # lies inside the 'intron' of the + strand pieces (NC_041549.1: 17,704 nt false intron);
                                     # the chained exons of the 11-20 kb hexacoral nad5 introns meet within 0 residues
ACROSS_ORIGIN = "gene spans the origin of the circular record"
TRANS_SPLICED = "trans-spliced"
HEADER_RE = re.compile(r"(pt|pt_cyano|mt|mt_alpha)\s*([\S]*)\s*;", re.I)   # mfannot:1180
# the annotation core's CDS note: reference protein (MITOS2 refseq89m featureProt) and the residues it aligned
CORE_REF_RE = re.compile(r"ref (\S+) score (\d+) residues (\d+)-(\d+)/(\d+)")
STUB_MIN_AA = 30                     # reference residues the core's CDS leaves unaligned at one end that make the
                                     # gene an intron candidate (an exon the intron-less core could not reach)


@dataclass
class LibProt:
    species: str
    gene: str
    seq: str


@dataclass
class Orf:
    strand: int
    start: int                        # transcript coords (1-based) on its strand
    end: int
    evalue: float = 1.0
    prot_id: str = ""
    q_min: int = 0                    # ORF aa range covered by the hit's HSPs
    q_max: int = 0


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


TIMEOUT = 900                          # seconds per external program call (--timeout)


def run_tool(cmd: List[str], what: str) -> str:
    """Run an external program; on timeout or failure return '' and warn, so one hard
    sequence cannot stall or kill a batch."""
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=TIMEOUT).stdout
    except subprocess.TimeoutExpired:
        print(f"[mito_intron] {what} timed out after {TIMEOUT}s - skipped", file=sys.stderr)
    except OSError as exc:
        print(f"[mito_intron] {what} failed: {exc}", file=sys.stderr)
    return ""


# ── step 0: library ──────────────────────────────────────────────────────────
def iter_files(paths: List[str], exts: set) -> List[Path]:
    out: List[Path] = []
    for p in map(Path, paths):
        if p.is_dir():
            out += sorted(x for x in p.rglob("*") if x.suffix.lower() in exts and x.is_file())
        elif p.is_file():
            out.append(p)
    return out


def norm_species(s: str) -> str:
    return " ".join(s.replace("_", " ").split()).lower()


def load_library(paths: List[str]) -> Dict[str, LibProt]:
    """MFannot-format protein FASTA ('>Homo sapiens mt cox1; 513 aa')."""
    lib: Dict[str, LibProt] = {}
    for f in iter_files(paths, {".faa", ".fa", ".fasta", ".pep"}):
        for r in SeqIO.parse(str(f), "fasta"):
            m = HEADER_RE.search(r.description)
            if not m or not m.group(2):
                continue
            lib[f"L{len(lib):05d}"] = LibProt(norm_species(r.description[:m.start()]), m.group(2),
                                              str(r.seq).upper().rstrip("*"))
    return lib


def load_ref_gb(paths: List[str], gencode: int) -> Dict[str, LibProt]:
    """Extra homologs from close-relative GenBank records (MitoFinder's reference idea): CDS
    translations of the standard mitochondrial protein genes, gene names normalised to the
    MFannot library spelling so they group with it; CDS with internal stops are skipped."""
    canon = {p.stem.lower(): p.stem for p in DEFAULT_LIB.glob("*.faa")}
    by_product = {v: k for k, v in PRODUCTS.items()}
    lib: Dict[str, LibProt] = {}
    for f in iter_files(paths, GB_EXT):
        for rec in SeqIO.parse(str(f), "genbank"):
            sp = norm_species(rec.annotations.get("organism", ""))
            for ft in rec.features:
                if ft.type != "CDS":
                    continue
                g = norm_gene(ft.qualifiers.get("gene", [""])[0])
                if g not in PRODUCTS:
                    g = by_product.get(ft.qualifiers.get("product", [""])[0].strip().lower(), "")
                if g not in PRODUCTS:
                    continue
                prot = ft.qualifiers.get("translation", [""])[0]
                if not prot:
                    table = int(ft.qualifiers.get("transl_table", [gencode])[0])
                    nt = ft.location.extract(rec.seq)
                    prot = str(nt[: len(nt) - len(nt) % 3].translate(table=table))
                prot = prot.rstrip("*").upper()
                if "*" in prot or len(prot) < 50:
                    continue
                lib[f"R{len(lib):05d}"] = LibProt(sp, canon.get(g.lower(), g), prot)
    return lib


_REFSEQ: Dict[str, Dict[str, str]] = {}


def refseq_protein(refdir: Optional[Path], ref_id: str) -> str:
    """Protein ref_id ('NC_027615:cox1-a-1-3175-5620') of the MITOS2 reference data (refseq89m/featureProt) that
    the annotation core aligns; '' when the data are not found."""
    gene = ref_id.split(":", 1)[-1].split("-", 1)[0]
    if gene not in _REFSEQ:
        f = Path(refdir or "") / "featureProt" / f"{gene}.fas"
        _REFSEQ[gene] = ({r.id: str(r.seq).upper().rstrip("*") for r in SeqIO.parse(str(f), "fasta")}
                         if refdir and f.is_file() else {})
    return _REFSEQ[gene].get(ref_id, "")


def core_references(rec, refdir: Optional[Path]) -> Dict[str, dict]:
    """Reference protein of each standard protein gene the annotation core annotated, from its CDS note (CORE_REF_RE):
    the closest homolog of the core's blastx search over all metazoan RefSeq proteins, far closer to a bilaterian
    gene than MFannot's library (annelids: Homo sapiens). {library gene name: {id, seq, strand, score, stub}};
    stub = the core's CDS leaves >= STUB_MIN_AA residues of the reference unaligned at one end."""
    canon = {p.stem.lower(): p.stem for p in DEFAULT_LIB.glob("*.faa")}
    out: Dict[str, dict] = {}
    for f in rec.features:
        m = CORE_REF_RE.search(" ".join(f.qualifiers.get("note", []))) if f.type == "CDS" else None
        g = norm_gene(f.qualifiers.get("gene", [""])[0])
        if not m or g not in PRODUCTS:
            continue
        g = canon.get(g.lower(), g)
        score, a, b, n = (int(x) for x in m.group(2, 3, 4, 5))
        seq = refseq_protein(refdir, m.group(1))
        if seq and (g not in out or score > out[g]["score"]):
            out[g] = {"id": m.group(1), "seq": seq, "strand": f.location.strand, "score": score,
                      "stub": a - 1 >= STUB_MIN_AA or n - b >= STUB_MIN_AA}
    return out


# ── step 1: pre-screen (flip + blastp + SelectProteinForExo) ─────────────────
def find_orfs(seq: str, gencode: int, min_aa: int = MIN_ORF_AA) -> List[Tuple[Orf, str]]:
    """flip -m -l <min_aa>: stop-to-stop segments, ORF from the first start codon."""
    table = CodonTable.unambiguous_dna_by_id[gencode]
    starts, stops = set(table.start_codons), set(table.stop_codons)
    out: List[Tuple[Orf, str]] = []
    for strand, s in ((1, seq), (-1, str(Seq(seq).reverse_complement()))):
        L = len(s)
        for frame in range(3):
            seg = frame
            for i in range(frame, L - 2, 3):
                last = i + 6 > L
                if s[i:i + 3] in stops or last:
                    end = i if s[i:i + 3] in stops else i + 3
                    j = next((x for x in range(seg, end, 3) if s[x:x + 3] in starts), None)
                    if j is not None and (end - j) // 3 >= min_aa:
                        prot = "M" + str(Seq(s[j + 3:end]).translate(table=gencode))
                        out.append((Orf(strand, j + 1, end), prot))
                    seg = i + 3
    return out


def prescreen(seq: str, lib: Dict[str, LibProt], lib_db: Path, gencode: int, tools: Dict[str, Optional[str]],
              tmp: Path, exclude_species: str = "", threads: int = 1) -> Dict[str, dict]:
    """Candidate intron-containing genes: {gene: {'prot': lib id, 'strand': 1|-1|2}}."""
    orfs = find_orfs(seq, gencode)
    if not orfs:
        return {}
    q = tmp / "orfs.faa"
    q.write_text("".join(f">O{i}\n{p}\n" for i, (_, p) in enumerate(orfs)))
    res = run_tool([str(tools["blastp"]), "-query", str(q), "-db", str(lib_db), "-matrix", "BLOSUM62",
                    "-outfmt", "6 qseqid sseqid evalue pident qstart qend", "-num_threads", str(threads)], "blastp")
    hits: Dict[str, Dict[str, List[Tuple[float, float, int, int]]]] = defaultdict(dict)
    for line in res.splitlines():
        qid, sid, ev_s, pid_s, qs, qe = line.split("\t")
        if exclude_species and lib[sid].species == exclude_species:
            continue
        hits[qid].setdefault(sid, []).append((float(ev_s), float(pid_s), int(qs), int(qe)))
    genes: Dict[str, List[Orf]] = defaultdict(list)
    for qid, by_subject in hits.items():
        base = orfs[int(qid[1:])][0]
        ranked = sorted(by_subject.items(), key=lambda kv: min(h[0] for h in kv[1]))
        accepted: List[Tuple[int, int]] = []
        for sid, hsps in ranked:
            ev, pid0 = min(h[0] for h in hsps), hsps[0][1]
            if ev > BLAST_EVALUE or (pid0 < 50 and ev > SHORT_HSP_EVALUE):
                break                                         # mfannot:1160 / 1204-1205
            good = [h for h in hsps if h[0] <= BLAST_EVALUE] or hsps[:1]
            fr = (min(h[2] for h in good), max(h[3] for h in good))
            if any(not (fr[1] < a or fr[0] > b) for a, b in accepted):
                continue                                      # overlapping frame (mfannot:1226-1233)
            accepted.append(fr)
            genes[lib[sid].gene].append(Orf(base.strand, base.start, base.end, ev, sid, fr[0], fr[1]))
    cands: Dict[str, dict] = {}
    for gene, os_ in genes.items():
        if len(os_) < 2:
            continue
        prot = min(os_, key=lambda o: o.evalue).prot_id
        strands = {o.strand for o in os_}
        if len(strands) == 2:
            cands[gene] = {"prot": prot, "strand": 2}
            continue
        segs = sorted((o.start + 3 * o.q_min, o.start + 3 * o.q_max) for o in os_)
        if any(abs(segs[i][1] - segs[i + 1][0]) >= MIN_INTRON for i in range(len(segs) - 1)):
            cands[gene] = {"prot": prot, "strand": os_[0].strand}       # mfannot:1447
    return cands


# ── step 2: Exonerate + report grouping ──────────────────────────────────────
@dataclass
class Report:
    qab: int
    qae: int
    tab: int
    tae: int
    score: int
    strand: int
    exons: List[Tuple[int, int]]


def run_exonerate(exe: str, prot_faa: Path, genome_fna: Path, gencode: int, maxintron: int) -> List[Report]:
    cmd = [exe, "--model", "protein2genome", "--showalignment", "no", "--showvulgar", "no",
           "--showtargetgff", "yes", "--minintron", str(MIN_INTRON), "--maxintron", str(maxintron),
           "--intronpenalty", "-35", "--geneticcode", str(gencode),
           *(["--splice3", str(SPLICE_PSSM), "--splice5", str(SPLICE_PSSM)] if SPLICE_PSSM.exists() else []),
           # %tS = target strand (%g, the splice-derived gene orientation, is "." without introns)
           "--ryo", ">>RYO\t%qab\t%qae\t%tab\t%tae\t%s\t%tS\n", str(prot_faa), str(genome_fna)]
    out = run_tool(cmd, "exonerate")
    reps: List[Report] = []
    exons: List[Tuple[int, int]] = []
    for line in out.splitlines():
        if line.startswith(">>RYO"):
            _, qab, qae, tab, tae, score, strand = line.split("\t")
            if exons:
                reps.append(Report(int(qab), int(qae), int(tab), int(tae), int(score),
                                   1 if strand == "+" else -1, sorted(exons)))
            exons = []
        elif "\texon\t" in line:
            c = line.split("\t")
            exons.append((int(c[3]), int(c[4])))
    return reps


def chain_step(a: List[Report], b: List[Report], L: int = 0, circular: bool = False) -> Optional[int]:
    """Distance (nt, transcript direction) from the end of chain a to the start of chain b, a before b in the
    protein; None when b lies upstream of a on a linear sequence (not collinear). On a circular record the
    distance runs across the origin (a gene whose exons lie on both sides of the origin)."""
    if a[-1].strand != b[0].strand:
        return abs(a[-1].tae - b[0].tab)              # exons on both strands: reported as trans-spliced later
    dt = (b[0].tab - a[-1].tae) * a[-1].strand
    if dt >= -CHAIN_OVERLAP:
        return max(dt, 0)
    return dt + L if circular and L else None


def group_reports(reps: List[Report], qlen: int, L: int = 0, circular: bool = False) -> List[List[Report]]:
    """SortAndGroupC4Report + DefineBestSol (mfannot): chain reports into gene models. Deviation: two
    alignments are chained only in protein order along the transcript (collinear; MFannot takes the genome
    distance in either direction), across the origin on circular records."""
    groups = [[r] for r in reps if r.score >= MIN_REPORT_SCORE]

    def ov(a: List[Report], b: List[Report]) -> float:
        """Query overlap allowed between chain a and the chain b after it: MFannot's 5 % of the protein, WRAP_OVERLAP
        across the origin of a circular record (both pieces then run into the intron, as exonerate cannot splice
        across the origin; the gene is aligned again on a rotated copy)."""
        wraps = circular and a[-1].strand == b[0].strand and (b[0].tab - a[-1].tae) * a[-1].strand < -CHAIN_OVERLAP
        return (WRAP_OVERLAP if wraps else 5) * qlen / 100

    def best_sol(i: int) -> Optional[int]:
        ps, pe = groups[i][0].qab, groups[i][-1].qae
        best, dmin = None, None
        for j in range(len(groups) - 1, -1, -1):
            if j == i:
                continue
            ps2, pe2 = groups[j][0].qab, groups[j][-1].qae
            if ps < ps2:
                if not (ps2 > pe or pe - ps2 < ov(groups[i], groups[j])):
                    continue
                d = chain_step(groups[i], groups[j], L, circular)
            elif ps > ps2:
                if not (ps > pe2 or pe2 - ps < ov(groups[j], groups[i])):
                    continue
                d = chain_step(groups[j], groups[i], L, circular)
            else:
                continue
            if d is not None and (dmin is None or d < dmin):
                best, dmin = j, d
        return best

    i = 0
    while i < len(groups):
        j = best_sol(i)
        if j is None or best_sol(j) != i:
            i += 1
            continue
        a, b = groups[i], groups[j]
        groups[i] = a + b if a[0].qab < b[0].qab else b + a
        del groups[j]
        if j < i:
            i -= 1
    groups.sort(key=lambda g: (-sum(abs(r.qae - r.qab) + 1 for r in g), len(g)))
    return groups


def wrap_steps(reps: List[Report]) -> List[int]:
    """Indices i where reps[i] -> reps[i + 1] (protein order) steps backwards along the transcript."""
    rs = sorted(reps, key=lambda r: r.qab)
    return [i for i in range(len(rs) - 1) if (rs[i + 1].tab - rs[i].tae) * rs[i].strand < -CHAIN_OVERLAP]


def build_model(prot_id: str, gene: str, reps: List[Report], qlen: int, circular: bool = False
                ) -> Tuple[Optional[Model], str]:
    """Chained reports -> gene model, or (None, reason). Reports out of protein order along the transcript are
    a gene across the origin of a circular record (ACROSS_ORIGIN: modelled on a rotated copy) or, on a linear
    record, the pieces of a gene cut by its ends (or paralogous hits): not modelled."""
    strands = {r.strand for r in reps}
    if len(strands) != 1:
        return None, f"exons on both strands ({TRANS_SPLICED}?) - not modelled"
    back = wrap_steps(reps)
    if back:
        if circular and len(back) == 1:
            return None, ACROSS_ORIGIN
        return None, ("alignments out of protein order along the sequence (gene cut by the ends of a linear "
                      "record, or paralogous hits) - not modelled")
    exons: List[Tuple[int, int]] = []
    for e in sorted(e for r in reps for e in r.exons):
        if exons and e[0] - exons[-1][1] - 1 < MIN_INTRON:   # gap between chained reports, not an intron
            exons[-1] = (exons[-1][0], max(exons[-1][1], e[1]))
        else:
            exons.append(e)
    if len(exons) < 2:
        return None, ""
    longest = max(exons[i + 1][0] - exons[i][1] - 1 for i in range(len(exons) - 1))
    if longest > MAX_CHAIN_INTRON:
        return None, (f"chained alignments {longest} nt apart, longer than any known mt intron ({MAX_CHAIN_INTRON} nt) "
                      "- not modelled")
    rs = sorted(reps, key=lambda r: r.qab)
    for x, y in zip(rs, rs[1:]):
        if (y.tab - x.tae) * x.strand > MAX_INTRON and y.qab - x.qae >= CHAIN_HOLE_AA:
            return None, (f"chained alignments {(y.tab - x.tae) * x.strand} nt apart leave residues {x.qae + 1}-{y.qab} "
                          f"of the homolog unaligned: an exon elsewhere ({TRANS_SPLICED} gene?) - not modelled")
    coding = sum(b - a + 1 for a, b in exons)
    if coding > 3 * qlen + 33:            # MitoFinder geneChecker_fasta_gaps.py:150 length sanity check
        return None, (f"chained exons ({coding} nt) longer than the homolog allows ({3 * qlen + 33} nt) "
                      "- likely joins pieces of another gene; not modelled")
    covered = sum(abs(r.qae - r.qab) for r in reps)
    return Model(prot_id, gene, strands.pop(), exons, sum(r.score for r in reps), covered / max(qlen, 1),
                 qlen - max(r.qae for r in reps), min(r.qab for r in reps)), ""


def gene_modules(prot_id: str, gene: str, reps: List[Report], qlen: int, circular: bool = False) -> List[Model]:
    """Cis-spliced pieces of a gene whose chained alignments lie on both strands or leave an exon elsewhere (trans-
    spliced placozoan cox1): runs of collinear same-strand alignments in protein order, each modelled with its own
    introns. Pieces with an intron are returned, their ends open where another piece continues the protein."""
    rs = sorted(reps, key=lambda r: r.qab)
    runs = [[rs[0]]]
    for x, y in zip(rs, rs[1:]):
        step = (y.tab - x.tae) * x.strand
        if y.strand != x.strand or step < -CHAIN_OVERLAP or (step > MAX_INTRON and y.qab - x.qae >= CHAIN_HOLE_AA):
            runs.append([y])
        else:
            runs[-1].append(y)
    out = []
    for i, run in enumerate(runs):
        m, _ = build_model(prot_id, gene, run, qlen, circular)
        if m is not None and len(m.exons) > 1:
            m.open5, m.open3 = i > 0, i < len(runs) - 1
            m.piece = (f"{TRANS_SPLICED} gene: piece {i + 1} of {len(runs)} (homolog residues {run[0].qab + 1}-"
                       f"{run[-1].qae}), the other pieces not joined")
            out.append(m)
    return out


# ── step 2: intron type with Infernal + Rfam ────────────────────────────────
SUBTYPE_CM_URL = "https://api.github.com/repos/LaraSellesVidal/Group1IntronDatabase/contents/subtype_CMs"
SUBTYPE_CM_DIR = Path.home() / ".cache" / "mito_intron" / "selles_vidal_2025_subtype_CMs"
SUBTYPE_REF = "Selles Vidal et al. 2025, doi:10.1093/bioadv/vbaf020"


def subtype_cms() -> List[Path]:
    """Group I subtype covariance models (IA1..IE3) of Selles Vidal et al. 2025, trained on 42,459
    group I introns. Not redistributed (no licence stated): downloaded once into a local cache;
    if the download fails the bundled Rfam models are used alone."""
    if not SUBTYPE_CM_DIR.exists() or not list(SUBTYPE_CM_DIR.glob("*.cm")):
        try:
            import json
            import urllib.request
            SUBTYPE_CM_DIR.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(SUBTYPE_CM_URL, timeout=30) as r:
                items = json.load(r)
            for it in items:
                if it["name"].endswith(".cm"):
                    urllib.request.urlretrieve(it["download_url"], SUBTYPE_CM_DIR / it["name"])
        except Exception as exc:                       # offline: carry on with Rfam only
            print(f"[mito_intron] group I subtype models not available ({exc}); using Rfam only",
                  file=sys.stderr)
    return sorted(SUBTYPE_CM_DIR.glob("*.cm"))


def intron_types(cmsearch: str, cms: List[Path], seqs: Dict[str, str], evalue: float, tmp: Path
                 ) -> Dict[str, Tuple[str, float, str]]:
    """Best covariance-model hit per intron: (group 'I'/'II', E-value, model name).
    Models: Rfam Intron_gpI (RF00028), Intron_gpII (RF00029, domains V-VI), the seven group II
    D1-D4 models (RF01998, RF01999, RF02001, RF02003-RF02005, RF02012), and the group I subtype
    models of Selles Vidal et al. 2025 when available."""
    cms = [c for c in cms if c.exists()]
    if not seqs or not cms:
        return {}
    cm = tmp / "all_models.cm"
    cm.write_text("".join(c.read_text() for c in cms))
    best: Dict[str, Tuple[str, float, str]] = {}

    def search(todo: Dict[str, str], fast: bool) -> None:
        fna, tbl = tmp / "introns.fna", tmp / "introns.tbl"
        fna.write_text("".join(f">{k}\n{v}\n" for k, v in todo.items()))
        tbl.unlink(missing_ok=True)
        # --mxsize 128 MB: a long group I subtype model hit took 266 MB (NC_027103.1 nad5 intron); same models and
        # E-values within 10 % at 128 (whole run: 293 -> 183 MB fast pass, 977 -> 237 MB full pass); 64 loses hits
        # --cpu 0 (no worker threads, each holds its own DP matrices): same hits and time, 187 -> 165 MB on NC_027103.1
        run_tool([cmsearch, "--noali", "--toponly", *(["--rfam"] if fast else []), "--cpu", "0", "--mxsize", "128",
                  "--tblout", str(tbl), str(cm), str(fna)], "cmsearch")
        for line in tbl.read_text().splitlines() if tbl.exists() else []:
            if line.startswith("#"):
                continue
            c = line.split()
            name, model, strand, e = c[0], c[2], c[9], float(c[15])
            if strand != "+" or e > evalue:
                continue
            kind = "II" if (model.endswith("gpII") or model.startswith("group-II")) else "I"
            if name not in best or e < best[name][1]:
                best[name] = (kind, e, model)

    # two passes: Infernal's fast --rfam filters first; full-sensitivity search only for introns still
    # untyped (benchmark: same 147/165 NCBI introns typed, same group for all, ~2.4x faster)
    search(seqs, fast=True)
    rest = {k: v for k, v in seqs.items() if k not in best}
    if rest:
        search(rest, fast=False)
    return best


def type_label(kind: str, ev: float, model: str) -> str:
    if not kind:
        return "intron, type not identified"
    if re.fullmatch(r"I[A-E]\d?", model):
        return f"group I intron, subtype {model} ({SUBTYPE_REF} model, E={ev:.2g})"
    acc = {"Intron_gpI": "RF00028", "Intron_gpII": "RF00029"}.get(model, model)
    return f"group {kind} intron (Infernal/Rfam {acc}, E={ev:.2g})"


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
        fs, fe = int(ft.location.start) + 1, int(ft.location.end)
        if overlap((fs, fe), g_span) >= 30:
            continue
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
        fs, fe = int(ft.location.start) + 1, int(ft.location.end)
        if overlap((fs, fe), g_span) >= 30:          # the gene being re-modelled (or a nested ORF)
            continue
        start_tx = fs if m.strand == 1 else L + 1 - fe
        if start_tx > end_tx:
            limit = min(limit, start_tx - 1)
    return limit


# ── write-back ───────────────────────────────────────────────────────────────
def spans_record_ends(m: Model, L: int, circular: bool, maxintron: int) -> str:
    """Reason not to model a gene whose exons touch both ends of a linear record: the two pieces of a gene cut by a
    rotation, chained across a gap longer than exonerate's own --maxintron (stress case 2026-10: 14,901 nt on a
    16,560-nt record). A short intron in a fragment that runs off both ends (LC670763 cox1, join(<1..235,1101..>1523))
    is kept. On a circular record the model is collinear (build_model), so its exons at both ends are one gene with a
    long intron (corallimorpharian nad5); a gene across the origin is modelled on a rotated copy."""
    if circular or not (m.exons[0][0] <= 50 and m.exons[-1][1] >= L - 50):
        return ""
    gap = max(m.exons[i + 1][0] - m.exons[i][1] - 1 for i in range(len(m.exons) - 1))
    return (f"gene cut by both ends of the linear record, pieces {gap} nt apart (rotated genome?) - not modelled"
            if gap > maxintron else "")


def to_transcript(m: Model, L: int) -> List[Tuple[int, int]]:
    if m.strand == 1:
        return list(m.exons)
    return sorted((L + 1 - b, L + 1 - a) for a, b in m.exons)


def to_genome(exons: List[Tuple[int, int]], strand: int, L: int) -> List[Tuple[int, int]]:
    return list(exons) if strand == 1 else sorted((L + 1 - b, L + 1 - a) for a, b in exons)


def overlap(a: Tuple[int, int], b: Tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


# Standard products of the mitochondrial protein genes (NCBI mitogenome submissions)
PRODUCTS = {
    "cox1": "cytochrome c oxidase subunit I", "cox2": "cytochrome c oxidase subunit II",
    "cox3": "cytochrome c oxidase subunit III", "cob": "cytochrome b", "cytb": "cytochrome b",
    "atp6": "ATP synthase F0 subunit 6", "atp8": "ATP synthase F0 subunit 8", "atp9": "ATP synthase F0 subunit 9",
    "nad1": "NADH dehydrogenase subunit 1", "nad2": "NADH dehydrogenase subunit 2",
    "nad3": "NADH dehydrogenase subunit 3", "nad4": "NADH dehydrogenase subunit 4",
    "nad4l": "NADH dehydrogenase subunit 4L", "nad5": "NADH dehydrogenase subunit 5",
    "nad6": "NADH dehydrogenase subunit 6",
}
ALIASES = {"coi": "cox1", "coii": "cox2", "coiii": "cox3", "cytb": "cob", "nd1": "nad1", "nd2": "nad2",
           "nd3": "nad3", "nd4": "nad4", "nd4l": "nad4l", "nd5": "nad5", "nd6": "nad6"}


def norm_gene(gene: str) -> str:
    g = gene.lower().replace("-", "").replace("_", "")
    return ALIASES.get(g, g)


def product_for(gene: str) -> str:
    return PRODUCTS.get(norm_gene(gene), f"{gene} protein")


def conflicting_genes(rec, gene: str, strand: int, exons: List[Tuple[int, int]]) -> List[str]:
    """Existing CDS (e.g. the core's) of a *different* standard mitochondrial protein gene that the
    model's exons overlap by >= 30 nt (intron-encoded ORFs such as heg are not conflicts)."""
    out = []
    for f in rec.features:
        if f.type != "CDS" or f.location.strand != strand or not f.qualifiers.get("gene"):
            continue
        other = f.qualifiers["gene"][0]
        if norm_gene(other) == norm_gene(gene) or norm_gene(other) not in PRODUCTS:
            continue                      # same gene, or intron-encoded ORF (heg/LAGLIDADG/orf...) that may overlap exons
        parts = [(int(p.start) + 1, int(p.end)) for p in f.location.parts]
        if sum(overlap(x, y) for x in exons for y in parts) >= 30:
            out.append(other)
    return out


def repair_chain(rec, prot_id: str, gene: str, groups: List[List[Report]], qlen: int
                 ) -> Optional[Tuple[Model, List[Tuple[int, int]]]]:
    """Best chain minus its exons on another gene's CDS, chained with the next report that lies on the
    side the dropped exons were on. NC_006902.1 ND5: exonerate aligns nad5 residues 135-600 as an exon on
    ND4 (nad4/nad5 are homologous) + the true exon 2, and residues 1-239 as the true exon 1 in a separate
    report; their query overlap (> 5 %) keeps group_reports from chaining them. Returns (model, dropped exons)."""
    best = groups[0]
    strand = best[0].strand
    kept = [Report(r.qab, r.qae, r.tab, r.tae, r.score, r.strand,
                   [e for e in r.exons if not conflicting_genes(rec, gene, strand, [e])]) for r in best]
    lo = min(e[0] for r in best for e in r.exons)
    k_ex = [e for r in kept for e in r.exons]
    if not k_ex or len(k_ex) == sum(len(r.exons) for r in best):
        return None
    k_lo, k_hi = min(a for a, _ in k_ex), max(b for _, b in k_ex)
    dropped_low = lo < k_lo                            # dropped exon(s) before the kept ones in genome order
    dropped = [e for r in best for e in r.exons if e not in k_ex]
    for g in groups[1:6]:
        g_ex = [e for r in g for e in r.exons]
        if g[0].strand != strand or conflicting_genes(rec, gene, strand, g_ex):
            continue
        if (dropped_low and max(b for _, b in g_ex) < k_lo) or (not dropped_low and min(a for a, _ in g_ex) > k_hi):
            m, _ = build_model(prot_id, gene, g + [r for r in kept if r.exons], qlen)
            if m is not None:
                return m, dropped
    return None


def realign_masked(rec, seq: str, prot_id: str, gene: str, groups: List[List[Report]], qlen: int, pf: Path,
                   tmp: Path, a, tools: Dict[str, Optional[str]]) -> Optional[Model]:
    """repair_chain, then exonerate again with the dropped exons masked (N) and --maxintron long enough for
    the repaired model, so both splice sites of the new junction come from one spliced alignment (the
    repaired junction keeps the acceptor exonerate chose for the wrong donor: 8 nt off in NC_006902.1)."""
    rep = repair_chain(rec, prot_id, gene, groups, qlen)
    if rep is None:
        return None
    m, dropped = rep
    masked = list(seq)
    for x, y in dropped:
        masked[x - 1:y] = "N" * (y - x + 1)
    fna = tmp / "masked.fna"
    fna.write_text(">g\n" + "".join(masked) + "\n")
    longest = max(m.exons[i + 1][0] - m.exons[i][1] for i in range(len(m.exons) - 1))
    reps = run_exonerate(str(tools["exonerate"]), pf, fna, a.genetic_code, max(a.maxintron, longest + 200))
    for g in group_reports(reps, qlen)[:3]:
        alt, _ = build_model(prot_id, gene, g, qlen)
        if alt is not None and alt.strand == m.strand and not conflicting_genes(rec, gene, alt.strand, alt.exons) \
                and overlap(alt.span, m.span) > 0:
            return alt
    return m


def apply_model(rec, m: Model, exons: List[Tuple[int, int]], intron_info: List[dict], gencode: int,
                partial5: bool = False, partial3: bool = False) -> str:
    """Replace the gene/CDS features that overlap the exons; add intron features. Returns gene name."""
    strand = m.strand

    def shared(ft) -> int:
        return sum(overlap((int(p.start) + 1, int(p.end)), x) for p in ft.location.parts for x in exons)
    # features that share >= 30 nt with the exons; a gene feature only for a hit CDS of that name (or without CDS):
    # a gene modelled earlier whose intron hosts this one (hexacoral nad5 around cox1) keeps its gene feature
    hit_cds = [ft for ft in rec.features if ft.type in ("CDS", "mRNA") and ft.location.strand == strand
               and shared(ft) >= 30]
    hit_names = {ft.qualifiers.get("gene", [""])[0] for ft in hit_cds}
    cds_names = {ft.qualifiers.get("gene", [""])[0] for ft in rec.features if ft.type == "CDS"}
    hit, keep = [], []
    for ft in rec.features:
        name = ft.qualifiers.get("gene", [""])[0]
        if any(ft is x for x in hit_cds) or (
                ft.type == "gene" and ft.location.strand == strand and (name in hit_names or name not in cds_names)
                and shared(ft) >= 30):
            hit.append(ft)
        else:
            keep.append(ft)
    names = [f.qualifiers.get("gene", [""])[0] for f in hit if f.type == "CDS" and f.qualifiers.get("gene")]
    name = names[0] if names else m.gene
    products = [f.qualifiers["product"][0] for f in hit if f.type == "CDS" and f.qualifiers.get("product")]
    product = products[0] if products else product_for(name)
    # partial ends (gene runs off the sequence): 5' end is the low coordinate on +, the high one on -
    lo_open = partial5 if strand == 1 else partial3
    hi_open = partial3 if strand == 1 else partial5

    def pos(i: int, j: int, a: int, b: int):
        st = BeforePosition(a - 1) if (i == 0 and lo_open) else ExactPosition(a - 1)
        en = AfterPosition(b) if (j == len(exons) - 1 and hi_open) else ExactPosition(b)
        return FeatureLocation(st, en, strand=strand)
    locs = [pos(i, i, a, b) for i, (a, b) in enumerate(exons)]
    if strand == -1:
        locs = locs[::-1]
    cds_loc = CompoundLocation(locs) if len(locs) > 1 else locs[0]
    gene_loc = FeatureLocation(BeforePosition(exons[0][0] - 1) if lo_open else ExactPosition(exons[0][0] - 1),
                               AfterPosition(exons[-1][1]) if hi_open else ExactPosition(exons[-1][1]), strand=strand)
    note = "exon-intron structure predicted by protein alignment (Exonerate, MFannot method)"
    nt = cds_loc.extract(rec.seq)
    prot = str(nt[: len(nt) - len(nt) % 3].translate(table=gencode, to_stop=True))
    quals: Dict[str, List[str]] = {"gene": [name], "codon_start": ["1"], "transl_table": [str(gencode)],
                                   "product": [product], "note": [note]}
    for k, v in ({} if partial3 else stop_qualifiers(rec, cds_loc, gencode)).items():
        quals[k] = quals.get(k, []) + v
    quals["translation"] = [prot]
    new = [SeqFeature(gene_loc, type="gene", qualifiers={"gene": [name]}),
           SeqFeature(cds_loc, type="CDS", qualifiers=quals)]
    for n, info in enumerate(intron_info, 1):
        new.append(SeqFeature(FeatureLocation(info["start"] - 1, info["end"], strand=strand), type="intron",
                              qualifiers={"gene": [name], "number": [str(n)], "note": [info["note"]]}))
    rec.features = sorted(keep + new, key=lambda f: (int(f.location.start), f.type != "source", f.type != "gene"))
    return name


UNTYPED_MIN_SPLICE = 14     # untyped introns: type-I splice score needed to be written (benchmark: true >=14.5, false <=12.5)
UNTYPED_MIN_COVERAGE = 0.8  # untyped introns: share of the homolog the model must align (true 0.98, false 0.65)


def evidence(s: str, exons: List[Tuple[int, int]], kinds: List[str], splice: List[float], gencode: int,
             partial5: bool = False) -> Tuple[str, str]:
    """Confidence of a gene model, from evidence independent of existing annotations.
    high: every intron typed by Rfam and the spliced CDS has no internal stop codon;
    medium: untyped introns with splice score >= UNTYPED_MIN_SPLICE; low (not written):
    internal stop codon, invalid start codon, or untyped intron with a weak splice site."""
    table = CodonTable.unambiguous_dna_by_id[gencode]
    prot = _translate(s, exons, gencode)
    if "*" in prot.rstrip("*"):
        return "low", "internal stop codon in spliced CDS"
    if not partial5 and s[exons[0][0] - 1:exons[0][0] + 2] not in table.start_codons:
        return "low", f"no start codon ({s[exons[0][0] - 1:exons[0][0] + 2]})"
    weak = [i + 1 for i, (k, sc) in enumerate(zip(kinds, splice)) if not k and sc < UNTYPED_MIN_SPLICE]
    if weak:
        return "low", f"intron {','.join(map(str, weak))} untyped with weak splice site (< {UNTYPED_MIN_SPLICE})"
    if all(kinds):
        return "high", "all introns typed (Rfam), spliced CDS without internal stop"
    return "medium", f"untyped intron(s) with splice score >= {UNTYPED_MIN_SPLICE}, spliced CDS without internal stop"


# ── genes across the origin of a circular record ─────────────────────────────
def move_location(loc, delta: int, L: int):
    """Location shifted by delta (mod L) on a circular sequence of length L; a part that crosses the origin is split
    (join, transcript order kept), parts that meet at the old origin are merged. Fuzzy ends (< / >) are kept."""
    def pos(old, v: int):
        return type(old)(v) if isinstance(old, (BeforePosition, AfterPosition)) else ExactPosition(v)
    parts = []
    for p in loc.parts:
        s, e = int(p.start), int(p.end)
        if e - s >= L:                                  # whole record (source)
            parts.append(FeatureLocation(p.start, p.end, strand=p.strand))
            continue
        s2 = (s + delta) % L
        e2 = s2 + e - s
        if e2 <= L:
            parts.append(FeatureLocation(pos(p.start, s2), pos(p.end, e2), strand=p.strand))
            continue
        a = FeatureLocation(pos(p.start, s2), ExactPosition(L), strand=p.strand)
        b = FeatureLocation(ExactPosition(0), pos(p.end, e2 - L), strand=p.strand)
        parts += [b, a] if p.strand == -1 else [a, b]
    merged = [parts[0]]
    for q in parts[1:]:
        r = merged[-1]
        if (q.strand != -1 and int(r.end) == int(q.start)) or (q.strand == -1 and int(q.end) == int(r.start)):
            merged[-1] = (FeatureLocation(r.start, q.end, strand=r.strand) if q.strand != -1
                          else FeatureLocation(q.start, r.end, strand=r.strand))
        else:
            merged.append(q)
    return merged[0] if len(merged) == 1 else CompoundLocation(merged, operator=getattr(loc, "operator", "join"))


def rotated_copy(rec, start: int):
    """Copy of a circular record that starts at 1-based position start; one feature per feature of rec, same order."""
    from Bio.SeqRecord import SeqRecord
    L, d = len(rec.seq), start - 1
    new = SeqRecord(rec.seq[d:] + rec.seq[:d], id=rec.id, name=rec.name, description=rec.description,
                    annotations=dict(rec.annotations))
    new.features = [SeqFeature(move_location(f.location, -d, L), type=f.type,
                               qualifiers={k: list(v) for k, v in f.qualifiers.items()}) for f in rec.features]
    return new


def origin_rotation(rec, chain: List[Report]) -> Optional[int]:
    """Start (1-based) of a rotated copy on which a gene across the origin is contiguous: inside the stretch its
    chained alignments leave free between the pieces on both sides of the origin, at the position not covered by
    any feature nearest the middle of that stretch."""
    L = len(rec.seq)
    rs = sorted(chain, key=lambda r: r.qab)
    i = wrap_steps(rs)[0]
    before = [(min(e[0] for e in r.exons), max(e[1] for e in r.exons)) for r in rs[:i + 1]]
    after = [(min(e[0] for e in r.exons), max(e[1] for e in r.exons)) for r in rs[i + 1:]]
    low, high = (after, before) if rs[0].strand == 1 else (before, after)
    lo, hi = max(e for _, e in low), min(s for s, _ in high)
    if hi - lo < 2:
        return None
    covered = bytearray(L + 2)
    for f in rec.features:
        if f.type != "source":
            for p in f.location.parts:
                covered[int(p.start) + 1:int(p.end) + 1] = b"\x01" * (int(p.end) - int(p.start))
    mid = (lo + hi) // 2
    free = [p for p in range(lo + 1, hi) if not covered[p]]
    return min(free, key=lambda p: abs(p - mid)) if free else mid


def move_pos(p: int, delta: int, L: int) -> int:
    return (p - 1 + delta) % L + 1


def model_across_origin(rec, gene: str, cand: dict, chain: List[Report], lib, lib_db, a, tools, tmp) -> List[dict]:
    """A gene whose exons lie on both sides of the origin of a circular record (cox1 exon 2 at the start of a
    record rotated to start at the core's cox1 piece, exon 1 at its end): modelled on a rotated copy of the record,
    its features moved back (joins across the origin) and the replaced features removed from rec."""
    start = origin_rotation(rec, chain)
    if start is None:                   # no free stretch between the pieces: not a gene across the origin
        return [{"record": rec.id, "gene": gene, "written": "no",
                 "note": "alignments out of protein order that overlap on the sequence (paralogous hits?) - not modelled"}]
    L, d = len(rec.seq), start - 1
    rot = rotated_copy(rec, start)
    before = list(rot.features)
    rows, _ = process_record(rot, lib, lib_db, a, tools, tmp, cands={gene: cand})
    after_ids = {id(f) for f in rot.features}
    before_ids = {id(f) for f in before}
    gone = {i for i, f in enumerate(before) if id(f) not in after_ids}
    added = [SeqFeature(move_location(f.location, d, L), type=f.type,
                        qualifiers=move_position_qualifiers(f.qualifiers, L, -d))     # /transl_except moves too
             for f in rot.features if id(f) not in before_ids]
    # an intron across the origin would be a join, which NCBI rejects (table2asn SEQ_FEAT.MultiIntervalIntron): it is
    # given by the CDS join only (as NCBI's NC_027106.1 nad5)
    added = [f for f in added if not (f.type == "intron" and len(f.location.parts) > 1)]
    if gone or added:
        rec.features = sorted([f for i, f in enumerate(rec.features) if i not in gone] + added,
                              key=lambda f: (int(f.location.start), f.type != "source", f.type != "gene"))

    def span(x: str) -> str:
        p, q = (int(v) for v in x.split(".."))
        p2, q2 = move_pos(p, d, L), move_pos(q, d, L)
        return f"{p2}..{q2}" if p2 <= q2 else f"{p2}..{L},1..{q2}"
    for r in rows:
        if r.get("cds"):
            r["cds"] = ",".join(span(x) for x in r["cds"].split(","))
        for k in ("intron_start", "intron_end"):
            if r.get(k):
                r[k] = move_pos(int(r[k]), d, L)
        r["evidence"] = (r.get("evidence", "") + f"; {ACROSS_ORIGIN}: modelled on a copy starting at {start}").lstrip("; ")
        if r.get("intron_start") and r["intron_start"] > r["intron_end"] and r.get("written") == "yes":
            r["evidence"] += "; intron across the origin: written in the CDS join only (no intron feature)"
    return rows


# ── driver ───────────────────────────────────────────────────────────────────
def annotate_record(rec, lib: Dict[str, LibProt], lib_db: Path, a, tools: Dict[str, Optional[str]], tmp: Path
                    ) -> List[dict]:
    """process_record, then each gene found across the origin of the circular record on a rotated copy."""
    rows, across = process_record(rec, lib, lib_db, a, tools, tmp)
    for gene, cand, chain in across:
        rows += model_across_origin(rec, gene, cand, chain, lib, lib_db, a, tools, tmp)
    return rows


def process_record(rec, lib: Dict[str, LibProt], lib_db: Path, a, tools: Dict[str, Optional[str]],
                   tmp: Path, cands: Optional[Dict[str, dict]] = None
                   ) -> Tuple[List[dict], List[Tuple[str, dict, List[Report]]]]:
    """Rows of the report, and the genes found across the origin of a circular record (gene, candidate, chain)
    for model_across_origin. cands: model only these genes (default: prescreen + core references)."""
    seq = str(rec.seq).upper()
    L = len(seq)
    circular = rec.annotations.get("topology") == "circular"
    fna = tmp / "genome.fna"
    fna.write_text(f">g\n{seq}\n")
    own = norm_species(rec.annotations.get("organism", "")) if a.exclude_own_species else ""
    core = core_references(rec, a.refdir)
    lib = {**lib, **{f"C:{r['id']}": LibProt(r["id"], g, r["seq"]) for g, r in core.items()}}
    if cands is None:
        cands = prescreen(seq, lib, lib_db, a.genetic_code, tools, tmp, own, a.threads)
        # a gene the core annotated only in part is a candidate as well, its reference protein the Exonerate query
        # (the pre-screen's MFannot protein stays the query of its own candidates: a core reference can be a flawed
        # RefSeq protein, NC_008827.1 cox1 with WRWA inserted at the intron, which moved 6 coral junctions by 12 nt)
        cands.update({g: {"prot": f"C:{r['id']}", "strand": r["strand"], "core_only": True}
                      for g, r in core.items() if r["stub"] and g not in cands})
    rc_seq = str(Seq(seq).reverse_complement())
    rows: List[dict] = []
    models: List[Model] = []
    across: List[Tuple[str, dict, List[Report]]] = []
    for gene, c in sorted(cands.items()):
        prot = lib[c["prot"]].seq
        pf = tmp / "prot.faa"
        pf.write_text(f">{c['prot']}\n{prot}\n")
        groups = group_reports(run_exonerate(str(tools["exonerate"]), pf, fna, a.genetic_code, a.maxintron), len(prot),
                               L, circular)
        m, why = build_model(c["prot"], gene, groups[0], len(prot), circular) if groups else (None, "")
        if why == ACROSS_ORIGIN:
            across.append((gene, c, groups[0]))
            continue
        # nad4/nad5 (and other paralogues) are homologous: when the best chain overlaps another gene's CDS,
        # the next chain that does not is used (NC_006902.1 ND5); otherwise the best one is reported as low.
        if m is not None and conflicting_genes(rec, m.gene, m.strand, m.exons):
            m = next((alt for alt in (build_model(c["prot"], gene, g, len(prot), circular)[0] for g in groups[1:6])
                      if alt is not None and not conflicting_genes(rec, alt.gene, alt.strand, alt.exons)),
                     None) or realign_masked(rec, seq, c["prot"], gene, groups, len(prot), pf, tmp, a, tools) or m
        if m is not None:
            found = [m]
        else:                                    # a trans-spliced gene: its cis-spliced pieces
            found = gene_modules(c["prot"], gene, groups[0], len(prot), circular) if TRANS_SPLICED in why else []
        if not found:
            if why:
                rows.append({"record": rec.id, "gene": gene, "written": "no", "note": why})
            continue
        for m in found:
            m.core_only = bool(c.get("core_only"))
            # two models must not share coding sequence; spans may nest (a hexacoral nad5 intron hosts the cox1 gene)
            if any(overlap(x, y) for k in models if k.strand == m.strand for x in m.exons for y in k.exons):
                continue
            why = spans_record_ends(m, L, circular, a.maxintron) or \
                copied_exons(seq if m.strand == 1 else rc_seq, to_transcript(m, L))
            if why:
                rows.append({"record": rec.id, "gene": gene, "written": "no", "note": why})
                continue
            models.append(m)
    raw = {}                                   # every intron of every model, typed in one cmsearch call
    for k, m in enumerate(models):
        use = seq if m.strand == 1 else rc_seq
        tx = to_transcript(m, L)
        raw.update({f"m{k}i{j}": use[tx[j][1]:tx[j + 1][0] - 1] for j in range(len(tx) - 1)})
    typed = intron_types(tools["cmsearch"], a.cm_files, raw, a.evalue, tmp) if tools["cmsearch"] and raw else {}
    for k, m in enumerate(models):
        use = seq if m.strand == 1 else rc_seq
        tx = to_transcript(m, L)
        tinfo: List[Any] = [typed.get(f"m{k}i{j}", ("", None, "")) for j in range(len(tx) - 1)]
        for j, (kind, ev, model) in enumerate(tinfo):
            # a weak group I subtype hit is no type when the intron has near-perfect group II splice sites
            # (Endomyzostoma OR637230.1 cox1: IA1 E=0.0013, GUGUG...AC; NCBI's and the gold placement are group II)
            if kind == "I" and ev > WEAK_SUBTYPE_E and re.fullmatch(r"I[A-E]\d?", model) and \
                    max(score_type_ii(use, tx[j][1] + 1, tx[j + 1][0] - 1, x) for x in ADJUST_RANGE) >= TYPE_II_STRONG:
                tinfo[j] = ("", None, "")
        kinds = [t[0] for t in tinfo]
        refs = [lib[m.query].seq] + [v.seq for k_, v in lib.items() if v.gene == m.gene and k_ != m.query
                                     and not (a.exclude_own_species and v.species == own)][:9]
        # boundary shifts are judged against the model's homolog, the core's reference of the gene (the closest one)
        # and three more of the library: median gain (aa_score)
        near = [core[m.gene]["seq"]] if m.gene in core and core[m.gene]["seq"] != refs[0] else []
        scorer = homolog_scorer(refs[0], *near, *refs[1:4])
        adj, notes, sscores = adjust_introns(use, tx, kinds, a.genetic_code, scorer)
        n_exons = len(adj)
        adj, tinfo, notes, sscores = add_mini_exons(use, adj, tinfo, notes, sscores, a.genetic_code, refs)
        kinds = [t[0] for t in tinfo]
        if len(adj) > n_exons:
            # the introns around an inserted mini-exon get the boundary step again: exonerate had spread the mini-exon's
            # residues over the neighbouring junction (Plakinastrella NC_010217.1 cox1: exon 1 ran 9 nt into the intron)
            redo = {j for j, n in enumerate(notes) if "; mini-exon " in n}
            adj, notes2, sc2 = adjust_introns(use, adj, kinds, a.genetic_code, scorer, only=redo, w5=MINI_REFINE)
            for j in redo:
                mini = adj[j + 1] if j + 1 in redo else adj[j]          # the exon between two redone introns
                notes[j] = notes2[j] + re.sub(r"mini-exon \d+ nt", f"mini-exon {mini[1] - mini[0] + 1} nt",
                                              notes[j][notes[j].index("; mini-exon "):])
                sscores[j] = sc2[j]
        z_start, z_end = core_ends(rec, m, adj, L)
        z_start, z_end = (None if m.open5 else z_start), (None if m.open3 else z_end)
        started, partial5 = (adj, True) if m.open5 else \
            extend_to_start(use, adj, a.genetic_code, prev_feature_end(rec, m, adj, L), m.qmiss_n)
        fallback, partial3 = (started, True) if m.open3 else \
            extend_to_stop(use, started, a.genetic_code, next_feature_start(rec, m, adj, L), m.qmiss_c, rule=a.stop_rule)
        # the existing ends only when they are a start / stop: the native core ends a CDS at its last aligned codon
        # (NC_008164.1 ND5 ...GTA, the T stop 1 nt further -> table2asn NoStop; mito_curate P3 skips joins)
        if z_start is not None and use[z_start - 1:z_start + 2] not in CodonTable.unambiguous_dna_by_id[a.genetic_code].start_codons:
            z_start = None
        cand = list(adj)
        cand[0] = (z_start if z_start is not None else fallback[0][0], cand[0][1])
        if z_end is not None and not ends_with_stop(use, cand[:-1] + [(cand[-1][0], z_end)], a.genetic_code):
            z_end = None
        cand[-1] = (cand[-1][0], z_end if z_end is not None else fallback[-1][1])
        partial5 = partial5 and z_start is None
        partial3 = partial3 and z_end is None
        if valid_exons(cand):
            adj = cand
        else:
            z_start = z_end = None
            adj = fallback
        if not valid_exons(adj):
            rows.append({"record": rec.id, "gene": m.gene, "note": "invalid exon order (gene spans origin?) - not modelled"})
            continue
        ends = ("5'/3' end from existing CDS: " + ", ".join(x for x, v in (("start", z_start), ("end", z_end)) if v)
                if z_start or z_end else "")
        g_exons = to_genome(adj, m.strand, L)
        info: List[Dict[str, Any]] = []
        for j in range(len(adj) - 1):
            s_tx, e_tx = adj[j][1] + 1, adj[j + 1][0] - 1
            gs, ge = (s_tx, e_tx) if m.strand == 1 else (L + 1 - e_tx, L + 1 - s_tx)
            kind, ev, model = tinfo[j]
            label = type_label(kind, ev, model)
            hint = ("UNCONFIRMED: long intron from chained alignments, check manually"
                    if not kind and ge - gs + 1 > a.maxintron else "")       # report only, never a feature note
            info.append({"start": gs, "end": ge, "note": f"{label}; {notes[j]}; MFannot method", "hint": hint})
        if ends:
            info[0]["note"] += f"; {ends}"
        confidence, reason = evidence(use, adj, kinds, sscores, a.genetic_code, partial5)
        if m.core_only and confidence != "low" and not all(kinds):
            # exonerate against a close homolog skips divergent or frameshifted stretches with 'introns' >= 142 nt
            # (sponge NC_006894.1 nad2 / nad5): a gene the ORF pre-screen did not flag needs introns typed by Rfam
            confidence, reason = "low", "untyped intron in a gene found only from the core's partial CDS"
        if confidence != "low" and not all(kinds) and m.qcov < UNTYPED_MIN_COVERAGE:
            # Hexamermis NC_008828.1 atp6: an untyped 'intron' between parts of two gene copies, homolog 65 % aligned
            confidence, reason = "low", (f"untyped intron(s) and only {m.qcov:.0%} of the homolog aligned "
                                         f"(< {UNTYPED_MIN_COVERAGE:.0%})")
        clash = conflicting_genes(rec, m.gene, m.strand, g_exons)
        if clash:
            confidence, reason = "low", f"model overlaps existing {','.join(sorted(set(clash)))} CDS (different gene)"
        hom = lib[m.query]
        name = (apply_model(rec, m, g_exons, info, a.genetic_code, partial5, partial3)
                if confidence != "low" else m.gene)
        for n, x in enumerate(info, 1):
            rows.append({"record": rec.id, "gene": name, "strand": "+" if m.strand == 1 else "-",
                         "cds": ",".join(f"{s}..{e}" for s, e in g_exons), "intron": n,
                         "intron_start": x["start"], "intron_end": x["end"], "intron_len": x["end"] - x["start"] + 1,
                         "confidence": confidence, "written": "no" if confidence == "low" else "yes",
                         "evidence": reason + (f"; {x['hint']}" if x.get("hint") else "") + (f"; {m.piece}" if m.piece else ""),
                         "note": x["note"],
                         "homolog": f"{hom.species} {hom.gene}",
                         "homolog_coverage": f"{m.qcov:.2f}", "exonerate_score": m.score})
    return rows, across


_JOB: Dict[str, Any] = {}               # shared state for process_file (set in main; inherited by fork workers)


def process_file(f: Path) -> List[dict]:
    """Annotate every record of one GenBank file; a failing record is reported and left unchanged."""
    lib, lib_db, a, tools, out = _JOB["lib"], _JOB["lib_db"], _JOB["a"], _JOB["tools"], _JOB["out"]
    rows: List[dict] = []
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        try:
            recs = list(SeqIO.parse(str(f), "genbank"))
        except Exception as exc:
            return [{"record": f.name, "written": "no", "note": f"cannot read GenBank file: {exc}"}]
        for rec in recs:
            try:
                rows += annotate_record(rec, lib, lib_db, a, tools, tmp)
            except Exception as exc:              # keep going: record is written back unchanged
                rows.append({"record": rec.id, "written": "no", "note": f"error, record left unchanged: {exc!r}"})
        src_root = Path(a.input) if Path(a.input).is_dir() else f.parent
        dest = out / f.relative_to(src_root)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if any(r.get("written") == "yes" for r in rows):
            write_gb(recs, dest)
        elif dest.resolve() != f.resolve():
            shutil.copy2(f, dest)
    return rows


def find_tool(name: str, bin_dir: Optional[str]) -> Optional[str]:
    return mito_env.find(name, bin_dir)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="GenBank file or directory (e.g. MitoM 04.final_gb)")
    p.add_argument("-o", "--out", required=True, help="Output directory (may equal input dir to update in place)")
    p.add_argument("--lib", action="append", default=[],
                   help=f"Protein library in MFannot header format (default {DEFAULT_LIB}); repeatable")
    p.add_argument("--ref_gb", action="append", default=[],
                   help="GenBank file/dir of close relatives: their CDS translations are added to the "
                        "homolog library (helps lineages poorly covered by MFannot's library, e.g. annelids)")
    p.add_argument("--genetic_code", type=int, default=5)
    p.add_argument("--maxintron", type=int, default=MAX_INTRON, help=f"Exonerate --maxintron (MFannot default {MAX_INTRON})")
    p.add_argument("--evalue", type=float, default=0.01, help="cmsearch E-value cut-off for intron type")
    p.add_argument("--cm", default=str(DEFAULT_CM), help="Rfam group I/II intron covariance models")
    p.add_argument("--no_subtype_cms", action="store_true",
                   help="Do not use the group I subtype models of Selles Vidal et al. 2025 (Rfam only)")
    p.add_argument("--bin", default=None, help="Directory with exonerate, blastp, makeblastdb, cmsearch")
    p.add_argument("--refdir", default=None,
                   help="MITOS2 reference data (refseq89m) holding the reference proteins named in the annotation "
                        "core's CDS notes, used as closest homologs; default from mito_env")
    p.add_argument("--threads", type=int, default=1, help="Threads for blastp")
    p.add_argument("--jobs", type=int, default=1, help="GenBank files processed in parallel")
    p.add_argument("--timeout", type=int, default=900, help="Seconds allowed per external program call")
    p.add_argument("--stop_rule", choices=["homolog", "first"], default="first",
                   help="3' end when no existing CDS end is available: 'homolog' = stop where the homolog ends "
                        "(complete or T/TA stop there); 'first' = first in-frame complete stop codon")
    p.add_argument("--exclude_own_species", action="store_true",
                   help="Ignore library proteins of the record's own species (benchmarking)")
    a = p.parse_args(argv)

    tools = {t: find_tool(t, a.bin) for t in ("exonerate", "blastp", "makeblastdb", "cmsearch")}
    missing = [t for t in ("exonerate", "blastp", "makeblastdb") if not tools[t]]
    if missing:
        print(f"[mito_intron] missing {', '.join(missing)} "
              "(conda install -c bioconda exonerate blast infernal; or --bin)", file=sys.stderr)
        return 2
    if not tools["cmsearch"]:
        print("[mito_intron] cmsearch not found: intron type will not be identified", file=sys.stderr)
    a.cm_files = [Path(a.cm)] + ([] if a.no_subtype_cms else subtype_cms())
    rd = Path(a.refdir or mito_env.mitos_refdir(None) or "").expanduser()
    rd = rd / "refseq89m" if (rd / "refseq89m").is_dir() else rd
    a.refdir = rd if (rd / "featureProt").is_dir() else None
    lib = load_library(a.lib or [str(DEFAULT_LIB)])
    if a.ref_gb:
        ref = load_ref_gb(a.ref_gb, a.genetic_code)
        lib.update(ref)
        print(f"[mito_intron] + {len(ref)} reference proteins from --ref_gb", file=sys.stderr)
    files = iter_files([a.input], GB_EXT)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows: List[dict] = []
    global TIMEOUT
    TIMEOUT = a.timeout
    with tempfile.TemporaryDirectory() as td:
        lib_faa = Path(td) / "lib.faa"
        lib_faa.write_text("".join(f">{k}\n{v.seq}\n" for k, v in lib.items()))
        subprocess.run([str(tools["makeblastdb"]), "-in", str(lib_faa), "-dbtype", "prot", "-out", str(Path(td) / "lib")],
                       capture_output=True, check=True)
        _JOB.update(lib=lib, lib_db=Path(td) / "lib", a=a, tools=tools, out=out)
        if a.jobs > 1 and len(files) > 1:
            import multiprocessing as mp
            with mp.get_context("fork").Pool(a.jobs) as pool:      # workers inherit _JOB (no pickling)
                for r in pool.imap(process_file, files):
                    rows += r
        else:
            for f in files:
                rows += process_file(f)
    cols = ["record", "gene", "strand", "cds", "intron", "intron_start", "intron_end", "intron_len",
            "confidence", "written", "evidence", "note", "homolog", "homolog_coverage", "exonerate_score"]
    write_tsv(out / "intron_report.tsv", cols, rows)
    w_ = [r for r in rows if r.get("written") == "yes"]
    low = [r for r in rows if r.get("written") == "no"]
    print(f"[mito_intron] {len(files)} file(s), library {len(lib)} proteins; "
          f"{len(w_)} intron(s) written in {len({(r['record'], r['gene']) for r in w_})} gene(s); "
          f"{len(low)} low-confidence not written -> {out / 'intron_report.tsv'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
