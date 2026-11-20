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
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from Bio import SeqIO
import mito_env  # noqa: E402  (finds external programs: config, PATH, conda environments)
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
import mito_common  # noqa: E402
from mito_common import cli_guard, count_seq_records, prepare_input, require_records  # noqa: E402
import mito_intron_exon  # noqa: E402
import mito_intron_rules  # noqa: E402,F401  (tests patch its functions)
from mito_intron_rules import (  # noqa: E402, F401  (re-exported: aa_score, refine_junction, REFINE_MIN_GAIN, TYPE_II_MIN)
    MAX_INTRON,
    MIN_INTRON,
    MIN_REPORT_SCORE,
    REFINE_MIN_GAIN,
    TYPE_II_MIN,
    _translate,
    aa_score,
    add_mini_exons,
    adjust_introns,
    ADJUST_RANGE,
    copied_exons,
    core_ends,
    ends_with_stop,
    extend_to_start,
    extend_to_stop,
    homolog_scorer,
    junction_continuity,
    splicing_gain,
    Model,
    next_feature_start,
    overlap,
    prev_feature_end,
    refine_junction,
    score_type_ii,
    to_genome,
    valid_exons,
)
from mito_common import run_tool as run_program  # noqa: E402  (CC0: logged, e.g. by `mitom run --debug`)
MINI_REFINE = 12                     # refine_junction's 5' window for the introns around an inserted mini-exon
WEAK_SUBTYPE_E, TYPE_II_STRONG = 1e-3, 55   # a group I subtype hit weaker than this yields to group II sites this strong
# MFannot defaults (PirModels/Option.pir, mfannot:284)
MIN_ORF_AA, BLAST_EVALUE, SHORT_HSP_EVALUE = 40, 1e-10, 1e-11
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


TIMEOUT = 900                          # seconds per external program call: default of --timeout (passed down)


_FAILED: set = set()                   # programs already reported as failed (one warning each)


def try_program(cmd: List[str], what: str, timeout: int = TIMEOUT) -> str:
    """Run an external program and return its stdout; on timeout or failure return '' and warn (once per program), so
    one hard sequence cannot stall or kill a batch: the step then finds no intron there (mito_common.run_tool: the
    programs whose failure stops the step)."""
    try:
        p = run_program(cmd, timeout=timeout)
        if p.returncode != 0 and what not in _FAILED:     # (2026-10-31 review: failures were silent: no introns)
            _FAILED.add(what)
            print(f"[mito_intron] WARNING: {what} failed (exit {p.returncode}): "
                  f"{(p.stderr or '').strip()[-200:]} - its results are missing", file=sys.stderr)
        return "" if p.returncode else p.stdout          # a failed program's partial output is not parsed (Kimi 1.13)
    except subprocess.TimeoutExpired:
        print(f"[mito_intron] {what} timed out after {timeout}s - skipped", file=sys.stderr)
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


def internal_stops(rec, gencode: int) -> List[Tuple[str, int, str]]:
    """Standard protein genes whose CDS, as annotated, has an internal stop codon: (library gene name, strand,
    translation with stops). An intron read as coding sequence shows as stops in frame (an annotation engine without
    introns, an imported annotation); so does a frameshift. Unknown genes and pseudo CDS are skipped."""
    canon = {p.stem.lower(): p.stem for p in DEFAULT_LIB.glob("*.faa")}
    out = []
    for f in rec.features:
        g = norm_gene(f.qualifiers.get("gene", [""])[0]) if f.type == "CDS" else ""
        if g not in PRODUCTS or "pseudo" in f.qualifiers:
            continue
        cs = int(f.qualifiers.get("codon_start", ["1"])[0]) - 1
        nt = str(f.location.extract(rec.seq)).upper()[cs:]
        prot = str(Seq(nt[: len(nt) - len(nt) % 3]).translate(table=gencode))
        if "*" in prot.rstrip("*"):
            out.append((canon.get(g.lower(), g), f.location.strand or 1, prot))
    return out


def closest_lib_protein(lib: Dict[str, "LibProt"], gene: str, prot: str, n: int = 60) -> Optional[str]:
    """The library protein of `gene` that aligns best to prot (stops as X; first n of the gene, BLOSUM62 local)."""
    ids = [k for k, v in lib.items() if v.gene == gene][:n]
    if not ids:
        return None
    score = homolog_scorer(*[lib[k].seq for k in ids])(prot.replace("*", "X"))
    return max(zip(score, ids))[1]


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
              tmp: Path, exclude_species: str = "", threads: int = 1, timeout: int = TIMEOUT) -> Dict[str, dict]:
    """Candidate intron-containing genes: {gene: {'prot': lib id, 'strand': 1|-1|2}}."""
    orfs = find_orfs(seq, gencode)
    if not orfs:
        return {}
    q = tmp / "orfs.faa"
    q.write_text("".join(f">O{i}\n{p}\n" for i, (_, p) in enumerate(orfs)))
    res = try_program([str(tools["blastp"]), "-query", str(q), "-db", str(lib_db), "-matrix", "BLOSUM62",
                       "-outfmt", "6 qseqid sseqid evalue pident qstart qend", "-num_threads", str(threads)], "blastp",
                      timeout)
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


def run_exonerate(exe: str, prot_faa: Path, genome_fna: Path, gencode: int, maxintron: int,
                  timeout: int = TIMEOUT) -> List[Report]:
    from mito_core import exonerate_code       # codes exonerate 2.4 lacks (24, 33): the 64-letter table (D2, 2026-11-05)
    cmd = [exe, "--model", "protein2genome", "--showalignment", "no", "--showvulgar", "no",
           "--showtargetgff", "yes", "--minintron", str(MIN_INTRON), "--maxintron", str(maxintron),
           "--intronpenalty", "-35", "--geneticcode", exonerate_code(int(gencode)),
           *(["--splice3", str(SPLICE_PSSM), "--splice5", str(SPLICE_PSSM)] if SPLICE_PSSM.exists() else []),
           # %tS = target strand (%g, the splice-derived gene orientation, is "." without introns)
           "--ryo", ">>RYO\t%qab\t%qae\t%tab\t%tae\t%s\t%tS\n", str(prot_faa), str(genome_fna)]
    out = try_program(cmd, "exonerate", timeout)
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
    failed = SUBTYPE_CM_DIR / ".failed"                    # offline / rate limited: try again after a day, not every
    recent = failed.is_file() and time.time() - failed.stat().st_mtime < 86400   # run (code review, 2026-11-06)
    if not (SUBTYPE_CM_DIR / ".complete").is_file() and not recent:   # interrupted: fetched again (Kimi 1.14)
        try:
            import json
            import urllib.request
            SUBTYPE_CM_DIR.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(SUBTYPE_CM_URL, timeout=30) as r:
                items = json.load(r)
            for it in items:
                name = Path(str(it["name"])).name      # a name from the web never leaves the folder (N7)
                if name == it["name"] and name.endswith(".cm"):
                    part = SUBTYPE_CM_DIR / f".{name}.part"
                    with urllib.request.urlopen(it["download_url"], timeout=60) as r, open(part, "wb") as h:
                        shutil.copyfileobj(r, h)
                    part.replace(SUBTYPE_CM_DIR / name)
            (SUBTYPE_CM_DIR / ".complete").touch()
            failed.unlink(missing_ok=True)
        except Exception as exc:                       # offline: carry on with Rfam only
            print(f"[mito_intron] group I subtype models not available ({exc}); using Rfam only (next try in a "
                  f"day; delete {failed} to try now)", file=sys.stderr)
            try:
                failed.touch()
            except OSError:
                pass
    # only a complete set: an interrupted download's first models alone would type introns on part of the
    # subtypes while the message says Rfam only (Bugbot on PR #4, 2026-11-06)
    return sorted(SUBTYPE_CM_DIR.glob("*.cm")) if (SUBTYPE_CM_DIR / ".complete").is_file() else []


def intron_types(cmsearch: str, cms: List[Path], seqs: Dict[str, str], evalue: float, tmp: Path,
                 timeout: int = TIMEOUT) -> Dict[str, Tuple[str, float, str]]:
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
        try_program([cmsearch, "--noali", "--toponly", *(["--rfam"] if fast else []), "--cpu", "0", "--mxsize", "128",
                     "--tblout", str(tbl), str(cm), str(fna)], "cmsearch", timeout)
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


def conflicting_cds(rec, gene: str, strand: int, exons: List[Tuple[int, int]]) -> list:
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
            out.append(f)
    return out


def short_exon_model(rec, L: int, gene: str, c: dict, lib: Dict[str, LibProt], grp: List[Report], m: Optional[Model],
                     why: str, a, own: str) -> Optional[Model]:
    """The chain of Exonerate alignments plus the short divergent exons it does not align (mito_intron_exon: homolog
    vote against a shuffled null, issue #15). A model with the new exons, or None when the chain has no hole of
    >= 12 residues or no significant locus. Only single-strand chains in protein order are completed."""
    if why or len({r.strand for r in grp}) != 1 or wrap_steps(grp):
        return None
    strand = grp[0].strand
    seq = str(rec.seq).upper()
    use = seq if strand == 1 else str(Seq(seq).reverse_complement())
    prot = lib[c["prot"]].seq
    homs = [v.seq for k_, v in lib.items() if v.gene == gene and k_ != c["prot"] and not k_.startswith("C:")
            and not (a.exclude_own_species and v.species == own)]
    blocked = []
    for f in rec.features:
        if f.type in ("CDS", "rRNA", "tRNA") and not (f.qualifiers.get("gene") and
                                                       norm_gene(f.qualifiers["gene"][0]) == norm_gene(gene)):
            for p_ in f.location.parts:
                blocked.append((int(p_.start), int(p_.end)) if strand == 1 else (L - int(p_.end), L - int(p_.start)))

    def tx(e):
        return e if strand == 1 else (L + 1 - e[1], L + 1 - e[0])
    reps = sorted(((r.qab, r.qae, min(tx(e)[0] for e in r.exons), max(tx(e)[1] for e in r.exons)) for r in grp))
    own_exons = sorted(tx(e) for r in grp for e in r.exons)
    new, note = mito_intron_exon.complete_chain(use, own_exons, reps, prot, homs, a.genetic_code, blocked, a.maxintron,
                                                MIN_INTRON)
    if not note:
        return None
    exons = sorted(tx(e) for e in new)
    first, last = min(e[0] for e in own_exons), max(e[1] for e in own_exons)
    n_hole = any(e[1] < first for e in new)                      # in transcript order
    c_hole = any(e[0] > last for e in new)
    qn, qc = (m.qmiss_n, m.qmiss_c) if m else (min(r.qab for r in grp), len(prot) - max(r.qae for r in grp))
    out = Model(c["prot"], gene, strand, exons, sum(r.score for r in grp), 1.0, 0 if c_hole else qc,
                0 if n_hole else qn)
    return out


def conflicting_genes(rec, gene: str, strand: int, exons: List[Tuple[int, int]]) -> List[str]:
    """Gene names of conflicting_cds."""
    return [f.qualifiers["gene"][0] for f in conflicting_cds(rec, gene, strand, exons)]


WEAK_CONFLICT_RATIO = 3          # a model whose exonerate score is >= 3 x the score of the CDS it overlaps wins over it
WEAK_CONFLICT_INSIDE = 0.8       # ... when >= 80 % of that CDS lies inside the model's exons (the same DNA, explained twice)


def outweighed_conflicts(rec, m, exons: List[Tuple[int, int]]) -> Optional[list]:
    """The conflicting CDS of model m when every one is a weak call of the core inside m's exons (note 'ref ... score N',
    N * WEAK_CONFLICT_RATIO <= m.score, WEAK_CONFLICT_INSIDE of its length in the exons); None when a conflict stands.
    MW592986.1 atp6: score 528, 161 + 68 aa on exons 16007..16488 / 17095..17314, against an ATP8 of score 84 inside
    exon 2 (the real atp8 is elsewhere) - the strong continuous gene is better evidence than the weak small one.
    A CDS without a core score (the submitter's annotation) is never outweighed."""
    found = conflicting_cds(rec, m.gene, m.strand, exons)
    if not found:
        return None
    for f in found:
        sc = re.search(r"\bscore (\d+(?:\.\d+)?)", " ".join(f.qualifiers.get("note", [])))
        if not sc or float(sc.group(1)) * WEAK_CONFLICT_RATIO > m.score:
            return None
        parts = [(int(p.start) + 1, int(p.end)) for p in f.location.parts]
        inside = sum(overlap(x, y) for x in exons for y in parts)
        if inside < WEAK_CONFLICT_INSIDE * sum(b - a + 1 for a, b in parts):
            return None
    return found


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
    reps = run_exonerate(str(tools["exonerate"]), pf, fna, a.genetic_code, max(a.maxintron, longest + 200), getattr(a, "timeout", TIMEOUT))
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
                and shared(ft) >= 30) or (
                ft.type == "intron" and ft.location.strand == strand and (name in hit_names or not name)
                and exons[0][0] - 1 <= int(ft.location.start) and int(ft.location.end) <= exons[-1][1]):
            hit.append(ft)              # an intron of the replaced gene (record annotated already): no second copy
        else:
            keep.append(ft)
    names = [f.qualifiers.get("gene", [""])[0] for f in hit if f.type == "CDS" and f.qualifiers.get("gene")]
    names = [x for x in names if norm_gene(x) == norm_gene(m.gene)] or names      # (not the name of a replaced weak call)
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
    if not partial5 and str(nt[:3]).upper() in CodonTable.unambiguous_dna_by_id[gencode].start_codons:
        prot = "M" + prot[1:]                     # NCBI: M for any start codon, as mito_curate (Kimi 3.3, 2026-11-05)
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


INTRON_NOTE = "intron across the origin of the record"   # CDS note of a dropped intron (see model_across_origin)
UNTYPED_MIN_SPLICE = 14     # untyped introns: type-I splice score needed to be written (benchmark: true >=14.5, false <=12.5)
UNTYPED_MIN_COVERAGE = 0.8  # untyped introns: share of the homolog the model must align (true 0.98, false 0.65)
# --untyped_rule continuity: an intron removes DNA only, so the homologs of the gene run on across its junction
JUNCTION_SPAN = 0.5         # more than half of the homologs align on both sides of the junction
JUNCTION_MAX_AA = 1         # homolog residues missing / added at the junction (median): 1 = the codon the junction splits
SPLICE_GAIN_MIN = 20        # BLOSUM62 points (~9 bits) splicing must add over reading on in place (random sequence: < 0)
CONT_WINDOW = 30            # nt each end of such a junction may move to make the homologs continuous
WRITE_ALL = bool(os.environ.get("MITOM_INTRON_WRITE_ALL"))   # benchmarks only: write every model without stop / clash


def continuous(c: Tuple[float, float, float]) -> bool:
    span, miss, extra = c
    return span > JUNCTION_SPAN and miss <= JUNCTION_MAX_AA and extra <= JUNCTION_MAX_AA


def continuity_text(c: Tuple[float, float, float]) -> str:
    return f"{c[0]:.0%} of homologs span it, {c[1]:g} aa missing, {c[2]:g} aa added"


def evidence(s: str, exons: List[Tuple[int, int]], kinds: List[str], splice: List[float], gencode: int,
             partial5: bool = False, cont: Optional[List[Tuple[float, float, float]]] = None,
             rule: str = "splice", gains: Optional[List[float]] = None) -> Tuple[str, str]:
    """Confidence of a gene model, from evidence independent of existing annotations.
    high: every intron typed by Rfam and the spliced CDS has no internal stop codon;
    medium: the untyped introns pass `rule`; low (not written): internal stop codon, invalid start codon, or an
    untyped intron that fails `rule`.
    rule 'splice' (before 2026-11-09): an untyped intron needs a type-I splice score >= UNTYPED_MIN_SPLICE;
    rule 'continuity': an untyped intron needs (1) the homologs continuous across its junction (cont,
    junction_continuity: most homologs span it, none of their residues missing and none added, JUNCTION_*) - an intron
    removes DNA only, whatever its type (degenerate group II introns of annelids and bryozoans have no group I splice
    sites) - and (2) to be needed (gains, splicing_gain >= SPLICE_GAIN_MIN): read on in place without splicing, the gene
    does not cover its homolog as well."""
    table = CodonTable.unambiguous_dna_by_id[gencode]
    prot = _translate(s, exons, gencode)
    if "*" in prot.rstrip("*"):
        return "low", "internal stop codon in spliced CDS"
    if not partial5 and s[exons[0][0] - 1:exons[0][0] + 2] not in table.start_codons:
        return "low", f"no start codon ({s[exons[0][0] - 1:exons[0][0] + 2]})"
    if all(kinds):
        return "high", "all introns typed (Rfam), spliced CDS without internal stop"
    if WRITE_ALL:
        return "medium", "untyped intron(s) (MITOM_INTRON_WRITE_ALL: no evidence rule)"
    if rule == "splice":
        weak = [i + 1 for i, (k, sc) in enumerate(zip(kinds, splice)) if not k and sc < UNTYPED_MIN_SPLICE]
        if weak:
            return "low", f"intron {','.join(map(str, weak))} untyped with weak splice site (< {UNTYPED_MIN_SPLICE})"
        return "medium", f"untyped intron(s) with splice score >= {UNTYPED_MIN_SPLICE}, spliced CDS without internal stop"
    broken = [f"{i + 1} ({continuity_text(c)})" for i, (k, c) in enumerate(zip(kinds, cont or [])) if not k and not continuous(c)]
    if broken:
        return "low", f"untyped intron {', '.join(broken)}: homologs not continuous across the junction"
    idle = [f"{i + 1} (gain {g:g})" for i, (k, g) in enumerate(zip(kinds, gains or [])) if not k and g < SPLICE_GAIN_MIN]
    if idle:
        return "low", (f"untyped intron {', '.join(idle)} not needed: the gene reads on in place without splicing "
                       f"(homolog score gain of splicing < {SPLICE_GAIN_MIN})")
    return "medium", ("untyped intron(s) placed by protein homology (homologs continuous across each junction; "
                      "type not identified), spliced CDS without internal stop")


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
    # given by the CDS join only (as NCBI's NC_027106.1 nad5). Its description stays on the CDS (INTRON_NOTE): when
    # the record is rotated to the gene's first exon (mitom_orient.reorient_after_introns) the intron lies inside the
    # record again and gets its feature back (MJ20260820041 CHY cox1 group II intron: lost before)
    crossing = [f for f in added if f.type == "intron" and len(f.location.parts) > 1]
    for f in crossing:
        cds = next((c for c in added if c.type == "CDS" and c.qualifiers.get("gene") == f.qualifiers.get("gene")), None)
        if cds is not None:
            cds.qualifiers.setdefault("note", []).append(
                f"{INTRON_NOTE} {f.qualifiers.get('number', ['?'])[0]}: {'; '.join(f.qualifiers.get('note', []))}")
    added = [f for f in added if f not in crossing]
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
    return rows + stop_check_rows(rec, a.genetic_code)


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
        cands = prescreen(seq, lib, lib_db, a.genetic_code, tools, tmp, own, a.threads, getattr(a, "timeout", TIMEOUT))
        # a gene the core annotated only in part is a candidate as well, its reference protein the Exonerate query
        # (the pre-screen's MFannot protein stays the query of its own candidates: a core reference can be a flawed
        # RefSeq protein, NC_008827.1 cox1 with WRWA inserted at the intron, which moved 6 coral junctions by 12 nt)
        cands.update({g: {"prot": f"C:{r['id']}", "strand": r["strand"], "core_only": True}
                      for g, r in core.items() if r["stub"] and g not in cands})
        # and a gene whose CDS has an internal stop codon (an intron read as coding sequence), its query the core's
        # reference or the closest library protein
        for g, strand, prot_ in internal_stops(rec, a.genetic_code):
            pid = f"C:{core[g]['id']}" if g in core else closest_lib_protein(lib, g, prot_)
            if pid and g not in cands:
                cands[g] = {"prot": pid, "strand": strand, "core_only": g in core}
    rc_seq = str(Seq(seq).reverse_complement())
    rows: List[dict] = []
    models: List[Model] = []
    across: List[Tuple[str, dict, List[Report]]] = []
    for gene, c in sorted(cands.items()):
        prot = lib[c["prot"]].seq
        pf = tmp / "prot.faa"
        pf.write_text(f">{c['prot']}\n{prot}\n")
        groups = group_reports(run_exonerate(str(tools["exonerate"]), pf, fna, a.genetic_code, a.maxintron, getattr(a, "timeout", TIMEOUT)), len(prot),
                               L, circular)
        m, why = build_model(c["prot"], gene, groups[0], len(prot), circular) if groups else (None, "")
        if why == ACROSS_ORIGIN:
            across.append((gene, c, groups[0]))
            continue
        if groups and getattr(a, "short_exons", False) and m is None and not why:
            m2 = short_exon_model(rec, L, gene, c, lib, groups[0], m, why if m is None else "", a, own)
            if m2 is not None:
                m, why = m2, ""
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
    typed = intron_types(tools["cmsearch"], a.cm_files, raw, a.evalue, tmp, getattr(a, "timeout", TIMEOUT)) if tools["cmsearch"] and raw else {}
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
        homs = [*near, *refs]
        def cont_of(ex):
            pr = _translate(use, ex, a.genetic_code)
            return [junction_continuity(pr, homs, sum(b - a_ + 1 for a_, b in ex[:j + 1]) // 3) for j in range(len(ex) - 1)]
        cont = cont_of(adj)
        if getattr(a, "untyped_rule", "continuity") == "continuity":
            # an untyped junction where the homologs are not continuous is searched for the placement that makes them
            # continuous (both ends independently, frame kept); none -> the model stays as it is (and is not written)
            for j, k in enumerate(kinds):
                if k or continuous(cont[j]):
                    continue
                new, gain = refine_junction(use, adj, j, a.genetic_code, scorer, w5=CONT_WINDOW, w3=CONT_WINDOW)
                if gain and continuous(cont_of(new)[j]):
                    adj = new
                    cont = cont_of(adj)
                    notes[j] += "; junction moved to where the homologs are continuous"
        need = [splicing_gain(use, adj, j, a.genetic_code, scorer) for j in range(len(adj) - 1)]
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
        for x, (sp, mi, ex), g in zip(info, cont, need):
            x.update(span=sp, missing=mi, extra=ex, gain=g)
        confidence, reason = evidence(use, adj, kinds, sscores, a.genetic_code, partial5, cont,
                                      getattr(a, "untyped_rule", "continuity"), need)
        if m.core_only and confidence != "low" and not all(kinds) and getattr(a, "untyped_rule", "continuity") == "splice" and not WRITE_ALL:
            # exonerate against a close homolog skips divergent or frameshifted stretches with 'introns' >= 142 nt
            # (sponge NC_006894.1 nad2 / nad5): a gene the ORF pre-screen did not flag needs introns typed by Rfam
            confidence, reason = "low", "untyped intron in a gene found only from the core's partial CDS"
        if confidence != "low" and not all(kinds) and m.qcov < UNTYPED_MIN_COVERAGE and not WRITE_ALL and \
                getattr(a, "untyped_rule", "continuity") == "splice":
            # (rule 'continuity' needs no coverage: a divergent homolog aligns in part; a false junction shows as
            # discontinuity, which the coverage only stood in for)
            # Hexamermis NC_008828.1 atp6: an untyped 'intron' between parts of two gene copies, homolog 65 % aligned
            confidence, reason = "low", (f"untyped intron(s) and only {m.qcov:.0%} of the homolog aligned "
                                         f"(< {UNTYPED_MIN_COVERAGE:.0%})")
        clash = conflicting_genes(rec, m.gene, m.strand, g_exons)
        if clash:
            weak = outweighed_conflicts(rec, m, g_exons) if confidence != "low" else None
            if weak is not None:             # apply_model replaces them: the strong gene explains that DNA
                reason += (f"; replaces the weak {','.join(sorted(set(clash)))} call "
                           f"({', '.join(f.qualifiers['note'][0] for f in weak)})")
            else:
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
                         "homolog_coverage": f"{m.qcov:.2f}", "exonerate_score": m.score,
                         "junction_span": f"{x['span']:.2f}", "junction_missing_aa": x["missing"],
                         "junction_extra_aa": x["extra"], "splice_score": sscores[n - 1],
                         "splicing_gain": f"{x['gain']:g}"})
    return rows, across


def stop_check_rows(rec, gencode: int) -> List[dict]:
    """After the intron step: genes whose CDS still has an internal stop codon (an intron not found, or a frameshift)."""
    return [{"record": rec.id, "gene": g, "strand": "+" if st == 1 else "-", "written": "no",
             "note": "CDS has internal stop codon(s) after the intron step: wrong frame or location, a frameshift, "
                     "or an intron not annotated - check this gene (the curation after this step may still correct "
                     "the frame: curation_report.tsv, rule P7; the final word is the NCBI check)"}
            for g, st, _ in internal_stops(rec, gencode)]


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
        src_root = Path(_JOB["src"]) if Path(_JOB["src"]).is_dir() else f.parent   # the prepared input (Kimi 1.4, 2026-11-05)
        dest = out / f.relative_to(src_root)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if any(r.get("written") == "yes" for r in rows):
            write_gb(recs, dest)
        elif dest.resolve() != f.resolve():
            shutil.copy2(f, dest)
    return rows


def find_tool(name: str, bin_dir: Optional[str]) -> Optional[str]:
    return mito_env.find(name, bin_dir)


def build_parser() -> argparse.ArgumentParser:
    """The command line (main); options() gives the same options without one."""
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
    p.add_argument("--untyped_rule", choices=["continuity", "splice"], default="continuity",
                   help="when an intron Rfam cannot type is written: continuity = the homologs of the gene are "
                        "continuous across its junction (an intron removes DNA only); splice = it has a type-I splice "
                        "score >= %d (before 2026-11-09)" % UNTYPED_MIN_SPLICE)
    p.add_argument("--short_exons", action="store_true",
                   help="EXPERIMENTAL: search short divergent exons Exonerate leaves unaligned, for genes with no "
                        "model (homolog vote against a shuffled null, issue #15); gained nothing on the NCBI set")
    p.add_argument("--stop_rule", choices=["homolog", "first"], default="first",
                   help="3' end when no existing CDS end is available: 'homolog' = stop where the homolog ends "
                        "(complete or T/TA stop there); 'first' = first in-frame complete stop codon")
    p.add_argument("--exclude_own_species", action="store_true",
                   help="Ignore library proteins of the record's own species (benchmarking)")
    return p


def options(**given) -> argparse.Namespace:
    """Options for run(): input and out are required, e.g. options(input="gb/", out="introns/", genetic_code=4)."""
    return mito_common.options(build_parser(), **given)


def run(a: argparse.Namespace) -> int:
    """Introns de novo (MFannot method) in the GenBank file / folder a.input -> a.out (intron_report.tsv). Returns the
    exit code (1: exonerate / BLAST missing)."""
    tools = {t: find_tool(t, a.bin) for t in ("exonerate", "blastp", "makeblastdb", "cmsearch")}
    missing = [t for t in ("exonerate", "blastp", "makeblastdb") if not tools[t]]
    if missing:
        print(f"[mito_intron] missing {', '.join(missing)} "
              "(conda install -c bioconda exonerate blast infernal; or --bin)", file=sys.stderr)
        return 1
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
    src_in = prepare_input(a.input)
    files = iter_files([src_in], GB_EXT)
    require_records(len(files) and count_seq_records(src_in), "GenBank records", a.input)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows: List[dict] = []
    with tempfile.TemporaryDirectory() as td:
        lib_faa = Path(td) / "lib.faa"
        lib_faa.write_text("".join(f">{k}\n{v.seq}\n" for k, v in lib.items()))
        run_program([str(tools["makeblastdb"]), "-in", str(lib_faa), "-dbtype", "prot", "-out", str(Path(td) / "lib")],
                    check=True)
        _JOB.update(lib=lib, lib_db=Path(td) / "lib", a=a, tools=tools, out=out, src=src_in)
        if a.jobs > 1 and len(files) > 1:
            import multiprocessing as mp
            with mp.get_context("fork").Pool(a.jobs) as pool:      # workers inherit _JOB (no pickling)
                for r in pool.imap(process_file, files):
                    rows += r
        else:
            for f in files:
                rows += process_file(f)
    cols = ["record", "gene", "strand", "cds", "intron", "intron_start", "intron_end", "intron_len",
            "confidence", "written", "evidence", "note", "homolog", "homolog_coverage", "exonerate_score",
            "junction_span", "junction_missing_aa", "junction_extra_aa", "splice_score", "splicing_gain"]
    write_tsv(out / "intron_report.tsv", cols, rows)
    w_ = [r for r in rows if r.get("written") == "yes"]
    low = [r for r in rows if r.get("written") == "no"]
    print(f"[mito_intron] {len(files)} file(s), library {len(lib)} proteins; "
          f"{len(w_)} intron(s) written in {len({(r['record'], r['gene']) for r in w_})} gene(s); "
          f"{len(low)} low-confidence not written -> {out / 'intron_report.tsv'}", file=sys.stderr)
    return 0


def main(argv=None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    sys.exit(cli_guard(main))
