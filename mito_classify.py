#!/usr/bin/env python3
"""mito_classify.py - phylum and genetic code of a mitochondrial genome without a taxonomy (FASTA input, GenBank
records without a lineage): --clade auto of `mitom run` uses it, so that no --clade / --genetic_code is needed.

Method (a few seconds per genome):
  1. blastx of the sequence (code 5: AGR = Ser, ATA = Met, TGA = Trp - no early stops in any animal mitochondrial
     code) against 8 conserved proteins of the MITOS2 reference data refseq89m (cox1 cox2 cox3 cob atp6 nad1 nad4
     nad5; ~8,000 genomes of all animal phyla).
  2. The phylum: the vote of the best hits (top 5 per gene, weighted by bit score relative to the gene's best hit),
     the lineage of each reference genome from models/ref_taxa.tsv.gz (scripts/build_ref_taxa.py). The genetic
     code: the vote of the CDS /transl_table of the best hits of that phylum (so a tunicate gets 13, not the
     vertebrate 2).
  3. Codon check: in the alignments to the closest references, the reference amino acid at the genome's AGA AGG
     ATA AAA TGA TAA codons; every animal mitochondrial code (2 4 5 9 13 14 21 24 33) is scored by how often the
     reference residue is the amino acid the code gives that codon; codes reading AGA / AGG as stop are left out
     when they occur inside the aligned genes. A code that fits clearly better than the one of the hits is used
     (the genome's own codons over its relatives' label; a WARNING in the run).
  4. Phylum support below MIN_SUPPORT (no close relative in the reference: small phyla): phylum uncertain (clade
     Animal), genetic code from the codons.
  $MITOM_EXCLUDE_REF (accessions, as for mito_core.py): those reference genomes are not counted (evaluation).

Command line (test / inspection):  python mito_classify.py genome.fasta [more.fasta ...] [--threads 8]
Output: one line per sequence: id, phylum, genetic code, support, codon check, closest reference.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from Bio.Data import CodonTable
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import mito_env  # noqa: E402
from mito_taxa import METAZOAN_PHYLA, PHYLUM_SYNONYMS  # noqa: E402

REF_TAXA = HERE / "models" / "ref_taxa.tsv.gz"
GENES = ("cox1", "cox2", "cox3", "cob", "atp6", "nad1", "nad4", "nad5")
TOP_PER_GENE = 5
MIN_BITS = 50.0
ANIMAL_MT_CODES = (2, 4, 5, 9, 13, 14, 21, 24, 33)
CHECK_CODONS = ("AGA", "AGG", "ATA", "AAA", "TGA", "TAA")
MIN_CODON_OBS = 5           # a codon seen fewer times in the alignments is not scored
CODON_MARGIN = 0.15         # a code fits "clearly better" when its mean agreement is this much higher
NOISE_FOLD = 3              # a codon must be seen > NOISE_FOLD x the internal TAG (noise) count to be scored
AGR_SENSE = 0.001          # AGA + AGG per aligned codon above which they are not stop codons
MIN_SUPPORT = 0.7           # phylum vote below this: no close relative in the reference (small phyla) - the phylum is
                            # a guess (24 / 59 right in the leave-relatives-out test, 80 / 81 above it): clade Animal,
                            # and the genetic code from the codons (69 / 71 right) instead of the hits


@dataclass
class Classification:
    record: str
    phylum: str = ""                 # '' = no decision (too few hits)
    genetic_code: int = 0
    lineage: str = ""                # the lineage shared by the best hits of the phylum
    support: float = 0.0             # weight of the phylum among all hits (0-1)
    closest: str = ""                # best reference: accession species
    codon_code: int = 0              # code that fits the codons best (0 = codons not informative)
    codon_check: str = ""            # consistent / conflict: ... / not informative
    n_genes: int = 0                 # genes with hits
    note: str = ""
    phylum_guess: str = ""           # the phylum vote when its support is below MIN_SUPPORT (phylum then '')
    code_from: str = "closest reference genomes"
    codon_scores: Dict[int, float] = field(default_factory=dict)


@lru_cache(maxsize=1)
def ref_taxa() -> Dict[str, Tuple[str, str]]:
    """accession -> (genetic code, lineage) of the reference genomes."""
    out: Dict[str, Tuple[str, str]] = {}
    with gzip.open(REF_TAXA, "rt") as fh:
        next(fh)
        for line in fh:
            acc, code, lin = line.rstrip("\n").split("\t")
            out[acc] = (code, lin)
    return out


@lru_cache(maxsize=1)
def _codes_by_taxon() -> Dict[str, Counter]:
    out: Dict[str, Counter] = defaultdict(Counter)
    for code, lin in ref_taxa().values():
        if code:
            for t in lin.split(";"):
                out[t][int(code)] += 1
    return out


def code_of_lineage(lineage: str) -> int:
    """The genetic code of the reference genomes of the deepest taxon of `lineage` that has any (a GenBank record
    with a taxonomy but no CDS /transl_table: Tunicata 13, Vertebrata 2 ...); 0 when none or no phylum."""
    idx = _codes_by_taxon()
    tax = [x for x in lineage.split(";") if x]
    ph = phylum_of(lineage)
    names = [t for t in tax if t == ph or PHYLUM_SYNONYMS.get(t) == ph]
    if not names:                                 # only the phylum and the taxa below it (not Metazoa, Bilateria)
        return 0
    for t in reversed(tax[tax.index(names[0]):]):
        if t in idx:
            return idx[t].most_common(1)[0][0]
    return 0


def phylum_of(lineage: str) -> str:
    tax = lineage.split(";")
    return next((t for t in tax if t in METAZOAN_PHYLA), "") or next((PHYLUM_SYNONYMS[t] for t in tax
                                                                      if t in PHYLUM_SYNONYMS), "")


def _blastx(blastx: str, query: Path, db: Path, threads: int) -> List[List[str]]:
    p = subprocess.run([blastx, "-query", str(query), "-db", str(db), "-query_gencode", "5", "-evalue", "1e-5",
                        "-seg", "no", "-max_target_seqs", "100", "-num_threads", str(max(1, threads)),
                        "-outfmt", "6 sseqid bitscore qstart qend qseq sseq stitle"],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"blastx failed on {db.name}: {p.stderr.strip()[-300:]}")
    return [c.split("\t") for c in p.stdout.splitlines() if c]


def _codons(genome: str, qstart: int, qend: int, qseq: str, sseq: str) -> Iterable[Tuple[str, str]]:
    """(genome codon, reference residue) of each aligned column of a blastx HSP."""
    if qstart <= qend:
        nt, pos = genome, qstart - 1
    else:                                       # minus strand: walk the reverse complement
        nt, pos = str(Seq(genome).reverse_complement()), len(genome) - qstart
    for q, s in zip(qseq, sseq):
        if q == "-":
            continue
        codon = nt[pos:pos + 3]
        pos += 3
        if s != "-" and len(codon) == 3:
            yield codon, s


def codon_scores(obs: Dict[str, Counter], n_codons: int = 0) -> Dict[int, float]:
    """Mean over the informative codons of the share of reference residues that are the code's amino acid. A code that
    reads AGA / AGG as stop (2) is left out when they occur inside the aligned genes (>= AGR_SENSE per codon:
    lancelets 3.0 / 1000, vertebrates 0-0.1 / 1000). A codon counts only when it is seen clearly more often than TAG
    (a stop in every code: its count inside the alignments is the noise of sequencing errors and frameshifts)."""
    scores = {}
    obs = defaultdict(Counter, obs)
    noise = sum(obs["TAG"].values())          # TAG is a stop in every animal code: inside genes = errors / frameshifts
    agr = sum(obs["AGA"].values()) + sum(obs["AGG"].values())
    agr_sense = n_codons and agr / n_codons >= AGR_SENSE and agr > NOISE_FOLD * noise + 2
    for code in ANIMAL_MT_CODES:
        t = CodonTable.unambiguous_dna_by_id[code]
        if agr_sense and "AGA" in t.stop_codons:
            continue
        vals = []
        for c in CHECK_CODONS:
            n = sum(obs[c].values())
            if n >= MIN_CODON_OBS + NOISE_FOLD * noise:     # (2026-10-31 stress test: random 10 % substitutions put
                                                            # AGR / TAA into vertebrate genes: codes 13 / 33 won)
                aa = t.forward_table.get(c, "*")
                vals.append(obs[c][aa] / n if aa != "*" else 0.0)
        if vals:
            scores[code] = sum(vals) / len(vals)
    return scores


def classify_seq(rid: str, genome: str, refdir: Path, blastx: str, threads: int = 4,
                 exclude: Iterable[str] = ()) -> Classification:
    """Phylum and genetic code of one genome (see the module docstring); `exclude`: taxon names (a genus, a family
    ...) whose reference genomes are not counted (evaluation: the genome and its relatives are in the reference)."""
    exclude = set(exclude)
    exclude_acc = {x.split(".")[0] for x in re.split(r"[,\s]+", os.environ.get("MITOM_EXCLUDE_REF", "")) if x}
    res = Classification(rid)
    genome = genome.upper()
    taxa = ref_taxa()
    with tempfile.TemporaryDirectory(prefix="mito_classify.") as td:
        q = Path(td) / "q.fa"
        q.write_text(f">q\n{genome}\n")
        dbs = [refdir / "featureProt" / f"{g}.fas" for g in GENES]
        with cf.ThreadPoolExecutor(max_workers=min(len(dbs), max(1, threads))) as ex:
            per_gene = list(ex.map(lambda d: _blastx(blastx, q, d, 1), dbs))
    votes: Counter = Counter()
    best_hits: List[Tuple[float, str, str, List[str]]] = []      # (weight, phylum, accession, blast columns)
    for rows in per_gene:
        seen, kept = set(), []
        for c in sorted(rows, key=lambda c: -float(c[1])):
            acc = c[0].split(":")[0]
            species = c[6].split(" ", 1)[1] if " " in c[6] else ""
            if float(c[1]) < MIN_BITS or acc in seen or acc not in taxa or acc in exclude_acc:
                continue
            if exclude and (species.split(" ")[0] in exclude or exclude & set(taxa[acc][1].split(";"))):
                continue
            seen.add(acc)
            kept.append((c, acc))
            if len(kept) == TOP_PER_GENE:
                break
        if not kept:
            continue
        res.n_genes += 1
        top = float(kept[0][0][1])
        for c, acc in kept:
            w = float(c[1]) / top
            ph = phylum_of(taxa[acc][1]) or "unknown"
            votes[ph] += w
            best_hits.append((w * float(c[1]), ph, acc, c))
    if not votes or res.n_genes < 2:
        res.note = f"too few protein hits ({res.n_genes} gene(s))"
        return res
    phylum, w = votes.most_common(1)[0]
    res.support = round(w / sum(votes.values()), 3)
    if phylum == "unknown":
        res.note = "the best hits have no phylum in their lineage"
        return res
    res.phylum = phylum
    mine = sorted((h for h in best_hits if h[1] == phylum), key=lambda h: -h[0])
    code_votes: Counter = Counter()
    for wt, _, acc, _ in mine:
        if taxa[acc][0]:
            code_votes[int(taxa[acc][0])] += wt
    res.genetic_code = code_votes.most_common(1)[0][0] if code_votes else 0
    lins = [taxa[acc][1].split(";") for _, _, acc, _ in mine[:5]]
    shared = []
    for parts in zip(*lins):
        if len(set(parts)) != 1:
            break
        shared.append(parts[0])
    res.lineage = ";".join(shared)
    top = mine[0][3]
    res.closest = f"{top[0].split(':')[0]} {top[6].split(' ', 1)[1] if ' ' in top[6] else ''}".strip()
    # codon check on the alignments to the references of the chosen code
    # (reference residues are homologous amino acids whatever the reference's own code: all hits of the phylum)
    obs: Dict[str, Counter] = defaultdict(Counter)
    n_codons = 0
    for _, _, acc, c in mine:
        for codon, aa in _codons(genome, int(c[2]), int(c[3]), c[4], c[5]):
            n_codons += 1
            if codon in CHECK_CODONS or codon == "TAG":
                obs[codon][aa] += 1
    res.codon_scores = {k: round(v, 3) for k, v in codon_scores(obs, n_codons).items()}
    if res.codon_scores and res.genetic_code not in res.codon_scores:      # the hits' code reads AGR as stop
        res.codon_code = max(res.codon_scores, key=res.codon_scores.get)
        res.codon_check = (f"conflict: AGA / AGG inside the genes, not stop codons as in code {res.genetic_code}; "
                           f"code {res.codon_code} fits the codons")
    elif not res.codon_scores:
        res.codon_check = "not informative"
    else:
        best = max(res.codon_scores, key=lambda k: (res.codon_scores[k], k == res.genetic_code))
        res.codon_code = best
        if best == res.genetic_code or res.codon_scores[best] - res.codon_scores[res.genetic_code] < CODON_MARGIN:
            res.codon_check = "consistent"
        else:
            res.codon_check = (f"conflict: code {best} fits the codons better "
                               f"({res.codon_scores[best]:.2f} vs {res.codon_scores[res.genetic_code]:.2f})")
    if res.codon_check.startswith("conflict"):          # the genome's own codons over the relatives' label
        res.note = f"genetic code {res.codon_code} from the codons, not {res.genetic_code} of the closest genomes"
        res.genetic_code, res.code_from = res.codon_code, "codons"
    if res.support < MIN_SUPPORT:
        res.note = (f"no close relative in the reference data: phylum uncertain (best guess {res.phylum}, support "
                    f"{res.support})")
        res.phylum_guess, res.phylum = res.phylum, ""
        if res.codon_code:
            res.genetic_code, res.code_from = res.codon_code, "codons"
    return res


def classify_records(records: List[Tuple[str, str]], refdir: Optional[str] = None, threads: int = 8,
                     blastx: Optional[str] = None) -> List[Classification]:
    """(id, sequence) records -> their Classification, in parallel."""
    rd = Path(refdir or mito_env.mitos_refdir(None) or "").expanduser()
    if (rd / "refseq89m").is_dir():
        rd = rd / "refseq89m"
    if not (rd / "featureProt").is_dir():
        raise FileNotFoundError(f"reference data refseq89m not found ({rd}): mitom setup")
    bx = blastx or mito_env.find("blastx")
    if not bx:
        raise FileNotFoundError("blastx not found (mitom doctor)")
    per = max(1, min(len(GENES), threads))
    with cf.ThreadPoolExecutor(max_workers=max(1, threads // per)) as ex:
        return list(ex.map(lambda r: classify_seq(r[0], r[1], rd, bx, per), records))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", nargs="+")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--refdir", default=None)
    a = ap.parse_args()
    from Bio import SeqIO
    recs = [(r.id, str(r.seq)) for f in a.fasta for r in SeqIO.parse(f, "fasta")]
    print("record\tphylum\tgenetic_code\tsupport\tcodon_check\tclosest\tlineage\tnote")
    for c in classify_records(recs, a.refdir, a.threads):
        print(f"{c.record}\t{c.phylum or '-'}\t{c.genetic_code or '-'}\t{c.support}\t{c.codon_check}\t{c.closest}\t"
              f"{c.lineage}\t{c.note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
