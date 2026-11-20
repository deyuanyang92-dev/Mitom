#!/usr/bin/env python3
"""mito_pcg_hmm.py - protein genes the annotation missed, found with a profile HMM of the nearest relatives.

Short, fast-evolving genes (atp8 above all, also atp6, nad6, nad4l) are missed when no close relative is among the
reference proteins: a pairwise search of a ~50-codon protein against distant single proteins scores too low (e.g.
atp8 of the tube worms Serpulidae, which NCBI itself annotates in only 2 of 9 records). A profile HMM built from the
gene in the record's nearest relatives knows which positions are conserved and finds such genes (as HMMER
searches with gene profiles do in MITOS / MFannot-style annotation).

For each standard protein gene absent from a record:
  1. references: the gene's proteins in models/pcg_refs.faa.gz (scripts/build_pcg_refs.py: RefSeq and INSDC
     animal mitogenomes, with lineages) nearest the record's lineage (most leading taxonomy names shared; among
     equal ones the most similar in a conserved gene of the record, cox1 else cob ...: 3-mer overlap), up to --k
     (30); profiles of the nearest 10, 30 and 100 (MAFFT, hmmbuild) - which one finds a ~50-codon gene depends
     on the references, so all three are tried and each E-value is multiplied by 3 (Bonferroni);
  2. targets: the open reading frames in the sequence no exon covers (mito_orf.candidates: tandem repeats masked;
     >= 20 codons), translated with the record's code; hmmsearch, HMMER's own inclusion threshold E <= 0.01 (for
     the record: the expected number of chance hits among its ORFs) times the prior odds of the gene among the 30
     nearest annotated relatives, min(1, p / (1 - p)) with p = (relatives with the gene + 0.5) / (relatives + 1) -
     a gene the lineage has lost (atp8 of most flatworms and nematodes) needs far stronger evidence; and >= 40 % of
     the model aligned (as the rRNA HMM step, mito_rrna_hmm.MIN_COVER);
  3. the gene: the frame of the best hit, from the start codon nearest the start of the alignment (upstream first)
     to its stop codon - mitochondrial genes overlap their neighbours by a few codons (atp8 / atp6), so up to
     OVERLAP nt of the stop may lie on the next gene; a frame that runs on into a tRNA ends there (tRNA
     punctuation, as mito_orf). The start is refined later by the curation (mito_curate.py).
Written as gene + CDS with the core's name and product, /note "found by MitoM mito_pcg_hmm: ..."; every gene
searched is a row of pcg_hmm_report.tsv (found or not, E-value, references used).

Usage: python3 mito_pcg_hmm.py search -i <GenBank folder or file> [--refs models/pcg_refs.faa.gz] [--k 30]
       [--genes atp8,atp6] [--exclude_acc A,B] [--exclude_genus Hydroides] [--report pcg_hmm_report.tsv]
(`mitom run` / `annotate` / `import`: step 6, after the MITOS2 rescue; --pcg_hmm no: off). Measured:
benchmark/results/pcg_hmm_2026-11-10/README.md.

HMMER: Eddy 2011, PLoS Comput Biol 7:e1002195; MAFFT: Katoh & Standley 2013, Mol Biol Evol 30:772.
"""
from __future__ import annotations

import argparse
import gzip
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, ExactPosition, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import mito_orf
from mito_common import cli_guard, iter_gb_files, require_path, run_tool, stop_qualifiers, write_gb, write_tsv
from mito_names import CORE_NAMES, PCG, gene_name

HERE = Path(__file__).resolve().parent
REFS = HERE / "models" / "pcg_refs.faa.gz"
K = 30                      # references for the prior
SIZES = (10, 30, 100)       # profiles of the 10 / 30 / 100 nearest references: a ~50-codon gene is found by one
                            # or another depending on the references (serpulid atp8); E x 3 for the three tries
E_MAX = 0.01                # HMMER's inclusion threshold
MIN_COVER = 0.4             # share of the model aligned (as mito_rrna_hmm.MIN_COVER)
MIN_CODONS = 20             # shortest ORF searched (atp8 of some phyla: ~30 codons)
OVERLAP = 30                # nt of the stop codon's ORF end that may lie on the next gene (atp8 / atp6: 7-10)
NOTE = "found by MitoM mito_pcg_hmm"
COLUMNS = ["record", "gene", "status", "start", "end", "strand", "codons", "evalue", "model_cover", "references",
           "nearest_reference", "shared_lineage", "prior", "note"]


def lineage_of(rec: SeqRecord) -> List[str]:
    return [t.replace(" ", "_") for t in rec.annotations.get("taxonomy", [])] + \
        [rec.annotations.get("organism", "").replace(" ", "_")]


def load_refs(path: Path = REFS) -> Dict[str, List[Tuple[str, List[str], str]]]:
    """gene -> [(accession, lineage, protein)]; also "_records" -> [(accession, lineage, genes)] for the prior."""
    out: Dict[str, list] = defaultdict(list)
    recs: Dict[str, Tuple[List[str], set]] = {}
    with gzip.open(path, "rt") as fh:
        for r in SeqIO.parse(fh, "fasta"):
            acc, gene, _code, lin = r.description.split("|", 3)
            lineage = lin.split(";")
            out[gene].append((acc, lineage, str(r.seq)))
            recs.setdefault(acc, (lineage, set()))[1].add(gene)
    out["_records"] = [(acc, lin, genes) for acc, (lin, genes) in recs.items() if len(genes) >= 10]
    return out


def prior(refs: Dict[str, list], gene: str, lineage: List[str], exclude_acc: set, exclude_genus: set,
          k: int = K) -> Tuple[float, int, int]:
    """(p, with, of): the share of the k nearest annotated records (>= 10 genes in the library) that have the gene,
    smoothed (with + 0.5) / (of + 1). Genes lost in a lineage (atp8 of most flatworms and nematodes) have p near 0."""
    near = nearest(refs.get("_records", []), lineage, k, exclude_acc, exclude_genus)
    have = sum(gene in genes for _, _, genes, _ in near)
    return (have + 0.5) / (len(near) + 1), have, len(near)


def shared(a: Sequence[str], b: Sequence[str]) -> int:
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


MARKERS = ("cox1", "cob", "nad1", "nad5", "nad4", "cox3")   # conserved genes that rank references of equal lineage


def kmers(prot: str, k: int = 3) -> set:
    return {prot[i: i + k] for i in range(len(prot) - k + 1)}


def similarity(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if a and b else 0.0


def nearest(refs: List[Tuple[str, List[str], str]], lineage: List[str], k: int = K, exclude_acc: set = frozenset(),
            exclude_genus: set = frozenset(), marker: Optional[Tuple[set, Dict[str, set]]] = None
            ) -> List[Tuple[str, List[str], str, int]]:
    """The k references sharing the longest leading part of the lineage; among equal ones the most similar in a
    conserved gene of the record (marker: its 3-mers and those of the references' same gene - Sabellidae before
    Siboglinidae for a serpulid, both in Sabellida), then accession order. A record without taxonomy (FASTA input):
    by the marker alone."""
    pool = [(acc, lin, prot, shared(lineage, lin)) for acc, lin, prot in refs
            if acc not in exclude_acc and acc.split(".")[0] not in exclude_acc
            and (lin[-1].split("_")[0] if lin else "") not in exclude_genus]
    if marker:
        q, of = marker
        pool.sort(key=lambda x: (-x[3], -similarity(q, of.get(x[0], set())), x[0]))
    else:
        pool.sort(key=lambda x: (-x[3], x[0]))
    return pool[:k]


def marker_of(rec: SeqRecord, refs: Dict[str, list], code: int) -> Optional[Tuple[set, Dict[str, set]]]:
    """(3-mers of the record's first marker gene, {reference accession: 3-mers of its same gene})."""
    for g in MARKERS:
        cds = [f for f in rec.features if f.type == "CDS" and gene_name(f.qualifiers, "CDS") == g]
        if not cds or g not in refs:
            continue
        try:
            nt = cds[0].extract(rec.seq)
            prot = str(nt[: len(nt) // 3 * 3].translate(table=code)).rstrip("*")
        except Exception:
            continue
        cache = refs.setdefault(f"_kmers_{g}", {})
        if not cache:
            cache.update({acc: kmers(p) for acc, _, p in refs[g]})
        return kmers(prot), cache
    return None


def build_hmm(name: str, prots: List[str], td: Path, bin_dir: Optional[str]) -> Path:
    fa, aln, hmm = td / f"{name}.faa", td / f"{name}.aln", td / f"{name}.hmm"
    fa.write_text("".join(f">r{i}\n{p}\n" for i, p in enumerate(prots)))
    mafft = mito_orf.find_tool("mafft", bin_dir)
    hmmbuild = mito_orf.find_tool("hmmbuild", bin_dir)
    if not (mafft and hmmbuild):
        raise RuntimeError("mafft / hmmbuild not found (conda install -c bioconda mafft hmmer)")
    if len(prots) > 1:
        p = run_tool([mafft, "--quiet", "--auto", "--amino", str(fa)], check=True)
        aln.write_text(p.stdout)
    else:
        aln.write_text(fa.read_text())
    run_tool([hmmbuild, "--amino", "-n", name, str(hmm), str(aln)], check=True)
    return hmm


def model_len(hmm: Path) -> int:
    for ln in hmm.read_text().splitlines():
        if ln.startswith("LENG"):
            return int(ln.split()[1])
    return 0


def hmmsearch(hmm: Path, targets: List[str], td: Path, bin_dir: Optional[str]) -> List[dict]:
    """Domain hits: target index, E-value of the sequence, model and target alignment ranges."""
    exe = mito_orf.find_tool("hmmsearch", bin_dir)
    if not exe:
        raise RuntimeError("hmmsearch not found (conda install -c bioconda hmmer)")
    fa, tbl = td / "t.faa", td / "t.dom"
    fa.write_text("".join(f">{i}\n{t}\n" for i, t in enumerate(targets)))
    run_tool([exe, "--max", "--noali", "--domtblout", str(tbl), "-E", "10", str(hmm), str(fa)], check=True)
    hits = []
    for ln in tbl.read_text().splitlines():
        if ln.startswith("#"):
            continue
        c = ln.split()
        hits.append({"target": int(c[0]), "evalue": float(c[6]), "hmm_from": int(c[15]), "hmm_to": int(c[16]),
                     "ali_from": int(c[17]), "ali_to": int(c[18])})
    return hits


def place(rec: SeqRecord, o: "mito_orf.Orf", hit: dict, code: int) -> Optional[Tuple[int, int, int, bool, str]]:
    """(start, end, strand, 5' open, tRNA end) of the gene in the frame of ORF candidate o: the start codon nearest
    the start of the alignment (upstream first, else up to 10 codons downstream, else the stretch's first codon with
    an open 5' end), the end at the stop codon of the frame (stop within OVERLAP nt of a covered stretch allowed) or
    where the frame runs into a tRNA."""
    n = len(rec.seq)
    seq = str(rec.seq).upper()
    s, e, p0 = o.seg
    strand = o.strand
    nt = seq[s:] + seq[: e - n] if e > n else seq[s:e]
    if strand == -1:
        from Bio.Seq import Seq
        nt = str(Seq(nt).reverse_complement())
    k = (len(nt) - 3) // 3                                 # codons before the closing stop
    starts = set(CodonTable.unambiguous_dna_by_id[code].start_codons)
    i0 = p0 + hit["ali_from"] - 1                           # codon of the segment where the alignment starts
    up = [j for j in range(i0, p0 - 1, -1) if nt[3 * j: 3 * j + 3] in starts]
    down = [j for j in range(i0 + 1, min(i0 + 11, k)) if nt[3 * j: 3 * j + 3] in starts]
    open5 = False
    if up:                                                 # (upstream: the start must not cut the aligned part)
        a = up[0]
    elif down:
        a = down[0]
    else:
        a, open5 = p0, True
    cov = mito_orf.covered_mask(rec)
    reading = [(s + i) % n if strand == 1 else (e - 1 - i) % n for i in range(len(nt))]
    first_cov = next((i for i in range(3 * a, 3 * k) if cov[reading[i]]), None)
    trna = ""
    if first_cov is None or 3 * k + 3 - first_cov <= OVERLAP:
        b = 3 * k + 3                                       # to the stop codon
    else:
        b = 3 * a + (first_cov - 3 * a) // 3 * 3          # runs into a gene / tRNA: last whole codon before it
        trna = o.trna_end or "next gene"
    if (b - 3 * a) // 3 < MIN_CODONS:
        return None
    lo, hi = (s + 3 * a, s + b) if strand == 1 else (e - b, e - 3 * a)
    if lo >= n:                                            # (the copy beyond the origin of a circular record)
        lo, hi = lo - n, hi - n
    return lo, hi, strand, open5, trna


def location(lo: int, hi: int, strand: int, open5: bool, open3: bool, n: int):
    lo_open = open5 if strand == 1 else open3
    hi_open = open3 if strand == 1 else open5
    lp = BeforePosition(lo % n) if lo_open else ExactPosition(lo % n)
    if hi <= n:
        return FeatureLocation(lp, AfterPosition(hi) if hi_open else ExactPosition(hi), strand=strand)
    parts = [FeatureLocation(lp, n, strand=strand),
             FeatureLocation(0, AfterPosition(hi - n) if hi_open else ExactPosition(hi - n), strand=strand)]
    return CompoundLocation(parts if strand == 1 else parts[::-1])


def missing_genes(rec: SeqRecord) -> List[str]:
    have = {gene_name(f.qualifiers, f.type) for f in rec.features if f.type in ("CDS", "gene")}
    return [g for g in sorted(PCG) if g not in have]


def overlaps(lo: int, hi: int, a: int, b: int, n: int) -> bool:
    """Whether the intervals lo..hi and a..b share a position on a circle of n nt (an end past n runs across the
    origin: place() gives hi > n for such a gene; 2026-11-20 review: it was compared as a line)."""
    return any(lo < b + shift and a + shift < hi for shift in (-n, 0, n))


def translation_of(rec: SeqRecord, cds: SeqFeature, code: int, partial5: bool) -> str:
    """/translation of a CDS made here ('' when it has an internal stop): codons of its location, a start codon of the
    table giving M unless the 5' end is open, the terminal stop dropped (a T / TA stop has no complete codon)."""
    nt = str(cds.extract(rec.seq)).upper()
    table = CodonTable.unambiguous_dna_by_id[code]
    nt = nt[: len(nt) - len(nt) % 3]
    prot = "".join(table.forward_table.get(nt[i:i + 3], "*" if nt[i:i + 3] in table.stop_codons else "X")
                   for i in range(0, len(nt), 3))
    prot = prot[:-1] if prot.endswith("*") else prot
    if not prot or "*" in prot:
        return ""
    if not partial5 and nt[:3] in table.start_codons:
        prot = "M" + prot[1:]
    return prot


def search_record(rec: SeqRecord, refs: Dict[str, list], default_code: int, genes: Optional[Sequence[str]] = None,
                  k: int = K, exclude_acc: set = frozenset(), exclude_genus: set = frozenset(),
                  bin_dir: Optional[str] = None, write: bool = True) -> List[dict]:
    """Search the record's missing genes (or `genes`); add the found ones (write); report rows."""
    code = mito_orf.code_of(rec, default_code)
    todo = [g for g in (genes or missing_genes(rec)) if g in refs]
    rows = []
    if not todo:
        return rows
    cands = [o for o in mito_orf.candidates(rec, code, MIN_CODONS)]
    lineage = lineage_of(rec)
    marker = marker_of(rec, refs, code)
    taken: List[Tuple[int, int]] = []
    with tempfile.TemporaryDirectory() as tdn:
        td = Path(tdn)
        for g in todo:
            near = nearest(refs[g], lineage, max(SIZES), exclude_acc | {rec.id}, exclude_genus, marker)
            base = {"record": rec.id, "gene": g, "references": len(near),
                    "nearest_reference": near[0][0] if near else "",
                    "shared_lineage": ";".join(lineage[: near[0][3]][-2:]) if near else ""}
            if not near or not cands:
                rows.append({**base, "status": "not found", "note": "no references" if not near else "no ORF"})
                continue
            pr, have, of = prior(refs, g, lineage, exclude_acc | {rec.id}, exclude_genus, k)
            e_max = E_MAX * min(1.0, pr / (1 - pr))      # prior odds: a gene the relatives lack needs more evidence
            base["prior"] = f"{have}/{of}"
            sizes = sorted({min(n, len(near)) for n in SIZES})
            raw = []                                       # one profile per size; E x the profiles tried (Bonferroni)
            for n_ref in sizes:
                hmm = build_hmm(f"{g}_{n_ref}", [p for _, _, p, _ in near[:n_ref]], td, bin_dir)
                leng_n = model_len(hmm)
                for h in hmmsearch(hmm, [o.piece for o in cands], td, bin_dir):
                    h.update(evalue=h["evalue"] * len(sizes), leng=leng_n, n_ref=n_ref)
                    raw.append(h)
            hits = [h for h in raw if h["evalue"] <= e_max and (h["hmm_to"] - h["hmm_from"] + 1) >= MIN_COVER * h["leng"]]
            hits.sort(key=lambda h: h["evalue"])
            done = False
            for h in hits:
                o = cands[h["target"]]
                pl = place(rec, o, h, code)
                if pl is None:
                    continue
                lo, hi, strand, open5, trna = pl
                if any(overlaps(lo, hi, a, b, len(rec.seq)) for a, b in taken):
                    continue
                cover = (h["hmm_to"] - h["hmm_from"] + 1) / h["leng"]
                note = (f"{NOTE}: profile HMM of the {h['n_ref']} nearest references ({base['shared_lineage']}), "
                        f"hmmsearch E={h['evalue']:.0e} (x {len(sizes)} profiles tried), {100 * cover:.0f}% of the "
                        f"model; check")
                if trna:
                    note += f"; no stop codon before {trna}: 3' end open"
                if write:
                    sym, product = CORE_NAMES[g]
                    loc = location(lo, hi, strand, open5, bool(trna), len(rec.seq))
                    quals = {"gene": [sym], "codon_start": ["1"], "transl_table": [str(code)], "product": [product],
                             "note": [note]}
                    for q, v in stop_qualifiers(rec, loc, code).items():
                        quals[q] = quals.get(q, []) + v
                    cds = SeqFeature(loc, type="CDS", qualifiers=quals)
                    tr = translation_of(rec, cds, code, open5)
                    if tr:
                        cds.qualifiers["translation"] = [tr]          # (run alone, or with --curate no, the record keeps it)
                    rec.features.append(SeqFeature(loc, type="gene", qualifiers={"gene": [sym]}))
                    rec.features.append(cds)
                taken.append((lo, hi))
                rows.append({**base, "status": "found", "start": lo + 1, "end": hi, "strand": "+" if strand == 1 else "-",
                             "codons": (hi - lo) // 3, "evalue": f"{h['evalue']:.1e}", "model_cover": f"{cover:.2f}",
                             "note": note})
                done = True
                break
            if not done:
                best = min((h["evalue"] for h in raw), default=None)
                note = f"best hmmsearch E={best:.0e}" if best is not None else "no hit"
                if best is not None and best <= E_MAX and e_max < E_MAX:
                    note += (f": below 0.01 but {g} is in only {have} of the {of} nearest annotated relatives "
                             f"(threshold {e_max:.0e}) - check")
                rows.append({**base, "status": "not found", "note": note})
    if write and any(r["status"] == "found" for r in rows):
        rec.features.sort(key=lambda f: (f.type != "source", int(f.location.parts[0].start) if f.location else 0,
                                          f.type != "gene"))
    return rows


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("search", help="search the missing protein genes of GenBank records, in place")
    s.add_argument("-i", "--input", required=True, help="GenBank file or folder (written in place)")
    s.add_argument("--refs", default=str(REFS))
    s.add_argument("--genetic_code", type=int, default=5, help="for records without /transl_table (default 5)")
    s.add_argument("--k", type=int, default=K, help=f"annotated relatives that give the prior of a gene being present (default {K}); the profiles use "
                                                       f"the {SIZES} nearest references")
    s.add_argument("--genes", default="", help="comma-separated genes to search (default: the missing ones)")
    s.add_argument("--exclude_acc", default="", help="accessions left out of the references (benchmarks; also "
                   "$MITOM_EXCLUDE_REF)")
    s.add_argument("--exclude_genus", default="", help="genera left out of the references (benchmarks; also "
                   "$MITOM_EXCLUDE_GENUS)")
    s.add_argument("--report", default=None, help="report TSV (default <input folder>/pcg_hmm_report.tsv)")
    s.add_argument("--bin", default=None, help="folder with mafft, hmmbuild, hmmsearch")
    return p


def main(argv: Sequence[str] = None) -> int:
    a = build_parser().parse_args(argv)
    src = require_path(a.input)
    files = [src] if src.is_file() else iter_gb_files(src)
    refs = load_refs(Path(a.refs))
    genes = [g for g in a.genes.split(",") if g] or None
    import os                                       # (benchmarks: the record and its genus left out of the references)
    ex_acc = {x for x in (a.exclude_acc + "," + os.environ.get("MITOM_EXCLUDE_REF", "")).split(",") if x}
    ex_gen = {x for x in (a.exclude_genus + "," + os.environ.get("MITOM_EXCLUDE_GENUS", "")).split(",") if x}
    rows = []
    for f in files:
        recs = list(SeqIO.parse(str(f), "genbank"))
        changed = False
        for rec in recs:
            try:
                rr = search_record(rec, refs, a.genetic_code, genes, a.k, ex_acc, ex_gen, a.bin)
            except Exception as exc:                        # one record's failure is a row, not the end of the run
                rr = [{"record": rec.id, "gene": "", "status": "failed", "note": f"{type(exc).__name__}: {exc}"}]
            rows += rr
            changed |= any(r["status"] == "found" for r in rr)
        if changed:
            write_gb(recs, f)
    report = Path(a.report) if a.report else (src.parent if src.is_file() else src) / "pcg_hmm_report.tsv"
    write_tsv(report, COLUMNS, rows)
    n = sum(r["status"] == "found" for r in rows)
    print(f"[mito_pcg_hmm] {len(rows)} missing protein genes searched, {n} found -> {report}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(cli_guard(main))
