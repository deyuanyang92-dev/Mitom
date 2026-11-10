#!/usr/bin/env python3
"""mito_rna.py - transcript evidence for an annotated mitochondrial genome from RNA-seq reads.

Animal mitochondrial genomes are transcribed as long polycistronic precursors that are cut at the tRNAs (tRNA
punctuation, Ojala et al. 1981); mature mRNAs get poly(A) tails, which complete the incomplete stop codons T / TA.
Coverage is therefore continuous (the precursors fill the gaps) and cannot show where genes end. What can:
  - relative depth: the mean depth of a gene or ORF over the median of the 13 protein genes (Selifanova et al. 2023
    judged the Polydora ORFans transcribed by depth like nad4l / nad2's);
  - poly(A) sites: reads whose clipped 3' end begins with a run of A (or whose clipped 5' end on the reverse strand
    ends with T): the 3' end of a mature transcript; sites at tRNA 3' ends are CCA tails and are not counted; a site
    whose last INTERNAL_WINDOW nt (transcript sense) hold an A-run >= INTERNAL_RUN or >= INTERNAL_A A's is oligo-dT
    priming on a genomic A-run inside a transcript (internal priming), reported but not counted as a 3' end;
  - the library type: stranded or not, from how reads on the protein genes lie relative to the genes.
Reads are mapped to this genome only, mapping quality >= MIN_MAPQ (nuclear copies of mitochondrial DNA map with low
quality to both). Mapper (--mapper): bwa mem (default; best placement, poly(A) tails and depth when the RNA-seq
individual differs from the assembled one), minimap2 -ax sr (the default before 2026-11-10; loses reads at a few %
divergence), bowtie2 --very-sensitive-local. Measured: benchmark/results/rna_mappers_2026-11-10/README.md.

Writes rna_features.tsv (one row per gene / ORF / non-coding region: mean depth, depth relative to the median of
the 13 protein genes, the share of its positions transcribed at >= 10 % of that median - a partly transcribed ORF
shows here -, "transcribed like the genes" when its depth is at least half that of the least-expressed protein gene,
sense share for stranded libraries, poly(A) reads at its 3' end) and rna_polya_sites.tsv (sites with >=
MIN_POLYA_READS reads, the feature they end). A reverse-stranded library (dUTP, read 1 antisense) is recognised and
read accordingly. Measured: benchmark/results/rnaseq_2026-11-10/README.md.

Usage: python3 mito_rna.py -g annotated.gb -1 R1.fq.gz [-2 R2.fq.gz] -o outdir [--threads 8]
       python3 mito_rna.py -g annotated.gb --sam aligned.sam -o outdir          (reads already mapped to it)
"""
from __future__ import annotations

import argparse
import re
import statistics
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, TextIO, Tuple

import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from mito_common import cli_guard, run_tool, write_tsv
from mito_names import PCG, gene_name

MIN_MAPQ = 20               # reads of nuclear copies (NUMTs) map with low quality
POLYA_MIN = 6               # A's at the start of a 3' clip that make a poly(A) tail
MIN_POLYA_READS = 3         # reads for a reported poly(A) site
SITE_WINDOW = 3             # poly(A) ends within this many nt are one site
END_WINDOW = 10             # a site this close to a feature's 3' end belongs to it
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
COVERED = 0.1              # a position is transcribed at >= 10 % of the protein genes' median depth
MIN_MEDIAN_DEPTH = 10      # below this protein-gene median depth the expression calls are left out (counting noise)
INTERNAL_WINDOW = 10        # internal priming: the last 10 nt before a poly(A) site (transcript sense) ...
INTERNAL_RUN = 7            # ... hold an A-run of >= 7 (true 3' ends of 67 annelid genomes: 1 of 1,047) ...
INTERNAL_A = 9              # ... or >= 9 A's (none of the 1,047)
MAPPERS = ("bwa", "minimap2", "bowtie2")
POLYA_FEW = 1e-4            # poly(A) reads below this share of the reads: no 3'-end evidence in this library
FEATURE_COLUMNS = ["record", "feature", "type", "start", "end", "strand", "length", "mean_depth", "relative_depth",
                   "covered_share", "transcribed_like_genes", "sense_share", "polya_reads_at_3prime_end"]
SITE_COLUMNS = ["record", "position", "strand", "reads", "feature_ending_here", "distance_to_feature_end",
                "internal_priming"]


def find_tool(name: str, bin_dir: Optional[str] = None) -> Optional[str]:
    import shutil
    for d in (bin_dir, str(Path(sys.executable).parent)):
        if d and (Path(d) / name).exists():
            return str(Path(d) / name)
    return shutil.which(name)


def map_reads(ref_fa: Path, r1: str, r2: Optional[str], threads: int, bin_dir: Optional[str],
              mapper: str = "bwa") -> subprocess.Popen:
    """Reads -> the mitogenome with `mapper` (bwa mem -Y | minimap2 -ax sr -Y | bowtie2 --very-sensitive-local);
    SAM on stdout. A missing mapper falls back to the next of MAPPERS, with a warning."""
    order = [mapper] + [m for m in MAPPERS if m != mapper]
    for m in order:
        exe = find_tool({"bwa": "bwa", "minimap2": "minimap2", "bowtie2": "bowtie2"}[m], bin_dir)
        if not exe:
            continue
        if m != mapper:
            print(f"[mito_rna] WARNING {mapper} not found: mapping with {m}", file=sys.stderr)
        reads = [r1] + ([r2] if r2 else [])
        if m == "bwa":
            run_tool([exe, "index", str(ref_fa)], check=True)
            cmd = [exe, "mem", "-Y", "-t", str(threads), str(ref_fa)] + reads
        elif m == "minimap2":
            cmd = [exe, "-ax", "sr", "-Y", "-t", str(threads), str(ref_fa)] + reads
        else:
            build = find_tool("bowtie2-build", bin_dir) or "bowtie2-build"
            run_tool([build, "-q", str(ref_fa), str(ref_fa)], check=True)
            cmd = [exe, "--very-sensitive-local", "-p", str(threads), "-x", str(ref_fa)] + \
                (["-1", r1, "-2", r2] if r2 else ["-U", r1])
        return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    raise RuntimeError("no mapper found (bwa, minimap2 or bowtie2; conda install -c bioconda bwa)")


def internal_priming(seq: str, pos: int, strand: int) -> bool:
    """A poly(A) site at 0-based `pos` whose last INTERNAL_WINDOW nt (transcript sense) are A-rich: oligo-dT primed
    on a genomic A-run inside a transcript, not a transcript end."""
    n = len(seq)
    if strand == 1:
        w = "".join(seq[(pos - INTERNAL_WINDOW + 1 + k) % n] for k in range(INTERNAL_WINDOW)).upper()
    else:
        w = "".join(seq[(pos + INTERNAL_WINDOW - 1 - k) % n] for k in range(INTERNAL_WINDOW)).upper()
        w = w.translate(str.maketrans("ACGT", "TGCA"))
    longest = max((len(x) for x in re.findall("A+", w)), default=0)
    return longest >= INTERNAL_RUN or w.count("A") >= INTERNAL_A


def parse_sam(lines: Iterable[str], n: int) -> dict:
    """Depth by fragment strand, poly(A) 3' ends, sense counts per position. The fragment strand of a pair is read 1's
    strand (read 2 reversed) - for a stranded library it is (or is opposite to) the transcript's strand."""
    depth = {1: np.zeros(n, dtype=np.int64), -1: np.zeros(n, dtype=np.int64)}
    polya: Counter = Counter()
    used = low = 0
    for ln in lines:
        if ln.startswith("@"):
            continue
        c = ln.split("\t", 11)
        flag = int(c[1])
        if flag & 0x904:                                    # unmapped, secondary, supplementary
            continue
        if int(c[4]) < MIN_MAPQ:
            low += 1
            continue
        used += 1
        pos = int(c[3]) - 1
        rev = bool(flag & 0x10)
        frag_strand = (-1 if rev else 1) * (-1 if flag & 0x80 else 1)
        ops = CIGAR.findall(c[5])
        ref = pos
        for ln_, op in ops:
            k = int(ln_)
            if op in "M=XD":
                a, b = ref, ref + k
                if b <= n:
                    depth[frag_strand][a:b] += 1
                else:                                        # (circular: past the end)
                    depth[frag_strand][a % n:] += 1
                    depth[frag_strand][: (b - n)] += 1
                ref += k
            elif op == "N":
                ref += k
        seq = c[9]
        # poly(A): SAM holds the read as on the reference's plus strand, whatever the read's own direction, so a
        # tail shows as a clipped A-run after the alignment (a plus-strand transcript's 3' end) or a clipped T-run
        # before it (a minus-strand transcript's); the library's strandedness does not matter here.
        if ops and ops[-1][1] == "S" and seq[len(seq) - int(ops[-1][0]):].startswith("A" * POLYA_MIN):
            polya[((ref - 1) % n, 1)] += 1
        if ops and ops[0][1] == "S" and seq[: int(ops[0][0])].endswith("T" * POLYA_MIN):
            polya[(pos % n, -1)] += 1
    return {"depth": depth, "polya": polya, "used": used, "low_mapq": low}


def features_of(rec: SeqRecord) -> List[Tuple[str, str, int, int, int]]:
    """(name, type, start, end, strand) of genes, rRNAs, tRNAs, ORFs and non-coding regions."""
    out = []
    for f in rec.features:
        if f.type not in ("CDS", "rRNA", "tRNA", "misc_feature") or f.location is None:
            continue
        name = f.qualifiers.get("gene", f.qualifiers.get("product", [""]))[0]
        if f.type == "misc_feature":
            note = " ".join(f.qualifiers.get("note", []))
            if not note.startswith("non-coding region"):
                continue
            name = note.split(",")[0].replace("non-coding region ", "")
        parts = list(f.location.parts)
        out.append((name, f.type, min(int(p.start) for p in parts), max(int(p.end) for p in parts),
                    f.location.strand or 1))
    return out


def polya_sites(polya: Counter) -> List[Tuple[int, int, int]]:
    """(position, strand, reads): poly(A) ends within SITE_WINDOW joined at their most used position."""
    out = []
    for strand in (1, -1):
        pos = sorted((p, c) for (p, s), c in polya.items() if s == strand)
        group: List[Tuple[int, int]] = []
        for p, c in pos + [(10 ** 12, 0)]:
            if group and p - group[-1][0] > SITE_WINDOW:
                top = max(group, key=lambda x: x[1])[0]
                out.append((top, strand, sum(x[1] for x in group)))
                group = []
            group.append((p, c))
    return [x for x in out if x[2] >= MIN_POLYA_READS]


def analyse(rec: SeqRecord, sam: Iterable[str]) -> Tuple[List[dict], List[dict], dict]:
    n = len(rec.seq)
    res = parse_sam(sam, n)
    tot = res["depth"][1] + res["depth"][-1]
    feats = features_of(rec)
    pcg = [f for f in feats if f[1] == "CDS" and gene_name({"gene": [f[0]]}, "CDS") in PCG]
    gene_depths = [float(tot[s:e].mean()) for _, _, s, e, _ in pcg]
    med = statistics.median(gene_depths) if pcg else 0.0
    lowest = min(gene_depths) if pcg else 0.0
    # library: share of depth on the protein genes lying on the gene's strand (0.5 unstranded; near 0: reverse
    # stranded, read 1 antisense as in dUTP libraries - the strands are swapped)
    same = sum(int(res["depth"][st][s:e].sum()) for _, _, s, e, st in pcg)
    allp = sum(int(tot[s:e].sum()) for _, _, s, e, _ in pcg) or 1
    sense = same / allp
    stranded = abs(sense - 0.5) > 0.3
    reverse = stranded and sense < 0.5
    if reverse:
        res["depth"] = {1: res["depth"][-1], -1: res["depth"][1]}
    seq = str(rec.seq)
    sites = [x for x in polya_sites(res["polya"]) if not internal_priming(seq, x[0], x[1])]
    primed = [x for x in polya_sites(res["polya"]) if internal_priming(seq, x[0], x[1])]
    tRNA_ends = {(e - 1 if st == 1 else s, st) for nm, t, s, e, st in feats if t == "tRNA"}
    rows = []
    for nm, t, s, e, st in feats:
        d = float(tot[s:e].mean()) if e > s else 0.0
        end = e - 1 if st == 1 else s
        pa = sum(c for p, sst, c in sites if (sst == st or not stranded) and abs(p - end) <= END_WINDOW)
        ss = float(res["depth"][st][s:e].sum()) / max(1.0, float(tot[s:e].sum()))
        cov = float((tot[s:e] >= COVERED * med).mean()) if e > s and med else 0.0
        rows.append({"record": rec.id, "feature": nm, "type": t, "start": s + 1, "end": e,
                     "strand": "+" if st == 1 else "-", "length": e - s, "mean_depth": f"{d:.1f}",
                     "relative_depth": f"{d / med:.2f}" if med else "", "covered_share": f"{cov:.2f}",
                     "transcribed_like_genes": ("" if med < MIN_MEDIAN_DEPTH else
                                                "yes" if pcg and d >= 0.5 * lowest else "no"),
                     "sense_share": f"{ss:.2f}" if stranded else "", "polya_reads_at_3prime_end": pa})
    srows = []
    for p, st, c in sites + primed:
        if (p, st) in tRNA_ends or any(abs(p - x) <= 1 and sx == st for x, sx in tRNA_ends):
            continue                                     # (CCA-tailed tRNA 3' end)
        ends = [(abs(p - (e - 1 if fs == 1 else s)), nm) for nm, t, s, e, fs in feats
                if t != "misc_feature" and (fs == st or not stranded)]
        dist, nm = min(ends) if ends else (None, "")
        srows.append({"record": rec.id, "position": p + 1, "strand": "+" if st == 1 else "-", "reads": c,
                      "feature_ending_here": nm if dist is not None and dist <= END_WINDOW else "",
                      "distance_to_feature_end": dist if dist is not None else "",
                      "internal_priming": "yes" if (p, st, c) in primed else ""})
    info = {"reads_used": res["used"], "reads_low_mapq": res["low_mapq"], "pcg_median_depth": med,
            "sense_share_on_genes": sense, "stranded": stranded, "reverse": reverse,
            "polya_sites": sum(1 for r in srows if not r["internal_priming"]),
            "internal_priming_sites": sum(1 for r in srows if r["internal_priming"]),
            "polya_reads": sum(res["polya"].values())}
    return rows, srows, info


def main(argv: Sequence[str] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-g", "--genbank", required=True, help="annotated mitogenome (one record)")
    p.add_argument("-1", "--reads1", default=None)
    p.add_argument("-2", "--reads2", default=None)
    p.add_argument("--sam", default=None, help="SAM of reads already mapped to this record ('-' = stdin)")
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--bin", default=None)
    p.add_argument("--mapper", choices=MAPPERS, default="bwa",
                   help="bwa mem (default), minimap2 -ax sr (the default before 2026-11-10), bowtie2 local")
    a = p.parse_args(argv)
    rec = SeqIO.read(a.genbank, "genbank")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.sam:
        fh: TextIO = sys.stdin if a.sam == "-" else open(a.sam)
        rows, sites, info = analyse(rec, fh)
    else:
        if not a.reads1:
            raise SystemExit("give reads (-1 [-2]) or --sam")
        ref = out / "ref.fa"
        ref.write_text(f">{rec.id}\n{rec.seq}\n")
        proc = map_reads(ref, a.reads1, a.reads2, a.threads, a.bin, a.mapper)
        rows, sites, info = analyse(rec, proc.stdout)
        proc.wait()
    write_tsv(out / "rna_features.tsv", FEATURE_COLUMNS, rows)
    write_tsv(out / "rna_polya_sites.tsv", SITE_COLUMNS, sites)
    lib = ("reverse-stranded (read 1 antisense)" if info["reverse"] else "stranded") if info["stranded"] \
        else "unstranded"
    print(f"[mito_rna] {info['reads_used']} reads (MAPQ >= {MIN_MAPQ}; {info['reads_low_mapq']} below), {lib} "
          f"(sense share {info['sense_share_on_genes']:.2f}), protein-gene median depth "
          f"{info['pcg_median_depth']:.1f}, {info['polya_reads']} poly(A) reads, {info['polya_sites']} sites (+ {info['internal_priming_sites']} internal "
          f"priming) -> {out}",
          file=sys.stderr)
    if info["pcg_median_depth"] < MIN_MEDIAN_DEPTH:
        print(f"[mito_rna] WARNING protein-gene median depth {info['pcg_median_depth']:.1f} < {MIN_MEDIAN_DEPTH}: too "
              "few reads for expression calls (transcribed_like_genes left empty)", file=sys.stderr)
    if info["polya_reads"] < POLYA_FEW * max(1, info["reads_used"]):
        print(f"[mito_rna] WARNING {info['polya_reads']} poly(A) reads: no 3'-end evidence in this library (random "
              "primers / total RNA, short reads, or tails trimmed, e.g. fastp --trim_poly_x)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(cli_guard(main))
