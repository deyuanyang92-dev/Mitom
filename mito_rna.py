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
sense share for stranded libraries, poly(A) reads at its 3' end), rna_depth.tsv.gz (depth per position and strand)
and rna_polya_sites.tsv (sites with >=
MIN_POLYA_READS reads, the feature they end). A reverse-stranded library (dUTP, read 1 antisense) is recognised and
read accordingly. Measured: benchmark/results/rnaseq_2026-11-10/README.md.

Refinement (--refine, on a finished output folder; no mapping): rna_refined.gb and rna_end_corrections.tsv. A protein
gene's 3' end is CONFIRMED by >= MIN_POLYA_READS poly(A) reads within END_WINDOW of it (note added). It is MOVED to a
poly(A) site inside the gene's last CORRECT_MAX_NT nt when the annotated end is not confirmed (Travisia COX3: a site
6 codons inside, 11 reads, would have cut the conserved ...YWWGS end), the site has >= CORRECT_MIN_READS reads, is not
internal priming, and leaves an in-frame incomplete stop T / TA (completed by polyadenylation;
/transl_except, translation and a note with the old end). ORFs get their expression as a note.

Usage: python3 mito_rna.py -g annotated.gb -1 R1.fq.gz [-2 R2.fq.gz] -o outdir [--threads 8]
       python3 mito_rna.py -g annotated.gb --sam aligned.sam -o outdir          (reads already mapped to it)
       python3 mito_rna.py -g annotated.gb -o outdir --refine [--run SRR...]    (evidence -> annotation)
"""
from __future__ import annotations

import argparse
import gzip
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

from mito_common import (POLYA_NOTE, abbreviated_stop, cli_guard, codon_translate, run_tool, transl_except_term,
                         write_gb, write_tsv)
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
CORRECT_MIN_READS = 10      # poly(A) reads for a site to move a gene's 3' end
CORRECT_MAX_NT = 90         # ... lying inside the gene's last 90 nt
UTR_MAX = 30                # a gene with a complete stop: a site up to 30 nt beyond (no feature between) confirms it
RRNA_MAX_NT = 30            # rRNA 3' end moved to a site within 30 nt ...
RRNA_SHARE = 0.9            # ... holding >= 90 % of the poly(A) reads there, when the annotated end has none
DUPLICATED = 0.5            # share of a feature's 25-mers found elsewhere in the record: reads map ambiguously
CORRECTION_COLUMNS = ["record", "gene", "strand", "old_start", "old_end", "new_start", "new_end", "old_aa", "new_aa",
                      "incomplete_stop", "site_reads", "reads_at_old_end", "kind"]
DIVERGED = 0.02             # reads differ from the genome by more than this: another individual / population (or
                            # species); fast-evolving genes then get fewer reads - their low depth is not low expression
POLYA_FEW = 1e-4            # poly(A) reads below this share of the reads: no 3'-end evidence in this library
FEATURE_COLUMNS = ["record", "feature", "type", "start", "end", "strand", "length", "mean_depth", "relative_depth",
                   "covered_share", "transcribed_like_genes", "sense_share", "polya_reads_at_3prime_end",
                   "duplicated_share"]
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


def parse_sam(lines: Iterable[str], n: int, ref_name: Optional[str] = None) -> dict:
    """Depth by fragment strand, poly(A) 3' ends, sense counts per position. The fragment strand of a pair is read 1's
    strand (read 2 reversed) - for a stranded library it is (or is opposite to) the transcript's strand."""
    depth = {1: np.zeros(n, dtype=np.int64), -1: np.zeros(n, dtype=np.int64)}
    polya: Counter = Counter()
    used = low = 0
    edits = aligned = 0
    n_sq = 0
    for ln in lines:
        if ln.startswith("@"):
            n_sq += ln.startswith("@SQ")
            continue
        c = ln.split("\t", 11)
        flag = int(c[1])
        if flag & 0x904:                                    # unmapped, secondary, supplementary
            continue
        if ref_name is not None and n_sq != 1 and c[2] != ref_name:   # another reference sequence in the BAM: not this genome
            continue
        if int(c[4]) < MIN_MAPQ:
            low += 1
            continue
        used += 1
        nm = re.search(r"\tNM:i:(\d+)", ln)
        if nm:
            edits += int(nm.group(1))
            aligned += sum(int(k) for k, op in CIGAR.findall(c[5]) if op in "MI=X")
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
    return {"depth": depth, "polya": polya, "used": used, "low_mapq": low,
            "edit_rate": edits / aligned if aligned else 0.0}


def _feature_rows(rec: SeqRecord) -> list:
    """(feature, name, type) of genes, rRNAs, tRNAs, ORFs and non-coding regions, in record order."""
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
        out.append((f, name, f.type))
    return out


def features_of(rec: SeqRecord) -> List[Tuple[str, str, int, int, int]]:
    """(name, type, start, end, strand) of genes, rRNAs, tRNAs, ORFs and non-coding regions (bounding span)."""
    out = []
    for f, name, t in _feature_rows(rec):
        parts = list(f.location.parts)
        out.append((name, t, min(int(p.start) for p in parts), max(int(p.end) for p in parts),
                    f.location.strand or 1))
    return out


def feature_parts(rec: SeqRecord) -> List[Tuple[List[Tuple[int, int]], int]]:
    """For each row of features_of: ([(start, end) of every part, 0-based], 0-based 3' end base in transcript sense).
    A joined feature is measured over its parts, not the span between them (introns, the origin of a circular record)."""
    out = []
    for f, _, _ in _feature_rows(rec):
        parts = [(int(p.start), int(p.end)) for p in f.location.parts]
        st = f.location.strand or 1
        last = list(f.location.parts)[-1]            # parts are in transcript order
        out.append((parts, int(last.end) - 1 if st == 1 else int(last.start)))
    return out


def _over(arr, parts) -> np.ndarray:
    return np.concatenate([arr[a:b] for a, b in parts]) if parts else arr[:0]


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


def _kmer_counts(seq: str, k: int = 25) -> Counter:
    """25-mers of both strands of the (circular) record."""
    s2 = seq + seq[: k - 1]
    rc = s2[::-1].translate(str.maketrans("ACGT", "TGCA"))
    return Counter(x[i: i + k] for x in (s2, rc) for i in range(len(seq)))


def _duplicated_share(seg: str, kcount: Counter, k: int = 25) -> float:
    """Share of a feature's 25-mers that occur more than once in the record (either strand): there reads get mapping
    quality 0 and are not counted, so the depth says nothing."""
    km = [seg[i: i + k] for i in range(len(seg) - k + 1)]
    return sum(1 for x in km if kcount[x] > 1) / len(km) if km else 0.0


def analyse(rec: SeqRecord, sam: Iterable[str], ref: Optional[str] = None) -> Tuple[List[dict], List[dict], dict]:
    n = len(rec.seq)
    res = parse_sam(sam, n, ref)
    tot = res["depth"][1] + res["depth"][-1]
    feats = features_of(rec)
    geom = feature_parts(rec)
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
    kcount = _kmer_counts(str(rec.seq).upper())
    useq = str(rec.seq).upper()
    for (nm, t, s, e, st), (parts, end) in zip(feats, geom):
        tt = _over(tot, parts)
        d = float(tt.mean()) if tt.size else 0.0
        dup = _duplicated_share("".join(useq[a:b] for a, b in parts), kcount)
        pa = sum(c for p, sst, c in sites if (sst == st) and abs(p - end) <= END_WINDOW)
        ss = float(_over(res["depth"][st], parts).sum()) / max(1.0, float(tt.sum()))
        cov = float((tt >= COVERED * med).mean()) if tt.size and med else 0.0
        rows.append({"record": rec.id, "feature": nm, "type": t, "start": s + 1, "end": e,
                     "strand": "+" if st == 1 else "-", "length": e - s, "mean_depth": f"{d:.1f}",
                     "relative_depth": f"{d / med:.2f}" if med else "", "covered_share": f"{cov:.2f}",
                     "transcribed_like_genes": ("" if med < MIN_MEDIAN_DEPTH or dup > DUPLICATED else
                                                "yes" if pcg and d >= 0.5 * lowest else "no"),
                     "sense_share": f"{ss:.2f}" if stranded else "", "polya_reads_at_3prime_end": pa,
                     "duplicated_share": f"{dup:.2f}"})
    srows = []
    for p, st, c in sites + primed:
        if (p, st) in tRNA_ends or any(abs(p - x) <= 1 and sx == st for x, sx in tRNA_ends):
            continue                                     # (CCA-tailed tRNA 3' end)
        ends = [(abs(p - g[1]), nm) for (nm, t, s, e, fs), g in zip(feats, geom)
                if t != "misc_feature" and fs == st]
        dist, nm = min(ends) if ends else (None, "")
        srows.append({"record": rec.id, "position": p + 1, "strand": "+" if st == 1 else "-", "reads": c,
                      "feature_ending_here": nm if dist is not None and dist <= END_WINDOW else "",
                      "distance_to_feature_end": dist if dist is not None else "",
                      "internal_priming": "yes" if (p, st, c) in primed else ""})
    info = {"reads_used": res["used"], "reads_low_mapq": res["low_mapq"], "pcg_median_depth": med,
            "sense_share_on_genes": sense, "stranded": stranded, "reverse": reverse,
            "polya_sites": sum(1 for r in srows if not r["internal_priming"]),
            "internal_priming_sites": sum(1 for r in srows if r["internal_priming"]),
            "polya_reads": sum(res["polya"].values()), "edit_rate": res["edit_rate"],
            "depth": res["depth"]}
    return rows, srows, info


def _site_reads(sites: List[dict], pos1: int, strand: str, stranded: bool, window: int) -> int:
    return sum(int(x["reads"]) for x in sites if not x.get("internal_priming")
               and (x["strand"] == strand) and abs(int(x["position"]) - pos1) <= window)


def _stop_codons(code: int) -> set:
    from Bio.Data import CodonTable
    return set(CodonTable.unambiguous_dna_by_id[code].stop_codons)


def refine(rec: SeqRecord, features: List[dict], sites: List[dict], stranded: bool, run: str = "") -> List[dict]:
    """Write the RNA evidence into rec (in place): confirmed / moved 3' ends of protein genes and rRNAs, ORF
    expression, duplicated genes. -> corrections (see the module docstring for the rules)."""
    from Bio.SeqFeature import FeatureLocation
    src = f" ({run})" if run else ""
    expr = {(r["feature"], int(r["start"])): r for r in features}
    out = []
    n = len(rec.seq)
    genes = {f.qualifiers.get("gene", [""])[0]: f for f in rec.features if f.type == "gene"}
    starts = {(int(f.location.start) if (f.location.strand or 1) == 1 else int(f.location.end) - 1,
               f.location.strand or 1)
              for f in rec.features if f.type in ("CDS", "tRNA", "rRNA") and f.location is not None}
    seq = str(rec.seq).upper()
    kcount = _kmer_counts(seq)
    usable = [x for x in sites if not x.get("internal_priming")]

    def free_between(a: int, b: int, st: int) -> bool:      # no feature starts on strand st in [a, b) (0-based)
        return not any(st2 == st and a <= p < b for p, st2 in starts)

    def move(f, name, loc, nt, reads, at_end, kind):
        tail = len(nt) % 3
        s0, e0 = int(f.location.start), int(f.location.end)
        old = (s0 + 1, e0, (e0 - s0) // 3)
        f.location = loc
        if name in genes:
            genes[name].location = loc
        notes = [x for x in f.qualifiers.get("note", []) if x != "stop codon not determined"]
        if f.type == "CDS":
            code = int(f.qualifiers.get("transl_table", ["5"])[0])
            f.qualifiers["transl_except"] = [transl_except_term(int(loc.start) + 1, int(loc.end), loc.strand, tail)]
            f.qualifiers["translation"] = [codon_translate(nt, code).rstrip("*")]
            notes.append(POLYA_NOTE)
        f.qualifiers["note"] = notes + [f"3' end from {reads} poly(A) RNA-seq reads{src} (was {old[0]}..{old[1]})"]
        out.append({"record": rec.id, "gene": name, "strand": "+" if loc.strand == 1 else "-", "old_start": old[0],
                    "old_end": old[1], "new_start": int(loc.start) + 1, "new_end": int(loc.end),
                    "old_aa": old[2] if f.type == "CDS" else "", "new_aa": len(nt) // 3 if f.type == "CDS" else "",
                    "incomplete_stop": nt[-tail:] if f.type == "CDS" and tail else "", "site_reads": reads,
                    "reads_at_old_end": at_end, "kind": kind})

    for f in rec.features:
        if f.type not in ("CDS", "rRNA") or f.location is None or len(f.location.parts) > 1:
            continue
        name = f.qualifiers.get("gene", f.qualifiers.get("product", [""]))[0]
        st = f.location.strand or 1
        s0, e0 = int(f.location.start), int(f.location.end)
        L = e0 - s0
        sign = "+" if st == 1 else "-"
        end1 = e0 if st == 1 else s0 + 1                                  # 1-based 3'-end position
        at_end = _site_reads(usable, end1, sign, stranded, END_WINDOW)
        notes = f.qualifiers.setdefault("note", [])
        r = expr.get((name, s0 + 1))
        dup = _duplicated_share(seq[s0:e0], kcount)
        if dup > DUPLICATED:
            notes.append(f"RNA-seq{src}: {dup:.0%} of this sequence occurs elsewhere in the "
                         "record, reads there map ambiguously: depth not measurable")
        elif name.lower().startswith("orf") and r and r["relative_depth"]:
            notes.append(f"RNA-seq{src}: depth {r['relative_depth']} x the protein genes' median, "
                         f"{float(r['covered_share']):.0%} of its length transcribed")
        cand = []                                    # (k = transcript length up to the site, reads)
        for x in usable:
            if x["strand"] != sign:      # a poly(A) site's strand is the transcript's, stranded library or not
                continue
            pos = int(x["position"])                                       # 1-based last base before the tail
            cand.append((pos - s0 if st == 1 else e0 - pos + 1, int(x["reads"])))
        done = False
        if f.type == "CDS":
            code = int(f.qualifiers.get("transl_table", ["5"])[0])
            nt0 = str(f.location.extract(rec.seq)).upper()
            has_stop = (len(nt0) % 3 == 0 and nt0[-3:] in _stop_codons(code)) or abbreviated_stop(nt0, code)
            if not has_stop:              # no stop codon: a site beyond the end does not confirm it
                at_end = sum(reads for k, reads in cand if L - END_WINDOW <= k <= L)
        if f.type == "CDS" and at_end < MIN_POLYA_READS:                    # a confirmed end is never moved
            best = None
            for k, reads in cand:
                shorten = max(3, L - CORRECT_MAX_NT) <= k <= L - 1
                extend = not has_stop and L + 1 <= k <= L + CORRECT_MAX_NT
                if not (shorten or extend) or reads < CORRECT_MIN_READS or reads <= at_end:
                    continue
                if st == 1 and s0 + k > n or st == -1 and e0 - k < 0:
                    continue
                loc = FeatureLocation(s0, s0 + k, st) if st == 1 else FeatureLocation(e0 - k, e0, st)
                nt = str(loc.extract(rec.seq)).upper()
                if not abbreviated_stop(nt, code) or "*" in codon_translate(nt, code):
                    continue
                if extend and not free_between(e0 if st == 1 else e0 - k, s0 + k if st == 1 else s0, st):
                    continue
                if best is None or reads > best[1]:
                    best = (loc, reads, nt, "extended" if extend else "shortened")
            if best:
                move(f, name, best[0], best[2], best[1], at_end, best[3])
                done = True
            elif has_stop and len(nt0) % 3 == 0:                          # complete stop: a short 3' UTR
                utr = [(k - L, reads) for k, reads in cand if END_WINDOW < k - L <= UTR_MAX and reads >= MIN_POLYA_READS
                       and free_between(e0 if st == 1 else s0 - (k - L), e0 + (k - L) if st == 1 else s0, st)]
                if utr:
                    d, reads = max(utr, key=lambda t: t[1])
                    notes.append(f"3' end confirmed by {reads} poly(A) RNA-seq reads{src}: the transcript ends {d} nt "
                                 "after the stop codon")
                    done = True
        if f.type == "rRNA" and at_end < MIN_POLYA_READS:
            near = [(k, reads) for k, reads in cand if abs(k - L) <= RRNA_MAX_NT and k > 0]
            tot = sum(reads for _, reads in near)
            if near:
                k, reads = max(near, key=lambda t: t[1])
                ok = reads >= CORRECT_MIN_READS and reads >= RRNA_SHARE * tot and \
                    (k <= L or free_between(e0 if st == 1 else e0 - k, s0 + k if st == 1 else s0, st))
                if ok and 0 < (s0 + k if st == 1 else e0 - k) < n:
                    loc = FeatureLocation(s0, s0 + k, st) if st == 1 else FeatureLocation(e0 - k, e0, st)
                    move(f, name, loc, str(loc.extract(rec.seq)), reads, at_end, "rRNA 3' end")
                    done = True
        if not done and at_end >= MIN_POLYA_READS:
            notes.append(f"3' end confirmed by {at_end} poly(A) RNA-seq reads{src}")
        if not f.qualifiers.get("note"):
            f.qualifiers.pop("note", None)
    return out


def main_refine(a) -> int:
    import csv
    rec = SeqIO.read(a.genbank, "genbank")
    out = Path(a.out)
    rd = lambda n: list(csv.DictReader(open(out / n), delimiter="\t"))
    features, sites = rd("rna_features.tsv"), rd("rna_polya_sites.tsv")
    stranded = any(r["sense_share"] for r in features)
    corr = refine(rec, features, sites, stranded, a.run or "")
    write_gb(rec, out / "rna_refined.gb")
    write_tsv(out / "rna_end_corrections.tsv", CORRECTION_COLUMNS, corr)
    conf = sum(1 for f in rec.features
               if f.type == "CDS" and any("confirmed by" in n for n in f.qualifiers.get("note", [])))
    print(f"[mito_rna] refine {rec.id}: {conf} protein-gene 3' ends confirmed, {len(corr)} moved "
          + "; ".join(f"{c['gene']} {c['old_aa']}->{c['new_aa']} codons" for c in corr) + f" -> {out}/rna_refined.gb",
          file=sys.stderr)
    return 0


def main(argv: Sequence[str] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-g", "--genbank", required=True, help="annotated mitogenome (one record)")
    p.add_argument("-1", "--reads1", default=None)
    p.add_argument("-2", "--reads2", default=None)
    p.add_argument("--sam", default=None, help="SAM of reads already mapped to this record ('-' = stdin)")
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--bin", default=None)
    p.add_argument("--refine", action="store_true",
                   help="write the evidence of a finished output folder (-o) into the annotation: rna_refined.gb")
    p.add_argument("--run", default=None, help="run accession quoted in the notes (--refine)")
    p.add_argument("--mapper", choices=MAPPERS, default="bwa",
                   help="bwa mem (default), minimap2 -ax sr (the default before 2026-11-10), bowtie2 local")
    a = p.parse_args(argv)
    if a.refine:
        return main_refine(a)
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
    with gzip.open(out / "rna_depth.tsv.gz", "wt") as fh:             # per position, transcript strand when stranded
        fh.write("position\tplus\tminus\n")
        d = info["depth"]
        for i in range(len(rec.seq)):
            fh.write(f"{i + 1}\t{d[1][i]}\t{d[-1][i]}\n")
    lib = ("reverse-stranded (read 1 antisense)" if info["reverse"] else "stranded") if info["stranded"] \
        else "unstranded"
    print(f"[mito_rna] {info['reads_used']} reads (MAPQ >= {MIN_MAPQ}; {info['reads_low_mapq']} below), {lib} "
          f"(sense share {info['sense_share_on_genes']:.2f}), protein-gene median depth "
          f"{info['pcg_median_depth']:.1f}, edit rate {info['edit_rate']:.4f}, {info['polya_reads']} poly(A) "
          f"reads, {info['polya_sites']} sites (+ {info['internal_priming_sites']} internal priming) -> {out}",
          file=sys.stderr)
    if info["pcg_median_depth"] < MIN_MEDIAN_DEPTH:
        print(f"[mito_rna] WARNING protein-gene median depth {info['pcg_median_depth']:.1f} < {MIN_MEDIAN_DEPTH}: too "
              "few reads for expression calls (transcribed_like_genes left empty)", file=sys.stderr)
    if info["edit_rate"] > DIVERGED:
        print(f"[mito_rna] WARNING the reads differ from this genome by {100 * info['edit_rate']:.1f} % (another "
              "individual, population or species): fast-evolving genes get fewer reads, so their low relative depth "
              "does not show low expression", file=sys.stderr)
    if info["polya_reads"] < POLYA_FEW * max(1, info["reads_used"]):
        print(f"[mito_rna] WARNING {info['polya_reads']} poly(A) reads: no 3'-end evidence in this library (random "
              "primers / total RNA, short reads, or tails trimmed, e.g. fastp --trim_poly_x)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(cli_guard(main))
