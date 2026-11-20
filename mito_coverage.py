#!/usr/bin/env python3
"""mito_coverage.py - read depth of a mitogenome (or any small reference): compute it, summarise it, draw it.

  mitom coverage -i ref.gb|ref.fa --fq1 R1.fq.gz [--fq2 R2.fq.gz] -o OUT [--plot flat,circos,mitom | none]
  mitom coverage -i ref.gb --depth_file depth.tsv -o OUT --plot flat,circos        (draw from `samtools depth` output)

Computing (benchmark/results/coverage_2026-11-19/README.md): the reads are mapped with --mapper (bwa mem, bwa-mem2,
minimap2 -x sr, bowtie2 --sensitive --dovetail; all SAM on a pipe) and the reference span of every primary alignment is
added to a difference array - no BAM, no sort, no `samtools depth`, nothing on disk but the index. A circular
reference is mapped as sequence + its first ORIGIN_EXT nt and the extension is folded back, so reads across the origin
count (GenBank topology, `topology=circular` in the FASTA header, or --circular). Depth is per base, zeros included.

Output (prefix = --prefix, default the reference file name):
  <prefix>.depth.tsv       record, position (1-based), depth - the format of `samtools depth -a`
  <prefix>.coverage.tsv    one row per record: length, mean / median / min / max depth, % bases >= 1x / 10x / 0.2 x mean /
                           100x, mapped_reads (with reads), topology, gc_pct; --depth_gz writes <prefix>.depth.tsv.gz
  <prefix>.<record>.flat.png|pdf     flat map: depth along the sequence on a log or linear axis, genes below (GenBank),
                                     bases with no reads marked (the view of GetOrganelle and of using_bwa_draw_coverage)
  <record>.png|pdf                   circular maps (--plot circos | mitom) with the depth ring, drawn by mito_plot.py
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

import mito_env
import mito_plot
from mito_common import cli_guard, command_line, iter_gb_files, prepare_input, STDERR_MARK, STDERR_TAIL, TOOL_LOG

MAPPERS = ("bwa", "bwa-mem2", "minimap2", "bowtie2")
PLOTS = ("flat", "circos", "mitom")
ORIGIN_EXT = mito_plot.ORIGIN_EXT
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
SKIP_FLAGS = 0x904                       # unmapped, secondary, supplementary
FILE_SUFFIXES = (".gb", ".gbk", ".gbf", ".gbff", ".genbank", ".fa", ".fas", ".fasta", ".fna")


# ---------------------------------------------------------------- reference

def read_reference(path: str) -> List[SeqRecord]:
    """Records of a GenBank or FASTA file / folder of GenBank files; circular = topology annotation, else a
    `topology=circular` word in the FASTA description."""
    src = Path(prepare_input(path))
    if src.is_dir():
        files = iter_gb_files(src)
    else:
        files = [src]
    recs: List[SeqRecord] = []
    for f in files:
        with open(f) as h:
            first = next((ln for ln in h if ln.strip()), "")
        fmt = "fasta" if first.startswith(">") else "genbank"
        for r in SeqIO.parse(str(f), fmt):
            if fmt == "fasta" and "topology=circular" in r.description.replace(" ", ""):
                r.annotations["topology"] = "circular"
            recs.append(r)
    if not recs:
        raise ValueError(f"no records in {path}")
    ids = [r.id for r in recs]
    if len(set(ids)) != len(ids):
        raise ValueError("record IDs are not unique: " + ", ".join(sorted({i for i in ids if ids.count(i) > 1})))
    return recs


def is_circular(rec: SeqRecord, force: Optional[bool]) -> bool:
    return force if force is not None else rec.annotations.get("topology") == "circular"


# ---------------------------------------------------------------- mapping and counting

def _find(name: str) -> str:
    exe = mito_env.find(name)
    if not exe:
        raise RuntimeError(f"{name} not found (put it on PATH or run `mitom doctor`)")
    return exe


def _tool(cmd: Sequence[str], log) -> None:
    rc = subprocess.run(list(cmd), stdout=log, stderr=log).returncode
    if rc:
        raise RuntimeError(f"{command_line(list(cmd))} failed (exit {rc})")


def mapper_command(mapper: str, fa: Path, reads: Sequence[str], threads: int, log) -> List[str]:
    """Index `fa` and return the command that writes SAM to stdout (every read has a line or is omitted)."""
    t = str(threads)
    if mapper == "bwa":
        exe = _find("bwa")
        _tool([exe, "index", str(fa)], log)
        return [exe, "mem", "-t", t, str(fa), *reads]
    if mapper == "bwa-mem2":
        exe = _find("bwa-mem2")
        _tool([exe, "index", str(fa)], log)
        return [exe, "mem", "-t", t, str(fa), *reads]
    if mapper == "minimap2":
        return [_find("minimap2"), "-ax", "sr", "--secondary=no", "-t", t, str(fa), *reads]
    if mapper == "bowtie2":
        exe = _find("bowtie2")
        build = mito_env.find("bowtie2-build") or str(Path(exe).with_name("bowtie2-build"))
        _tool([build, "-q", str(fa), str(fa)], log)
        mates = ["-1", reads[0], "-2", reads[1]] if len(reads) == 2 else ["-U", reads[0]]
        return [exe, "--sensitive", "--no-discordant", "--dovetail", "--no-unal", "-p", t, "-x", str(fa), *mates]
    raise ValueError(f"unknown mapper {mapper!r} ({', '.join(MAPPERS)})")


def count_sam(lines, ref_len: Dict[str, int], min_mapq: int = 0, reads: Optional[Dict[str, int]] = None
              ) -> Dict[str, np.ndarray]:
    """Per-base depth from SAM lines: the reference span (CIGAR M D N = X) of every mapped primary alignment; reads
    (given: filled) = the number of those alignments per record (a pair counts twice)."""
    diff = {k: np.zeros(n + 1, dtype=np.int64) for k, n in ref_len.items()}
    for line in lines:
        if line[0] == "@":
            continue
        f = line.split("\t", 7)
        if int(f[1]) & SKIP_FLAGS or int(f[4]) < min_mapq:
            continue
        d = diff.get(f[2])
        if d is None:
            continue
        start = int(f[3]) - 1
        cig = f[5]
        span = int(cig[:-1]) if cig[-1] == "M" and cig[:-1].isdigit() else \
            sum(int(k) for k, op in CIGAR.findall(cig) if op in "MDN=X")
        d[start] += 1
        d[min(start + span, len(d) - 1)] -= 1
        if reads is not None:
            reads[f[2]] = reads.get(f[2], 0) + 1
    return {k: np.cumsum(d[:-1]) for k, d in diff.items()}


def depth_from_reads(recs: Sequence[SeqRecord], reads: Sequence[str], mapper: str = "bwa", threads: int = 1,
                     circular: Optional[bool] = None, min_mapq: int = 0, keep_log: Optional[Path] = None,
                     reads_out: Optional[Dict[str, int]] = None) -> Dict[str, np.ndarray]:
    """{record id: per-base depth} of `recs` from one FASTQ or a pair. A circular record is mapped with its first
    ORIGIN_EXT nt appended and the extension is folded back onto the start. reads_out (given: filled): mapped
    alignments per record."""
    if not 1 <= len(reads) <= 2:
        raise ValueError(f"1 or 2 FASTQ files expected, got {len(reads)}")
    ext: Dict[str, int] = {}
    with tempfile.TemporaryDirectory(prefix="mito_cov_") as tmp:
        fa, logf = Path(tmp) / "ref.fa", Path(tmp) / "mapper.log"
        with open(fa, "w") as h:
            for r in recs:
                s = str(r.seq)
                ext[r.id] = min(ORIGIN_EXT, len(s) - 1) if is_circular(r, circular) else 0
                h.write(f">{r.id}\n{s + s[:ext[r.id]]}\n")
        ref_len = {r.id: len(r.seq) + ext[r.id] for r in recs}
        t0 = time.time()
        with open(logf, "w") as log:
            cmd = mapper_command(mapper, fa, list(map(str, reads)), threads, log)
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log, text=True)
            cov = count_sam(p.stdout, ref_len, min_mapq, reads_out)
            rc = p.wait()
        tail = logf.read_text().strip().splitlines()[-STDERR_TAIL:]
        TOOL_LOG.debug(f"[CMD rc={rc} {time.time() - t0:.1f}s] {command_line(cmd)}" +
                       "".join(f"\n{STDERR_MARK}{ln}" for ln in tail))
        if rc:
            raise RuntimeError(f"{mapper} failed (exit {rc}): " + "\n".join(tail))
        if keep_log:
            keep_log.write_text(logf.read_text())
    out = {}
    for r in recs:
        n, e = len(r.seq), ext[r.id]
        d = cov[r.id][:n].copy()
        if e:
            d[:e] += cov[r.id][n:]
        out[r.id] = d
    return out


def read_depth_file(path: str, recs: Sequence[SeqRecord]) -> Dict[str, np.ndarray]:
    """`samtools depth` (record, position, depth; positions without a line = 0) -> {record id: depth}."""
    out = {r.id: np.zeros(len(r.seq), dtype=np.int64) for r in recs}
    with open(path) as h:
        for line in h:
            f = line.split("\t")
            if len(f) < 3 or f[0] not in out:
                continue
            pos = int(f[1])
            if 1 <= pos <= len(out[f[0]]):
                out[f[0]][pos - 1] = int(f[2])
    return out


# ---------------------------------------------------------------- summary and tables

def summarise(depth: np.ndarray) -> Dict[str, float]:
    mean = float(depth.mean()) if len(depth) else 0.0
    cut = 0.2 * mean
    return {"length": len(depth), "mean": mean, "median": float(np.median(depth)), "min": int(depth.min()),
            "max": int(depth.max()), "pct_ge1": 100 * float((depth >= 1).mean()),
            "pct_ge10": 100 * float((depth >= 10).mean()),
            "pct_ge100": 100 * float((depth >= 100).mean()), "pct_ge_0.2mean": 100 * float((depth >= cut).mean())}


def gc_pct(rec: SeqRecord) -> float:
    seq = str(rec.seq).upper()
    acgt = sum(seq.count(c) for c in "ACGT")
    return 100 * (seq.count("G") + seq.count("C")) / acgt if acgt else 0.0


def write_tables(depths: Dict[str, np.ndarray], prefix: Path, recs: Optional[Sequence[SeqRecord]] = None,
                 mapped: Optional[Dict[str, int]] = None, gz: bool = False) -> List[Path]:
    """<prefix>.depth.tsv (.gz: gz) and <prefix>.coverage.tsv. recs adds topology and GC %, mapped (reads per record)
    adds mapped_reads (issue #12); columns after the first eight keep their order for old readers."""
    import gzip
    dpath = Path(f"{prefix}.depth.tsv" + (".gz" if gz else ""))
    spath = Path(f"{prefix}.coverage.tsv")
    with (gzip.open(dpath, "wt") if gz else open(dpath, "w")) as h:
        for rid, d in depths.items():
            h.write("".join(f"{rid}\t{i}\t{v}\n" for i, v in enumerate(d.tolist(), 1)))
    cols = ["length", "mean", "median", "min", "max", "pct_ge1", "pct_ge10", "pct_ge_0.2mean", "pct_ge100"]
    by_id = {r.id: r for r in recs or []}
    extra = (["mapped_reads"] if mapped is not None else []) + (["topology", "gc_pct"] if by_id else [])
    with open(spath, "w") as h:
        h.write("record\t" + "\t".join(cols + extra) + "\n")
        for rid, d in depths.items():
            s = summarise(d)
            cells = [f"{s[c]:.2f}" if isinstance(s[c], float) else str(s[c]) for c in cols]
            if mapped is not None:
                cells.append(str(mapped.get(rid, 0)))
            if by_id:
                r = by_id.get(rid)
                cells += [("circular" if is_circular(r, None) else "linear") if r is not None else "",
                          f"{gc_pct(r):.2f}" if r is not None else ""]
            h.write(rid + "\t" + "\t".join(cells) + "\n")
    return [dpath, spath]


# ---------------------------------------------------------------- flat map

def _runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    """(start, end) of the runs of True in `mask` (end exclusive)."""
    if not mask.any():
        return []
    m = np.concatenate(([0], mask.astype(np.int8), [0]))
    ch = np.flatnonzero(np.diff(m))
    return list(zip(ch[::2].tolist(), ch[1::2].tolist()))


def flat_figure(rec: SeqRecord, depth: np.ndarray, log: bool = False, window: int = 0) -> Figure:
    """Depth along the sequence (x: position, y: depth; linear, log with log=True) with the gene track below the curve and
    uncovered bases marked; `window` > 1 smooths the curve with a running mean of that many bases."""
    n = len(depth)
    s = summarise(depth)
    fig = Figure(figsize=(12, 4.8 if any(mito_plot.category(f) for f in rec.features) else 3.6))
    genes = [(f, mito_plot.category(f)) for f in rec.features if mito_plot.category(f)]
    if genes:
        ax, gx = fig.subplots(2, 1, sharex=True, height_ratios=[5, 1], gridspec_kw={"hspace": 0.05})
    else:
        ax, gx = fig.subplots(1, 1), None
    y = depth.astype(float)
    if window > 1:
        y = np.convolve(y, np.ones(window) / window, mode="same")
    x = np.arange(1, n + 1)
    ytop = max(s["max"], 1) * 1.15
    if log:
        ax.set_yscale("symlog", linthresh=1)
        ytop = max(s["max"], 1) * 2
    ax.fill_between(x, 0, y, step="mid", color=mito_plot.DEPTH, alpha=0.55, linewidth=0)
    ax.plot(x, y, color=mito_plot.INK, linewidth=0.6, drawstyle="steps-mid")
    ax.axhline(s["mean"], color=mito_plot.SKEW_POS, linestyle="--", linewidth=0.9)
    ax.text(n * 1.003, s["mean"], f"mean {s['mean']:,.0f}x", va="center", fontsize=8, color=mito_plot.SKEW_POS,
            clip_on=False)
    for a, b in _runs(depth == 0):                       # bases with no read
        ax.axvspan(a + 0.5, b + 0.5, color="#d03030", alpha=0.35, linewidth=0)
    ax.set_xlim(0.5, n + 0.5)
    ax.set_ylim(0, ytop)
    ax.set_ylabel("depth (x, log)" if log else "depth (x)")
    ax.grid(True, axis="y", linestyle=":", alpha=0.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_title(f"{rec.id}  {n:,} bp   mean {s['mean']:,.1f}x   median {s['median']:,.0f}x   "
                 f"min {s['min']}x   max {s['max']:,}x   >=1x {s['pct_ge1']:.1f}%   "
                 f"{'circular' if rec.annotations.get('topology') == 'circular' else 'linear'}",
                 fontsize=10, loc="left")
    if gx is not None:
        gx.set_ylim(-1.2, 1.2)
        gx.set_yticks([])
        for sp in ("top", "right", "left"):
            gx.spines[sp].set_visible(False)
        for f, cat in genes:
            strand = f.location.strand or 1
            for part in f.location.parts:
                a, b = int(part.start), int(part.end)
                gx.add_patch(Rectangle((a + 0.5, 0.05 if strand > 0 else -1.0), b - a, 0.95,
                                       color=mito_plot.COLORS[cat], linewidth=0))
            lab = mito_plot.label(f, cat)
            if lab and cat != "tRNA":
                a = int(f.location.start)
                gx.text((a + int(f.location.end)) / 2, 0.5 if strand > 0 else -0.5, lab, ha="center", va="center",
                        fontsize=6.5, color="white" if cat in ("COX", "ATP", "rRNA", "NADH") else mito_plot.INK)
        gx.set_xlabel("position (bp)")
    else:
        ax.set_xlabel("position (bp)")
    return fig


def save_figure(fig: Figure, stem: Path, formats: Sequence[str]) -> List[Path]:
    out = []
    with matplotlib.rc_context({"pdf.fonttype": 42}):
        for fmt in formats:
            p = Path(f"{stem}.{fmt}")
            fig.savefig(p, dpi=300, facecolor="white", bbox_inches="tight")
            out.append(p)
    return out


# ---------------------------------------------------------------- command line

def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="reference: GenBank / FASTA file, or a folder of GenBank files")
    p.add_argument("-o", "--out", required=True, metavar="OUT", help="output folder")
    p.add_argument("--fq1", help="reads (FASTQ, .gz allowed)")
    p.add_argument("--fq2", help="mates of --fq1")
    p.add_argument("--depth_file", help="an existing `samtools depth` table (record, position, depth) instead of reads")
    p.add_argument("--mapper", choices=MAPPERS, default="bwa", help="read mapper (default bwa)")
    p.add_argument("--threads", type=int, default=4, help="mapper threads (default 4)")
    p.add_argument("--min_mapq", type=int, default=0, help="count only alignments with this MAPQ or more (default 0)")
    p.add_argument("--circular", action="store_true", default=None,
                   help="treat every record as circular (default: GenBank topology / FASTA 'topology=circular')")
    p.add_argument("--plot", default="flat",
                   help="comma-separated: flat (depth along the sequence), circos / mitom (circular maps, GenBank "
                        "input), or none = tables only (default flat)")
    p.add_argument("--log_scale", action="store_true", help="flat map with a log depth axis (default linear)")
    p.add_argument("--smooth", type=int, default=0, metavar="BP", help="flat map: running mean over BP bases")
    p.add_argument("--formats", default="png,pdf", help="comma-separated (default png,pdf)")
    p.add_argument("--prefix", help="prefix of the output files (default: the reference file name)")
    p.add_argument("--depth_gz", action="store_true", help="write <prefix>.depth.tsv.gz (default: plain .tsv)")
    a = p.parse_args(argv)
    if bool(a.fq1) == bool(a.depth_file):
        p.error("give reads (--fq1 [--fq2]) or --depth_file")
    if a.fq2 and not a.fq1:
        p.error("--fq2 needs --fq1")
    plots = [] if a.plot.strip() in ("", "none") else [x.strip() for x in a.plot.split(",")]
    bad = [x for x in plots if x not in PLOTS]
    if bad:
        p.error(f"unknown --plot {bad} (choose from {', '.join(PLOTS)} or none)")
    formats = [x.strip() for x in a.formats.split(",") if x.strip()]
    mapped: Optional[Dict[str, int]] = None if a.depth_file else {}
    try:
        recs = read_reference(a.input)
        if a.depth_file:
            depths = read_depth_file(a.depth_file, recs)
        else:
            depths = depth_from_reads(recs, [x for x in (a.fq1, a.fq2) if x], a.mapper, a.threads, a.circular,
                                      a.min_mapq, reads_out=mapped)
    except (RuntimeError, ValueError, OSError) as exc:
        print(f"[mito_coverage] {exc}", file=sys.stderr)
        return 1
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    prefix = out / (a.prefix or re.sub(r"\.(gb\w*|genbank|fa\w*|fna)$", "", Path(a.input).name))
    written = write_tables(depths, prefix, recs, mapped, a.depth_gz)
    used: set = set()
    for rec in recs:
        d = depths[rec.id]
        if "flat" in plots:
            safe = re.sub(r"[^\w.-]", "_", rec.id)
            stem = Path(f"{prefix}.flat" if len(recs) == 1 else f"{prefix}.{safe}.flat")
            written += save_figure(flat_figure(rec, d, a.log_scale, a.smooth), stem, formats)
        for style in (s for s in plots if s != "flat"):
            if not any(mito_plot.category(f) for f in rec.features):
                print(f"[mito_coverage] {rec.id}: no gene annotation, no {style} map (use --plot flat)",
                      file=sys.stderr)
                continue
            sub = out / style if len([s for s in plots if s != "flat"]) > 1 else out
            written += mito_plot.plot_record(rec, sub, d, formats, style, used=used if sub == out else set())
    for w in written:
        print(w)
    return 0


if __name__ == "__main__":
    sys.exit(cli_guard(main))
