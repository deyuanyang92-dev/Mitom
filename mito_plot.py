#!/usr/bin/env python3
"""mito_plot.py - circular map of an annotated mitochondrial genome (MitoM), drawn with matplotlib.

Replaces the MitoZ / circos map. From the outside in:
  genes       + strand outside the backbone circle, - strand inside, coloured by class (COX, NADH, ATP, CYTB,
              rRNA, tRNA; other = D-loop, misc_feature, rep_origin, repeat_region and CDS of other names, e.g.
              intron ORFs). Protein genes and rRNAs are labelled by /gene, tRNAs by their one-letter code (L1 / L2,
              S1 / S2 kept); crowded labels are fanned out with leader lines. Every part of a join is drawn (a gene
              across the origin, a frameshift join, the exons of a gene with an intron, linked by a line).
  GC content  deviation of each window from the genome mean (WINDOW nt windows every STEP nt)
  GC skew     (G - C) / (G + C) of the same windows
  depth       only with reads: minimap2 -x sr -c when it is installed (about 8x faster), else bwa mem; the
              reference span of every primary alignment; a circular record is mapped as sequence + its first
              ORIGIN_EXT nt and positions past the end are folded back onto the start
  centre      organism, length, topology, accession; a kb ruler runs inside the innermost ring
Position 1 is at the top and positions run clockwise; a linear record leaves a gap between its end and start.

Usage:
  python3 mito_plot.py -i <GenBank file or directory> -o <out_dir> [--fq1 R1.fq.gz [--fq2 R2.fq.gz]]
                       [--threads N] [--formats png,pdf]
Every record of every *.gb / *.gbk / *.gbf file (directory: top level only) becomes <out_dir>/<record id>.png
(300 dpi) and .pdf; reads are used only when the input holds one record.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import matplotlib
import numpy as np
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.textpath import TextPath

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402

WINDOW, STEP, ORIGIN_EXT = 100, 20, 1000
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
# five hues that stay apart pairwise (any two classes can be neighbours: all-pairs normal-vision dE >= 15 in the
# dataviz reference palette), dark / light neutral for tRNA / other; every gene also carries its name
COLORS = {"COX": "#e34948", "NADH": "#eda100", "ATP": "#008300", "CYTB": "#1baf7a", "rRNA": "#2a78d6",
          "tRNA": "#3d3c39", "other": "#c3c2b7"}
LEGEND = {"COX": "cytochrome c oxidase (COX)", "NADH": "NADH dehydrogenase (ND)", "ATP": "ATP synthase (ATP)",
          "CYTB": "cytochrome b (CYTB)", "rRNA": "rRNA", "tRNA": "tRNA", "other": "other"}
GC_HI, GC_LO, SKEW_POS, SKEW_NEG, DEPTH = "#52514e", "#b5b3ab", "#6250d6", "#eb6834", "#e87ba4"
TINT = 0.7                                   # opacity of the large skew / depth areas
INK, MUTED = "#0b0b0b", "#898781"
OTHER_TYPES = ("D-loop", "misc_feature", "rep_origin", "repeat_region")
AA = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H",
      "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
      "Tyr": "Y", "Val": "V"}
SIZE, R_MAX = 10.0, 1.36                     # figure inches; axes radius in data units (backbone radius 1)
BLOCK, TRACK, TRACK_GAP = 0.075, 0.12, 0.03  # gene block height on each side; ring height; space between rings
LABEL_GAP, FS_GENE, FS_TRNA = 0.045, 9.5, 8.0
LINEAR_GAP = np.radians(4)                   # between the end and the start of a linear record


def gene_name(f: SeqFeature) -> str:
    return (f.qualifiers.get("gene") or [""])[0]


def category(f: SeqFeature) -> Optional[str]:
    """Colour class of a feature; None for features that are not drawn (source, gene, intron ...)."""
    if f.type in ("tRNA", "rRNA"):
        return f.type
    if f.type == "CDS":
        g = gene_name(f).upper()
        for cat, pat in (("COX", r"COX?[123]|CO(I|II|III)"), ("NADH", r"NA?D[1-6]L?"), ("ATP", r"ATP[68]"),
                         ("CYTB", r"CYTB|COB|CYB")):
            if re.fullmatch(pat, g):
                return cat
        return "other"
    return "other" if f.type in OTHER_TYPES else None


def label(f: SeqFeature, cat: str) -> str:
    """/gene; for a tRNA its one-letter code (trnL1 -> L1, /product tRNA-Phe -> F); CR for an unnamed control region."""
    g = gene_name(f)
    if cat == "tRNA":
        m = re.match(r"trn([A-Z])(\d?)", g, re.I)
        if m and not re.match(r"trna[-_ ]?[A-Z]{3}", g, re.I):        # trnA / TRNA, not tRNA-Leu
            return m.group(1).upper() + m.group(2)
        m = re.search(r"tRNA-([A-Za-z]{3})", " ".join([g] + f.qualifiers.get("product", [])))
        return AA.get(m.group(1).capitalize(), "") if m else ""
    if not g and (f.type == "D-loop" or "control region" in " ".join(f.qualifiers.get("note", [])).lower()):
        return "CR"
    return g


def gc_windows(seq: str, circular: bool, window: int = WINDOW, step: int = STEP
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Centre (0-based), GC fraction of the A/C/G/T bases and GC skew (G - C) / (G + C) of `window` nt windows
    every `step` nt. A circular sequence has ceil(L / step) windows, the last ones wrapping across the origin; a
    linear one has full windows only (one over the whole sequence when it is shorter than `window`)."""
    s = np.frombuffer(seq.upper().encode(), dtype=np.uint8)
    n, w = len(s), min(window, len(s))
    if circular:
        s = np.concatenate([s, s[:w - 1]])
        starts = np.arange(0, n, step)
    else:
        starts = np.arange(0, n - w + 1, step)

    def count(mask: np.ndarray) -> np.ndarray:
        cs = np.concatenate([[0], np.cumsum(mask)])
        return cs[starts + w] - cs[starts]
    g, c = count(s == ord("G")), count(s == ord("C"))
    acgt = count(np.isin(s, np.frombuffer(b"ACGT", dtype=np.uint8)))
    gc = np.divide(g + c, acgt, out=np.zeros(len(starts)), where=acgt > 0)
    skew = np.divide(g - c, g + c, out=np.zeros(len(starts)), where=g + c > 0)
    return starts + w / 2, gc, skew


def theta(pos, length: int, circular: bool):
    """Angle (radians, clockwise from the top) of 0-based position(s) `pos`."""
    pos = np.asarray(pos, dtype=float)
    if circular:
        return 2 * np.pi * pos / length
    return LINEAR_GAP / 2 + (2 * np.pi - LINEAR_GAP) * pos / length


def span(parts: Sequence[Tuple[int, int]], length: int, circular: bool) -> Tuple[int, int]:
    """(start, length in nt) of the shortest arc covering all `parts`; on a circular record it may run across the
    origin (start + length > record length)."""
    ps = sorted(parts)
    start, end = ps[0][0], max(e for _, e in ps)
    if not circular:
        return start, end - start
    best, reach = (start + length - end, 0), ps[0][1]          # widest empty arc: across the origin, or a gap
    for i in range(1, len(ps)):
        best = max(best, (ps[i][0] - reach, i))
        reach = max(reach, ps[i][1])
    gap, i = best
    return (ps[i][0], length - gap) if i else (start, end - start)


def spread(targets: Sequence[float], sizes: Sequence[float]) -> List[float]:
    """Positions of labels of `sizes` near the sorted `targets` without overlap: a run of crowded labels becomes one
    block centred on the mean of its targets (least squares, order kept)."""
    blocks: List[List[int]] = []
    firsts: List[float] = []                                   # first edge of each block
    for i, t in enumerate(targets):
        blocks.append([i])
        firsts.append(t - sizes[i] / 2)
        while len(blocks) > 1 and firsts[-2] + sum(sizes[j] for j in blocks[-2]) > firsts[-1]:
            last = blocks.pop()
            blocks[-1] += last
            firsts.pop()
            offs = np.cumsum([0.0] + [sizes[j] for j in blocks[-1]])
            firsts[-1] = float(np.mean([targets[j] - offs[k] - sizes[j] / 2 for k, j in enumerate(blocks[-1])]))
    out = [0.0] * len(targets)
    for block, first in zip(blocks, firsts):
        for j in block:
            out[j] = first + sizes[j] / 2
            first += sizes[j]
    return out


def fan(angles: Sequence[float], sizes: Sequence[float]) -> List[float]:
    """spread() around the circle (radians), starting after the widest empty arc between labels."""
    n = len(angles)
    if n < 2:
        return list(angles)
    order = sorted(range(n), key=lambda i: angles[i])
    gaps = [(angles[order[(k + 1) % n]] - angles[order[k]]) % (2 * np.pi) for k in range(n)]
    cut = (int(np.argmax(gaps)) + 1) % n
    order = order[cut:] + order[:cut]
    base = angles[order[0]]
    placed = spread([base + (angles[i] - base) % (2 * np.pi) for i in order], [sizes[i] for i in order])
    out = [0.0] * n
    for i, p in zip(order, placed):
        out[i] = p
    return out


def _text_width(s: str, fontsize: float) -> float:
    """Width of `s` in points."""
    return TextPath((0, 0), s, size=fontsize).get_extents().width


def _align(phi: float, outward: bool) -> Tuple[str, str]:
    """ha / va that make upright text at angle phi grow away from its anchor, outwards or towards the centre."""
    x, y = np.sin(phi), np.cos(phi)
    if not outward:
        x, y = -x, -y
    ha = "left" if x > 0.3 else "right" if x < -0.3 else "center"
    return ha, "bottom" if y > 0.3 else "top" if y < -0.3 else "center"


def _label_text(ax, phi: float, r: float, s: str, fontsize: float, outward: bool, unit: float) -> None:
    """Label from (phi, r) outwards or inwards: names along the radius (never upside down), codes of 1-2
    characters upright, centred on the radius (so the room fan() gives them is symmetric)."""
    if len(s) <= 2:
        half = (_text_width(s, fontsize) * abs(np.sin(phi)) + fontsize * abs(np.cos(phi))) / 2 / unit
        ax.text(phi, r + half if outward else r - half, s, fontsize=fontsize, color=INK, ha="center", va="center")
        return
    a = 90 - np.degrees(phi)
    right = np.cos(np.radians(a)) > -1e-9
    ha = ("left" if right else "right") if outward else ("right" if right else "left")
    ax.text(phi, r, s, fontsize=fontsize, color=INK, rotation=a if right else a + 180, rotation_mode="anchor",
            ha=ha, va="center")


def _draw_genes(ax, rec: SeqRecord, circular: bool, unit: float) -> Tuple[List[str], float]:
    """Gene blocks and labels; returns the classes drawn and the innermost radius the labels reach."""
    n = len(rec.seq)
    cats, labels = [], []                                      # labels: angle, text, font size, outside
    for f in rec.features:
        cat = category(f) if f.location is not None else None
        if not cat:
            continue
        cats.append(cat)
        parts = [(int(p.start), int(p.end), p.strand) for p in f.location.parts]
        for s, e, strand in parts:
            bottom = 1.0 if strand == 1 else 1.0 - BLOCK if strand == -1 else 1.0 - BLOCK / 2
            ax.bar(theta(s, n, circular), BLOCK, width=theta(e, n, circular) - theta(s, n, circular), bottom=bottom,
                   align="edge", color=COLORS[cat], edgecolor="white", linewidth=0.5, zorder=3)
        start, length = span([(s, e) for s, e, _ in parts], n, circular)
        strand = f.location.strand
        if len(parts) > 1:                                     # line through introns / across the origin
            r = 1.0 + BLOCK / 2 if strand == 1 else 1.0 - BLOCK / 2 if strand == -1 else 1.0
            t = theta(np.linspace(start, start + length, 100), n, circular)
            ax.plot(t, np.full(len(t), r), color=COLORS[cat], lw=0.8, zorder=2)
        text = label(f, cat)
        if text:
            labels.append((float(theta((start + length / 2) % n, n, circular)), text,
                           FS_TRNA if cat == "tRNA" else FS_GENE, strand != -1))
    reach = 1.0 - BLOCK
    for outside in (True, False):
        side = [lb for lb in labels if lb[3] == outside]
        sign = 1 if outside else -1
        edge = 1.0 + sign * BLOCK
        r_text = edge + sign * LABEL_GAP
        # angular room: one line of radial text; an upright code its width or height, whichever is larger
        sizes = [max(1.2 * fs, _text_width(s, fs) if len(s) <= 2 else 0) / unit / r_text for _, s, fs, _ in side]
        for (t, text, fs, _), p in zip(side, fan([lb[0] for lb in side], sizes)):
            ax.plot([t, t, p], [edge, edge + sign * LABEL_GAP * 0.4, r_text - sign * 0.008], color=MUTED, lw=0.4)
            _label_text(ax, p, r_text, text, fs, outside, unit)
            if not outside:
                reach = min(reach, r_text - max(_text_width(text, fs), 1.2 * fs) / unit)
    return cats, reach


def _closed(x: np.ndarray, length: int, circular: bool, *ys: np.ndarray) -> List[np.ndarray]:
    """Window centres and values with the last window repeated before the first, so a circular ring is closed."""
    if not circular:
        return [x, *ys]
    return [np.concatenate([[x[-1] - length], x])] + [np.concatenate([[y[-1]], y]) for y in ys]


def _ring(ax, t: np.ndarray, v: np.ndarray, r0: float, colors: Tuple[str, str], base: np.ndarray,
          alpha: float = 1.0) -> None:
    """Values v (-1..1) as areas above / below the baseline circle r0, in two colours."""
    ax.fill_between(t, r0, r0 + np.clip(v, 0, None) * TRACK / 2, color=colors[0], alpha=alpha, lw=0)
    ax.fill_between(t, r0 + np.clip(v, None, 0) * TRACK / 2, r0, color=colors[1], alpha=alpha, lw=0)
    ax.plot(base, np.full(len(base), r0), color=MUTED, lw=0.3)


def _ruler(ax, r: float, length: int, circular: bool, base: np.ndarray) -> None:
    """Circle with position ticks (about 6-12), labels inside it; a linear record's ends are marked start / end."""
    step = next((s for s in (100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000) if length / s <= 12), 100000)
    ax.plot(base, np.full(len(base), r), color=MUTED, lw=0.5)
    ticks = [(p, f"{p / 1000:g} kb" if step >= 1000 else f"{p} bp") for p in range(step, length, step)]
    ticks = ([(0, "0")] if circular else [(0, "start"), (length, "end")]) + ticks
    for pos, text in ticks:
        t = float(theta(pos, length, circular))
        ax.plot([t, t], [r, r - 0.015], color=MUTED, lw=0.5)
        ha, va = _align(t, False)
        if not circular and pos in (0, length):
            ha = "left" if pos == 0 else "right"
        ax.text(t, r - 0.025, text, fontsize=7, color=MUTED, ha=ha, va=va)


def plot_record(rec: SeqRecord, out_dir: Path, depth: Optional[np.ndarray] = None,
                formats: Sequence[str] = ("png", "pdf")) -> List[Path]:
    """Draw the map of `rec` (depth: per-position read depth, e.g. from depth_from_reads) and write
    <out_dir>/<record id>.<format> for each format (png at 300 dpi); returns the written paths."""
    n = len(rec.seq)
    circular = rec.annotations.get("topology") == "circular"
    if depth is not None and len(depth) != n:
        raise ValueError(f"{rec.id}: depth has {len(depth)} values for {n} bp")
    fig = Figure(figsize=(SIZE, SIZE))
    FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1), projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, R_MAX)
    ax.set_axis_off()
    unit = 72 * SIZE / 2 / R_MAX                                  # points per data unit
    base = theta(np.linspace(0, n, 721), n, circular)            # angles of a full ring

    ax.plot(base, np.ones(len(base)), color=MUTED, lw=0.8, zorder=1)
    if not circular:                                              # end caps of a linear record
        for t in (theta(0, n, circular), theta(n, n, circular)):
            ax.plot([t, t], [1 - BLOCK - 0.02, 1 + BLOCK + 0.02], color=INK, lw=1.0)
    cats, reach = _draw_genes(ax, rec, circular, unit)

    x, gc, skew = gc_windows(str(rec.seq), circular)
    seq = str(rec.seq).upper()
    acgt = sum(seq.count(b) for b in "ACGT")
    mean_gc = (seq.count("G") + seq.count("C")) / acgt if acgt else 0.0
    dev = gc - mean_gc
    x, dev, skew = _closed(x, n, circular, dev, skew)
    t = theta(x, n, circular)
    r = reach - TRACK_GAP - TRACK / 2                            # first ring: centre radius
    _ring(ax, t, dev / max(np.abs(dev).max(), 1e-9), r, (GC_HI, GC_LO), base)
    r -= TRACK + TRACK_GAP
    _ring(ax, t, skew / max(np.abs(skew).max(), 1e-9), r, (SKEW_POS, SKEW_NEG), base, TINT)
    handles = [Patch(color=GC_HI, label=f"GC content above mean ({100 * mean_gc:.1f} %)"),
               Patch(color=GC_LO, label="GC content below mean"),
               Patch(facecolor=SKEW_POS, alpha=TINT, lw=0, label="GC skew + (G > C)"),
               Patch(facecolor=SKEW_NEG, alpha=TINT, lw=0, label="GC skew - (C > G)")]
    if depth is not None:
        r -= TRACK + TRACK_GAP
        bottom = r - TRACK / 2
        d = np.asarray(depth, dtype=float)
        edges = np.linspace(0, n, min(n, 2000) + 1).astype(int)
        binned = np.add.reduceat(d, edges[:-1]) / np.diff(edges)
        xb, binned = _closed((edges[:-1] + edges[1:]) / 2, n, circular, binned)
        top = max(d.max(), 1.0)
        ax.fill_between(theta(xb, n, circular), bottom, bottom + binned / top * TRACK, color=DEPTH, alpha=TINT, lw=0)
        ax.plot(base, np.full(len(base), bottom), color=MUTED, lw=0.3)
        ax.plot(base, np.full(len(base), bottom + d.mean() / top * TRACK), color=INK, lw=0.5, ls=(0, (3, 2)))
        handles += [Patch(facecolor=DEPTH, alpha=TINT, lw=0, label=f"read depth (max {d.max():,.0f}x)"),
                    Line2D([], [], color=INK, lw=0.5, ls=(0, (3, 2)), label=f"mean depth {d.mean():,.0f}x")]
    r_ruler = r - TRACK / 2 - 0.025
    _ruler(ax, r_ruler, n, circular, base)

    organism = rec.annotations.get("organism")
    organism = organism if organism and organism != "." else ""
    lines = [(s, 12.0, "italic") for s in textwrap.wrap(organism or rec.id, 22)]
    lines += [(f"{n:,} bp, {'circular' if circular else 'linear'}", 9.5, "normal")]
    lines += [(rec.id, 9.5, "normal")] if organism else []
    y = sum(1.4 * fs for _, fs, _ in lines) / 2                  # points above the centre
    for s, fs, style in lines:
        y -= 0.7 * fs
        ax.text(0.5, 0.5 + y / (72 * SIZE), s, transform=ax.transAxes, fontsize=fs, style=style, color=INK,
                ha="center", va="center")
        y -= 0.7 * fs

    genes = [Patch(color=COLORS[c], label=LEGEND[c]) for c in COLORS if c in cats]
    if genes:
        fig.legend(handles=genes, loc="lower left", frameon=False, fontsize=8, title_fontsize=8, alignment="left",
                   title="genes: + strand outside, - strand inside")
    fig.legend(handles=handles, loc="lower right", frameon=False, fontsize=8, title_fontsize=8, alignment="left",
               title=f"rings (GC: {WINDOW} nt windows, step {STEP} nt)")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = re.sub(r"[^\w.-]", "_", rec.id)
    written = []
    with matplotlib.rc_context({"pdf.fonttype": 42}):           # editable text in the PDF
        for fmt in formats:
            path = out_dir / f"{stem}.{fmt}"
            fig.savefig(path, dpi=300, facecolor="white")
            written.append(path)
    return written


def depth_from_reads(rec: SeqRecord, reads: Sequence, threads: int = 1) -> np.ndarray:
    """Per-position read depth of `rec` from one FASTQ file or a pair (.gz allowed), from the reference span of every
    primary alignment: minimap2 -x sr -c when it is installed (PAF; -c for exact ends; about 8x faster), else bwa mem
    (SAM, CIGAR M / D / N / = / X; bwa is in mitotools, minimap2 is not: its k8 pulls sysroot_linux-64, +0.5 GB on
    Linux). A circular record is mapped as sequence + its first ORIGIN_EXT nt and positions past its end are folded
    back onto the start, so reads across the origin count. RuntimeError when neither is found or the aligner fails."""
    mm2, bwa = mito_env.find("minimap2"), mito_env.find("bwa")
    if not (mm2 or bwa):
        raise RuntimeError("neither minimap2 nor bwa found (needed for the depth ring): install one (bioconda) or "
                           "put it on PATH")
    if not 1 <= len(reads) <= 2:
        raise ValueError(f"depth_from_reads: 1 or 2 FASTQ files, got {len(reads)}")
    seq = str(rec.seq)
    n = len(seq)
    circular = rec.annotations.get("topology") == "circular"
    ref = seq + seq[:ORIGIN_EXT] if circular else seq
    diff = np.zeros(len(ref) + 1, dtype=np.int64)
    with tempfile.TemporaryDirectory(prefix="mito_plot_") as tmp:
        fa, log = Path(tmp) / "ref.fa", Path(tmp) / "aligner.log"
        fa.write_text(f">ref\n{ref}\n")
        with open(log, "w") as err:
            if mm2:
                p = subprocess.Popen([mm2, "-x", "sr", "-c", "--secondary=no", "-t", str(threads), str(fa),
                                      *map(str, reads)], stdout=subprocess.PIPE, stderr=err, text=True)
                for line in p.stdout:                             # PAF: mapped reads only
                    f = line.rstrip("\n").split("\t")
                    if "tp:A:P" in f[12:]:
                        diff[int(f[7])] += 1
                        diff[int(f[8])] -= 1
                rc = p.wait()
            else:
                rc = subprocess.run([bwa, "index", str(fa)], stdout=err, stderr=err).returncode
                if not rc:
                    p = subprocess.Popen([bwa, "mem", "-t", str(threads), str(fa), *map(str, reads)],
                                         stdout=subprocess.PIPE, stderr=err, text=True)
                    for line in p.stdout:                         # SAM: unmapped reads have a line too
                        f = line.split("\t", 6)
                        if line[0] == "@" or int(f[1]) & 0x904:   # unmapped, secondary, supplementary
                            continue
                        start = int(f[3]) - 1
                        diff[start] += 1
                        diff[start + sum(int(k) for k, op in CIGAR.findall(f[5]) if op in "MDN=X")] -= 1
                    rc = p.wait()
        if rc:
            raise RuntimeError(f"{Path(mm2 or bwa).name} failed (exit {rc}): {log.read_text().strip()[-1000:]}")
    cov = np.cumsum(diff[:-1])
    depth = cov[:n].copy()
    if circular:
        depth[:len(ref) - n] += cov[n:]
    return depth


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="GenBank file, or a directory of *.gb / *.gbk / *.gbf files")
    p.add_argument("-o", "--out_dir", required=True)
    p.add_argument("--fq1", help="reads for the depth ring (FASTQ, .gz allowed); only with a single record")
    p.add_argument("--fq2", help="mates of --fq1")
    p.add_argument("--threads", type=int, default=1, help="aligner threads (minimap2, else bwa mem; default 1)")
    p.add_argument("--formats", default="png,pdf", help="comma-separated (default png,pdf)")
    a = p.parse_args(argv)
    src = Path(a.input)
    if not src.exists():
        p.error(f"{src} not found")
    files = sorted(f for f in src.iterdir() if f.suffix in (".gb", ".gbk", ".gbf")) if src.is_dir() else [src]
    recs = [r for f in files for r in SeqIO.parse(str(f), "genbank")]
    if not recs:
        p.error(f"no GenBank records in {src}")
    if a.fq2 and not a.fq1:
        p.error("--fq2 needs --fq1")
    reads = [x for x in (a.fq1, a.fq2) if x]
    if reads and len(recs) != 1:
        p.error(f"--fq1 / --fq2 need an input with one record ({len(recs)} found)")
    formats = [x.strip() for x in a.formats.split(",") if x.strip()]
    for rec in recs:
        try:
            depth = depth_from_reads(rec, reads, a.threads) if reads else None
        except RuntimeError as exc:
            print(f"[mito_plot] {rec.id}: {exc}", file=sys.stderr)
            return 1
        for path in plot_record(rec, Path(a.out_dir), depth, formats):
            print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
