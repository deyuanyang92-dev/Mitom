#!/usr/bin/env python3
"""mito_plot.py - circular map of an annotated mitochondrial genome (MitoM), drawn with matplotlib.

Two styles. circos (default, --style circos), after mitoPlotR (K. Onthank, github.com/KirtOnthank/mitoPlotR, MIT
licence): a framed lane per strand (+ outside, - inside) with blocks coloured by type (protein gene, tRNA, rRNA,
other), gene names along the circle outside, tRNA codes, a position axis, GC content, read depth in 25 nt bins
scaled to the 99th percentile with the mean as a dashed line, the organism in the centre.

mitom (--style mitom; the default until 2026-10-30), from the outside in:
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
  mitom plot -i <GenBank file or folder> -o <out_dir> [--fq1 R1.fq.gz [--fq2 R2.fq.gz]] [--threads N]
             [--formats png,pdf] [--style circos|mitom]         (= python3 mito_plot.py ...)
Every record of the GenBank file, or of every GenBank file in the folder and its subfolders (.gb .gbk .gbf .genbank
.gbff), becomes <out_dir>/<record id>.png (300 dpi) and .pdf; reads are used only when the input holds one record.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
import textwrap
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib
import numpy as np
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath, text_to_path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402
import mito_names  # noqa: E402
from mito_common import command_line, iter_gb_files, old_spelling, STDERR_MARK, STDERR_TAIL, TOOL_LOG  # noqa: E402
from mito_common import cli_guard, InputError, prepare_input, require_records  # noqa: E402
from mito_ncr import MIN_LEN as MIN_NONCODING, extent, is_ncr, noncoding, span  # noqa: E402,F401

WINDOW, STEP, ORIGIN_EXT = 100, 20, 1000
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
# five hues that stay apart pairwise (any two classes can be neighbours: all-pairs normal-vision dE >= 15 in the
# dataviz reference palette), dark / light neutral for tRNA / other; every gene also carries its name
COLORS = {"COX": "#e34948", "NADH": "#eda100", "ATP": "#008300", "CYTB": "#1baf7a", "rRNA": "#2a78d6",
          "tRNA": "#3d3c39", "other": "#c3c2b7"}
LEGEND = {"COX": "cytochrome c oxidase (COX)", "NADH": "NADH dehydrogenase (ND)", "ATP": "ATP synthase (ATP)",
          "CYTB": "cytochrome b (CYTB)", "rRNA": "rRNA", "tRNA": "tRNA", "other": "other"}
GC_HI, GC_LO, SKEW_POS, SKEW_NEG, DEPTH = "#52514e", "#b5b3ab", "#6250d6", "#eb6834", "#e87ba4"
NONCODING = "#e1e0d9"                         # unannotated stretch: lighter than 'other' (annotated control region)
TINT = 0.7                                   # opacity of the large skew / depth areas
INK, MUTED = "#0b0b0b", "#898781"
OTHER_TYPES = ("D-loop", "misc_feature", "rep_origin", "repeat_region")
AA = mito_names.AA3
SIZE, R_MAX = 10.0, 1.36                     # figure inches; axes radius in data units (backbone radius 1)
BLOCK, TRACK, TRACK_GAP = 0.075, 0.12, 0.03  # gene block height on each side; ring height; space between rings
LABEL_GAP, FS_GENE, FS_TRNA = 0.045, 9.5, 8.0
LINEAR_GAP = np.radians(4)                   # between the end and the start of a linear record
TICK_ROOM = 0.4                              # no position tick closer than this many steps to the end (its label
                                             # would cover 0 / end: 16 kb of a 16,063 bp record)
MIN_INTRON = 142                             # a gap between the parts of a joined gene drawn as an intron when it has
                                             # no intron feature (= mito_intron.MIN_INTRON, MFannot's --minintron;
                                             # shorter gaps are frameshift joins, e.g. bird nad3)
HATCH = "////"


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
    if is_ncr(f):                                # MitoM's NCR features: drawn as non-coding regions
        return None
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


def introns(rec: SeqRecord) -> List[Tuple[int, int, int, str]]:
    """(start, end, strand, label) of every intron: the intron features, and the gaps of >= MIN_INTRON nt between
    consecutive parts of a joined CDS / rRNA / tRNA on one strand that no intron feature covers (records of tools
    that write the join only). end > record length for an intron across the origin of a circular record; label
    '<gene> intron[ k] <length> bp' (k when the gene has several)."""
    n = len(rec.seq)
    circular = rec.annotations.get("topology") == "circular"
    found: List[Tuple[int, int, int, str]] = []
    for f in rec.features:
        if f.type == "intron" and f.location is not None:
            s, length = extent(f.location, n, circular)
            found.append((s, s + length, f.location.strand or 1, gene_name(f)))
    explicit = list(found)
    for f in rec.features:
        if f.type not in ("CDS", "rRNA", "tRNA") or f.location is None or len(f.location.parts) < 2:
            continue
        parts = list(f.location.parts)                            # transcript order
        for a, b in zip(parts, parts[1:]):
            if a.strand != b.strand:                              # trans-spliced
                continue
            left, right = (a, b) if a.strand != -1 else (b, a)    # the two parts in genome order
            s, e = int(left.end), int(right.start)
            if e >= s:
                length = e - s
            elif circular and int(right.end) <= int(left.start):  # the join runs across the origin
                length = n - s + e
            else:                                                 # parts overlap: a frameshift join
                continue
            if length < MIN_INTRON:
                continue
            if any(_covered(x, y, s, s + length, n) >= 0.5 * length for x, y, _, _ in explicit):
                continue                                          # given by an intron feature (one nested in the
                #                                                   gap, e.g. a hexacoral cox1 intron inside the nad5
                #                                                   intron, does not count: eval_introns.py's rule)
            found.append((s, s + length, a.strand or 1, gene_name(f)))
    count: Dict[str, int] = {}
    for _, _, _, g in found:
        count[g] = count.get(g, 0) + 1
    seen: Dict[str, int] = {}
    out = []
    for s, e, strand, g in sorted(found):
        seen[g] = seen.get(g, 0) + 1
        k = f" {seen[g]}" if count[g] > 1 else ""
        out.append((s, e, strand, f"{g + ' ' if g else ''}intron{k} {e - s:,} bp"))
    return out


def _covered(s1: int, e1: int, s2: int, e2: int, n: int) -> int:
    """nt of the arc [s2, e2) that the arc [s1, e1) covers on a circle of n (ends may exceed n)."""
    return max(max(0, min(e1 + k, e2) - max(s1 + k, s2)) for k in (0, n, -n))


def theta(pos, length: int, circular: bool):
    """Angle (radians, clockwise from the top) of 0-based position(s) `pos`."""
    pos = np.asarray(pos, dtype=float)
    if circular:
        return 2 * np.pi * pos / length
    return LINEAR_GAP / 2 + (2 * np.pi - LINEAR_GAP) * pos / length


def name_pos(parts: Sequence[Tuple[int, int, Optional[int]]], start: int, length: int, n: int) -> float:
    """Where a feature's name goes: the middle of its extent, or of its largest part when it has an intron (the
    middle of a 10-kb intron is far from both exons); a gene merely cut by the origin keeps the middle."""
    if len(parts) > 1 and sum(e - s for s, e, _ in parts) < length:
        s, e, _ = max(parts, key=lambda p: p[1] - p[0])
        return (s + e) / 2
    return (start + length / 2) % n


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


def _label_text(ax, phi: float, r: float, s: str, fontsize: float, outward: bool, unit: float,
                color: str = INK) -> None:
    """Label from (phi, r) outwards or inwards: names along the radius (never upside down), codes of 1-2
    characters upright, centred on the radius (so the room fan() gives them is symmetric)."""
    if len(s) <= 2:
        half = (_text_width(s, fontsize) * abs(np.sin(phi)) + fontsize * abs(np.cos(phi))) / 2 / unit
        ax.text(phi, r + half if outward else r - half, s, fontsize=fontsize, color=color, ha="center", va="center")
        return
    a = 90 - np.degrees(phi)
    right = np.cos(np.radians(a)) > -1e-9
    ha = ("left" if right else "right") if outward else ("right" if right else "left")
    ax.text(phi, r, s, fontsize=fontsize, color=color, rotation=a if right else a + 180, rotation_mode="anchor",
            ha=ha, va="center")


def _draw_genes(ax, rec: SeqRecord, circular: bool, unit: float, min_noncoding: int = MIN_NONCODING
                ) -> Tuple[List[str], float]:
    """Gene blocks, introns (hatched), unannotated stretches and labels; returns the classes drawn ('intron',
    'noncoding' included) and the innermost radius the labels reach."""
    n = len(rec.seq)
    cats, labels = [], []                                      # labels: angle, text, font size, outside, colour
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
        start, length = extent(f.location, n, circular)
        strand = f.location.strand
        if len(parts) > 1:                                     # line through introns / across the origin
            r = 1.0 + BLOCK / 2 if strand == 1 else 1.0 - BLOCK / 2 if strand == -1 else 1.0
            t = theta(np.linspace(start, start + length, 100), n, circular)
            ax.plot(t, np.full(len(t), r), color=COLORS[cat], lw=0.8, zorder=2)
        text = label(f, cat)
        if text:
            labels.append((float(theta(name_pos(parts, start, length, n), n, circular)), text,
                           FS_TRNA if cat == "tRNA" else FS_GENE, strand != -1, INK))
    for s, e, strand, text in introns(rec):                    # hatched in the gene's lane, labelled with length
        cats.append("intron")
        ax.bar(theta(s, n, circular), BLOCK, width=theta(e, n, circular) - theta(s, n, circular),
               bottom=1.0 if strand != -1 else 1.0 - BLOCK, align="edge", facecolor="white", edgecolor=MUTED,
               hatch=HATCH, linewidth=0, zorder=2.5)
        # on the other side of the circle than the gene's name (both would sit at the middle of the gene)
        labels.append((float(theta(((s + e) / 2) % n, n, circular)), text, FS_TRNA, strand == -1, MUTED))
    for k, (s, e) in enumerate(noncoding(rec, min_noncoding), 1):   # across the backbone (strand unknown);
        #                                                            numbered as in the GenBank (mito_ncr.add_ncr)
        cats.append("noncoding")
        ax.bar(theta(s, n, circular), BLOCK, width=theta(e, n, circular) - theta(s, n, circular),
               bottom=1.0 - BLOCK / 2, align="edge", color=NONCODING, linewidth=0, zorder=2.5)
        labels.append((float(theta(((s + e) / 2) % n, n, circular)), f"NCR{k} {e - s:,} bp", FS_TRNA, True, MUTED))
    reach = 1.0 - BLOCK
    for outside in (True, False):
        side = [lb for lb in labels if lb[3] == outside]
        sign = 1 if outside else -1
        edge = 1.0 + sign * BLOCK
        r_text = edge + sign * LABEL_GAP
        # angular room: one line of radial text; an upright code its width or height, whichever is larger
        sizes = [max(1.2 * fs, _text_width(s, fs) if len(s) <= 2 else 0) / unit / r_text for _, s, fs, _, _ in side]
        for (t, text, fs, _, color), p in zip(side, fan([lb[0] for lb in side], sizes)):
            ax.plot([t, t, p], [edge, edge + sign * LABEL_GAP * 0.4, r_text - sign * 0.008], color=MUTED, lw=0.4)
            _label_text(ax, p, r_text, text, fs, outside, unit, color)
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
    ticks = [(p, f"{p / 1000:g} kb" if step >= 1000 else f"{p} bp") for p in range(step, length, step)
             if length - p >= TICK_ROOM * step]                 # (not on top of 0 / end)
    ticks = ([(0, "0")] if circular else [(0, "start"), (length, "end")]) + ticks
    for pos, text in ticks:
        t = float(theta(pos, length, circular))
        ax.plot([t, t], [r, r - 0.015], color=MUTED, lw=0.5)
        ha, va = _align(t, False)
        if not circular and pos in (0, length):
            ha = "left" if pos == 0 else "right"
        ax.text(t, r - 0.025, text, fontsize=7, color=MUTED, ha=ha, va=va)


STYLES = ("circos", "mitom")


def plot_record(rec: SeqRecord, out_dir: Path, depth: Optional[np.ndarray] = None,
                formats: Sequence[str] = ("png", "pdf"), style: str = "circos",
                min_noncoding: int = MIN_NONCODING, used: Optional[set] = None) -> List[Path]:
    """Draw the map of `rec` (depth: per-position read depth, e.g. from depth_from_reads) in `style` - circos (strand
    lanes, after mitoPlotR) or mitom (genes on both sides of one circle, GC content and GC skew rings) - and write
    <out_dir>/<record id>.<format> for each format (png at 300 dpi); returns the written paths. Introns are hatched
    and labelled with their length, unannotated stretches of >= min_noncoding nt (0: none) shaded and labelled
    'NCR <length>'."""
    n = len(rec.seq)
    if depth is not None and len(depth) != n:
        raise ValueError(f"{rec.id}: depth has {len(depth)} values for {n} bp")
    if style not in STYLES:
        raise ValueError(f"unknown map style {style!r} ({', '.join(STYLES)})")
    fig = (_circos_figure if style == "circos" else _mitom_figure)(rec, depth, min_noncoding)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = base = re.sub(r"[^\w.-]", "_", rec.id)
    if used is not None:                          # records with one ID (two input files): two maps, not one overwritten
        k = 2
        while stem in used:
            stem, k = f"{base}_{k}", k + 1
        used.add(stem)
    written = []
    with matplotlib.rc_context({"pdf.fonttype": 42}):           # editable text in the PDF
        for fmt in formats:
            path = out_dir / f"{stem}.{fmt}"
            fig.savefig(path, dpi=300, facecolor="white")
            written.append(path)
    return written


def _polar(n: int, circular: bool):
    """Figure with one polar axes over it (position 1 at the top, clockwise; backbone radius 1), points per data
    unit and the angles of a full ring."""
    fig = Figure(figsize=(SIZE, SIZE))
    FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1), projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, R_MAX)
    ax.set_axis_off()
    return fig, ax, 72 * SIZE / 2 / R_MAX, theta(np.linspace(0, n, 721), n, circular)


def _mitom_figure(rec: SeqRecord, depth: Optional[np.ndarray], min_noncoding: int = MIN_NONCODING) -> Figure:
    """The mitom style: genes on both sides of one circle, GC content and GC skew rings, read depth, a kb ruler."""
    n = len(rec.seq)
    circular = rec.annotations.get("topology") == "circular"
    fig, ax, unit, base = _polar(n, circular)

    ax.plot(base, np.ones(len(base)), color=MUTED, lw=0.8, zorder=1)
    if not circular:                                              # end caps of a linear record
        for t in (theta(0, n, circular), theta(n, n, circular)):
            ax.plot([t, t], [1 - BLOCK - 0.02, 1 + BLOCK + 0.02], color=INK, lw=1.0)
    cats, reach = _draw_genes(ax, rec, circular, unit, min_noncoding)

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
    genes += [Patch(facecolor="white", edgecolor=MUTED, hatch=HATCH, lw=0, label="intron")] if "intron" in cats else []
    genes += [Patch(color=NONCODING, label=f"non-coding region (NCR, unannotated, >= {min_noncoding} nt)")] \
        if "noncoding" in cats else []
    if genes:
        fig.legend(handles=genes, loc="lower left", frameon=False, fontsize=8, title_fontsize=8, alignment="left",
                   title="genes: + strand outside, - strand inside")
    fig.legend(handles=handles, loc="lower right", frameon=False, fontsize=8, title_fontsize=8, alignment="left",
               title=f"rings (GC: {WINDOW} nt windows, step {STEP} nt)")
    return fig


# ── circos style: after mitoPlotR (K. Onthank, MIT licence, github.com/KirtOnthank/mitoPlotR) ────────────────────
# Its layout - a framed lane per strand (+ outside, - inside) with feature blocks coloured by type, names written
# along the circle outside, a position axis, binned read depth, the species in the centre - with its green / blue /
# red replaced by the validated hues below (its red / green pair is too close for deuteranopia: dE 7.8, these 9.2),
# and tRNA codes, GC content and the mean depth added.
CIRCOS_COLORS = {"gene": "#1baf7a", "tRNA": "#2a78d6", "rRNA": "#eb6834", "other": "#c3c2b7"}
CIRCOS_LEGEND = {"gene": "protein gene", "tRNA": "tRNA", "rRNA": "rRNA", "other": "other (control region ...)"}
LANE, LANE_GAP, EDGE = 0.07, 0.012, 0.4       # strand lane height; between the lanes; outline width (pt)
R_PLUS = 1.0                                  # + strand lane: R_PLUS .. R_PLUS + LANE
R_MINUS = R_PLUS - LANE_GAP - LANE            # - strand lane: R_MINUS .. R_MINUS + LANE
FS_NAME, FS_CODE, FS_TICK = 10.0, 7.5, 7.5
DEPTH_BIN, DEPTH_CAP = 25, 0.99               # read depth: mean of 25 nt bins, scaled to the 99th percentile
DEPTH_INK, MEAN_LINE = "#52514e", "#e34948"


def _advance(s: str, fontsize: float) -> float:
    """Advance width of `s` in points (letters placed one by one along a circle)."""
    return text_to_path.get_text_width_height_descent(s, FontProperties(size=fontsize), ismath=False)[0] if s else 0.0


def _curved_text(ax, phi: float, r: float, s: str, fontsize: float, unit: float, outward: bool = True,
                 color: str = INK) -> float:
    """`s` written along the circle of radius r and centred on angle phi, as circos' bending text: left to right on
    the upper half (letters upright) and on the lower half (turned), the letters outside r (outward) or inside it.
    Returns its angular width."""
    upper = np.cos(phi) >= 0
    total = _advance(s, fontsize)
    va = ("bottom" if upper else "top") if outward else ("top" if upper else "bottom")
    for i, ch in enumerate(s):
        if ch.isspace():
            continue
        x = (_advance(s[:i], fontsize) + _advance(s[:i + 1], fontsize)) / 2 - total / 2   # from the label centre
        a = phi + x / (r * unit) if upper else phi - x / (r * unit)
        ax.text(a, r, ch, fontsize=fontsize, color=color, rotation=-np.degrees(a) + (0 if upper else 180),
                rotation_mode="anchor", ha="center", va=va)
    return total / (r * unit)


def _circos_type(cat: str) -> str:
    return "gene" if cat in ("COX", "NADH", "ATP", "CYTB") else cat


def _circos_features(ax, rec: SeqRecord, circular: bool, unit: float, min_noncoding: int = MIN_NONCODING
                     ) -> List[str]:
    """Strand lanes with the feature blocks, introns (hatched), unannotated stretches, tRNA codes and curved names
    outside; returns the types drawn ('intron', 'noncoding' included)."""
    n = len(rec.seq)
    base = theta(np.linspace(0, n, 721), n, circular)
    for r0 in (R_PLUS, R_MINUS):                                  # the lanes: white, framed
        ax.fill_between(base, r0, r0 + LANE, color="white", zorder=1)
        for r in (r0, r0 + LANE):
            ax.plot(base, np.full(len(base), r), color=INK, lw=EDGE, zorder=2)
    if not circular:                                              # a linear record's ends close the lanes
        for t in (theta(0, n, circular), theta(n, n, circular)):
            for r0 in (R_PLUS, R_MINUS):
                ax.plot([t, t], [r0, r0 + LANE], color=INK, lw=EDGE, zorder=2)
    types, names, codes = [], [], []                              # names: (angle, text, size, colour); codes
    for f in rec.features:
        cat = category(f) if f.location is not None else None
        if not cat:
            continue
        typ = _circos_type(cat)
        types.append(typ)
        parts = [(int(p.start), int(p.end), p.strand) for p in f.location.parts]
        for s, e, strand in parts:
            bottom, height = (R_PLUS, LANE) if strand == 1 else (R_MINUS, LANE) if strand == -1 else \
                (R_MINUS, 2 * LANE + LANE_GAP)                    # no strand (control region): both lanes
            ax.bar(theta(s, n, circular), height, width=theta(e, n, circular) - theta(s, n, circular), bottom=bottom,
                   align="edge", color=CIRCOS_COLORS[typ], edgecolor=INK, linewidth=EDGE, zorder=3)
        start, length = extent(f.location, n, circular)
        if len(parts) > 1:                                        # intron / origin: the parts joined by a line
            strand = f.location.strand
            r = R_PLUS + LANE / 2 if strand == 1 else R_MINUS + LANE / 2 if strand == -1 else R_PLUS - LANE_GAP / 2
            t = theta(np.linspace(start, start + length, 100), n, circular)
            ax.plot(t, np.full(len(t), r), color=INK, lw=EDGE, zorder=2)
        text = label(f, cat)
        if text:
            phi = float(theta(name_pos(parts, start, length, n), n, circular))
            if typ == "tRNA":
                codes.append((phi, text))
            else:
                names.append((phi, text, FS_NAME, INK))
    for s, e, strand, text in introns(rec):                       # hatched in the gene's lane, labelled
        types.append("intron")
        ax.bar(theta(s, n, circular), LANE, width=theta(e, n, circular) - theta(s, n, circular),
               bottom=R_PLUS if strand != -1 else R_MINUS, align="edge", facecolor="white", edgecolor=MUTED,
               hatch=HATCH, linewidth=0, zorder=2.5)
        names.append((float(theta(((s + e) / 2) % n, n, circular)), text, FS_CODE, MUTED))
    for k, (s, e) in enumerate(noncoding(rec, min_noncoding), 1):   # both lanes (strand unknown)
        types.append("noncoding")
        ax.bar(theta(s, n, circular), 2 * LANE + LANE_GAP, width=theta(e, n, circular) - theta(s, n, circular),
               bottom=R_MINUS, align="edge", color=NONCODING, linewidth=0, zorder=2.5)
        names.append((float(theta(((s + e) / 2) % n, n, circular)), f"NCR{k} {e - s:,} bp", FS_CODE, MUTED))
    top = R_PLUS + LANE
    r_code = top + 0.014                                          # tRNA codes: upright, spread, with leader lines
    sizes = [(max(_advance(c, FS_CODE), FS_CODE) + 2) / unit / r_code for _, c in codes]
    for (t, c), p in zip(codes, fan([t for t, _ in codes], sizes)):
        ax.plot([t, t, p], [top, top + 0.004, r_code - 0.002], color=MUTED, lw=0.4, zorder=2)
        half = (_advance(c, FS_CODE) * abs(np.sin(p)) + FS_CODE * abs(np.cos(p))) / 2 / unit
        ax.text(p, r_code + half, c, fontsize=FS_CODE, color=INK, ha="center", va="center")
    r_name = r_code + 1.2 * FS_CODE / unit + 0.012                # names: outside the codes, along the circle
    sizes = [(_advance(x, fs) + 6) / unit / r_name for _, x, fs, _ in names]
    for (t, x, fs, color), p in zip(names, fan([t for t, _, _, _ in names], sizes)):
        if abs((p - t + np.pi) % (2 * np.pi) - np.pi) * r_name * unit > 2:   # moved: a leader line to the gene
            ax.plot([t, t, p], [top, r_code - 0.004, r_name - 0.004], color=MUTED, lw=0.4, zorder=2)
        _curved_text(ax, p, r_name, x, fs, unit, color=color)
    return types


def _circos_axis(ax, r: float, n: int, circular: bool, unit: float) -> None:
    """Position axis (mitoPlotR's circos.axis): a circle, about 8 ticks pointing inwards, numbers along the circle."""
    step = next((s for s in (500, 1000, 2000, 5000, 10000, 20000, 50000) if n / s <= 9), 100000)
    base = theta(np.linspace(0, n, 721), n, circular)
    ax.plot(base, np.full(len(base), r), color=INK, lw=0.8)
    for pos in [p for p in range(0, n, step) if p == 0 or n - p >= TICK_ROOM * step] + ([] if circular else [n]):
        t = float(theta(pos, n, circular))
        ax.plot([t, t], [r, r - 0.018], color=INK, lw=0.8)
        _curved_text(ax, t, r - 0.026, str(pos), FS_TICK, unit, outward=False)


def _circos_figure(rec: SeqRecord, depth: Optional[np.ndarray], min_noncoding: int = MIN_NONCODING) -> Figure:
    """The circos style (after mitoPlotR): strand lanes, names along the circle, position axis, GC content and
    read depth inside, the species in the centre."""
    n = len(rec.seq)
    circular = rec.annotations.get("topology") == "circular"
    fig, ax, unit, base = _polar(n, circular)
    types = _circos_features(ax, rec, circular, unit, min_noncoding)
    r_axis = R_MINUS - 0.022
    _circos_axis(ax, r_axis, n, circular, unit)

    x, gc, _ = gc_windows(str(rec.seq), circular)                 # GC content: area over the track's baseline
    x, gc = _closed(x, n, circular, gc)
    lo, hi = float(gc.min()), float(gc.max())
    r0, h = r_axis - 0.075 - TRACK, TRACK
    ax.fill_between(theta(x, n, circular), r0, r0 + (gc - lo) / max(hi - lo, 1e-9) * h, color=GC_LO, lw=0)
    ax.plot(base, np.full(len(base), r0), color=MUTED, lw=0.4)
    handles = [Patch(facecolor=GC_LO, lw=0, label=f"GC content {100 * lo:.0f}-{100 * hi:.0f} % ({WINDOW} nt windows)")]
    inner = r0
    if depth is not None:                                         # read depth: bars of DEPTH_BIN nt, capped
        d = np.asarray(depth, dtype=float)
        edges = np.arange(0, n + DEPTH_BIN, DEPTH_BIN).clip(max=n)
        edges = edges[np.r_[True, np.diff(edges) > 0]]
        binned = np.add.reduceat(d, edges[:-1]) / np.diff(edges)
        cap = max(float(np.quantile(binned, DEPTH_CAP)), 1.0)
        r0 = inner - 0.03 - 2 * TRACK
        top = r0 + np.minimum(binned / cap, 1) * 2 * TRACK
        ax.fill_between(theta(edges, n, circular), r0, np.append(top, top[-1]), step="post", color=DEPTH_INK, lw=0)
        ax.plot(base, np.full(len(base), r0 + min(d.mean() / cap, 1) * 2 * TRACK), color=MEAN_LINE, lw=0.8,
                ls=(0, (4, 3)))
        handles += [Patch(facecolor=DEPTH_INK, lw=0, label=f"read depth ({DEPTH_BIN} nt bins, top = {cap:,.0f}x)"),
                    Line2D([], [], color=MEAN_LINE, lw=0.8, ls=(0, (4, 3)), label=f"mean depth {d.mean():,.0f}x")]
        inner = r0

    organism = rec.annotations.get("organism")
    organism = organism if organism and organism != "." else ""
    lines = [(s, 17.0, "italic", INK) for s in textwrap.wrap(organism or rec.id, 20)]
    lines += [("mitochondrial genome", 11.0, "normal", INK),
              (f"{n:,} bp, {'circular' if circular else 'linear'}", 11.0, "normal", INK)]
    lines += [(rec.id, 9.0, "normal", MUTED)] if organism else []
    y = sum(1.45 * fs for _, fs, _, _ in lines) / 2              # points above the centre
    for s, fs, style, color in lines:
        y -= 0.725 * fs
        ax.text(0.5, 0.5 + y / (72 * SIZE), s, transform=ax.transAxes, fontsize=fs, style=style, color=color,
                ha="center", va="center")
        y -= 0.725 * fs
    keys = [t for t in CIRCOS_COLORS if t in types]
    if keys:
        fig.legend(handles=[Patch(facecolor=CIRCOS_COLORS[t], edgecolor=INK, lw=EDGE, label=CIRCOS_LEGEND[t])
                            for t in keys]
                   + ([Patch(facecolor="white", edgecolor=MUTED, hatch=HATCH, lw=0, label="intron")]
                      if "intron" in types else [])
                   + ([Patch(color=NONCODING, label=f"non-coding region (NCR, unannotated, >= {min_noncoding} nt)")]
                      if "noncoding" in types else []), loc="lower left", frameon=False, fontsize=8, title_fontsize=8,
                   alignment="left", title="+ strand: outer lane, \u2212 strand: inner lane")
    fig.legend(handles=handles, loc="lower right", frameon=False, fontsize=8, alignment="left")
    return fig


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
        t0 = time.time()
        with open(log, "w") as err:
            if mm2:
                cmd = [mm2, "-x", "sr", "-c", "--secondary=no", "-t", str(threads), str(fa), *map(str, reads)]
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, text=True)
                for line in p.stdout:                             # PAF: mapped reads only
                    f = line.rstrip("\n").split("\t")
                    if "tp:A:P" in f[12:]:
                        diff[int(f[7])] += 1
                        diff[int(f[8])] -= 1
                rc = p.wait()
            else:
                cmd = [bwa, "index", str(fa)]
                rc = subprocess.run(cmd, stdout=err, stderr=err).returncode
                if not rc:
                    cmd = [bwa, "mem", "-t", str(threads), str(fa), *map(str, reads)]
                    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, text=True)
                    for line in p.stdout:                         # SAM: unmapped reads have a line too
                        f = line.split("\t", 6)
                        if line[0] == "@" or int(f[1]) & 0x904:   # unmapped, secondary, supplementary
                            continue
                        start = int(f[3]) - 1
                        diff[start] += 1
                        diff[start + sum(int(k) for k, op in CIGAR.findall(f[5]) if op in "MDN=X")] -= 1
                    rc = p.wait()
        tail = log.read_text().strip().splitlines()[-STDERR_TAIL:] if rc else []
        TOOL_LOG.debug(f"[CMD rc={rc} {time.time() - t0:.1f}s] {command_line(cmd)}" + "".join(f"\n{STDERR_MARK}{ln}"
                                                                                       for ln in tail))
        if rc:
            raise RuntimeError(f"{Path(mm2 or bwa).name} failed (exit {rc}): {log.read_text().strip()[-1000:]}")
    cov = np.cumsum(diff[:-1])
    depth = cov[:n].copy()
    if circular:
        depth[:len(ref) - n] += cov[n:]
    return depth


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="GenBank file, or a folder of GenBank files (subfolders included)")
    p.add_argument("-o", "--out", dest="out_dir", metavar="OUT", required=True, help="folder for the maps")
    old_spelling(p, "--out_dir", "out_dir")
    p.add_argument("--fq1", help="reads for the depth ring (FASTQ, .gz allowed); only with a single record")
    p.add_argument("--fq2", help="mates of --fq1")
    p.add_argument("--threads", type=int, default=1, help="aligner threads (minimap2, else bwa mem; default 1)")
    p.add_argument("--formats", default="png,pdf", help="comma-separated (default png,pdf)")
    p.add_argument("--style", choices=STYLES, default="circos",
                   help="circos (default, after mitoPlotR): a framed lane per strand, names along the circle, a "
                        "position axis, GC content, read depth; mitom: genes on both sides of one circle, GC content "
                        "and GC skew rings")
    p.add_argument("--min_noncoding", type=int, default=MIN_NONCODING,
                   help=f"shortest unannotated stretch shaded and labelled 'NCR <length>' (default {MIN_NONCODING} "
                        "nt; 0: none)")
    a = p.parse_args(argv)
    src = Path(prepare_input(a.input))
    recs = [r for f in iter_gb_files(src) for r in SeqIO.parse(str(f), "genbank")]
    if not recs:
        p.error(f"no GenBank records in {src}")
    if a.fq2 and not a.fq1:
        p.error("--fq2 needs --fq1")
    reads = [x for x in (a.fq1, a.fq2) if x]
    if reads and len(recs) != 1:
        p.error(f"--fq1 / --fq2 need an input with one record ({len(recs)} found)")
    formats = [x.strip() for x in a.formats.split(",") if x.strip()]
    used_stems: set = set()
    for rec in recs:
        try:
            depth = depth_from_reads(rec, reads, a.threads) if reads else None
        except RuntimeError as exc:
            print(f"[mito_plot] {rec.id}: {exc}", file=sys.stderr)
            return 1
        for path in plot_record(rec, Path(a.out_dir), depth, formats, a.style, a.min_noncoding, used=used_stems):
            print(path)
    return 0


if __name__ == "__main__":
    sys.exit(cli_guard(main))
