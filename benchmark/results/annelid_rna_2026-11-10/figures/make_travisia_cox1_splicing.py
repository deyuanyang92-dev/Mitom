#!/usr/bin/env python3
"""Figure: how splicing of the Travisia forbesii COX1 group II intron is judged from RNA-seq (OZ479805, ERR17229609)
against DNA of the same specimen (10x Illumina ERR17229605).

Coverage alone cannot show splicing - unspliced precursor RNA covers the intron too. Splicing is shown by junction
reads: reads whose sequence runs from the last base of exon 1 straight into the first base of exon 2 (counted as the
exact 30-nt exon|exon sequence in the reads, travisia_dna_check.py). DNA has the intron and no such reads.
Counts below: the same 4 M RNA read pairs as depth_rna.tsv; all 10x mitochondrial reads for DNA.

Usage: python3 make_travisia_cox1_splicing.py   (writes travisia_cox1_splicing.png / .svg next to this script)
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Rectangle

HERE = Path(__file__).parent
RNA, DNA = "#2a78d6", "#eb6834"          # categorical slots 1 / 2 (dataviz reference palette)
INK, MUTED, GRID = "#2b2b29", "#6b6a63", "#e4e3dd"
EX1, INTRON, EX2 = (1, 273), (274, 1525), (1526, 2800)
COUNTS = {"rna": {"spliced": 5799, "ex_in": 2991, "in_ex": 4064},
          "dna": {"spliced": 0, "ex_in": 141, "in_ex": 98}}
for f in ("Noto Sans CJK SC", "Noto Sans CJK JP", "Noto Serif CJK SC", "Noto Serif CJK JP"):
    if any(f == x.name for x in matplotlib.font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break


def depth(name):
    xs, ys = [], []
    for ln in open(HERE / f"depth_{name}.tsv"):
        _, p, d = ln.split("\t")
        xs.append(int(p))
        ys.append(int(d))
    return xs, ys


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def coverage(ax, name, color, title):
    xs, ys = depth(name)
    ax.fill_between(xs, ys, color=color, alpha=0.18, linewidth=0)
    ax.plot(xs, ys, color=color, linewidth=2)
    top = max(ys) * 1.75
    ax.set_ylim(0, top)
    ax.set_xlim(0, 2900)
    ax.axvspan(*INTRON, color=GRID, alpha=0.5, linewidth=0)
    c = COUNTS[name]
    intron_depth = sorted(ys[INTRON[0]:INTRON[1]])[(INTRON[1] - INTRON[0]) // 2]
    ly = top * 0.36 if intron_depth < 0.25 * top else top * 0.03      # above a low intron line, else at the bottom
    # sashimi arc: exon 1 end -> exon 2 start, labelled with the spliced junction reads
    mid, width = (EX1[1] + EX2[0]) / 2, EX2[0] - EX1[1]
    ax.add_patch(Arc((mid, max(ys) * 1.05), width, max(ys) * 0.75, theta1=0, theta2=180, color=color,
                     linewidth=2 if c["spliced"] else 1.2, linestyle="-" if c["spliced"] else (0, (4, 3))))
    lab = (f"剪接连接 reads / spliced exon|exon: {c['spliced']:,}" if c["spliced"]
           else "剪接连接 reads / spliced exon|exon: 0")
    ax.text(mid, max(ys) * 1.6, lab, ha="center", va="center", fontsize=10, color=INK, weight="bold")
    ax.text(EX1[1] + 15, ly, f"← 外显子|内含子\n   exon|intron {c['ex_in']:,}", ha="left", va="bottom",
            fontsize=8.5, color=MUTED)
    ax.text(EX2[0] - 15, ly, f"内含子|外显子 →\nintron|exon {c['in_ex']:,}   ", ha="right", va="bottom",
            fontsize=8.5, color=MUTED)
    ax.set_ylabel("深度 / depth", color=MUTED, fontsize=9)
    ax.set_title(title, loc="left", fontsize=11, color=INK)
    style(ax)


def gene_model(ax):
    ax.set_xlim(0, 2900)
    ax.set_ylim(0, 4.9)
    ax.axis("off")
    y = 3.3
    ax.plot([INTRON[0], INTRON[1]], [y, y], color=MUTED, linewidth=1.5)
    for s, e, lab in ((*EX1, "外显子 1 / exon 1"), (*EX2, "外显子 2 / exon 2")):
        ax.add_patch(Rectangle((s, y - 0.28), e - s, 0.56, color=INK, linewidth=0))
        ax.text((s + e) / 2, y + 0.55, lab, ha="center", va="bottom", fontsize=9, color=INK)
    ax.text(sum(INTRON) / 2, y + 0.18, "II 型内含子 / group II intron (1,252 nt)", ha="center", fontsize=9,
            color=MUTED)
    ax.text(0, 4.9, "COX1, Travisia forbesii OZ479805 (DToL wuTraForb1)", fontsize=11, color=INK, va="top")
    # three read types
    rows = [(2.25, "剪接后的 read：外显子 1 末端直接接外显子 2\n→ 只有剪接过的 RNA 才有", "spliced"),
            (1.45, "未剪接的 read：跨外显子|内含子边界\n（前体 RNA 或 DNA）", "unspliced"),
            (0.65, "内含子内部的 read：前体 RNA 和 DNA 都有，不能说明剪接", "inside")]
    for yy, txt, kind in rows:
        if kind == "spliced":
            ax.plot([EX1[1] - 150, EX1[1]], [yy, yy], color=RNA, linewidth=4, solid_capstyle="round")
            ax.plot([EX2[0], EX2[0] + 150], [yy, yy], color=RNA, linewidth=4, solid_capstyle="round")
            ax.plot([EX1[1], EX2[0]], [yy, yy], color=RNA, linewidth=0.8, linestyle=(0, (2, 3)))
        elif kind == "unspliced":
            ax.plot([EX1[1] - 150, EX1[1] + 150], [yy, yy], color=MUTED, linewidth=4, solid_capstyle="round")
            ax.plot([EX2[0] - 150, EX2[0] + 150], [yy, yy], color=MUTED, linewidth=4, solid_capstyle="round")
        else:
            ax.plot([800, 1100], [yy, yy], color=MUTED, linewidth=4, solid_capstyle="round")
        ax.text(1760 if kind != "inside" else 1180, yy, txt, va="center", fontsize=8.5, color=INK, wrap=True)


def main():
    fig = plt.figure(figsize=(11, 8.6), dpi=150)
    gs = fig.add_gridspec(3, 1, height_ratios=[1.25, 1, 1], hspace=0.45)
    gene_model(fig.add_subplot(gs[0]))
    coverage(fig.add_subplot(gs[1]), "rna", RNA, "转录组 RNA-seq（ERR17229609，前 400 万对）")
    ax = fig.add_subplot(gs[2])
    coverage(ax, "dna", DNA, "同一个体的 DNA（10x Illumina ERR17229605）")
    ax.set_xlabel("线粒体基因组位置 / position in OZ479805 (nt)", color=MUTED, fontsize=9)
    fig.text(0.01, 0.005, "判断剪接看跨剪接点的 reads，不看覆盖度：RNA 有 5,799 条外显子|外显子 reads，DNA 为 0 → 内含子在 DNA 中，"
             "RNA 中被剪掉。Splicing = junction reads, not coverage.", fontsize=8.5, color=MUTED)
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"travisia_cox1_splicing.{ext}", bbox_inches="tight", facecolor="white")


if __name__ == "__main__":
    main()
