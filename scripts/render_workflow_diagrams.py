#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Render the bilingual workflow diagrams (SVG and Mermaid) of `mitom run` (mitom_run.RUN_STAGES) into docs/flowcharts/.
"""

from __future__ import annotations

import html
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "docs" / "flowcharts"

FONT_STACK = "Arial, Helvetica, sans-serif"


Node = Tuple[str, Sequence[str]]


def _text_lines(lines: Sequence[str], x: int, y: int, size: int = 16) -> str:
    tspans: List[str] = []
    for i, line in enumerate(lines):
        dy = "0" if i == 0 else "1.35em"
        tspans.append(
            f'<tspan x="{x}" dy="{dy}">{html.escape(line)}</tspan>'
        )
    return (
        f'<text x="{x}" y="{y}" font-size="{size}" '
        f'text-anchor="middle" dominant-baseline="middle">'
        + "".join(tspans)
        + "</text>"
    )


def render_svg(title: str, subtitle: str, nodes: Iterable[Node], out_path: Path) -> None:
    nodes = list(nodes)
    width = 980
    node_x = 140
    node_w = 700
    node_h = 72
    gap = 28
    top = 112
    height = top + len(nodes) * node_h + (len(nodes) - 1) * gap + 72
    center_x = node_x + node_w // 2

    parts: List[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-label="{html.escape(title)}">',
        "<defs>",
        "<style>",
        f"text {{ font-family: {FONT_STACK}; fill: #17202a; }}",
        ".title { font-weight: 700; }",
        ".subtitle { fill: #566573; }",
        ".node { fill: #f8fbff; stroke: #2f6f9f; stroke-width: 2; rx: 8; }",
        ".node-accent { fill: #eef7f1; stroke: #2f7d50; stroke-width: 2; rx: 8; }",
        ".arrow { stroke: #586069; stroke-width: 2.4; fill: none; marker-end: url(#arrowhead); }",
        "</style>",
        '<marker id="arrowhead" markerWidth="10" markerHeight="10" refX="8" refY="3" '
        'orient="auto" markerUnits="strokeWidth">',
        '<path d="M0,0 L0,6 L9,3 z" fill="#586069" />',
        "</marker>",
        "</defs>",
        '<rect x="0" y="0" width="980" height="100%" fill="#ffffff" />',
        f'<text class="title" x="{width // 2}" y="42" font-size="26" '
        f'text-anchor="middle">{html.escape(title)}</text>',
        f'<text class="subtitle" x="{width // 2}" y="72" font-size="15" '
        f'text-anchor="middle">{html.escape(subtitle)}</text>',
    ]

    for i, (_node_id, lines) in enumerate(nodes):
        y = top + i * (node_h + gap)
        css = "node-accent" if i in (0, len(nodes) - 1) else "node"
        parts.append(
            f'<rect class="{css}" x="{node_x}" y="{y}" width="{node_w}" height="{node_h}" />'
        )
        parts.append(_text_lines(lines, center_x, y + node_h // 2 - 8))
        if i < len(nodes) - 1:
            y1 = y + node_h
            y2 = y + node_h + gap - 6
            parts.append(f'<path class="arrow" d="M {center_x} {y1 + 5} L {center_x} {y2}" />')

    parts.append("</svg>")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(parts) + "\n", encoding="utf-8")


# The 11 stages of `mitom run` (mitom_run.RUN_STAGES; `mitom annotate` skips 3, 4, 9, 11). Keep in step with the
# step table of README.md section 3.
ZH_NODES: List[Node] = [
    ("input", ["1 读入与去重", "同一基因组只留一条；qc/dedup_report.tsv"]),
    ("annotate", ["2 第一轮注释", "13 个蛋白基因、2 个 rRNA、22 个 tRNA；核对环状 / 线状"]),
    ("reorient", ["3 统一方向与起点", "正链，从 cox1（--which_gene_first）开始"]),
    ("round2", ["4 转动后核查，必要时第二轮注释", "--if_reannotation_2 auto：只重注受转动影响的记录"]),
    ("write", ["5 写 GenBank", "写回原记录的物种、文献、标本信息"]),
    ("add_missing", ["6 补漏", "12S / 16S（rRNA HMM）、内含子；--rescue 时 MITOS2"]),
    ("curate", ["7 校正", "按 NCBI 规范修整端点和基因名；标出非编码区"]),
    ("ncbi", ["8 NCBI 检查与提交文件", "table2asn → submission/"]),
    ("compare", ["9 与原注释比较", "qc/original_vs_mitom.tsv（GenBank 输入）"]),
    ("maps", ["10 基因组圈图", "annotation/plots/"]),
    ("output", ["11 整理输出", "summary.tsv、annotation/、submission/、qc/、input/、logs/"]),
]
EN_NODES: List[Node] = [
    ("input", ["1 Read input, remove duplicates", "one record per genome; qc/dedup_report.tsv"]),
    ("annotate", ["2 Round-1 annotation", "13 protein genes, 2 rRNAs, 22 tRNAs; circular or linear checked"]),
    ("reorient", ["3 Re-orient", "plus strand, start at cox1 (--which_gene_first)"]),
    ("round2", ["4 Check after re-orientation (round 2)", "--if_reannotation_2 auto: only records the shift broke"]),
    ("write", ["5 Write GenBank", "original header (organism, references, voucher) copied back"]),
    ("add_missing", ["6 Add what the core missed", "12S / 16S (rRNA HMM), introns; MITOS2 with --rescue"]),
    ("curate", ["7 Curate", "gene ends and names to NCBI rules; non-coding regions"]),
    ("ncbi", ["8 NCBI check and submission files", "table2asn -> submission/"]),
    ("compare", ["9 Compare with the original annotation", "qc/original_vs_mitom.tsv (GenBank input)"]),
    ("maps", ["10 Genome maps", "annotation/plots/"]),
    ("output", ["11 Organise the output", "summary.tsv, annotation/, submission/, qc/, input/, logs/"]),
]


def write_mermaid(nodes: Sequence[Node], out_path: Path) -> None:
    ids = [chr(ord("A") + i) for i in range(len(nodes))]
    lines = ["flowchart TD"] + [f'    {i}["{t}<br/>{d}"]' for i, (_, (t, d)) in zip(ids, nodes)]
    lines += ["", "    " + " --> ".join(ids)]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    render_svg("MitoM 工作流程（mitom run）", "11 步，与终端的 [Step k/11] 一致", ZH_NODES,
               OUT_DIR / "mitom_workflow.zh.svg")
    render_svg("MitoM workflow (mitom run)", "11 steps, as [Step k/11] on the terminal", EN_NODES,
               OUT_DIR / "mitom_workflow.en.svg")
    write_mermaid(ZH_NODES, OUT_DIR / "mitom_workflow.zh.mmd")
    write_mermaid(EN_NODES, OUT_DIR / "mitom_workflow.en.mmd")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
