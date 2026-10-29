"""mito_common.py - small helpers shared by the MitoM modules: GenBank file suffixes, safe file names, report tables,
file copies."""
from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import List


def safe_name(s: str) -> str:
    s = str(s).strip().replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", s)
    return s or "NA"


# ONE list of GenBank file suffixes for auto_clade and for input detection (.gbff included)
GB_SUFFIXES = (".gb", ".gbk", ".gbf", ".genbank", ".gbff")


# reports written only when there is something in them
def write_nonempty_list(out_path: Path, items: List[str], comment: str = "") -> None:
    unique = sorted(set(items))
    if not unique:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    if comment:
        lines.append(comment + "\n")
    for item in unique:
        lines.append(item + "\n")
    out_path.write_text("".join(lines), encoding="utf-8")


def write_nonempty_table(out_path: Path, header: str, rows: List[str]) -> None:
    if not rows:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write(r + "\n")


def safe_copy(src: Path, dst: Path) -> None:
    """Copy *src* -> *dst*; src and dst being the same file is not an error (e.g. --out_root = the input folder)."""
    try:
        shutil.copy2(str(src), str(dst))
    except shutil.SameFileError:
        pass
