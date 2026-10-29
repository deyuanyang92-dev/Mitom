"""mito_common.py - helpers shared by the MitoM modules: GenBank files (suffixes, finding them, atomic writes), FASTA
reading, safe file names, report tables, file copies, CDS stop qualifiers and codon translation."""
from __future__ import annotations

import os
import re
import shutil
import tempfile
from pathlib import Path
from typing import Dict, Iterator, List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

_UMASK = os.umask(0)                 # the process umask, read once (os.umask can only be read by setting it)
os.umask(_UMASK)


def safe_name(s: str) -> str:
    s = str(s).strip().replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", s)
    return s or "NA"


# ONE list of GenBank file suffixes for every MitoM module (.gbff included)
GB_SUFFIXES = (".gb", ".gbk", ".gbf", ".genbank", ".gbff")


def iter_gb_files(path) -> List[Path]:
    """GenBank files of an input: the file itself, or every file with a GenBank suffix (GB_SUFFIXES, any case) below a
    folder, sorted; [] for a path that does not exist."""
    p = Path(path)
    if p.is_file():
        return [p]
    return sorted(x for x in p.rglob("*") if x.is_file() and x.suffix.lower() in GB_SUFFIXES)


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


def write_gb(recs, path) -> None:
    """Write GenBank records atomically: a temporary file in the same folder replaces `path` only when complete, so a
    failure (or a killed run) never leaves a half-written file; permissions as the file had (a new file: umask)."""
    from Bio import SeqIO
    path = Path(path)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as h:
            SeqIO.write(recs, h, "genbank")
        if path.exists():
            shutil.copymode(path, tmp)
        else:
            os.chmod(tmp, 0o666 & ~_UMASK)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def read_fasta(path) -> Iterator[SeqRecord]:
    """FASTA records (id = first word, description = whole header), streamed, without Bio.SeqIO (it imports numpy,
    ~27 MB, resident while the annotation programs run)."""
    head, parts = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if head is not None:
                    yield SeqRecord(Seq("".join(parts)), id=(head.split() or [""])[0], description=head)
                head, parts = line[1:].strip(), []
            elif head is not None:
                parts.append("".join(line.split()))
    if head is not None:
        yield SeqRecord(Seq("".join(parts)), id=(head.split() or [""])[0], description=head)


def topology_from_header(description: str) -> str:
    """linear when the FASTA header says topology=linear, else circular (the annotation core's rule)."""
    return "linear" if re.search(r"topology\s*=\s*linear", description or "", re.I) else "circular"


def contig_record_id(prefix: str, seq_id: str, k: int, n_seqs: int) -> str:
    """Record name of sequence k of a sample as the core gives it: the prefix for one sequence, else the sequence's
    own (cleaned, <= 16 characters) name."""
    return prefix if n_seqs == 1 else (re.sub(r"[^A-Za-z0-9_.-]", "_", seq_id)[:16] or f"{prefix}_{k + 1}")


POLYA_NOTE = "TAA stop codon is completed by the addition of 3' A residues to the mRNA"


def transl_except_term(start: int, end: int, strand: int, tail: int) -> str:
    """/transl_except of an abbreviated stop (T / TA: `tail` nt, completed by polyadenylation) at the 3' end of a CDS
    whose 3'-most part is start..end (1-based) on `strand`."""
    if strand == -1:
        pos = f"complement({start}..{start + tail - 1})" if tail > 1 else f"complement({start})"
    else:
        pos = f"{end - tail + 1}..{end}" if tail > 1 else f"{end}"
    return f"(pos:{pos},aa:TERM)"


def stop_qualifiers(rec, loc, code: int) -> Dict[str, List[str]]:
    """Incomplete stop (T / TA completed by polyadenylation): /transl_except and the standard note."""
    nt = str(loc.extract(rec.seq)).upper()
    tail = len(nt) % 3
    if tail == 0 or nt[-tail:] not in ("T", "TA"):
        return {}
    last = loc.parts[-1]
    return {"transl_except": [transl_except_term(int(last.start) + 1, int(last.end), last.strand, tail)],
            "note": [POLYA_NOTE]}


def codon_translate(nt, code: int) -> str:
    """Translation of the whole codons of nt (a str or Seq) with genetic code `code`."""
    return str(Seq(str(nt)[: len(nt) - len(nt) % 3]).translate(table=code))


def safe_copy(src: Path, dst: Path) -> None:
    """Copy *src* -> *dst*; src and dst being the same file is not an error (e.g. --out_root = the input folder)."""
    try:
        shutil.copy2(str(src), str(dst))
    except shutil.SameFileError:
        pass
