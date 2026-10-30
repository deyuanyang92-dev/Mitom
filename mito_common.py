"""mito_common.py - helpers shared by the MitoM modules: GenBank files (suffixes, finding them, atomic writes), FASTA
reading, safe file names, report tables, file copies, CDS stop qualifiers and codon translation, the annotation
core's gene record (CoreFeature), running external programs (run_tool)."""
from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import shlex
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Tuple

from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

_UMASK = os.umask(0)                 # the process umask, read once (os.umask can only be read by setting it)
os.umask(_UMASK)


def old_spelling(parser, name: str, dest: str, **kw) -> None:
    """An option's old spelling: still accepted (scripts written for earlier versions keep working), not shown in -h."""
    parser.add_argument(name, dest=dest, default=argparse.SUPPRESS, help=argparse.SUPPRESS, **kw)


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


def write_tsv(path, columns, rows) -> None:
    """A report table: tab-separated, LF line ends, the header line always (also without rows). Rows are dicts (a
    missing key is written empty, a key not in columns is left out) or sequences in column order; fields with a tab,
    a newline or a double quote are quoted as the csv module quotes them."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        for r in rows:
            w.writerow([r.get(c, "") for c in columns] if isinstance(r, dict) else r)


def write_table(out_path: Path, header: str, rows: List[str]) -> None:
    """A report table from tab-joined lines: the header line always, LF line ends."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(x + "\n" for x in [header, *rows]), encoding="utf-8")


REPORTS = "reports"                  # the folder of a run's reports (run, import, annotate)


def reports_dir(out_root: Path) -> Path:
    """<out_root>/reports, created: every step of a run writes its report there."""
    d = Path(out_root) / REPORTS
    d.mkdir(parents=True, exist_ok=True)
    return d


# problem lists: written only when there is something in them
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


# ── the annotation core's genes (mito_core.py; its tRNA step in mito_trna.py) ─────────────────────────────────────
@dataclass
class CoreFeature:
    kind: str              # CDS / rRNA / tRNA
    name: str              # MITOS2-style: cox1, rrnL, trnL1
    start: int             # 0-based half-open
    end: int
    strand: int
    note: str = ""
    partial5: bool = False
    anticodon: str = ""
    acpos: int = -1            # tRNA: 0-based position of the anticodon within the tRNA (MiTFi)
    model: str = ""            # tRNA: the covariance model that found it
    score: float = 0.0


def positions(f: CoreFeature, n: int) -> set:
    return {i % n for i in range(f.start, f.end)}


# ── locations on a circular sequence: rotation, reverse complement ─────────────────────────────────────────────────
def rc_location(loc, seq_len: int):
    """loc on the reverse complement of a sequence of length seq_len; partial ends stay partial and change sides
    (a 5'-partial plus-strand start <s becomes a 5'-partial minus-strand end >)."""
    if isinstance(loc, CompoundLocation):
        return CompoundLocation(
            [rc_location(p, seq_len) for p in loc.parts][::-1], operator=loc.operator)
    start, end = int(loc.start), int(loc.end)
    new_strand  = None if loc.strand is None else (-loc.strand)
    new_start = BeforePosition(seq_len - end) if isinstance(loc.end, AfterPosition) else seq_len - end
    new_end = AfterPosition(seq_len - start) if isinstance(loc.start, BeforePosition) else seq_len - start
    return FeatureLocation(new_start, new_end, strand=new_strand)


# qualifiers that name positions of the sequence: (pos:<location>,aa:...[,seq:...])
POSITION_QUALIFIERS = ("transl_except", "anticodon")


def _split_pos(text: str):
    """'(pos:LOC,rest' -> ('(pos:', LOC, ',rest') with LOC up to the comma outside its parentheses; None if no pos:."""
    i = text.find("pos:")
    if i < 0:
        return None
    depth, j = 0, i + 4
    while j < len(text):
        ch = text[j]
        if ch == "(":
            depth += 1
        elif ch == ")":
            if depth == 0:
                break
            depth -= 1
        elif ch == "," and depth == 0:
            break
        j += 1
    return text[:i + 4], text[i + 4:j], text[j:]


def move_positions(text: str, seq_len: int, shift: int = 0, rc: bool = False) -> str:
    """A /transl_except or /anticodon value with its positions moved as the sequence is: reverse-complemented first
    (rc), then rotated so that 0-based position `shift` becomes the start (as shift_location). Ranges cut by the new
    origin become join(a..L,1..b); ranges that meet across the old origin are merged. Unparsable values unchanged."""
    parts = _split_pos(text)
    if parts is None:
        return text
    head, loc, rest = parts
    comp = loc.startswith("complement(") and loc.endswith(")")
    inner = loc[len("complement("):-1] if comp else loc
    if inner.startswith("join(") and inner.endswith(")"):
        inner = inner[len("join("):-1]
    try:
        ranges = [tuple(int(x) for x in r.split("..")) if ".." in r else (int(r), int(r)) for r in inner.split(",")]
    except ValueError:
        return text
    if rc:
        ranges = [(seq_len - b + 1, seq_len - a + 1) for a, b in reversed(ranges)]
        comp = not comp
    out: List[Tuple[int, int]] = []
    for a, b in ranges:
        a2, b2 = (a - 1 - shift) % seq_len + 1, (b - 1 - shift) % seq_len + 1
        pieces = [(a2, b2)] if a2 <= b2 else [(a2, seq_len), (1, b2)]
        for p in pieces:
            if out and out[-1][1] + 1 == p[0]:
                out[-1] = (out[-1][0], p[1])
            else:
                out.append(p)
    body = ",".join(f"{a}..{b}" if a != b else f"{a}" for a, b in out)
    body = f"join({body})" if len(out) > 1 else body
    return head + (f"complement({body})" if comp else body) + rest


def move_position_qualifiers(qualifiers: Dict[str, List[str]], seq_len: int, shift: int = 0,
                             rc: bool = False) -> Dict[str, List[str]]:
    """A copy of a feature's qualifiers whose position qualifiers (POSITION_QUALIFIERS) follow the sequence
    (move_positions); the other qualifiers unchanged."""
    q = dict(qualifiers)
    for k in POSITION_QUALIFIERS:
        if k in q:
            q[k] = [move_positions(str(v), seq_len, shift, rc) for v in q[k]]
    return q


def reverse_complement_record(record: "SeqRecord") -> "SeqRecord":
    seq_len = len(record.seq)
    new_rec = SeqRecord(record.seq.reverse_complement(), id=record.id,
                        name=record.name, description=record.description)
    new_rec.annotations = dict(record.annotations)
    new_rec.dbxrefs     = list(getattr(record, "dbxrefs", []))
    new_features: List["SeqFeature"] = []
    for feat in record.features:
        try:
            new_loc = rc_location(feat.location, seq_len)
            new_features.append(SeqFeature(location=new_loc, type=feat.type,
                                            id=getattr(feat, "id", "<unknown>"),
                                            qualifiers=move_position_qualifiers(feat.qualifiers, seq_len, rc=True)))
        except Exception:
            new_features.append(feat)
    new_rec.features = new_features
    return new_rec


def merge_abutting_parts(parts: List, operator: str = "join"):
    """Consecutive parts of the same strand that abut (a join across the old origin, [s, n) + [0, e), after the
    shift) become one interval; parts are in transcript order, so on the minus strand the later part lies
    before the earlier one. The outer ends keep their (fuzzy) positions. One part left: a plain FeatureLocation."""
    merged = [parts[0]]
    for part in parts[1:]:
        prev = merged[-1]
        if part.strand == prev.strand and part.strand == -1 and int(part.end) == int(prev.start):
            merged[-1] = FeatureLocation(part.start, prev.end, strand=prev.strand)
        elif part.strand == prev.strand and part.strand != -1 and int(prev.end) == int(part.start):
            merged[-1] = FeatureLocation(prev.start, part.end, strand=prev.strand)
        else:
            merged.append(part)
    return merged[0] if len(merged) == 1 else CompoundLocation(merged, operator=operator)


def shift_location(location, shift: int, seq_len: int):
    """Shift all feature coordinates by -shift (mod seq_len). Nested joins are flattened, parts that meet at the old
    origin merged; a whole-sequence feature (source) stays whole; a zero-length point stays a point; partial ends
    (< >) stay partial (trans-spliced placozoan cox1 pieces, 2026-10-29)."""
    if isinstance(location, CompoundLocation):
        flat_parts: List = []
        for p in location.parts:
            shifted = shift_location(p, shift, seq_len)
            if isinstance(shifted, CompoundLocation):
                flat_parts.extend(shifted.parts)
            else:
                flat_parts.append(shifted)
        return merge_abutting_parts(flat_parts, location.operator)

    orig_start = int(location.start)
    orig_end   = int(location.end)
    strand     = location.strand
    length     = orig_end - orig_start

    if length == seq_len:
        return FeatureLocation(0, seq_len, strand=strand)

    start = (orig_start - shift) % seq_len
    end   = (orig_end   - shift) % seq_len

    if end == 0 and length > 0:
        end = seq_len
    start, end = (type(old)(new) if isinstance(old, (BeforePosition, AfterPosition)) else new
                  for old, new in ((location.start, start), (location.end, end)))

    if start <= end:
        return FeatureLocation(start, end, strand=strand)

    return CompoundLocation([FeatureLocation(start, seq_len, strand=strand),
                              FeatureLocation(0, end, strand=strand)])


def feature_start_pos(feat: "SeqFeature") -> int:
    loc = feat.location
    if isinstance(loc, CompoundLocation) and loc.parts:
        return int(loc.parts[0].start)
    return int(loc.start)


POLYA_NOTE = "TAA stop codon is completed by the addition of 3' A residues to the mRNA"


def transl_except_term(start: int, end: int, strand: int, tail: int) -> str:
    """/transl_except of an abbreviated stop (T / TA: `tail` nt, completed by polyadenylation) at the 3' end of a CDS
    whose 3'-most part is start..end (1-based) on `strand`."""
    if strand == -1:
        pos = f"complement({start}..{start + tail - 1})" if tail > 1 else f"complement({start})"
    else:
        pos = f"{end - tail + 1}..{end}" if tail > 1 else f"{end}"
    return f"(pos:{pos},aa:TERM)"


def abbreviated_stop(nt: str) -> bool:
    """The CDS sequence nt ends with an incomplete stop codon, T or TA after its last whole codon (mitochondrial
    transcripts complete it by polyadenylation)."""
    tail = len(nt) % 3
    return bool(tail) and nt[-tail:] in ("T", "TA")


def ends_with_stop(nt: str, stops) -> bool:
    """The CDS sequence nt ends with a stop codon of `stops` in frame, or with an abbreviated stop (T / TA)."""
    return abbreviated_stop(nt) if len(nt) % 3 else nt[-3:] in stops


def stop_qualifiers(rec, loc, code: int) -> Dict[str, List[str]]:
    """Incomplete stop (T / TA completed by polyadenylation): /transl_except and the standard note."""
    nt = str(loc.extract(rec.seq)).upper()
    tail = len(nt) % 3
    if not abbreviated_stop(nt):
        return {}
    last = loc.parts[-1]
    return {"transl_except": [transl_except_term(int(last.start) + 1, int(last.end), last.strand, tail)],
            "note": [POLYA_NOTE]}


def codon_translate(nt, code: int) -> str:
    """Translation of the whole codons of nt (a str or Seq) with genetic code `code`."""
    return str(Seq(str(nt)[: len(nt) - len(nt) % 3]).translate(table=code))


# ── external programs ─────────────────────────────────────────────────────────────────────────────────────────────
TOOL_LOG = logging.getLogger("mitom.tools")
STDERR_TAIL = 20                     # lines of a failed program's stderr in the log ...
DEBUG_STDERR = 200                   # ... of any program's stderr in debug mode
STDERR_MARK = "    | "               # each under the [CMD] line, after this (a log line of the program itself stays
                                     # inside the [CMD] entry for `mitom debug`)


def debug_mode() -> bool:
    """`mitom ... --debug` (MITOM_DEBUG=1, set for every program MitoM starts): log what each program did."""
    return bool(os.environ.get("MITOM_DEBUG"))


if debug_mode() and not TOOL_LOG.handlers:      # a step run by `mitom ... --debug`: the programs it runs, on its stderr
    _h = logging.StreamHandler()                 # (the workflow keeps that in the step's log; its own run log gets
    _h.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))   # them through its root logger:
    TOOL_LOG.addHandler(_h)                      # mitom_log.setup_root_logger removes this handler)
    TOOL_LOG.setLevel(logging.DEBUG)
    TOOL_LOG.propagate = False


def command_line(cmd) -> str:
    """The command as one shell-quoted line (logs, reproduction)."""
    return " ".join(shlex.quote(str(x)) for x in cmd)


def run_tool(cmd, *, cwd=None, env=None, timeout=None, check: bool = False, capture: bool = True,
             stdout=None, stderr=None) -> subprocess.CompletedProcess:
    """Run an external program: the one way the MitoM modules run programs. Output is captured as text unless
    capture=False (it then goes where this process's output goes) or stdout / stderr are given. Logged at DEBUG level
    on the 'mitom.tools' logger: '[CMD rc=<exit> <seconds>s] <command line>' (and the folder), then the last
    STDERR_TAIL lines of stderr when the program failed (DEBUG_STDERR lines of any program in debug mode), each line
    after STDERR_MARK. Exceptions
    as subprocess.run (OSError: program not found, TimeoutExpired; CalledProcessError when check=True), each logged
    before it is raised."""
    args = [str(x) for x in cmd]
    kw = {"cwd": str(cwd) if cwd else None, "env": env, "timeout": timeout, "text": True, "errors": "replace"}
    if stdout is not None or stderr is not None:
        kw.update(stdout=stdout, stderr=stderr)
    elif capture:
        kw.update(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    where = f"  (in {cwd})" if cwd else ""
    t0 = time.time()
    try:
        p = subprocess.run(args, **kw)
    except (OSError, subprocess.SubprocessError) as exc:
        TOOL_LOG.debug(f"[CMD failed {time.time() - t0:.1f}s] {command_line(args)}{where}: {exc!r}")
        raise
    err = p.stderr.strip().splitlines() if (p.returncode or debug_mode()) and isinstance(p.stderr, str) else []
    TOOL_LOG.debug(f"[CMD rc={p.returncode} {time.time() - t0:.1f}s] {command_line(args)}{where}"
                   + "".join(f"\n{STDERR_MARK}{ln}" for ln in err[-(DEBUG_STDERR if debug_mode() else STDERR_TAIL):]))
    if check and p.returncode:
        raise subprocess.CalledProcessError(p.returncode, args, p.stdout, p.stderr)
    return p


def run_or_raise(cmd, **kw) -> subprocess.CompletedProcess:
    """run_tool for the annotation core's programs: RuntimeError('<program> failed (<exit code>): <end of stderr>')
    when the program fails."""
    p = run_tool(cmd, **kw)
    if p.returncode != 0:
        raise RuntimeError(f"{Path(cmd[0]).name} failed ({p.returncode}): {p.stderr.strip()[-400:]}")
    return p


def safe_copy(src: Path, dst: Path) -> None:
    """Copy *src* -> *dst*; src and dst being the same file is not an error (e.g. --out_root = the input folder)."""
    try:
        shutil.copy2(str(src), str(dst))
    except shutil.SameFileError:
        pass
