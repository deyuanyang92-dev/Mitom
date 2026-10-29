#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""mitom_pipeline.py - the MitoM workflow: `mitom run`, `import`, `annotate`, `reorient`, `transfer_metadata`,
`replace_organism` and `dedup` (`mitom <command> -h` lists the options; batch_mitoz.py is the old name).

run      GenBank input: records deduplicated (RefSeq NC_ kept), the original annotation saved for the comparison;
         FASTA input: one sample per record -> round-1 annotation (MitoM's core mito_core.py, or --engine mitoz /
         mitos2 via mito_engines.py) -> re-orientation (plus strand, cox1 first) -> round 2 only for the records
         whose annotation the shift breaks (--if_reannotation_2 auto) -> final GenBank with the original header and
         /organism -> MITOS2 rescue (optional), rRNA HMM fallback, introns (mito_intron.py de novo, then
         mito_intron_ref.py site transfer) -> boundary rules and NCBI names (mito_curate.py) -> NCBI submission files
         (ncbi_validate.py, table2asn) -> original vs MitoM (mito_compare.py) -> genome maps (mito_plot.py) -> the
         output layout of OUTPUT_LAYOUT (README.txt in the output folder).
import   the same steps on existing MitoZ / MITOS2 results or any annotated GenBank, without annotating again.
annotate one round of annotation and the steps after it, in the work layout (04.final_gb/, 05.curated/ ...).

The annotation runs --max_tasks samples at a time with --threads each, in threads or with Snakemake
(--scheduler_backend snakemake / auto; Snakemake is optional and not installed with MitoM); everything after the
annotation runs in this process.

Biopython: Cock et al. (2009) Bioinformatics doi:10.1093/bioinformatics/btp163
Snakemake: Mölder et al. (2021) F1000Res doi:10.12688/f1000research.29032.2
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import difflib
import functools
import json
import logging
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, TypeVar

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
except Exception:
    SeqIO = None  # type: ignore

# ─────────────────────────────────────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────────────────────────────────────
_LOG_LOCK = threading.Lock()


class _LockedHandler(logging.Handler):
    def __init__(self, inner: logging.Handler):
        super().__init__()
        self._inner = inner
        self.setFormatter(inner.formatter)
        self.setLevel(inner.level)

    def emit(self, record: logging.LogRecord) -> None:
        with _LOG_LOCK:
            self._inner.emit(record)


def setup_root_logger(log_file: Path, console_level: int = logging.INFO, mode: str = "a") -> None:
    """Everything goes to log_file (mode "w": a new log for this run); the console gets console_level and above
    (run: WARNING, its steps are printed by Steps / say())."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    for h in list(root.handlers):
        root.removeHandler(h)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh = logging.FileHandler(str(log_file), mode=mode, encoding="utf-8")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter("    ! %(message)s") if console_level > logging.INFO else fmt)
    sh.setLevel(console_level)
    root.addHandler(_LockedHandler(fh))
    root.addHandler(_LockedHandler(sh))


# ── what the user sees on the terminal during `run` (GetOrganelle style: numbered steps, counts, time, files) ──
def say(text: str = "") -> None:
    """One line on the terminal and in the run log."""
    print(text, flush=True)
    if text.strip():
        logging.info(text.strip())


def duration(sec: float) -> str:
    sec = int(round(sec))
    if sec < 60:
        return f"{sec} s"
    if sec < 3600:
        return f"{sec // 60} min {sec % 60:02d} s"
    return f"{sec // 3600} h {sec % 3600 // 60:02d} min"


class Steps:
    """Numbered steps: start() prints '[Step n/N] title', info() an indented line, done() the result and the time."""

    def __init__(self, total: int):
        self.total, self.n, self.t0, self.run_t0 = total, 0, time.time(), time.time()

    def start(self, title: str) -> None:
        self.n += 1
        self.t0 = time.time()
        say("")
        say(f"[Step {self.n}/{self.total}] {title}")

    def info(self, text: str) -> None:
        say(f"    {text}")

    def done(self, text: str = "") -> None:
        say(f"    -> {text}  ({duration(time.time() - self.t0)})" if text else f"    -> done ({duration(time.time() - self.t0)})")


class Progress:
    """'    label: 120/636 (19%), 5 min 02 s so far, about 18 min left' - at most every `every` seconds, and at the end."""

    def __init__(self, label: str, total: int, every: float = 30.0):
        self.label, self.total, self.every = label, max(1, total), every
        self.t0 = self.last = time.time()
        self.shown = -1
        self.lock = threading.Lock()

    def update(self, done: int) -> None:
        with self.lock:
            now = time.time()
            if done == self.shown or (done < self.total and now - self.last < self.every):
                return
            self.last, self.shown = now, done
        el = now - self.t0
        left = f", about {duration(el / done * (self.total - done))} left" if 0 < done < self.total else ""
        say(f"    {self.label}: {done}/{self.total} ({100 * done // self.total}%), {duration(el)} so far{left}")


_STAGE_PROGRESS: Optional[Progress] = None          # set by run_annotation_stage; updated per finished sample
STAGE_LABEL = {"01.anno_round1": "round-1 annotation", "03.anno_round2_reoriented": "round-2 annotation"}


def watch_files(pattern_dir: Path, pattern: str, progress: Progress, stop: threading.Event) -> threading.Thread:
    """Background thread: progress.update(number of files matching pattern) every 10 s until stop is set."""
    def loop() -> None:
        while not stop.wait(10):
            progress.update(len(list(pattern_dir.glob(pattern))))
    t = threading.Thread(target=loop, daemon=True)
    t.start()
    return t


def get_sample_logger(name: str, log_file: Path) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    lg = logging.getLogger(name)
    lg.setLevel(logging.INFO)
    lg.propagate = False
    for h in list(lg.handlers):
        try:
            h.close()
        except Exception:
            pass
        lg.removeHandler(h)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh = logging.FileHandler(str(log_file), encoding="utf-8")
    fh.setFormatter(fmt)
    lg.addHandler(fh)
    return lg


def close_sample_logger(name: str) -> None:
    lg = logging.getLogger(name)
    for h in list(lg.handlers):
        try:
            h.close()
        except Exception:
            pass
        lg.removeHandler(h)


# ─────────────────────────────────────────────────────────────────────────────
# Generic helpers
# ─────────────────────────────────────────────────────────────────────────────
def tail_text(s: str, n_lines: int = 120) -> str:
    return "\n".join((s or "").splitlines()[-n_lines:])


def safe_name(s: str) -> str:
    s = str(s).strip().replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", s)
    return s or "NA"


def header_first_token(desc_or_id: str) -> str:
    s = (desc_or_id or "").strip()
    if s.startswith(">"):
        s = s[1:]
    return s.split()[0] if s else "NA"


def which_or_none(exe: str) -> Optional[str]:
    try:
        return shutil.which(exe)
    except Exception:
        return None



MITO_CORE = Path(__file__).resolve().parent / "mito_core.py"      # MitoM's annotation core
MITO_ENGINES = Path(__file__).resolve().parent / "mito_engines.py"  # MitoZ / MITOS2 as engines, imported results
ENGINE_LABEL = {"native": "MitoM core", "mitoz": "MitoZ", "mitos2": "MITOS2", "input": "annotation of the input"}
# --clade choices: taxonomic groups that set the genetic code and the curation rules of the group (names kept
# from earlier MitoM versions so that existing commands keep working)
BASE_CLADES = ("Chordata", "Arthropoda", "Echinodermata", "Annelida-segmented-worms", "Bryozoa", "Mollusca",
                "Nematoda", "Nemertea-ribbon-worms", "Porifera-sponges")
# genetic code per clade when --genetic_code is not given
CLADE_CODES = {"Chordata": 2, "Echinodermata": 9, "Porifera-sponges": 4}
# further groups that need --genetic_code; "Animal" = any other phylum
EXTRA_CLADES = ("Cnidaria", "Platyhelminthes-flatworms", "Chaetognatha-arrow-worms", "Animal")
ALL_CLADES = BASE_CLADES + EXTRA_CLADES
# --clade auto: phylum (from the GenBank taxonomy) -> clade; any other metazoan phylum -> "Animal"
PHYLUM_CLADE = {"Chordata": "Chordata", "Arthropoda": "Arthropoda", "Echinodermata": "Echinodermata",
                "Annelida": "Annelida-segmented-worms", "Bryozoa": "Bryozoa", "Mollusca": "Mollusca",
                "Nematoda": "Nematoda", "Nemertea": "Nemertea-ribbon-worms", "Porifera": "Porifera-sponges",
                "Cnidaria": "Cnidaria", "Platyhelminthes": "Platyhelminthes-flatworms",
                "Chaetognatha": "Chaetognatha-arrow-worms"}
CLADE_PHYLUM = {clade: phylum for phylum, clade in PHYLUM_CLADE.items()}      # inverse: explicit --clade -> phylum
METAZOAN_PHYLA = set(PHYLUM_CLADE) | {
    "Hemichordata", "Xenacoelomorpha", "Ctenophora", "Placozoa", "Tardigrada", "Onychophora", "Nematomorpha",
    "Priapulida", "Kinorhyncha", "Loricifera", "Brachiopoda", "Phoronida", "Entoprocta", "Cycliophora", "Rotifera",
    "Acanthocephala", "Gastrotricha", "Gnathostomulida", "Orthonectida", "Dicyemida", "Micrognathozoa"}
# usual NCBI mitochondrial code when the input has no CDS transl_table (benchmark/mito_sample.py DEFAULT_CODE)
PHYLUM_CODE = {"Chordata": 2, "Echinodermata": 9, "Hemichordata": 9, "Platyhelminthes": 9, "Cnidaria": 4,
               "Porifera": 4, "Placozoa": 4, "Ctenophora": 4}
# class / subphylum names that stand for a phylum when the lineage has no phylum token (Geneious, MitoFish and
# hand-edited headers: "Vertebrata", "Hexapoda", "Gastropoda"...); used only when METAZOAN_PHYLA matches nothing
PHYLUM_SYNONYMS = {
    **dict.fromkeys(("Vertebrata", "Craniata", "Tunicata", "Cephalochordata"), "Chordata"),
    **dict.fromkeys(("Hexapoda", "Insecta", "Crustacea", "Chelicerata", "Myriapoda"), "Arthropoda"),
    **dict.fromkeys(("Gastropoda", "Bivalvia", "Cephalopoda", "Polyplacophora", "Scaphopoda"), "Mollusca"),
    **dict.fromkeys(("Polychaeta", "Clitellata", "Sipuncula"), "Annelida"),
    **dict.fromkeys(("Anthozoa", "Hydrozoa", "Scyphozoa", "Cubozoa"), "Cnidaria"),
    **dict.fromkeys(("Demospongiae", "Hexactinellida", "Calcarea", "Homoscleromorpha"), "Porifera"),
    **dict.fromkeys(("Turbellaria", "Trematoda", "Cestoda", "Monogenea"), "Platyhelminthes"),
    **dict.fromkeys(("Chromadorea", "Enoplea"), "Nematoda"),
}
# ONE list of GenBank file suffixes for auto_clade and for input detection (.gbff included)
GB_SUFFIXES = (".gb", ".gbk", ".gbf", ".genbank", ".gbff")


def auto_clade(in_path: Path, genetic_code: Optional[int] = None) -> Tuple[Optional[str], Optional[int], Optional[str], str]:
    """--clade auto: (clade, genetic code, phylum, '') from the taxonomy and CDS transl_table of the GenBank
    input; (None, None, None, error) when it cannot be decided: several phyla, a taxonomy that names no known
    phylum, CDS transl_tables that differ between records (unless `genetic_code` is given), unreadable file.
    Class names in PHYLUM_SYNONYMS stand for their phylum when the lineage has no phylum token. Input with no
    taxonomy at all (FASTA, or GenBank without a lineage) keeps the previous default clade with a WARNING
    message, as before auto existed (2026-10-21 check); records without taxonomy do not block the
    known phylum of the others. A non-empty taxonomy that matches nothing is an ERROR, never a fallback."""
    in_path = Path(in_path)
    files = [in_path] if in_path.is_file() else sorted(p for p in in_path.rglob("*") if p.is_file())
    files = [p for p in files if p.suffix.lower() in GB_SUFFIXES]
    need = "give --clade (and --genetic_code), e.g. --clade Animal --genetic_code 5"
    phyla: Dict[str, int] = {}
    codes: Dict[int, int] = {}              # genetic code -> number of RECORDS (not CDS)
    unmatched: Dict[str, int] = {}          # lineage of records whose taxonomy names no known phylum
    for f in files:
        try:
            for rec in SeqIO.parse(str(f), "genbank"):
                tax = rec.annotations.get("taxonomy")
                tax = [str(t) for t in tax] if isinstance(tax, list) else []
                p = next((t for t in tax if t in METAZOAN_PHYLA), "") \
                    or next((PHYLUM_SYNONYMS[t] for t in tax if t in PHYLUM_SYNONYMS), "")
                if p:
                    phyla[p] = phyla.get(p, 0) + 1
                elif tax:
                    lineage = "; ".join(tax[-3:])
                    unmatched[lineage] = unmatched.get(lineage, 0) + 1
                rec_codes = [str(ft.qualifiers["transl_table"][0]) for ft in rec.features
                             if ft.type == "CDS" and "transl_table" in ft.qualifiers]
                if rec_codes:
                    try:
                        c = int(max(dict.fromkeys(rec_codes), key=rec_codes.count))
                    except ValueError:
                        continue
                    codes[c] = codes.get(c, 0) + 1
        except Exception as exc:                          # truncated / not a GenBank file: name it
            return None, None, None, f"--clade auto: cannot read {f}: {exc}; fix the file or {need}"
    if unmatched:
        got = "; ".join(f"'{lin}' ({n})" for lin, n in sorted(unmatched.items())[:5])
        return None, None, None, (f"--clade auto: the taxonomy of {sum(unmatched.values())} record(s) names no known "
                                  f"metazoan phylum (lineage ends {got}); {need}")
    if len(phyla) > 1:
        got = ", ".join(f"{p} ({n})" for p, n in sorted(phyla.items()))
        return None, None, None, (f"--clade auto: the input must be one metazoan phylum, found {got}; "
                                  f"run each phylum separately or {need}")
    if genetic_code is None and len(codes) > 1:
        got = ", ".join(f"{c}: {n} records" for c, n in sorted(codes.items()))
        return None, None, None, (f"--clade auto: the CDS transl_table differs between records ({got}); "
                                  f"run the groups separately or give --genetic_code")
    if not phyla:
        fallback = "Annelida-segmented-worms"
        code = genetic_code or (next(iter(codes)) if codes else 5)
        return fallback, code, None, (
            f"WARNING: --clade auto found no taxonomy in the input (FASTA or GenBank without a lineage); "
            f"using the previous default --clade {fallback} --genetic_code {code}. "
            f"For other taxa {need}")
    phylum = next(iter(phyla))
    code = genetic_code or (next(iter(codes)) if codes else PHYLUM_CODE.get(phylum, 5))
    return PHYLUM_CLADE.get(phylum, "Animal"), code, phylum, ""



_FASTA_EXTS = {".fa", ".fasta", ".fna", ".fas"}
_GB_EXTS    = GB_SUFFIXES


def _is_fasta_file(p: Path) -> bool:
    return p.suffix.lower() in _FASTA_EXTS


def _is_gb_file(p: Path) -> bool:
    return p.suffix.lower() in _GB_EXTS


def _sniff_input_type(in_path: Path) -> str:
    if in_path.is_file():
        if _is_fasta_file(in_path):
            return "fasta"
        if _is_gb_file(in_path):
            return "genbank"
        try:
            with in_path.open(encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    line = line.strip()
                    if line:
                        if line.startswith(">"):
                            return "fasta"
                        if line.startswith("LOCUS"):
                            return "genbank"
                        break
        except Exception:
            pass
        return "fasta"
    if in_path.is_dir():
        has_fasta = any(_is_fasta_file(p) for p in in_path.rglob("*") if p.is_file())
        has_gb    = any(_is_gb_file(p)    for p in in_path.rglob("*") if p.is_file())
        if has_gb and not has_fasta:
            return "genbank"
        if has_fasta and not has_gb:
            return "fasta"
        if has_gb and has_fasta:
            return "mixed"
        return "fasta"
    raise FileNotFoundError(f"Input not found: {in_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Output helpers — only write non-empty content
# ─────────────────────────────────────────────────────────────────────────────
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


def int_mem_mb(x, minimum: int = 128) -> int:
    try:
        v = int(math.ceil(float(x)))
    except Exception:
        v = minimum
    return max(minimum, v)


def positive_int_arg(x: str) -> int:
    try:
        v = int(x)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(f"expected positive integer, got {x!r}") from exc
    if v < 1:
        raise argparse.ArgumentTypeError(f"expected positive integer, got {x!r}")
    return v


def id_maxlen_arg(x: str) -> int:
    v = positive_int_arg(x)
    if v < len("S00001"):
        raise argparse.ArgumentTypeError(
            f"--id_maxlen must be >= {len('S00001')} to keep internal IDs unique")
    return v


def _safe_copy(src: Path, dst: Path) -> None:
    """Copy *src* -> *dst*; src and dst being the same file is not an error (e.g. --out_root = the input folder)."""
    try:
        shutil.copy2(str(src), str(dst))
    except shutil.SameFileError:
        pass


# ─────────────────────────────────────────────────────────────────────────────
# GenBank record accession helpers
# ─────────────────────────────────────────────────────────────────────────────
def gb_record_accession(rec: "SeqRecord") -> str:
    ver  = rec.annotations.get("sequence_version", "")
    accs = rec.annotations.get("accessions", [])
    acc  = accs[0] if accs else ""
    if not acc:
        acc = rec.name if rec.name and rec.name not in (".", "<unknown>") else ""
    if not acc:
        acc = rec.id  if rec.id  and rec.id  not in (".", "<unknown>") else ""
    if ver and acc and not acc.endswith(f".{ver}"):
        versioned = f"{acc}.{ver}"
    else:
        versioned = acc
    return safe_name(versioned or rec.id or "unknown")


@dataclass
class GbRecordMeta:
    acc_versioned: str
    acc_base:      str
    version_int:   int
    seq_hash:      str
    seq_len:       int
    topology:      str
    source_file:   Path
    record:        object


def _gb_version_int(rec: "SeqRecord") -> int:
    ver = rec.annotations.get("sequence_version", "")
    try:
        return int(ver)
    except (TypeError, ValueError):
        pass
    m = re.search(r"\.(\d+)$", rec.id or "")
    return int(m.group(1)) if m else 0


REFSEQ_ACC_RE = re.compile(r"^[A-Z]{2}_\d")
REFSEQ_SOURCE_RE = re.compile(
    r"(?:identical to|derived from)\s+((?:[A-Z]{1,6}\d{5,9}(?:\.\d+)?(?:\s*,\s*|\s+and\s+)?)+)", re.I)


def is_refseq_acc(acc: str) -> bool:
    """RefSeq accessions have a two-letter prefix and an underscore (NC_, NR_, ...)."""
    return bool(REFSEQ_ACC_RE.match(acc or ""))


def refseq_source_accessions(rec: "SeqRecord") -> List[str]:
    """INSDC accessions a RefSeq record was made from, read from its COMMENT
    ('The reference sequence is identical to MN123456' / 'was derived from ...')."""
    comment = " ".join(str(rec.annotations.get("comment", "")).split())
    out: List[str] = []
    for m in REFSEQ_SOURCE_RE.finditer(comment):
        out += [a.split(".")[0] for a in re.findall(r"[A-Z]{1,6}\d{5,9}(?:\.\d+)?", m.group(1))]
    return out


def dedup_gb_records(
    gb_files: List[Path],
    dedup_report_path: Optional[Path] = None,
) -> Tuple[List[GbRecordMeta], Dict[str, Path], Dict]:
    import hashlib

    stats: Dict = {
        "total_in": 0, "kept": 0,
        "dup_exact": 0, "dup_version": 0, "dup_seq": 0, "dup_refseq": 0,
        "parse_errors": 0, "no_sequence": 0,
        "dup_exact_list": [], "dup_version_list": [], "dup_seq_list": [], "dup_refseq_list": [],
        "no_sequence_list": [],
    }

    all_metas:               List[GbRecordMeta] = []
    metadata_only_acc_to_gb: Dict[str, Path]    = {}
    dropped: List[Tuple[GbRecordMeta, str, str]] = []     # (record, reason, accession kept instead)

    for gb_path in gb_files:
        try:
            records = list(SeqIO.parse(str(gb_path), "genbank"))
        except Exception as exc:
            logging.warning(f"[DEDUP] Cannot parse {gb_path}: {exc}")
            stats["parse_errors"] += 1
            continue

        for rec in records:
            stats["total_in"] += 1
            acc_ver  = gb_record_accession(rec)
            acc_base = _normalize_acc(acc_ver) or acc_ver
            ver_int  = _gb_version_int(rec)

            try:
                seq_str = str(rec.seq).upper().replace(" ", "")
                has_seq = bool(seq_str)
            except Exception:
                seq_str, has_seq = "", False

            if not has_seq:
                stats["no_sequence"] += 1
                stats["no_sequence_list"].append(f"{acc_ver} ({gb_path.name})")
                metadata_only_acc_to_gb.setdefault(safe_name(acc_ver), gb_path.resolve())
                continue

            seq_hash = hashlib.md5(seq_str.encode()).hexdigest()
            topo     = rec.annotations.get("topology", "circular").lower()
            if topo not in ("circular", "linear"):
                topo = "circular"

            all_metas.append(GbRecordMeta(
                acc_versioned=acc_ver, acc_base=acc_base,
                version_int=ver_int, seq_hash=seq_hash,
                seq_len=len(seq_str), topology=topo,
                source_file=gb_path, record=rec,
            ))

    by_acc: Dict[str, GbRecordMeta] = {}
    for meta in all_metas:
        key = meta.acc_base
        if key not in by_acc:
            by_acc[key] = meta
        else:
            prev = by_acc[key]
            if meta.version_int > prev.version_int:
                kind = "dup_exact" if prev.acc_versioned == meta.acc_versioned else "dup_version"
                stats[kind] += 1
                stats[f"{kind}_list"].append(f"{prev.acc_versioned} superseded by {meta.acc_versioned}")
                dropped.append((prev, "older version" if kind == "dup_version" else "same accession listed twice",
                                meta.acc_versioned))
                by_acc[key] = meta
            else:
                kind = "dup_version" if meta.version_int < prev.version_int else "dup_exact"
                stats["dup_exact"] += 1
                stats["dup_exact_list"].append(
                    f"{meta.acc_versioned} (dup of {prev.acc_versioned})")
                dropped.append((meta, "older version" if kind == "dup_version" else "same accession listed twice",
                                prev.acc_versioned))

    # RefSeq (NC_...) records are curated copies of INSDC submissions: keep the RefSeq
    # record and drop the submission it was made from (named in its COMMENT).
    superseded: Dict[str, str] = {}
    for meta in by_acc.values():
        if is_refseq_acc(meta.acc_base):
            for src in refseq_source_accessions(meta.record):
                superseded.setdefault(src, meta.acc_versioned)
    for key in [k for k in by_acc if k in superseded]:
        stats["dup_refseq"] += 1
        stats["dup_refseq_list"].append(
            f"{by_acc[key].acc_versioned} superseded by RefSeq {superseded[key]} (derived from it)")
        dropped.append((by_acc[key], "INSDC submission that RefSeq record was derived from (its COMMENT)",
                        superseded[key]))
        del by_acc[key]

    seq_hash_to_acc: Dict[str, str] = {}
    kept_records: List[GbRecordMeta] = []
    # RefSeq first, so that among sequence-identical records the curated one is kept
    for meta in sorted(by_acc.values(), key=lambda m: not is_refseq_acc(m.acc_base)):
        if meta.seq_hash in seq_hash_to_acc:
            stats["dup_seq"] += 1
            stats["dup_seq_list"].append(
                f"{meta.acc_versioned} (seq identical to {seq_hash_to_acc[meta.seq_hash]}, "
                f"len={meta.seq_len}, dropped)")
            logging.warning(
                f"[DEDUP_SEQ] '{meta.acc_versioned}' has identical sequence to "
                f"'{seq_hash_to_acc[meta.seq_hash]}' (len={meta.seq_len}). Keeping only first (RefSeq preferred).")
            dropped.append((meta, "identical sequence", seq_hash_to_acc[meta.seq_hash]))
        else:
            seq_hash_to_acc[meta.seq_hash] = meta.acc_versioned
            kept_records.append(meta)

    stats["kept"] = len(kept_records)

    if dedup_report_path is not None and (kept_records or stats["total_in"] > 0):
        write_dedup_report(dedup_report_path, kept_records, dropped, stats)

    return kept_records, metadata_only_acc_to_gb, stats


def write_dedup_report(path: Path, kept: List[GbRecordMeta],
                       dropped: List[Tuple[GbRecordMeta, str, str]], stats: Dict) -> None:
    """dedup_report.tsv: one row per input record. status kept / dropped; reason = why a record was dropped;
    related = the accession kept instead (dropped rows) or the accessions this record replaced (kept rows)."""
    replaced: Dict[str, List[str]] = {}
    for meta, _reason, by in dropped:
        replaced.setdefault(by, []).append(meta.acc_versioned)

    def organism(meta: GbRecordMeta) -> str:
        return str(getattr(meta.record, "annotations", {}).get("organism", ""))

    rows = []
    for meta in kept:
        rows.append((meta.acc_versioned, "kept", "", ",".join(replaced.get(meta.acc_versioned, [])),
                     organism(meta), meta.seq_len, meta.topology, meta.source_file.name))
    for meta, reason, by in dropped:
        rows.append((meta.acc_versioned, "dropped", reason, by,
                     organism(meta), meta.seq_len, meta.topology, meta.source_file.name))
    rows.sort(key=lambda r: (r[4], r[0]))
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# GB deduplication report: one row per record with a sequence (kept / dropped); records without ORIGIN\n",
        "# and parse errors are counted in the total line only. Rules: (1) same accession -> newest version;\n",
        "# (2) RefSeq NC_ record kept, the INSDC submission it was derived from (named in its COMMENT) dropped;\n",
        "# (3) identical sequence -> one copy kept (RefSeq preferred). Different records of the same species\n",
        "# are NOT merged. related = accession kept instead (dropped) / accessions it replaced (kept).\n",
        f"# total_in={stats['total_in']}  kept={stats['kept']}  dropped={len(dropped)}  "
        f"no_seq(metadata_only)={stats['no_sequence']}  parse_errors={stats['parse_errors']}\n",
        "accession\tstatus\treason\trelated\torganism\tseq_len\ttopology\tsource_file\n",
    ]
    lines += ["\t".join(str(x) for x in r) + "\n" for r in rows]
    path.write_text("".join(lines), encoding="utf-8")


def convert_gb_to_fasta_for_annotation(
    gb_files: List[Path],
    out_dir: Path,
    dedup_report_path: Optional[Path] = None,
    orig_dir: Optional[Path] = None,
) -> Tuple[List[Path], Dict[str, Path]]:
    if SeqIO is None:
        raise RuntimeError("Biopython not installed: pip install biopython")
    out_dir.mkdir(parents=True, exist_ok=True)

    kept_records, metadata_only_acc_to_gb, stats = dedup_gb_records(
        gb_files, dedup_report_path=dedup_report_path)

    n_total  = stats["total_in"]
    n_kept   = stats["kept"]
    n_no_seq = stats["no_sequence"]
    n_dup    = n_total - n_kept - n_no_seq

    logging.info(
        f"[GB_DEDUP] total={n_total}  with_seq={n_kept}  "
        f"no_seq(metadata_only)={n_no_seq}  dropped_dup={n_dup}  "
        f"(dup_exact={stats['dup_exact']}, dup_version={stats['dup_version']}, "
        f"dup_seq={stats['dup_seq']}, dup_refseq={stats['dup_refseq']}, parse_errors={stats['parse_errors']})"
    )
    if n_no_seq > 0:
        logging.info(
            f"[GB_NO_SEQ] {n_no_seq} record(s) have no sequence (ORIGIN absent) "
            f"— skipped for annotation, metadata source still tracked.")
    for msg in stats["dup_seq_list"]:
        logging.warning(f"  [DUP_SEQ]     {msg}")

    fasta_paths:    List[Path]      = []
    acc_to_orig_gb: Dict[str, Path] = {}

    for meta in kept_records:
        fa_name = safe_name(meta.acc_versioned)
        out_fa  = out_dir / f"{fa_name}.fasta"
        out_rec = SeqRecord(
            meta.record.seq, id=fa_name, name=fa_name,
            description=f"topology={meta.topology}",
        )
        with out_fa.open("w", encoding="utf-8") as fh:
            SeqIO.write([out_rec], fh, "fasta")
        fasta_paths.append(out_fa.resolve())
        acc_to_orig_gb[fa_name] = meta.source_file.resolve()

    for acc_key, src_path in metadata_only_acc_to_gb.items():
        acc_to_orig_gb.setdefault(acc_key, src_path)

    if orig_dir is not None:
        # the input records unchanged, one file each, for the before / after comparison (run: original_annotation/)
        if orig_dir.is_dir():
            shutil.rmtree(orig_dir)                    # a previous run's records
        orig_dir.mkdir(parents=True, exist_ok=True)
        for meta in kept_records:
            try:
                SeqIO.write([meta.record], str(orig_dir / f"{safe_name(meta.acc_versioned)}.gb"), "genbank")
            except Exception as exc:               # noqa: BLE001  (a record Biopython cannot write back: compared later as missing)
                logging.warning(f"[ORIGINAL_GB] {meta.acc_versioned}: not written to {orig_dir.name}/: {exc!r}")
    LAST_DEDUP_STATS.clear()
    LAST_DEDUP_STATS.update(stats)
    return fasta_paths, acc_to_orig_gb


LAST_DEDUP_STATS: Dict = {}        # stats of the last convert_gb_to_fasta_for_annotation (run: step 1 line)


def iter_inputs(
    in_path: Path,
    gb_to_fasta_dir: Optional[Path] = None,
    orig_dir: Optional[Path] = None,
) -> Tuple[List[Path], str, Dict[str, Path]]:
    itype = _sniff_input_type(in_path)
    logging.info(f"[INPUT] Detected input type: {itype}  ({in_path})")

    if itype == "fasta":
        if in_path.is_file():
            return [in_path.resolve()], "fasta", {}
        files = sorted(p.resolve() for p in in_path.rglob("*")
                       if p.is_file() and _is_fasta_file(p))
        if not files:
            raise FileNotFoundError(f"No FASTA files in: {in_path}")
        return files, "fasta", {}

    if itype in ("genbank", "mixed"):
        if in_path.is_file():
            gb_files = [in_path.resolve()]
        else:
            gb_files = sorted(p.resolve() for p in in_path.rglob("*")
                              if p.is_file() and _is_gb_file(p))
        if not gb_files:
            raise FileNotFoundError(f"No GenBank files in: {in_path}")
        if itype == "mixed":
            logging.warning(
                "[INPUT] Mixed FASTA+GB directory — only GB files will be used.")
        if gb_to_fasta_dir is None:
            raise RuntimeError("gb_to_fasta_dir required when input is GenBank")
        dedup_report = gb_to_fasta_dir / "dedup_report.tsv"
        fasta_paths, acc_to_orig_gb = convert_gb_to_fasta_for_annotation(
            gb_files, gb_to_fasta_dir, dedup_report_path=dedup_report, orig_dir=orig_dir)
        if not fasta_paths:
            raise RuntimeError(f"No records extracted from GenBank input: {in_path}")
        logging.info(
            f"[INPUT] Converted {len(gb_files)} GB file(s) -> "
            f"{len(fasta_paths)} unique FASTA record(s)")
        return fasta_paths, "converted_from_gb", acc_to_orig_gb

    raise FileNotFoundError(f"Input not found or unrecognised: {in_path}")


MIN_CIRCULAR_LEN = 10000    # a FASTA record shorter than this, or one of several in a file, is a contig: linear


def infer_topology_from_description(desc: str) -> str:
    m = re.search(r"topology\s*=\s*(circular|linear)", (desc or "").lower())
    return m.group(1) if m else "circular"


def _is_topology_tag(s: str) -> bool:
    return bool(re.match(r"^topology\s*=\s*(circular|linear)$", (s or "").strip(), re.I))


def make_internal_id(idx: int, max_len: int) -> str:
    base = f"S{idx:05d}"
    return base[:max_len] if len(base) > max_len else base


# ─────────────────────────────────────────────────────────────────────────────
# SampleInfo / RunSummary
# ─────────────────────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class SampleInfo:
    sample_name:    str
    orig_token:     str
    internal_id:    str
    topology:       str
    source_fasta:   Path
    internal_fasta: Path
    orig_gb_path:   Optional[Path] = None


@dataclass
class RunSummary:
    annotation_failed:        List[str] = field(default_factory=list)
    final_write_failed:       List[str] = field(default_factory=list)
    no_target_gene:           List[str] = field(default_factory=list)
    metadata_transferred:     List[str] = field(default_factory=list)
    metadata_compat_match:    List[str] = field(default_factory=list)
    metadata_key_mismatch:    List[str] = field(default_factory=list)
    metadata_skipped:         List[str] = field(default_factory=list)
    metadata_transfer_failed: List[str] = field(default_factory=list)
    reorient_failed:          List[str] = field(default_factory=list)       # reorient workflow (nothing annotated)

    def all_problem_samples(self) -> List[str]:
        return sorted(set(
            self.annotation_failed + self.final_write_failed
            + self.no_target_gene  + self.metadata_transfer_failed
            + self.metadata_key_mismatch                           # header matched only by a compatible accession
            + self.reorient_failed
        ))

    def log(self) -> None:
        logging.info(
            "[SUMMARY] annotation_failed=%d  final_write_failed=%d  no_target_gene=%d"
            "  reorient_failed=%d",
            len(self.annotation_failed), len(self.final_write_failed), len(self.no_target_gene),
            len(self.reorient_failed))
        logging.info(
            "[SUMMARY] metadata_transferred=%d  compat_match=%d  key_mismatch=%d  "
            "skipped=%d  transfer_failed=%d",
            len(self.metadata_transferred), len(self.metadata_compat_match),
            len(self.metadata_key_mismatch), len(self.metadata_skipped),
            len(self.metadata_transfer_failed))


# ─────────────────────────────────────────────────────────────────────────────
# FASTA helpers
# ─────────────────────────────────────────────────────────────────────────────
def write_one_record_fasta(record: "SeqRecord", out_fa: Path, seq_id: str, topology: str) -> None:
    rec = record[:]
    rec.id = rec.name = seq_id
    rec.description = f"topology={topology}"
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with out_fa.open("w", encoding="utf-8") as fh:
        SeqIO.write([rec], fh, "fasta")


def split_or_sanitize_fastas(
    in_files: List[Path],
    out_dir: Path,
    if_split: bool,
    id_maxlen: int,
    acc_to_orig_gb: Optional[Dict[str, Path]] = None,
) -> Tuple[List[SampleInfo], Path]:
    if SeqIO is None:
        raise RuntimeError("Biopython not installed: pip install biopython")
    out_dir.mkdir(parents=True, exist_ok=True)
    if acc_to_orig_gb is None:
        acc_to_orig_gb = {}
    sample_counts: Dict[str, int] = {}
    _seen_sample_names: set = set()
    samples: List[SampleInfo] = []
    idx = 0
    id_map_lines = [
        "internal_id\tsample_name\torig_token\ttopology\tsource_fasta\tinternal_fasta\torig_gb_path\n"
    ]
    for fa in in_files:
        try:
            records = list(SeqIO.parse(str(fa), "fasta"))
        except Exception as exc:
            logging.warning(f"[SPLIT] Cannot parse {fa}: {exc}")
            continue
        if not records:
            logging.warning(f"[SPLIT] No records in {fa}, skipping.")
            continue
        if (not if_split) and len(records) > 1:
            raise ValueError(f"FASTA has >1 record: {fa}. Use --if_split_fasta yes.")
        for rec in records:
            idx += 1
            desc = (rec.description or "").strip()
            orig_token = header_first_token(rec.id) if _is_topology_tag(desc) \
                         else header_first_token(desc)
            base_name = safe_name(orig_token)
            sample_counts[base_name] = sample_counts.get(base_name, 0) + 1
            cnt = sample_counts[base_name]
            cand_name = base_name if cnt == 1 else f"{base_name}_{cnt}"
            while cand_name in _seen_sample_names:
                sample_counts[base_name] += 1
                cand_name = f"{base_name}_{sample_counts[base_name]}"
            _seen_sample_names.add(cand_name)
            sample_name = cand_name
            topology    = infer_topology_from_description(desc)
            if topology == "circular" and not re.search(r"topology\s*=\s*circular", desc.lower()) \
                    and (len(records) > 1 or len(rec.seq) < MIN_CIRCULAR_LEN):
                # contigs of an assembly: rotating / joining across a false origin cuts real genes
                topology = "linear"
                logging.warning(f"[TOPOLOGY] {orig_token}: {len(records)} records / {len(rec.seq)} nt — treated "
                                "as linear; add topology=circular to the header to override")
            internal_id = make_internal_id(idx, id_maxlen)
            out_fa      = out_dir / f"{internal_id}.fasta"
            write_one_record_fasta(rec, out_fa, seq_id=internal_id, topology=topology)
            orig_gb_path = acc_to_orig_gb.get(safe_name(rec.id)) \
                        or acc_to_orig_gb.get(orig_token)
            si = SampleInfo(
                sample_name=sample_name, orig_token=orig_token,
                internal_id=internal_id, topology=topology,
                source_fasta=fa.resolve(), internal_fasta=out_fa.resolve(),
                orig_gb_path=orig_gb_path,
            )
            samples.append(si)
            id_map_lines.append(
                f"{internal_id}\t{sample_name}\t{orig_token}\t{topology}\t"
                f"{si.source_fasta}\t{si.internal_fasta}\t{orig_gb_path or 'NA'}\n"
            )
    id_map_path = out_dir / "id_map.tsv"
    id_map_path.write_text("".join(id_map_lines), encoding="utf-8")
    return samples, id_map_path


def parse_suffix_fq(s: str) -> Tuple[str, str]:
    parts = [p.strip() for p in (s or "").split(",") if p.strip()]
    if len(parts) != 2:
        raise ValueError("--suffix_fq needs two comma-separated suffixes")
    return parts[0], parts[1]


def find_fq_pair_for_sample(
    sample_name: str, infq: Path, fq_position: str, suffix1: str, suffix2: str,
) -> Tuple[Optional[Path], Optional[Path]]:
    sn = str(sample_name)
    if fq_position == "flat":
        f1 = infq / f"{sn}{suffix1}"
        f2 = infq / f"{sn}{suffix2}"
        return (f1.resolve() if f1.exists() else None), (f2.resolve() if f2.exists() else None)
    if fq_position == "sample_dir":
        d = infq / sn
        if not d.exists():
            return None, None
        f1 = d / f"{sn}{suffix1}"
        f2 = d / f"{sn}{suffix2}"
        if f1.exists() and f2.exists():
            return f1.resolve(), f2.resolve()
        c1 = sorted(d.glob(f"*{suffix1}"))
        c2 = sorted(d.glob(f"*{suffix2}"))
        return (c1[0].resolve() if c1 else None), (c2[0].resolve() if c2 else None)
    raise ValueError(f"Unknown --fq_position: {fq_position}")


def run_cmd_capture(
    cmd: List[str], cwd: Optional[Path], stdout_path: Path, stderr_path: Path,
) -> Tuple[int, float]:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    p  = subprocess.run(cmd, cwd=str(cwd) if cwd else None,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    dt = time.time() - t0
    if p.stdout:
        stdout_path.write_text(p.stdout, encoding="utf-8", errors="ignore")
    if p.stderr:
        stderr_path.write_text(p.stderr, encoding="utf-8", errors="ignore")
    return int(p.returncode), float(dt)


# ─────────────────────────────────────────────────────────────────────────────
# GenBank helpers
# ─────────────────────────────────────────────────────────────────────────────
def collect_genbank_files(src: Path) -> List[Path]:
    if src.is_file():
        return [src.resolve()]
    if src.is_dir():
        return sorted(
            p.resolve() for p in src.rglob("*")
            if p.is_file() and _is_gb_file(p)
        )
    return []


def merge_gb_files(gb_files: List[Path], merged_out: Path) -> int:
    """The files' text one after the other (each checked by parsing it): a Biopython round trip would write the
    value-less /pseudo of mito_curate P10 as /pseudo=""."""
    merged_out.parent.mkdir(parents=True, exist_ok=True)
    n, chunks = 0, []
    for gb in gb_files:
        try:
            k = sum(1 for _ in SeqIO.parse(str(gb), "genbank"))
        except Exception as exc:
            logging.warning(f"[MERGE_SKIP] {gb}: {exc}")
            continue
        if k:
            text = gb.read_text()
            chunks.append(text if text.endswith("\n") else text + "\n")
            n += k
    if chunks:
        merged_out.write_text("".join(chunks))
    return n


def refresh_summary(out_root: Path, gb_dir: Path) -> None:
    """Recount PCG / tRNA / rRNA of annotation-summary.tsv from the final GenBank files (after rescue, rRNA
    HMM fallback and curation); rows without a final file are kept unchanged."""
    f = out_root / "annotation-summary.tsv"
    if not f.exists():
        return
    finals = {p.stem: p for p in Path(gb_dir).rglob("*.gbf")}
    lines = f.read_text(encoding="utf-8").splitlines()
    out = lines[:1]
    for line in lines[1:]:
        c = line.split("\t")
        if c and c[0] in finals:
            c[1:4] = [str(x) for x in counts_from_genbank(finals[c[0]])]
        out.append("\t".join(c))
    f.write_text("\n".join(out) + "\n", encoding="utf-8")


def counts_from_genbank(gb_path: Path) -> Tuple[int, int, int]:
    pcg = trna = rrna = 0
    for rec in SeqIO.parse(str(gb_path), "genbank"):
        for feat in rec.features:
            if feat.type == "CDS":
                pcg  += 1
            elif feat.type == "tRNA":
                trna += 1
            elif feat.type == "rRNA":
                rrna += 1
    return pcg, trna, rrna


_LOCUS_NAME_MAXLEN = 16


def rewrite_genbank_record_ids(gb_in: Path, gb_out: Path, new_base_id: str) -> Tuple[int, int, int]:
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    if not recs:
        raise ValueError(f"No records in {gb_in}")
    pcg = trna = rrna = 0
    out_recs: List["SeqRecord"] = []
    for i, r in enumerate(recs, start=1):
        rid = new_base_id if i == 1 else f"{new_base_id}_{i}"
        locus = rid[:_LOCUS_NAME_MAXLEN]                  # a new record, not r[:] (which copies every feature)
        if len(rid) > _LOCUS_NAME_MAXLEN:
            logging.info(
                f"[LOCUS_TRUNCATE] '{rid}' ({len(rid)} chars) truncated to '{locus}' "
                f"in LOCUS line (GenBank limit={_LOCUS_NAME_MAXLEN}). "
                f"ACCESSION/VERSION unaffected.")
        new_annotations = dict(r.annotations)
        new_annotations["accessions"] = [rid]
        new_annotations.pop("sequence_version", None)
        rr = SeqRecord(
            r.seq,
            id=rid,
            name=locus,
            description=r.description or "",
            dbxrefs=list(getattr(r, "dbxrefs", [])),
            annotations=new_annotations,
            features=list(r.features),
        )
        out_recs.append(rr)
        for feat in r.features:
            if feat.type == "CDS":
                pcg  += 1
            elif feat.type == "tRNA":
                trna += 1
            elif feat.type == "rRNA":
                rrna += 1
    gb_out.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_recs, str(gb_out), "genbank")
    return pcg, trna, rrna


def gb_to_internal_fasta(gb_in: Path, fasta_out: Path, internal_id: str) -> None:
    """Convert a GenBank file to FASTA for Phase-2 re-annotation input: all records, the second contig as
    <internal_id>_2 and so on."""
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    if not recs:
        raise ValueError(f"No records in {gb_in}")
    out_recs = []
    for i, rec in enumerate(recs):
        topo   = rec.annotations.get("topology", "circular")
        seq_id = internal_id if i == 0 else f"{internal_id}_{i + 1}"
        out_recs.append(SeqRecord(rec.seq, id=seq_id, name=seq_id,
                                  description=f"topology={topo}"))
    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    with fasta_out.open("w", encoding="utf-8") as fh:
        SeqIO.write(out_recs, fh, "fasta")


# ─────────────────────────────────────────────────────────────────────────────
# GenBank reorientation
# ─────────────────────────────────────────────────────────────────────────────
def normalize_gene_key(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (name or "").strip().lower())


_GENE_ALIAS_GROUPS = [
    ["cox1", "coi", "coxi", "co1", "cytochromecoxidasesubuniti", "cytochromecoxidasesubunit1"],
    ["cox2", "coii", "coxii", "co2", "cytochromecoxidasesubunitii", "cytochromecoxidasesubunit2"],
    ["cox3", "coiii", "coxiii", "co3", "cytochromecoxidasesubunitiii", "cytochromecoxidasesubunit3"],
    ["cob", "cytb", "cyb", "cytochromeb"],
    ["nad1", "nd1", "nadh1", "nadhdehydrogenasesubunit1"],
    ["atp6", "atpase6", "atpsynthasef0subunit6", "atpsynthasesubunit6"],
]


def gene_aliases(gene: str) -> List[str]:
    g = normalize_gene_key(gene)
    for group in _GENE_ALIAS_GROUPS:
        if g in group:
            return list(group)
    return [g]


# Records oriented by a fallback gene because the first gene was absent: {record id: gene}
ORIENTATION_FALLBACK: Dict[str, str] = {}


def orientation_genes(args: argparse.Namespace) -> str:
    """'cox1' plus --fallback_genes, as the comma list reorient_one_gb tries in order."""
    fb = [g.strip() for g in (getattr(args, "fallback_genes", "") or "").split(",")
          if g.strip() and g.strip().lower() != "none"]
    return ",".join([args.which_gene_first] + [g for g in fb if g != args.which_gene_first])


def find_orientation_feature(record: "SeqRecord", genes: str) -> Tuple[Optional["SeqFeature"], str]:
    """First gene of the comma list present in the record; when the gene is annotated in
    several pieces (e.g. exons around an intron), the 5'-most piece in its own orientation."""
    for gene in [g.strip() for g in genes.split(",") if g.strip()]:
        aliases = set(gene_aliases(gene))
        hits = [f for f in record.features if f.type in {"gene", "CDS", "tRNA", "rRNA"}
                and any(normalize_gene_key(n) in aliases for n in feature_gene_names(f))]
        if not hits:
            continue
        strand = hits[0].location.strand
        same = [f for f in hits if f.location.strand == strand]
        if strand == -1:
            return max(same, key=lambda f: int(f.location.end)), gene
        return min(same, key=lambda f: int(f.location.start)), gene
    return None, ""


def feature_gene_names(feature: "SeqFeature") -> List[str]:
    names: List[str] = []
    for k in ("gene", "product", "label"):
        for v in feature.qualifiers.get(k, []):
            names.append(str(v))
    return names


def _rc_location(loc, seq_len: int):
    if isinstance(loc, CompoundLocation):
        return CompoundLocation(
            [_rc_location(p, seq_len) for p in loc.parts][::-1], operator=loc.operator)
    start, end = int(loc.start), int(loc.end)
    new_strand  = None if loc.strand is None else (-loc.strand)
    return FeatureLocation(seq_len - end, seq_len - start, strand=new_strand)


def reverse_complement_record(record: "SeqRecord") -> "SeqRecord":
    seq_len = len(record.seq)
    new_rec = SeqRecord(record.seq.reverse_complement(), id=record.id,
                        name=record.name, description=record.description)
    new_rec.annotations = dict(record.annotations)
    new_rec.dbxrefs     = list(getattr(record, "dbxrefs", []))
    new_features: List["SeqFeature"] = []
    for feat in record.features:
        try:
            new_loc = _rc_location(feat.location, seq_len)
            new_features.append(SeqFeature(location=new_loc, type=feat.type,
                                            id=getattr(feat, "id", "<unknown>"),
                                            qualifiers=dict(feat.qualifiers)))
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


_LINEAR_SKIP_SEEN: set = set()      # sequences already reported by [LINEAR_REORIENT_SKIP] (run: step 3, then the final write)


def rotate_record(record: "SeqRecord", start_pos: int) -> "SeqRecord":
    seq_len = len(record.seq)
    topology = record.annotations.get("topology", "circular")
    if str(topology).lower() == "linear" and start_pos != 0:
        key = hash(str(record.seq))
        (logging.info if key in _LINEAR_SKIP_SEEN else logging.warning)(
            f"[LINEAR_REORIENT_SKIP] record={record.id}: target start={start_pos}; "
            "linear topology cannot be circularly rotated, sequence left unchanged.")
        _LINEAR_SKIP_SEEN.add(key)
        return SeqRecord(                     # moving a linear contig's tail before its head would corrupt it
            record.seq,
            id=record.id, name=record.name, description=record.description,
            dbxrefs=list(getattr(record, "dbxrefs", [])),
            annotations=dict(record.annotations),
            features=list(record.features),
        )
    new_seq = record.seq[start_pos:] + record.seq[:start_pos]
    new_rec = SeqRecord(
        Seq(str(new_seq)),
        id=record.id, name=record.name, description=record.description,
        dbxrefs=list(getattr(record, "dbxrefs", [])),
        annotations={**record.annotations, "topology": topology},
    )
    new_features: List["SeqFeature"] = []
    for feat in record.features:
        try:
            new_loc = shift_location(feat.location, start_pos, seq_len)
            new_features.append(SeqFeature(location=new_loc, type=feat.type,
                                            id=getattr(feat, "id", "<unknown>"),
                                            qualifiers=dict(feat.qualifiers)))
        except Exception:
            new_features.append(feat)
    new_features.sort(key=lambda f: feature_start_pos(f))
    new_rec.features = new_features
    return new_rec


# --if_reannotation_2 auto (user 2026-10-21: re-annotating every record after reorientation wastes resources)
# ATP8 and ATP6 are NOT in the list: absent in whole phyla (nematodes, flatworms, many bivalves, anthozoans...),
# so their absence says nothing about a cut gene (review B finding 4).
STANDARD_GENES = ("COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
                  "s-rRNA", "l-rRNA")
REANNOT_EDGE = 50          # nt from the ends of the round-1 sequence where a gene may have been cut
TRNA_LEN = (55, 95)        # complete tRNA genes are this long (nt); outside = cut or merged
RRNA_MIN = {"s": 600, "l": 900}   # shorter s-rRNA / l-rRNA (nt) at the sequence ends = cut


def edge_gene_truncated(f: "SeqFeature", seq: "Seq", total: int, code: int = 5) -> bool:
    """A CDS / rRNA / tRNA feature within REANNOT_EDGE nt of an end of the round-1 sequence (length `total`):
    does it look cut there? (A complete gene that merely starts at 1 or ends at the last base is not: NCBI
    vertebrate records start with trnF.) CDS: length not a multiple of 3 (except a T / TA stop tail, the
    mitochondrial incomplete stop codon), or no start codon where its 5' end is the sequence end, or no stop
    codon in the feature nor in frame within REANNOT_EDGE nt behind it where its 3' end is the sequence end
    (codon tables of `code`); tRNA: outside 55-95 nt; rRNA: s-rRNA < 600 nt, l-rRNA < 900 nt (unknown rRNA
    name: < 600). A join across the origin is an origin-spanning feature, judged like any interior one."""
    loc = f.location
    wraps = len(loc.parts) > 1 and any(int(p.end) == total for p in loc.parts) \
        and any(int(p.start) == 0 for p in loc.parts)       # join across the origin: annotated as such, its ends are inside
    left, right = (False, False) if wraps else (int(loc.start) < REANNOT_EDGE, int(loc.end) > total - REANNOT_EDGE)
    if not (left or right or wraps):
        return False
    n = len(loc)
    if f.type == "tRNA":
        return not TRNA_LEN[0] <= n <= TRNA_LEN[1]
    if f.type == "rRNA":
        name = " ".join(f.qualifiers.get("gene", []) + f.qualifiers.get("product", [])).lower()
        large = any(k in name for k in ("l-rrna", "rrnl", "16s", "large", "lsu"))
        return n < RRNA_MIN["l" if large else "s"]
    if left and right:                         # spans the whole sequence: cannot judge, re-annotate
        return True
    from Bio.Data import CodonTable
    table = CodonTable.unambiguous_dna_by_id.get(int(code), CodonTable.unambiguous_dna_by_id[5])
    nt = str(f.extract(seq)).upper()
    rem = len(nt) % 3
    tail = (rem == 1 and nt.endswith("T")) or (rem == 2 and nt.endswith("TA"))
    if rem and not tail:
        return True
    minus = loc.strand == -1
    five_at_edge = right if minus else left       # which end of the gene is a sequence end
    three_at_edge = left if minus else right
    if five_at_edge and nt[:3] not in table.start_codons:
        return True
    if three_at_edge and not tail and nt[-3:] not in table.stop_codons:
        # the core annotates a CDS up to its last aligned codon (the stop is added later): a stop in frame within
        # REANNOT_EDGE nt behind it, or a T / TA tail ending exactly at the sequence end, makes it complete
        if minus:
            beyond = str(seq[max(0, int(loc.start) - REANNOT_EDGE):int(loc.start)].reverse_complement()).upper()
            reaches_end = int(loc.start) <= REANNOT_EDGE
        else:
            beyond = str(seq[int(loc.end):int(loc.end) + REANNOT_EDGE]).upper()
            reaches_end = int(loc.end) + REANNOT_EDGE >= total
        codons = [beyond[i:i + 3] for i in range(0, len(beyond) - 2, 3)]
        if any(c in table.stop_codons for c in codons):
            return False
        return not (reaches_end and beyond[len(codons) * 3:] in ("T", "TA"))
    return False


def select_reannotation(samples: List["SampleInfo"], results1: Dict[str, Dict], out_root: Path,
                        code: int = 5) -> List["SampleInfo"]:
    """--if_reannotation_2 auto: the samples whose reoriented round-1 annotation fails reorient_check()
    (genetic code `code`); every decision with its reasons goes to reannotation_decision.tsv."""
    chosen, rows = [], []
    for si in samples:
        iid = si.internal_id
        reor = out_root / "01.anno_round1" / iid / "reoriented" / f"{iid}.reoriented.gbf"
        best = results1.get(iid, {}).get("best_gb")
        try:
            probs = reorient_check(Path(best), reor, code) if best else ["no round-1 GenBank"]
        except Exception as exc:                     # any doubt -> re-annotate
            probs = [f"check failed: {exc}"]
        if probs:
            chosen.append(si)
        rows.append(f"{si.sample_name}\t{'reannotate' if probs else 'shifted_round1'}\t{'; '.join(probs)}")
    write_nonempty_table(out_root / "reannotation_decision.tsv", "sample_name\tdecision\treasons", rows)
    logging.info(f"[PHASE2] auto: {len(chosen)} of {len(samples)} sample(s) re-annotated, "
                 f"{len(samples) - len(chosen)} use the shifted round-1 annotation (reannotation_decision.tsv)")
    return chosen


def reorient_check(round1_gb: Path, reor_gb: Path, code: int = 5) -> List[str]:
    """Reasons to re-annotate a reoriented record; [] = the shifted round-1 annotation can be used as is.
    Checks: (1) a gene within REANNOT_EDGE nt of the ends of the round-1 sequence that looks cut there
    (edge_gene_truncated, genetic code `code`; the core annotates the sequence as linear, so a gene across the
    original origin is cut or missed) - only when the sequence moved; (2) a standard protein gene / rRNA
    (STANDARD_GENES; not ATP8 / ATP6) missing in round 1; (3) every feature's sequence identical
    before and after the shift, same features; (4) no feature across the new origin. A linear record is never rotated
    (only turned to the plus strand): its end genes are cut by the real sequence ends and a missing gene is missing
    from the sequence, so (1) and (2) do not apply - a second round on the same sequence finds nothing new (2026-10
    benchmark, 15 linear annelid records: the same true tRNAs, 2 more false ones from the re-run)."""
    probs: List[str] = []
    r1 = {r.id: r for r in SeqIO.parse(str(round1_gb), "genbank")}
    for ro in SeqIO.parse(str(reor_gb), "genbank"):
        a = r1.get(ro.id)
        if a is None:
            probs.append(f"{ro.id}: no round-1 record")
            continue
        L = len(a.seq)
        genes = [f for f in a.features if f.type in ("CDS", "rRNA", "tRNA")]
        moved = str(ro.seq) != str(a.seq)
        linear = "linear" in (a.annotations.get("topology"), ro.annotations.get("topology"))
        if moved and not linear:
            edge = sorted({f.qualifiers.get("gene", [f.type])[0] for f in genes
                           if edge_gene_truncated(f, a.seq, L, code)})
            if edge:
                probs.append(f"{ro.id}: gene(s) at the end of the round-1 sequence, look cut: {', '.join(edge)}")
        have = {f.qualifiers.get("gene", [""])[0] for f in genes}
        miss = [g for g in STANDARD_GENES if g not in have]
        if miss and not linear:
            probs.append(f"{ro.id}: standard gene(s) missing in round 1: {', '.join(miss)}")

        def seqs(rec: "SeqRecord") -> Dict[Tuple[str, str], List[str]]:
            out: Dict[Tuple[str, str], List[str]] = {}
            for f in rec.features:
                if f.type in ("CDS", "rRNA", "tRNA"):
                    out.setdefault((f.type, f.qualifiers.get("gene", [""])[0]), []).append(str(f.extract(rec.seq)))
            return {k: sorted(v) for k, v in out.items()}
        before, after = seqs(a), seqs(ro)
        diff = sorted(k[1] or k[0] for k in set(before) | set(after) if before.get(k) != after.get(k))
        if diff:
            probs.append(f"{ro.id}: gene sequence differs after the shift: {', '.join(diff)}")
        n = len(ro.seq)
        cross = sorted({f.qualifiers.get("gene", [f.type])[0] for f in ro.features if f.type != "source"
                        and len(f.location.parts) > 1 and any(int(p.end) == n for p in f.location.parts)
                        and any(int(p.start) == 0 for p in f.location.parts)})
        if cross:
            probs.append(f"{ro.id}: feature(s) across the new origin: {', '.join(cross)}")
    return probs


def feature_start_pos(feat: "SeqFeature") -> int:
    loc = feat.location
    if isinstance(loc, CompoundLocation) and loc.parts:
        return int(loc.parts[0].start)
    return int(loc.start)


def reorient_one_gb(
    gb_in: Path, gb_out: Path, gb_no_gene: Path,
    which_gene_first: str, force_forward: bool,
) -> Tuple[int, int, List[str]]:
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    if not recs:
        raise ValueError(f"No GenBank records: {gb_in}")
    ok: List["SeqRecord"] = []
    no_gene: List["SeqRecord"] = []
    no_gene_ids: List[str] = []
    for rec in recs:
        try:
            feat, used = find_orientation_feature(rec, which_gene_first)
            if feat is None:
                no_gene.append(rec)
                no_gene_ids.append(rec.id)
                continue
            first = which_gene_first.split(",")[0].strip()
            if used != first:
                ORIENTATION_FALLBACK[rec.id] = used
                logging.info(f"[REORIENT_FALLBACK] record={rec.id}: '{first}' not found, "
                             f"sequence starts at '{used}'")
            work_rec = rec
            if force_forward and feat.location.strand == -1:
                work_rec = reverse_complement_record(rec)
                feat, _ = find_orientation_feature(work_rec, used)
                if feat is None:
                    no_gene.append(rec)
                    no_gene_ids.append(rec.id)
                    continue
            ok.append(rotate_record(work_rec, feature_start_pos(feat)))
        except Exception as exc:
            logging.warning(f"[REORIENT_REC_FAIL] record={rec.id}: {exc}")
            no_gene.append(rec)
            no_gene_ids.append(rec.id)
    gb_out.parent.mkdir(parents=True, exist_ok=True)
    gb_no_gene.parent.mkdir(parents=True, exist_ok=True)
    if ok:
        SeqIO.write(ok, str(gb_out), "genbank")
    if no_gene:
        SeqIO.write(no_gene, str(gb_no_gene), "genbank")
    return len(ok), len(no_gene), no_gene_ids


# ─────────────────────────────────────────────────────────────────────────────
# Metadata transfer engine
# ─────────────────────────────────────────────────────────────────────────────
_REFSEQ_PREFIX_PAT = re.compile(r"^([A-Z]{2}_)")
# only accession-shaped IDs lose a trailing version: local sample names such as B64_2 are distinct samples
_NON_REFSEQ_VERSION_PAT = re.compile(r"^[A-Z]{1,4}\d{5,}[._]\d+$", re.I)
_LOCUS_BPCOUNT_RE  = re.compile(r"(LOCUS\s+\S+\s+)\d+(\s+[a-zA-Z]+)")      # LOCUS bp count


def _clean_acc(acc: Optional[str]) -> Optional[str]:
    if not acc:
        return None
    return acc.strip().rstrip(";").strip().replace(" ", "") or None


def _normalize_acc(acc: Optional[str]) -> Optional[str]:
    s = _clean_acc(acc)
    if not s:
        return None
    m = _REFSEQ_PREFIX_PAT.match(s)
    if m:
        prefix = m.group(1)
        body   = re.sub(r"[._]\d+$", "", s[len(prefix):])
        s      = prefix + body
    else:
        if _NON_REFSEQ_VERSION_PAT.match(s):
            s = re.sub(r"[._]\d+$", "", s)
    return s.rstrip(";").strip() or None


@functools.lru_cache(maxsize=4096)
def _acc_variants(acc: Optional[str]) -> Tuple[str, ...]:
    """All meaningful accession variants for lookup matching (cached; a tuple, so the cached value cannot be
    changed by a caller)."""
    clean = _clean_acc(acc)
    if not clean:
        return ()
    norm = _normalize_acc(clean)
    seen: List[str] = []

    def _add(x: Optional[str]) -> None:
        if x and x not in seen:
            seen.append(x)

    _add(clean)
    _add(norm)
    m = re.search(r"[._](\d+)$", clean)
    if m and norm and norm != clean:             # a version was removed (B64_2 must not gain B64 as an alias)
        ver = m.group(1)
        _add(f"{norm}.{ver}")
        _add(f"{norm}_{ver}")
    return tuple(seen)


def _split_gb_records(text: str) -> List[str]:
    """Split a multi-record GenBank text into individual record strings (no empty pseudo-record from newlines after
    the last //)."""
    records: List[str] = []
    buf: List[str] = []
    for line in text.splitlines(keepends=True):
        buf.append(line)
        if line.strip() == "//":
            chunk = "".join(buf)
            if chunk.strip():          # skip empty / whitespace-only chunks
                records.append(chunk)
            buf = []
    if buf:
        chunk = "".join(buf)
        if chunk.strip():
            records.append(chunk)
    return records


def _split_meta_body(record: str) -> Tuple[str, str]:
    lines = record.splitlines(keepends=True)
    meta: List[str] = []
    body: List[str] = []
    in_body = False
    for line in lines:
        if not in_body and (line.startswith("FEATURES") or line.startswith("ORIGIN")):
            in_body = True
        (body if in_body else meta).append(line)
    return "".join(meta), "".join(body)


def _get_record_accession(record: str) -> Optional[str]:
    """Accession from GB text: VERSION first, then ACCESSION, then LOCUS (Biopython truncates the LOCUS name to 16
    characters)."""
    for prefix in ("VERSION", "ACCESSION", "LOCUS"):
        for line in record.splitlines():
            if line.startswith(prefix):
                parts = line.split()
                if len(parts) >= 2:
                    return _clean_acc(parts[1])
    return None


def _meta_is_meaningful(meta: str) -> bool:
    if not meta.strip():
        return False
    lines = [ln for ln in meta.splitlines() if ln.strip()]
    if not lines or not lines[0].startswith("LOCUS"):
        return False
    useful = ("DEFINITION", "ACCESSION", "VERSION", "KEYWORDS",
              "SOURCE", "REFERENCE", "COMMENT", "DBLINK", "ORGANISM")
    return any(any(ln.startswith(p) for p in useful) for ln in lines[1:])


def build_metadata_index(src_text: str) -> Tuple[Dict, Dict, Dict]:
    exact_dict: Dict[str, str] = {}
    clean_dict: Dict[str, str] = {}
    norm_dict:  Dict[str, Tuple[str, str]] = {}
    for rec in _split_gb_records(src_text):
        raw_acc = _get_record_accession(rec)
        if not raw_acc:
            continue
        meta, _ = _split_meta_body(rec)
        if not _meta_is_meaningful(meta):
            continue
        for v in _acc_variants(raw_acc):
            exact_dict[v] = meta
        cl = _clean_acc(raw_acc)
        if cl:
            clean_dict.setdefault(cl, meta)
        norm = _normalize_acc(raw_acc)
        if norm:
            norm_dict.setdefault(norm, (raw_acc, meta))
    return exact_dict, clean_dict, norm_dict


def find_metadata_from_index(acc: str, exact_dict, clean_dict, norm_dict,
                              ) -> Tuple[Optional[str], Optional[str]]:
    for v in _acc_variants(acc):
        if v in exact_dict:
            mtype = "exact" if v == acc else f"clean '{acc}'->'{v}'"
            return exact_dict[v], mtype
    cl = _clean_acc(acc)
    if cl and cl in clean_dict:
        return clean_dict[cl], f"clean '{acc}'->'{cl}'"
    norm = _normalize_acc(acc)
    if norm and norm in norm_dict:
        orig, meta = norm_dict[norm]
        return meta, f"compat '{orig}'<->'{acc}'"
    return None, None


def _similarity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    return difflib.SequenceMatcher(None, a.upper(), b.upper()).ratio()


def diagnose_unmatched(acc: str, exact_dict: Dict) -> Tuple[List[str], List[Tuple[int, str]]]:
    hints: List[str] = []
    candidates: List[Tuple[int, str]] = []
    acc_cl   = _clean_acc(acc)     or acc
    acc_norm = _normalize_acc(acc) or acc
    acc_up   = acc_cl.upper()
    seen: set = set()
    for k in exact_dict:
        ck = _clean_acc(k)
        if not ck or ck in seen:
            continue
        seen.add(ck)
        orig_norm = _normalize_acc(ck) or ck
        score = 0
        reason: Optional[str] = None
        if ck.upper() == acc_up and ck != acc_cl:
            score = 90
            reason = f"[case mismatch] orig='{k}' new='{acc}'"
        elif orig_norm == acc_norm and ck != acc_cl:
            score = 88
            reason = f"[version/style mismatch] orig='{k}' new='{acc}'"
        else:
            sim = _similarity(acc_norm, orig_norm)
            if sim > 0.78:
                score = int(sim * 55)
        if score > 0:
            if reason:
                hints.append(reason)
            candidates.append((score, k))
    if not hints:
        hints.append(f"[no diagnosis] no similar accession found for '{acc}'")
    candidates.sort(key=lambda x: -x[0])
    return hints[:5], candidates[:5]


def interactive_resolve(acc: str, candidates: List[Tuple[int, str]],
                        exact_dict: Dict, auto_accept: bool = False,
                        ) -> Tuple[Optional[str], Optional[str]]:
    if not candidates:
        return None, None
    best_score, best_acc = candidates[0]
    if auto_accept:
        if best_score >= 60:
            logging.info(f"[AUTO-ACCEPT] '{acc}' -> '{best_acc}' (score={best_score})")
            return exact_dict.get(best_acc), best_acc
        logging.warning(f"[AUTO-SKIP] '{acc}' top score {best_score} < 60")
        return None, None
    if not sys.stdin.isatty():
        logging.warning(
            f"[NON-INTERACTIVE] stdin is not a TTY; skipping interactive resolution "
            f"for '{acc}'. Top candidate: '{best_acc}' (score={best_score}). "
            f"Use --auto_accept_metadata yes to resolve automatically.")
        return None, None
    print(f"\n⚠ '{acc}' not matched. Candidates:")
    for i, (score, orig_acc) in enumerate(candidates, start=1):
        print(f"  [{i}] {orig_acc:<30} score={score}")
    print("  [s] skip   [q] quit")
    while True:
        try:
            choice = input("choose: ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            print()
            return None, None
        if choice == "q":
            sys.exit(0)
        if choice == "s":
            return None, None
        if choice.isdigit():
            idx = int(choice) - 1
            if 0 <= idx < len(candidates):
                chosen = candidates[idx][1]
                return exact_dict.get(chosen), chosen
        print("  invalid, retry.")


def replace_metadata_in_text(
    new_gb_text: str, exact_dict, clean_dict, norm_dict,
    smart: bool = False, auto_accept: bool = False,
) -> Tuple[str, Dict]:
    """Replace GB header sections in *new_gb_text* using the metadata index."""
    stats: Dict = {
        "total_new": 0, "no_accession": 0,
        "exact": 0, "compat": 0, "smart_resolved": 0, "not_found": 0,
        "exact_list": [], "compat_list": [], "smart_list": [], "not_found_list": [],
    }
    out_records: List[str] = []
    for rec in _split_gb_records(new_gb_text):
        stats["total_new"] += 1
        acc = _get_record_accession(rec)
        if not acc:
            stats["no_accession"] += 1
            out_records.append(rec)
            continue
        meta, match_type = find_metadata_from_index(acc, exact_dict, clean_dict, norm_dict)
        if meta is not None:
            _, body = _split_meta_body(rec)
            # the bp count of the source metadata's LOCUS line: the length of the re-annotated sequence
            new_locus_m = re.search(r"LOCUS\s+\S+\s+(\d+)\s+[a-zA-Z]+", rec)
            if new_locus_m:
                real_len = new_locus_m.group(1)
                meta = _LOCUS_BPCOUNT_RE.sub(rf"\g<1>{real_len}\g<2>", meta, count=1)
            out_records.append(meta + body if body else rec)
            if match_type == "exact":
                stats["exact"] += 1
                stats["exact_list"].append(acc)
            else:
                stats["compat"] += 1
                stats["compat_list"].append((acc, match_type))
            continue
        hints, candidates = diagnose_unmatched(acc, exact_dict)
        if smart and candidates:
            resolved_meta, chosen_acc = interactive_resolve(
                acc, candidates, exact_dict, auto_accept=auto_accept)
            if resolved_meta is not None:
                _, body = _split_meta_body(rec)
                out_records.append(resolved_meta + body if body else rec)
                stats["smart_resolved"] += 1
                stats["smart_list"].append((acc, chosen_acc))
                continue
        stats["not_found"] += 1
        stats["not_found_list"].append((acc, hints))
        out_records.append(rec)
    return "".join(out_records), stats


def transfer_gb_metadata_file(
    metadata_src_gb: Path, new_body_gb: Path, output_gb: Path,
    smart: bool = False, auto_accept: bool = False,
) -> Tuple[bool, int, int, List[str], str, List[str]]:
    if not metadata_src_gb.exists():
        return False, 0, 0, [], f"metadata source not found: {metadata_src_gb}", []
    if not new_body_gb.exists():
        return False, 0, 0, [], f"body source not found: {new_body_gb}", []
    src_text = metadata_src_gb.read_text(encoding="utf-8", errors="replace")
    new_text = new_body_gb.read_text(encoding="utf-8", errors="replace")
    exact_dict, clean_dict, norm_dict = build_metadata_index(src_text)
    if not exact_dict and not norm_dict:
        return False, 0, 0, [], "metadata source has no usable metadata", []
    replaced_text, stats = replace_metadata_in_text(
        new_text, exact_dict, clean_dict, norm_dict, smart=smart, auto_accept=auto_accept)
    output_gb.parent.mkdir(parents=True, exist_ok=True)
    output_gb.write_text(replaced_text, encoding="utf-8")
    n_replaced  = stats["exact"] + stats["compat"] + stats["smart_resolved"]
    unmatched   = [x[0] for x in stats["not_found_list"]]
    compat_keys = [x[0] for x in stats["compat_list"]]
    return True, n_replaced, stats["not_found"], unmatched, "ok", compat_keys


# ─────────────────────────────────────────────────────────────────────────────
# Organism / SOURCE replacement helpers
# ─────────────────────────────────────────────────────────────────────────────
def build_organism_index(gb_files: List[Path]) -> Dict[str, Tuple[str, str]]:
    """Accession -> (organism_name, source_line_text) lookup, under every accession variant (_acc_variants, which
    include the normalised form)."""
    index: Dict[str, Tuple[str, str]] = {}
    for gb_path in gb_files:
        try:
            records = list(SeqIO.parse(str(gb_path), "genbank"))
        except Exception as exc:
            logging.warning(f"[ORG_INDEX] Cannot parse {gb_path}: {exc}")
            continue
        for rec in records:
            acc = gb_record_accession(rec)

            organism: str = rec.annotations.get("organism", "").strip()
            if not organism:
                for feat in rec.features:
                    if feat.type == "source":
                        vals = feat.qualifiers.get("organism", [])
                        if vals:
                            organism = str(vals[0]).strip()
                        break
            if not organism:
                organism = rec.annotations.get("source", "").strip()

            source_text: str = rec.annotations.get("source", "").strip()

            if not organism and not source_text:
                continue

            entry: Tuple[str, str] = (organism, source_text)

            for variant in _acc_variants(acc):
                if variant in index and index[variant] != entry:      # first wins, the conflict is logged
                    logging.warning(
                        f"[ORG_INDEX_CONFLICT] {variant}: keeping first "
                        f"{index[variant]!r}, ignoring {entry!r} from {gb_path.name}")
                else:
                    index.setdefault(variant, entry)

    return index


def lookup_organism(
    acc: str, org_index: Dict[str, Tuple[str, str]]
) -> Optional[Tuple[str, str]]:
    """Return (organism, source_text) for *acc* using exact → normalised fallback."""
    for v in _acc_variants(acc):
        if v in org_index:
            return org_index[v]
    norm = _normalize_acc(acc)
    if norm and norm in org_index:
        return org_index[norm]
    return None


def replace_organism_in_record(
    rec: "SeqRecord",
    new_organism: str,
    new_source: Optional[str] = None,
) -> "SeqRecord":
    """A copy of *rec* with the organism name replaced throughout: rec.annotations["organism"] when new_organism is
    non-empty, rec.annotations["source"] only when new_source is given and non-empty, /organism= on every source
    feature (added when absent)."""
    new_annotations = dict(rec.annotations)
    if new_organism:
        new_annotations["organism"] = new_organism
    if new_source:
        new_annotations["source"] = new_source

    new_features: List["SeqFeature"] = []
    for feat in rec.features:
        new_quals = dict(feat.qualifiers)
        if feat.type == "source" and new_organism:
            new_quals["organism"] = [new_organism]
        new_features.append(
            SeqFeature(
                location=feat.location,
                type=feat.type,
                id=getattr(feat, "id", "<unknown>"),
                qualifiers=new_quals,
            )
        )

    new_rec = SeqRecord(
        rec.seq,
        id=rec.id, name=rec.name, description=rec.description,
        dbxrefs=list(getattr(rec, "dbxrefs", [])),
        annotations=new_annotations,
        features=new_features,
    )
    return new_rec


def replace_organism_in_gb_file(
    gb_in:     Path,
    gb_out:    Path,
    org_index: Dict[str, Tuple[str, str]],
) -> Tuple[int, int, int]:
    """Replace /organism= (and optionally SOURCE) in every record of *gb_in*."""
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    if not recs:
        raise ValueError(f"No records in {gb_in}")
    out_recs: List["SeqRecord"] = []
    n_replaced = n_not_found = 0
    for rec in recs:
        acc    = gb_record_accession(rec)
        result = lookup_organism(acc, org_index)
        if result is not None:
            organism, source_text = result
            out_recs.append(
                replace_organism_in_record(rec, organism,
                                           new_source=source_text if source_text else None)
            )
            n_replaced += 1
            logging.info(
                f"[REPLACE_ORG] {acc}: /organism= → '{organism}'"
                + (f"  SOURCE → '{source_text}'" if source_text else "")
            )
        else:
            out_recs.append(rec)
            n_not_found += 1
            logging.warning(
                f"[REPLACE_ORG_MISS] {acc}: no entry in organism index — unchanged.")
    gb_out.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_recs, str(gb_out), "genbank")
    return n_replaced, n_not_found, len(recs)


def _sync_organism_qualifier(gb_in: Path, gb_out: Path) -> int:
    """Sync /organism= in FEATURES source to match rec.annotations["organism"] (added when absent)."""
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    out_recs: List["SeqRecord"] = []
    n_fixed = 0
    for rec in recs:
        organism = rec.annotations.get("organism", "").strip()
        if not organism:
            out_recs.append(rec)
            continue
        new_features: List["SeqFeature"] = []
        changed = False
        for feat in rec.features:
            new_quals = dict(feat.qualifiers)
            if feat.type == "source":
                current = str(new_quals["organism"][0]).strip() \
                          if new_quals.get("organism") else ""
                if current != organism:
                    new_quals["organism"] = [organism]
                    changed = True
            new_features.append(
                SeqFeature(
                    location=feat.location,
                    type=feat.type,
                    id=getattr(feat, "id", "<unknown>"),
                    qualifiers=new_quals,
                )
            )
        new_rec = SeqRecord(
            rec.seq,
            id=rec.id, name=rec.name, description=rec.description,
            dbxrefs=list(getattr(rec, "dbxrefs", [])),
            annotations=dict(rec.annotations),
            features=new_features,
        )
        out_recs.append(new_rec)
        if changed:
            n_fixed += 1
    gb_out.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_recs, str(gb_out), "genbank")
    return n_fixed


def gb_file_has_meaningful_metadata(gb_path: Path) -> bool:
    if not gb_path.is_file():
        return False
    try:
        txt = gb_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return False
    return any(_meta_is_meaningful(_split_meta_body(r)[0]) for r in _split_gb_records(txt))


# ─────────────────────────────────────────────────────────────────────────────
# Environment preflight
# ─────────────────────────────────────────────────────────────────────────────
def check_environment(args: argparse.Namespace) -> List[str]:
    problems: List[str] = []
    if SeqIO is None:
        raise RuntimeError("Biopython not installed. Run: pip install biopython")
    try:
        import Bio
        logging.info(f"[ENV] Biopython {getattr(Bio, '__version__', '?')} ✔")
    except Exception:
        pass
    major, minor = sys.version_info.major, sys.version_info.minor
    if (major, minor) < (3, 8):
        raise RuntimeError(f"Python 3.8+ required, found {major}.{minor}")
    logging.info(f"[ENV] Python {major}.{minor} ✔")

    if getattr(args, "cmd", None) in ("annotate", "run"):
        import mito_env
        for prog in ("blastx", "exonerate", "cmsearch", "java"):
            found = mito_env.find(prog, None)
            if found:
                logging.info(f"[ENV] {prog} ✔  {found}")
            else:
                problems.append(f"[ENV] MitoM needs {prog} (install.sh; mito_env.py doctor)")
        if not mito_env.mitos_refdir(getattr(args, "mitos_refdir", None)):
            problems.append("[ENV] MitoM needs the MITOS2 reference data refseq89m (--mitos_refdir; install.sh)")

    smk_bin = getattr(args, "snakemake_bin", None)
    backend = getattr(args, "scheduler_backend", "threads")
    if backend in ("snakemake", "auto"):
        if snakemake_available(smk_bin):
            smk_exe = smk_bin or which_or_none("snakemake") or "snakemake"
            try:
                sv = subprocess.run([smk_exe, "--version"], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, text=True, timeout=10)
                smk_ver = (sv.stdout or sv.stderr or "").strip().splitlines()[0] \
                          if sv.returncode == 0 else "?"
                logging.info(f"[ENV] Snakemake {smk_ver} ✔")
            except Exception:
                logging.info("[ENV] Snakemake available (version check failed)")
        else:
            if backend == "snakemake":
                raise RuntimeError("[ENV] --scheduler_backend snakemake requested but snakemake not found.")
            logging.info("[ENV] Snakemake not found, will use threads backend.")

    if hasattr(args, "input") and args.input:
        in_p = Path(args.input)
        if not in_p.exists():
            raise RuntimeError(f"[ENV] Input path does not exist: {in_p}")
        logging.info(f"[ENV] Input path exists ✔  ({in_p})")

    if hasattr(args, "infq") and args.infq:
        fq_p = Path(args.infq)
        if not fq_p.exists():
            problems.append(f"[ENV] --infq path does not exist: {fq_p}")
        else:
            logging.info(f"[ENV] FASTQ dir exists ✔  ({fq_p})")

    return problems


def preflight(args: argparse.Namespace) -> None:
    problems = check_environment(args)
    if problems:
        for p in problems:
            logging.warning(p)
        logging.warning(f"[ENV] {len(problems)} warning(s) above — pipeline continues but may fail.")


# ─────────────────────────────────────────────────────────────────────────────
# Annotation worker (mito_core.py)
# ─────────────────────────────────────────────────────────────────────────────
def sample_lineage(si: "SampleInfo") -> str:
    """NCBI lineage of a sample from its original GenBank record ('' for FASTA input): the MitoM-ML tRNA caller of
    mito_core.py uses the clade models on this lineage (full lineage 98.5 % vs phylum only 97.9 % tRNA recall,
    projects/invert_trna)."""
    if not si.orig_gb_path:
        return ""
    from Bio import SeqIO
    try:
        for rec in SeqIO.parse(str(si.orig_gb_path), "genbank"):
            tax = [t for t in rec.annotations.get("taxonomy", []) if t]
            if tax:
                return ";".join(tax)
    except (OSError, ValueError):
        pass
    return ""


def native_core_cmd(fasta_file: Path, workdir: Path, internal_id: str, args, lineage: str = "") -> List[str]:
    """mito_core.py annotate command for one sample (--core native); `args` is the Namespace or the worker's dict."""
    g = (lambda k, d=None: getattr(args, k, d)) if not isinstance(args, dict) else (lambda k, d=None: args.get(k, d))
    cmd = [sys.executable, "-W", "ignore", str(MITO_CORE), "annotate", "-i", str(fasta_file), "-o", str(workdir),
           "--prefix", internal_id, "--genetic_code", str(int(g("genetic_code", 5))), "--threads", str(int(g("threads", 1)))]
    if g("mitos_refdir"):
        cmd += ["--refdir", str(g("mitos_refdir"))]
    if g("mitos_bin"):
        cmd += ["--mitos_bin", str(g("mitos_bin"))]
    if g("phylum"):
        cmd += ["--phylum", str(g("phylum"))]
    if g("rrna_hmm_lib"):
        cmd += ["--rrna_lib", str(g("rrna_hmm_lib"))]
    if g("rrna_method"):
        cmd += ["--rrna_method", str(g("rrna_method"))]
    if lineage:
        cmd += ["--lineage", lineage]
    return cmd


def engine_cmd(fasta_file: Path, workdir: Path, internal_id: str, args, lineage: str = "") -> List[str]:
    """Command that writes <workdir>/<internal_id>.gbf with the annotation of --engine: MitoM's core (native), or
    MitoZ / MITOS2 run by mito_engines.py; `args` is the Namespace or a dict."""
    g = (lambda k, d=None: getattr(args, k, d)) if not isinstance(args, dict) else (lambda k, d=None: args.get(k, d))
    engine = g("engine") or "native"
    if engine == "native":
        return native_core_cmd(fasta_file, workdir, internal_id, args, lineage)
    cmd = [sys.executable, "-W", "ignore", str(MITO_ENGINES), engine, "-i", str(fasta_file), "-o", str(workdir),
           "--prefix", internal_id, "--genetic_code", str(int(g("genetic_code", 5))), "--threads", str(int(g("threads", 1)))]
    if engine == "mitoz":
        cmd += ["--clade", str(g("clade"))]
        cmd += ["--mitoz_bin", str(g("mitoz_bin"))] if g("mitoz_bin") else []
        cmd += ["--species_name", str(g("species_name"))] if g("species_name") else []
    else:
        cmd += ["--refdir", str(g("mitos_refdir"))] if g("mitos_refdir") else []
        cmd += ["--mitos_bin", str(g("mitos_bin"))] if g("mitos_bin") else []
    return cmd


def write_input_annotation(si: SampleInfo, dest: Path) -> str:
    """`mitom import`: the input GenBank record of this sample (annotated by MitoZ, MITOS2 or another program) as its
    round-1 annotation, named like the core's record. Returns an error message, '' when written."""
    if not si.orig_gb_path:
        return f"{si.sample_name}: no input GenBank record"
    recs = list(SeqIO.parse(str(si.orig_gb_path), "genbank"))
    rec = next((r for r in recs if si.orig_token in (r.id, r.name, safe_name(r.id))), recs[0] if len(recs) == 1 else None)
    if rec is None:
        return f"{si.sample_name}: record {si.orig_token} not found in {si.orig_gb_path}"
    if not any(f.type in ("CDS", "rRNA", "tRNA") for f in rec.features):
        return f"{si.sample_name}: the input record has no gene annotation"
    out = SeqRecord(rec.seq, id=si.internal_id, name=si.internal_id[:16], description="", features=list(rec.features))
    out.annotations.update(molecule_type="DNA", topology=si.topology,
                           data_file_division=rec.annotations.get("data_file_division", "INV"))
    SeqIO.write(out, str(dest), "genbank")
    return ""


def annotate_one(
    si: SampleInfo, fasta_file: Path, out_root: Path, phase_dirname: str,
    args: argparse.Namespace, total: int = 0, done_counter: Optional[List] = None,
) -> Tuple[str, bool, Path, Optional[Path]]:
    """Annotate one sample (engine_cmd, or the input record for `mitom import`): (internal_id, ok, workdir, the
    <workdir>/<internal_id>.gbf the engine wrote or None)."""
    internal_id = si.internal_id
    workdir     = out_root / phase_dirname / internal_id
    workdir.mkdir(parents=True, exist_ok=True)
    log_name = f"sample.{phase_dirname}.{internal_id}"
    lg = get_sample_logger(log_name, workdir / "sample.log")
    lg.info(f"[SAMPLE] internal_id={internal_id} sample_name={si.sample_name}")
    lg.info(f"[INPUT_FASTA] {fasta_file}")
    logs_dir    = workdir / "logs"
    stdout_path = logs_dir / "core.stdout.log"
    stderr_path = logs_dir / "core.stderr.log"
    if getattr(args, "engine", "native") == "input":          # mitom import: no annotation program runs
        t0 = time.time()
        logs_dir.mkdir(parents=True, exist_ok=True)
        err = write_input_annotation(si, workdir / f"{internal_id}.gbf")
        stderr_path.write_text(err + "\n" if err else "", encoding="utf-8")
        rc, dt = (1 if err else 0), time.time() - t0
    else:
        cmd = engine_cmd(fasta_file, workdir, internal_id, args, sample_lineage(si))
        lg.info("[CMD] " + " ".join(cmd))
        rc, dt = run_cmd_capture(cmd, cwd=workdir, stdout_path=stdout_path, stderr_path=stderr_path)
    if rc != 0:
        (logs_dir / "core.returncode.txt").write_text(str(rc) + "\n", encoding="utf-8")
    (logs_dir / "core.runtime_sec.txt").write_text(f"{dt:.2f}\n", encoding="utf-8")

    if done_counter is not None:
        lock: threading.Lock = done_counter[1]
        with lock:
            done_counter[0] += 1
            done_n = done_counter[0]
        logging.info(f"[PROGRESS] {done_n}/{total} '{si.sample_name}' done (rc={rc}, {dt:.0f}s)")
        if _STAGE_PROGRESS is not None:
            _STAGE_PROGRESS.update(done_n)

    if rc != 0:
        err_tail = tail_text(stderr_path.read_text(encoding="utf-8", errors="ignore") if stderr_path.exists() else "")
        lg.error(f"[FAIL] annotation core returned {rc}")
        lg.error("---- stderr tail ----\n" + err_tail)
        last = next((x.strip() for x in reversed(err_tail.splitlines()) if x.strip()), "no stderr")
        logging.warning(f"[FAIL] {si.sample_name}: annotation core exit {rc}: {last}  (sample log: {workdir / 'sample.log'})")
        close_sample_logger(log_name)
        return internal_id, False, workdir, None

    own = workdir / f"{internal_id}.gbf"                 # what the core and mito_engines.py write (exit 0 = written)
    close_sample_logger(log_name)
    return internal_id, True, workdir, own if own.is_file() else None


# ─────────────────────────────────────────────────────────────────────────────
# Snakemake scheduler
# ─────────────────────────────────────────────────────────────────────────────
def snakemake_available(snakemake_bin: Optional[str]) -> bool:
    if snakemake_bin:
        return Path(snakemake_bin).exists() or which_or_none(snakemake_bin) is not None
    return which_or_none("snakemake") is not None


@dataclass
class AnnotTask:
    internal_id:   str
    sample_name:   str
    fasta_file:    str
    phase_dirname: str
    out_root:      str
    lineage:       str = ""
    cmd:           List[str] = field(default_factory=list)      # engine_cmd(): what the worker runs


def render_snakemake_worker_script(out_dir: Path) -> Path:
    """Generate the self-contained annotation worker script (runs the task's engine command, engine_cmd()). The
    worker annotates only; finalize_genbank_output and every later step run in the main process."""
    fp   = out_dir / "_worker.py"
    code = r'''#!/usr/bin/env python3
"""Self-contained annotation worker invoked by Snakemake: runs the task's command (the annotation engine).

This script handles annotation only.  Finalization steps
(ID-rewrite, reorientation, metadata overlay, organism sync) are
executed in the main pipeline process after all workers finish.
"""
import json, subprocess, sys, time
from pathlib import Path

def tail_text(s, n=120): return "\n".join((s or "").splitlines()[-n:])

def main():
    task = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    internal_id   = task["internal_id"]
    phase_dirname = task["phase_dirname"]
    out_root      = Path(task["out_root"]).resolve()
    fasta_file    = Path(task["fasta_file"]).resolve()
    workdir       = out_root / phase_dirname / internal_id
    workdir.mkdir(parents=True, exist_ok=True)
    logs_dir = workdir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    cmd = task["cmd"]                               # engine_cmd() of the main process

    t0 = time.time()
    p  = subprocess.run(cmd, cwd=str(workdir),
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    dt = time.time() - t0

    if p.stdout:
        (logs_dir / "core.stdout.log").write_text(
            p.stdout, encoding="utf-8", errors="ignore")
    if p.stderr:
        (logs_dir / "core.stderr.log").write_text(
            p.stderr, encoding="utf-8", errors="ignore")
    if p.returncode != 0:
        (logs_dir / "core.returncode.txt").write_text(
            str(p.returncode) + "\n", encoding="utf-8")
    (logs_dir / "core.runtime_sec.txt").write_text(
        f"{dt:.2f}\n", encoding="utf-8")

    result = {
        "internal_id": internal_id,
        "sample_name": task.get("sample_name", internal_id),
        "ok": False, "workdir": str(workdir), "best_gb": None,
        "returncode": int(p.returncode),
    }
    if p.returncode == 0:
        own = workdir / (internal_id + ".gbf")          # what the engine writes (exit 0 = written)
        result.update({"ok": own.is_file(), "best_gb": str(own) if own.is_file() else None})
    else:
        result["stderr_tail"] = tail_text(p.stderr or "")

    Path(sys.argv[1]).with_suffix(".result.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    Path(sys.argv[1]).with_suffix(".done").write_text("done\n", encoding="utf-8")

if __name__ == "__main__":
    main()
'''
    fp.write_text(code, encoding="utf-8")
    return fp


def render_snakemake_snakefile(
    workdir: Path, task_jsons: List[Path], worker_script: Path,
    python_exe: str = sys.executable or "python3",
    job_threads: int = 1,
    job_mem_mb: int = 128,
) -> Path:
    snakefile = workdir / "Snakefile"
    rule_lines: List[str] = []
    all_done:   List[str] = []
    for tp in task_jsons:
        done_s  = str(tp.with_suffix(".done"))
        rule_id = re.sub(r"[^A-Za-z0-9_]", "_", tp.stem)
        all_done.append(done_s)
        cmdline = " ".join(shlex.quote(x) for x in [python_exe, str(worker_script), str(tp)])
        # threads / mem_mb per job: what the engine gets (--threads), so that Snakemake schedules the real load
        rule_lines.append(f"""
rule task_{rule_id}:
    input:  task={json.dumps(str(tp), ensure_ascii=False)}
    output: done={json.dumps(done_s, ensure_ascii=False)}
    threads: {max(1, int(job_threads))}
    resources: mem_mb={int_mem_mb(job_mem_mb)}
    shell:  {json.dumps(cmdline, ensure_ascii=False)}
""")
    header = f"rule all:\n    input: {json.dumps(all_done, ensure_ascii=False)}\n"
    snakefile.write_text(header + "\n".join(rule_lines), encoding="utf-8")
    return snakefile


def run_snakemake_for_annotation_tasks(
    *, tasks: List[AnnotTask], out_root: Path, phase_name: str, args: argparse.Namespace,
) -> Dict[str, Dict]:
    snk_work = out_root / f"_snakemake_{phase_name}"
    snk_work.mkdir(parents=True, exist_ok=True)
    worker_script = render_snakemake_worker_script(snk_work)
    task_jsons: List[Path] = []
    for t in tasks:
        tp = snk_work / f"{t.internal_id}.json"
        tp.write_text(json.dumps(asdict(t), indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
        task_jsons.append(tp)
    per_job_mem_mb = max(128, int_mem_mb(args.smk_mem_mb) // max(1, int(args.max_tasks)))
    snakefile = render_snakemake_snakefile(
        snk_work, task_jsons, worker_script,
        python_exe=sys.executable or "python3",
        job_threads=int(args.threads),
        job_mem_mb=per_job_mem_mb)
    smk_bin   = args.snakemake_bin or which_or_none("snakemake")
    if not smk_bin:
        raise FileNotFoundError("snakemake executable not found")
    cmd = [
        str(smk_bin), "--snakefile", str(snakefile), "--directory", str(snk_work),
        "--cores", str(max(1, int(args.smk_cores))),
        "--jobs",  str(max(1, int(args.max_tasks))),
        "--resources", f"mem_mb={int_mem_mb(args.smk_mem_mb)}",
        "--rerun-incomplete", "--keep-going", "--printshellcmds",
        "--latency-wait", str(int(args.smk_latency_wait)),
    ]
    smk_extra: List[str] = []
    for token in (getattr(args, "smk_extra", None) or []):
        smk_extra.extend(shlex.split(str(token)))
    if smk_extra:
        logging.info("[SNAKEMAKE_EXTRA] " + " ".join(smk_extra))
        cmd.extend(smk_extra)

    def _run(c: List[str]) -> "subprocess.CompletedProcess[str]":
        logging.info("[SNAKEMAKE] " + " ".join(c))
        return subprocess.run(c, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    global _STAGE_PROGRESS
    stop = threading.Event()
    if _STAGE_PROGRESS is not None:
        watch_files(snk_work, "*.result.json", _STAGE_PROGRESS, stop)
    try:
        p = _run(cmd)
        if p.returncode != 0:
            combined = (p.stdout or "") + (p.stderr or "")
            if "LockException" in combined or "Directory cannot be locked" in combined:
                logging.warning(f"[SNAKEMAKE_LOCK] Stale lock in {snk_work} — attempting --unlock…")
                pu = subprocess.run(
                    [str(smk_bin), "--snakefile", str(snakefile),
                     "--directory", str(snk_work), "--unlock"],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if pu.returncode == 0:
                    logging.warning("[SNAKEMAKE_LOCK] Unlock OK — retrying…")
                    p = _run(cmd)
                else:
                    logging.error(f"[SNAKEMAKE_UNLOCK_FAIL] rc={pu.returncode}")
    finally:
        stop.set()
    if _STAGE_PROGRESS is not None:                            # the watcher polls every 10 s: print the final count
        _STAGE_PROGRESS.update(sum(1 for tp in task_jsons if tp.with_suffix(".result.json").exists()))

    (snk_work / "snakemake.stdout.log").write_text(p.stdout or "", encoding="utf-8",
                                                    errors="ignore")
    (snk_work / "snakemake.stderr.log").write_text(p.stderr or "", encoding="utf-8",
                                                    errors="ignore")
    if p.returncode != 0:
        raise RuntimeError(
            f"Snakemake failed for phase '{phase_name}'\n"
            f"stdout:\n{tail_text(p.stdout or '')}\n\nstderr:\n{tail_text(p.stderr or '')}")

    results: Dict[str, Dict] = {}
    for tp in task_jsons:
        rp = tp.with_suffix(".result.json")
        if rp.exists():
            obj = json.loads(rp.read_text(encoding="utf-8"))
            results[obj["internal_id"]] = obj
            if not obj.get("ok"):                              # the worker's sample.log holds the details
                last = next((x.strip() for x in reversed(str(obj.get("stderr_tail") or "").splitlines())
                             if x.strip()), "no stderr")
                logging.warning(f"[FAIL] {obj.get('sample_name', obj['internal_id'])}: annotation core exit "
                                f"{obj.get('returncode', '?')}: {last}  (sample log: {obj.get('workdir')}/sample.log)")
        else:
            iid = tp.stem
            logging.warning(f"[SNAKEMAKE_RESULT_MISSING] {rp} not found for task {iid}")
            results[iid] = {
                "internal_id": iid, "sample_name": iid, "ok": False,
                "workdir": str(out_root / phase_name / iid), "best_gb": None,
            }
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Unified annotation stage dispatcher
# ─────────────────────────────────────────────────────────────────────────────
def run_annotation_stage(
    *, samples: List[SampleInfo], phase_dirname: str, out_root: Path,
    args: argparse.Namespace, fasta_selector,
) -> Dict[str, Dict]:
    backend = getattr(args, "scheduler_backend", "threads")
    if getattr(args, "engine", "native") == "input":            # mitom import: nothing to run
        backend = "threads"
    if backend == "auto":
        # Snakemake schedules nothing when every sample runs at once, and its process costs ~120 MB (a third of the
        # peak memory of a one-record run): threads then
        backend = "snakemake" if snakemake_available(getattr(args, "snakemake_bin", None)) \
                  and len(samples) > int(getattr(args, "max_tasks", 1) or 1) else "threads"
    logging.info(f"[ANNOT_STAGE] phase={phase_dirname}  backend={backend}  n={len(samples)}")
    global _STAGE_PROGRESS
    _STAGE_PROGRESS = Progress(STAGE_LABEL.get(phase_dirname, phase_dirname), len(samples))

    if backend == "snakemake":
        tasks: List[AnnotTask] = []
        for si in samples:
            tasks.append(AnnotTask(
                internal_id=si.internal_id, sample_name=si.sample_name,
                fasta_file=str(fasta_selector(si)), phase_dirname=phase_dirname,
                out_root=str(out_root), lineage=sample_lineage(si),
                cmd=engine_cmd(Path(fasta_selector(si)).resolve(), (out_root / phase_dirname / si.internal_id).resolve(),
                               si.internal_id, args, sample_lineage(si))))
        return run_snakemake_for_annotation_tasks(
            tasks=tasks, out_root=out_root, phase_name=phase_dirname, args=args)

    results: Dict[str, Dict] = {}
    total = len(samples)
    done_counter: List = [0, threading.Lock()]

    def _task(si: SampleInfo) -> Tuple:
        return annotate_one(si, fasta_selector(si), out_root, phase_dirname, args,
                            total=total, done_counter=done_counter)

    with cf.ThreadPoolExecutor(max_workers=int(args.max_tasks)) as ex:
        futs = {ex.submit(_task, si): si for si in samples}
        for fut in cf.as_completed(futs):
            si = futs[fut]
            try:
                internal_id, ok, workdir, best_gb = fut.result()
            except Exception as exc:
                logging.error(f"[TASK_EXCEPTION] {si.sample_name}: {exc}")
                internal_id, ok, workdir, best_gb = si.internal_id, False, out_root / phase_dirname / si.internal_id, None
            results[internal_id] = {
                "internal_id": internal_id, "sample_name": si.sample_name,
                "ok": ok, "workdir": str(workdir), "best_gb": str(best_gb) if best_gb else None,
            }
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Final GenBank output
# ─────────────────────────────────────────────────────────────────────────────
def starts_with_gene(gb: Path, which_gene_first: str, force_forward: bool) -> bool:
    """Does a record of `gb` begin with its orientation gene (the first of which_gene_first present; on the plus
    strand when force_forward), i.e. no gene feature lies before it?"""
    for rec in SeqIO.parse(str(gb), "genbank"):
        feat, _ = find_orientation_feature(rec, which_gene_first)
        if feat is None or (force_forward and feat.location.strand == -1):
            continue
        genes = [f for f in rec.features if f.type in ("CDS", "rRNA", "tRNA")]
        if genes and feature_start_pos(feat) <= min(feature_start_pos(f) for f in genes):
            return True
    return False


def finalize_genbank_output(
    *,
    sample_name:       str,
    body_src_gb:       Path,
    metadata_src_gb:   Optional[Path],
    dest_gb:           Path,
    which_gene_first:  str,
    force_forward:     bool,
    preserve_metadata: bool,
    smart_metadata:    bool = False,
    auto_accept:       bool = False,
) -> Tuple[int, int, int, str, bool]:
    """Produce the final, publication-grade GenBank file for one sample.

    Steps
    -----
    1. ID rewrite      — replace internal S0xxxx IDs with sample_name.
    2. Reorientation   — guarantee the record starts at which_gene_first.
    3. Metadata overlay — (when preserve_metadata) replace GB header section.
    3.5 Organism sync  — fix /organism= inside the FEATURES source feature.
    4. Recount         — recount PCG / tRNA / rRNA from the definitive file.

    Every step writes a draft next to dest_gb that replaces it only at the end: a failure never leaves a half-final
    file at the final path.
    """
    pcg = trna = rrna = 0
    orientation_ok = False

    td = Path(tempfile.mkdtemp(prefix="mitom_final_"))
    dest_gb.parent.mkdir(parents=True, exist_ok=True)
    draft_gb = dest_gb.with_name(f".{dest_gb.name}.tmp")
    committed = False
    try:
        if draft_gb.exists():
            draft_gb.unlink()

        # ── Step 1: ID rewrite ────────────────────────────────────────
        id_gb = td / "id_rewrite.gbf"
        try:
            pcg, trna, rrna = rewrite_genbank_record_ids(body_src_gb, id_gb, sample_name)
        except Exception as exc:
            raise RuntimeError(f"ID rewrite failed for {sample_name}: {exc}") from exc

        # ── Step 2: Final guaranteed reorientation ────────────────────
        reor_gb    = td / "reoriented.gbf"
        no_gene_gb = td / "no_gene.gbf"
        try:
            n_ok, n_no, _ = reorient_one_gb(
                id_gb, reor_gb, no_gene_gb, which_gene_first, force_forward)
            if n_ok > 0 and reor_gb.exists():
                orientation_ok = True
                if n_no > 0 and no_gene_gb.exists():
                    merged_reor = td / "merged_reoriented.gbf"
                    merge_gb_files([reor_gb, no_gene_gb], merged_reor)
                    working_gb = merged_reor
                    logging.info(
                        f"[FINAL_MULTICONTIG] {sample_name}: merged {n_ok} reoriented "
                        f"+ {n_no} un-reoriented contig(s) into final GB.")
                else:
                    working_gb = reor_gb
            else:
                working_gb = id_gb
                logging.info(
                    f"[FINAL_NO_GENE] {sample_name}: '{which_gene_first}' not found "
                    f"in chosen_gb; orientation unchanged.")
        except Exception as exc:
            logging.warning(f"[FINAL_REORIENT_FAIL] {sample_name}: {exc} — using ID-rewritten GB")
            working_gb = id_gb

        # ── Step 3: Optional metadata overlay (header section, text-level) ──
        metadata_mode = "rewritten_no_metadata"

        if preserve_metadata and metadata_src_gb and metadata_src_gb.exists():
            if gb_file_has_meaningful_metadata(metadata_src_gb):
                meta_out = td / "with_metadata.gbf"
                ok, n_replaced, n_unmatched, unmatched, reason, compat = transfer_gb_metadata_file(
                    metadata_src_gb=metadata_src_gb,
                    new_body_gb=working_gb,
                    output_gb=meta_out,
                    smart=smart_metadata,
                    auto_accept=auto_accept,
                )
                if ok and n_replaced > 0:
                    _safe_copy(meta_out, draft_gb)
                    metadata_mode = "metadata_transferred"
                    if compat:
                        metadata_mode += "+compat"
                    if unmatched:
                        logging.warning(
                            f"[META_UNMATCHED] {sample_name}: {len(unmatched)} record(s) "
                            f"unmatched: {unmatched}")
                elif ok and n_replaced == 0:
                    _safe_copy(working_gb, draft_gb)
                    metadata_mode = "metadata_no_match"
                    logging.info(
                        f"[META_NO_MATCH] {sample_name}: 0 records matched in "
                        f"{metadata_src_gb.name}; using ID-rewritten GB.")
                else:
                    _safe_copy(working_gb, draft_gb)
                    metadata_mode = "metadata_failed"
                    logging.warning(f"[META_FAIL] {sample_name}: {reason}")
            else:
                _safe_copy(working_gb, draft_gb)
                metadata_mode = "metadata_src_not_usable"
        else:
            _safe_copy(working_gb, draft_gb)

        # ── Step 3.5: Sync /organism= in FEATURES source ──
        try:
            org_fixed = td / "org_fixed.gbf"
            applied = False

            if metadata_src_gb and metadata_src_gb.exists():
                org_index = build_organism_index([metadata_src_gb])
                if org_index:
                    n_rep, n_nf, n_tot = replace_organism_in_gb_file(draft_gb, org_fixed, org_index)
                    if n_rep > 0:
                        _safe_copy(org_fixed, draft_gb)
                        applied = True
                        logging.info(
                            f"[ORG_SYNC] {sample_name}: /organism= restored in "
                            f"{n_rep}/{n_tot} record(s) from original GB"
                            + (f" ({n_nf} not matched)" if n_nf else ""))

            if not applied:
                n_synced = _sync_organism_qualifier(draft_gb, org_fixed)
                if n_synced > 0:
                    _safe_copy(org_fixed, draft_gb)
                    logging.info(
                        f"[ORG_SYNC] {sample_name}: /organism= synced from header "
                        f"in {n_synced} record(s)")

        except Exception as exc:
            logging.warning(
                f"[ORG_SYNC_FAIL] {sample_name}: organism qualifier sync failed: {exc}")

        # ── Step 4: recount from definitive file ──────────────────────
        try:
            pcg, trna, rrna = counts_from_genbank(draft_gb)
            # cox1_first from the file itself: a linear record is never rotated (rotate_record), whatever
            # reorient_one_gb counted
            orientation_ok = starts_with_gene(draft_gb, which_gene_first, force_forward)
        except Exception:
            pass
        os.replace(str(draft_gb), str(dest_gb))
        committed = True

    finally:
        if not committed:
            try:
                draft_gb.unlink()
            except Exception:
                pass
        shutil.rmtree(td, ignore_errors=True)

    return pcg, trna, rrna, metadata_mode, orientation_ok


# ─────────────────────────────────────────────────────────────────────────────
# Workflows
# ─────────────────────────────────────────────────────────────────────────────
KNOWN_LIMITATION = ("[KNOWN_LIMITATION] the annotation core does not annotate introns (group I/II); they are added "
                    "afterwards by mito_intron.py (MFannot method) and mito_intron_ref.py (intron sites of "
                    "annotated relatives; --intron_fix) when exonerate is available - check intron_report.tsv "
                    "and intron_ref_report.tsv.")


def workflow_annotate(args: argparse.Namespace) -> int:
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "pipeline.log")
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))
    preflight(args)

    in_path     = Path(args.input).resolve()
    split_dir   = out_root / "00.internal_fasta"
    gb_conv_dir = out_root / "00.gb_converted_fasta"

    in_files, input_type, acc_to_orig_gb = iter_inputs(in_path, gb_to_fasta_dir=gb_conv_dir)
    logging.info(f"[INPUT_TYPE] {input_type}  ({len(in_files)} file(s))")
    samples, id_map = split_or_sanitize_fastas(
        in_files, split_dir, args.if_split_fasta == "yes", int(args.id_maxlen),
        acc_to_orig_gb=acc_to_orig_gb)
    logging.info(f"[ID_MAP] {id_map}  (total: {len(samples)})")

    results = run_annotation_stage(
        samples=samples, phase_dirname="01.anno_round1",
        out_root=out_root, args=args,
        fasta_selector=lambda si: si.internal_fasta)

    summary    = RunSummary()
    final_dir  = out_root / "04.final_gb"
    final_dir.mkdir(parents=True, exist_ok=True)
    summary_rows: List[str] = []

    for si in samples:
        r       = results.get(si.internal_id, {})
        best_gb = Path(r["best_gb"]) if r.get("best_gb") else None
        if not r.get("ok") or not best_gb or not best_gb.exists():
            summary.annotation_failed.append(si.sample_name)
            continue
        dest = final_dir / f"{si.sample_name}.gbf"
        try:
            pcg, trna, rrna, mode, orient_ok = finalize_genbank_output(
                sample_name=si.sample_name, body_src_gb=best_gb,
                metadata_src_gb=si.orig_gb_path, dest_gb=dest,
                which_gene_first=orientation_genes(args),
                force_forward=args.force_forward,
                preserve_metadata=(args.preserve_phase1_metadata == "yes")
                    if hasattr(args, "preserve_phase1_metadata") else False,
            )
            summary_rows.append(f"{si.sample_name}\t{pcg}\t{trna}\t{rrna}\t{orientation_text(orient_ok, si)}\t{mode}")
        except Exception as exc:
            logging.error(f"[FINAL_WRITE_FAIL] {si.sample_name}: {exc}")
            summary.final_write_failed.append(si.sample_name)

    write_nonempty_list(out_root / "annotation_failed.txt",
                        summary.annotation_failed, "# annotation core non-zero exit")
    write_nonempty_list(out_root / "final_write_failed.txt",
                        summary.final_write_failed, "# Failed during final write")
    write_nonempty_table(
        out_root / "annotation-summary.tsv",
        "sample_name\tPCG\ttRNA\trRNA\tcox1_orientation_ok\tmetadata_mode",
        summary_rows)
    run_rescue(args, final_dir, out_root)
    run_rrna_hmm(args, final_dir, out_root)
    run_intron_fix(args, final_dir, out_root)
    reorient_after_introns(args, final_dir)
    merge_dir = run_curate_and_submit(args, final_dir, out_root)
    refresh_summary(out_root, merge_dir)
    all_final_gbs = sorted(                      # rglob: the no_<gene>_gene/ records are merged too
        p for p in merge_dir.rglob("*.gbf")
        if p.is_file() and p.stat().st_size > 0)
    if all_final_gbs:
        merge_gb_files(all_final_gbs, out_root / "all.final.gbf")
    logging.info(f"[PLOT] {run_plots(args, merge_dir, out_root)}")
    (out_root / "run_params.json").write_text(
        json.dumps(vars(args), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    summary.log()
    logging.warning(KNOWN_LIMITATION)
    write_orientation_fallback(out_root, samples)
    if args.keep_intermediate == "no":
        cleanup_intermediates(out_root, samples, summary.all_problem_samples())
    logging.info("[DONE]")
    return 0


def workflow_reorient(args: argparse.Namespace) -> int:
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "pipeline.log")
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))

    in_path  = Path(args.input).resolve()
    gb_files = collect_genbank_files(in_path)
    if not gb_files:
        logging.error("[ERROR] No GB files found.")
        return 2

    out_dir     = out_root / "reoriented"
    no_gene_dir = out_root / f"no_{args.which_gene_first}_gene"
    out_dir.mkdir(parents=True, exist_ok=True)
    no_gene_dir.mkdir(parents=True, exist_ok=True)

    summary      = RunSummary()
    summary_rows: List[str] = []

    for gb_in in gb_files:
        sample     = safe_name(gb_in.stem)
        gb_out     = out_dir     / f"{sample}.gbf"
        gb_no_gene = no_gene_dir / f"{sample}.gbf"
        try:
            n_ok, n_no, no_gene_ids = reorient_one_gb(
                gb_in, gb_out, gb_no_gene, orientation_genes(args), args.force_forward)
            if n_no > 0:
                logging.warning(f"[NO_TARGET_GENE] '{sample}': "
                                 f"{n_no} record(s) missing '{args.which_gene_first}'")
                summary.no_target_gene.append(sample)
            if n_ok > 0 and gb_out.exists():
                pcg, trna, rrna = counts_from_genbank(gb_out)
                summary_rows.append(f"{sample}\t{pcg}\t{trna}\t{rrna}\tyes")
            elif n_ok == 0 and n_no > 0:
                _safe_copy(gb_in, gb_out)
                pcg, trna, rrna = counts_from_genbank(gb_out)
                summary_rows.append(f"{sample}\t{pcg}\t{trna}\t{rrna}\tno")
                logging.warning(f"[FALLBACK_COPY] '{sample}': original copied unchanged")
        except Exception as exc:
            logging.error(f"[FAIL] {sample}: {exc}")
            summary.reorient_failed.append(sample)

    write_nonempty_list(out_root / "no_target_gene_samples.txt",
                        summary.no_target_gene, f"# samples without {args.which_gene_first}")
    write_nonempty_list(out_root / "failed_samples.txt",
                        summary.all_problem_samples(), "# all problematic samples")
    write_nonempty_table(out_root / "annotation-summary.tsv",
                         "sample_name\tPCG\ttRNA\trRNA\tcox1_orientation_ok",
                         summary_rows)
    summary.log()
    logging.info("[DONE]")
    return 0


def orientation_text(orient_ok: bool, si: "SampleInfo") -> str:
    """cox1_orientation_ok column: yes / no / no (linear, not rotated)."""
    return "yes" if orient_ok else "no (linear, not rotated)" if si.topology == "linear" else "no"


def write_orientation_fallback(out_root: Path, samples: List["SampleInfo"]) -> None:
    """oriented_by_fallback.txt: samples whose first gene was absent and that start at a fallback gene."""
    names = {si.internal_id: si.sample_name for si in samples}
    rows = sorted({f"{names.get(k.split('.')[0], k)}\t{g}" for k, g in ORIENTATION_FALLBACK.items()})
    write_nonempty_table(out_root / "oriented_by_fallback.txt", "sample_name\tstart_gene", rows)


def announce(msg: str) -> None:
    """A decision made before the run's log file exists (resolve_tools): log it and show it on the console."""
    logging.info(msg)
    print(msg, file=sys.stderr)


def resolve_tools(args: argparse.Namespace) -> str:
    """Fill unset tool locations from mito_env (config file, PATH, conda environments) and check the clade.
    Returns an error message, '' when fine."""
    if str(Path(__file__).resolve().parent) not in sys.path:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
    import mito_env
    if args.clade == "auto":
        clade, code, phylum, err = auto_clade(Path(args.input), getattr(args, "genetic_code", None))
        if err and not err.startswith("WARNING"):
            return f"ERROR: {err}"
        if err:
            print(f"\n{'=' * 68}\n  [CLADE_DEFAULT] {err}\n{'=' * 68}\n", file=sys.stderr)
        args.clade = clade
        if getattr(args, "genetic_code", None) is None:
            args.genetic_code = code
        if hasattr(args, "rrna_hmm") and not getattr(args, "phylum", None):
            args.phylum = phylum
        announce(f"[CLADE] auto: {phylum or 'no taxonomy'} -> --clade {args.clade} --genetic_code {args.genetic_code}")
    if args.clade not in ALL_CLADES:
        return (f"ERROR: --clade '{args.clade}' is not a supported clade. Choose one of:\n  " + ", ".join(ALL_CLADES)
                + "\n  (Animal = all-animal protein database, for any phylum not listed.)")
    if hasattr(args, "rrna_hmm") and not getattr(args, "phylum", None) and args.clade != "Animal":
        args.phylum = CLADE_PHYLUM.get(args.clade)       # explicit --clade: its phylum selects the rRNA HMM library
    if getattr(args, "genetic_code", 0) is None:
        if args.clade in EXTRA_CLADES:
            return (f"ERROR: --clade {args.clade} needs --genetic_code (e.g. 4 Cnidaria, 9 Platyhelminthes, "
                    "5 most invertebrates; NCBI table of your taxon)")
        args.genetic_code = CLADE_CODES.get(args.clade, 5)
    if not args.intron_bin:
        args.intron_bin = mito_env.bin_with("exonerate", "cmsearch", "blastp", "makeblastdb", "tblastn", "blastn")
    args.mitos_bin = mito_env.mitos_bin(args.mitos_bin)
    args.mitos_refdir = mito_env.mitos_refdir(args.mitos_refdir)
    engine = getattr(args, "engine", "native")
    if engine == "mitoz":
        import mito_engines
        if args.clade not in BASE_CLADES:
            return (f"ERROR: --engine mitoz: MitoZ has no protein database for --clade {args.clade} (MitoZ clades: "
                    f"{', '.join(BASE_CLADES)}); use --engine native or mitos2")
        if not mito_env.find("mitoz", getattr(args, "mitoz_bin", None)):
            return f"ERROR: --engine mitoz: MitoZ not found (install: {mito_engines.MITOZ_INSTALL}; or give --mitoz_bin)"
    if engine == "mitos2" and not (args.mitos_bin and args.mitos_refdir):
        import mito_engines
        return ("ERROR: --engine mitos2: MITOS2 (runmitos) or its reference data refseq89m not found (install: "
                f"{mito_engines.MITOS2_INSTALL}; reference data: mitom setup; or give --mitos_bin / --mitos_refdir)")
    if engine == "input":                                  # mitom import: nothing can annotate again
        args.if_reannotation_2 = "no"
    return ""


def prepare_import(args: argparse.Namespace) -> str:
    """mitom import: existing results (mito_engines.import_results) -> <out>/00.imported/imported.gb, which the run
    workflow then reads with engine 'input' (no annotation program runs, round 2 off). Returns an error, '' when fine."""
    import mito_engines
    try:
        recs, notes = mito_engines.import_results(Path(args.input), [Path(x) for x in args.fasta], args.genetic_code,
                                                  args.topology)
    except (ValueError, OSError) as e:
        return f"ERROR: {e}"
    if not recs:
        return f"ERROR: no annotated GenBank record and no MITOS2 result (result.bed) found in {args.input}"
    out = Path(args.out_root).resolve() / "00.imported"
    out.mkdir(parents=True, exist_ok=True)
    SeqIO.write(recs, str(out / "imported.gb"), "genbank")
    for note in notes:
        announce(f"[IMPORT] {note}")
    announce(f"[IMPORT] {len(recs)} annotated record(s) from {args.input}")
    args.input, args.engine = str(out / "imported.gb"), "input"
    return ""


def run_curate_and_submit(args: argparse.Namespace, final_dir: Path, out_root: Path) -> Path:
    """MitoM final steps: mito_curate.py (boundary rules, NCBI gene names) -> 05.curated/, then the NCBI
    submission inputs per record in 06.ncbi_submission/: <id>.fsa + <id>.tbl (BankIt), <id>.sqn + <id>.val
    (table2asn; --sbt adds the submitter template) and submission_summary.tsv. Returns the folder of the
    final GenBank files (05.curated, or final_dir with --curate no)."""
    cur = run_curate(args, final_dir, out_root)
    if cur != final_dir:
        run_submission(args, cur, out_root)
    return cur


def aggressive_start_genes(clade: Optional[str]) -> str:
    """mito_curate.py --aggressive_start_genes for a clade: nad5 / nad4L start at the most upstream in-frame
    start codon in NCBI annelid records (7 / 7, annelid benchmark 2026-10) but not in other phyla (173 metazoa
    records: 17 right, 52 wrong)."""
    return "ND5,ND4L" if (clade or "").lower().startswith("annelida") else ""


def abbreviated_overlap(clade: Optional[str]) -> List[str]:
    """mito_curate.py --abbreviated_overlap for the Chordata clade: vertebrate RefSeq records end a gene whose complete
    stop would run 1-2 nt into the next tRNA / gene with T / TA (leave-out benchmark: +14 exact ends in Chordata;
    other phyla keep the complete stop)."""
    return ["--abbreviated_overlap"] if (clade or "").lower().startswith("chordata") else []


def no_upstream_atg(clade: Optional[str]) -> List[str]:
    """mito_curate.py --no_upstream_atg for nematodes and arthropods, whose genes often start with TTG / ATT / ATA
    (moving such starts to an upstream ATG: 6 right, 14 wrong in the leave-out benchmark; other phyla +37)."""
    return ["--no_upstream_atg"] if (clade or "").lower().startswith(("nematoda", "arthropoda")) else []


def run_curate(args: argparse.Namespace, final_dir: Path, out_root: Path) -> Path:
    """mito_curate.py on final_dir -> 05.curated/ (returned); final_dir when --curate no or curation failed."""
    if getattr(args, "curate", "yes") == "no":
        return final_dir
    here = Path(__file__).resolve().parent
    cur = out_root / "05.curated"
    p = subprocess.run([sys.executable, "-W", "ignore", str(here / "mito_curate.py"), "-i", str(final_dir), "-o", str(cur),
                        "--genetic_code", str(int(args.genetic_code)), "--no_validate",
                        "--aggressive_start_genes", aggressive_start_genes(getattr(args, "clade", None)),
                        *abbreviated_overlap(getattr(args, "clade", None)),
                        *no_upstream_atg(getattr(args, "clade", None))],
                       capture_output=True, text=True)
    if p.returncode:
        logging.error(f"[CURATE] failed (exit {p.returncode}): {p.stderr[-300:]} - final files left uncurated")
        return final_dir
    logging.info(f"[CURATE] {(p.stdout.strip().splitlines() or [''])[-1]}  -> {cur}")
    return cur


def run_submission(args: argparse.Namespace, cur: Path, out_root: Path) -> None:
    """NCBI table2asn validation of every final record and the submission files -> 06.ncbi_submission/."""
    here = Path(__file__).resolve().parent
    if str(here) not in sys.path:
        sys.path.insert(0, str(here))
    import mito_env
    import ncbi_validate as nv
    t2a = mito_env.find("table2asn")
    if not t2a:
        logging.info("[SUBMISSION] table2asn not found: no .sqn / validation (install.sh installs it)")
        return
    sub = out_root / "06.ncbi_submission"
    sub.mkdir(parents=True, exist_ok=True)
    summary, msgs = [], []
    for gbf in sorted(cur.rglob("*.gbf")):
        for rec in SeqIO.parse(str(gbf), "genbank"):
            rows = nv.validate_record(rec, sub / "_work", t2a, getattr(args, "sbt", None))
            d = sub / "_work" / re.sub(r"[^\w.]+", "_", rec.id)
            for ext in (".fsa", ".tbl", ".sqn", ".val"):
                if (d / f"rec{ext}").exists():
                    shutil.copy2(d / f"rec{ext}", sub / f"{d.name}{ext}")
            n = {sev: sum(r["severity"] == sev for r in rows) for sev in ("ERROR", "WARNING")}
            summary.append(f"{rec.id}\t{n['ERROR']}\t{n['WARNING']}\t"
                           + "; ".join(sorted({r["code"] for r in rows if r["severity"] == "ERROR"})))
            msgs += [f"{r['record']}\t{r['severity']}\t{r['code']}\t{r['message']}" for r in rows]
    shutil.rmtree(sub / "_work", ignore_errors=True)
    write_nonempty_table(sub / "submission_summary.tsv", "record\tERROR\tWARNING\terror_codes", summary)
    write_nonempty_table(sub / "validation_messages.tsv", "record\tseverity\tcode\tmessage", msgs)
    bad = [x.split("\t")[0] for x in summary if x.split("\t")[1] != "0"]
    logging.info(f"[SUBMISSION] {len(summary)} record(s) -> {sub}; NCBI ERROR in {len(bad)}: {bad[:10]}"
                 + (" - see validation_messages.tsv" if bad else ""))


F = TypeVar("F", bound=Callable[..., Any])


def optional_step(tag: str) -> Callable[[F], F]:
    """An optional post-annotation step that must never abort the run: any exception is logged as
    '[TAG] failed: ...' and the workflow goes on (intron fix, curation, submission files, merge)."""
    def deco(fn: F) -> F:
        @functools.wraps(fn)
        def wrapper(*a: Any, **k: Any) -> Any:
            try:
                return fn(*a, **k)
            except Exception as exc:
                logging.error(f"[{tag}] failed: {type(exc).__name__}: {exc} - step skipped, the run continues")
                return None
        return wrapper  # type: ignore[return-value]
    return deco


@optional_step("RESCUE")
def run_rescue(args: argparse.Namespace, final_dir: Path, out_root: Path) -> None:
    """Add genes the core missed from a MITOS2 annotation (mito_rescue.py, --rescue yes; MITOS2 is optional), in
    place; report rescue_report.tsv."""
    if getattr(args, "rescue", "no") == "no":
        return
    if args.rescue == "auto":
        # the native core retries other loci / blastx placement itself; MITOS2 took 5-7 min and ~1 GB per record
        # lacking a gene (mostly atp8 of nematodes / flatworms, which have none) for 2 genes in 24 records
        logging.info("[RESCUE] skipped: native core (--rescue yes runs MITOS2 for records missing a gene)")
        return
    if not (args.mitos_bin and args.mitos_refdir):
        msg = "[RESCUE] skipped: MITOS2 not found (--mitos_bin / --mitos_refdir)"
        logging.error(msg) if args.rescue == "yes" else logging.info(msg)
        return
    import mito_rescue
    rows = mito_rescue.rescue_dir(final_dir, int(args.genetic_code), jobs=int(getattr(args, "max_tasks", 1) or 1),
                                  mitos_bin=args.mitos_bin, refdir=args.mitos_refdir)
    write_nonempty_table(out_root / "rescue_report.tsv", "file\trecord\tgene\tlocation\tstatus\thint",
                         ["\t".join([r[k] for k in ("file", "record", "gene", "location", "status")] + [r.get("hint", "")])
                          for r in rows])
    added = [f"{r['record']}:{r['gene']}" for r in rows if r["status"] == "added"]
    failed = [r["record"] for r in rows if r["status"] != "added"]
    logging.info(f"[RESCUE] {len(added)} gene(s) added from MITOS2 {added}" + (f"; MITOS2 failed: {failed}" if failed else ""))


@optional_step("RRNA_HMM")
def run_rrna_hmm(args: argparse.Namespace, final_dir: Path, out_root: Path) -> None:
    """Last-resort 12S/16S search with per-phylum profile HMMs for rRNAs still missing after the cores and
    the MITOS2 rescue (mito_rrna_hmm.py), in place; report rrna_hmm_report.tsv."""
    mode = getattr(args, "rrna_hmm", "no")
    if mode == "no":
        return
    import mito_env
    import mito_rrna_hmm
    lib = Path(args.rrna_hmm_lib) if getattr(args, "rrna_hmm_lib", None) else mito_rrna_hmm.LIB
    if not (mito_rrna_hmm.read_manifest(lib) and mito_env.find("nhmmer")):
        msg = f"[RRNA_HMM] skipped: no HMM library in {lib} or nhmmer not found"
        logging.error(msg) if mode == "yes" else logging.info(msg)
        return
    rows = mito_rrna_hmm.search_dir(final_dir, lib, getattr(args, "phylum", None))
    write_nonempty_table(out_root / "rrna_hmm_report.tsv",
                         "file\trecord\tphylum\tgene\tlocation\tstatus\tcoverage\tnote",
                         ["\t".join([r[k] for k in ("file", "record", "phylum", "gene", "location", "status")]
                                    + [str(r.get(k, "")) for k in ("coverage", "note")]) for r in rows])
    added = [f"{r['record']}:{r['gene']}" for r in rows if str(r["status"]).startswith("added")]
    logging.info(f"[RRNA_HMM] {len(added)} rRNA(s) added from per-phylum HMMs {added}")
    n_skip = sum(str(r["status"]).startswith(("skipped: no library", "skipped: no model")) for r in rows)
    if n_skip:
        msg = (f"[RRNA_HMM] WARNING: {n_skip} missing rRNA(s) were not searched because no HMM model was chosen for "
               f"their phylum (phylum unknown, or no model for it); give the phylum with --phylum, e.g. --phylum Mollusca "
               f"(see rrna_hmm_report.tsv)")
        logging.warning(msg)
        print(msg, file=sys.stderr)



def run_intron_fix(args: argparse.Namespace, final_dir: Path, out_root: Path) -> str:
    """Add introns to final GenBank files in place with mito_intron.py (de novo; protein library MFannot's, or
    --intron_lib), then mito_intron_ref.py (intron sites of annotated relatives, run_intron_ref); returns the line for
    the step summary."""
    if args.intron_fix == "no":
        return "skipped (--intron_fix no)"
    script = Path(__file__).resolve().with_name("mito_intron.py")
    exe = (str(Path(args.intron_bin) / "exonerate") if args.intron_bin else shutil.which("exonerate")) or ""
    if args.intron_bin and not (Path(args.intron_bin) / "blastp").exists():
        exe = ""
    if not script.exists() or not Path(exe).exists():          # on the console: introns would be missing silently
        msg = "[INTRON_FIX] skipped: exonerate / blastp not found (install.sh installs them; or --intron_bin)"
        logging.error(msg) if args.intron_fix == "yes" else logging.warning(msg)
        return "skipped: exonerate / blastp not found (--intron_bin)"
    cmd = [sys.executable, str(script), "-i", str(final_dir), "-o", str(final_dir),
           "--genetic_code", str(int(args.genetic_code)), "--threads", str(int(args.threads)),
           "--jobs", str(int(getattr(args, "max_tasks", 1) or 1))]
    for lib in args.intron_lib:
        cmd += ["--lib", lib]
    for ref in args.intron_ref_gb:
        cmd += ["--ref_gb", ref]
    if args.intron_bin:
        cmd += ["--bin", args.intron_bin]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    report = final_dir / "intron_report.tsv"
    if report.exists():
        shutil.move(str(report), str(out_root / "intron_report.tsv"))
    last = (proc.stdout.strip().splitlines() or [f"failed (exit {proc.returncode}), see intron_fix.log"])[-1]
    logging.info(f"[INTRON_FIX] {last}")
    ref, ref_log = run_intron_ref(args, final_dir, out_root)
    (out_root / "intron_fix.log").write_text(proc.stdout + proc.stderr + ref_log, encoding="utf-8")
    return (last.replace("[mito_intron] ", "").split(" -> ")[0] + f"; {ref}; per intron: intron_report.tsv, "
            "intron_ref_report.tsv")


def run_intron_ref(args: argparse.Namespace, final_dir: Path, out_root: Path) -> Tuple[str, str]:
    """mito_intron_ref.py on the de novo output, in place: the intron sites of the library models/intron_sites.*
    (193 NCBI records, verified gold standard 2026-10-29) re-place de novo boundaries and add introns de novo misses
    (mini-exons, trans-spliced pieces, rRNA introns). --phylum gives the lineage of records without a taxonomy.
    Returns (summary, log text)."""
    script = Path(__file__).resolve().with_name("mito_intron_ref.py")
    cmd = [sys.executable, str(script), "annotate", "-i", str(final_dir), "-o", str(final_dir),
           "--genetic_code", str(int(args.genetic_code)), "--jobs", str(int(getattr(args, "max_tasks", 1) or 1))]
    if args.intron_bin:
        cmd += ["--bin", args.intron_bin]
    if getattr(args, "phylum", None):
        cmd += ["--lineage", f"Eukaryota;Metazoa;{args.phylum}"]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    report = final_dir / "intron_ref_report.tsv"
    if report.exists():
        shutil.move(str(report), str(out_root / "intron_ref_report.tsv"))
    last = (proc.stderr.strip().splitlines() or [f"failed (exit {proc.returncode}), see intron_fix.log"])[-1]
    if proc.returncode:
        logging.warning(f"[INTRON_REF] {last}")
        return f"site transfer failed (exit {proc.returncode}), see intron_fix.log", proc.stdout + proc.stderr
    logging.info(f"[INTRON_REF] {last}")
    return "site transfer: " + last.replace("[mito_intron_ref] ", "").split(" -> ")[0], proc.stdout + proc.stderr


def reorient_after_introns(args: argparse.Namespace, final_dir: Path) -> int:
    """The intron step can turn the orientation gene into exons across the origin: MitoM rotated the record to the
    core's piece of cox1, and when that piece is exon 2, exon 1 and the intron lie at the record's end (annelid cox1
    group II introns). Rotate such circular records again so the gene's 5'-most exon starts the sequence. Returns the
    number of records rotated."""
    genes, n = orientation_genes(args), 0
    for gbf in sorted(final_dir.rglob("*.gbf")):
        recs, changed = list(SeqIO.parse(str(gbf), "genbank")), False
        for i, rec in enumerate(recs):
            feat, _ = find_orientation_feature(rec, genes)
            if rec.annotations.get("topology") != "circular" or feat is None or feat.location.strand == -1:
                continue
            pos = feature_start_pos(feat)
            if pos:
                recs[i], changed = rotate_record(rec, pos), True
                n += 1
                logging.info(f"[REORIENT_INTRON] {rec.id}: {feat.qualifiers.get('gene', ['?'])[0]} starts at {pos + 1} "
                             "after the intron step: rotated to its first exon")
        if changed:
            SeqIO.write(recs, str(gbf), "genbank")
    return n


INTERMEDIATE_DIRS = (
    "00.internal_fasta", "00.gb_converted_fasta", "00.reoriented_fasta_internal",
    "01.anno_round1", "03.anno_round2_reoriented", "04.final_body_merged", "00.imported",
)


def cleanup_intermediates(out_root: Path, samples: List[SampleInfo],
                          problem_samples: List[str]) -> None:
    """Delete internal-ID intermediates once final output is written.

    Kept: 04.final_gb/, all.final.gbf, *.tsv / *.txt reports, pipeline.log,
    run_params.json.  dedup_report.tsv is moved to out_root; core logs of
    problem samples are copied to failed_logs/ before deletion. `run` then
    arranges what is kept into its final layout (organize_outputs).
    """
    dedup = out_root / "00.gb_converted_fasta" / "dedup_report.tsv"
    if dedup.is_file():
        shutil.move(str(dedup), str(out_root / "dedup_report.tsv"))
    problems = set(problem_samples)
    for si in samples:
        if si.sample_name not in problems:
            continue
        for phase in ("01.anno_round1", "03.anno_round2_reoriented"):
            wd = out_root / phase / si.internal_id
            for log in [wd / "sample.log", *sorted((wd / "logs").glob("*.log"))]:
                if log.is_file():
                    dst = out_root / "failed_logs" / f"{si.sample_name}.{phase}.{log.name}"
                    dst.parent.mkdir(exist_ok=True)
                    shutil.copy2(log, dst)
    freed = 0
    for name in INTERMEDIATE_DIRS:
        d = out_root / name
        if d.is_dir():
            freed += sum(p.stat().st_size for p in d.rglob("*") if p.is_file())
            shutil.rmtree(d, ignore_errors=True)
    logging.info(f"[CLEANUP] intermediates removed, {freed / 1e6:.1f} MB freed "
                 f"(use --keep_intermediate yes to keep them)")


def run_phase2(phase2_samples: List["SampleInfo"], args: argparse.Namespace, out_root: Path,
               reor_fasta_dir: Path) -> Dict[str, Dict]:
    """Round 2: annotate the reoriented sequences of `phase2_samples`. Nothing to do (every sample passed the
    --if_reannotation_2 auto check) -> no annotation stage at all: Snakemake would otherwise get an empty rule."""
    logging.info(f"[PHASE2] re-annotation: {len(phase2_samples)} samples")
    if not phase2_samples:
        return {}
    return run_annotation_stage(
        samples=phase2_samples, phase_dirname="03.anno_round2_reoriented",
        out_root=out_root, args=args,
        fasta_selector=lambda si: reor_fasta_dir / f"{si.internal_id}.fasta")


def fallback_src_tag(mode: str, iid: str, reannotated_ids: set) -> str:
    """annotation-summary body_source of a sample whose final body is the shifted round-1 file: 'shifted_round1'
    when --if_reannotation_2 auto decided that round 2 is not needed, 'phase1b_reoriented' when round 2 was
    meant to run but failed (or is off)."""
    return "shifted_round1" if mode == "auto" and iid not in reannotated_ids else "phase1b_reoriented"



def count_rows(table: Path, column: Optional[str] = None, value: Any = None) -> int:
    """Data rows of a TSV report (0 when missing); with column, rows where it equals value / startswith value /
    value(cell) is true."""
    if not table.is_file():
        return 0
    lines = [x for x in table.read_text(encoding="utf-8", errors="ignore").splitlines() if x and not x.startswith("#")]
    if len(lines) < 2:
        return 0
    if column is None:
        return len(lines) - 1
    head = lines[0].split("\t")
    if column not in head:
        return 0
    i = head.index(column)
    cells = [r.split("\t")[i] if len(r.split("\t")) > i else "" for r in lines[1:]]
    if callable(value):
        return sum(bool(value(c)) for c in cells)
    return sum(c == value or str(c).startswith(str(value)) for c in cells)


# the output folder after `run` (organize_outputs); README.txt in the folder repeats this
OUTPUT_LAYOUT = [
    ("final_annotation/", "MitoM annotation, one GenBank file per record (<accession>.gb) - the result"),
    ("final_annotation_all.gb", "the same records in one file"),
    ("original_annotation/", "the input records unchanged (GenBank input only), same file names"),
    ("original_vs_mitom.tsv", "original vs MitoM annotation, one row per gene"),
    ("summary.tsv", "one row per record: gene counts, orientation, round, NCBI errors"),
    ("ncbi_submission/", "NCBI submission files per record (.fsa .tbl .sqn) and table2asn results"),
    ("plots/", "circular genome map per record (.png, .pdf): genes, GC content, GC skew, read depth with --infq"),
    ("reports/", "one report per step (deduplication, round-2 decisions, curation, ...)"),
    ("run.log", "everything the run did, with all commands"),
    ("work/", "intermediate files (only with --keep_intermediate yes)"),
]
REPORT_NOTES = {
    "dedup_report.tsv": "deduplication: one row per input record with a sequence, kept / dropped and why",
    "reannotation_decision.tsv": "after re-orientation: round 2 or not, and why",
    "curation_report.tsv": "every boundary / name change: rule, old and new location, reason",
    "rescue_report.tsv": "genes added from MITOS2",
    "rrna_hmm_report.tsv": "rRNAs searched / added by the per-phylum rRNA HMMs",
    "intron_report.tsv": "introns found",
    "intron_ref_report.tsv": "introns placed from the intron sites of annotated relatives (models/intron_sites)",
    "metadata_transfer_log.tsv": "where each record's GenBank header came from",
    "oriented_by_fallback.txt": "records oriented by another gene because cox1 was not found",
    "all_problem_samples.txt": "records with any problem during the run",
    "run_params.json": "all options of this run",
}


REPORT_FILES = tuple(REPORT_NOTES) + (
    "intron_fix.log", "annotation_failed.txt", "final_write_failed.txt", "reorient_failed.txt",
    "metadata_transfer_failed.txt", "metadata_key_mismatch.txt", "no_target_gene_samples.txt", "failed_samples.txt",
    "annotation-summary.tsv", "pipeline.log",
)                                                   # top-level files `run` writes; organize_outputs moves these to reports/
WORK_DIRS = INTERMEDIATE_DIRS + ("04.final_gb", "05.curated", "_snakemake_01.anno_round1",
                                 "_snakemake_03.anno_round2_reoriented")   # folders of the run: work/ or deleted


def _plot_file(job: Tuple[Path, Path, List[str], int]) -> Tuple[int, int, List[str]]:
    """run_plots worker, one final GenBank file: (maps written, maps with read depth, messages for the log)."""
    import mito_plot
    gbf, out, reads, threads = job
    n = depth_n = 0
    msgs: List[str] = []
    for rec in SeqIO.parse(str(gbf), "genbank"):
        depth = None
        if reads:
            try:
                depth = mito_plot.depth_from_reads(rec, reads, threads=threads)
                depth_n += 1
            except Exception as exc:                       # noqa: BLE001  (a figure must not end the run)
                msgs.append(f"[PLOT] {rec.id}: no depth ring ({exc})")
        try:
            mito_plot.plot_record(rec, out, depth=depth)
            n += 1
        except Exception as exc:                           # noqa: BLE001
            msgs.append(f"[PLOT] {rec.id}: {exc!r}")
    return n, depth_n, msgs


def run_plots(args: argparse.Namespace, gb_dir: Path, out_root: Path) -> str:
    """mito_plot.py: one circular map per final record -> out_root/plots/<record>.png / .pdf (--max_tasks processes,
    ~1.7 s per record); with --infq the reads of the record's sample add a sequencing-depth ring (minimap2 when installed, else bwa mem; --threads
    each). Returns the line for the step summary; never raises."""
    if getattr(args, "plot", "yes") == "no":
        return "skipped (--plot no)"
    try:
        import mito_plot  # noqa: F401  (matplotlib)
    except Exception as exc:                               # noqa: BLE001  (a figure must not end the run)
        logging.error(f"[PLOT] mito_plot.py not usable: {exc!r}")
        return "skipped: mito_plot.py not usable (see run.log)"
    infq = Path(args.infq).resolve() if getattr(args, "infq", None) else None
    s1, s2 = parse_suffix_fq(args.suffix_fq) if infq else ("", "")
    jobs = [(gbf, out_root / "plots",
             [str(x) for x in find_fq_pair_for_sample(gbf.stem, infq, args.fq_position, s1, s2) if x] if infq else [],
             int(args.threads)) for gbf in sorted(gb_dir.rglob("*.gbf"))]
    n = depth_n = 0
    try:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(max(1, min(len(jobs), int(getattr(args, "max_tasks", 1) or 1)))) as pool:
            for k, d, msgs in pool.imap(_plot_file, jobs):          # matplotlib is not thread-safe: processes
                n, depth_n = n + k, depth_n + d
                for m in msgs:
                    logging.warning(m)
    except Exception as exc:                               # noqa: BLE001
        logging.error(f"[PLOT] failed: {exc!r}")
    return f"{n} map(s) in plots/" + (f", {depth_n} with read depth" if depth_n else "")


def write_summary(out_root: Path) -> None:
    """summary.tsv from annotation-summary.tsv (internal names) + NCBI error counts, with readable columns."""
    src = out_root / "annotation-summary.tsv"
    if not src.is_file():
        return
    errs: Dict[str, str] = {}
    sub = out_root / "ncbi_submission" / "submission_summary.tsv"
    if sub.is_file():
        for line in sub.read_text(encoding="utf-8").splitlines()[1:]:
            c = line.split("\t")
            if len(c) > 1:
                errs[c[0]] = c[1]
    rows = ["record\tprotein_genes\ttRNAs\trRNAs\tcox1_first\tannotation_used\toriginal_header_restored\tncbi_errors"]
    for line in src.read_text(encoding="utf-8").splitlines()[1:]:
        c = line.split("\t") + [""] * 7
        name, pcg, trn, rrn, ori, body, meta = c[:7]
        used = ("round 2" if body.startswith("phase2") else
                "round 1, shifted" if body.startswith(("shifted_round1", "phase1b")) else
                "round 1, not re-oriented" if body.startswith("phase1_raw") else body)
        if "partial_no_gene_merged" in body:
            used += " + extra contigs"
        rows.append("\t".join([name, pcg, trn, rrn, ori, used, "yes" if "transferred" in meta else "no",
                               errs.get(name, "")]))
    (out_root / "summary.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    src.unlink()


def organize_outputs(out_root: Path, final_gb_dir: Path, keep: bool = False) -> None:
    """Put the run's output into the layout of OUTPUT_LAYOUT (final_annotation/, reports/ ...) and write
    README.txt. Intermediate folders go to work/ (keep) or are deleted. Only the files and folders MitoM writes
    are touched (the output folder may hold the user's own files); the results of a previous run into the same
    folder are replaced."""
    def remove(p: Path) -> None:
        if p.is_dir() and not p.is_symlink():
            shutil.rmtree(p)
        elif p.exists() or p.is_symlink():
            p.unlink()

    def replace(src: Path, dst: Path) -> None:                # a previous run's copy must not swallow the new one
        remove(dst)
        shutil.move(str(src), str(dst))

    for name in ("final_annotation", "final_annotation_all.gb", "ncbi_submission", "reports", "summary.tsv"):
        remove(out_root / name)                                # results of a previous run into the same folder
    fin = out_root / "final_annotation"
    fin.mkdir()
    for f in sorted(final_gb_dir.rglob("*.gbf")):
        shutil.move(str(f), str(fin / f"{f.stem}.gb"))
    reports = out_root / "reports"
    reports.mkdir()
    kept_dedup = out_root / "00.gb_converted_fasta" / "dedup_report.tsv"     # --keep_intermediate yes
    if kept_dedup.is_file():
        shutil.copy2(kept_dedup, reports / "dedup_report.tsv")
    for f in final_gb_dir.glob("*.tsv"):                      # curation_report.tsv
        shutil.move(str(f), str(reports / f.name))
    if (out_root / "all.final.gbf").exists():
        shutil.move(str(out_root / "all.final.gbf"), str(out_root / "final_annotation_all.gb"))
    if (out_root / "06.ncbi_submission").is_dir():
        shutil.move(str(out_root / "06.ncbi_submission"), str(out_root / "ncbi_submission"))
    write_summary(out_root)
    for name in REPORT_FILES:                                  # MitoM's own files only: the folder may hold the user's
        if (out_root / name).is_file():
            replace(out_root / name, reports / name)
    if (out_root / "failed_logs").is_dir():
        replace(out_root / "failed_logs", reports / "failed_logs")
    for name in WORK_DIRS:
        d = out_root / name
        if not d.is_dir():
            continue
        if keep:
            (out_root / "work").mkdir(exist_ok=True)
            replace(d, out_root / "work" / name)
        else:
            shutil.rmtree(d, ignore_errors=True)
    lines = ["MitoM output", "=" * 60, ""]
    lines += [f"{n:28s} {w}" for n, w in OUTPUT_LAYOUT if (out_root / n.rstrip("/")).exists()]
    lines += ["", "reports/", "-" * 60]
    lines += [f"  {n:28s} {w}" for n, w in REPORT_NOTES.items() if (reports / n).exists()]
    others = sorted(f.name for f in reports.iterdir() if f.name not in REPORT_NOTES)
    if others:
        lines += [f"  {n}" for n in others]
    lines += ["", "summary.tsv columns", "-" * 60,
              "  annotation_used: round 1, shifted = round-1 annotation moved with the sequence (passed the check);",
              "                   round 2 = annotated again after re-orientation; not re-oriented = no orientation gene",
              "  original_header_restored: the GenBank header of the input record (organism, references) was copied back",
              "  ncbi_errors: NCBI table2asn ERROR count (0 = ready for submission)",
              "", "original_vs_mitom.tsv result column", "-" * 60,
              "  same / ends_differ / name_differs / strand_differs / only_original / only_mitom / same_place",
              "  locations are on the original record's sequence; difference_nt = larger shift of the two ends"]
    (out_root / "README.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def check_transl_tables(orig_dir: Path, code: int) -> List[str]:
    """[CODE] WARNING per input GenBank record whose CDS /transl_table differs from the genetic code in use
    (a wrong --genetic_code gives internal stops without any other sign). Returns the records' names."""
    bad: List[str] = []
    if not orig_dir.is_dir():
        return bad
    for gb in sorted(orig_dir.glob("*.gb")):
        try:
            for rec in SeqIO.parse(str(gb), "genbank"):
                tables = sorted({str(f.qualifiers["transl_table"][0]) for f in rec.features
                                 if f.type == "CDS" and "transl_table" in f.qualifiers})
                if tables and tables != [str(code)]:
                    bad.append(rec.id)
                    logging.warning(f"[CODE] {rec.id}: record transl_table={'/'.join(tables)} differs from "
                                    f"--genetic_code {code} — check")
        except Exception as exc:                   # noqa: BLE001  (a check only)
            logging.info(f"[CODE] {gb.name}: not checked: {exc!r}")
    return bad


def workflow_run(args: argparse.Namespace) -> int:
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "run.log", console_level=logging.WARNING, mode="w")   # one run, one log
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))
    preflight(args)
    steps = Steps(11)
    say("")
    say("MitoM - mitochondrial genome annotation")
    say(f"    input        : {Path(args.input).resolve()}")
    say(f"    output       : {out_root}")
    say(f"    clade        : {args.clade}    genetic code: {args.genetic_code}")
    say(f"    parallel     : {args.max_tasks} records at a time x {args.threads} threads")
    say(f"    annotation   : {ENGINE_LABEL[getattr(args, 'engine', 'native')]}")
    say(f"    round 2      : {args.if_reannotation_2}  (auto = only records whose annotation is broken by the cox1 shift)")
    say(f"    full log     : {out_root / 'run.log'}")

    in_path     = Path(args.input).resolve()
    split_dir   = out_root / "00.internal_fasta"
    gb_conv_dir = out_root / "00.gb_converted_fasta"

    steps.start("Read input and remove duplicate records")
    in_files, input_type, acc_to_orig_gb = iter_inputs(in_path, gb_to_fasta_dir=gb_conv_dir,
                                                       orig_dir=out_root / "original_annotation")
    logging.info(f"[INPUT_TYPE] {input_type}  ({len(in_files)} file(s))")
    samples, id_map = split_or_sanitize_fastas(
        in_files, split_dir, args.if_split_fasta == "yes", int(args.id_maxlen),
        acc_to_orig_gb=acc_to_orig_gb)
    total = len(samples)
    logging.info(f"[ID_MAP] {id_map}  (total: {total})")
    if input_type == "converted_from_gb":
        st = LAST_DEDUP_STATS
        steps.info(f"GenBank records read: {st.get('total_in', '?')}")
        dropped = [(st.get(k, 0), why) for k, why in (
            ("dup_refseq", "GenBank record that a RefSeq NC_ record was made from (NC_ kept)"),
            ("dup_version", "older version of the same accession"), ("dup_exact", "same accession twice"),
            ("dup_seq", "identical sequence"), ("no_sequence", "no sequence"))]
        for n, why in dropped:
            if n:
                steps.info(f"removed {n}: {why}")
        steps.info("the original annotation of every kept record is saved for comparison (original_annotation/)")
        steps.done(f"{total} records to annotate; reasons per record: reports/dedup_report.tsv")
    else:
        steps.done(f"{total} FASTA sequences to annotate")

    summary = RunSummary()
    check_transl_tables(out_root / "original_annotation", int(args.genetic_code))

    steps.start(f"Round-1 annotation ({ENGINE_LABEL[getattr(args, 'engine', 'native')]}: protein genes, rRNAs, tRNAs)")
    results1 = run_annotation_stage(
        samples=samples, phase_dirname="01.anno_round1",
        out_root=out_root, args=args,
        fasta_selector=lambda si: si.internal_fasta)
    for si in samples:
        if not results1.get(si.internal_id, {}).get("ok"):
            summary.annotation_failed.append(si.sample_name)
    steps.done(f"{total - len(summary.annotation_failed)} annotated, {len(summary.annotation_failed)} failed"
               + (f" ({', '.join(summary.annotation_failed[:5])}{' ...' if len(summary.annotation_failed) > 5 else ''})"
                  if summary.annotation_failed else ""))

    steps.start(f"Re-orient: plus strand, sequence starts with {args.which_gene_first} "
                f"(if missing: {args.fallback_genes})")

    logging.info(
        f"[PHASE1b] reorient gene='{args.which_gene_first}' force_forward={args.force_forward}")
    reor_fasta_dir  = out_root / "00.reoriented_fasta_internal"
    reor_fasta_dir.mkdir(parents=True, exist_ok=True)
    reor_ok_ids:    List[str]        = []
    no_gene_iids:   set              = set()
    partial_no_gene_gb: Dict[str, Path] = {}     # contigs without the gene: merged back into the final record

    for si in samples:
        r1      = results1.get(si.internal_id, {})
        best_gb = Path(r1["best_gb"]) if r1.get("best_gb") else None
        if not best_gb or not best_gb.exists():
            continue
        reor_dir   = out_root / "01.anno_round1" / si.internal_id / "reoriented"
        reor_dir.mkdir(parents=True, exist_ok=True)
        gb_reor    = reor_dir / f"{si.internal_id}.reoriented.gbf"
        gb_no_gene = reor_dir / f"{si.internal_id}.no_{args.which_gene_first}.gbf"
        try:
            n_ok, n_no, no_gene_ids = reorient_one_gb(
                best_gb, gb_reor, gb_no_gene, orientation_genes(args), args.force_forward)
            if n_ok == 0 and n_no > 0:
                logging.warning(
                    f"[NO_TARGET_GENE] '{si.sample_name}' ({si.internal_id}): "
                    f"ALL {n_no} record(s) missing '{args.which_gene_first}' — skipping Phase-2")
                summary.no_target_gene.append(si.sample_name)
                no_gene_iids.add(si.internal_id)
            elif n_no > 0:
                logging.info(
                    f"[PARTIAL_NO_GENE] '{si.sample_name}' ({si.internal_id}): "
                    f"{n_no} minor contig(s) lack '{args.which_gene_first}'; "
                    f"primary reoriented OK — proceeding to Phase-2.")
                if gb_no_gene.exists():
                    partial_no_gene_gb[si.internal_id] = gb_no_gene
            if n_ok > 0 and gb_reor.exists():
                reor_ok_ids.append(si.internal_id)
                gb_to_internal_fasta(gb_reor, reor_fasta_dir / f"{si.internal_id}.fasta",
                                     si.internal_id)
        except Exception as exc:
            logging.error(f"[REORIENT_FAIL] {si.sample_name} ({si.internal_id}): {exc}")
            summary.reorient_failed.append(si.sample_name)      # an error, not a record without the gene
            no_gene_iids.add(si.internal_id)

    fallback_used = sum(1 for si in samples if si.internal_id in reor_ok_ids
                        and any(k.startswith(si.internal_id) for k in ORIENTATION_FALLBACK))
    steps.done(f"{len(reor_ok_ids)} re-oriented"
               + (f" ({fallback_used} by a fallback gene)" if fallback_used else "")
               + (f", {len(no_gene_iids)} without any orientation gene (kept as annotated)" if no_gene_iids else ""))

    steps.start("Check the annotation after re-orientation; annotate again where it is offset (round 2)")
    do_phase2 = args.if_reannotation_2 in ("yes", "auto")
    results2:  Dict[str, Dict] = {}
    reannotated_ids: set = set()
    if do_phase2:
        phase2_samples = [si for si in samples if si.internal_id in reor_ok_ids]
        n_checked = len(phase2_samples)
        if args.if_reannotation_2 == "auto":
            phase2_samples = select_reannotation(phase2_samples, results1, out_root,
                                                 int(getattr(args, "genetic_code", None) or 5))
            steps.info(f"checked {n_checked}: {len(phase2_samples)} need round 2 (gene cut at the old ends, gene "
                       f"changed by the shift, gene across the new start, standard gene missing); "
                       f"{n_checked - len(phase2_samples)} keep the shifted round-1 annotation")
            steps.info("decision and reason per record: reports/reannotation_decision.tsv")
        else:
            steps.info(f"--if_reannotation_2 yes: all {n_checked} re-oriented records are annotated again")
        results2 = run_phase2(phase2_samples, args, out_root, reor_fasta_dir)
        reannotated_ids = {si.internal_id for si in phase2_samples}
        n_fail2 = 0
        for si in phase2_samples:
            if not results2.get(si.internal_id, {}).get("ok"):
                n_fail2 += 1
                logging.warning(
                    f"[PHASE2_FAIL_FALLBACK] '{si.sample_name}': "
                    f"phase2 failed, fallback to phase1b reoriented GB")
        steps.done(f"round 2: {len(phase2_samples) - n_fail2} annotated"
                   + (f", {n_fail2} failed (they keep the shifted round-1 annotation)" if n_fail2 else ""))
    else:
        steps.done("--if_reannotation_2 no: skipped")

    steps.start("Write the annotation; copy the original GenBank header (organism, references ...) back")

    final_dir         = out_root / "04.final_gb"
    no_gene_final_dir = final_dir / f"no_{args.which_gene_first}_gene"
    final_dir.mkdir(parents=True, exist_ok=True)
    no_gene_final_dir.mkdir(parents=True, exist_ok=True)

    summary_rows:    List[str] = []
    meta_log_rows:   List[str] = []

    for si in samples:
        sample      = si.sample_name
        iid         = si.internal_id
        is_no_gene  = iid in no_gene_iids
        r1          = results1.get(iid, {})
        metadata_src = si.orig_gb_path

        chosen_gb: Optional[Path] = None
        src_tag   = "phase1_raw"

        if do_phase2 and not is_no_gene:
            r2 = results2.get(iid, {})
            if r2.get("best_gb") and Path(r2["best_gb"]).exists():
                chosen_gb = Path(r2["best_gb"])
                src_tag   = "phase2"

        if chosen_gb is None:
            gb_reor = (out_root / "01.anno_round1" / iid / "reoriented"
                       / f"{iid}.reoriented.gbf")
            if gb_reor.exists():
                chosen_gb = gb_reor
                src_tag   = fallback_src_tag(args.if_reannotation_2 if do_phase2 else "no", iid, reannotated_ids)

        if chosen_gb is None and r1.get("best_gb") and Path(r1["best_gb"]).exists():
            chosen_gb = Path(r1["best_gb"])

        if chosen_gb is None or not chosen_gb.exists():
            if sample not in summary.annotation_failed:
                summary.final_write_failed.append(sample)
            continue

        body_for_final = chosen_gb
        extra_no_gene_gb = partial_no_gene_gb.get(iid)
        if extra_no_gene_gb and extra_no_gene_gb.exists():
            chosen_is_phase1_raw = False
            if r1.get("best_gb"):
                try:
                    chosen_is_phase1_raw = chosen_gb.resolve() == Path(r1["best_gb"]).resolve()
                except Exception:
                    chosen_is_phase1_raw = False
            if not chosen_is_phase1_raw:
                merged_body_dir = out_root / "04.final_body_merged"
                merged_body_dir.mkdir(parents=True, exist_ok=True)
                merged_body = merged_body_dir / f"{iid}.with_no_gene.gbf"
                merge_gb_files([chosen_gb, extra_no_gene_gb], merged_body)
                body_for_final = merged_body
                src_tag += "+partial_no_gene_merged"
                logging.info(
                    f"[FINAL_PARTIAL_NO_GENE_MERGE] '{sample}': merged "
                    f"{extra_no_gene_gb.name} back into final body.")

        dest = (no_gene_final_dir if is_no_gene else final_dir) / f"{sample}.gbf"

        try:
            pcg, trna, rrna, meta_mode, orient_ok = finalize_genbank_output(
                sample_name=sample,
                body_src_gb=body_for_final,                          # with the contigs without the gene
                metadata_src_gb=metadata_src,
                dest_gb=dest,
                which_gene_first=orientation_genes(args),
                force_forward=args.force_forward,
                preserve_metadata=(args.preserve_phase1_metadata == "yes"),
                smart_metadata=(args.smart_metadata == "yes"),
                auto_accept=(args.auto_accept_metadata == "yes"),
            )
            summary_rows.append(
                f"{sample}\t{pcg}\t{trna}\t{rrna}\t{orientation_text(orient_ok, si)}\t{src_tag}\t{meta_mode}")
            meta_log_rows.append(
                f"{sample}\t{iid}\t{metadata_src or 'NA'}\t{body_for_final}\t{dest}\t{meta_mode}")

            if "transferred" in meta_mode:
                summary.metadata_transferred.append(sample)
                if "compat" in meta_mode:
                    summary.metadata_compat_match.append(sample)
                    summary.metadata_key_mismatch.append(sample)
            elif meta_mode in ("metadata_no_match", "metadata_src_not_usable"):
                summary.metadata_skipped.append(sample)
            elif meta_mode == "metadata_failed":
                summary.metadata_transfer_failed.append(sample)

        except Exception as exc:
            logging.error(f"[FINAL_WRITE_FAIL] {sample}: {exc}")
            summary.final_write_failed.append(sample)

    write_nonempty_table(
        out_root / "annotation-summary.tsv",
        "sample_name\tPCG\ttRNA\trRNA\tcox1_orientation_ok\tbody_source\tmetadata_mode",
        summary_rows)
    steps.done(f"{len(summary_rows)} records written"
               + (f", original header restored in {len(summary.metadata_transferred)}" if summary.metadata_transferred else "")
               + (f", {len(summary.final_write_failed)} failed" if summary.final_write_failed else ""))
    write_nonempty_table(
        out_root / "metadata_transfer_log.tsv",
        "sample_name\tinternal_id\tmetadata_src\tbody_src\tdest\tmode",
        meta_log_rows)
    write_nonempty_list(out_root / "annotation_failed.txt",
                        summary.annotation_failed, "# annotation core returned non-zero exit")
    write_nonempty_list(out_root / f"no_{args.which_gene_first}_samples.txt",
                        summary.no_target_gene,
                        f"# samples where '{args.which_gene_first}' was not found")
    write_nonempty_list(out_root / "final_write_failed.txt",
                        summary.final_write_failed, "# failed during final GB write")
    write_nonempty_list(out_root / "reorient_failed.txt",
                        summary.reorient_failed, "# failed during reorientation")
    write_nonempty_list(out_root / "metadata_transfer_failed.txt",
                        summary.metadata_transfer_failed, "# metadata transfer failed")
    write_nonempty_list(out_root / "metadata_key_mismatch.txt",
                        summary.metadata_key_mismatch,
                        "# accession key mismatch during metadata transfer")
    write_nonempty_list(out_root / "all_problem_samples.txt",
                        summary.all_problem_samples(), "# all problematic samples")
    (out_root / "run_params.json").write_text(
        json.dumps(vars(args), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    steps.start("Add what the core missed: MITOS2 rescue (optional), 12S/16S by rRNA HMM, introns")
    run_rescue(args, final_dir, out_root)
    steps.info(f"genes added by MITOS2 rescue : {count_rows(out_root / 'rescue_report.tsv', 'status', 'added')}")
    run_rrna_hmm(args, final_dir, out_root)
    steps.info(f"rRNAs added by rRNA HMM      : {count_rows(out_root / 'rrna_hmm_report.tsv', 'status', 'added')}")
    steps.info(f"introns                      : {run_intron_fix(args, final_dir, out_root)}")
    n_rot = reorient_after_introns(args, final_dir)
    if n_rot:
        steps.info(f"re-oriented to the first exon : {n_rot} record(s) whose {args.which_gene_first} gained an exon upstream")
    steps.done()

    steps.start("Boundary rules for protein genes and rRNAs, NCBI gene names")
    merge_dir = run_curate(args, final_dir, out_root)
    curate_skipped = ("skipped (--curate no)" if args.curate == "no" else
                      "curation failed (see run.log): the annotation before curation is used")
    steps.done(f"{count_rows(merge_dir / 'curation_report.tsv')} changes, each with rule and reason in "
               f"reports/curation_report.tsv" if merge_dir != final_dir else curate_skipped)

    steps.start("NCBI check (table2asn) and submission files")
    if merge_dir != final_dir:
        run_submission(args, merge_dir, out_root)
        sub_sum = out_root / "06.ncbi_submission" / "submission_summary.tsv"
        n_rec = count_rows(sub_sum)
        n_err = count_rows(sub_sum, "ERROR", lambda v: v != "0")
        steps.done(f"{n_rec} records checked: {n_rec - n_err} without NCBI ERROR, {n_err} with ERROR "
                   f"(ncbi_submission/validation_messages.tsv)" if n_rec else
                   ("validation produced no records (see run.log)" if sub_sum.parent.is_dir() else
                    "table2asn not found: skipped"))
    else:
        steps.done(curate_skipped)
    refresh_summary(out_root, merge_dir)
    all_final_gbs = sorted(                      # rglob: the no_<gene>_gene/ records are merged too
        p for p in merge_dir.rglob("*.gbf")
        if p.is_file() and p.stat().st_size > 0
    )
    if all_final_gbs:
        n_merged = merge_gb_files(all_final_gbs, out_root / "all.final.gbf")
        logging.info(f"[MERGED_ALL] {out_root / 'all.final.gbf'}  ({n_merged} records)")
    else:
        logging.warning("[MERGED_ALL] No final GB files to merge.")

    n_ok = len(summary_rows)
    logging.info(f"[FINAL_DIR] {final_dir}")
    logging.info(f"[OUTPUT_COUNT] {n_ok}/{total}")
    logging.info(f"[METADATA_MODE] preserve={args.preserve_phase1_metadata}")
    summary.log()
    logging.info(KNOWN_LIMITATION)
    write_orientation_fallback(out_root, samples)
    if args.keep_intermediate == "no":
        cleanup_intermediates(out_root, samples, summary.all_problem_samples())

    steps.start("Compare with the original annotation")
    orig_dir = out_root / "original_annotation"
    if input_type == "converted_from_gb" and orig_dir.is_dir():
        try:
            import mito_compare
            n = mito_compare.compare_dirs(orig_dir, merge_dir, out_root / "original_vs_mitom.tsv",
                                          int(args.genetic_code))
            steps.done(mito_compare.summary_text(n) + "; one row per gene: original_vs_mitom.tsv")
        except Exception as exc:                   # noqa: BLE001  (a report: must not end the run after annotation)
            logging.error(f"[COMPARE] failed: {exc!r}")
            steps.done("comparison failed (see run.log); the annotation itself is complete")
    else:
        steps.done("FASTA input has no original annotation: skipped")

    steps.start("Genome maps (circular plots)")
    steps.done(run_plots(args, merge_dir, out_root))

    steps.start("Organise the output")
    organize_outputs(out_root, merge_dir, keep=args.keep_intermediate == "yes")
    for name, what in OUTPUT_LAYOUT:
        if (out_root / name).exists():
            steps.info(f"{name:28s} {what}")
    steps.done(f"explained in {out_root / 'README.txt'}")
    say("")
    say(f"MitoM finished in {duration(time.time() - steps.run_t0)}: {out_root}")
    logging.info("[DONE]")
    if total and len(summary.annotation_failed) == total:
        logging.warning(f"[FAIL] every record failed in the annotation core ({total}): see reports/failed_logs/")
        return 1
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# Standalone transfer_metadata workflow
# ─────────────────────────────────────────────────────────────────────────────
def workflow_transfer_metadata(args: argparse.Namespace) -> int:
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "transfer_metadata.log")
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))

    metadata_files = collect_genbank_files(Path(args.metadata).resolve())
    body_files     = collect_genbank_files(Path(args.new_body).resolve())
    if not metadata_files:
        logging.error(f"[ERROR] No metadata GB files: {args.metadata}")
        return 2
    if not body_files:
        logging.error(f"[ERROR] No body GB files: {args.new_body}")
        return 2

    out_dir = out_root / "transferred_gb"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_rows: List[str] = []

    def _process(body_file: Path, meta_file: Optional[Path], out_file: Path) -> str:
        if args.if_replace_metadata == "no":
            if args.overwrite_new_body == "no":
                _safe_copy(body_file, out_file)
            log_rows.append(f"{body_file}\t{meta_file or 'NA'}\t{out_file}\tskipped_by_user\t0\t0")
            return "skip"
        if meta_file is None:
            if args.overwrite_new_body == "no":
                _safe_copy(body_file, out_file)
            log_rows.append(f"{body_file}\tNA\t{out_file}\tno_match\t0\t0")
            return "no_match"
        ok, n_rep, n_un, unmatched, reason, compat = transfer_gb_metadata_file(
            meta_file, body_file, out_file,
            smart=(args.smart_metadata == "yes"), auto_accept=(args.auto_accept_metadata == "yes"))
        if not ok and args.overwrite_new_body == "no":
            _safe_copy(body_file, out_file)
        log_rows.append(
            f"{body_file}\t{meta_file}\t{out_file}\t"
            f"{'ok' if ok else 'no_usable_metadata'}\t{n_rep}\t{n_un}")
        return "ok" if ok else "skip"

    is_single = (len(metadata_files) == 1 and len(body_files) == 1
                 and Path(args.metadata).is_file() and Path(args.new_body).is_file())
    if is_single:
        out_name = getattr(args, "output_name", None) or body_files[0].name
        out_file = body_files[0] if args.overwrite_new_body == "yes" else (out_dir / out_name)
        _process(body_files[0], metadata_files[0], out_file)
        write_nonempty_table(out_root / "transfer_metadata.tsv",
                             "body_file\tmetadata_file\toutput_file\tstatus\treplaced\tunmatched",
                             log_rows)
        return 0

    meta_by_stem: Dict[str, Path] = {}           # the same stem in two folders is an error, never a silent choice
    dup_stems: Dict[str, List[Path]] = {}
    for mf in metadata_files:
        if mf.stem in meta_by_stem:
            dup_stems.setdefault(mf.stem, [meta_by_stem[mf.stem]]).append(mf)
        else:
            meta_by_stem[mf.stem] = mf
    if dup_stems:
        for stem, paths in sorted(dup_stems.items()):
            logging.error(
                f"[DUP_METADATA_STEM] '{stem}' appears in multiple metadata files: "
                + ", ".join(str(p) for p in paths))
        return 2
    n_ok = n_skip = n_fail = 0
    for bf in body_files:
        matched  = meta_by_stem.get(bf.stem)
        out_file = bf if args.overwrite_new_body == "yes" else (out_dir / bf.name)
        try:
            status = _process(bf, matched, out_file)
            if status == "ok":
                n_ok += 1
            else:
                n_skip += 1
        except Exception:
            if args.overwrite_new_body == "no":
                _safe_copy(bf, out_file)
            n_fail += 1
            log_rows.append(f"{bf}\t{matched or 'NA'}\t{out_file}\tfail\t0\t0")
    write_nonempty_table(out_root / "transfer_metadata.tsv",
                         "body_file\tmetadata_file\toutput_file\tstatus\treplaced\tunmatched",
                         log_rows)
    logging.info(f"[SUMMARY] transferred={n_ok}  skipped={n_skip}  failed={n_fail}")
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# Standalone replace_organism workflow
# ─────────────────────────────────────────────────────────────────────────────
def workflow_replace_organism(args: argparse.Namespace) -> int:
    """Replace /organism= qualifiers in FEATURES source blocks (`run` does this in step 3.5 of
    finalize_genbank_output)."""
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "replace_organism.log")
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))

    meta_files = collect_genbank_files(Path(args.metadata).resolve())
    body_files = collect_genbank_files(Path(args.new_body).resolve())
    if not meta_files:
        logging.error(f"[ERROR] No metadata GB files: {args.metadata}")
        return 2
    if not body_files:
        logging.error(f"[ERROR] No body GB files: {args.new_body}")
        return 2

    logging.info(f"[ORG_INDEX] Building organism+source index from {len(meta_files)} file(s)…")
    org_index = build_organism_index(meta_files)
    logging.info(f"[ORG_INDEX] {len(org_index)} accession variant(s) indexed.")
    if not org_index:
        logging.error("[ERROR] Organism index is empty — check metadata GB files.")
        return 2

    out_dir = out_root / "replaced_organism"
    out_dir.mkdir(parents=True, exist_ok=True)

    log_rows: List[str] = []
    n_ok = n_miss = n_fail = 0

    for bf in body_files:
        out_file = bf if args.overwrite_new_body == "yes" else (out_dir / bf.name)
        try:
            n_rep, n_nf, n_tot = replace_organism_in_gb_file(bf, out_file, org_index)
            log_rows.append(f"{bf.name}\t{out_file}\t{n_rep}\t{n_nf}\t{n_tot}\t")         # note column empty
            if n_rep > 0:
                n_ok += 1
                logging.info(
                    f"[DONE] {bf.name}: {n_rep}/{n_tot} record(s) updated, "
                    f"{n_nf} unchanged.")
            else:
                n_miss += 1
                logging.warning(
                    f"[NO_MATCH] {bf.name}: 0 records matched — organism/source not replaced.")
        except Exception as exc:
            n_fail += 1
            logging.error(f"[FAIL] {bf.name}: {exc}")
            log_rows.append(f"{bf.name}\t{out_file}\t0\t0\t0\tFAIL: {exc}")
            _safe_copy(bf, out_file)

    write_nonempty_table(
        out_root / "replace_organism.tsv",
        "source_file\toutput_file\treplaced\tnot_found\ttotal\tnote",
        log_rows)
    logging.info(
        f"[SUMMARY] files_updated={n_ok}  no_match={n_miss}  failed={n_fail}")
    logging.info("[DONE]")
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# Argument parser
# ─────────────────────────────────────────────────────────────────────────────
_CLADE_DEFAULT = "auto"


def default_threads() -> int:
    """Threads per annotation job when --threads is not given: min(8, CPUs of this machine)."""
    return min(8, os.cpu_count() or 8)


def default_max_tasks() -> int:
    """Concurrent jobs when --max_tasks is not given: CPUs // default threads (at least 1)."""
    return max(1, (os.cpu_count() or 8) // default_threads())


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="mitom" if os.environ.get("MITOM_LAUNCHER") else "mitom_pipeline.py",   # sub-commands show `mitom run`
        description=(
            "MitoM - annotation workflow for animal mitochondrial genomes (MitoM's own annotation core, mito_core.py)\n\n"
            "  run               - full workflow: annotate -> re-orient -> round 2 where needed -> final GenBank ->\n"
            "                      rRNA HMM fallback, introns, curation, NCBI submission files, comparison, genome maps\n"
            "  import            - the same steps on existing MitoZ / MITOS2 results or any annotated GenBank\n"
            "  annotate          - one round of annotation and the steps after it (work layout: 04.final_gb/ ...)\n"
            "  reorient          - re-orient GenBank files to start at a target gene\n"
            "  transfer_metadata - copy GenBank headers from the original records\n"
            "  replace_organism  - restore /organism= of the FEATURES source from the original records\n"
            "  dedup             - deduplicate GenBank records (RefSeq NC_ kept)\n\n"
            "The annotation runs in threads or with Snakemake (--scheduler_backend; cluster: --smk_extra '--profile\n"
            "slurm'); everything after it runs in this process. `mitom <command> -h` lists the options.\n\n"
            "Citation:\n"
            "  Snakemake:  Mölder et al. (2021) F1000Res doi:10.12688/f1000research.29032.2\n"
            "  Biopython:  Cock et al. (2009) Bioinformatics doi:10.1093/bioinformatics/btp163\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("-i", "--input",  required=True,
                        help="Input FASTA file/dir or GenBank file/dir.")
        sp.add_argument("-o", "--out_root", default="mitom_out",
                        help="Output root directory (default: mitom_out).")
        sp.add_argument("--if_split_fasta", choices=["yes", "no"], default="yes")
        sp.add_argument("--id_maxlen", type=id_maxlen_arg, default=15,
                        help="Max length of internal sample IDs (default: 15).")
        sp.add_argument("--which_gene_first", default="cox1",
                        help="Gene used as sequence start (default: cox1).")
        sp.add_argument("--fallback_genes", default="cox2,cox3,cob,nad1",
                        help="Genes tried in order when --which_gene_first is absent "
                             "(default: cox2,cox3,cob,nad1; 'none' = leave such records unrotated). "
                             "Records oriented this way are listed in oriented_by_fallback.txt.")
        sp.add_argument("--force_forward", action="store_true", default=True,
                        help="RC the sequence when target gene is on minus strand (default: on).")
        sp.add_argument("--no_force_forward", action="store_false", dest="force_forward")

    def add_annotation(sp: argparse.ArgumentParser) -> None:
        g = sp.add_argument_group("Annotation")
        g.add_argument("--core", choices=["native"], default="native", help=argparse.SUPPRESS)   # old scripts
        g.add_argument("--mitos_bin", default=None,
                       help="MITOS2 environment bin (runmitos). Default: found automatically (mito_env.py).")
        g.add_argument("--mitos_refdir", default=None,
                       help="Directory containing MITOS2 refseq89m. Default: found automatically.")
        g.add_argument("--plot", choices=["yes", "no"], default="yes",
                       help="Circular genome map per record (mito_plot.py, matplotlib) -> plots/<record>.png and .pdf: "
                            "genes by strand and type, GC content, GC skew; with --infq also the read depth (minimap2 when installed, else bwa mem). "
                            "Default yes (about 2 s per record, --max_tasks in parallel).")
        g.add_argument("--curate", choices=["yes", "no"], default="yes",
                       help="yes (default): mito_curate.py boundary rules and NCBI gene names -> 05.curated/, and "
                            "NCBI submission files (.fsa/.tbl/.sqn/.val) -> 06.ncbi_submission/")
        g.add_argument("--sbt", default=None,
                       help="NCBI submitter template (template.sbt from submit.ncbi.nlm.nih.gov/genbank/template) "
                            "for the .sqn files")
        g.add_argument("--rescue", choices=["auto", "yes", "no"], default="auto",
                       help="Add protein genes / rRNAs the core missed from a MITOS2 annotation (mito_rescue.py; "
                            "nothing annotated is changed; needs the optional MITOS2 installation). auto (default) / "
                            "no = not run (the core retries other loci itself); yes = run for records missing a gene.")
        g.add_argument("--rrna_hmm", choices=["auto", "yes", "no"], default="auto",
                       help="Last-resort 12S/16S search for rRNAs still missing after the cores and the rescue, "
                            "with per-phylum profile HMMs built from NCBI rRNAs (mito_rrna_hmm.py; only missing "
                            "rRNAs are searched). auto = run when the library and nhmmer are found (default).")
        g.add_argument("--rrna_hmm_lib", default=None,
                       help="Per-phylum rRNA HMM library (mito_rrna_hmm.py build). Default: models/rrna_hmm.")
        g.add_argument("--rrna_method", choices=["hmm", "cm"], default="hmm",
                       help="12S/16S search of the core: hmm (default, faster) = per-phylum profile HMMs, "
                            "MITOS2's covariance models where no model or hit; cm = MITOS2's covariance models only. "
                            "Same accuracy (2026-10-28, 450 leave-out records: found 98.5 / 98.5 %%, both ends exact "
                            "47.4 / 47.6 %%); hmm takes seconds and ~0.1 GB, cm ~10 s and ~0.5 GB per gene.")
        g.add_argument("--phylum", default=None,
                       help="Phylum of all input records for --rrna_hmm (default: the phylum found by --clade auto, "
                            "or the phylum of an explicit --clade; with --clade Animal and FASTA input give it here, "
                            "otherwise the record's own taxonomy is used when it has one).")
        g.add_argument("--intron_fix", choices=["auto", "yes", "no"], default="auto",
                       help="Add introns with mito_intron.py (MFannot method: Exonerate "
                            "protein2genome + MFannot boundary rules + Rfam intron type), then with "
                            "mito_intron_ref.py (intron sites of annotated relatives, models/intron_sites). "
                            "auto = run when exonerate is found (default).")
        g.add_argument("--intron_bin", default=None,
                       help="Directory containing exonerate, cmsearch, blastp, tblastn and blastn (e.g. a conda env "
                            "bin).")
        g.add_argument("--intron_lib", action="append", default=[],
                       help="Protein library for intron models in MFannot header format "
                            "('>Species mt cox1; ...'), repeatable. Default: MFannot's curated "
                            "collections bundled in models/mfannot_proteins.")
        g.add_argument("--intron_ref_gb", action="append", default=[],
                       help="GenBank file/dir of close relatives used as extra homologs for intron models "
                            "(repeatable).")
        g.add_argument("--keep_intermediate", choices=["yes", "no"], default="no",
                       help="Keep the intermediate folders (00.*, 01.*, 03.*, 04.*, 05.*, _snakemake_*; run: moved to "
                            "work/). Default no: deleted after the final output is written; core logs of problem "
                            "samples are always kept (run: reports/failed_logs/, annotate: failed_logs/).")
        g.add_argument("--threads",      type=positive_int_arg, default=default_threads(),
                       help="Threads per annotation job (default: min(8, number of CPUs) = %(default)s).")
        g.add_argument("--max_tasks",    type=positive_int_arg, default=default_max_tasks(),
                       help="Max concurrent annotation jobs (default: CPUs // threads = %(default)s).")
        g.add_argument("--clade",        default=_CLADE_DEFAULT,
                       help="Taxonomic group (sets the genetic code and group rules). auto (default): from the taxonomy of the GenBank "
                            "input (one phylum per run; class names such as Vertebrata, Gastropoda, Insecta count "
                            "for their phylum; a group below or Animal). A taxonomy that names no "
                            "known phylum, several phyla, or CDS transl_tables that differ between records are "
                            "errors; input with NO taxonomy at all (FASTA, GenBank without a lineage) falls back to "
                            "Annelida-segmented-worms with a WARNING - give the clade explicitly for any other "
                            f"taxon: {', '.join(ALL_CLADES)}.")
        g.add_argument("--genetic_code", type=int, default=None,
                       help="Genetic code (default with --clade auto: the CDS transl_table of the input, else the "
                            "phylum's usual code; otherwise the code of the clade - 2 Chordata, 9 Echinodermata, "
                            "4 Porifera, else 5; required for Cnidaria, Platyhelminthes, Chaetognatha, Animal).")
        g.add_argument("--species_name", default=None, help=argparse.SUPPRESS)   # hidden; passed to MitoZ (--engine mitoz)

        fq = sp.add_argument_group("FASTQ (optional - read depth ring of the genome maps, --plot)")
        fq.add_argument("--infq", default=None,
                        help="Directory with the paired FASTQ files of the samples (matched by sample name).")
        fq.add_argument("--fq_position", choices=["flat", "sample_dir"], default="flat")
        fq.add_argument("--suffix_fq", default="_1.clean.fq.gz,_2.clean.fq.gz")

        sch = sp.add_argument_group("Scheduler")
        sch.add_argument("--scheduler_backend",
                         choices=["auto", "snakemake", "threads"], default="auto",
                         help=(
                             "Dispatch backend (default: auto).\n"
                             "  auto      — use Snakemake if available, else threads\n"
                             "  snakemake — always use Snakemake\n"
                             "  threads   — always use ThreadPoolExecutor\n"
                             "Annotation is parallelised; finalization runs\n"
                             "in the main process regardless of backend."
                         ))
        sch.add_argument("--snakemake_bin",    default=None)
        sch.add_argument("--smk_cores",        type=positive_int_arg, default=64)
        sch.add_argument("--smk_mem_mb",       type=positive_int_arg, default=256000)
        sch.add_argument("--smk_latency_wait", type=positive_int_arg, default=60)
        sch.add_argument("--smk_extra",        action="append", default=[],
                         metavar="FLAG",
                         help="Extra Snakemake flags appended verbatim.\n"
                              "Example: --smk_extra '--profile slurm'")

    def add_engine(sp: argparse.ArgumentParser) -> None:
        g = sp.add_argument_group("Annotation engine")
        g.add_argument("--engine", choices=["native", "mitoz", "mitos2"], default="native",
                       help="Program that finds the genes: native (default) = MitoM's own core (mito_core.py); mitoz = "
                            "MitoZ 3.6, mitos2 = MITOS2 2.1.10 (not installed with MitoM - install them yourself, see "
                            "mito_engines.py; found on PATH or in a conda environment). Every later MitoM step runs on "
                            "their annotation: reorientation, round 2, curation, NCBI files, genome maps.")
        g.add_argument("--mitoz_bin", default=None,
                       help="Folder with the mitoz program (--engine mitoz). Default: found automatically.")

    def add_metadata(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--preserve_phase1_metadata", choices=["yes", "no"], default="yes",
                        help="Restore original GenBank header metadata into final files "
                             "(default: yes; requires GenBank input).")
        sp.add_argument("--smart_metadata",       choices=["yes", "no"], default="no",
                        help="Enable fuzzy/interactive accession matching (default: no).")
        sp.add_argument("--auto_accept_metadata", choices=["yes", "no"], default="no",
                        help="Auto-accept best fuzzy match (requires --smart_metadata yes).")

    # ── annotate ──────────────────────────────────────────────────────────
    a = sub.add_parser("annotate", help="Batch annotation only (one round of the annotation engine).")
    add_common(a)
    add_annotation(a)
    add_engine(a)
    a.add_argument("--preserve_phase1_metadata", choices=["yes", "no"], default="no")

    # ── reorient ──────────────────────────────────────────────────────────
    r = sub.add_parser("reorient", help="Reorient GenBank files to start at a target gene.")
    add_common(r)

    # ── run ───────────────────────────────────────────────────────────────
    run_p = sub.add_parser(
        "run",
        help="Full pipeline: annotate → reorient → re-annotate → guaranteed final reorientation.")
    add_common(run_p)
    add_annotation(run_p)
    add_engine(run_p)
    run_p.add_argument("--if_reannotation_2", choices=["yes", "auto", "no"], default="auto",
                       help="Phase-2 re-annotation of the reoriented sequences. auto (default): round 1 runs in "
                            "full and only samples that fail the check after the cox1 shift (gene at the original "
                            "sequence ends, standard gene missing, gene changed by the shift, gene across the new "
                            "origin) are re-annotated, decisions in reannotation_decision.tsv; yes: every sample; "
                            "no: never.")
    add_metadata(run_p)

    # ── import ────────────────────────────────────────────────────────────
    im = sub.add_parser(
        "import",
        help="Existing MitoZ / MITOS2 results or any annotated GenBank -> the MitoM steps without annotating again: "
             "reorientation, curation and NCBI names, NCBI submission files, genome maps (-i: a result folder, "
             "a .gbf / .gb file, or a folder of them).")
    add_common(im)
    add_annotation(im)
    add_metadata(im)
    im.add_argument("--fasta", action="append", default=[],
                    help="Genome FASTA of MITOS2 results (sequence named as in result.bed; repeatable). Default: "
                         "FASTA files in the MITOS2 result folder or the folder above it.")
    im.add_argument("--topology", choices=["auto", "circular", "linear"], default="auto",
                    help="auto (default): MITOS2 from result.gff, GenBank from the LOCUS line, MitoZ output circular "
                         "(MitoZ writes linear for every record; give linear for contigs).")

    # ── transfer_metadata ─────────────────────────────────────────────────
    tm = sub.add_parser("transfer_metadata",
                        help="Standalone GenBank header transfer between two files/dirs.")
    tm.add_argument("-m", "--metadata", required=True,
                    help="Original GenBank file/dir (metadata source).")
    tm.add_argument("-n", "--new_body", required=True,
                    help="Re-annotated GenBank file/dir (body source).")
    tm.add_argument("-o", "--out_root", default="transfer_metadata_out")
    tm.add_argument("--if_replace_metadata",  choices=["yes", "no"], default="yes")
    tm.add_argument("--overwrite_new_body",   choices=["yes", "no"], default="no")
    tm.add_argument("--output_name", default=None,
                    help="Output filename (single-file mode only).")
    tm.add_argument("--smart_metadata",       choices=["yes", "no"], default="no")
    tm.add_argument("--auto_accept_metadata", choices=["yes", "no"], default="no")

    # ── replace_organism ──────────────────────────────────────────────────
    ro = sub.add_parser(
        "replace_organism",
        help=(
            "Restore /organism= and SOURCE text in FEATURES source blocks.\n"
            "Fixes placeholder organism names (e.g. 'Test sp.') that survive\n"
            "transfer_metadata because they live in the body, not the header.\n"
            "In 'run' mode this is applied automatically as Step 3.5."
        ),
    )
    ro.add_argument("-m", "--metadata", required=True,
                    help="Original GenBank file/dir (organism/source name source).")
    ro.add_argument("-n", "--new_body", required=True,
                    help="Re-annotated GenBank file/dir whose /organism= needs fixing.")
    ro.add_argument("-o", "--out_root", default="replace_organism_out",
                    help="Output root directory (default: replace_organism_out).")
    ro.add_argument("--overwrite_new_body", choices=["yes", "no"], default="no",
                    help="Overwrite the input files in-place (default: no).")

    # ── dedup ─────────────────────────────────────────────────────────────
    dd = sub.add_parser(
        "dedup",
        help="Deduplicate GenBank records (same accession / older version / RefSeq copy of an "
             "INSDC submission / identical sequence; RefSeq NC_ kept) and write GenBank + FASTA.")
    dd.add_argument("-i", "--input", required=True, help="GenBank file or directory.")
    dd.add_argument("-o", "--out_root", default="dedup_out", help="Output directory (default: dedup_out).")

    return p


def workflow_dedup(args: argparse.Namespace) -> int:
    """Standalone deduplication: dedup.gb, dedup.fasta, dedup_report.tsv."""
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "dedup.log")
    files = collect_genbank_files(Path(args.input).resolve())
    kept, _, stats = dedup_gb_records(files, dedup_report_path=out_root / "dedup_report.tsv")
    recs = [m.record for m in kept]
    SeqIO.write(recs, str(out_root / "dedup.gb"), "genbank")
    SeqIO.write(recs, str(out_root / "dedup.fasta"), "fasta")
    logging.info(f"[DEDUP] in={stats['total_in']} kept={stats['kept']} no_seq={stats['no_sequence']} "
                 f"dup_exact={stats['dup_exact']} dup_version={stats['dup_version']} "
                 f"dup_refseq={stats['dup_refseq']} dup_seq={stats['dup_seq']} -> {out_root}")
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────
def main() -> int:
    if SeqIO is None:
        print("ERROR: Biopython not installed.  pip install biopython", file=sys.stderr)
        return 2

    args = build_parser().parse_args()

    if args.cmd == "import":
        err = prepare_import(args)
        if err:
            print(err, file=sys.stderr)
            return 2
    if args.cmd in ("annotate", "run", "import"):
        err = resolve_tools(args)
        if err:
            print(err, file=sys.stderr)
            return 2

    if getattr(args, "suffix_fq", None):
        parse_suffix_fq(args.suffix_fq)

    if (getattr(args, "auto_accept_metadata", "no") == "yes"
            and getattr(args, "smart_metadata", "no") != "yes"):
        print("ERROR: --auto_accept_metadata yes requires --smart_metadata yes", file=sys.stderr)
        return 2


    dispatch = {
        "annotate":          workflow_annotate,
        "reorient":          workflow_reorient,
        "run":               workflow_run,
        "import":            workflow_run,
        "transfer_metadata": workflow_transfer_metadata,
        "replace_organism":  workflow_replace_organism,
        "dedup":             workflow_dedup,
    }
    fn = dispatch.get(args.cmd)
    if fn is None:
        print(f"Unknown command: {args.cmd}", file=sys.stderr)
        return 2
    return fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
