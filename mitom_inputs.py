"""mitom_inputs.py - what goes into a MitoM run: the input type (FASTA / GenBank, file or folder), --clade auto from
the GenBank taxonomy, GenBank deduplication (dedup_report.tsv; RefSeq NC_ kept) and conversion to FASTA, one sample per
FASTA record with an internal ID (id_map.tsv), the FASTQ pair of a sample (--infq), and the CDS transl_table check."""
from __future__ import annotations

import dataclasses
import gzip
import io
import logging
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mito_common import GB_SUFFIXES, InputError, safe_name, topology_from_header, truncated_genbank
from mito_taxa import METAZOAN_PHYLA, PHYLUM_CLADE, PHYLUM_CODE, PHYLUM_SYNONYMS
from mitom_metadata import _normalize_acc, gb_record_accession


def header_first_token(desc_or_id: str) -> str:
    s = (desc_or_id or "").strip()
    if s.startswith(">"):
        s = s[1:]
    return s.split()[0] if s else "NA"


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
EMBL_EXTS = {".embl", ".ebl", ".ena"}
_GB_EXTS    = GB_SUFFIXES


def _plain_suffix(p: Path) -> str:
    """The suffix of a file name without a trailing .gz (genome.fasta.gz -> .fasta)."""
    return Path(p.name[:-3]).suffix.lower() if p.name.lower().endswith(".gz") else p.suffix.lower()


NORMALIZED_DIR = "00.input"


def normalize_input(in_path: Path, out_root: Path) -> Tuple[Path, List[str]]:
    """The input as the steps read it: gzip-compressed files (.gz) unpacked and text that is not UTF-8 (Latin-1
    headers of Windows programs) converted, into out_root/00.input/ (2026-10-31: both made the run end with 0
    records). Returns (the path to read, notes); the input itself when nothing needs to change."""
    in_path = Path(in_path)
    if in_path.is_file():
        files, base = [in_path], in_path.parent
    elif in_path.is_dir():
        files, base = sorted(p for p in in_path.rglob("*") if p.is_file()
                             and _plain_suffix(p) in (set(_FASTA_EXTS) | set(_GB_EXTS) | EMBL_EXTS)), in_path
    else:
        return in_path, []

    def needs(p: Path) -> str:
        with p.open("rb") as fh:
            head = fh.read(3)
        if head[:2] == b"\x1f\x8b":
            return "gzip"
        if head == b"\xef\xbb\xbf":                 # UTF-8 with a byte-order mark (Windows Notepad): 0 records
            return "bom"
        with p.open("rb") as fh:
            if fh.read(5) == b"ID   ":                   # EMBL / ENA flat file: written as GenBank
                return "embl"
        try:
            p.read_bytes().decode("utf-8")
        except UnicodeDecodeError:
            return "latin1"
        return ""

    todo = {p: needs(p) for p in files}
    if not any(todo.values()):
        return in_path, []
    dest = Path(out_root) / NORMALIZED_DIR
    if dest.is_dir():
        shutil.rmtree(dest)
    notes = []
    for p, what in todo.items():
        rel = p.relative_to(base)
        t = dest / (rel.parent / rel.name[:-3] if what == "gzip" and rel.name.lower().endswith(".gz") else rel)
        if what == "embl":
            t = t.with_suffix(".gb")
        t.parent.mkdir(parents=True, exist_ok=True)
        data = gzip.open(p, "rb").read() if what == "gzip" else p.read_bytes()
        try:
            text, latin = data.decode("utf-8-sig"), False
        except UnicodeDecodeError:
            text, latin = data.decode("latin-1"), True
        if what == "embl" or (what == "gzip" and text.startswith("ID   ")):
            t = t.with_suffix(".gb") if what == "embl" else t
            recs = list(SeqIO.parse(io.StringIO(text), "embl"))
            for r in recs:                               # (GenBank needs molecule_type)
                r.annotations.setdefault("molecule_type", "DNA")
            with t.open("w", encoding="utf-8") as fh:
                SeqIO.write(recs, fh, "genbank")
        else:
            t.write_text(text, encoding="utf-8")
        if what:
            notes.append(f"{p.name}: " + "; ".join(x for x in ("unpacked" if what == "gzip" else "",
                                                               "EMBL written as GenBank" if what == "embl" else "",
                                                               "byte-order mark removed" if what == "bom" else "",
                                                               "not UTF-8 text: read as Latin-1" if latin else "") if x))
    if in_path.is_file():
        name = in_path.name[:-3] if todo[in_path] == "gzip" and in_path.name.lower().endswith(".gz") else in_path.name
        return (dest / Path(name).with_suffix(".gb").name if todo[in_path] == "embl" else dest / name), notes
    return dest, notes


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


# a FASTA record takes part in the accession rules only when its first word is an accession (MK123456.1, NC_012345.1):
# local names such as NODE_1 or sample_3 recur in unrelated assemblies
_FASTA_ACC_RE = re.compile(r"^(?:[A-Z]{2}_\d{6,9}|[A-Z]{1,6}\d{5,9})(?:\.\d+)?$")


def record_format(path: Path) -> str:
    """'fasta' or 'genbank' for one input file (suffix, else its first line)."""
    return "fasta" if _is_fasta_file(path) else "genbank" if _is_gb_file(path) else _sniff_input_type(path)


def collect_sequence_files(src: Path) -> List[Path]:
    """The GenBank and FASTA files of an input (the file itself, or every such file below a folder), sorted."""
    if src.is_file():
        return [src.resolve()]
    if src.is_dir():
        return sorted(p.resolve() for p in src.rglob("*") if p.is_file() and (_is_gb_file(p) or _is_fasta_file(p)))
    return []


def dedup_gb_records(
    gb_files: List[Path],
    dedup_report_path: Optional[Path] = None,
) -> Tuple[List[GbRecordMeta], Dict[str, Path], Dict]:
    """One record per genome from GenBank and / or FASTA files: (1) same accession -> newest version; (2) a RefSeq
    NC_ record replaces the INSDC submission named in its COMMENT (GenBank only); (3) identical sequence -> one copy
    (RefSeq first). FASTA records take part in (1) only when their first header word is an accession."""
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
        fmt = record_format(gb_path)
        try:
            records = list(SeqIO.parse(str(gb_path), fmt))
        except Exception as exc:
            logging.warning(f"[DEDUP] Cannot parse {gb_path}: {exc}")
            stats["parse_errors"] += 1
            continue
        cut = set(truncated_genbank(gb_path, records)) if fmt == "genbank" else set()

        for i, rec in enumerate(records):
            stats["total_in"] += 1
            if i in cut:
                logging.warning(f"[DEDUP] {rec.id} ({gb_path.name}): the record is cut off (sequence shorter than its "
                                f"LOCUS length, or no closing //): not annotated - download it again")
                stats["parse_errors"] += 1
                continue
            if fmt == "fasta":
                acc_ver  = rec.id
                acc_base = (_normalize_acc(acc_ver) or acc_ver) if _FASTA_ACC_RE.match(acc_ver) else f"{gb_path}#{i}"
            else:
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
            topo     = (topology_from_header(rec.description) if fmt == "fasta"
                        else rec.annotations.get("topology", "circular").lower())
            if topo not in ("circular", "linear"):
                topo = "circular"

            all_metas.append(GbRecordMeta(
                acc_versioned=acc_ver, acc_base=acc_base,
                version_int=ver_int, seq_hash=seq_hash,
                seq_len=len(seq_str), topology=topo,
                source_file=gb_path, record=rec,
            ))

    by_acc: Dict[str, GbRecordMeta] = {}
    names_used = {m.acc_versioned for m in all_metas}
    for meta in all_metas:
        key = meta.acc_base
        prev = by_acc.get(key)
        if prev is not None and prev.version_int == meta.version_int and prev.seq_hash != meta.seq_hash:
            # same name and version, another sequence: not a duplicate (records without ACCESSION fall back to the
            # LOCUS name - contig1 in two files, YZZ20231_1 / _2 read as versions; 2026-10-31: the second was dropped)
            if meta.acc_versioned == prev.acc_versioned:
                n = 2
                while f"{meta.acc_versioned}_{n}" in names_used:
                    n += 1
                logging.warning(f"[DEDUP] {meta.acc_versioned} ({meta.source_file.name}): another sequence under the "
                                f"name of {prev.acc_versioned} ({prev.source_file.name}): kept as "
                                f"{meta.acc_versioned}_{n}")
                meta = dataclasses.replace(meta, acc_versioned=f"{meta.acc_versioned}_{n}")
                names_used.add(meta.acc_versioned)
            k = 2
            while f"{key}#{k}" in by_acc:
                k += 1
            by_acc[f"{key}#{k}"] = meta
            continue
        if key not in by_acc:
            by_acc[key] = meta
        else:
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
    dedup_report: Optional[Path] = None,
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
        dedup_report = dedup_report or gb_to_fasta_dir / "dedup_report.tsv"
        fasta_paths, acc_to_orig_gb = convert_gb_to_fasta_for_annotation(
            gb_files, gb_to_fasta_dir, dedup_report_path=dedup_report, orig_dir=orig_dir)
        if not fasta_paths:
            raise InputError(f"no usable GenBank record in {in_path} (empty, cut off, without sequence, or unreadable: see the warnings)")
        logging.info(
            f"[INPUT] Converted {len(gb_files)} GB file(s) -> "
            f"{len(fasta_paths)} unique FASTA record(s)")
        return fasta_paths, "converted_from_gb", acc_to_orig_gb

    raise FileNotFoundError(f"Input not found or unrecognised: {in_path}")


MIN_CIRCULAR_LEN = 10000    # an unlabelled FASTA record shorter than this is a contig: linear
UNLABELLED = "FASTA without topology=: "   # start of the topology_source of an unlabelled record


MIN_OVERLAP = 50   # nt the end of a sequence must repeat its start to count as an assembly's circular overlap
# words of a FASTA header that tell the topology when there is no topology=: GetOrganelle '(circular)', NCBI
# 'complete genome' / 'partial genome' / 'nearly complete'; checked in this order (partial first)
HEADER_HINTS = ((r"\b(partial (genome|sequence)|nearly complete|incomplete|linear)\b", "linear"),
                (r"(?<!non-)(?<!not )\bcircular\b|\bcircularized\b|\bcomplete genome\b", "circular"))


def header_hint(desc: str) -> Tuple[Optional[str], str]:
    """(topology, the words) a FASTA header's text gives (HEADER_HINTS); (None, '') when it gives none."""
    low = (desc or "").lower()
    for pat, topo in HEADER_HINTS:
        m = re.search(pat, low)
        if m:
            return topo, m.group(0)
    return None, ""


def terminal_overlap(seq: str, min_len: int = MIN_OVERLAP) -> int:
    """Length of the longest end of `seq` that repeats its start (>= min_len nt, at most a quarter of the sequence,
    leaving >= MIN_CIRCULAR_LEN nt, not itself a tandem repeat; 0: none) -
    the overlap an assembler leaves on a circular molecule (MITObim circules.py, MitoHiFi, NOVOPlasty trim it).
    Prefix function of start + separator + end, linear in the length."""
    s = seq.upper()
    m = len(s) // 4                                           # an overlap is a small part of the molecule
    if m < min_len:
        return 0
    k = _prefix_function(s[:m] + "\x00" + s[-m:])[-1]
    if k < min_len or len(s) - k < MIN_CIRCULAR_LEN:
        return 0
    if k - _prefix_function(s[:k])[-1] < k // 2:             # the overlap is a tandem repeat itself (period < k/2):
        return 0                                              # both ends in one repeat array, not a closed circle
    return k


def _prefix_function(t: str) -> List[int]:
    pi = [0] * len(t)
    for i in range(1, len(t)):
        k = pi[i - 1]
        while k and t[i] != t[k]:
            k = pi[k - 1]
        if t[i] == t[k]:
            k += 1
        pi[i] = k
    return pi


def infer_topology_from_description(desc: str) -> str:
    m = re.search(r"topology\s*=\s*(circular|linear)", (desc or "").lower())
    return m.group(1) if m else "circular"


def _is_topology_tag(s: str) -> bool:
    return bool(re.match(r"^topology\s*=\s*(circular|linear)$", (s or "").strip(), re.I))


def make_internal_id(idx: int, max_len: int) -> str:
    base = f"S{idx:05d}"
    return base[:max_len] if len(base) > max_len else base


@dataclass(frozen=True)
class SampleInfo:
    sample_name:    str
    orig_token:     str
    internal_id:    str
    topology:       str
    source_fasta:   Path
    internal_fasta: Path
    orig_gb_path:   Optional[Path] = None
    topology_source: str = ""      # where the label came from (reports/topology_decision.tsv)
    label:          str = ""       # the input's topology (topology: what round 1 annotates; --topology circular /
                                   # linear and linear-mtDNA lineages can differ, mitom_topology)


def write_one_record_fasta(record: "SeqRecord", out_fa: Path, seq_id: str, topology: str) -> None:
    rec = record[:]
    rec.id = rec.name = seq_id
    rec.description = f"topology={topology}"
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with out_fa.open("w", encoding="utf-8") as fh:
        SeqIO.write([rec], fh, "fasta")


_IUPAC = set("ACGTRYKMSWBDHVN")
MIN_ACGT = 500        # a sequence with fewer A/C/G/T than this (empty, all N) is not annotated


def clean_sequence(rec: "SeqRecord") -> Tuple["SeqRecord", str]:
    """The sequence as the programs of the core read it: alignment gaps ('-', '.'), stops ('*'), spaces and digits
    removed, U -> T, letters that are not IUPAC nucleotides -> N (2026-10-31: a '-' made exonerate fail on every
    protein gene). Returns (record, warning, '' when nothing changed); a sequence with fewer than MIN_ACGT A/C/G/T
    comes back empty (skipped, with the warning)."""
    try:
        raw = str(rec.seq)
    except Exception:                                    # (Biopython UndefinedSequenceError: CONTIG / WGS record)
        return SeqRecord(Seq(""), id=rec.id, name=rec.name, description=rec.description), \
            "no sequence (a CONTIG / WGS record): not annotated"
    s = re.sub(r"[-.*\s0-9]", "", raw.upper()).replace("U", "T")
    bad = sum(1 for c in s if c not in _IUPAC)
    if bad:
        s = "".join(c if c in _IUPAC else "N" for c in s)
    acgt = sum(s.count(c) for c in "ACGT")
    notes = []
    letters = re.sub(r"[^A-Z]", "", raw.upper())
    nt = sum(letters.count(c) for c in "ACGTUN")
    if letters and nt < 0.8 * len(letters):          # protein letters are IUPAC codes too (A G K M V ...): 2026-10-31
        return rec[:0], (f"only {100 * nt / len(letters):.0f} % A/C/G/T/N: not a nucleotide sequence (protein?), "
                         f"not annotated")
    if len(s) != len(raw):
        notes.append(f"{len(raw) - len(s)} gap / stop / space characters removed")
    if bad:
        notes.append(f"{bad} letters that are not nucleotides read as N")
    if acgt < MIN_ACGT:
        notes.append(f"only {acgt} A/C/G/T nucleotides (empty or all N): not annotated")
        s = ""
    if not notes:                        # (case kept: only gaps, odd letters and empty sequences are changed)
        return rec, ""
    out = rec[:0]
    out.seq = Seq(s)
    out.id, out.name, out.description = rec.id, rec.name, rec.description
    return out, "; ".join(notes)


def split_or_sanitize_fastas(
    in_files: List[Path],
    out_dir: Path,
    if_split: bool,
    id_maxlen: int,
    acc_to_orig_gb: Optional[Dict[str, Path]] = None,
    topology: Optional[str] = None,
) -> Tuple[List[SampleInfo], Path]:
    """One single-sequence FASTA per sample (out_dir/<internal_id>.fasta, header topology=...). The topology: `topology`
    (--topology circular / linear) when given, else the label - the GenBank LOCUS line, a FASTA header's
    topology=circular|linear, else circular for one sequence of >= MIN_CIRCULAR_LEN nt and linear for several
    sequences in a file or a shorter one."""
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
            rec, why = clean_sequence(rec)
            if why:
                logging.warning(f"[SEQ] {rec.id}: {why}")
            if not len(rec.seq):
                continue
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
            if sample_name != base_name:
                logging.warning(f"[SPLIT] {fa.name}: a second sequence named '{base_name}': annotated as '{sample_name}'")
            orig_gb_path = acc_to_orig_gb.get(safe_name(rec.id)) \
                        or acc_to_orig_gb.get(orig_token)
            topo     = infer_topology_from_description(desc)
            tagged   = re.search(r"topology\s*=\s*(circular|linear)", desc.lower())
            source   = "GenBank LOCUS line" if orig_gb_path and tagged else "FASTA header" if tagged else \
                       UNLABELLED + "decided from the sequence"
            hint, words = (None, "") if tagged else header_hint(desc)
            if hint:                                  # assembler / NCBI words: a label the ends can still overturn
                topo, tagged = hint, True
                source = f"FASTA header text '{words}'"
            k = terminal_overlap(str(rec.seq)) if not orig_gb_path and topo == "circular" else 0
            if k:                                     # an assembler's overlap: circular, the repeated end removed
                logging.warning(f"[TOPOLOGY] {orig_token}: its last {k} nt repeat its first {k} nt (assembly overlap "
                                f"of a circular molecule): removed, {len(rec.seq) - k} nt annotated as circular")
                rec = rec[:len(rec.seq) - k]
                topo, tagged = "circular", True
                source = f"terminal overlap of {k} nt (removed)" + ("" if source.startswith(UNLABELLED)
                                                                     else f"; {source}")
            if topo == "circular" and not tagged and len(rec.seq) < MIN_CIRCULAR_LEN:
                # too short for a whole mitogenome: a contig (rotating / joining across a false origin cuts genes)
                topo = "linear"
                source = UNLABELLED + "shorter than 10 kb"
            # (an unlabelled sequence of >= 10 kb - one in the file or several - is annotated circular first and then
            # decided from its ends and gene set, mitom_topology.decide: assemblies often hold several mitogenomes)
            label = topo
            if topology in ("circular", "linear"):
                topo = topology
            internal_id = make_internal_id(idx, id_maxlen)
            out_fa      = out_dir / f"{internal_id}.fasta"
            write_one_record_fasta(rec, out_fa, seq_id=internal_id, topology=topo)
            si = SampleInfo(
                sample_name=sample_name, orig_token=orig_token,
                internal_id=internal_id, topology=topo,
                source_fasta=fa.resolve(), internal_fasta=out_fa.resolve(),
                orig_gb_path=orig_gb_path, topology_source=source, label=label,
            )
            samples.append(si)
            id_map_lines.append(
                f"{internal_id}\t{sample_name}\t{orig_token}\t{topo}\t"
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


def collect_genbank_files(src: Path) -> List[Path]:
    if src.is_file():
        return [src.resolve()]
    if src.is_dir():
        return sorted(
            p.resolve() for p in src.rglob("*")
            if p.is_file() and _is_gb_file(p)
        )
    return []


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
