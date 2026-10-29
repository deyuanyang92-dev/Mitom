#!/usr/bin/env python3
"""mito_check.py - audit the existing annotation of mitochondrial genome GenBank records.

A report only: it shows what is wrong with an existing annotation (e.g. records downloaded from NCBI, or a
final MitoM result); it does not change any file. `mitom run` re-annotates every record regardless.
An ERROR means the annotation is unreliable; a WARNING is a smaller problem or legitimate variation.

ERROR   no_annotation      no CDS feature at all
        missing_pcg        one of the 13 protein genes absent (cox1-3, cob, atp6, atp8, nad1-6, nad4l)
        missing_rrna       rrnS (12S) or rrnL (16S) absent
        duplicate_pcg      the same protein gene annotated twice (split genes COX1a/COX1b are a warning)
        internal_stop      stop codon inside a CDS (translated with its transl_table)
        frame              CDS length not a multiple of 3 without partial marker or transl_except
        pcg_length         PCG shorter than 60 % or longer than 160 % of the usual length
        pcg_overlap        two protein/rRNA genes overlap by more than 100 nt
        trna_count         fewer tRNAs than --min_trna (22 in most bilaterians)
WARNING no_transl_table    CDS without /transl_table
        start_codon        first codon is not a start codon of the table
        no_stop            last codon is not a stop and no partial marker / transl_except
        split_gene         gene annotated in two parts (COX1a / COX1b)
        trna_set           amino acids without a tRNA, or extra copies
        orientation        cox1 not on the plus strand, or not the first gene (the run re-orients)
        unknown_cds        a CDS that is not one of the 13 genes (ORF ...)

Usage:
  mitom check -i genomes.gb -o check_out [--genetic_code 5] [--min_trna 22]
Outputs in -o:  check_report.tsv (one row per record)  check_details.tsv (one row per issue)
                check_flagged.txt (records with an ERROR)
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from mito_common import iter_gb_files, write_tsv  # noqa: E402
from mito_names import PCG, gene_name, trna_name  # noqa: E402

# usual protein gene lengths (nt, incl. stop) in bilaterian mitogenomes; 60 % / 160 % are the error limits
TYPICAL_LEN = {"cox1": 1536, "cox2": 690, "cox3": 780, "cob": 1140, "nad1": 930, "nad2": 1000, "nad3": 350,
               "nad4": 1350, "nad4l": 290, "nad5": 1750, "nad6": 480, "atp6": 700, "atp8": 160}
LEN_LOW, LEN_HIGH = 0.6, 1.6
MAX_OVERLAP = 100
AA20 = set("ACDEFGHIKLMNPQRSTVWY")
TWO_COPIES = {"L", "S"}


@dataclass
class Issue:
    severity: str      # ERROR / WARNING
    code: str
    feature: str       # gene name or ''
    location: str
    message: str


def _loc(f: SeqFeature) -> str:
    return str(f.location)


def _qual(f: SeqFeature, k: str) -> str:
    return (f.qualifiers.get(k, [""]) or [""])[0]


def _partial(f: SeqFeature) -> Tuple[bool, bool]:
    """(5' partial, 3' partial) from fuzzy positions, in biological orientation."""
    s = str(f.location)
    left, right = s.startswith("<") or "[<" in s, ">" in s
    return (right, left) if f.location.strand == -1 else (left, right)


def cds_table(f: SeqFeature, default: Optional[int]) -> Optional[int]:
    t = _qual(f, "transl_table")
    if t.isdigit():
        return int(t)
    return default


def record_code(rec: SeqRecord, default: Optional[int]) -> Optional[int]:
    """Genetic code of a record: the most common /transl_table of its CDS, else default."""
    c = Counter(_qual(f, "transl_table") for f in rec.features if f.type == "CDS" and _qual(f, "transl_table").isdigit())
    return int(c.most_common(1)[0][0]) if c else default


def check_cds(rec: SeqRecord, f: SeqFeature, name: str, code: Optional[int]) -> List[Issue]:
    out: List[Issue] = []
    label = name or _qual(f, "gene") or _qual(f, "product") or "CDS"
    table_id = cds_table(f, code)
    if not _qual(f, "transl_table"):
        out.append(Issue("WARNING", "no_transl_table", label, _loc(f), "CDS without /transl_table"))
    table = None                                   # without a genetic code only the frame and the length are checked
    if table_id is not None:
        try:
            table = CodonTable.unambiguous_dna_by_id[table_id]
        except KeyError:
            out.append(Issue("WARNING", "no_transl_table", label, _loc(f), f"unknown transl_table {table_id}"))
    try:
        nt = str(f.extract(rec.seq)).upper()
    except Exception as exc:  # noqa: BLE001  (location outside the sequence etc.)
        out.append(Issue("ERROR", "frame", label, _loc(f), f"cannot extract the CDS: {exc}"))
        return out
    p5, p3 = _partial(f)
    cs = _qual(f, "codon_start")
    offset = int(cs) - 1 if cs in ("1", "2", "3") else 0
    nt = nt[offset:]
    abbreviated = bool(f.qualifiers.get("transl_except"))      # T / TA stop completed by polyadenylation
    if len(nt) % 3 and not (p5 or p3 or abbreviated):
        out.append(Issue("ERROR", "frame", label, _loc(f), f"length {len(nt) + offset} is not a multiple of 3"))
    codons = [nt[i:i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]
    if not codons:
        out.append(Issue("ERROR", "frame", label, _loc(f), "empty CDS"))
        return out
    if table is not None:
        stops = set(table.stop_codons)
        internal = [i + 1 for i, c in enumerate(codons[:-1]) if c in stops]
        if internal:
            out.append(Issue("ERROR", "internal_stop", label, _loc(f),
                             f"{len(internal)} internal stop codon(s) at codon {internal[:5]} (table {table_id})"))
        if not p5 and codons[0] not in table.start_codons and "N" not in codons[0]:
            out.append(Issue("WARNING", "start_codon", label, _loc(f), f"starts with {codons[0]} (table {table_id})"))
        last = codons[-1] if len(nt) % 3 == 0 else nt[-(len(nt) % 3):]
        if not p3 and last not in stops and not abbreviated:
            out.append(Issue("WARNING", "no_stop", label, _loc(f), f"ends with {last}, no stop codon / transl_except"))
    if name in TYPICAL_LEN:
        n = len(nt) + offset
        lo, hi = TYPICAL_LEN[name] * LEN_LOW, TYPICAL_LEN[name] * LEN_HIGH
        if n < lo:
            out.append(Issue("ERROR", "pcg_length", label, _loc(f), f"{n} nt, shorter than {lo:.0f} (usual {TYPICAL_LEN[name]})"))
        elif n > hi:
            out.append(Issue("ERROR", "pcg_length", label, _loc(f), f"{n} nt, longer than {hi:.0f} (usual {TYPICAL_LEN[name]})"))
    return out


def _segments(f: SeqFeature) -> List[Tuple[int, int]]:
    return [(int(p.start), int(p.end)) for p in f.location.parts]


def overlap_nt(a: SeqFeature, b: SeqFeature) -> int:
    return sum(max(0, min(e1, e2) - max(s1, s2)) for s1, e1 in _segments(a) for s2, e2 in _segments(b))


def check_record(rec: SeqRecord, code: Optional[int] = None, min_trna: int = 22) -> List[Issue]:
    """All issues of one record (ERROR first)."""
    issues: List[Issue] = []
    cds = [f for f in rec.features if f.type == "CDS"]
    rrna = [f for f in rec.features if f.type == "rRNA"]
    trna = [f for f in rec.features if f.type == "tRNA"]
    if not cds:
        issues.append(Issue("ERROR", "no_annotation", "", "", "no CDS feature"))
        return issues
    code = record_code(rec, code)

    named: Dict[str, List[SeqFeature]] = {}
    for f in cds:
        g = gene_name(f.qualifiers, "CDS")
        if g:
            named.setdefault(g, []).append(f)
        else:
            issues.append(Issue("WARNING", "unknown_cds", _qual(f, "gene") or _qual(f, "product"), _loc(f),
                                "CDS that is not one of the 13 protein genes"))
        issues += check_cds(rec, f, g, code)
    for g in sorted(PCG):
        if g not in named:
            issues.append(Issue("ERROR", "missing_pcg", g, "", "protein gene not annotated"))
    for g, fs in named.items():
        if len(fs) > 1:
            raw = {_qual(f, "gene").lower() for f in fs}
            if len(raw) == len(fs) and all(re.search(r"[abc]$", r) for r in raw):
                issues.append(Issue("WARNING", "split_gene", g, "; ".join(_loc(f) for f in fs),
                                    f"annotated in {len(fs)} parts ({', '.join(sorted(raw))})"))
            else:
                issues.append(Issue("ERROR", "duplicate_pcg", g, "; ".join(_loc(f) for f in fs),
                                    f"annotated {len(fs)} times"))

    rr = {gene_name(f.qualifiers, "rRNA") for f in rrna}
    for g in ("rrnS", "rrnL"):
        if g not in rr:
            issues.append(Issue("ERROR", "missing_rrna", g, "", "rRNA not annotated"))

    big = [(gene_name(f.qualifiers, f.type) or _qual(f, "gene"), f) for f in cds + rrna]
    for i, (na, a) in enumerate(big):
        for nb, b in big[i + 1:]:
            o = overlap_nt(a, b)
            if o > MAX_OVERLAP:
                issues.append(Issue("ERROR", "pcg_overlap", f"{na}/{nb}", f"{_loc(a)} ; {_loc(b)}", f"overlap {o} nt"))

    if len(trna) < min_trna:
        issues.append(Issue("ERROR", "trna_count", "", "", f"{len(trna)} tRNA(s), fewer than {min_trna}"))
    aas = Counter(trna_name(f.qualifiers)[0] for f in trna)
    missing = sorted(a for a in AA20 if aas.get(a, 0) < (2 if a in TWO_COPIES else 1))
    extra = sorted(a for a, n in aas.items() if a and n > (2 if a in TWO_COPIES else 1))
    if missing or extra:
        issues.append(Issue("WARNING", "trna_set", "", "",
                            (f"no tRNA for {','.join(missing)}" if missing else "")
                            + ("; " if missing and extra else "") + (f"extra copies of {','.join(extra)}" if extra else "")))

    cox1 = named.get("cox1")
    if cox1:
        f = cox1[0]
        first = min(cds + rrna + trna, key=lambda x: int(x.location.start))
        if f.location.strand == -1:
            issues.append(Issue("WARNING", "orientation", "cox1", _loc(f), "cox1 on the minus strand"))
        elif first is not f and first.type != "tRNA":
            issues.append(Issue("WARNING", "orientation", "cox1", _loc(f),
                                f"cox1 is not the first gene ({_qual(first, 'gene') or first.type} starts at {int(first.location.start) + 1})"))
    issues.sort(key=lambda i: (i.severity != "ERROR", i.code, i.feature))
    return issues


def counts(rec: SeqRecord) -> Tuple[int, int, int]:
    return (sum(f.type == "CDS" for f in rec.features), sum(f.type == "tRNA" for f in rec.features),
            sum(f.type == "rRNA" for f in rec.features))


def gb_files(path: Path) -> List[Path]:
    return iter_gb_files(path)


def check_files(paths: Iterable[Path], code: Optional[int] = None, min_trna: int = 22) -> List[dict]:
    rows = []
    for p in paths:
        for rec in SeqIO.parse(str(p), "genbank"):
            issues = check_record(rec, code, min_trna)
            errs = [i for i in issues if i.severity == "ERROR"]
            warns = [i for i in issues if i.severity == "WARNING"]
            pcg, trn, rrn = counts(rec)
            rows.append({
                "record": rec.id, "organism": rec.annotations.get("organism", ""),
                "status": "flagged" if errs else "ok", "n_error": len(errs), "n_warning": len(warns),
                "errors": "; ".join(f"{i.code}:{i.feature}" if i.feature else i.code for i in errs),
                "warnings": "; ".join(f"{i.code}:{i.feature}" if i.feature else i.code for i in warns),
                "PCG": pcg, "tRNA": trn, "rRNA": rrn, "length": len(rec.seq), "file": p.name,
                "_issues": issues,
            })
    return rows


def write_reports(rows: List[dict], out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    cols = ["record", "organism", "status", "n_error", "n_warning", "errors", "warnings", "PCG", "tRNA", "rRNA",
            "length", "file"]
    write_tsv(out / "check_report.tsv", cols, rows)
    write_tsv(out / "check_details.tsv", ["record", "severity", "code", "feature", "location", "message"],
              [[r["record"], i.severity, i.code, i.feature, i.location, i.message] for r in rows for i in r["_issues"]])
    (out / "check_flagged.txt").write_text("".join(f"{r['record']}\n" for r in rows if r["status"] == "flagged"),
                                           encoding="utf-8")


def summary_line(rows: List[dict]) -> str:
    n = len(rows)
    flagged = sum(r["status"] == "flagged" for r in rows)
    codes = Counter(i.code for r in rows for i in r["_issues"] if i.severity == "ERROR")
    top = ", ".join(f"{c} {k}" for k, c in codes.most_common())
    return f"[CHECK] {n} record(s): {n - flagged} ok, {flagged} flagged for re-annotation" + (f" ({top})" if top else "")


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="GenBank file or directory.")
    ap.add_argument("-o", "--out_dir", default="check_out", help="Output directory (default: check_out).")
    ap.add_argument("--genetic_code", type=int, default=None,
                    help="Translation table for CDS without /transl_table (default: the record's own tables).")
    ap.add_argument("--min_trna", type=int, default=22, help="Fewer tRNAs than this is an ERROR (default 22; 0 = off).")
    args = ap.parse_args(argv)
    files = gb_files(Path(args.input))
    if not files:
        print(f"no GenBank files in {args.input}", file=sys.stderr)
        return 2
    rows = check_files(files, args.genetic_code, args.min_trna)
    write_reports(rows, Path(args.out_dir))
    print(summary_line(rows), file=sys.stderr)
    print(f"reports: {Path(args.out_dir) / 'check_report.tsv'}, check_details.tsv, check_flagged.txt", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
