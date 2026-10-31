#!/usr/bin/env python3
"""ncbi_validate.py - check GenBank annotation against NCBI's own submission validator.

Each record is converted to the NCBI submission inputs (.fsa with source modifiers and
5-column feature table .tbl) and validated with table2asn -V vb, the program NCBI uses
when it processes a submission. Messages are summarised per severity and code.
Submitter/publication messages (no template) are ignored: only sequence and feature
problems (SEQ_FEAT, SEQ_INST, SEQ_DESCR on the record) count.

Usage:
  python3 ncbi_validate.py -i gb_file_or_dir -o out_dir [--table2asn path] [--focus intron]
Outputs: out_dir/<record>/{.fsa,.tbl,.val}, out_dir/validation.tsv (one row per message),
and a summary on stdout.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.SeqFeature import AfterPosition, BeforePosition

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402  (finds external programs: config, PATH, conda environments)
from mito_common import iter_gb_files, run_tool, write_tsv  # noqa: E402
from mito_common import cli_guard, InputError, prepare_input, require_records  # noqa: E402

SKIP_QUALS = {"translation", "db_xref", "protein_id", "locus_tag_old"}
INTERNAL_NOTE = re.compile(r"check .* manually|mito_\w+\.py")   # e.g. mito_rescue.py / mito_intron.py notes
IGNORE_CODES =re.compile(r"SEQ_DESCR\.(NoPubFound|NoSourceDescriptor|MissingLineage|NoTaxonID|BadOrgMod)|"
                          r"GENERIC\.MissingPubRequirement|SEQ_PKG\.|SEQ_DESCR\.TitleNotAppropriateForSet")
MSG_RE = re.compile(r"^(INFO|WARNING|ERROR|REJECT|FATAL):\s+valid\s+\[([A-Z_]+\.[A-Za-z0-9_]+)\]\s+(.*)$", re.I)


def feature_lines(rec) -> List[str]:
    out = [f">Feature {rec.id}"]
    for f in rec.features:
        if f.type == "source":
            continue
        parts = list(f.location.parts)               # Biopython keeps minus-strand parts in transcript order
        for i, p in enumerate(parts):
            s, e = int(p.start) + 1, int(p.end)
            lo = "<" if isinstance(p.start, BeforePosition) else ""
            hi = ">" if isinstance(p.end, AfterPosition) else ""
            if p.strand == -1:      # '<' always on the first column, '>' on the second (table2asn: '>150<TAB>40' = Bad interval)
                a, b = ("<" if hi else "") + str(e), (">" if lo else "") + str(s)
            else:
                a, b = f"{lo}{s}", f"{hi}{e}"
            out.append(f"{a}\t{b}" + (f"\t{f.type}" if i == 0 else ""))
        for k, vals in f.qualifiers.items():
            if k in SKIP_QUALS:
                continue
            for v in vals:
                if k == "note" and INTERNAL_NOTE.search(str(v)):     # tool hints for the curator, not for the database
                    continue
                out.append(f"\t\t\t{k}\t{' '.join(str(v).split())}")
    return out


def fasta_defline(rec, gencode: Optional[int], nuclear: bool = False) -> str:
    org = rec.annotations.get("organism", "unknown organism")
    topo = rec.annotations.get("topology", "linear")
    mods = [f"[organism={org}]"] + ([] if nuclear else ["[location=mitochondrion]"]) + [f"[topology={topo}]"]
    if topo == "circular":                       # else SEQ_INST.CompleteCircleProblem (table2asn 2026-10-21)
        mods.append("[completeness=complete]")
    if gencode:
        mods.append(f"[{'gcode' if nuclear else 'mgcode'}={gencode}]")
    return f">{rec.id} " + " ".join(mods)


def record_gencode(rec) -> Optional[int]:
    codes = [int(f.qualifiers["transl_table"][0]) for f in rec.features
             if f.type == "CDS" and f.qualifiers.get("transl_table", [""])[0].isdigit()]
    return max(set(codes), key=codes.count) if codes else None


# notes of the features the intron steps write: mito_intron.py ('... MFannot method'), mito_intron_ref.py (CDS: 'known
# intron sites (...; mito_intron_ref)', introns: 'site transferred from <accession>'), older runs ('mito_intron')
INTRON_STEP_NOTE = re.compile(r"MFannot method|mito_intron|site transferred from")


def our_intervals(rec, tag=INTRON_STEP_NOTE) -> List[str]:
    """'a-b' strings of every interval of features whose /note matches tag (a regex or a text; default: the features
    written by MitoM's intron steps)."""
    rx = re.compile(re.escape(tag)) if isinstance(tag, str) else tag
    out = []
    for f in rec.features:
        if any(rx.search(n) for n in f.qualifiers.get("note", [])):
            out += [f"{int(p.start) + 1}-{int(p.end)}" for p in f.location.parts]
    return out


def record_dir(out: Path, rec_id: str) -> Path:
    """Folder of one record's submission files (rec.fsa, rec.tbl, rec.sqn, rec.val) below `out`."""
    return out / re.sub(r"[^\w.]+", "_", rec_id)


TABLE2ASN_TIMEOUT = 600      # s; one mitogenome takes < 1 s


def validate_record(rec, out: Path, table2asn: str, sbt: Optional[str] = None,
                    nuclear: bool = False) -> List[Dict[str, str]]:
    """Write rec as NCBI submission inputs (rec.fsa + rec.tbl) and run table2asn on them (rec.sqn, rec.val);
    sbt = submitter template from https://submit.ncbi.nlm.nih.gov/genbank/template/submission/ (optional)."""
    d = record_dir(out, rec.id)
    d.mkdir(parents=True, exist_ok=True)
    base = d / "rec"
    (d / "rec.fsa").write_text(fasta_defline(rec, record_gencode(rec), nuclear) + "\n" + str(rec.seq) + "\n")
    (d / "rec.tbl").write_text("\n".join(feature_lines(rec)) + "\n")
    rows = []
    for old in (d / "rec.val", d / "rec.sqn"):       # a failing table2asn must not leave the last run's messages
        old.unlink(missing_ok=True)
    try:
        run_tool([table2asn, "-i", str(d / "rec.fsa"), "-V", "vb", "-o", str(base) + ".sqn",
                  *(["-t", sbt] if sbt else [])], timeout=TABLE2ASN_TIMEOUT)
    except subprocess.TimeoutExpired:                  # never wait for ever (table2asn 1.28 can hang)
        return [{"record": rec.id, "severity": "ERROR", "code": "MITOM.Table2asnTimeout", "ours": "no",
                 "message": f"table2asn did not finish in {TABLE2ASN_TIMEOUT} s: not validated"}]
    val = d / "rec.val"
    mine = our_intervals(rec)
    for line in val.read_text(errors="ignore").splitlines() if val.exists() else []:
        m = MSG_RE.match(line.strip())
        if not m or IGNORE_CODES.search(m.group(2)):
            continue
        msg = m.group(3)
        rows.append({"record": rec.id, "severity": m.group(1).upper(), "code": m.group(2),
                     "ours": "yes" if any(re.search(rf"\b{re.escape(iv)}\b", msg) for iv in mine) else "no",
                     "message": msg[:300]})
    return rows


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--table2asn", default=None, help="table2asn executable (default: found by mito_env.py)")
    p.add_argument("--focus", default="", help="Only count messages whose text matches this regex (e.g. intron|CDS)")
    p.add_argument("--only_ours", action="store_true",
                   help="Only count messages on features written by MitoM's intron steps (mito_intron.py, mito_intron_ref.py)")
    p.add_argument("--nuclear", action="store_true", help="nuclear sequences (mitom nuc): no location, [gcode=]")
    a = p.parse_args(argv)
    exe = mito_env.find("table2asn", a.table2asn)
    if not exe:
        print("[ncbi_validate] table2asn not found (conda install -c bioconda table2asn)", file=sys.stderr)
        return 1
    src = Path(prepare_input(a.input))
    files = iter_gb_files(src) if src.is_dir() else [src]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, str]] = []
    n = 0
    for f in files:
        for rec in SeqIO.parse(str(f), "genbank"):
            n += 1
            rows += validate_record(rec, out, exe, nuclear=a.nuclear)
    require_records(n, "GenBank records", a.input)
    write_tsv(out / "validation.tsv", ["record", "severity", "code", "ours", "message"], rows)
    focus = re.compile(a.focus, re.I) if a.focus else None
    sel = [r for r in rows if (not focus or focus.search(r["message"])) and (not a.only_ours or r["ours"] == "yes")]
    sev = Counter(r["severity"] for r in sel)
    print(f"[ncbi_validate] {n} record(s); messages: " + ", ".join(f"{k} {v}" for k, v in sev.most_common()),
          file=sys.stderr)
    for (s, c), k in Counter((r["severity"], r["code"]) for r in sel if r["severity"] in ("ERROR", "REJECT", "FATAL")).most_common(15):
        print(f"   {s:6s} {c:40s} {k}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(cli_guard(main))
