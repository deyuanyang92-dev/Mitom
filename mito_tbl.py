#!/usr/bin/env python3
"""mito_tbl.py - GenBank (.gb .gbk .gbf .gbff .genbank) -> NCBI submission files: 5-column feature table (.tbl) and
FASTA with source modifiers (.fsa), the pair BankIt and table2asn take (`mitom gb2tbl`).

Input: one or more GenBank files or folders (every GenBank file below a folder), any mix; MitoM's final_annotation,
MitoZ / MITOS2 / Geneious output, NCBI downloads.
Output (-o OUT):
  --split all (default)  OUT/<name>.tbl + OUT/<name>.fsa with every record (one batch submission); <name>: --name,
                         else the stem of a single input file, else 'submission'
  --split file           one pair per input file: OUT/<file stem>.tbl / .fsa
  --split record         one pair per record:     OUT/<SeqID>.tbl / .fsa
  gb2tbl_report.tsv      every record: its SeqID, source file, output pair, organism, length, topology, location,
                         genetic code, gene counts, warnings (and the NCBI ERROR / WARNING counts with --validate)
Records:
  - the same record in two inputs (same ID and sequence) is written once; different records with one ID are renamed
    ID_2, ID_3 ... (NCBI needs unique SeqIDs in a submission); an ID is cleaned to NCBI's characters
    (A-Z a-z 0-9 _ - . : * #) and cut to 50
  - the FASTA header carries [organism=], [location=mitochondrion] / [location=chloroplast] (from /organelle; none
    for a nuclear record), [topology=], [completeness=complete] for a circular record, [mgcode=] / [gcode=] (the CDS
    transl_table), [taxid=] (/db_xref="taxon:..."), and the source qualifiers NCBI takes as modifiers: isolate, strain, specimen_voucher, country /
    geo_loc_name, collection_date, lat_lon, collected_by, identified_by, haplotype, host, sex, tissue_type, dev_stage,
    note ... (--source_modifiers no: only the first ones)
  - the feature table: every feature but source, partial ends < >, joins, qualifiers (not /translation, /db_xref,
    /protein_id, MitoM's internal notes) - the table `mitom run` writes to ncbi_submission/
--validate: table2asn (found by mito_env.py) on every pair -> <pair>.sqn (the file to upload; --sbt adds the submitter
template from https://submit.ncbi.nlm.nih.gov/genbank/template/submission/) and <pair>.val; NCBI ERROR / WARNING per
record in the report.

Usage:
  mitom gb2tbl -i final_annotation/ -o ncbi                      # all records -> ncbi/final_annotation.tbl + .fsa
  mitom gb2tbl -i a.gbf b.gb c.gbk -o ncbi --name batch1          # three files -> ncbi/batch1.tbl + .fsa
  mitom gb2tbl -i genomes/ -o ncbi --split record --validate --sbt template.sbt
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402
import ncbi_validate as nv  # noqa: E402  (the feature table and the FASTA header of the NCBI submission files)
from mito_common import GB_SUFFIXES, iter_gb_files, run_tool, write_tsv  # noqa: E402

MAX_ID = 50
# GenBank source qualifier -> table2asn FASTA modifier (table2asn 1.28 checked; /country is NCBI's geo_loc_name now)
MODIFIERS = [("isolate", "isolate"), ("strain", "strain"), ("specimen_voucher", "specimen-voucher"),
             ("culture_collection", "culture-collection"), ("bio_material", "bio-material"),
             ("geo_loc_name", "geo_loc_name"), ("country", "geo_loc_name"), ("collection_date", "collection-date"),
             ("lat_lon", "lat-lon"), ("collected_by", "collected-by"), ("identified_by", "identified-by"),
             ("haplotype", "haplotype"), ("host", "host"), ("sex", "sex"), ("tissue_type", "tissue-type"),
             ("dev_stage", "dev-stage"), ("cultivar", "cultivar"), ("breed", "breed"), ("ecotype", "ecotype"),
             ("note", "note")]
COLUMNS = ["seq_id", "record", "source_file", "output", "organism", "length", "topology", "location", "genetic_code",
           "genes", "CDS", "tRNA", "rRNA", "NCBI_ERROR", "NCBI_WARNING", "warnings"]


def clean_id(rid: str) -> str:
    """A SeqID NCBI accepts: A-Z a-z 0-9 _ - . : * # only, at most MAX_ID characters."""
    return re.sub(r"[^A-Za-z0-9_.:*#-]", "_", rid or "record")[:MAX_ID] or "record"


def source_of(rec) -> Dict[str, List[str]]:
    return next((f.qualifiers for f in rec.features if f.type == "source"), {})


def location(rec) -> Tuple[str, str]:
    """(FASTA location modifier or '', genetic code modifier: mgcode for an organelle, gcode for a nuclear record)."""
    org = " ".join(source_of(rec).get("organelle", [])).lower()
    if "mitochondr" in org:
        return "mitochondrion", "mgcode"
    if "chloroplast" in org or "plastid" in org:
        return "chloroplast", "gcode"
    return "", "gcode"


def defline(rec, seq_id: str, modifiers: bool) -> str:
    """FASTA header: [organism=] [location=] [topology=] [completeness=] [mgcode= / gcode=] + source modifiers."""
    loc, code_mod = location(rec)
    src = source_of(rec)
    org = (rec.annotations.get("organism") or " ".join(src.get("organism", [])) or "unknown organism").strip()
    topo = rec.annotations.get("topology", "linear")
    mods = [f"[organism={org}]"] + ([f"[location={loc}]"] if loc else []) + [f"[topology={topo}]"]
    if topo == "circular":                        # else SEQ_INST.CompleteCircleProblem (as ncbi_validate.fasta_defline)
        mods.append("[completeness=complete]")
    code = nv.record_gencode(rec)
    if code:
        mods.append(f"[{code_mod}={code}]")
    taxid = next((x.split(":", 1)[1] for x in src.get("db_xref", []) if x.startswith("taxon:")), "")
    if taxid.isdigit():                           # (else table2asn warns SEQ_DESCR.NoTaxonID)
        mods.append(f"[taxid={taxid}]")
    if modifiers:
        seen = set()
        for qual, mod in MODIFIERS:
            for v in src.get(qual, []):
                v = " ".join(str(v).replace("[", "(").replace("]", ")").split())
                if v and (mod, v) not in seen:
                    seen.add((mod, v))
                    mods.append(f"[{mod}={v}]")
    return f">{seq_id} " + " ".join(mods)


def counts(rec) -> Counter:
    return Counter(f.type for f in rec.features if f.type in ("gene", "CDS", "tRNA", "rRNA"))


def collect(paths: List[Path]) -> Tuple[List[Tuple[Path, object]], List[str]]:
    """(file, record) of every GenBank record of the inputs, files in the order given (a folder: sorted)."""
    out, missing = [], []
    for p in paths:
        files = iter_gb_files(p)
        if not files:
            missing.append(str(p))
        for f in files:
            out += [(f, r) for r in SeqIO.parse(str(f), "genbank")]
    return out, missing


def assign_ids(items) -> List[Tuple[Path, object, str, List[str]]]:
    """(file, record, SeqID, warnings); the same record twice (ID + sequence) is kept once, other ID clashes renamed."""
    kept, used, seqs = [], Counter(), {}
    for f, rec in items:
        warn = []
        raw = rec.id if rec.id and rec.id not in (".", "<unknown id>") else (rec.name or f.stem)
        sid = clean_id(raw)
        if sid != raw:
            warn.append(f"ID {raw!r} written as {sid!r}")
        if seqs.get(sid) == str(rec.seq).upper():
            continue                                   # the same record in another input: once
        if sid in used:
            used[sid] += 1
            new = clean_id(f"{sid[:MAX_ID - 4]}_{used[sid]}")
            warn.append(f"ID {sid} also used by another record: written as {new}")
            sid = new
        used[sid] += 1
        seqs[sid] = str(rec.seq).upper()
        if not (rec.annotations.get("organism") or source_of(rec).get("organism")):
            warn.append("no organism: [organism=unknown organism] - set it before submitting")
        if not any(x.type != "source" for x in rec.features):
            warn.append("no features: the table has the header line only")
        kept.append((f, rec, sid, warn))
    return kept


def pair(out: Path, name: str, ext: str) -> Path:
    """OUT/<name><ext> (not Path.with_suffix: a name such as NC_000931.1 keeps its version)."""
    return out / f"{name}{ext}"


def write_pair(out: Path, name: str, items, modifiers: bool) -> None:
    """<name>.tbl + <name>.fsa with the records of items (file, record, SeqID, warnings)."""
    tbl, fsa = [], []
    for _, rec, sid, _ in items:
        lines = nv.feature_lines(rec)
        lines[0] = f">Feature {sid}"
        tbl += lines
        fsa.append(defline(rec, sid, modifiers) + "\n" + str(rec.seq).upper())
    pair(out, name, ".tbl").write_text("\n".join(tbl) + "\n")
    pair(out, name, ".fsa").write_text("\n".join(fsa) + "\n")


def validate(out: Path, name: str, t2a: str, sbt: Optional[str]) -> Dict[str, Counter]:
    """table2asn on <name>.fsa (+ <name>.tbl) -> <name>.sqn, <name>.val; ERROR / WARNING per SeqID (submitter /
    publication messages left out, as in ncbi_validate.py)."""
    run_tool([t2a, "-i", str(pair(out, name, ".fsa")), "-V", "vb", "-o", str(pair(out, name, ".sqn")),
              *(["-t", sbt] if sbt else [])])
    per: Dict[str, Counter] = defaultdict(Counter)
    val = pair(out, name, ".val")
    for line in val.read_text(errors="ignore").splitlines() if val.is_file() else []:
        m = nv.MSG_RE.match(line.strip())
        if not m or nv.IGNORE_CODES.search(m.group(2)):
            continue
        sid = re.search(r"lcl\|([^\s:\]]+)", m.group(3))
        per[sid.group(1) if sid else ""][m.group(1).upper()] += 1
    return per


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="mitom gb2tbl", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", nargs="+", action="extend", required=True,
                   help="GenBank files and / or folders (repeatable: -i a.gbf b.gb -i more/)")
    p.add_argument("-o", "--out", required=True, help="output folder")
    p.add_argument("--split", choices=["all", "file", "record"], default="all",
                   help="all (default): one .tbl + .fsa with every record; file: one pair per input file; "
                        "record: one pair per record")
    p.add_argument("--name", default=None,
                   help="base name of the pair with --split all (default: the stem of a single input file, else "
                        "'submission')")
    p.add_argument("--source_modifiers", choices=["yes", "no"], default="yes",
                   help="copy the source qualifiers NCBI takes as FASTA modifiers (isolate, specimen_voucher, "
                        "country, collection_date ...) into the .fsa headers (default yes)")
    p.add_argument("--validate", action="store_true",
                   help="run NCBI table2asn on every pair: .sqn (the file to upload) and .val, errors in the report")
    p.add_argument("--sbt", default=None, help="NCBI submitter template (.sbt) for the .sqn files (with --validate)")
    p.add_argument("--table2asn", default=None, help="table2asn executable (default: found by mito_env.py)")
    a = p.parse_args(argv)
    inputs = [Path(x) for x in a.input]
    items, missing = collect(inputs)
    for m in missing:
        print(f"[gb2tbl] no GenBank file in {m} (suffixes {', '.join(GB_SUFFIXES)})", file=sys.stderr)
    if not items:
        return 2
    kept = assign_ids(items)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    groups: Dict[str, list] = defaultdict(list)
    for it in kept:
        if a.split == "record":
            groups[it[2]].append(it)
        elif a.split == "file":
            groups[clean_id(it[0].stem)].append(it)
        else:
            single = len(inputs) == 1 and inputs[0].is_file()
            groups[a.name or (clean_id(inputs[0].stem) if single else "submission")].append(it)
    t2a = mito_env.find("table2asn", a.table2asn) if a.validate else None
    if a.validate and not t2a:
        print("[gb2tbl] table2asn not found (mitom doctor; conda install -c bioconda table2asn): no .sqn / check",
              file=sys.stderr)
    rows = []
    for name, its in groups.items():
        write_pair(out, name, its, a.source_modifiers == "yes")
        val = validate(out, name, t2a, a.sbt) if t2a else {}
        for f, rec, sid, warn in its:
            c = counts(rec)
            loc, _ = location(rec)
            rows.append({"seq_id": sid, "record": rec.id, "source_file": str(f), "output": f"{name}.tbl",
                         "organism": rec.annotations.get("organism", ""), "length": len(rec.seq),
                         "topology": rec.annotations.get("topology", "linear"), "location": loc or "nucleus",
                         "genetic_code": nv.record_gencode(rec) or "", "genes": c["gene"], "CDS": c["CDS"],
                         "tRNA": c["tRNA"], "rRNA": c["rRNA"],
                         "NCBI_ERROR": val.get(sid, Counter())["ERROR"] if t2a else "",
                         "NCBI_WARNING": val.get(sid, Counter())["WARNING"] if t2a else "", "warnings": "; ".join(warn)})
    write_tsv(out / "gb2tbl_report.tsv", COLUMNS, rows)
    dropped = len(items) - len(kept)
    print(f"[gb2tbl] {len(kept)} record(s) from {len(items)} in {len({str(f) for f, _ in items})} file(s)"
          + (f" ({dropped} repeated record(s) written once)" if dropped else "")
          + f" -> {len(groups)} .tbl / .fsa pair(s) in {out}", file=sys.stderr)
    renamed = sum("also used" in r["warnings"] for r in rows)
    if renamed:
        print(f"[gb2tbl] {renamed} record(s) renamed: their ID was used by another record (gb2tbl_report.tsv)",
              file=sys.stderr)
    if t2a:
        e = sum(int(r["NCBI_ERROR"] or 0) for r in rows)
        print(f"[gb2tbl] table2asn: {e} ERROR(s) in {sum(1 for r in rows if r['NCBI_ERROR'])} record(s); .sqn files "
              "to upload" + ("" if a.sbt else " (add --sbt template.sbt for the submitter)"), file=sys.stderr)
        if e:
            print("[gb2tbl] the messages: <pair>.val; most annotation errors (missing stop codons, overlaps) are "
                  "fixed by `mitom curate -i <GenBank> -o curated` (or `mitom import` for MitoZ / MITOS2 results), "
                  "then gb2tbl again", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
