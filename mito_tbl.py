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
Cleaning (what other programs leave in a GenBank file and NCBI does not take; every change in the report):
  - qualifiers outside the INSDC feature table definition (Geneious /created_by /modified_by, SnapGene / ApE /label
    /ApEinfo_*, ...) are left out; a feature key outside it (Geneious 'ORF' ...) becomes misc_feature with the old key
    as a note
  - text is made ASCII (NCBI rejects anything else, SEQ_INST.NonAsciiAsn): accents dropped (Muller, Koln), mu -> u,
    typographic quotes and dashes plain; other characters removed
  - a feature running past the end of the sequence is reported; a site between two bases (100^101) is written as the
    2-nt interval around it (as NCBI's gbf2tbl)
  - a file that is not readable GenBank, and a record without sequence (CONTIG / WGS master records) are reported and
    skipped; the other inputs are written
  - genetic code: the CDS /transl_table; --genetic_code N for records without one (else, for a mitochondrion, the
    usual code of its phylum, reported)
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
import copy
import re
import sys
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402
import ncbi_validate as nv  # noqa: E402  (the feature table and the FASTA header of the NCBI submission files)
from mito_common import GB_SUFFIXES, iter_gb_files, run_tool, write_tsv  # noqa: E402
from mito_taxa import METAZOAN_PHYLA, PHYLUM_CODE, PHYLUM_SYNONYMS  # noqa: E402

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
# INSDC Feature Table Definition (v11.3): feature keys and qualifiers NCBI takes (+ /codon_recognized of older NCBI
# tRNA annotation, which table2asn still takes)
FEATURE_KEYS = {
    "assembly_gap", "C_region", "CDS", "centromere", "D-loop", "D_segment", "exon", "gap", "gene", "iDNA", "intron",
    "J_segment", "mat_peptide", "misc_binding", "misc_difference", "misc_feature", "misc_recomb", "misc_RNA",
    "misc_structure", "mobile_element", "modified_base", "mRNA", "ncRNA", "N_region", "old_sequence", "operon", "oriT",
    "polyA_site", "precursor_RNA", "prim_transcript", "primer_bind", "propeptide", "protein_bind", "regulatory",
    "repeat_region", "rep_origin", "rRNA", "S_region", "sig_peptide", "source", "stem_loop", "STS", "telomere", "tmRNA",
    "transit_peptide", "tRNA", "unsure", "V_region", "V_segment", "variation", "3'UTR", "5'UTR"}
QUALIFIERS = {
    "allele", "altitude", "anticodon", "artificial_location", "bio_material", "bound_moiety", "cell_line", "cell_type",
    "chromosome", "circular_RNA", "citation", "clone", "clone_lib", "codon_recognized", "codon_start", "collected_by",
    "collection_date",
    "compare", "country", "cultivar", "culture_collection", "db_xref", "dev_stage", "direction", "EC_number", "ecotype",
    "environmental_sample", "estimated_length", "exception", "experiment", "focus", "frequency", "function", "gap_type",
    "gene", "gene_synonym", "geo_loc_name", "germline", "haplogroup", "haplotype", "host", "identified_by", "inference",
    "isolate", "isolation_source", "lab_host", "lat_lon", "linkage_evidence", "locus_tag", "macronuclear", "map",
    "mating_type", "metagenome_source", "mobile_element_type", "mod_base", "mol_type", "ncRNA_class", "note", "number",
    "old_locus_tag", "operon", "organelle", "organism", "partial", "PCR_conditions", "PCR_primers", "phenotype",
    "plasmid", "pop_variant", "product", "protein_id", "proviral", "pseudo", "pseudogene", "rearranged",
    "recombination_class", "regulatory_class", "replace", "ribosomal_slippage", "rpt_family", "rpt_type",
    "rpt_unit_range", "rpt_unit_seq", "satellite", "segment", "serotype", "serovar", "sex", "specimen_voucher",
    "standard_name", "strain", "sub_clone", "submitter_seqid", "sub_species", "sub_strain", "tag_peptide", "tissue_lib",
    "tissue_type", "transgenic", "translation", "transl_except", "transl_table", "trans_splicing", "type_material",
    "variety"}
SUBPHYLUM_CODE = {"Tunicata": 13, "Cephalochordata": 5}          # the chordates without the vertebrate code 2
ASCII = {"\u00b5": "u", "\u03bc": "u", "\u00df": "ss", "\u00e6": "ae", "\u00c6": "AE", "\u0153": "oe", "\u00f8": "o",
         "\u00d8": "O", "\u0131": "i", "\u0142": "l", "\u0141": "L", "\u00b0": " degrees", "\u2018": "'", "\u2019": "'",
         "\u201c": '"', "\u201d": '"', "\u2013": "-", "\u2014": "-", "\u2212": "-", "\u00d7": "x", "\u00a0": " "}


def clean_id(rid: str) -> str:
    """A SeqID NCBI accepts: A-Z a-z 0-9 _ - . : * # only, at most MAX_ID characters."""
    return re.sub(r"[^A-Za-z0-9_.:*#-]", "_", rid or "record")[:MAX_ID] or "record"


def ascii_text(v: str) -> str:
    """NCBI takes ASCII only: accents dropped (u-umlaut -> u), ASCII stand-ins for mu, sharp s, quotes, dashes ...;
    other characters (CJK ...) removed."""
    v = "".join(ASCII.get(ch, ch) for ch in str(v))
    v = unicodedata.normalize("NFKD", v).encode("ascii", "ignore").decode("ascii")
    return " ".join(v.split())


def clean_record(rec) -> Tuple[object, List[str]]:
    """A copy of rec NCBI takes (see the module text), and what was changed."""
    rec = copy.deepcopy(rec)
    dropped, keys, texts, beyond, sites = Counter(), Counter(), 0, [], 0
    n = len(rec.seq)
    for f in rec.features:
        loc = f.location
        if loc is not None and len(loc.parts) == 1 and int(loc.start) == int(loc.end) and 0 < int(loc.start) < n:
            # a site between two bases (GenBank 100^101): the 2-nt interval around it, as NCBI's gbf2tbl writes it
            f.location = FeatureLocation(int(loc.start) - 1, int(loc.start) + 1, loc.strand)
            sites += 1
        if f.type not in FEATURE_KEYS:
            keys[f.type] += 1
            f.qualifiers.setdefault("note", []).insert(0, f.type)
            f.type = "misc_feature"
        for k in list(f.qualifiers):
            if k not in QUALIFIERS:
                dropped[k] += 1
                del f.qualifiers[k]
                continue
            vals = [str(v) for v in f.qualifiers[k]]
            new = [ascii_text(v) if not v.isascii() else v for v in vals]
            texts += sum(a != b for a, b in zip(vals, new))
            f.qualifiers[k] = new
        if f.location is not None and int(f.location.end) > n:
            beyond.append(f"{f.type} {int(f.location.start) + 1}..{int(f.location.end)}")
    for k in ("organism", "source"):
        if rec.annotations.get(k) and not str(rec.annotations[k]).isascii():
            rec.annotations[k] = ascii_text(rec.annotations[k])
            texts += 1
    warn = []
    if keys:
        warn.append("feature keys NCBI does not take written as misc_feature: "
                    + ", ".join(f"{k} x{v}" for k, v in keys.most_common()))
    if dropped:
        warn.append("qualifiers NCBI does not take left out: " + ", ".join(f"{k} x{v}" for k, v in dropped.most_common()))
    if texts:
        warn.append(f"{texts} text(s) made ASCII (NCBI takes no other characters)")
    if beyond:
        warn.append(f"feature(s) past the end of the {n} bp sequence: {', '.join(beyond[:3])} - check")
    if sites:
        warn.append(f"{sites} site(s) between two bases (100^101) written as the 2-nt interval around them")
    return rec, warn


def record_code(rec, override: Optional[int]) -> Tuple[Optional[int], str]:
    """(genetic code, warning): --genetic_code, else the CDS /transl_table, else for a mitochondrion the usual code
    of the phylum in the record's taxonomy (5 when unknown; tunicates 13, lancelets 5)."""
    if override:
        return override, ""
    code = nv.record_gencode(rec)
    if code or location(rec)[0] != "mitochondrion" or not any(f.type == "CDS" for f in rec.features):
        return code, ""
    tax = rec.annotations.get("taxonomy", []) or []
    phylum = next((t for t in tax if t in METAZOAN_PHYLA), "") or next((PHYLUM_SYNONYMS[t] for t in tax
                                                                       if t in PHYLUM_SYNONYMS), "")
    code = next((SUBPHYLUM_CODE[t] for t in tax if t in SUBPHYLUM_CODE), PHYLUM_CODE.get(phylum, 5))
    return code, (f"no /transl_table: genetic code {code} assumed ({phylum or 'phylum unknown'}); "
                  "give --genetic_code to set it")


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


def defline(rec, seq_id: str, modifiers: bool, code: Optional[int] = None) -> str:
    """FASTA header: [organism=] [location=] [topology=] [completeness=] [mgcode= / gcode=] + source modifiers;
    code: the genetic code (record_code), default the CDS /transl_table."""
    loc, code_mod = location(rec)
    src = source_of(rec)
    org = (rec.annotations.get("organism") or " ".join(src.get("organism", [])) or "unknown organism").strip()
    topo = rec.annotations.get("topology", "linear")
    mods = [f"[organism={org}]"] + ([f"[location={loc}]"] if loc else []) + [f"[topology={topo}]"]
    if topo == "circular":                        # else SEQ_INST.CompleteCircleProblem (as ncbi_validate.fasta_defline)
        mods.append("[completeness=complete]")
    code = code or nv.record_gencode(rec)
    if code:
        mods.append(f"[{code_mod}={code}]")
    taxid = next((x.split(":", 1)[1] for x in src.get("db_xref", []) if x.startswith("taxon:")), "")
    if taxid.isdigit():                           # (else table2asn warns SEQ_DESCR.NoTaxonID)
        mods.append(f"[taxid={taxid}]")
    if modifiers:
        seen = set()
        for qual, mod in MODIFIERS:
            for v in src.get(qual, []):
                v = ascii_text(str(v).replace("[", "(").replace("]", ")"))
                if v and (mod, v) not in seen:
                    seen.add((mod, v))
                    mods.append(f"[{mod}={v}]")
    return f">{seq_id} " + " ".join(mods)


def counts(rec) -> Counter:
    return Counter(f.type for f in rec.features if f.type in ("gene", "CDS", "tRNA", "rRNA"))


def collect(paths: List[Path]) -> Tuple[List[Tuple[Path, object]], List[str]]:
    """(file, record) of every GenBank record of the inputs, files in the order given (a folder: sorted), and the
    problems: an input without GenBank files, a file that is not readable GenBank, a record without sequence."""
    out, problems = [], []
    for p in paths:
        files = iter_gb_files(p)
        if not files:
            problems.append(f"no GenBank file in {p} (suffixes {', '.join(GB_SUFFIXES)})")
        for f in files:
            try:
                recs = list(SeqIO.parse(str(f), "genbank"))
            except Exception as exc:                  # noqa: BLE001  (a broken file must not stop the others)
                problems.append(f"{f}: not readable as GenBank ({exc}) - skipped")
                continue
            if not recs:
                problems.append(f"{f}: no GenBank record - skipped")
            for r in recs:
                try:
                    bytes(r.seq)
                except Exception:                     # noqa: BLE001  (Bio.Seq.UndefinedSequenceError: CONTIG / WGS)
                    problems.append(f"{f}: {r.id} has no sequence (a CONTIG or WGS master record?) - skipped")
                    continue
                out.append((f, r))
    return out, problems


def assign_ids(items) -> List[Tuple[Path, object, str, List[str]]]:
    """(file, record, SeqID, warnings); the same record twice (ID + sequence) is kept once, other ID clashes renamed."""
    kept, used, seqs = [], Counter(), {}
    for f, rec in items:
        rec, warn = clean_record(rec)
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


def write_pair(out: Path, name: str, items, modifiers: bool, override: Optional[int] = None) -> None:
    """<name>.tbl + <name>.fsa with the records of items (file, record, SeqID, warnings)."""
    tbl, fsa = [], []
    for _, rec, sid, _ in items:
        lines = nv.feature_lines(rec)
        lines[0] = f">Feature {sid}"
        tbl += lines
        fsa.append(defline(rec, sid, modifiers, record_code(rec, override)[0]) + "\n" + str(rec.seq).upper())
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
    p.add_argument("--genetic_code", type=int, default=None,
                   help="genetic code of every record (default: the CDS /transl_table; records without one: the usual "
                        "mitochondrial code of their phylum, reported)")
    p.add_argument("--validate", action="store_true",
                   help="run NCBI table2asn on every pair: .sqn (the file to upload) and .val, errors in the report")
    p.add_argument("--sbt", default=None, help="NCBI submitter template (.sbt) for the .sqn files (with --validate)")
    p.add_argument("--table2asn", default=None, help="table2asn executable (default: found by mito_env.py)")
    a = p.parse_args(argv)
    inputs = [Path(x) for x in a.input]
    items, problems = collect(inputs)
    for m in problems:
        print(f"[gb2tbl] {m}", file=sys.stderr)
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
        write_pair(out, name, its, a.source_modifiers == "yes", a.genetic_code)
        val = validate(out, name, t2a, a.sbt) if t2a else {}
        for f, rec, sid, warn in its:
            c = counts(rec)
            loc, _ = location(rec)
            code, code_warn = record_code(rec, a.genetic_code)
            warn = warn + ([code_warn] if code_warn else [])
            rows.append({"seq_id": sid, "record": rec.id, "source_file": str(f), "output": f"{name}.tbl",
                         "organism": rec.annotations.get("organism", ""), "length": len(rec.seq),
                         "topology": rec.annotations.get("topology", "linear"), "location": loc or "nucleus",
                         "genetic_code": code or "", "genes": c["gene"], "CDS": c["CDS"],
                         "tRNA": c["tRNA"], "rRNA": c["rRNA"],
                         "NCBI_ERROR": val.get(sid, Counter())["ERROR"] if t2a else "",
                         "NCBI_WARNING": val.get(sid, Counter())["WARNING"] if t2a else "", "warnings": "; ".join(warn)})
    rows += [{"seq_id": "", "record": "", "source_file": "", "output": "skipped", "warnings": m} for m in problems]
    write_tsv(out / "gb2tbl_report.tsv", COLUMNS, rows)
    dropped = len(items) - len(kept)
    print(f"[gb2tbl] {len(kept)} record(s) from {len(items)} in {len({str(f) for f, _ in items})} file(s)"
          + (f" ({dropped} repeated record(s) written once)" if dropped else "")
          + f" -> {len(groups)} .tbl / .fsa pair(s) in {out}", file=sys.stderr)
    cleaned = sum(any(w in r["warnings"] for w in ("left out", "misc_feature", "ASCII", "past the end", "two bases"))
                  for r in rows if r["output"] != "skipped")
    if cleaned:
        print(f"[gb2tbl] {cleaned} record(s) cleaned for NCBI (non-INSDC qualifiers / feature keys, non-ASCII text, "
              "features past the end, 100^101 sites): gb2tbl_report.tsv", file=sys.stderr)
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
