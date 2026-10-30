#!/usr/bin/env python3
"""mito_tbl.py - GenBank / EMBL flat files -> NCBI submission files: 5-column feature table (.tbl) and FASTA with
source modifiers (.fsa), the pair BankIt and table2asn take (`mitom gb2tbl`).

Input: one or more GenBank (.gb .gbk .gbf .gbff .genbank) or EMBL / ENA (.embl .ebl .ena) files or folders (every such
file below a folder), any mix, or '-' for standard input (efetch -format gb | mitom gb2tbl -i - ...); the format is read
from the file (LOCUS / ID line). MitoM's final_annotation, MitoZ / MITOS2 / Geneious output, NCBI and ENA downloads.
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
  - the FASTA header carries [organism=] (the source /organism), [location=mitochondrion] / [location=chloroplast]
    (from /organelle; none for a nuclear record), [topology=], [completeness=complete] for a circular record,
    [mgcode=] / [gcode=] (the CDS transl_table), [taxid=] (/db_xref="taxon:..."), the source qualifiers NCBI takes as
    modifiers: isolate, strain, specimen_voucher, country / geo_loc_name, collection_date, lat_lon, collected_by,
    identified_by, isolation_source, haplotype, host, sex, tissue_type, dev_stage, note ... (--source_modifiers no: only
    the first ones), then the DEFINITION as the title (as NCBI's gbf2fsa; not a placeholder such as MitoZ's
    'topology=linear'); the sequence in lines of 70
  - the feature table: every feature but source, partial ends < >, joins, qualifiers (not /translation, /db_xref,
    /protein_id, MitoM's internal notes) - the table `mitom run` writes to ncbi_submission/
Your source information and products (the tables of GB2sequin, https://chlorobox.mpimp-golm.mpg.de/GenBank2Sequin.html):
  --source_table FILE    one row per record, tab-separated (NCBI's .src) or CSV: the first column the record (its
                         SeqID, ID with or without version, LOCUS name, or the stem of its file), the others FASTA
                         modifiers by name (isolate, specimen_voucher, country, collection_date, lat_lon, organism,
                         location, topology, genetic_code, definition, bioproject, biosample ...; case, spaces, '_'
                         and '-' do not matter); other columns are left out and reported
  --modifiers TEXT       the same for every record, as table2asn -j: "[location=mitochondrion] [topology=circular]";
                         a --source_table value wins over --modifiers, both over the record's own
  --product_table FILE   two columns, gene and product: the /product of the CDS, tRNA and rRNA with that /gene (any
                         case). Products still missing are filled from the gene name for the mitochondrial genes
                         (cox1 -> cytochrome c oxidase subunit I, trnL1 -> tRNA-Leu, rrnL -> 16S ribosomal RNA): a CDS
                         without one is an NCBI ERROR (SEQ_FEAT.MissingProteinName)
Cleaning (what other programs leave in a GenBank file and NCBI does not take; every change in the report):
  - qualifiers outside the INSDC feature table definition (Geneious /created_by /modified_by, SnapGene / ApE /label
    /ApEinfo_*, ...) are left out; a feature key outside it (Geneious 'ORF' ...) becomes misc_feature with the old key
    as a note
  - text is made ASCII (NCBI rejects anything else, SEQ_INST.NonAsciiAsn): accents dropped (Muller, Koln), mu -> u,
    typographic quotes and dashes plain; other characters removed. A file is read as UTF-8, else Windows-1252 /
    Latin-1 (older Windows programs)
  - a feature running past the end of the sequence is reported; a site between two bases (100^101) is written as the
    2-nt interval around it (as NCBI's gbf2tbl)
  - a file that is not readable GenBank / EMBL, and a record without sequence (CONTIG / WGS master records) are
    reported and skipped; the other inputs are written
  - genetic code: the CDS /transl_table; --genetic_code N for records without one (else, for a mitochondrion, the
    usual code of its phylum, reported)
--validate: table2asn (found by mito_env.py) on every pair -> <pair>.sqn (the file to upload; --sbt adds the submitter
template from https://submit.ncbi.nlm.nih.gov/genbank/template/submission/) and <pair>.val; NCBI ERROR / WARNING per
record in the report.

Usage:
  mitom gb2tbl -i final_annotation/ -o ncbi                      # all records -> ncbi/final_annotation.tbl + .fsa
  mitom gb2tbl -i a.gbf b.gb c.embl -o ncbi --name batch1         # three files -> ncbi/batch1.tbl + .fsa
  mitom gb2tbl -i genomes/ -o ncbi --split record --validate --sbt template.sbt
  mitom gb2tbl -i genomes/ -o ncbi --source_table samples.tsv --modifiers "[topology=circular]"
"""
from __future__ import annotations

import argparse
import copy
import csv
import io
import re
import sys
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402
import ncbi_validate as nv  # noqa: E402  (the feature table and the FASTA header of the NCBI submission files)
from mito_common import GB_SUFFIXES, iter_gb_files, run_tool, write_tsv  # noqa: E402
from mito_names import AA1_TO_3, AA3, PRODUCTS, gene_name, name_key, trna_name  # noqa: E402
from mito_taxa import METAZOAN_PHYLA, PHYLUM_CODE, PHYLUM_SYNONYMS  # noqa: E402

EMBL_SUFFIXES = (".embl", ".ebl", ".ena")
FASTA_WIDTH = 70                                   # sequence line length of NCBI's FASTA (gbf2fsa, efetch)
MAX_ID = 50
# GenBank source qualifier -> table2asn FASTA modifier (table2asn 1.28 checked; /country is NCBI's geo_loc_name now)
MODIFIERS = [("isolate", "isolate"), ("strain", "strain"), ("specimen_voucher", "specimen-voucher"),
             ("culture_collection", "culture-collection"), ("bio_material", "bio-material"),
             ("geo_loc_name", "geo_loc_name"), ("country", "geo_loc_name"), ("collection_date", "collection-date"),
             ("lat_lon", "lat-lon"), ("altitude", "altitude"), ("collected_by", "collected-by"),
             ("identified_by", "identified-by"), ("isolation_source", "isolation-source"), ("haplotype", "haplotype"),
             ("haplogroup", "haplogroup"), ("host", "host"), ("lab_host", "lab-host"), ("sex", "sex"),
             ("tissue_type", "tissue-type"), ("tissue_lib", "tissue-lib"), ("dev_stage", "dev-stage"),
             ("cell_line", "cell-line"), ("cell_type", "cell-type"), ("clone", "clone"), ("pop_variant", "pop-variant"),
             ("sub_species", "sub-species"), ("variety", "variety"), ("cultivar", "cultivar"), ("breed", "breed"),
             ("ecotype", "ecotype"), ("type_material", "type-material"), ("note", "note")]
# the names a --source_table column or --modifiers may use: table2asn 1.28 FASTA modifiers (each checked 2026-10-30;
# 'voucher', 'title' or 'definition' are 'Unrecognized modifier' there) + the ones MitoM puts into the record
FASTA_MODS = {m for _, m in MODIFIERS} | {
    "completeness", "taxid", "moltype", "bioproject", "biosample", "sra", "environmental-sample", "metagenomic",
    "metagenome-source", "germline", "mating-type", "fwd-pcr-primer-name", "fwd-pcr-primer-seq", "rev-pcr-primer-name",
    "rev-pcr-primer-seq"}
MOD_KEY = {**{re.sub(r"[^a-z0-9]", "", m): m for m in FASTA_MODS},
           "country": "geo_loc_name", "voucher": "specimen-voucher",
           "organism": "organism", "location": "location", "topology": "topology", "definition": "definition",
           "mgcode": "code", "gcode": "code", "geneticcode": "code", "transltable": "code"}
LOCATIONS = {"mitochondrion": "mitochondrion", "chloroplast": "chloroplast", "plastid": "plastid", "genomic": "",
             "nucleus": "", "nuclear": ""}         # [location=] value -> /organelle ('' : a nuclear record)
COLUMNS = ["seq_id", "record", "source_file", "output", "organism", "length", "topology", "location", "genetic_code",
           "genes", "CDS", "tRNA", "rRNA", "NCBI_ERROR", "NCBI_WARNING", "warnings"]
# INSDC Feature Table Definition (v11.3): feature keys and qualifiers NCBI takes (+ /codon_recognized of older NCBI
# tRNA annotation, which table2asn still takes, and NCBI's source /breed)
FEATURE_KEYS = {
    "assembly_gap", "C_region", "CDS", "centromere", "D-loop", "D_segment", "exon", "gap", "gene", "iDNA", "intron",
    "J_segment", "mat_peptide", "misc_binding", "misc_difference", "misc_feature", "misc_recomb", "misc_RNA",
    "misc_structure", "mobile_element", "modified_base", "mRNA", "ncRNA", "N_region", "old_sequence", "operon", "oriT",
    "polyA_site", "precursor_RNA", "prim_transcript", "primer_bind", "propeptide", "protein_bind", "regulatory",
    "repeat_region", "rep_origin", "rRNA", "S_region", "sig_peptide", "source", "stem_loop", "STS", "telomere", "tmRNA",
    "transit_peptide", "tRNA", "unsure", "V_region", "V_segment", "variation", "3'UTR", "5'UTR"}
QUALIFIERS = {
    "allele", "altitude", "anticodon", "artificial_location", "bio_material", "bound_moiety", "breed", "cell_line",
    "cell_type", "chromosome", "circular_RNA", "citation", "clone", "clone_lib", "codon_recognized", "codon_start",
    "collected_by", "collection_date",
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


class Entry(NamedTuple):
    """A record as written: its input file, the cleaned record, SeqID, warnings, and what --genetic_code /
    --modifiers / --source_table gave (genetic code, other FASTA modifiers)."""
    file: Path
    rec: object
    sid: str
    warn: List[str]
    code: Optional[int]
    extra: Dict[str, str]


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
    if not str(rec.description or "").isascii():
        rec.description = ascii_text(rec.description)
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


def organism(rec) -> str:
    """The source /organism (an EMBL OS line also has the common name), else the ORGANISM line; '' for none."""
    org = " ".join(source_of(rec).get("organism", [])).strip() or str(rec.annotations.get("organism") or "").strip()
    return "" if org.lower().rstrip(".") in ("", "unknown") else org


def location(rec) -> Tuple[str, str]:
    """(FASTA location modifier or '', genetic code modifier: mgcode for an organelle, gcode for a nuclear record)."""
    org = " ".join(source_of(rec).get("organelle", [])).lower()
    if "mitochondr" in org:
        return "mitochondrion", "mgcode"
    if "chloroplast" in org or "plastid" in org:
        return "chloroplast", "gcode"
    return "", "gcode"


def title(rec) -> str:
    """The FASTA title (NCBI's definition line): the record's DEFINITION / DE line; '' for none or a placeholder
    ('.', MitoZ's 'topology=linear', the record's ID)."""
    d = " ".join(str(rec.description or "").split()).rstrip(".").strip()
    if not d or "=" in d or d in ("<unknown description>", rec.id, rec.name):
        return ""
    return ascii_text(d.replace("[", "(").replace("]", ")"))


def defline(rec, seq_id: str, modifiers: bool, code: Optional[int] = None,
            extra: Optional[Dict[str, str]] = None) -> str:
    """FASTA header: [organism=] [location=] [topology=] [completeness=] [mgcode= / gcode=] + source modifiers + extra
    (--modifiers / --source_table, in place of the record's own), then the title; code: the genetic code
    (record_code), default the CDS /transl_table."""
    extra = extra or {}
    loc, code_mod = location(rec)
    src = source_of(rec)
    org = organism(rec) or "unknown organism"
    topo = rec.annotations.get("topology", "linear")
    mods = [f"[organism={org}]"] + ([f"[location={loc}]"] if loc else []) + [f"[topology={topo}]"]
    if topo == "circular" and "completeness" not in extra:    # else SEQ_INST.CompleteCircleProblem (as fasta_defline)
        mods.append("[completeness=complete]")
    code = code or nv.record_gencode(rec)
    if code:
        mods.append(f"[{code_mod}={code}]")
    taxid = next((x.split(":", 1)[1] for x in src.get("db_xref", []) if x.startswith("taxon:")), "")
    if taxid.isdigit() and "taxid" not in extra:  # (else table2asn warns SEQ_DESCR.NoTaxonID)
        mods.append(f"[taxid={taxid}]")
    if modifiers:
        seen = set()
        for qual, mod in MODIFIERS:
            if mod in extra:
                continue
            for v in src.get(qual, []):
                v = ascii_text(str(v).replace("[", "(").replace("]", ")"))
                if v and (mod, v) not in seen:
                    seen.add((mod, v))
                    mods.append(f"[{mod}={v}]")
    mods += [f"[{m}={v}]" for m, v in extra.items()]
    t = title(rec)
    return f">{seq_id} " + " ".join(mods) + (f" {t}" if t else "")


def counts(rec) -> Counter:
    return Counter(f.type for f in rec.features if f.type in ("gene", "CDS", "tRNA", "rRNA"))


def read_text(path: Path) -> str:
    """A file's text ('-': standard input): UTF-8, else Windows-1252 / Latin-1 (older Windows programs)."""
    if str(path) == "-":
        return sys.stdin.read()
    raw = Path(path).read_bytes()
    for enc in ("utf-8-sig", "cp1252"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            pass
    return raw.decode("latin-1")


def collect(paths: List[Path]) -> Tuple[List[Tuple[Path, object]], List[str]]:
    """(file, record) of every GenBank / EMBL record of the inputs, files in the order given (a folder: sorted; '-':
    standard input, file 'stdin'), and the problems: an input without such files, a file that is not readable, a record
    without sequence."""
    out, problems = [], []
    suffixes = GB_SUFFIXES + EMBL_SUFFIXES
    for p in paths:
        files = [p] if str(p) == "-" else iter_gb_files(p, suffixes)
        if not files:
            problems.append(f"no GenBank / EMBL file in {p} (suffixes {', '.join(suffixes)})")
        for f in files:
            name = Path("stdin") if str(f) == "-" else f
            try:
                text = read_text(f)
                recs = list(SeqIO.parse(io.StringIO(text), "embl" if re.match(r"\s*ID   ", text) else "genbank"))
            except Exception as exc:                  # noqa: BLE001  (a broken file must not stop the others)
                problems.append(f"{name}: not readable as GenBank / EMBL ({exc}) - skipped")
                continue
            if not recs:
                problems.append(f"{name}: no GenBank / EMBL record - skipped")
            for r in recs:
                try:
                    empty = not bytes(r.seq)
                except Exception:                     # noqa: BLE001  (Bio.Seq.UndefinedSequenceError: CONTIG / WGS)
                    empty = True
                if empty:
                    problems.append(f"{name}: {r.id} has no sequence (a CONTIG or WGS master record?) - skipped")
                    continue
                out.append((name, r))
    return out, problems


def assign_ids(items) -> List[Tuple[Path, object, str, List[str]]]:
    """(file, cleaned record, SeqID, warnings); the same record twice (ID + sequence) is kept once, other ID clashes
    renamed."""
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
        kept.append((f, rec, sid, warn))
    return kept


def mod_name(name: str) -> str:
    """The FASTA modifier a --source_table column / --modifiers name stands for ('Specimen voucher' ->
    specimen-voucher; 'code' for the genetic code); '' when none."""
    return MOD_KEY.get(re.sub(r"[^a-z0-9]", "", str(name).lower()), "")


def parse_mods(text: Optional[str]) -> Tuple[Dict[str, str], List[str]]:
    """--modifiers "[organism=Xus yus] [location=mitochondrion]" -> ({modifier: value}, problems)."""
    mods, notes = {}, []
    pairs = re.findall(r"\[\s*([^=\]]+?)\s*=\s*([^\]]*)\]", text or "")
    if text and not pairs:
        notes.append(f"--modifiers {text!r}: write them as [name=value]")
    for k, v in pairs:
        if mod_name(k):
            mods[mod_name(k)] = v.strip()
        else:
            notes.append(f"--modifiers: {k} is not an NCBI source modifier - left out")
    return mods, notes


def read_rows(path: Path) -> List[List[str]]:
    """Rows of a tab-separated (NCBI .src) or comma-separated table; blank lines left out."""
    lines = [x for x in read_text(path).splitlines() if x.strip()]
    delim = "\t" if lines and "\t" in lines[0] else ","
    return [[c.strip() for c in row] for row in csv.reader(lines, delimiter=delim)]


def source_table(path: Path) -> Tuple[Dict[str, Tuple[str, Dict[str, str]]], List[str]]:
    """--source_table: {record key (lower case): (key as written, {modifier: value})}, and what was left out."""
    rows = read_rows(path)
    head = rows[0][1:] if rows else []
    names = [mod_name(h) for h in head]
    other = [h for h, n in zip(head, names) if h and not n]
    notes = [f"{path.name}: columns that are not NCBI source modifiers left out: {', '.join(other)}"] if other else []
    table: Dict[str, Tuple[str, Dict[str, str]]] = {}
    for r in rows[1:]:
        if not r or not r[0]:
            continue
        if r[0].lower() in table:
            notes.append(f"{path.name}: {r[0]} in two rows - the first used")
            continue
        table[r[0].lower()] = (r[0], {n: v for n, v in zip(names, r[1:]) if n and v})
    return table, notes


def table_key(table, f: Path, rec, sid: str) -> Optional[str]:
    """The --source_table row of a record: by SeqID, ID, ID without version, LOCUS name or file stem."""
    for k in (sid, rec.id, rec.id.split(".")[0], rec.name, f.stem):
        if k and k.lower() in table:
            return k.lower()
    return None


def product_table(path: Path) -> Dict[str, str]:
    """--product_table: {gene name key: product}; a first row 'gene ...' is the header."""
    rows = read_rows(path)
    if rows and rows[0][0].lower() == "gene":
        rows = rows[1:]
    return {name_key(r[0]): r[1] for r in rows if len(r) > 1 and r[0] and r[1]}


def apply_source(rec, mods: Dict[str, str]) -> Tuple[Optional[int], Dict[str, str], List[str]]:
    """--modifiers / --source_table values on rec: organism, location, topology and definition go into the record (the
    report and the genetic code read them there); returns the genetic code given, the other modifiers, problems."""
    code, extra, warn = None, {}, []
    if not mods:
        return code, extra, warn
    src = next((f for f in rec.features if f.type == "source"), None)
    if src is None:
        src = SeqFeature(FeatureLocation(0, len(rec.seq), 1), type="source")
        rec.features.insert(0, src)
    for mod, v in mods.items():
        v = ascii_text(v.replace("[", "(").replace("]", ")"))
        if not v:
            continue
        if mod == "organism":
            old = organism(rec)
            if old and old != v:                       # its taxid and its name in the definition go with it
                if v not in str(rec.description or ""):
                    rec.description = str(rec.description or "").replace(old, v)
                src.qualifiers["db_xref"] = [x for x in src.qualifiers.get("db_xref", []) if not x.startswith("taxon:")]
            src.qualifiers["organism"] = [v]
            rec.annotations["organism"] = v
        elif mod == "location":
            if v.lower() not in LOCATIONS:
                warn.append(f"location {v!r} left out: it is the place in the cell (mitochondrion, chloroplast, "
                            "genomic); the place of collection is geo_loc_name (country)")
            elif LOCATIONS[v.lower()]:
                src.qualifiers["organelle"] = [LOCATIONS[v.lower()]]
            else:
                src.qualifiers.pop("organelle", None)
        elif mod == "topology":
            if v.lower() in ("circular", "linear"):
                rec.annotations["topology"] = v.lower()
            else:
                warn.append(f"topology {v!r} left out (circular or linear)")
        elif mod == "code":
            if v.isdigit():
                code = int(v)
            else:
                warn.append(f"genetic code {v!r} left out (a number, e.g. 5)")
        elif mod == "definition":
            rec.description = v
        else:
            extra[mod] = v
    return code, extra, warn


def fill_products(rec, table: Dict[str, str]) -> List[str]:
    """/product of the CDS, tRNA and rRNA: from --product_table by /gene; where still missing, the name of a
    mitochondrial gene (a CDS without one is NCBI's ERROR SEQ_FEAT.MissingProteinName). What was changed."""
    edited, filled, missing = 0, 0, Counter()
    for f in rec.features:
        if f.type not in ("CDS", "tRNA", "rRNA"):
            continue
        q = f.qualifiers
        new = table.get(name_key(q.get("gene", [""])[0])) if table else None
        if new:
            edited += q.get("product") != [ascii_text(new)]
            q["product"] = [ascii_text(new)]
            continue
        if q.get("product"):
            continue
        if f.type == "tRNA":
            aa = trna_name(q)[0]
            m = re.search(r"aa:([A-Z][a-z]{2})", " ".join(q.get("anticodon", [])))   # (pos:...,aa:Leu,seq:taa)
            prod = f"tRNA-{AA1_TO_3[aa]}" if aa else f"tRNA-{m.group(1)}" if m and m.group(1) in AA3 else ""
        else:
            prod = PRODUCTS.get(gene_name(q, f.type, atp9=True), "")
        if prod:
            q["product"] = [prod]
            filled += 1
        else:
            missing[f.type] += 1
    warn = [f"{edited} product(s) set from the product table"] if edited else []
    if filled:
        warn.append(f"{filled} missing product(s) filled from the gene name")
    if missing:
        warn.append("no product and the gene not recognised: " + ", ".join(f"{k} x{v}" for k, v in missing.items())
                    + (" (a CDS without one is an NCBI ERROR: --product_table)" if missing["CDS"] else ""))
    return warn


def finish(f: Path, rec, sid: str, warn: List[str], mods: Dict[str, str], products: Dict[str, str],
           gcode: Optional[int]) -> Entry:
    """A record with the source values and products given, and the checks of what NCBI needs."""
    code, extra, w = apply_source(rec, mods)
    warn = warn + w + fill_products(rec, products)
    if not organism(rec):
        warn.append("no organism: [organism=unknown organism] - set it (--source_table / --modifiers) before submitting")
    if not any(x.type != "source" for x in rec.features):
        warn.append("no features: the table has the header line only")
    d = " ".join(str(rec.description or "").split())
    if d and not title(rec) and d not in (".", "<unknown description>"):
        warn.append(f"DEFINITION {d!r} left out of the FASTA title (not a description)")
    return Entry(f, rec, sid, warn, code or gcode, extra)


def pair(out: Path, name: str, ext: str) -> Path:
    """OUT/<name><ext> (not Path.with_suffix: a name such as NC_000931.1 keeps its version)."""
    return out / f"{name}{ext}"


def write_pair(out: Path, name: str, entries: List[Entry], modifiers: bool) -> None:
    """<name>.tbl + <name>.fsa with the records of entries."""
    tbl, fsa = [], []
    for e in entries:
        lines = nv.feature_lines(e.rec)
        lines[0] = f">Feature {e.sid}"
        tbl += lines
        seq = str(e.rec.seq).upper()
        fsa.append(defline(e.rec, e.sid, modifiers, record_code(e.rec, e.code)[0], e.extra) + "\n"
                   + "\n".join(seq[i:i + FASTA_WIDTH] for i in range(0, len(seq), FASTA_WIDTH)))
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
                   help="GenBank / EMBL files and / or folders, '-' for standard input (repeatable: -i a.gbf b.embl "
                        "-i more/)")
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
    p.add_argument("--source_table", default=None,
                   help="source information per record: tab-separated (NCBI .src) or CSV; first column the record "
                        "(SeqID, ID, LOCUS name or file stem), the others FASTA modifiers by name (isolate, "
                        "specimen_voucher, country, collection_date, lat_lon, organism, topology, definition ...)")
    p.add_argument("--modifiers", default=None,
                   help="FASTA modifiers for every record, as table2asn -j: \"[location=mitochondrion] "
                        "[topology=circular]\" (a --source_table value wins)")
    p.add_argument("--product_table", default=None,
                   help="gene and product (two columns): the /product of the CDS / tRNA / rRNA with that /gene "
                        "(missing products of mitochondrial genes are filled from the gene name anyway)")
    p.add_argument("--genetic_code", type=int, default=None,
                   help="genetic code of every record (default: the CDS /transl_table; records without one: the usual "
                        "mitochondrial code of their phylum, reported)")
    p.add_argument("--validate", action="store_true",
                   help="run NCBI table2asn on every pair: .sqn (the file to upload) and .val, errors in the report")
    p.add_argument("--sbt", default=None, help="NCBI submitter template (.sbt) for the .sqn files (with --validate)")
    p.add_argument("--table2asn", default=None, help="table2asn executable (default: found by mito_env.py)")
    a = p.parse_args(argv)
    for opt, path in (("--source_table", a.source_table), ("--product_table", a.product_table)):
        if path and not Path(path).is_file():
            print(f"[gb2tbl] {opt} {path}: no such file", file=sys.stderr)
            return 2
    inputs = [Path(x) for x in a.input]
    items, problems = collect(inputs)
    for m in problems:
        print(f"[gb2tbl] {m}", file=sys.stderr)
    if not items:
        return 2
    kept = assign_ids(items)
    glob, notes = parse_mods(a.modifiers)
    table, table_notes = source_table(Path(a.source_table)) if a.source_table else ({}, [])
    products = product_table(Path(a.product_table)) if a.product_table else {}
    entries, used = [], set()
    for f, rec, sid, warn in kept:
        key = table_key(table, f, rec, sid) if table else None
        if table and key is None:
            warn = warn + ["no row in the source table"]
        used.add(key)
        entries.append(finish(f, rec, sid, warn, {**glob, **(table[key][1] if key else {})}, products,
                              a.genetic_code))
    notes += table_notes + [f"{Path(a.source_table).name}: row {table[k][0]} matches no record"
                            for k in table if k not in used]
    for m in notes:
        print(f"[gb2tbl] {m}", file=sys.stderr)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    groups: Dict[str, List[Entry]] = defaultdict(list)
    for e in entries:
        if a.split == "record":
            groups[e.sid].append(e)
        elif a.split == "file":
            groups[clean_id(e.file.stem)].append(e)
        else:
            single = len(inputs) == 1 and inputs[0].is_file()
            groups[a.name or (clean_id(inputs[0].stem) if single else "submission")].append(e)
    t2a = mito_env.find("table2asn", a.table2asn) if a.validate else None
    if a.validate and not t2a:
        print("[gb2tbl] table2asn not found (mitom doctor; conda install -c bioconda table2asn): no .sqn / check",
              file=sys.stderr)
    rows = []
    for name, its in groups.items():
        write_pair(out, name, its, a.source_modifiers == "yes")
        val = validate(out, name, t2a, a.sbt) if t2a else {}
        for e in its:
            c = counts(e.rec)
            loc, _ = location(e.rec)
            code, code_warn = record_code(e.rec, e.code)
            warn = e.warn + ([code_warn] if code_warn else [])
            rows.append({"seq_id": e.sid, "record": e.rec.id, "source_file": str(e.file), "output": f"{name}.tbl",
                         "organism": organism(e.rec), "length": len(e.rec.seq),
                         "topology": e.rec.annotations.get("topology", "linear"), "location": loc or "nucleus",
                         "genetic_code": code or "", "genes": c["gene"], "CDS": c["CDS"],
                         "tRNA": c["tRNA"], "rRNA": c["rRNA"],
                         "NCBI_ERROR": val.get(e.sid, Counter())["ERROR"] if t2a else "",
                         "NCBI_WARNING": val.get(e.sid, Counter())["WARNING"] if t2a else "",
                         "warnings": "; ".join(warn)})
    rows += [{"seq_id": "", "record": "", "source_file": "", "output": "skipped", "warnings": m} for m in problems]
    rows += [{"seq_id": "", "record": "", "source_file": "", "output": "note", "warnings": m} for m in notes]
    write_tsv(out / "gb2tbl_report.tsv", COLUMNS, rows)
    dropped = len(items) - len(kept)
    print(f"[gb2tbl] {len(kept)} record(s) from {len(items)} in {len({str(f) for f, _ in items})} file(s)"
          + (f" ({dropped} repeated record(s) written once)" if dropped else "")
          + f" -> {len(groups)} .tbl / .fsa pair(s) in {out}", file=sys.stderr)
    if table:
        print(f"[gb2tbl] --source_table: {sum(k is not None for k in used)} row(s) used, "
              f"{sum('no row in the source table' in e.warn for e in entries)} record(s) without a row",
              file=sys.stderr)
    cleaned = sum(any(w in r["warnings"] for w in ("left out", "misc_feature", "ASCII", "past the end", "two bases",
                                                   "product"))
                  for r in rows if r["output"] not in ("skipped", "note"))
    if cleaned:
        print(f"[gb2tbl] {cleaned} record(s) cleaned for NCBI (non-INSDC qualifiers / feature keys, non-ASCII text, "
              "features past the end, 100^101 sites, missing products): gb2tbl_report.tsv", file=sys.stderr)
    renamed = sum("also used" in r["warnings"] for r in rows)
    if renamed:
        print(f"[gb2tbl] {renamed} record(s) renamed: their ID was used by another record (gb2tbl_report.tsv)",
              file=sys.stderr)
    if t2a:
        e = sum(int(r["NCBI_ERROR"] or 0) for r in rows if r["output"] not in ("skipped", "note"))
        print(f"[gb2tbl] table2asn: {e} ERROR(s) in {sum(1 for r in rows if r.get('NCBI_ERROR'))} record(s); .sqn "
              "files to upload" + ("" if a.sbt else " (add --sbt template.sbt for the submitter)"), file=sys.stderr)
        if e:
            print("[gb2tbl] the messages: <pair>.val; most annotation errors (missing stop codons, overlaps) are "
                  "fixed by `mitom curate -i <GenBank> -o curated` (or `mitom import` for MitoZ / MITOS2 results), "
                  "then gb2tbl again", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
