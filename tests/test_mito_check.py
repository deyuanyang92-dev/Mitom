"""Tests for mito_check.py (annotation audit report) and the run defaults.

Run: python3 -m pytest -q tests/test_mito_check.py
"""
import copy
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mitom_pipeline as bm  # noqa: E402
import mito_check as mc  # noqa: E402


def platynereis():
    return SeqIO.read(str(HERE / "data" / "NC_000931.gb"), "genbank")


def codes(issues, severity=None):
    return sorted({i.code for i in issues if severity is None or i.severity == severity})


def test_good_record_has_no_error():
    issues = mc.check_record(platynereis(), None, 22)
    assert codes(issues, "ERROR") == [], [(i.code, i.feature, i.message) for i in issues if i.severity == "ERROR"]


def test_no_annotation():
    rec = platynereis()
    rec.features = [f for f in rec.features if f.type == "source"]
    issues = mc.check_record(rec)
    assert [i.code for i in issues] == ["no_annotation"]


def test_missing_and_duplicate_genes_and_trna_count():
    rec = platynereis()
    atp8 = next(f for f in rec.features if f.type == "CDS" and mc.gene_name(f.qualifiers, "CDS") == "atp8")
    cox2 = next(f for f in rec.features if f.type == "CDS" and mc.gene_name(f.qualifiers, "CDS") == "cox2")
    rrnl = next(f for f in rec.features if f.type == "rRNA" and mc.gene_name(f.qualifiers, "rRNA") == "rrnL")
    trna = next(f for f in rec.features if f.type == "tRNA")
    rec.features.remove(atp8)
    rec.features.remove(rrnl)
    rec.features.remove(trna)
    rec.features.append(copy.deepcopy(cox2))
    issues = mc.check_record(rec)
    errs = {(i.code, i.feature) for i in issues if i.severity == "ERROR"}
    assert ("missing_pcg", "atp8") in errs and ("missing_rrna", "rrnL") in errs
    assert ("duplicate_pcg", "cox2") in errs and ("trna_count", "") in errs
    assert ("trna_count", "") not in {(i.code, i.feature) for i in mc.check_record(rec, min_trna=0)}


def test_internal_stop_frame_and_abbreviated_stop():
    rec = platynereis()
    cox1 = next(f for f in rec.features if f.type == "CDS" and mc.gene_name(f.qualifiers, "CDS") == "cox1")
    s, e = int(cox1.location.start), int(cox1.location.end)
    seq = str(rec.seq)
    # a TAA 30 codons in -> internal stop
    rec.seq = type(rec.seq)(seq[:s + 90] + "TAA" + seq[s + 93:])
    issues = mc.check_record(rec)
    assert ("internal_stop", "cox1") in {(i.code, i.feature) for i in issues if i.severity == "ERROR"}
    # cox1 is 1534 nt (abbreviated TA stop + transl_except): two nt cut -> 1532, no transl_except -> frame
    # error; with transl_except -> no error
    rec = platynereis()
    cox1 = next(f for f in rec.features if f.type == "CDS" and mc.gene_name(f.qualifiers, "CDS") == "cox1")
    cox1.location = FeatureLocation(s, e - 2, cox1.location.strand)
    cox1.qualifiers.pop("transl_except", None)
    assert ("frame", "cox1") in {(i.code, i.feature) for i in mc.check_record(rec) if i.severity == "ERROR"}
    cox1.qualifiers["transl_except"] = ["(pos:complement(1..2),aa:TERM)"]
    assert "frame" not in codes(mc.check_record(rec), "ERROR")


def test_overlap_length_and_orientation_warning():
    rec = platynereis()
    cox1 = next(f for f in rec.features if f.type == "CDS" and mc.gene_name(f.qualifiers, "CDS") == "cox1")
    cox2 = next(f for f in rec.features if f.type == "CDS" and mc.gene_name(f.qualifiers, "CDS") == "cox2")
    e = int(cox1.location.end)
    cox2.location = FeatureLocation(e - 200, e - 200 + len(cox2), 1)          # 200 nt into cox1
    issues = mc.check_record(rec)
    assert "pcg_overlap" in codes(issues, "ERROR")
    rec = platynereis()
    atp8 = next(f for f in rec.features if f.type == "CDS" and mc.gene_name(f.qualifiers, "CDS") == "atp8")
    atp8.location = FeatureLocation(int(atp8.location.start), int(atp8.location.start) + 60, 1)
    assert ("pcg_length", "atp8") in {(i.code, i.feature) for i in mc.check_record(rec)}
    rec = platynereis()
    cox1 = next(f for f in rec.features if f.type == "CDS" and mc.gene_name(f.qualifiers, "CDS") == "cox1")
    cox1.location = FeatureLocation(int(cox1.location.start), int(cox1.location.end), -1)
    w = [i for i in mc.check_record(rec) if i.code == "orientation"]
    assert w and w[0].severity == "WARNING" and "minus" in w[0].message


def test_check_files_reports_and_cli(tmp_path):
    rec_ok = platynereis()
    rec_bad = platynereis()
    rec_bad.id = "XX000001.1"
    rec_bad.name = "XX000001"
    rec_bad.annotations["accessions"] = ["XX000001"]
    rec_bad.features = [f for f in rec_bad.features if f.type == "source"]
    gb = tmp_path / "in.gb"
    SeqIO.write([rec_ok, rec_bad], str(gb), "genbank")
    out = tmp_path / "out"
    r = subprocess.run([sys.executable, str(HERE.parent / "mito_check.py"), "-i", str(gb), "-o", str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 0 and "2 record(s): 1 ok, 1 flagged" in r.stderr, r.stdout + r.stderr   # summary: stderr
    rows = [l.split("\t") for l in (out / "check_report.tsv").read_text().splitlines()]
    assert rows[0][:3] == ["record", "organism", "status"]
    status = {r[0]: r[2] for r in rows[1:]}
    assert status == {"NC_000931.1": "ok", "XX000001.1": "flagged"}
    assert (out / "check_flagged.txt").read_text().split() == ["XX000001.1"]
    details = (out / "check_details.tsv").read_text()
    assert "XX000001.1\tERROR\tno_annotation" in details


def test_run_parser_defaults():
    p = bm.build_parser()
    a = p.parse_args(["run", "-i", "x.gb", "-o", "o"])
    assert a.if_reannotation_2 == "auto" and not hasattr(a, "existing_annotation")


def test_length_and_frame_checked_without_transl_table():
    """Review 2026-10: a CDS without /transl_table (and no --genetic_code) still gets the frame and length checks;
    only the codon-table checks (internal stop, start, stop) need the code."""
    rec = platynereis()
    for f in rec.features:
        if f.type == "CDS":
            f.qualifiers.pop("transl_table", None)
            if f.qualifiers.get("gene", [""])[0].lower() == "atp8":
                s = int(f.location.start)
                f.location = FeatureLocation(s, s + 61, f.location.strand)
                f.qualifiers.pop("transl_except", None)
    issues = mc.check_record(rec, None, 22)
    assert "pcg_length" in codes(issues, "ERROR") and "frame" in codes(issues, "ERROR"), codes(issues)
    assert "no_transl_table" in codes(issues, "WARNING") and "internal_stop" not in codes(issues)
