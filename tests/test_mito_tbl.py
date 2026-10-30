"""mitom gb2tbl (mito_tbl.py, 2026-10-30): GenBank -> NCBI feature table (.tbl) + FASTA with source modifiers (.fsa);
all records in one pair, one per input file or one per record; repeated records, ID clashes, modifiers, validation.

Run: python3 -m pytest -q tests/test_mito_tbl.py   (the --validate test needs table2asn)
"""
import copy
import csv
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_env  # noqa: E402
import mito_tbl  # noqa: E402

DATA = HERE / "data"


def report(out):
    return list(csv.DictReader(open(out / "gb2tbl_report.tsv"), delimiter="\t"))


def headers(path):
    return [line for line in path.read_text().splitlines() if line.startswith(">")]


def test_all_records_in_one_pair_named_smartly(tmp_path):
    assert mito_tbl.main(["-i", str(DATA / "NC_000933.gb"), "-o", str(tmp_path / "one")]) == 0
    assert sorted(p.name for p in (tmp_path / "one").iterdir()) == ["NC_000933.fsa", "NC_000933.tbl", "gb2tbl_report.tsv"]
    assert mito_tbl.main(["-i", str(DATA / "NC_000933.gb"), str(DATA / "NC_000931.gb"), "-i",
                          str(DATA / "KY753833_mitoz.gbf"), "-o", str(tmp_path / "all")]) == 0
    tbl = (tmp_path / "all" / "submission.tbl").read_text()
    assert [x for x in tbl.splitlines() if x.startswith(">Feature")] == \
        [">Feature NC_000933.1", ">Feature NC_000931.1", ">Feature KY753833"]
    assert [h.split()[0] for h in headers(tmp_path / "all" / "submission.fsa")] == [">NC_000933.1", ">NC_000931.1", ">KY753833"]
    assert mito_tbl.main(["-i", str(DATA / "NC_000933.gb"), str(DATA / "NC_000931.gb"), "-o", str(tmp_path / "n"),
                          "--name", "batch1"]) == 0
    assert (tmp_path / "n" / "batch1.tbl").is_file()


def test_split_by_file_and_by_record_keeps_versions(tmp_path):
    src = tmp_path / "in"
    src.mkdir()
    (src / "x.gb").write_text((DATA / "NC_000933.gb").read_text())
    (src / "y.gbf").write_text((DATA / "NC_000931.gb").read_text())
    assert mito_tbl.main(["-i", str(src), "-o", str(tmp_path / "f"), "--split", "file"]) == 0
    assert {p.name for p in (tmp_path / "f").glob("*.tbl")} == {"x.tbl", "y.tbl"}
    assert mito_tbl.main(["-i", str(src), "-o", str(tmp_path / "r"), "--split", "record"]) == 0
    assert {p.name for p in (tmp_path / "r").glob("*.tbl")} == {"NC_000933.1.tbl", "NC_000931.1.tbl"}


def test_repeated_record_once_and_id_clash_renamed(tmp_path):
    other = SeqIO.read(str(DATA / "NC_000931.gb"), "genbank")
    other.id = "NC_000933.1"                                      # another genome under the same ID
    SeqIO.write(other, str(tmp_path / "clash.gbf"), "genbank")
    assert mito_tbl.main(["-i", str(DATA / "NC_000933.gb"), str(DATA / "NC_000933.gb"), str(tmp_path / "clash.gbf"),
                          "-o", str(tmp_path / "o")]) == 0
    rows = report(tmp_path / "o")
    assert [r["seq_id"] for r in rows] == ["NC_000933.1", "NC_000933.1_2"]       # the repeated record written once
    assert "also used by another record" in rows[1]["warnings"]
    assert ">Feature NC_000933.1_2" in (tmp_path / "o" / "submission.tbl").read_text()


def test_fasta_header_modifiers(tmp_path):
    rec = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    src = next(f for f in rec.features if f.type == "source")
    src.qualifiers.update(isolate=["MS-1"], specimen_voucher=["QD [001]"], country=["China: Qingdao"],
                          collection_date=["2024-05-01"])
    h = mito_tbl.defline(rec, "S1", True)
    for mod in ("[organism=Metridium senile]", "[location=mitochondrion]", "[topology=circular]",
                "[completeness=complete]", "[mgcode=4]", "[taxid=6116]", "[isolate=MS-1]",
                "[specimen-voucher=QD (001)]", "[geo_loc_name=China: Qingdao]", "[collection-date=2024-05-01]",
                "[note=white color morphe]"):
        assert mod in h, mod
    assert h.startswith(">S1 ") and "isolate" not in mito_tbl.defline(rec, "S1", False)
    nuc = copy.deepcopy(rec)
    next(f for f in nuc.features if f.type == "source").qualifiers.pop("organelle")
    assert "[location=" not in mito_tbl.defline(nuc, "S1", True) and "[gcode=4]" in mito_tbl.defline(nuc, "S1", True)
    assert mito_tbl.clean_id("my sample/1 (A)") == "my_sample_1__A_" and len(mito_tbl.clean_id("x" * 80)) == 50


def test_no_genbank_input(tmp_path):
    assert mito_tbl.main(["-i", str(tmp_path / "nothing"), "-o", str(tmp_path / "o")]) == 2


@pytest.mark.skipif(not mito_env.find("table2asn"), reason="table2asn not found")
def test_validate_makes_sqn_and_counts_errors(tmp_path):
    assert mito_tbl.main(["-i", str(DATA / "NC_000933.gb"), str(DATA / "KY753833_mitoz.gbf"), "-o", str(tmp_path / "v"),
                          "--validate"]) == 0
    assert (tmp_path / "v" / "submission.sqn").is_file()
    errs = {r["seq_id"]: int(r["NCBI_ERROR"]) for r in report(tmp_path / "v")}
    assert errs["NC_000933.1"] == 0 and errs["KY753833"] > 0                      # raw MitoZ: missing stop codons


def test_mitom_command(tmp_path):
    p = subprocess.run([sys.executable, str(HERE.parent / "mitom"), "gb2tbl", "-i", str(DATA / "NC_000931.gb"), "-o",
                        str(tmp_path / "c")], capture_output=True, text=True)
    assert p.returncode == 0 and (tmp_path / "c" / "NC_000931.tbl").is_file() and "[gb2tbl] 1 record(s)" in p.stderr
