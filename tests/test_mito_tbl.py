"""mitom gb2tbl (mito_tbl.py, 2026-10-30): GenBank / EMBL -> NCBI feature table (.tbl) + FASTA with source modifiers
(.fsa); all records in one pair, one per input file or one per record; repeated records, ID clashes, modifiers, the
source table, --modifiers, products, titles, cleaning, validation.

Run: python3 -m pytest -q tests/test_mito_tbl.py   (the --validate test needs table2asn)
"""
import copy
import csv
import io
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_env  # noqa: E402
import mito_names  # noqa: E402
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


def test_cleaning_for_ncbi(tmp_path):
    """What other programs leave in GenBank files and NCBI does not take: editor qualifiers, a non-INSDC feature key,
    non-ASCII text, a feature past the end - cleaned and reported (stress tests 2026-10-30)."""
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    rec = SeqIO.read(str(DATA / "NC_000931.gb"), "genbank")
    cds = next(f for f in rec.features if f.type == "CDS")
    cds.qualifiers.update(created_by=["User"], label=["COX1"], ApEinfo_fwdcolor=["#85dae9"],
                          note=["collected by M\u00fcller \u6d4b\u8bd5 5 \u00b5m, \u201cK\u00f6ln\u201d"])
    rec.features.append(SeqFeature(FeatureLocation(8000, 8300, 1), type="ORF", qualifiers={}))
    rec.features.append(SeqFeature(FeatureLocation(len(rec) - 50, len(rec) + 40, 1), type="misc_feature"))
    clean, warn = mito_tbl.clean_record(rec)
    c = next(f for f in clean.features if f.type == "CDS")
    assert not {"created_by", "label", "ApEinfo_fwdcolor"} & set(c.qualifiers)
    assert c.qualifiers["note"] == ['collected by Muller 5 um, "Koln"']
    assert any(f.type == "misc_feature" and f.qualifiers.get("note") == ["ORF"] for f in clean.features)
    text = "; ".join(warn)
    assert "ORF x1" in text and "created_by x1" in text and "made ASCII" in text and "past the end" in text
    rec.features.append(SeqFeature(FeatureLocation(100, 100, 1), type="misc_feature", qualifiers={"note": ["site"]}))
    site = next(f for f in mito_tbl.clean_record(rec)[0].features if f.qualifiers.get("note") == ["site"])
    assert (int(site.location.start), int(site.location.end)) == (99, 101)             # 100^101 -> 100..101
    trna = next(f for f in rec.features if f.type == "tRNA")
    trna.qualifiers["codon_recognized"] = ["UUA"]                  # older NCBI tRNA qualifier, table2asn takes it
    assert "codon_recognized" in next(f for f in mito_tbl.clean_record(rec)[0].features if f.type == "tRNA").qualifiers


def test_broken_inputs_are_reported_and_skipped(tmp_path):
    (tmp_path / "half.gb").write_text((DATA / "NC_000931.gb").read_text()[:6000])
    (tmp_path / "contig.gb").write_text(
        "LOCUS       XX000001                 500 bp    DNA     linear   INV 01-JAN-2020\n"
        "DEFINITION  contig-only record.\nACCESSION   XX000001\nVERSION     XX000001.1\n"
        "FEATURES             Location/Qualifiers\n     source          1..500\n"
        "                     /organism=\"Testus testus\"\nCONTIG      join(XX000002.1:1..500)\n//\n")
    (tmp_path / "good.gb").write_text((DATA / "NC_000933.gb").read_text())
    assert mito_tbl.main(["-i", str(tmp_path), "-o", str(tmp_path / "o")]) == 0          # the good file is written
    rows = report(tmp_path / "o")
    assert [r["seq_id"] for r in rows if r["output"] != "skipped"] == ["NC_000933.1"]
    skipped = " ".join(r["warnings"] for r in rows if r["output"] == "skipped")
    assert "not readable as GenBank" in skipped and "has no sequence" in skipped
    assert mito_tbl.main(["-i", str(tmp_path / "half.gb"), "-o", str(tmp_path / "x")]) == 2


def test_genetic_code_from_option_or_taxonomy():
    rec = SeqIO.read(str(DATA / "NC_000931.gb"), "genbank")
    for f in rec.features:
        f.qualifiers.pop("transl_table", None)
    code, warn = mito_tbl.record_code(rec, None)
    assert code == 5 and "assumed (Annelida)" in warn and "[mgcode=5]" in mito_tbl.defline(rec, "X", True, code)
    assert mito_tbl.record_code(rec, 9) == (9, "")
    vert = SeqIO.read(str(DATA / "NC_000931.gb"), "genbank")
    for f in vert.features:
        f.qualifiers.pop("transl_table", None)
    vert.annotations["taxonomy"] = ["Eukaryota", "Metazoa", "Chordata", "Vertebrata"]
    assert mito_tbl.record_code(vert, None)[0] == 2
    vert.annotations["taxonomy"] = ["Eukaryota", "Metazoa", "Chordata", "Tunicata", "Ascidiacea"]
    assert mito_tbl.record_code(vert, None)[0] == 13


def test_embl_stdin_and_latin1_inputs(tmp_path, monkeypatch):
    """GB2sequin reads EMBL / ENA files, NCBI's gbf2tbl reads a pipe; older Windows programs write Latin-1."""
    src = tmp_path / "in"
    src.mkdir()
    SeqIO.write(SeqIO.read(str(DATA / "NC_000931.gb"), "genbank"), str(src / "pd.embl"), "embl")
    assert mito_tbl.main(["-i", str(src / "pd.embl"), "-o", str(tmp_path / "e")]) == 0
    assert mito_tbl.main(["-i", str(DATA / "NC_000931.gb"), "-o", str(tmp_path / "g")]) == 0
    assert (tmp_path / "e" / "pd.tbl").read_text() == (tmp_path / "g" / "NC_000931.tbl").read_text()
    assert mito_tbl.main(["-i", str(src), "-o", str(tmp_path / "d")]) == 0               # a folder: .embl found
    assert (tmp_path / "d" / "submission.tbl").is_file()
    monkeypatch.setattr(sys, "stdin", io.StringIO((DATA / "NC_000933.gb").read_text()))
    assert mito_tbl.main(["-i", "-", "-o", str(tmp_path / "s"), "--split", "file"]) == 0
    assert ">Feature NC_000933.1" in (tmp_path / "s" / "stdin.tbl").read_text()
    text = (DATA / "NC_000931.gb").read_text().replace(
        '/organism="Platynereis dumerilii"', '/organism="Platynereis dumerilii"\n                     '
        '/collected_by="J. M\u00fcller"', 1)
    (tmp_path / "l1.gb").write_bytes(text.encode("latin-1"))
    assert mito_tbl.main(["-i", str(tmp_path / "l1.gb"), "-o", str(tmp_path / "l")]) == 0
    fsa = (tmp_path / "l" / "l1.fsa").read_text().splitlines()
    assert "[collected-by=J. Muller]" in fsa[0] and max(len(x) for x in fsa[1:]) == 70       # lines of 70, as NCBI
    (tmp_path / "w.gb").write_bytes(text.replace("J. M\u00fcller", "J. \u201cM\u00fcller\u201d").encode("cp1252"))
    assert mito_tbl.main(["-i", str(tmp_path / "w.gb"), "-o", str(tmp_path / "w")]) == 0      # Windows-1252 quotes
    assert '[collected-by=J. "Muller"]' in headers(tmp_path / "w" / "w.fsa")[0]


def test_title_from_definition():
    rec = SeqIO.read(str(DATA / "NC_000931.gb"), "genbank")
    assert mito_tbl.defline(rec, "X", True).endswith("[taxid=6359] Platynereis dumerilii mitochondrion, complete genome")
    for junk in ("topology=linear.", ".", "NC_000931.1", ""):                    # MitoZ writes 'topology=linear'
        rec.description = junk
        assert mito_tbl.title(rec) == ""
    rec.description = "Xus yus [sample 1] mitochondrion, complete genome."
    assert mito_tbl.title(rec) == "Xus yus (sample 1) mitochondrion, complete genome"


def test_source_table_and_modifiers(tmp_path):
    """GB2sequin's source modifier table / table2asn's .src and -j: per record and for all, the table winning."""
    (tmp_path / "samples.tsv").write_text(
        "Sequence_ID\tIsolate\tCountry\tVoucher\tLocation\tFamily\tcollection_date\tDefinition\n"
        "NC_000933\tMS-1\tChina: Qingdao\tQD [1]\tQingdao\tMetridiidae\t2024-05-01\t\n"
        "NC_000931.1\tPD-7\t\t\tmitochondrion\t\t\tPlatynereis dumerilii isolate PD-7 mitochondrion, complete genome\n"
        "nothere\tX\t\t\t\t\t\t\n")
    assert mito_tbl.main(["-i", str(DATA / "NC_000933.gb"), str(DATA / "NC_000931.gb"), "-o", str(tmp_path / "o"),
                          "--source_table", str(tmp_path / "samples.tsv"),
                          "--modifiers", "[isolate=ALL] [sex=female] [colour=red]"]) == 0
    h = {x.split()[0]: x for x in headers(tmp_path / "o" / "submission.fsa")}
    a, b = h[">NC_000933.1"], h[">NC_000931.1"]
    assert "[isolate=MS-1]" in a and "[isolate=ALL]" not in a and "[sex=female]" in a and "[sex=female]" in b
    for mod in ("[geo_loc_name=China: Qingdao]", "[specimen-voucher=QD (1)]", "[collection-date=2024-05-01]",
                "[location=mitochondrion]"):
        assert mod in a, mod
    assert "Metridiidae" not in a and "Qingdao]" not in a.split("geo_loc_name")[0]
    assert "[isolate=PD-7]" in b and b.endswith("] Platynereis dumerilii isolate PD-7 mitochondrion, complete genome")
    text = " ".join(r["warnings"] for r in report(tmp_path / "o"))
    for msg in ("location 'Qingdao' left out", "not NCBI source modifiers left out: Family", "row nothere matches no "
                "record", "colour is not an NCBI source modifier"):
        assert msg in text, msg
    (tmp_path / "org.csv").write_text("id,organism,topology,genetic_code\nNC_000931,Platynereis massiliensis,linear,9\n")
    assert mito_tbl.main(["-i", str(DATA / "NC_000931.gb"), "-o", str(tmp_path / "m"), "--source_table",
                          str(tmp_path / "org.csv")]) == 0
    h = headers(tmp_path / "m" / "NC_000931.fsa")[0]
    assert "[organism=Platynereis massiliensis]" in h and "[taxid=" not in h and "[mgcode=9]" in h
    assert "[topology=linear]" in h and "[completeness=" not in h
    assert h.endswith("] Platynereis massiliensis mitochondrion, complete genome")
    assert report(tmp_path / "m")[0]["organism"] == "Platynereis massiliensis"


def no_products(path):
    """NC_000931 as an editor export with gene names and no /product: tRNA-Leu -> /gene="trnL", small / large subunit
    rRNA -> rrnS / rrnL."""
    rec = SeqIO.read(str(DATA / "NC_000931.gb"), "genbank")
    for f in rec.features:
        if f.type in ("CDS", "tRNA", "rRNA"):
            prod = f.qualifiers.pop("product")[0]
            if f.type == "tRNA":
                f.qualifiers["gene"] = ["trn" + mito_names.AA3[prod[5:8]]]
            elif f.type == "rRNA":
                f.qualifiers["gene"] = ["rrnS" if prod.startswith("small") else "rrnL"]
    SeqIO.write(rec, str(path), "genbank")


def test_products_filled_or_from_table(tmp_path):
    no_products(tmp_path / "noprod.gb")
    (tmp_path / "prod.tsv").write_text("gene\tproduct\ncox1\tcytochrome c oxidase subunit 1\n")
    assert mito_tbl.main(["-i", str(tmp_path / "noprod.gb"), "-o", str(tmp_path / "o"), "--product_table",
                          str(tmp_path / "prod.tsv")]) == 0
    tbl = (tmp_path / "o" / "noprod.tbl").read_text()
    for prod in ("cytochrome c oxidase subunit 1", "cytochrome c oxidase subunit II", "NADH dehydrogenase subunit 4L",
                 "cytochrome b", "tRNA-Asn", "tRNA-Leu", "12S ribosomal RNA", "16S ribosomal RNA"):
        assert f"\t\t\tproduct\t{prod}\n" in tbl, prod
    w = report(tmp_path / "o")[0]["warnings"]
    assert "1 product(s) set from the product table" in w and "36 missing product(s) filled" in w   # 12 + 22 + 2
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(Seq("A" * 300), id="T")
    rec.features = [SeqFeature(FeatureLocation(0, 60, 1), type="tRNA",
                               qualifiers={"anticodon": ["(pos:30..32,aa:Leu,seq:taa)"]}),
                    SeqFeature(FeatureLocation(70, 130, 1), type="tRNA"),
                    SeqFeature(FeatureLocation(140, 290, 1), type="CDS", qualifiers={"gene": ["orf150"]})]
    w = "; ".join(mito_tbl.fill_products(rec, {}))
    assert rec.features[0].qualifiers["product"] == ["tRNA-Leu"]                 # from the /anticodon
    assert "no product and the gene not recognised: tRNA x1, CDS x1 (a CDS without one is an NCBI ERROR" in w


def test_no_genbank_input(tmp_path):
    assert mito_tbl.main(["-i", str(tmp_path / "nothing"), "-o", str(tmp_path / "o")]) == 2


@pytest.mark.skipif(not mito_env.find("table2asn"), reason="table2asn not found")
def test_validate_makes_sqn_and_counts_errors(tmp_path):
    assert mito_tbl.main(["-i", str(DATA / "NC_000933.gb"), str(DATA / "KY753833_mitoz.gbf"), "-o", str(tmp_path / "v"),
                          "--validate"]) == 0
    assert (tmp_path / "v" / "submission.sqn").is_file()
    errs = {r["seq_id"]: int(r["NCBI_ERROR"]) for r in report(tmp_path / "v")}
    assert errs["NC_000933.1"] == 0 and errs["KY753833"] > 0                      # raw MitoZ: missing stop codons


@pytest.mark.skipif(not mito_env.find("table2asn"), reason="table2asn not found")
def test_validate_with_filled_products_and_a_broken_file(tmp_path):
    no_products(tmp_path / "noprod.gb")                     # CDS without /product: SEQ_FEAT.MissingProteinName
    (tmp_path / "half.gb").write_text((DATA / "NC_000931.gb").read_text()[:6000])
    assert mito_tbl.main(["-i", str(tmp_path / "noprod.gb"), str(tmp_path / "half.gb"), "-o", str(tmp_path / "v"),
                          "--validate", "--modifiers", "[isolation-source=marine sediment]"]) == 0
    rows = [r for r in report(tmp_path / "v") if r["output"] not in ("skipped", "note")]
    assert [(r["seq_id"], r["NCBI_ERROR"]) for r in rows] == [("NC_000931.1", "0")]


def test_mitom_command(tmp_path):
    p = subprocess.run([sys.executable, str(HERE.parent / "mitom"), "gb2tbl", "-i", str(DATA / "NC_000931.gb"), "-o",
                        str(tmp_path / "c")], capture_output=True, text=True)
    assert p.returncode == 0 and (tmp_path / "c" / "NC_000931.tbl").is_file() and "[gb2tbl] 1 record(s)" in p.stderr
