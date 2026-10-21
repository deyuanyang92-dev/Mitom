"""Tests for mito_fetch.py (NCBI download). E-utilities are replaced by fakes; no network.

Run: python3 -m pytest -q tests/test_mito_fetch.py
"""
import io
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_fetch as mf  # noqa: E402

GB = (HERE / "data" / "NC_000931.gb").read_text()


def test_build_query_taxon_and_filters():
    q = mf.build_query(6340, None, False)
    assert q.startswith("txid6340[Organism:exp] AND mitochondrion[Filter]")
    assert '"complete genome"[Title]' in q and "NOT partial[Title]" in q and "srcdb_refseq" not in q
    assert mf.build_query(6340, None, True).endswith(" AND srcdb_refseq[PROP]")
    assert mf.build_query(6340, None, False, 12000, 20000).endswith("12000:20000[SLEN]")
    assert mf.build_query(6340, None, False, 12000, 0).endswith("12000:999999999[SLEN]")


def test_build_query_raw_overrides_taxon():
    q = mf.build_query(None, "Nereididae[Organism] AND mitochondrion[Filter]", True)
    assert q == "(Nereididae[Organism] AND mitochondrion[Filter]) AND srcdb_refseq[PROP]"
    with pytest.raises(ValueError):
        mf.build_query(None, None, False)


def test_resolve_taxid_numeric_and_name(monkeypatch):
    assert mf.resolve_taxid("6340") == 6340
    calls = []

    def esearch(db, term, **kw):
        calls.append((db, term))
        return io.StringIO("")

    monkeypatch.setattr(mf.Entrez, "esearch", esearch)
    monkeypatch.setattr(mf.Entrez, "read", lambda h: {"IdList": ["6340"]})
    assert mf.resolve_taxid("Annelida") == 6340
    assert calls == [("taxonomy", '"Annelida"[Scientific Name]')]
    monkeypatch.setattr(mf.Entrez, "read", lambda h: {"IdList": []})
    with pytest.raises(SystemExit, match="not found"):
        mf.resolve_taxid("Nonexistus")
    monkeypatch.setattr(mf.Entrez, "read", lambda h: {"IdList": ["1", "2"]})
    with pytest.raises(SystemExit, match="give the taxid"):
        mf.resolve_taxid("Ambiguous")


def test_search_uids_pages(monkeypatch):
    pages = {0: ["A.1", "B.1"], 2: ["C.1"]}
    seen = []

    def esearch(db, term, retstart, retmax, idtype):
        seen.append(retstart)
        return pages[retstart]

    monkeypatch.setattr(mf.Entrez, "esearch", esearch)
    monkeypatch.setattr(mf.Entrez, "read", lambda ids: {"IdList": ids, "Count": "3"})
    assert mf.search_uids("x") == ["A.1", "B.1", "C.1"]
    assert seen == [0, 2]


def test_download_appends_and_counts(tmp_path, monkeypatch):
    asked = []

    def efetch(db, id, rettype, retmode):
        asked.append(id.split(","))
        return io.StringIO(GB * len(id.split(",")))

    monkeypatch.setattr(mf.Entrez, "efetch", efetch)
    monkeypatch.setattr(mf, "GB_BATCH", 2)
    gb = tmp_path / "x.gb"
    assert mf.download(["A", "B", "C"], gb) == 3
    assert asked == [["A", "B"], ["C"]]
    assert gb.read_text().count("\n//\n") == 3


def test_existing_accessions_and_report(tmp_path):
    gb = tmp_path / "x.gb"
    gb.write_text(GB)
    assert mf.existing_accessions(gb) == {"NC_000931.1": "NC_000931"}
    assert mf.existing_accessions(tmp_path / "missing.gb") == {}
    rep = tmp_path / "r.tsv"
    assert mf.write_report(gb, rep) == 1
    lines = rep.read_text().splitlines()
    assert lines[0].split("\t")[:4] == ["accession", "refseq", "organism", "length"]
    f = lines[1].split("\t")
    assert f[0] == "NC_000931.1" and f[1] == "yes" and f[2] == "Platynereis dumerilii" and f[4] == "circular"
    assert "Annelida" in f[6]


def test_retry_then_raise(monkeypatch):
    monkeypatch.setattr(mf.time, "sleep", lambda s: None)
    n = {"calls": 0}

    def flaky():
        n["calls"] += 1
        if n["calls"] < 3:
            raise ConnectionError("reset")
        return "ok"

    assert mf._retry(flaky, "x") == "ok" and n["calls"] == 3

    def dead():
        raise TimeoutError("never")

    with pytest.raises(TimeoutError):
        mf._retry(dead, "y")


def test_main_resume_skips_present_records(tmp_path, monkeypatch):
    """Second run only fetches the records that are not in the file; --max caps new ones."""
    monkeypatch.setattr(mf, "setup_entrez", lambda e, k: None)
    monkeypatch.setattr(mf, "search_uids", lambda term: ["NC_000931.1", "X1.1", "X2.1"])
    got = []
    monkeypatch.setattr(mf, "download", lambda uids, gb: got.append(uids) or 0)
    monkeypatch.setattr(mf, "write_report", lambda gb, rep: 0)
    out = tmp_path / "o"
    out.mkdir()
    (out / "ncbi.gb").write_text(GB)
    assert mf.main(["--query", "q", "-o", str(out), "--max", "1"]) == 0
    assert got == [["X1.1"]]
    assert (out / "query.txt").read_text().strip() == "(q)"


def test_mitom_lists_fetch():
    text = (HERE.parent / "mitom").read_text()
    assert '"fetch": ("mito_fetch.py", [])' in text


def test_download_retries_incomplete_batch(tmp_path, monkeypatch):
    """Review 2026-10: a truncated response (no terminating //) is not written; the batch is asked again."""
    monkeypatch.setattr(mf.time, "sleep", lambda s: None)
    calls = []

    def efetch(db, id, rettype, retmode):
        calls.append(id)
        return io.StringIO(GB[:-40] if len(calls) == 1 else GB)

    monkeypatch.setattr(mf.Entrez, "efetch", efetch)
    gb = tmp_path / "x.gb"
    assert mf.download(["A"], gb) == 1 and len(calls) == 2
    text = gb.read_text()
    assert text.count("LOCUS ") == 1 and text.count("\n//\n") == 1


def test_download_error_page_is_never_written(tmp_path, monkeypatch):
    monkeypatch.setattr(mf.time, "sleep", lambda s: None)
    monkeypatch.setattr(mf, "RETRIES", 2)
    monkeypatch.setattr(mf.Entrez, "efetch", lambda db, id, rettype, retmode: io.StringIO("Error: Invalid ID\n"))
    gb = tmp_path / "x.gb"
    with pytest.raises(RuntimeError):
        mf.download(["A"], gb)
    assert gb.read_text() == ""


def test_main_keeps_the_present_version_of_a_record(tmp_path, monkeypatch):
    """NC_000931.1 in the file, NCBI now serves NC_000931.2: not downloaded again (two versions of one genome)."""
    monkeypatch.setattr(mf, "setup_entrez", lambda e, k: None)
    monkeypatch.setattr(mf, "search_uids", lambda term: ["NC_000931.2", "X1.1"])
    got = []
    monkeypatch.setattr(mf, "download", lambda uids, gb: got.append(uids) or 0)
    monkeypatch.setattr(mf, "write_report", lambda gb, rep: 0)
    out = tmp_path / "o"
    out.mkdir()
    (out / "ncbi.gb").write_text(GB)
    assert mf.main(["--query", "q", "-o", str(out)]) == 0
    assert got == [["X1.1"]]


def test_main_without_matches_writes_an_empty_report(tmp_path, monkeypatch):
    monkeypatch.setattr(mf, "setup_entrez", lambda e, k: None)
    monkeypatch.setattr(mf, "search_uids", lambda term: [])
    out = tmp_path / "o"
    assert mf.main(["--query", "q", "-o", str(out)]) == 0
    assert (out / "fetch_report.tsv").read_text().startswith("accession\t")
