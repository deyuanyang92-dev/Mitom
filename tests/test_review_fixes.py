"""Review fixes 2026-10-22 (mitom_pipeline / mito_env): origin joins after the shift, core failures on the console,
native core without MitoZ / MITOS2 messages, cox1_first of linear records, run.log per run, wrong --genetic_code,
multi-contig FASTA topology, $YZZ_MITO_CONFIG for a new file.

Run: python3 -m pytest -q tests/test_review_fixes.py
"""
import logging
import os
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mitom_pipeline as bm  # noqa: E402
import mito_env  # noqa: E402


@pytest.fixture(autouse=True)
def restore_environ():
    old = os.environ.copy()
    yield
    os.environ.clear()
    os.environ.update(old)


def parse(tmp_path, *extra, inp=None):
    return bm.build_parser().parse_args(["run", "-i", str(inp or tmp_path), "-o", str(tmp_path / "out"), *extra])


# ---- shift_location: a join across the old origin becomes one interval ---------------------------------------
@pytest.mark.parametrize("strand", [1, -1])
def test_origin_join_shifted_past_the_cut_is_one_interval(strand):
    """join(16000..16560,1..300) on a 16560-nt circular record, shifted by 1000: one simple interval
    15000..15860 (1-based; 0-based [14999, 15860)), both strands."""
    parts = [FeatureLocation(15999, 16560, strand), FeatureLocation(0, 300, strand)]
    loc = CompoundLocation(parts[::-1] if strand == -1 else parts)          # Biopython: parts in transcript order
    new = bm.shift_location(loc, 1000, 16560)
    assert isinstance(new, FeatureLocation) and not isinstance(new, CompoundLocation)
    assert (int(new.start), int(new.end), new.strand) == (14999, 15860, strand)


def test_true_exons_stay_separate_and_outer_fuzzy_ends_are_kept():
    loc = CompoundLocation([FeatureLocation(100, 200, 1), FeatureLocation(300, 400, 1)])
    new = bm.shift_location(loc, 50, 16560)
    assert len(new.parts) == 2 and [(int(p.start), int(p.end)) for p in new.parts] == [(50, 150), (250, 350)]
    merged = bm.merge_abutting_parts([FeatureLocation(BeforePosition(10), 20, 1),
                                      FeatureLocation(20, AfterPosition(30), 1)])
    assert isinstance(merged.start, BeforePosition) and isinstance(merged.end, AfterPosition)
    assert (int(merged.start), int(merged.end)) == (10, 30)


def test_rotated_record_keeps_the_origin_gene_sequence():
    L = 16560
    seq = Seq("".join("ACGT"[i % 4] for i in range(L)))
    r = SeqRecord(seq, id="X", name="X", annotations={"molecule_type": "DNA", "topology": "circular"})
    for strand in (1, -1):
        parts = [FeatureLocation(15999, L, strand), FeatureLocation(0, 300, strand)]
        r.features.append(SeqFeature(CompoundLocation(parts[::-1] if strand == -1 else parts), type="CDS",
                                     qualifiers={"gene": ["nad5" if strand == -1 else "cox1"]}))
    rot = bm.rotate_record(r, 1000)
    for old, new in zip(r.features, rot.features):
        assert not isinstance(new.location, CompoundLocation)
        assert str(old.extract(r.seq)) == str(new.extract(rot.seq))


# ---- core failures reach the console and the exit code -------------------------------------------------------
def fake_core(tmp_path, rc, stderr="boom: no refseq89m"):
    p = tmp_path / "fake_core.py"
    p.write_text(f"import sys\nprint({stderr!r}, file=sys.stderr)\nsys.exit({rc})\n")
    return p


def test_core_failure_is_warned_on_the_root_logger(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(bm, "MITO_CORE", fake_core(tmp_path, 3))
    a = parse(tmp_path, "--clade", "Mollusca", "--genetic_code", "5")
    fa = tmp_path / "S00001.fasta"
    fa.write_text(">S00001 topology=circular\nACGT\n")
    si = bm.SampleInfo("sampleA", "sampleA", "S00001", "circular", fa, fa)
    with caplog.at_level(logging.WARNING, logger=""):
        iid, ok, *_ = bm.annotate_one(si, fa, tmp_path / "out", "01.anno_round1", a)
    assert not ok
    msgs = [r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING]
    assert any(m.startswith("[FAIL] sampleA: annotation core exit 3: boom: no refseq89m") for m in msgs), msgs
    assert not any("MitoZ" in m for m in msgs)


def test_run_returns_nonzero_when_every_record_failed(tmp_path, monkeypatch):
    monkeypatch.setattr(bm, "MITO_CORE", fake_core(tmp_path, 3))
    monkeypatch.setenv("YZZ_MITO_NO_SCAN", "1")
    inp = tmp_path / "in.fa"
    inp.write_text(">rec1 topology=circular\n" + "ACGT" * 4000 + "\n")
    a = parse(tmp_path, "--clade", "Mollusca", "--genetic_code", "5", "--rescue", "no", "--rrna_hmm", "no",
              "--intron_fix", "no", "--curate", "no", inp=inp)
    assert bm.resolve_tools(a) == ""
    assert bm.workflow_run(a) != 0
    log = (tmp_path / "out" / "run.log").read_text()
    assert "[FAIL] rec1: annotation core exit 3" in log and "[FAIL] every record failed" in log


# ---- resolve_tools under the native core --------------------------------------------------------------------
def test_native_core_does_not_complain_about_optional_mitos2(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(mito_env, "mitos_bin", lambda explicit=None: None)
    monkeypatch.setattr(mito_env, "mitos_refdir", lambda explicit=None: None)
    a = parse(tmp_path, "--clade", "Mollusca", "--genetic_code", "5")
    assert bm.resolve_tools(a) == ""
    err = capsys.readouterr().err
    assert "MitoZ" not in err and "MITOS2 not found" not in err


# ---- cox1_first of a linear record ---------------------------------------------------------------------------
def record(topology, cox1_at, L=3000):
    r = SeqRecord(Seq("".join("ACGT"[i % 4] for i in range(L))), id="A", name="A",
                  annotations={"molecule_type": "DNA", "topology": topology})
    r.features = [SeqFeature(FeatureLocation(100, 160, 1), type="tRNA", qualifiers={"gene": ["trnF"]}),
                  SeqFeature(FeatureLocation(cox1_at, cox1_at + 600, 1), type="CDS", qualifiers={"gene": ["COX1"]})]
    return r


def test_linear_record_is_not_rotated_and_summary_says_so(tmp_path, caplog):
    src = tmp_path / "A.gbf"
    SeqIO.write(record("linear", 1000), str(src), "genbank")
    bm._LINEAR_SKIP_SEEN.clear()
    with caplog.at_level(logging.INFO):
        pcg, trna, rrna, mode, ok = bm.finalize_genbank_output(
            sample_name="A", body_src_gb=src, metadata_src_gb=None, dest_gb=tmp_path / "out" / "A.gbf",
            which_gene_first="cox1", force_forward=True, preserve_metadata=False)
        bm.finalize_genbank_output(                                             # the run's second reorientation
            sample_name="A", body_src_gb=src, metadata_src_gb=None, dest_gb=tmp_path / "out" / "A2.gbf",
            which_gene_first="cox1", force_forward=True, preserve_metadata=False)
    assert ok is False
    si = bm.SampleInfo("A", "A", "S00001", "linear", src, src)
    assert bm.orientation_text(ok, si) == "no (linear, not rotated)"
    assert bm.orientation_text(False, bm.SampleInfo("A", "A", "S00001", "circular", src, src)) == "no"
    skips = [r for r in caplog.records if "[LINEAR_REORIENT_SKIP]" in r.getMessage()]
    assert len(skips) == 2 and [r.levelno for r in skips] == [logging.WARNING, logging.INFO]
    # a circular record is rotated: cox1 first
    SeqIO.write(record("circular", 1000), str(src), "genbank")
    *_, ok = bm.finalize_genbank_output(
        sample_name="A", body_src_gb=src, metadata_src_gb=None, dest_gb=tmp_path / "out" / "A3.gbf",
        which_gene_first="cox1", force_forward=True, preserve_metadata=False)
    assert ok is True
    # a linear record that already starts with cox1
    SeqIO.write(record("linear", 0), str(src), "genbank")
    *_, ok = bm.finalize_genbank_output(
        sample_name="A", body_src_gb=src, metadata_src_gb=None, dest_gb=tmp_path / "out" / "A4.gbf",
        which_gene_first="cox1", force_forward=True, preserve_metadata=False)
    assert ok is True


# ---- run.log is written anew for each run --------------------------------------------------------------------
def test_run_log_is_not_appended_across_runs(tmp_path):
    log = tmp_path / "run.log"
    for _ in range(2):
        bm.setup_root_logger(log, console_level=logging.WARNING, mode="w")
        logging.info("[PARAMS] x")
        logging.getLogger().handlers[0].flush()
    assert log.read_text().count("[PARAMS]") == 1
    bm.setup_root_logger(tmp_path / "pipeline.log")
    for h in logging.getLogger().handlers:
        logging.getLogger().removeHandler(h)


# ---- a wrong --genetic_code is said ----------------------------------------------------------------------------
def test_transl_table_of_the_input_differing_from_the_code_is_warned(tmp_path, caplog):
    orig = tmp_path / "original_annotation"
    orig.mkdir()
    r = SeqRecord(Seq("ATG" * 100), id="NC_1.1", name="NC_1", annotations={"molecule_type": "DNA"})
    r.features = [SeqFeature(FeatureLocation(0, 300, 1), type="CDS", qualifiers={"gene": ["COX1"], "transl_table": ["5"]})]
    SeqIO.write(r, str(orig / "NC_1.1.gb"), "genbank")
    with caplog.at_level(logging.WARNING):
        assert bm.check_transl_tables(orig, 2) == ["NC_1.1"]
        assert bm.check_transl_tables(orig, 5) == []
    assert any("[CODE] NC_1.1: record transl_table=5 differs from --genetic_code 2" in rec.getMessage()
               for rec in caplog.records)


# ---- multi-contig / short FASTA records are linear unless the header says circular -------------------------
def test_contigs_and_short_records_are_linear_unless_declared_circular(tmp_path, caplog):
    fa = tmp_path / "contigs.fa"
    fa.write_text(">c1\n" + "ACGT" * 4000 + "\n>c2\n" + "ACGT" * 500 + "\n>c3 topology=circular\n" + "ACGT" * 500 + "\n")
    one = tmp_path / "one.fa"
    one.write_text(">g1\n" + "ACGT" * 4000 + "\n")
    short = tmp_path / "short.fa"
    short.write_text(">g2\n" + "ACGT" * 500 + "\n")
    with caplog.at_level(logging.WARNING):
        samples, _ = bm.split_or_sanitize_fastas([fa, one, short], tmp_path / "split", True, 16)
    topo = {si.sample_name: si.topology for si in samples}
    assert topo == {"c1": "linear", "c2": "linear", "c3": "circular", "g1": "circular", "g2": "linear"}
    msgs = [r.getMessage() for r in caplog.records if r.getMessage().startswith("[TOPOLOGY]")]
    assert len(msgs) == 3 and "[TOPOLOGY] c1: 3 records / 16000 nt — treated as linear" in msgs[0]
    assert ">c1 topology=linear" in (tmp_path / "split" / "S00001.fasta").read_text().replace(">S00001", ">c1")


# ---- $YZZ_MITO_CONFIG names the file to write even before it exists ------------------------------------------
def test_env_config_path_is_used_before_the_file_exists(tmp_path, monkeypatch):
    new = tmp_path / "new" / "mito.json"
    monkeypatch.setenv("YZZ_MITO_CONFIG", str(new))
    assert mito_env.config_path() == new
    monkeypatch.delenv("YZZ_MITO_CONFIG")
    assert mito_env.config_path() != new


def test_rescue_auto_skipped_for_native_core(tmp_path, monkeypatch):
    """--rescue auto does not run MITOS2: under the native core it cost 5-7 min and ~1 GB per record lacking a gene
    (mostly atp8 of taxa without atp8) for 2 genes in 24 records."""
    import argparse
    import mitom_pipeline as bm
    import mito_rescue
    called = []
    monkeypatch.setattr(mito_rescue, "rescue_dir", lambda *a, **k: called.append(1) or [])
    a = argparse.Namespace(rescue="auto", core="native", mitos_bin="x", mitos_refdir="y", genetic_code=5, max_tasks=1)
    bm.run_rescue(a, tmp_path, tmp_path)
    assert not called
    a.rescue = "yes"
    bm.run_rescue(a, tmp_path, tmp_path)
    assert called
