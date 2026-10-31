"""Inputs found in the stress test of 2026-10-31: alignment gaps, empty / all-N sequences, gzip and Latin-1 files, and a
run with nothing to annotate (it ended 'finished', exit 0, with 0 records).

Run: python3 -m pytest -q tests/test_input_hardening.py
"""
import gzip
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mitom_inputs as mi  # noqa: E402


def test_gaps_stops_and_odd_letters_are_cleaned():
    r, why = mi.clean_sequence(SeqRecord(Seq("ACGT-" * 200 + "**" + "ACGUX" * 10), id="a", description="a x"))
    assert "-" not in str(r.seq) and "*" not in str(r.seq) and "U" not in str(r.seq)
    assert str(r.seq).endswith("ACGTN" * 10) and r.description == "a x"
    assert "removed" in why and "read as N" in why


def test_clean_sequence_leaves_a_good_sequence_alone():
    rec = SeqRecord(Seq("acgt" * 290 + "RYKMSWBDHV" * 4), id="a")
    r, why = mi.clean_sequence(rec)
    assert r is rec and why == ""                     # soft-masked case and IUPAC codes kept


def test_empty_and_all_n_sequences_are_skipped(tmp_path):
    fa = tmp_path / "x.fasta"
    fa.write_text(">empty\n>allN\n" + "N" * 5000 + "\n>good\n" + "ACGT" * 300 + "\n")
    samples, _ = mi.split_or_sanitize_fastas([fa], tmp_path / "split", True, 15)
    assert [s.sample_name for s in samples] == ["good"]


def test_gzip_and_latin1_inputs_are_normalised(tmp_path):
    src = tmp_path / "in"
    src.mkdir()
    with gzip.open(src / "a.fasta.gz", "wt") as fh:
        fh.write(">a\n" + "ACGT" * 300 + "\n")
    (src / "b.fasta").write_bytes(b">caf\xe9\n" + b"ACGT" * 300 + b"\n")
    (src / "c.fasta").write_text(">c\n" + "ACGT" * 300 + "\n")
    out, notes = mi.normalize_input(src, tmp_path / "out")
    assert out == tmp_path / "out" / mi.NORMALIZED_DIR
    assert sorted(p.name for p in out.iterdir()) == ["a.fasta", "b.fasta", "c.fasta"]
    assert next(SeqIO.parse(str(out / "b.fasta"), "fasta")).id == "café"
    assert len(notes) == 2
    one, notes1 = mi.normalize_input(src / "a.fasta.gz", tmp_path / "out2")
    assert one.name == "a.fasta" and notes1
    plain, notes2 = mi.normalize_input(src / "c.fasta", tmp_path / "out3")
    assert plain == src / "c.fasta" and notes2 == []        # nothing to change: the input itself


def test_run_with_nothing_to_annotate_fails(tmp_path):
    fa = tmp_path / "allN.fasta"
    fa.write_text(">allN\n" + "N" * 16000 + "\n")
    p = subprocess.run([sys.executable, str(HERE.parent / "mitom_pipeline.py"), "run", "-i", str(fa), "-o",
                        str(tmp_path / "out"), "--clade", "Chordata", "--genetic_code", "2", "--plot", "no"],
                       capture_output=True, text=True)
    assert p.returncode == 2 and "no sequence to annotate" in p.stderr


def test_missing_input_is_a_clean_error(tmp_path):
    p = subprocess.run([sys.executable, str(HERE.parent / "mitom_pipeline.py"), "run", "-i", str(tmp_path / "nope.fa"),
                        "-o", str(tmp_path / "out")], capture_output=True, text=True)
    assert p.returncode != 0 and "Traceback" not in p.stderr and "not found" in p.stderr


def test_one_run_per_output_folder(tmp_path):
    """2026-10-31 stress test: two runs into one folder both failed with half the records. The second is refused while
    the first runs; a lock left by a killed run is taken over."""
    import os
    import mitom_pipeline as bm
    out = tmp_path / "out"
    out.mkdir()
    (out / bm.LOCK_NAME).write_text(f"{os.getppid()}\n")          # a running process (pytest's parent)
    assert "another MitoM run" in bm.lock_output(out)
    (out / bm.LOCK_NAME).write_text("999999999\n")                 # no such process: stale
    assert bm.lock_output(out) == ""
    assert (out / bm.LOCK_NAME).read_text().strip() == str(os.getpid())


def test_output_folder_not_writable_is_a_clean_error(tmp_path):
    import mitom_pipeline as bm
    ro = tmp_path / "ro"
    ro.mkdir()
    ro.chmod(0o555)
    try:
        assert "cannot write into the output folder" in bm.lock_output(ro / "out")
    finally:
        ro.chmod(0o755)


def test_utf8_byte_order_mark_is_removed(tmp_path):
    """2026-10-31 review: a Notepad 'UTF-8 with BOM' FASTA gave 0 records."""
    fa = tmp_path / "bom.fasta"
    fa.write_bytes(b"\xef\xbb\xbf>s1\n" + b"ACGT" * 300 + b"\n")
    out, notes = mi.normalize_input(fa, tmp_path / "out")
    assert notes and [r.id for r in SeqIO.parse(str(out), "fasta")] == ["s1"]


def test_protein_fasta_is_not_annotated_as_dna():
    """2026-10-31 fuzz: a protein FASTA (A G K M V ... are IUPAC nucleotide codes too) ran as a genome."""
    r, why = mi.clean_sequence(SeqRecord(Seq("MKVLAAGIVGLLLAAQPAMA" * 100), id="p"))
    assert len(r.seq) == 0 and "not a nucleotide sequence" in why


def test_embl_input_is_written_as_genbank(tmp_path):
    rec = SeqRecord(Seq("ACGT" * 300), id="X1.1", name="X1", description="test")
    rec.annotations["molecule_type"] = "DNA"
    SeqIO.write(rec, str(tmp_path / "a.embl"), "embl")
    out, notes = mi.normalize_input(tmp_path / "a.embl", tmp_path / "out")
    assert out.suffix == ".gb" and "EMBL" in notes[0]
    assert len(SeqIO.read(str(out), "genbank").seq) == 1200


def test_a_previous_runs_work_folders_are_cleared(tmp_path):
    """2026-10-31 review: a stale 04.final_gb/*.gbf of an earlier run into the same -o was published as a result,
    and a stale reoriented body of another genome became a final record."""
    import mitom_output as mo
    for d in ("04.final_gb", "01.anno_round1/S00001/reoriented", "06.ncbi_submission", "00.input", "my_notes"):
        (tmp_path / d).mkdir(parents=True)
        (tmp_path / d / "x.gbf").write_text("LOCUS STALE\n//\n")
    mo.clear_previous_work(tmp_path, keep=("00.input",))
    assert sorted(p.name for p in tmp_path.iterdir()) == ["00.input", "my_notes"]


def _gb_text(name, n, seed):
    import io
    import random
    r = SeqRecord(Seq("".join(random.Random(seed).choice("ACGT") for _ in range(n))), id=f"{name}.1", name=name,
                  description="test")
    r.annotations.update(molecule_type="DNA", topology="circular", accessions=[name], sequence_version=1)
    h = io.StringIO()
    SeqIO.write(r, h, "genbank")
    return h.getvalue()


def test_a_cut_off_last_record_is_not_annotated(tmp_path):
    """2026-10-31 review: a download killed while writing left a record without its end; it was read as a complete
    circular genome (11,840 of 15,000 bp)."""
    full = _gb_text("MN000001", 3000, 1) + _gb_text("MN000002", 3000, 2)
    cut = full[: len(full) - 900]
    (tmp_path / "x.gb").write_text(cut)
    kept, _, stats = mi.dedup_gb_records([tmp_path / "x.gb"])
    assert [m.acc_versioned for m in kept] == ["MN000001.1"]


def test_fetch_resume_downloads_a_cut_off_record_again(tmp_path):
    import mito_fetch
    full = _gb_text("MN000001", 3000, 1) + _gb_text("MN000002", 3000, 2)
    (tmp_path / "x.gb").write_text(full[: len(full) - 900])
    assert list(mito_fetch.existing_accessions(tmp_path / "x.gb")) == ["MN000001.1"]


def test_contig_record_without_sequence_is_skipped():
    from Bio.Seq import UndefinedSequenceError  # noqa: F401
    r = SeqRecord(Seq(None, length=1000), id="CONTIG1")
    out, why = mi.clean_sequence(r)
    assert len(out.seq) == 0 and "CONTIG" in why


def test_gb2fasta_same_file_name_in_two_subfolders(tmp_path):
    for d, seed in (("a", 1), ("b", 2)):
        (tmp_path / "in" / d).mkdir(parents=True)
        (tmp_path / "in" / d / "mt.gb").write_text(_gb_text(f"MN00000{seed}", 1200, seed))
    p = subprocess.run([sys.executable, str(HERE.parent / "mitom_pipeline.py"), "gb2fasta", "-i", str(tmp_path / "in"),
                        "-o", str(tmp_path / "out")], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert sorted(x.name for x in (tmp_path / "out").iterdir()) == ["mt.fasta", "mt_2.fasta"]
