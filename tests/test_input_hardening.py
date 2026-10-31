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
    rec = SeqRecord(Seq("acgtRYKM" * 300), id="a")
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
