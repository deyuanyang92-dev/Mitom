"""Codex release-readiness review 2026-11-08 (Mitom_review/RELEASE_READINESS_2026-11-08.md): each confirmed bug.

Run: python3 -m pytest -q tests/test_review_2026_11_08.py
"""
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
import mito_engines as me  # noqa: E402
import mito_env  # noqa: E402


@pytest.mark.parametrize("bad", ["-c2", "-o/tmp/x", "-rrefseq63m", "-R/refs", "-iother.fa", "--cod 2", "--outd x",
                                 "--refd /x", "--lin", "--best -c2", "--json=p.json"])
def test_1_mitos_args_attached_or_abbreviated_own_option_refused(bad):
    """-c2 (attached value) and --cod (argparse abbreviation) reached runmitos and overrode MitoM's genetic code."""
    with pytest.raises(ValueError):
        me.mitos_extra_args(bad)


def test_1_mitos_args_other_options_pass():
    assert me.mitos_extra_args("--best --circrot 90 --cutoff 0.3 --ncbicode") == [
        "--best", "--circrot", "90", "--cutoff", "0.3", "--ncbicode"]


def test_2_runmitos_py_found_by_every_lookup(tmp_path):
    """bioconda mitos 2.1.9 has only runmitos.py: find() (doctor, debug, mito_trna) and mitos_bin() agree."""
    b = tmp_path / "bin"
    b.mkdir()
    (b / "runmitos.py").write_text("#!/bin/sh\n")
    (b / "runmitos.py").chmod(0o755)
    assert mito_env.find("runmitos", str(b)) == str(b / "runmitos.py")
    assert mito_env.mitos_bin(str(b)) == str(b)


def test_3_gb2tbl_reads_a_folder_of_gb_gz(tmp_path):
    """gb2tbl said 'no GenBank / EMBL file' for a folder holding only .gb.gz; output named without .gb."""
    src = tmp_path / "in"
    src.mkdir()
    (src / "NC_000931.gb.gz").write_bytes(gzip.compress((HERE / "data" / "NC_000931.gb").read_bytes()))
    (src / "readme.txt.gz").write_bytes(gzip.compress(b"x"))           # not GenBank: ignored
    out = tmp_path / "out"
    p = subprocess.run([sys.executable, str(ROOT / "mito_tbl.py"), "-i", str(src / "NC_000931.gb.gz"), "-o", str(out)],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert (out / "NC_000931.tbl").is_file() and (out / "NC_000931.fsa").is_file()
    out2 = tmp_path / "out2"
    p = subprocess.run([sys.executable, str(ROOT / "mito_tbl.py"), "-i", str(src), "-o", str(out2)],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert "COX1" in (out2 / "submission.tbl").read_text() or "cox1" in (out2 / "submission.tbl").read_text().lower()
