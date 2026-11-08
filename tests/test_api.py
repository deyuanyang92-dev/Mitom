"""In-process API of the leaf commands: <module>.run(<module>.options(...)) does what the command line does
(mito_common.options: defaults of the parser, required options checked)."""
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_curate  # noqa: E402
import mito_intron  # noqa: E402
import mito_intron_ref  # noqa: E402
import mito_nuc  # noqa: E402

DATA = HERE / "data"


def test_options_are_the_command_line_defaults():
    a = mito_curate.options(input="x", out="y")
    b = mito_curate.build_parser().parse_args(["-i", "x", "-o", "y"])
    assert vars(a) == vars(b)
    a = mito_intron.options(input="x", out="y")
    assert vars(a) == vars(mito_intron.build_parser().parse_args(["-i", "x", "-o", "y"]))
    a = mito_intron_ref.options("annotate", input="x", out="y", jobs=2)
    assert vars(a) == vars(mito_intron_ref.build_parser().parse_args(["annotate", "-i", "x", "-o", "y", "--jobs", "2"]))
    a = mito_nuc.options("annotate", input=Path("x"), out=Path("y"))
    assert vars(a) == vars(mito_nuc.build_parser().parse_args(["annotate", "-i", "x", "-o", "y"]))


def test_options_refuse_unknown_and_missing():
    with pytest.raises(TypeError, match="unknown option"):
        mito_curate.options(input="x", out="y", no_such=1)
    with pytest.raises(TypeError, match="missing out"):
        mito_curate.options(input="x")


def test_curate_in_process_equals_the_command(tmp_path):
    src = DATA / "OW387151_mitoz.gbf"
    assert mito_curate.run(mito_curate.options(input=str(src), out=str(tmp_path / "api"), no_validate=True)) == 0
    subprocess.run([sys.executable, "-W", "ignore", str(HERE.parent / "mito_curate.py"), "-i", str(src), "-o",
                    str(tmp_path / "cli"), "--no_validate"], check=True, capture_output=True)
    for name in ("OW387151_mitoz.gbf", "curation_report.tsv"):
        assert (tmp_path / "api" / name).read_text() == (tmp_path / "cli" / name).read_text()
