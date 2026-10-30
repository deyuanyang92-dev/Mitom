"""The installed MitoM has no benchmark/, knowledge/, packaging/ or tests/ folder (bioconda: packaging/bioconda/build.sh;
the Docker image leaves out the same folders, .dockerignore): every command and every top-level module must work
without them. The bioconda build script itself installs MitoM into a temporary $PREFIX here."""
import concurrent.futures as cf
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "tests" / "data"
COMMANDS = ("run", "import", "annotate", "reorient", "transfer_metadata", "replace_organism", "dedup", "setup",
            "curate", "intron", "trna", "rescue", "rrna-hmm", "validate", "fetch", "check", "compare", "nuc", "plot",
            "debug")


@pytest.fixture(scope="module")
def prefix(tmp_path_factory):
    """$PREFIX after build.sh, run the way conda-build runs it (in the source folder, $PREFIX/bin/python present)."""
    pre = tmp_path_factory.mktemp("prefix")
    (pre / "bin").mkdir()
    (pre / "bin" / "python").symlink_to(sys.executable)
    subprocess.run(["bash", str(REPO / "packaging" / "bioconda" / "build.sh")], cwd=REPO, check=True,
                   env=dict(os.environ, PREFIX=str(pre)))
    yield pre
    shutil.rmtree(pre, ignore_errors=True)                    # models/ is large: not left in the pytest temp folder


def mitom(pre: Path, *args: str, cwd=None) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    return subprocess.run([str(pre / "bin" / "mitom"), *args], capture_output=True, text=True, env=env, cwd=cwd,
                          timeout=300)


def test_installed_tree_has_no_repo_only_folders(prefix):
    dest = prefix / "share" / "mitom"
    assert (dest / "mitom").is_file() and (dest / "mitom_pipeline.py").is_file()
    assert not any((dest / d).exists() for d in ("benchmark", "knowledge", "packaging", "tests"))


def test_every_command_starts(prefix):
    with cf.ThreadPoolExecutor(6) as ex:
        res = dict(zip(COMMANDS, ex.map(lambda c: mitom(prefix, c, "-h"), COMMANDS)))
    bad = {c: (p.returncode, p.stderr[-300:]) for c, p in res.items() if p.returncode != 0 or "usage" not in p.stdout}
    assert not bad


def test_every_module_imports_from_the_install(prefix):
    """Each top-level module imports in a Python that sees only the installed folder (not this repository)."""
    dest = prefix / "share" / "mitom"
    mods = sorted(p.stem for p in dest.glob("*.py"))
    code = ("import importlib, sys\nfrom pathlib import Path\n"
            f"mods = {mods!r}\n"
            "for m in mods:\n"
            "    mod = importlib.import_module(m)\n"
            "    assert Path(mod.__file__).resolve().parent == Path.cwd().resolve(), (m, mod.__file__)\n"
            "assert not {'mito_bench_eval', 'mito_bench_run'} & set(sys.modules)\n"
            "print('ok', len(mods))\n")
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    p = subprocess.run([str(prefix / "bin" / "python"), "-c", code], cwd=dest, capture_output=True, text=True,
                       env=env, timeout=300)
    assert p.returncode == 0, p.stderr[-2000:]
    assert p.stdout.strip() == f"ok {len(mods)}"


def test_check_and_compare_run_from_the_install(prefix, tmp_path):
    gb = DATA / "NC_000931.gb"
    p = mitom(prefix, "check", "-i", str(gb), "-o", str(tmp_path / "check"))
    assert p.returncode == 0, p.stderr
    assert (tmp_path / "check" / "check_report.tsv").is_file()
    out = tmp_path / "cmp.tsv"
    p = mitom(prefix, "compare", "--original", str(gb), "--mitom", str(gb), "-o", str(out), "--genetic_code", "5")
    assert p.returncode == 0, p.stderr
    rows = [ln.split("\t") for ln in out.read_text().splitlines()[1:]]
    assert rows and {r[-1] for r in rows} <= {"same", "same_place"}         # same_place: Leu / Ser without anticodon
