"""mito_env.py (one place that finds programs), clade check, MITOS2 mode, run defaults."""
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import batch_mitoz as bm  # noqa: E402
import mito_env  # noqa: E402


def fake_exe(d: Path, name: str) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    p = d / name
    p.write_text("#!/bin/sh\nexit 0\n")
    p.chmod(0o755)
    return p


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    """Isolated config file; caches cleared before and after."""
    path = tmp_path / "mito.json"
    monkeypatch.setenv("YZZ_MITO_CONFIG", str(path))
    mito_env.config.cache_clear()
    yield path
    mito_env.config.cache_clear()


def test_config_tools_bin_used(cfg, tmp_path):
    fake = fake_exe(tmp_path / "tools", "exonerate-yzz-test")
    cfg.write_text('{"tools_bin": "%s"}' % fake.parent)
    assert mito_env.find("exonerate-yzz-test") == str(fake)


def test_explicit_wrong_path_is_not_replaced(cfg, tmp_path):
    # 2026-10-20: an explicit --mitos_bin that does not contain runmitos must fail, not silently use another MITOS2
    fake_exe(tmp_path / "tools", "exonerate-yzz-test")
    cfg.write_text('{"tools_bin": "%s"}' % (tmp_path / "tools"))
    assert mito_env.find("exonerate-yzz-test", str(tmp_path / "empty")) is None
    assert mito_env.bin_with("exonerate-yzz-test", explicit=str(tmp_path / "empty")) is None


def test_bin_with_needs_all_programs(cfg, tmp_path):
    d = tmp_path / "half"
    fake_exe(d, "prog-a-yzz")
    assert mito_env.bin_with("prog-a-yzz", "prog-b-yzz") is None
    fake_exe(d, "prog-b-yzz")
    cfg.write_text('{"tools_bin": "%s"}' % d)
    mito_env.config.cache_clear()
    assert mito_env.bin_with("prog-a-yzz", "prog-b-yzz") == str(d)


def test_wrong_clade_rejected_before_any_job(tmp_path):
    # 2026-10-20: an unknown --clade only failed per sample after setup (Cnidaria is now an extended clade)
    r = subprocess.run([sys.executable, str(HERE.parent / "batch_mitoz.py"), "run", "-i", str(tmp_path),
                        "-o", str(tmp_path / "out"), "--clade", "Cnidarians"], capture_output=True, text=True)
    assert r.returncode == 2 and "not a supported clade" in r.stderr and not (tmp_path / "out").exists()


def test_mitos_mode_is_lowercased_and_validated(capsys):
    """C18 (scripts/mitoz_fast.py:54): YZZ_MITOS_MODE=Fast was treated as 'sensitive' without a word."""
    assert mito_env.mitos_mode("Fast") == "fast" and mito_env.mitos_mode(" FULL ") == "full"
    assert mito_env.mitos_mode("sensitive") == "sensitive" and mito_env.mitos_mode(None) == "sensitive"
    assert mito_env.mitos_mode("") == "sensitive"
    assert capsys.readouterr().err == ""
    assert mito_env.mitos_mode("sensetive") == "sensitive"            # unknown -> warning + sensitive
    err = capsys.readouterr().err
    assert "YZZ_MITOS_MODE" in err and "sensetive" in err and "sensitive" in err


def test_mito_trna_rejects_unknown_mitos_mode(tmp_path):
    """C18: any unknown mode silently behaved as 'sensitive'."""
    import types
    import mito_trna
    with pytest.raises(ValueError, match="mitos_mode"):
        mito_trna.mitos2(tmp_path / "g.fa", 5, {"runmitos": "/x/runmitos"}, tmp_path,
                         types.SimpleNamespace(mitos_mode="sensetive", mitos_ref="r", mitos_refdir="d"))


def test_extended_clade_accepted_and_needs_code(tmp_path):
    assert "Animal" in bm.ALL_CLADES and "Cnidaria" in bm.ALL_CLADES
    r = subprocess.run([sys.executable, str(HERE.parent / "batch_mitoz.py"), "run", "-i", str(tmp_path), "-o", str(tmp_path / "o"),
                        "--clade", "Cnidaria"], capture_output=True, text=True)
    assert r.returncode != 0 and "--genetic_code" in (r.stdout + r.stderr)


def test_native_core_only_and_mitom_steps_on_by_default(capsys):
    """MitoZ removed (user 2026-10-28: only its annotation and circular map were wanted, both now MitoM's own):
    --core mitoz / hybrid are refused, --core native (old scripts) is accepted. Reorientation, round 2, metadata
    transfer, introns, curation and the genome maps are on by default; the MITOS2 rescue is off (auto)."""
    p = bm.build_parser()
    for core in ("mitoz", "hybrid"):
        with pytest.raises(SystemExit):
            p.parse_args(["run", "-i", "x", "-o", "y", "--core", core])
        assert "invalid choice" in capsys.readouterr().err
    a = p.parse_args(["run", "-i", "x", "-o", "y", "--clade", "Mollusca", "--core", "native"])
    assert bm.resolve_tools(a) == "" and (a.core, a.rescue, a.plot, a.curate) == ("native", "auto", "yes", "yes")
    assert (a.which_gene_first, a.if_reannotation_2, a.preserve_phase1_metadata, a.intron_fix) == ("cox1", "auto", "yes", "auto")
    assert not hasattr(a, "trna") and not hasattr(a, "circos") and not hasattr(a, "conda_env_mitoz")


@pytest.mark.skipif(not mito_env.find("table2asn"), reason="table2asn not available")
def test_run_ends_with_curated_and_submission_files(tmp_path):
    """User 2026-10-21: output must be submittable to NCBI. run = MitoM annotation (04.final_gb) -> curation with
    NCBI names (05.curated) -> per record .fsa/.tbl (BankIt) + .sqn/.val (table2asn) in 06.ncbi_submission."""
    import argparse
    final = tmp_path / "04.final_gb"
    final.mkdir()
    (final / "KY753833.gbf").write_bytes((HERE / "data" / "KY753833_mitoz.gbf").read_bytes())
    args = argparse.Namespace(curate="yes", genetic_code=5, sbt=None)
    merged_from = bm.run_curate_and_submit(args, final, tmp_path)
    assert merged_from == tmp_path / "05.curated"
    rec = SeqIO.read(str(merged_from / "KY753833.gbf"), "genbank")
    assert {f.qualifiers["gene"][0] for f in rec.features if f.type == "rRNA"} == {"rrnS", "rrnL"}
    sub = tmp_path / "06.ncbi_submission"
    rid = rec.id
    assert all((sub / f"{rid}{ext}").exists() for ext in (".fsa", ".tbl", ".sqn", ".val"))
    head = (sub / "submission_summary.tsv").read_text().splitlines()
    assert head[0].startswith("record\tERROR\tWARNING") and head[1].startswith(rid)


def test_doctor_has_required_biopython_row(capsys, monkeypatch):
    """B1: the doctor said `All required programs found` although the scripts could not import Biopython."""
    import Bio
    mito_env.doctor()
    out = capsys.readouterr().out
    assert any(line.lstrip().startswith("OK") and "Biopython" in line and Bio.__version__ in line
               for line in out.splitlines())
    monkeypatch.setitem(sys.modules, "Bio", None)               # import Bio -> ImportError
    assert mito_env.doctor() == 1
    out = capsys.readouterr().out
    assert any(line.lstrip().startswith("!!") and "Biopython" in line for line in out.splitlines())
    assert "pip install biopython" in out and "All required programs found" not in out


@pytest.mark.parametrize("cpus,threads,tasks", [(4, 4, 1), (32, 8, 4), (2, 2, 1), (8, 8, 1), (None, 8, 1)])
def test_default_threads_follow_the_machine(monkeypatch, cpus, threads, tasks):
    """B22: 16 threads x 16 jobs regardless of the machine; now threads = min(8, CPUs), jobs = CPUs // threads."""
    monkeypatch.setattr(bm.os, "cpu_count", lambda: cpus)
    a = bm.build_parser().parse_args(["run", "-i", "x", "-o", "y"])
    assert (a.threads, a.max_tasks) == (threads, tasks)
    a = bm.build_parser().parse_args(["run", "-i", "x", "-o", "y", "--threads", "3", "--max_tasks", "5"])
    assert (a.threads, a.max_tasks) == (3, 5)


def test_prog_name_follows_the_launcher(monkeypatch):
    monkeypatch.delenv("MITOM_LAUNCHER", raising=False)
    assert bm.build_parser().prog == "batch_mitoz.py"
    monkeypatch.setenv("MITOM_LAUNCHER", "1")
    assert bm.build_parser().prog == "mitom"
