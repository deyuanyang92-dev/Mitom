"""batch_mitoz.run_curate: the options passed to mito_curate.py depend on the clade."""
import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import batch_mitoz as bm  # noqa: E402


def test_aggressive_start_genes_only_for_annelids():
    assert bm.aggressive_start_genes("Annelida-segmented-worms") == "ND5,ND4L"
    assert bm.aggressive_start_genes("Mollusca") == "" and bm.aggressive_start_genes(None) == ""


def test_run_curate_passes_the_clade_option(tmp_path, monkeypatch):
    seen = []

    def fake_run(cmd, **kw):
        seen.append(cmd)
        (tmp_path / "05.curated").mkdir(exist_ok=True)
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(bm.subprocess, "run", fake_run)
    args = argparse.Namespace(curate="yes", genetic_code=5, clade="Annelida-segmented-worms")
    out = bm.run_curate(args, tmp_path / "04.final_gb", tmp_path)
    assert out == tmp_path / "05.curated"
    i = seen[0].index("--aggressive_start_genes")
    assert seen[0][i + 1] == "ND5,ND4L"
    args.clade = "Nematoda"
    bm.run_curate(args, tmp_path / "04.final_gb", tmp_path)
    assert seen[1][seen[1].index("--aggressive_start_genes") + 1] == ""


def test_native_core_is_the_default():
    """The native core (mito_core.py, one conda environment) is the default; MitoZ is only needed for --core hybrid."""
    import subprocess
    out = subprocess.run([sys.executable, str(HERE.parent / "batch_mitoz.py"), "run", "-h"], capture_output=True, text=True).stdout
    assert "--core {native,hybrid,mitoz}" in out.replace("\n", "").replace(" ", "") or "{native,hybrid,mitoz}" in out
    import mito_env
    assert "mitoz" not in mito_env.REQUIRED and "blastx" in mito_env.REQUIRED and "java" in mito_env.REQUIRED


def test_abbreviated_overlap_only_for_chordata():
    assert bm.abbreviated_overlap("Chordata") == ["--abbreviated_overlap"]
    assert bm.abbreviated_overlap("Annelida-segmented-worms") == [] and bm.abbreviated_overlap(None) == []


def test_no_upstream_atg_for_nematodes_and_arthropods():
    assert bm.no_upstream_atg("Nematoda") == ["--no_upstream_atg"] == bm.no_upstream_atg("Arthropoda")
    assert bm.no_upstream_atg("Mollusca") == [] and bm.no_upstream_atg(None) == []
