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
