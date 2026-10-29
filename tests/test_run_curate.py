"""mitom_pipeline.run_curate: the options passed to mito_curate.py depend on the clade."""
import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mitom_pipeline as bm  # noqa: E402


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
    assert out.gb_dir == tmp_path / "05.curated" and out.ok
    i = seen[0].index("--aggressive_start_genes")
    assert seen[0][i + 1] == "ND5,ND4L"
    args.clade = "Nematoda"
    bm.run_curate(args, tmp_path / "04.final_gb", tmp_path)
    assert seen[1][seen[1].index("--aggressive_start_genes") + 1] == ""


def test_native_core_by_default_engines_optional():
    """MitoM's own core (mito_core.py) annotates by default; MitoZ / MITOS2 are optional engines (--engine), never
    needed by MitoM itself (user 2026-10-28)."""
    import subprocess
    out = subprocess.run([sys.executable, str(HERE.parent / "mitom_pipeline.py"), "run", "-h"], capture_output=True, text=True).stdout
    assert "--core" not in out and "--engine {native,mitoz,mitos2}" in out and "--plot" in out
    assert bm.build_parser().parse_args(["run", "-i", "x", "-o", "y"]).engine == "native"
    import mito_env
    assert "mitoz" not in mito_env.TOOLS and "blastx" in mito_env.REQUIRED and "java" in mito_env.REQUIRED


def test_abbreviated_overlap_only_for_chordata():
    assert bm.abbreviated_overlap("Chordata") == ["--abbreviated_overlap"]
    assert bm.abbreviated_overlap("Annelida-segmented-worms") == [] and bm.abbreviated_overlap(None) == []


def test_no_upstream_atg_for_nematodes_and_arthropods():
    assert bm.no_upstream_atg("Nematoda") == ["--no_upstream_atg"] == bm.no_upstream_atg("Arthropoda")
    assert bm.no_upstream_atg("Mollusca") == [] and bm.no_upstream_atg(None) == []


def test_intron_step_runs_de_novo_then_site_transfer(tmp_path, monkeypatch):
    """2026-10-29: mito_intron.py (de novo) is followed by mito_intron_ref.py (intron sites of annotated relatives),
    in place; both reports go to the output root; --phylum is passed as the lineage of records without a taxonomy."""
    seen = []
    final = tmp_path / "04.final_gb"
    final.mkdir()
    bindir = tmp_path / "bin"
    bindir.mkdir()
    for name in ("exonerate", "blastp"):
        (bindir / name).write_text("")

    def fake_run(cmd, **kw):
        seen.append(cmd)
        if Path(cmd[1]).name == "mito_intron.py":
            (final / "intron_report.tsv").write_text("record\n")
            return subprocess.CompletedProcess(cmd, 0, "[mito_intron] 1 file(s); 2 intron(s) written -> x\n", "")
        (final / "intron_ref_report.tsv").write_text("record\n")
        return subprocess.CompletedProcess(cmd, 0, "", "[mito_intron_ref] 1 file(s), library 297 genes + 0 own; "
                                                        "1 gene model(s) written, 0 low-confidence not written -> y\n")

    monkeypatch.setattr(bm.subprocess, "run", fake_run)
    args = argparse.Namespace(intron_fix="auto", intron_bin=str(bindir), genetic_code=4, threads=2, max_tasks=3,
                              intron_lib=[], intron_ref_gb=[], phylum="Cnidaria")
    line = bm.run_intron_fix(args, final, tmp_path).summary
    assert [Path(c[1]).name for c in seen] == ["mito_intron.py", "mito_intron_ref.py"]
    ref = seen[1]
    assert ref[2] == "annotate" and ref[ref.index("-i") + 1] == ref[ref.index("-o") + 1] == str(final)
    assert ref[ref.index("--genetic_code") + 1] == "4" and ref[ref.index("--bin") + 1] == str(bindir)
    assert ref[ref.index("--lineage") + 1] == "Eukaryota;Metazoa;Cnidaria"
    assert (tmp_path / "reports" / "intron_report.tsv").exists() and (tmp_path / "reports" / "intron_ref_report.tsv").exists()
    assert "site transfer: 1 file(s)" in line and "1 gene model(s) written" in line
    args.phylum = None
    seen.clear()
    bm.run_intron_fix(args, final, tmp_path)
    assert "--lineage" not in seen[1]
