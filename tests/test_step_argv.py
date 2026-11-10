"""The pipeline starts mito_curate.py / mito_intron.py / mito_intron_ref.py / mito_orf.py as programs (crash isolation, captured
stderr summary lines, fork pools in the child). These tests keep the argument lists the pipeline builds in step with
the scripts' own parsers: an option renamed or retyped in a script fails here, not in a user's run (review 2026-11-08)."""
import argparse
import subprocess

import pytest

import mito_curate
import mito_intron
import mito_intron_ref
import mito_orf
import mitom_steps


def _script_args(cmd, script):
    """The arguments after the script path of a [python, (-W ignore), script, ...] command."""
    i = next(k for k, x in enumerate(cmd) if x.endswith(script))
    return cmd[i + 1:]


@pytest.mark.parametrize("clade", ["Annelida-segmented-worms", "Chordata", "Nematoda", "Mollusca", None])
def test_curate_command_parses(tmp_path, monkeypatch, clade):
    seen = []
    monkeypatch.setattr(subprocess, "run", lambda cmd, **kw: seen.append(cmd) or subprocess.CompletedProcess(cmd, 0, "", ""))
    args = argparse.Namespace(curate="yes", genetic_code=5, clade=clade, phylum=None)
    mitom_steps.run_curate(args, tmp_path / "in", tmp_path)
    a = mito_curate.build_parser().parse_args(_script_args(seen[0], "mito_curate.py"))
    assert a.genetic_code == 5 and str(a.input) == str(tmp_path / "in") and str(a.out) == str(tmp_path / "05.curated")


def test_intron_commands_parse(tmp_path):
    gb = tmp_path / "gb"
    denovo = mitom_steps.intron_denovo_cmd(gb, gb, 4, 2, 3, libs=["x.faa"], ref_gbs=["r.gb"], bin_dir="B")
    a = mito_intron.build_parser().parse_args(_script_args(denovo, "mito_intron.py"))
    assert a.genetic_code == 4 and a.threads == 2 and a.jobs == 3
    ref = mitom_steps.intron_ref_cmd(gb, gb, 4, 3, "B", "Eukaryota;Metazoa;Cnidaria")
    b = mito_intron_ref.build_parser().parse_args(_script_args(ref, "mito_intron_ref.py"))
    assert b.cmd == "annotate" and b.genetic_code == 4 and b.jobs == 3 and b.lineage == "Eukaryota;Metazoa;Cnidaria"


def test_orf_command_parses(tmp_path):
    cmd = mitom_steps.orf_cmd(tmp_path / "gb", tmp_path / "o", 4, "all", 3, "B", ["r1.gb", "r2"])
    a = mito_orf.build_parser().parse_args(_script_args(cmd, "mito_orf.py"))
    assert a.genetic_code == 4 and a.orf == "all" and a.threads == 3 and a.bin == "B" and a.ref == ["r1.gb", "r2"]


def test_orf_step_failure_leaves_the_records(tmp_path, monkeypatch):
    """mito_orf.py exits non-zero: the step reports it, the records are untouched, the run goes on."""
    gb = tmp_path / "gb"
    gb.mkdir()
    (gb / "a.gbf").write_text("LOCUS x\n//\n")
    monkeypatch.setattr(mitom_steps, "run_tool", lambda cmd, **kw: subprocess.CompletedProcess(cmd, 1, "", "boom"))
    r = mitom_steps.run_orf(argparse.Namespace(orf="yes", genetic_code=5, threads=1), gb, tmp_path)
    assert not r.ok and (gb / "a.gbf").read_text() == "LOCUS x\n//\n"


def test_pcg_hmm_command_parses(tmp_path):
    import mito_pcg_hmm
    cmd = mitom_steps.pcg_hmm_cmd(tmp_path / "gb", tmp_path / "r.tsv", 4, "B")
    a = mito_pcg_hmm.build_parser().parse_args(_script_args(cmd, "mito_pcg_hmm.py"))
    assert a.cmd == "search" and a.genetic_code == 4 and a.report == str(tmp_path / "r.tsv") and a.bin == "B"
