"""Structure rules of the refactoring (2026-11-08), checked on the source so they stay true:
no module globals rewritten at run time, no sys.path edits, GenBank written only by mito_common.write_gb, programs
started only by mito_common.run_tool (except the probes and pipes listed), mitom_pipeline.py only the entry point."""
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MODULES = {p.name: p.read_text() for p in sorted(REPO.glob("*.py"))}
# subprocess used directly: run_tool itself; version / git / Python probes; mito_plot's read-mapping pipes (Popen);
# the Snakemake worker script (a template run outside MitoM's modules)
RAW_SUBPROCESS = {"mito_common.py", "mito_env.py", "mitom_prepare.py", "mitom_main.py", "mito_plot.py",
                  "mitom_snakemake.py",
                  "mito_rna.py"}        # (mito_plot / mito_rna: the mapper's SAM is read as a stream, hundreds of
                                        # millions of reads - run_tool would hold it in memory)


def code_lines(text: str):
    return [ln for ln in text.splitlines() if not ln.lstrip().startswith("#")]


def test_no_global_statement():
    bad = [f"{n}: {ln.strip()}" for n, t in MODULES.items() for ln in code_lines(t) if re.match(r"\s+global\s", ln)]
    assert not bad, "pass the value (RunContext, a parameter) instead of a module global:\n" + "\n".join(bad)


def test_no_sys_path_edit():
    assert not [n for n, t in MODULES.items() if "sys.path.insert" in t or "sys.path.append" in t]


def test_genbank_written_by_write_gb():
    pat = re.compile(r"SeqIO\.write\([^\n]*,\s*[\"']genbank[\"']\)\s*$", re.M)   # format argument genbank
    assert [n for n, t in MODULES.items() if pat.search(t)] == ["mito_common.py"]


def test_programs_started_by_run_tool():
    pat = re.compile(r"subprocess\.(run|Popen|call|check_call|check_output)\(")
    bad = sorted(n for n, t in MODULES.items() if n not in RAW_SUBPROCESS and pat.search("\n".join(code_lines(t))))
    assert not bad, f"use mito_common.run_tool in {bad}"


def test_mitom_pipeline_is_only_the_entry_point():
    tree = MODULES["mitom_pipeline.py"]
    defs = re.findall(r"^def (\w+)", tree, re.M)
    assert defs == ["main"] and len(tree.splitlines()) < 250


def test_annotate_stages_are_run_stages_without_round2_compare_output():
    """`mitom annotate` = the run stages minus re-orientation, round 2, comparison and output layout; the terminal
    counts its own steps (start_run: Steps(len(stages))) and shows the round-2 line only when round 2 is a stage."""
    import inspect
    import mitom_run as mr
    assert set(mr.ANNOTATE_STAGES) <= set(mr.RUN_STAGES)
    assert set(mr.RUN_STAGES) - set(mr.ANNOTATE_STAGES) == {mr.stage_reorient, mr.stage_round2, mr.stage_compare,
                                                            mr.stage_output}
    src = inspect.getsource(mr.start_run)
    assert "Steps(len(stages))" in src and "stage_round2 in stages" in src and "hasattr" not in src
