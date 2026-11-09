"""mitom_output.organize_outputs: the final layout of `run`, the user's own files and a second run into the same folder."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mitom_output  # noqa: E402


def _fake_run(out: Path, errors: str) -> Path:
    (out / "05.curated").mkdir(exist_ok=True)
    (out / "05.curated" / "A.gbf").write_text("LOCUS A\n//\n")
    (out / "05.curated" / "curation_report.tsv").write_text("record\trule\n")
    (out / "06.ncbi_submission").mkdir(exist_ok=True)
    (out / "06.ncbi_submission" / "submission_summary.tsv").write_text(f"record\tERROR\nA\t{errors}\n")
    (out / "rescue_report.tsv").write_text("file\trecord\n")
    return out / "05.curated"


def test_organize_touches_only_mitom_files_and_can_run_twice(tmp_path):
    out = tmp_path
    (out / "my_input.gb").write_text("LOCUS X\n//\n")             # the user's own file in the output folder
    (out / "01.raw").mkdir()
    (out / "01.raw" / "keep.txt").write_text("mine")
    mitom_output.organize_outputs(out, _fake_run(out, "1"))
    assert (out / "final_annotation" / "A.gb").is_file() and not (out / "05.curated").exists()
    assert (out / "reports" / "curation_report.tsv").is_file() and (out / "reports" / "rescue_report.tsv").is_file()
    assert (out / "my_input.gb").is_file() and (out / "01.raw" / "keep.txt").read_text() == "mine"
    assert (out / "ncbi_submission" / "submission_summary.tsv").read_text().endswith("A\t1\n")
    assert (out / "README.txt").is_file()
    mitom_output.organize_outputs(out, _fake_run(out, "0"))                 # second run into the same folder
    assert (out / "ncbi_submission" / "submission_summary.tsv").read_text().endswith("A\t0\n")
    assert not (out / "ncbi_submission" / "06.ncbi_submission").exists()
    assert sorted(p.name for p in (out / "final_annotation").iterdir()) == ["A.gb"]
    assert not (out / "reports" / "reports").exists()


def test_organize_keep_moves_work_dirs(tmp_path):
    out = tmp_path
    final = _fake_run(out, "0")
    (out / "01.anno_round1").mkdir()
    (out / "01.anno_round1" / "x.log").write_text("")
    mitom_output.organize_outputs(out, final, keep=True)
    assert (out / "work" / "01.anno_round1" / "x.log").is_file() and not (out / "01.anno_round1").exists()
    assert (out / "work" / "05.curated").is_dir()


def test_records_without_the_orientation_gene_are_listed_in_reports(tmp_path):
    """2026-10-30 (integration audit): no_<gene>_samples.txt (records without any orientation gene) stayed at the top
    of the output folder; it is a report of the run like the others."""
    out = tmp_path
    (out / "no_cox1_samples.txt").write_text("# samples where 'cox1' was not found\nS1\n")
    mitom_output.organize_outputs(out, _fake_run(out, "0"), reports_named=("no_cox1_samples.txt",))
    assert (out / "reports" / "no_cox1_samples.txt").is_file() and not (out / "no_cox1_samples.txt").exists()


def _real_record(path: Path) -> None:
    from Bio import SeqIO
    rec = SeqIO.read(str(HERE / "data" / "NC_000931.gb"), "genbank")
    SeqIO.write([rec], str(path), "genbank")


def test_grouped_layout(tmp_path):
    """Issue #9: --output_layout grouped (default) = summary.tsv, README.txt and annotation/ (GenBank, genes.fna,
    proteins.faa, plots/) submission/<record>/ qc/ input/ logs/ at the top; the user's own files stay; a second run
    replaces the first."""
    import argparse
    out = tmp_path
    (out / "my_notes.txt").write_text("mine")
    for run in (1, 2):
        final = _fake_run(out, "0")
        _real_record(final / "A.gbf")
        for ext in ("fsa", "tbl", "sqn", "val"):
            (out / "06.ncbi_submission" / f"A.{ext}").write_text("x")
        (out / "original_annotation").mkdir(exist_ok=True)
        (out / "original_annotation" / "A.gb").write_text("LOCUS A\n//\n")
        (out / "plots").mkdir(exist_ok=True)
        (out / "plots" / "A.png").write_text("png")
        (out / "original_vs_mitom.tsv").write_text("record\n")
        (out / "all.final.gbf").write_text("LOCUS A\n//\n")
        (out / "reports").mkdir(exist_ok=True)
        (out / "reports" / "run_params.json").write_text("{}")
        (out / "reports" / "intron_fix.log").write_text("log")
        (out / "run.log").write_text(f"[ENV] MitoM run {run}\n")
        (out / "annotation-summary.tsv").write_text("sample_name\tPCG\nA\t13\n")
        mitom_output.organize_outputs(out, final)
        mitom_output.regroup_outputs(out)
    top = sorted(p.name for p in out.iterdir())
    assert top == ["README.txt", "annotation", "input", "logs", "my_notes.txt", "qc", "submission", "summary.tsv"]
    ann = out / "annotation"
    assert sorted(p.name for p in ann.iterdir()) == ["A.gb", "all_records.gb", "genes.fna", "plots", "proteins.faa"]
    genes = (ann / "genes.fna").read_text().split(">")[1:]
    assert len(genes) == 15 and genes[0].startswith("NC_000931.1|") and "|Platynereis_dumerilii\n" in genes[0]
    prots = (ann / "proteins.faa").read_text()
    assert prots.count(">") == 13 and "\nM" in prots
    assert sorted(p.name for p in (out / "submission" / "A").iterdir()) == ["A.fsa", "A.sqn", "A.tbl", "A.val"]
    assert (out / "submission" / "submission_summary.tsv").is_file()
    assert {"curation_report.tsv", "rescue_report.tsv", "original_vs_mitom.tsv"} <= {p.name for p in (out / "qc").iterdir()}
    assert sorted(p.name for p in (out / "logs").iterdir()) == ["intron_fix.log", "run.log", "run_params.json"]
    assert (out / "logs" / "run.log").read_text() == "[ENV] MitoM run 2\n" and (out / "input" / "A.gb").is_file()
    readme = (out / "README.txt").read_text()
    assert "annotation/" in readme and "qc/" in readme and "final_annotation" not in readme
    assert mitom_output._mitom_output(out)                       # a third run knows the folder is MitoM's
    removed = mitom_output.clear_previous_work(out)
    assert {"annotation/", "qc/", "logs/", "submission/", "input/"} <= set(removed) and (out / "my_notes.txt").is_file()
    a = argparse.Namespace(output_layout="grouped", cmd="run")
    assert mitom_output.where(a, "reports/curation_report.tsv") == "qc/curation_report.tsv"
    assert mitom_output.where(a, "reports/failed_logs/") == "logs/failed_logs/"
    assert mitom_output.where(argparse.Namespace(output_layout="classic"), "reports/x.tsv") == "reports/x.tsv"
