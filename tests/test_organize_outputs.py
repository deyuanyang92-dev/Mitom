"""mitom_pipeline.organize_outputs: the final layout of `run`, the user's own files and a second run into the same folder."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mitom_pipeline as bm  # noqa: E402


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
    bm.organize_outputs(out, _fake_run(out, "1"))
    assert (out / "final_annotation" / "A.gb").is_file() and not (out / "05.curated").exists()
    assert (out / "reports" / "curation_report.tsv").is_file() and (out / "reports" / "rescue_report.tsv").is_file()
    assert (out / "my_input.gb").is_file() and (out / "01.raw" / "keep.txt").read_text() == "mine"
    assert (out / "ncbi_submission" / "submission_summary.tsv").read_text().endswith("A\t1\n")
    assert (out / "README.txt").is_file()
    bm.organize_outputs(out, _fake_run(out, "0"))                 # second run into the same folder
    assert (out / "ncbi_submission" / "submission_summary.tsv").read_text().endswith("A\t0\n")
    assert not (out / "ncbi_submission" / "06.ncbi_submission").exists()
    assert sorted(p.name for p in (out / "final_annotation").iterdir()) == ["A.gb"]
    assert not (out / "reports" / "reports").exists()


def test_organize_keep_moves_work_dirs(tmp_path):
    out = tmp_path
    final = _fake_run(out, "0")
    (out / "01.anno_round1").mkdir()
    (out / "01.anno_round1" / "x.log").write_text("")
    bm.organize_outputs(out, final, keep=True)
    assert (out / "work" / "01.anno_round1" / "x.log").is_file() and not (out / "01.anno_round1").exists()
    assert (out / "work" / "05.curated").is_dir()


def test_records_without_the_orientation_gene_are_listed_in_reports(tmp_path):
    """2026-10-30 (integration audit): no_<gene>_samples.txt (records without any orientation gene) stayed at the top
    of the output folder; it is a report of the run like the others."""
    out = tmp_path
    (out / "no_cox1_samples.txt").write_text("# samples where 'cox1' was not found\nS1\n")
    bm.organize_outputs(out, _fake_run(out, "0"), reports_named=("no_cox1_samples.txt",))
    assert (out / "reports" / "no_cox1_samples.txt").is_file() and not (out / "no_cox1_samples.txt").exists()
