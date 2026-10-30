"""Debug function: stable log lines (steps, records, errors, exit), --debug, `mitom debug <out | log>`.

Run: python3 -m pytest -q tests/test_mitom_debug.py   (the real failing run takes a few seconds, no external program)
"""
import json
import logging
import os
import subprocess
import sys
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_debug  # noqa: E402
import mitom_pipeline as bm  # noqa: E402
from mitom_log import record_event, Steps  # noqa: E402

MITOM = HERE.parent / "mitom"


def _drop_handlers():
    for h in list(logging.getLogger().handlers):
        logging.getLogger().removeHandler(h)


# ── the log lines `mitom debug` reads ────────────────────────────────────────────────────────────────────────────
def test_steps_parts_failures_and_records_have_fixed_log_lines(caplog):
    class Res:
        summary, ok = "skipped: exonerate not found", False

    with caplog.at_level(logging.INFO):
        st = Steps(2)
        st.start("add_missing", "Add what the core missed")
        st.part("introns", lambda: Res())
        st.done("1 rRNA added", ok=True)
        st.start("curate", "Boundary rules")
        st.failed(RuntimeError("boom"))
        record_event("NC_1.1", "annotation failed (round-1 annotation, exit 1): boom")
    msgs = [r.getMessage() for r in caplog.records]
    assert "[STEP 1/2 add_missing] start: Add what the core missed" in msgs
    assert any(m.startswith("[STEP 1/2 add_missing] introns failed ") and m.endswith(": skipped: exonerate not found")
               for m in msgs)
    assert any(m.startswith("[STEP 1/2 add_missing] done ") and m.endswith("s: 1 rRNA added") for m in msgs)
    assert any(m.startswith("[STEP 2/2 curate] failed ") and m.endswith("s: RuntimeError: boom") for m in msgs)
    assert ("[RECORD NC_1.1] annotation failed (round-1 annotation, exit 1): boom", logging.WARNING) in \
        [(r.getMessage(), r.levelno) for r in caplog.records]


def test_an_error_goes_into_the_log_with_its_traceback_and_the_exit_code(tmp_path, monkeypatch, capsys):
    """An exception inside a workflow: the log gets [ERROR] with the traceback and '[DONE] exit=1', the console the
    error and `mitom debug OUT`; main returns 1 (before, the traceback only went to the terminal)."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    gb = tmp_path / "in.gb"
    SeqIO.write([SeqRecord(Seq("ACGT" * 10), id="X1.1", annotations={"molecule_type": "DNA"})], gb, "genbank")

    def broken(*a, **k):
        raise ValueError("broken input")
    monkeypatch.setattr(bm, "dedup_gb_records", broken)
    monkeypatch.setattr(sys, "argv", ["mitom_pipeline.py", "dedup", "-i", str(gb), "-o", str(tmp_path / "out")])
    try:
        assert bm.main() == 1
    finally:
        _drop_handlers()
    log = (tmp_path / "out" / "dedup.log").read_text()
    assert "[ERROR] ValueError: broken input" in log and "Traceback (most recent call last)" in log
    assert "[DONE] exit=1" in log and "[ENV] command: " in log          # every workflow's log starts with [ENV]
    assert f"mitom debug {tmp_path / 'out'}" in capsys.readouterr().err


# ── mitom debug on logs written by hand ──────────────────────────────────────────────────────────────────────────
def write_run(out: Path, lines: str) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    (out / "run.log").write_text(textwrap.dedent(lines).lstrip())
    return out


FAILED = """\
2026-10-30 10:00:00,001 [INFO] [ENV] MitoM 1.1.0  /opt/mitom, git abc1234
2026-10-30 10:00:00,002 [INFO] [ENV] command: mitom run -i in.fa -o out --clade Mollusca
2026-10-30 10:00:00,003 [INFO] [PARAMS]
{"cmd": "run", "input": "in.fa"}
2026-10-30 10:00:00,004 [WARNING] [ENV] MitoM needs table2asn for the NCBI submission files (install.sh; mitom doctor)
2026-10-30 10:00:00,005 [INFO] [STEP 1/11 input] start: Read input and remove duplicate records
2026-10-30 10:00:00,006 [INFO] [STEP 1/11 input] done 0.1s: 2 FASTA sequences to annotate
2026-10-30 10:00:00,007 [INFO] [STEP 2/11 annotate] start: Round-1 annotation
2026-10-30 10:00:02,000 [DEBUG] [CMD rc=2 1.5s] /env/bin/blastx -query g.fa -db cox1.fas  (in /work/S00001)
    | BLAST Database error: No alias or index file found
2026-10-30 10:00:02,001 [WARNING] [RECORD rec1] annotation failed (round-1 annotation, exit 1): RuntimeError: blastx failed (2)
2026-10-30 10:00:03,000 [INFO] [STEP 2/11 annotate] done 3.0s: 1 annotated, 1 failed (rec1)
2026-10-30 10:00:04,000 [WARNING] [RECORD rec2] NCBI ERROR: 2 (SEQ_FEAT.InternalStop) - ncbi_submission/validation_messages.tsv
2026-10-30 10:00:04,500 [WARNING] [INTRON_FIX] skipped: exonerate / blastp not found
2026-10-30 10:00:04,600 [WARNING] [INTRON_FIX] skipped: exonerate / blastp not found
2026-10-30 10:00:05,000 [WARNING] [FAIL] every record failed in the annotation core (2): see reports/failed_logs/
2026-10-30 10:00:05,001 [INFO] [DONE] exit=1 5 s
"""


def test_report_of_a_failed_run(tmp_path, capsys):
    out = write_run(tmp_path / "out", FAILED)
    logs = out / "reports" / "failed_logs"
    logs.mkdir(parents=True)
    (logs / "rec1.01.anno_round1.core.stderr.log").write_text(
        "2026-10-30 10:00:01,000 INFO [core] cox1: no blastx hit\nTraceback (most recent call last):\n"
        '  File "mito_core.py", line 1, in run\n    raise RuntimeError\n'
        "RuntimeError: blastx failed (2): BLAST Database error: No alias or index file found\n")
    assert mito_debug.main([str(out)]) == 1
    rep = capsys.readouterr().out
    assert "command : mitom run -i in.fa -o out --clade Mollusca" in rep and "MitoM   : 1.1.0" in rep
    assert "status  : FAILED (exit 1): every record failed in the annotation core" in rep
    assert " 2/11 annotate          3.0 s  done    1 annotated, 1 failed (rec1)" in rep
    assert ("first error in reports/failed_logs/rec1.01.anno_round1.core.stderr.log: RuntimeError: blastx failed (2): "
            "BLAST Database error") in rep
    assert "[exit 2, 1.5 s] /env/bin/blastx -query g.fa -db cox1.fas" in rep and "in /work/S00001" in rep
    assert "| BLAST Database error: No alias or index file found" in rep
    assert "table2asn (needed by the NCBI submission files): " in rep and "table2asn" in rep.split("Missing programs")[1]
    assert "NCBI ERRORs (1 records" in rep and "rec2: NCBI ERROR: 2 (SEQ_FEAT.InternalStop)" in rep
    assert "   2x  [INTRON_FIX] skipped: exonerate / blastp not found" in rep
    assert "the run in debug mode : mitom run -i in.fa -o " + str(out) + "_debug --clade Mollusca --debug" in rep
    assert "cd /work/S00001 &&\n      /env/bin/blastx -query g.fa -db cox1.fas" in rep


def test_report_of_a_run_stopped_by_an_error_killed_or_clean(tmp_path, capsys):
    stopped = write_run(tmp_path / "stopped", """\
        2026-10-30 10:00:00,005 [INFO] [STEP 6/11 add_missing] start: Add what the core missed
        2026-10-30 10:00:01,000 [ERROR] [STEP 6/11 add_missing] failed 1.0s: ValueError: bad location
        2026-10-30 10:00:01,001 [ERROR] [ERROR] ValueError: bad location
        Traceback (most recent call last):
          File "mitom_pipeline.py", line 9, in main
        ValueError: bad location
        2026-10-30 10:00:01,002 [INFO] [DONE] exit=1 1 s
        """)
    assert mito_debug.main([str(stopped / "run.log")]) == 1
    rep = capsys.readouterr().out
    assert "STOPPED ON AN ERROR: ValueError: bad location" in rep and "File \"mitom_pipeline.py\", line 9" in rep
    assert " 6/11 add_missing       1.0 s  failed" in rep
    killed = write_run(tmp_path / "killed", """\
        2026-10-30 10:00:00,005 [INFO] [STEP 2/11 annotate] start: Round-1 annotation
        2026-10-30 10:00:09,000 [INFO] [PROGRESS] 3/10 'rec3' done (rc=0, 9s)
        """)
    os.utime(killed / "run.log", (0, 0))
    assert mito_debug.main([str(killed)]) == 1
    assert "DID NOT FINISH in step 2/11 annotate (Round-1 annotation): the log ends at 2026-10-30 10:00:09" in \
        capsys.readouterr().out
    clean = write_run(tmp_path / "clean", """\
        2026-10-30 10:00:00,005 [INFO] [STEP 1/1 input] start: Read input
        2026-10-30 10:00:00,006 [INFO] [STEP 1/1 input] done 0.0s: 1 record
        2026-10-30 10:00:00,007 [INFO] [DONE] exit=0 0 s
        """)
    assert mito_debug.main([str(clean)]) == 0
    assert "status  : finished (exit 0)" in capsys.readouterr().out
    old = write_run(tmp_path / "old", """\
        2026-10-28 10:00:00,005 [INFO] [Step 1/11] Read input and remove duplicate records
        2026-10-28 10:02:15,006 [INFO] -> 1 records to annotate  (2 min 15 s)
        2026-10-28 10:02:15,007 [INFO] [DONE]
        """)                                                    # a log of MitoM 1.0.0
    assert mito_debug.main([str(old)]) == 0
    assert " 1/11 -               135.0 s  done    1 records to annotate" in capsys.readouterr().out
    old_failed = write_run(tmp_path / "old_failed", """\
        2026-10-28 10:00:00,007 [INFO] [DONE]
        2026-10-28 10:00:00,008 [WARNING] [FAIL] every record failed in the annotation core (1): see reports/failed_logs/
        """)                                                    # 1.0.0 wrote this after [DONE]
    assert mito_debug.main([str(old_failed)]) == 1
    assert "FAILED (exit 1): every record failed" in capsys.readouterr().out
    assert mito_debug.main([str(tmp_path / "nothing_here")]) == 2


# ── a real failing run, end to end through the `mitom` command ─────────────────────────────────────────────────
def failing_setup(tmp_path: Path) -> dict:
    """Programs that do nothing, a blastx that fails as BLAST does without its database, empty reference data."""
    b = tmp_path / "bin"
    b.mkdir()
    for prog in ("exonerate", "cmsearch", "cmalign", "nhmmer", "java", "blastp", "makeblastdb", "tblastn", "blastn",
                 "table2asn"):
        (b / prog).write_text("#!/bin/sh\nexit 0\n")
    (b / "blastx").write_text("#!/bin/sh\necho 'BLAST Database error: No alias or index file found' >&2\nexit 2\n")
    for p in b.iterdir():
        p.chmod(0o755)
    for d in ("featureProt", "ncRNA"):
        (tmp_path / "refs" / "refseq89m" / d).mkdir(parents=True)
    (tmp_path / "cfg.json").write_text(json.dumps({"tools_bin": str(b), "mitos_refdir": str(tmp_path / "refs")}))
    (tmp_path / "in.fa").write_text(">rec1 topology=circular\n" + "ACGTTGCA" * 400 + "\n")
    return dict(os.environ, YZZ_MITO_CONFIG=str(tmp_path / "cfg.json"), YZZ_MITO_NO_SCAN="1")


def mitom(*args, env, cwd):
    return subprocess.run([sys.executable, str(MITOM), *args], capture_output=True, text=True, env=env, cwd=cwd,
                          timeout=300)


def test_real_failing_run_is_explained_by_mitom_debug(tmp_path):
    env = failing_setup(tmp_path)
    opts = ["--clade", "Mollusca", "--genetic_code", "5", "--refdir", str(tmp_path / "refs"), "--rescue", "no",
            "--rrna_hmm", "no", "--intron_fix", "no", "--curate", "no", "--plot", "no", "--jobs", "1", "--threads", "1"]
    run = mitom("run", "-i", "in.fa", "-o", "out", *opts, env=env, cwd=tmp_path)
    assert run.returncode == 1, run.stderr[-2000:]
    log = (tmp_path / "out" / "run.log").read_text()
    assert "[ENV] command: mitom run -i in.fa -o out" in log and "[DONE] exit=1" in log
    assert "[RECORD rec1] annotation failed (round-1 annotation, exit 1): RuntimeError: blastx failed (2)" in log
    dbg = mitom("debug", "out", env=env, cwd=tmp_path)
    assert dbg.returncode == 1
    rep = dbg.stdout
    assert "status  : FAILED (exit 1): every record failed in the annotation core" in rep
    assert "first error in reports/failed_logs/rec1.01.anno_round1.core.stderr.log: RuntimeError: blastx failed (2): " \
           "BLAST Database error: No alias or index file found" in rep
    assert "mito_core.py annotate" in rep.split("Failed programs")[1]
    assert "the run in debug mode : mitom run -i in.fa -o out_debug" in rep and "--debug" in rep
    assert "(its folder was removed after the run: run with --debug to keep it)" in rep

    run = mitom("run", "-i", "in.fa", "-o", "dbg", *opts, "--debug", env=env, cwd=tmp_path)
    assert run.returncode == 1
    log = (tmp_path / "dbg" / "run.log").read_text()
    assert "[DEBUG]" in log and "[ENV] program blastx: " in log and "[ENV] package Bio " in log
    assert '"keep_intermediate": "yes"' in log and (tmp_path / "dbg" / "work" / "01.anno_round1").is_dir()
    rep = mitom("debug", "dbg", env=env, cwd=tmp_path).stdout
    assert "failed inside it: [CMD rc=2 " in rep and "/bin/blastx -query " in rep    # the core's own failed program
    assert "first error in work/01.anno_round1/S00001/logs/core.stderr.log: RuntimeError: blastx failed (2)" in rep
    assert "record rec1 alone   : mitom run -i dbg/work/00.internal_fasta/S00001.fasta -o dbg_rec1" in rep
    hand = rep.split("the first failed program by hand:")[1].split("record rec1")[0]
    assert "/dbg/work/01.anno_round1/S00001 &&" in hand
    again = subprocess.run(hand.strip().replace("&&\n", "&&"), shell=True, capture_output=True, text=True, env=env)
    assert again.returncode != 0 and "blastx failed (2)" in again.stderr        # the failure, reproduced by hand


def test_debug_flag_works_with_every_command_and_steps_log_their_programs(tmp_path):
    """`mitom <command> --debug`: the launcher sets MITOM_DEBUG=1 for everything the command starts and drops the flag
    for the scripts that do not know it; a step script in debug mode logs each program it runs on stderr (the workflow
    keeps that in its log)."""
    env = dict(os.environ)
    env.pop("MITOM_DEBUG", None)
    assert subprocess.run([sys.executable, str(MITOM), "curate", "-h", "--debug"], capture_output=True, text=True,
                          env=env, timeout=120).returncode == 0
    p = subprocess.run([sys.executable, str(MITOM), "debug", "-h"], capture_output=True, text=True, env=env, timeout=120)
    assert p.returncode == 0 and "reproduce" in p.stdout
    code = "import mito_common; mito_common.run_tool(['sh', '-c', 'echo note >&2'])"
    q = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=HERE.parent,
                       env=dict(env, MITOM_DEBUG="1"), timeout=120)
    assert "[DEBUG] [CMD rc=0 " in q.stderr and "sh -c 'echo note >&2'" in q.stderr and "    | note" in q.stderr
    q = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=HERE.parent, env=env, timeout=120)
    assert "[CMD" not in q.stderr                               # not in debug mode: nothing on stderr
