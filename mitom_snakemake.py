"""mitom_snakemake.py - the Snakemake backend of the annotation stage (--scheduler_backend snakemake / auto): one rule
per sample that runs the sample's engine command through a self-contained worker script. Snakemake is optional (not
installed with MitoM); everything after the annotation runs in the main process."""
from __future__ import annotations

import argparse
import json
import logging
import math
import re
import shlex
import shutil
import subprocess
import sys
import threading
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from mitom_log import Progress, tail_text, watch_files


def which_or_none(exe: str) -> Optional[str]:
    try:
        return shutil.which(exe)
    except Exception:
        return None


def int_mem_mb(x, minimum: int = 128) -> int:
    try:
        v = int(math.ceil(float(x)))
    except Exception:
        v = minimum
    return max(minimum, v)


def snakemake_available(snakemake_bin: Optional[str]) -> bool:
    if snakemake_bin:
        return Path(snakemake_bin).exists() or which_or_none(snakemake_bin) is not None
    return which_or_none("snakemake") is not None


@dataclass
class AnnotTask:
    internal_id:   str
    sample_name:   str
    fasta_file:    str
    phase_dirname: str
    out_root:      str
    lineage:       str = ""
    cmd:           List[str] = field(default_factory=list)      # engine_cmd(): what the worker runs


def render_snakemake_worker_script(out_dir: Path) -> Path:
    """Generate the self-contained annotation worker script (runs the task's engine command, engine_cmd()). The
    worker annotates only; finalize_genbank_output and every later step run in the main process."""
    fp   = out_dir / "_worker.py"
    code = r'''#!/usr/bin/env python3
"""Self-contained annotation worker invoked by Snakemake: runs the task's command (the annotation engine).

This script handles annotation only.  Finalization steps
(ID-rewrite, reorientation, metadata overlay, organism sync) are
executed in the main pipeline process after all workers finish.
"""
import json, subprocess, sys, time
from pathlib import Path

def tail_text(s, n=120): return "\n".join((s or "").splitlines()[-n:])

def main():
    task = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    internal_id   = task["internal_id"]
    phase_dirname = task["phase_dirname"]
    out_root      = Path(task["out_root"]).resolve()
    fasta_file    = Path(task["fasta_file"]).resolve()
    workdir       = out_root / phase_dirname / internal_id
    workdir.mkdir(parents=True, exist_ok=True)
    logs_dir = workdir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    cmd = task["cmd"]                               # engine_cmd() of the main process

    t0 = time.time()
    p  = subprocess.run(cmd, cwd=str(workdir),
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    dt = time.time() - t0

    if p.stdout:
        (logs_dir / "core.stdout.log").write_text(
            p.stdout, encoding="utf-8", errors="ignore")
    if p.stderr:
        (logs_dir / "core.stderr.log").write_text(
            p.stderr, encoding="utf-8", errors="ignore")
    if p.returncode != 0:
        (logs_dir / "core.returncode.txt").write_text(
            str(p.returncode) + "\n", encoding="utf-8")
    (logs_dir / "core.runtime_sec.txt").write_text(
        f"{dt:.2f}\n", encoding="utf-8")

    result = {
        "internal_id": internal_id,
        "sample_name": task.get("sample_name", internal_id),
        "ok": False, "workdir": str(workdir), "best_gb": None,
        "returncode": int(p.returncode),
    }
    if p.returncode == 0:
        own = workdir / (internal_id + ".gbf")          # what the engine writes (exit 0 = written)
        result.update({"ok": own.is_file(), "best_gb": str(own) if own.is_file() else None})
    else:
        result["stderr_tail"] = tail_text(p.stderr or "")

    Path(sys.argv[1]).with_suffix(".result.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    Path(sys.argv[1]).with_suffix(".done").write_text("done\n", encoding="utf-8")

if __name__ == "__main__":
    main()
'''
    fp.write_text(code, encoding="utf-8")
    return fp


def render_snakemake_snakefile(
    workdir: Path, task_jsons: List[Path], worker_script: Path,
    python_exe: str = sys.executable or "python3",
    job_threads: int = 1,
    job_mem_mb: int = 128,
) -> Path:
    snakefile = workdir / "Snakefile"
    rule_lines: List[str] = []
    all_done:   List[str] = []
    for tp in task_jsons:
        done_s  = str(tp.with_suffix(".done"))
        rule_id = re.sub(r"[^A-Za-z0-9_]", "_", tp.stem)
        all_done.append(done_s)
        cmdline = " ".join(shlex.quote(x) for x in [python_exe, str(worker_script), str(tp)])
        # threads / mem_mb per job: what the engine gets (--threads), so that Snakemake schedules the real load
        rule_lines.append(f"""
rule task_{rule_id}:
    input:  task={json.dumps(str(tp), ensure_ascii=False)}
    output: done={json.dumps(done_s, ensure_ascii=False)}
    threads: {max(1, int(job_threads))}
    resources: mem_mb={int_mem_mb(job_mem_mb)}
    shell:  {json.dumps(cmdline, ensure_ascii=False)}
""")
    header = f"rule all:\n    input: {json.dumps(all_done, ensure_ascii=False)}\n"
    snakefile.write_text(header + "\n".join(rule_lines), encoding="utf-8")
    return snakefile


def run_snakemake_for_annotation_tasks(
    *, tasks: List[AnnotTask], out_root: Path, phase_name: str, args: argparse.Namespace,
    progress: Optional[Progress] = None,
) -> Dict[str, Dict]:
    """Run the tasks with Snakemake; `progress` (the stage's) is updated from the result files every 10 s."""
    snk_work = out_root / f"_snakemake_{phase_name}"
    snk_work.mkdir(parents=True, exist_ok=True)
    worker_script = render_snakemake_worker_script(snk_work)
    task_jsons: List[Path] = []
    for t in tasks:
        tp = snk_work / f"{t.internal_id}.json"
        tp.write_text(json.dumps(asdict(t), indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
        task_jsons.append(tp)
    per_job_mem_mb = max(128, int_mem_mb(args.smk_mem_mb) // max(1, int(args.max_tasks)))
    snakefile = render_snakemake_snakefile(
        snk_work, task_jsons, worker_script,
        python_exe=sys.executable or "python3",
        job_threads=int(args.threads),
        job_mem_mb=per_job_mem_mb)
    smk_bin   = args.snakemake_bin or which_or_none("snakemake")
    if not smk_bin:
        raise FileNotFoundError("snakemake executable not found")
    cmd = [
        str(smk_bin), "--snakefile", str(snakefile), "--directory", str(snk_work),
        "--cores", str(max(1, int(args.smk_cores))),
        "--jobs",  str(max(1, int(args.max_tasks))),
        "--resources", f"mem_mb={int_mem_mb(args.smk_mem_mb)}",
        "--rerun-incomplete", "--keep-going", "--printshellcmds",
        "--latency-wait", str(int(args.smk_latency_wait)),
    ]
    smk_extra: List[str] = []
    for token in (getattr(args, "smk_extra", None) or []):
        smk_extra.extend(shlex.split(str(token)))
    if smk_extra:
        logging.info("[SNAKEMAKE_EXTRA] " + " ".join(smk_extra))
        cmd.extend(smk_extra)

    def _run(c: List[str]) -> "subprocess.CompletedProcess[str]":
        logging.info("[SNAKEMAKE] " + " ".join(c))
        return subprocess.run(c, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    stop = threading.Event()
    if progress is not None:
        watch_files(snk_work, "*.result.json", progress, stop)
    try:
        p = _run(cmd)
        if p.returncode != 0:
            combined = (p.stdout or "") + (p.stderr or "")
            if "LockException" in combined or "Directory cannot be locked" in combined:
                logging.warning(f"[SNAKEMAKE_LOCK] Stale lock in {snk_work} — attempting --unlock…")
                pu = subprocess.run(
                    [str(smk_bin), "--snakefile", str(snakefile),
                     "--directory", str(snk_work), "--unlock"],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if pu.returncode == 0:
                    logging.warning("[SNAKEMAKE_LOCK] Unlock OK — retrying…")
                    p = _run(cmd)
                else:
                    logging.error(f"[SNAKEMAKE_UNLOCK_FAIL] rc={pu.returncode}")
    finally:
        stop.set()
    if progress is not None:                                   # the watcher polls every 10 s: print the final count
        progress.update(sum(1 for tp in task_jsons if tp.with_suffix(".result.json").exists()))

    (snk_work / "snakemake.stdout.log").write_text(p.stdout or "", encoding="utf-8",
                                                    errors="ignore")
    (snk_work / "snakemake.stderr.log").write_text(p.stderr or "", encoding="utf-8",
                                                    errors="ignore")
    if p.returncode != 0:
        raise RuntimeError(
            f"Snakemake failed for phase '{phase_name}'\n"
            f"stdout:\n{tail_text(p.stdout or '')}\n\nstderr:\n{tail_text(p.stderr or '')}")

    results: Dict[str, Dict] = {}
    for tp in task_jsons:
        rp = tp.with_suffix(".result.json")
        if rp.exists():
            obj = json.loads(rp.read_text(encoding="utf-8"))
            results[obj["internal_id"]] = obj
            if not obj.get("ok"):                              # the worker's sample.log holds the details
                last = next((x.strip() for x in reversed(str(obj.get("stderr_tail") or "").splitlines())
                             if x.strip()), "no stderr")
                logging.warning(f"[FAIL] {obj.get('sample_name', obj['internal_id'])}: annotation core exit "
                                f"{obj.get('returncode', '?')}: {last}  (sample log: {obj.get('workdir')}/sample.log)")
        else:
            iid = tp.stem
            logging.warning(f"[SNAKEMAKE_RESULT_MISSING] {rp} not found for task {iid}")
            results[iid] = {
                "internal_id": iid, "sample_name": iid, "ok": False,
                "workdir": str(out_root / phase_name / iid), "best_gb": None,
            }
    return results
