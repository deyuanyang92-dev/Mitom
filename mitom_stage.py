"""mitom_stage.py - the annotation stage: the command of each sample's engine (MitoM's core mito_core.py, or MitoZ /
MITOS2 through mito_engines.py; the input record itself for `mitom import`), one job per sample with its own log, the
jobs in threads or with Snakemake (mitom_snakemake.py), progress on the terminal."""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import functools
import logging
import os
import sys
import threading
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from mito_common import run_tool, safe_name
from mitom_inputs import SampleInfo
from mitom_log import Progress, close_sample_logger, get_sample_logger, record_event, tail_text
from mitom_snakemake import AnnotTask, run_snakemake_for_annotation_tasks, snakemake_available


_STAGE_PROGRESS: Optional[Progress] = None          # set by run_annotation_stage; updated per finished sample
STAGE_LABEL = {"01.anno_round1": "round-1 annotation", "02.anno_topology": "round-1 annotation, topology changed",
               "03.anno_round2_reoriented": "round-2 annotation"}


MITO_CORE = Path(__file__).resolve().parent / "mito_core.py"      # MitoM's annotation core
MITO_ENGINES = Path(__file__).resolve().parent / "mito_engines.py"  # MitoZ / MITOS2 as engines, imported results
ENGINE_LABEL = {"native": "MitoM core", "mitoz": "MitoZ", "mitos2": "MITOS2", "input": "annotation of the input"}


def run_cmd_capture(
    cmd: List[str], cwd: Optional[Path], stdout_path: Path, stderr_path: Path,
) -> Tuple[int, float]:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    p  = run_tool(cmd, cwd=cwd)
    dt = time.time() - t0
    if p.stdout:
        stdout_path.write_text(p.stdout, encoding="utf-8", errors="ignore")
    if p.stderr:
        stderr_path.write_text(p.stderr, encoding="utf-8", errors="ignore")
    return int(p.returncode), float(dt)


def sample_lineage(si: "SampleInfo") -> str:
    """NCBI lineage of a sample from its original GenBank record ('' for FASTA input): the MitoM-ML tRNA caller of
    mito_core.py uses the clade models on this lineage (full lineage 98.5 % vs phylum only 97.9 % tRNA recall,
    projects/invert_trna)."""
    if not si.orig_gb_path:
        return ""
    try:
        lineages = _file_lineages(str(si.orig_gb_path), os.path.getmtime(si.orig_gb_path))
    except OSError:
        return ""
    # the sample's own record of a multi-record file (2026-10-31: every record took the first record's lineage)
    for key in (si.orig_token, si.sample_name, safe_name(si.orig_token)):
        if key in lineages:
            return lineages[key]
    return lineages.get("", "")


@functools.lru_cache(maxsize=64)
def _file_lineages(path: str, mtime: float = 0.0) -> Dict[str, str]:
    """{record id / name / accession.version: lineage} of the records of a GenBank file ('' = the lineage when the
    file has one record); read once per file (and modification time)."""
    from mitom_metadata import gb_record_accession
    out: Dict[str, str] = {}
    try:
        recs = list(SeqIO.parse(path, "genbank"))
    except (OSError, ValueError):
        return out
    for rec in recs:
        tax = ";".join(t for t in rec.annotations.get("taxonomy", []) if t)
        for key in (rec.id, rec.name, safe_name(rec.id), gb_record_accession(rec)):
            out.setdefault(key, tax)
    if len(recs) == 1:
        out[""] = out.get(recs[0].id, "")
    return out


def native_core_cmd(fasta_file: Path, workdir: Path, internal_id: str, args, lineage: str = "") -> List[str]:
    """mito_core.py annotate command for one sample (--core native); `args` is the Namespace or the worker's dict."""
    g = (lambda k, d=None: getattr(args, k, d)) if not isinstance(args, dict) else (lambda k, d=None: args.get(k, d))
    cmd = [sys.executable, "-W", "ignore", str(MITO_CORE), "annotate", "-i", str(fasta_file), "-o", str(workdir),
           "--prefix", internal_id, "--genetic_code", str(int(g("genetic_code", 5))), "--threads", str(int(g("threads", 1)))]
    if g("mitos_refdir"):
        cmd += ["--refdir", str(g("mitos_refdir"))]
    if g("mitos_bin"):
        cmd += ["--mitos_bin", str(g("mitos_bin"))]
    if g("phylum"):
        cmd += ["--phylum", str(g("phylum"))]
    if g("rrna_hmm_lib"):
        cmd += ["--rrna_lib", str(g("rrna_hmm_lib"))]
    if g("rrna_method"):
        cmd += ["--rrna_method", str(g("rrna_method"))]
    if lineage:
        cmd += ["--lineage", lineage]
    if g("topology", "auto") == "auto":                  # the evidence --topology auto decides with (mitom_topology.py)
        cmd += ["--check_ends"]
    return cmd


def engine_cmd(fasta_file: Path, workdir: Path, internal_id: str, args, lineage: str = "") -> List[str]:
    """Command that writes <workdir>/<internal_id>.gbf with the annotation of --engine: MitoM's core (native), or
    MitoZ / MITOS2 run by mito_engines.py; `args` is the Namespace or a dict."""
    g = (lambda k, d=None: getattr(args, k, d)) if not isinstance(args, dict) else (lambda k, d=None: args.get(k, d))
    engine = g("engine") or "native"
    if engine == "native":
        return native_core_cmd(fasta_file, workdir, internal_id, args, lineage)
    cmd = [sys.executable, "-W", "ignore", str(MITO_ENGINES), engine, "-i", str(fasta_file), "-o", str(workdir),
           "--prefix", internal_id, "--genetic_code", str(int(g("genetic_code", 5))), "--threads", str(int(g("threads", 1)))]
    if engine == "mitoz":
        cmd += ["--clade", str(g("clade"))]
        cmd += ["--mitoz_bin", str(g("mitoz_bin"))] if g("mitoz_bin") else []
        cmd += ["--species_name", str(g("species_name"))] if g("species_name") else []
    else:
        cmd += ["--refdir", str(g("mitos_refdir"))] if g("mitos_refdir") else []
        cmd += ["--mitos_bin", str(g("mitos_bin"))] if g("mitos_bin") else []
    return cmd


def write_input_annotation(si: SampleInfo, dest: Path) -> str:
    """`mitom import`: the input GenBank record of this sample (annotated by MitoZ, MITOS2 or another program) as its
    round-1 annotation, named like the core's record. Returns an error message, '' when written."""
    if not si.orig_gb_path:
        return f"{si.sample_name}: no input GenBank record"
    recs = list(SeqIO.parse(str(si.orig_gb_path), "genbank"))
    rec = next((r for r in recs if si.orig_token in (r.id, r.name, safe_name(r.id))), recs[0] if len(recs) == 1 else None)
    if rec is None:
        return f"{si.sample_name}: record {si.orig_token} not found in {si.orig_gb_path}"
    if not any(f.type in ("CDS", "rRNA", "tRNA") for f in rec.features):
        return f"{si.sample_name}: the input record has no gene annotation"
    out = SeqRecord(rec.seq, id=si.internal_id, name=si.internal_id[:16], description="", features=list(rec.features))
    out.annotations.update(molecule_type="DNA", topology=si.topology,
                           data_file_division=rec.annotations.get("data_file_division", "INV"))
    SeqIO.write(out, str(dest), "genbank")
    return ""


def annotate_one(
    si: SampleInfo, fasta_file: Path, out_root: Path, phase_dirname: str,
    args: argparse.Namespace, total: int = 0, done_counter: Optional[List] = None,
) -> Tuple[str, bool, Path, Optional[Path]]:
    """Annotate one sample (engine_cmd, or the input record for `mitom import`): (internal_id, ok, workdir, the
    <workdir>/<internal_id>.gbf the engine wrote or None)."""
    internal_id = si.internal_id
    workdir     = out_root / phase_dirname / internal_id
    workdir.mkdir(parents=True, exist_ok=True)
    log_name = f"sample.{phase_dirname}.{internal_id}"
    lg = get_sample_logger(log_name, workdir / "sample.log")
    lg.info(f"[SAMPLE] internal_id={internal_id} sample_name={si.sample_name}")
    lg.info(f"[INPUT_FASTA] {fasta_file}")
    logs_dir    = workdir / "logs"
    stdout_path = logs_dir / "core.stdout.log"
    stderr_path = logs_dir / "core.stderr.log"
    if getattr(args, "engine", "native") == "input":          # mitom import: no annotation program runs
        t0 = time.time()
        logs_dir.mkdir(parents=True, exist_ok=True)
        err = write_input_annotation(si, workdir / f"{internal_id}.gbf")
        stderr_path.write_text(err + "\n" if err else "", encoding="utf-8")
        rc, dt = (1 if err else 0), time.time() - t0
    else:
        cmd = engine_cmd(fasta_file, workdir, internal_id, args, sample_lineage(si))
        lg.info("[CMD] " + " ".join(cmd))
        rc, dt = run_cmd_capture(cmd, cwd=workdir, stdout_path=stdout_path, stderr_path=stderr_path)
    if rc != 0:
        (logs_dir / "core.returncode.txt").write_text(str(rc) + "\n", encoding="utf-8")
    (logs_dir / "core.runtime_sec.txt").write_text(f"{dt:.2f}\n", encoding="utf-8")

    if done_counter is not None:
        lock: threading.Lock = done_counter[1]
        with lock:
            done_counter[0] += 1
            done_n = done_counter[0]
        logging.info(f"[PROGRESS] {done_n}/{total} '{si.sample_name}' done (rc={rc}, {dt:.0f}s)")
        if _STAGE_PROGRESS is not None:
            _STAGE_PROGRESS.update(done_n)

    if rc != 0:
        err_tail = tail_text(stderr_path.read_text(encoding="utf-8", errors="ignore") if stderr_path.exists() else "")
        lg.error(f"[FAIL] annotation core returned {rc}")
        lg.error("---- stderr tail ----\n" + err_tail)
        last = next((x.strip() for x in reversed(err_tail.splitlines()) if x.strip()), "no stderr")
        record_event(si.sample_name, f"annotation failed ({STAGE_LABEL.get(phase_dirname, phase_dirname)}, exit {rc}): "
                                     f"{last}  (sample log: {workdir / 'sample.log'})")
        close_sample_logger(log_name)
        return internal_id, False, workdir, None

    own = workdir / f"{internal_id}.gbf"                 # what the core and mito_engines.py write (exit 0 = written)
    close_sample_logger(log_name)
    return internal_id, True, workdir, own if own.is_file() else None


def fit_threads(args: argparse.Namespace, backend: str, n: int) -> None:
    """--threads per record lowered so that the records annotated at once (--jobs, at most n) fit the cores (Snakemake:
    --smk_cores; threads: the CPUs): 16 records x 48 threads with 64 cores ran one record at a time (2026-10-31: 148
    records in 57 min instead of about 4)."""
    cores = int(getattr(args, "smk_cores", 0) or 0) if backend == "snakemake" else (os.cpu_count() or 8)
    jobs = max(1, min(int(getattr(args, "max_tasks", 1) or 1), n))
    fit = max(1, (cores or (os.cpu_count() or 8)) // jobs)
    if int(getattr(args, "threads", 1) or 1) > fit:
        logging.warning(f"[ANNOT_STAGE] --threads {args.threads} x {jobs} records at a time exceeds {cores} cores: "
                        f"{fit} thread(s) per record, so that {jobs} records run at once")
        args.threads = fit


def run_annotation_stage(
    *, samples: List[SampleInfo], phase_dirname: str, out_root: Path,
    args: argparse.Namespace, fasta_selector,
) -> Dict[str, Dict]:
    backend = getattr(args, "scheduler_backend", "threads")
    if getattr(args, "engine", "native") == "input":            # mitom import: nothing to run
        backend = "threads"
    if backend == "auto":
        # Snakemake schedules nothing when every sample runs at once, and its process costs ~120 MB (a third of the
        # peak memory of a one-record run): threads then
        backend = "snakemake" if snakemake_available(getattr(args, "snakemake_bin", None)) \
                  and len(samples) > int(getattr(args, "max_tasks", 1) or 1) else "threads"
    fit_threads(args, backend, len(samples))
    logging.info(f"[ANNOT_STAGE] phase={phase_dirname}  backend={backend}  n={len(samples)}")
    global _STAGE_PROGRESS
    _STAGE_PROGRESS = Progress(STAGE_LABEL.get(phase_dirname, phase_dirname), len(samples))

    if backend == "snakemake":
        tasks: List[AnnotTask] = []
        for si in samples:
            tasks.append(AnnotTask(
                internal_id=si.internal_id, sample_name=si.sample_name,
                fasta_file=str(fasta_selector(si)), phase_dirname=phase_dirname,
                out_root=str(out_root), lineage=sample_lineage(si),
                cmd=engine_cmd(Path(fasta_selector(si)).resolve(), (out_root / phase_dirname / si.internal_id).resolve(),
                               si.internal_id, args, sample_lineage(si))))
        return run_snakemake_for_annotation_tasks(
            tasks=tasks, out_root=out_root, phase_name=phase_dirname, args=args, progress=_STAGE_PROGRESS)

    results: Dict[str, Dict] = {}
    total = len(samples)
    done_counter: List = [0, threading.Lock()]

    def _task(si: SampleInfo) -> Tuple:
        return annotate_one(si, fasta_selector(si), out_root, phase_dirname, args,
                            total=total, done_counter=done_counter)

    with cf.ThreadPoolExecutor(max_workers=int(args.max_tasks)) as ex:
        futs = {ex.submit(_task, si): si for si in samples}
        for fut in cf.as_completed(futs):
            si = futs[fut]
            try:
                internal_id, ok, workdir, best_gb = fut.result()
            except Exception as exc:
                record_event(si.sample_name, f"annotation failed ({STAGE_LABEL.get(phase_dirname, phase_dirname)}): "
                                             f"{type(exc).__name__}: {exc}", logging.ERROR)
                internal_id, ok, workdir, best_gb = si.internal_id, False, out_root / phase_dirname / si.internal_id, None
            results[internal_id] = {
                "internal_id": internal_id, "sample_name": si.sample_name,
                "ok": ok, "workdir": str(workdir), "best_gb": str(best_gb) if best_gb else None,
            }
    return results


def run_phase2(phase2_samples: List["SampleInfo"], args: argparse.Namespace, out_root: Path,
               reor_fasta_dir: Path) -> Dict[str, Dict]:
    """Round 2: annotate the reoriented sequences of `phase2_samples`. Nothing to do (every sample passed the
    --if_reannotation_2 auto check) -> no annotation stage at all: Snakemake would otherwise get an empty rule."""
    logging.info(f"[PHASE2] re-annotation: {len(phase2_samples)} samples")
    if not phase2_samples:
        return {}
    return run_annotation_stage(
        samples=phase2_samples, phase_dirname="03.anno_round2_reoriented",
        out_root=out_root, args=args,
        fasta_selector=lambda si: reor_fasta_dir / f"{si.internal_id}.fasta")
