"""mitom_prepare.py - before a workflow runs: the environment check and the [ENV] lines of the log, the tools and
the clade (--clade auto), `mitom import` results read in, one run per output folder (.mitom.lock)."""
from __future__ import annotations

import argparse
from mitom_output import where
import logging
import os
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import List, Optional

import mitom_groups
from mitom_context import RunContext
from mito_common import command_line, write_gb
from mito_taxa import ALL_CLADES, BASE_CLADES, CLADE_CODES, CLADE_PHYLUM, EXTRA_CLADES
from mitom_inputs import auto_clade
from mitom_snakemake import snakemake_available, which_or_none


def check_environment(args: argparse.Namespace) -> List[str]:
    problems: List[str] = []
    try:
        import Bio
        logging.info(f"[ENV] Biopython {getattr(Bio, '__version__', '?')} ✔")
    except Exception:
        pass
    major, minor = sys.version_info.major, sys.version_info.minor
    if (major, minor) < (3, 8):
        raise RuntimeError(f"Python 3.8+ required, found {major}.{minor}")
    logging.info(f"[ENV] Python {major}.{minor} ✔")

    if getattr(args, "cmd", None) in ("annotate", "run"):
        import mito_env
        seen = set()
        for step in run_steps(args):                       # the programs of the steps this run does (mito_env.REQUIRES)
            for prog in mito_env.REQUIRES[step]:
                if prog in seen:
                    continue
                seen.add(prog)
                found = mito_env.find(prog, None)
                if found:
                    logging.info(f"[ENV] {prog} ✔  {found}")
                else:
                    problems.append(f"[ENV] MitoM needs {prog} for {mito_env.STEP_NAMES[step]} (install.sh; mitom doctor)")
        if not mito_env.mitos_refdir(getattr(args, "mitos_refdir", None)):
            problems.append("[ENV] MitoM needs the MITOS2 reference data refseq89m (--mitos_refdir; install.sh)")

    smk_bin = getattr(args, "snakemake_bin", None)
    backend = getattr(args, "scheduler_backend", "threads")
    if backend in ("snakemake", "auto"):
        if snakemake_available(smk_bin):
            smk_exe = smk_bin or which_or_none("snakemake") or "snakemake"
            try:
                sv = subprocess.run([smk_exe, "--version"], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, text=True, timeout=10)
                smk_ver = (sv.stdout or sv.stderr or "").strip().splitlines()[0] \
                          if sv.returncode == 0 else "?"
                logging.info(f"[ENV] Snakemake {smk_ver} ✔")
            except Exception:
                logging.info("[ENV] Snakemake available (version check failed)")
        else:
            if backend == "snakemake":
                raise RuntimeError("[ENV] --scheduler_backend snakemake requested but snakemake not found.")
            logging.info("[ENV] Snakemake not found, will use threads backend.")

    if hasattr(args, "input") and args.input:
        in_p = Path(args.input)
        if not in_p.exists():
            raise RuntimeError(f"[ENV] Input path does not exist: {in_p}")
        logging.info(f"[ENV] Input path exists ✔  ({in_p})")

    if hasattr(args, "infq") and args.infq:
        fq_p = Path(args.infq)
        if not fq_p.exists():
            problems.append(f"[ENV] --infq path does not exist: {fq_p}")
        else:
            logging.info(f"[ENV] FASTQ dir exists ✔  ({fq_p})")

    return problems


def git_state(folder: Path) -> str:
    """', git <commit>[, uncommitted changes]' for a MitoM folder that is a git checkout, else ''."""
    try:
        head = subprocess.run(["git", "-C", str(folder), "rev-parse", "--short", "HEAD"], capture_output=True,
                              text=True, timeout=5)
        if head.returncode:
            return ""
        dirty = subprocess.run(["git", "-C", str(folder), "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True, timeout=5).stdout.strip()
        return f", git {head.stdout.strip()}" + (", uncommitted changes" if dirty else "")
    except (OSError, subprocess.SubprocessError):
        return ""


def log_environment(args: argparse.Namespace) -> None:
    """[ENV] lines at the start of a run's log: MitoM (version, folder, git commit), the command as typed, the working
    folder, Python and the machine; with --debug also every program with its version, the Python packages, the config
    file, the reference and data folders, the environment variables, memory and free disk space."""
    import platform

    import mito_env
    here = Path(__file__).resolve().parent
    logging.info(f"[ENV] MitoM {os.environ.get('MITOM_VERSION', '?')}  {here}{git_state(here)}")
    cmd = (["mitom"] if os.environ.get("MITOM_LAUNCHER") else [sys.executable, sys.argv[0]]) + sys.argv[1:]
    logging.info(f"[ENV] command: {command_line(cmd)}")
    logging.info(f"[ENV] working folder: {Path.cwd()}")
    logging.info(f"[ENV] Python {platform.python_version()} ({sys.executable}); {platform.platform()}; "
                 f"{os.cpu_count()} CPUs")
    if not getattr(args, "debug", False):
        return
    for name in mito_env.TOOLS:
        exe = mito_env.find(name)
        logging.info(f"[ENV] program {name}: " + (f"{exe}  ({mito_env.program_version(exe, name) or 'version ?'})"
                                                   if exe else "not found"))
    for pkg in ("Bio", "numpy", "sklearn", "matplotlib"):
        try:
            logging.info(f"[ENV] package {pkg} {getattr(__import__(pkg), '__version__', '?')}")
        except ImportError:
            logging.info(f"[ENV] package {pkg} not installed")
    cfg = mito_env.config_path()
    logging.info(f"[ENV] config {cfg}: " + (cfg.read_text(encoding="utf-8").strip() if cfg.is_file() else "none"))
    logging.info(f"[ENV] MITOS2 reference data: {mito_env.mitos_refdir(getattr(args, 'mitos_refdir', None)) or 'not found'}")
    logging.info(f"[ENV] data folder (MitoM-ML models): {mito_env.data_dir()}")
    for k in sorted(os.environ):
        if k.startswith(("MITOM", "MITOS", "CONDA", "OMP_", "JAVA", "YZZ_")) or k in ("PATH", "LANG", "TMPDIR"):
            logging.info(f"[ENV] ${k}={os.environ[k]}")
    try:
        mem = {ln.split(":")[0]: ln.split()[1] for ln in Path("/proc/meminfo").read_text().splitlines()[:3]}
        logging.info(f"[ENV] memory {int(mem['MemTotal']) // 1024} MB, available {int(mem['MemAvailable']) // 1024} MB")
    except (OSError, KeyError, ValueError, IndexError):
        pass
    out = Path(getattr(args, "out_root", ".") or ".")
    try:
        logging.info(f"[ENV] free disk space in {out}: {shutil.disk_usage(out).free // 2**30} GB")
    except OSError:
        pass


def run_steps(args: argparse.Namespace) -> List[str]:
    """The steps of mito_env.REQUIRES this run does (the core only with --engine native)."""
    steps = ["core", "trna_ml"] if getattr(args, "engine", "native") == "native" else []
    steps += ["introns"] if getattr(args, "intron_fix", "no") != "no" else []
    steps += ["submission"] if getattr(args, "curate", "no") != "no" else []
    steps += ["rrna_hmm"] if getattr(args, "rrna_hmm", "no") != "no" else []
    return steps


def preflight(args: argparse.Namespace) -> None:
    problems = check_environment(args)
    if problems:
        for p in problems:
            logging.warning(p)
        logging.warning(f"[ENV] {len(problems)} warning(s) above — pipeline continues but may fail.")


def announce(msg: str) -> None:
    """A decision made before the run's log file exists (resolve_tools): log it and show it on the console."""
    logging.info(msg)
    print(msg, file=sys.stderr)


def resolve_tools(args: argparse.Namespace, ctx: Optional[RunContext] = None) -> str:
    """Fill unset tool locations from mito_env (config file, PATH, conda environments) and check the clade; the
    --clade auto plan goes to ctx.plan (the run writes reports/clade_decision.tsv from it).
    Returns an error message, '' when fine."""
    import mito_env
    if args.clade == "auto" and args.cmd == "import":      # results of other programs: one clade for the run
        clade, code, phylum, err = auto_clade(Path(args.input), getattr(args, "genetic_code", None))
        if err and not err.startswith("WARNING"):
            return f"ERROR: {err}"
        if err:
            print(f"\n{'=' * 68}\n  [CLADE_DEFAULT] {err}\n{'=' * 68}\n", file=sys.stderr)
        args.clade = clade
        if getattr(args, "genetic_code", None) is None:
            args.genetic_code = code
        if hasattr(args, "rrna_hmm") and not getattr(args, "phylum", None):
            args.phylum = phylum
        announce(f"[CLADE] auto: {phylum or 'no taxonomy'} -> --clade {args.clade} --genetic_code {args.genetic_code}")
    elif args.clade == "auto":                             # per record: taxonomy, else blastx (mitom_groups.py)
        plan = mitom_groups.decide(args)
        if plan.error:
            return f"ERROR: {plan.error}"
        if ctx is not None:
            ctx.plan = plan
        if plan.warnings:                                # one box: the first ones, all in clade_decision.tsv
            shown = plan.warnings[:5] + ([f"... {len(plan.warnings) - 5} more: {where(args, 'reports/clade_decision.tsv')} (column note)"]
                                         if len(plan.warnings) > 5 else [])
            print(f"\n{'=' * 68}\n  [CLADE] WARNING ({len(plan.warnings)}):\n" + "\n".join(f"  - {w}" for w in shown)
                  + f"\n{'=' * 68}\n", file=sys.stderr)
        groups = plan.groups()
        if len(groups) > 1:
            if args.cmd != "run":
                got = ", ".join(f"{mitom_groups.group_label(k)} ({len(v)})" for k, v in groups.items())
                return (f"ERROR: --clade auto: the input has several phyla / genetic codes ({got}); `mitom run` "
                        f"annotates them group by group, `mitom {args.cmd}` needs one: give --clade and --genetic_code")
            args.clade_groups = len(groups)              # main(): workflow_groups, one run per group
            return ""
        if not groups:
            return f"ERROR: no sequence records in {args.input} (empty file, or not GenBank / FASTA)"
        (clade, code, phylum), rows = next(iter(groups.items()))
        args.clade, args.genetic_code = clade, code
        if hasattr(args, "rrna_hmm") and not getattr(args, "phylum", None):
            args.phylum = phylum
        how = Counter(r.source for r in rows)
        announce(f"[CLADE] auto: {phylum} ({', '.join(f'{n} by {s}' for s, n in how.items())}) -> --clade {clade} "
                 f"--genetic_code {code}")
    if args.clade not in ALL_CLADES:
        return (f"ERROR: --clade '{args.clade}' is not a supported clade. Choose one of:\n  " + ", ".join(ALL_CLADES)
                + "\n  (Animal = all-animal protein database, for any phylum not listed.)")
    if hasattr(args, "rrna_hmm") and not getattr(args, "phylum", None) and args.clade != "Animal":
        args.phylum = CLADE_PHYLUM.get(args.clade)       # explicit --clade: its phylum selects the rRNA HMM library
    if getattr(args, "genetic_code", 0) is None:
        if args.clade in EXTRA_CLADES:
            return (f"ERROR: --clade {args.clade} needs --genetic_code (e.g. 4 Cnidaria, 9 Platyhelminthes, "
                    "5 most invertebrates; NCBI table of your taxon)")
        args.genetic_code = CLADE_CODES.get(args.clade, 5)
        if args.clade == "Chordata" and getattr(args, "input", None):
            # tunicates 13, lancelets 5: the GenBank records' own code, as --clade auto (GLM 3.2, 2026-11-05)
            clade, code, _, err = auto_clade(Path(args.input))
            if clade == "Chordata" and code and not err:
                args.genetic_code = code
            elif "differs between records" in err:
                # vertebrates 2 + tunicates 13 / lancelets 5 in one input: never one code for all (Cursor F1, 2026-11-06)
                return (f"ERROR: --clade Chordata: {err.replace('--clade auto: ', '')}. `mitom run --clade auto` "
                        f"annotates each code separately; or give --genetic_code")
    if not args.intron_bin:
        args.intron_bin = mito_env.bin_with("exonerate", "cmsearch", "blastp", "makeblastdb", "tblastn", "blastn")
    args.mitos_bin = mito_env.mitos_bin(args.mitos_bin)
    args.mitos_refdir = mito_env.mitos_refdir(args.mitos_refdir)
    engine = getattr(args, "engine", "native")
    if engine == "mitoz":
        import mito_engines
        if args.clade not in BASE_CLADES:
            return (f"ERROR: --engine mitoz: MitoZ has no protein database for --clade {args.clade} (MitoZ clades: "
                    f"{', '.join(BASE_CLADES)}); use --engine native or mitos2")
        if not mito_env.find("mitoz", getattr(args, "mitoz_bin", None)):
            return f"ERROR: --engine mitoz: MitoZ not found (install: {mito_engines.MITOZ_INSTALL}; or give --mitoz_bin)"
    if getattr(args, "mitos_args", ""):
        import mito_engines
        if engine != "mitos2":
            return "ERROR: --mitos_args is for --engine mitos2 (MITOS2 options)"
        try:
            mito_engines.mitos_extra_args(args.mitos_args)
        except ValueError as e:
            return f"ERROR: {e}"
    if engine == "mitos2" and not (args.mitos_bin and args.mitos_refdir):
        import mito_engines
        return ("ERROR: --engine mitos2: MITOS2 (runmitos) or its reference data refseq89m not found (install: "
                f"{mito_engines.MITOS2_INSTALL}; reference data: mitom setup; or give --mitos_bin / --mitos_refdir)")
    if engine == "input":                                  # mitom import: nothing can annotate again
        args.if_reannotation_2 = "no"
    return ""


def prepare_import(args: argparse.Namespace) -> str:
    """mitom import: existing results (mito_engines.import_results) -> <out>/00.imported/imported.gb, which the run
    workflow then reads with engine 'input' (no annotation program runs, round 2 off). Returns an error, '' when fine."""
    import mito_engines
    try:
        recs, notes = mito_engines.import_results(Path(args.input), [Path(x) for x in args.fasta], args.genetic_code,
                                                  args.topology)
    except (ValueError, OSError) as e:
        return f"ERROR: {e}"
    if not recs:
        return f"ERROR: no annotated GenBank record and no MITOS2 result (result.bed) found in {args.input}"
    out = Path(args.out_root).resolve() / "00.imported"
    out.mkdir(parents=True, exist_ok=True)
    write_gb(recs, out / "imported.gb")
    for note in notes:
        announce(f"[IMPORT] {note}")
    announce(f"[IMPORT] {len(recs)} annotated record(s) from {args.input}")
    args.input_given = args.input                  # protected from clear_previous_work (Codex review, 2026-11-06)
    args.input, args.engine = str(out / "imported.gb"), "input"
    return ""


LOCK_NAME = ".mitom.lock"


def lock_output(out_root: Path) -> str:
    """One run per output folder: out_root/.mitom.lock holds the pid of the run writing there (2026-10-31 stress test:
    two runs started into one folder both failed with half the records). A lock of a process that no longer runs
    (killed) is taken over. Removed at exit. Returns an error message, '' when the folder is ours."""
    import atexit
    try:
        out_root.mkdir(parents=True, exist_ok=True)
        (out_root / ".mitom.write_test").write_text("")
        (out_root / ".mitom.write_test").unlink()
    except OSError as exc:
        return f"ERROR: cannot write into the output folder {out_root}: {exc.strerror or exc} (give -o a folder you can write to)"
    lock = out_root / LOCK_NAME
    me = str(os.getpid())
    for _ in range(2):
        try:
            fd = os.open(str(lock), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            try:
                held = lock.read_text().split()
                other = int(held[0])
                if not _lock_holder_alive(other, held[1] if len(held) > 1 else ""):
                    raise ProcessLookupError(other)
            except (ValueError, IndexError, OSError):
                lock.unlink(missing_ok=True)           # stale (the run was killed): take it over
                continue
            return (f"ERROR: another MitoM run (pid {other}) is writing into {out_root}; wait for it, stop it, or "
                    f"use another -o folder (a stale lock: delete {lock})")
        with os.fdopen(fd, "w") as fh:
            fh.write(f"{me} {_process_start(os.getpid())}".rstrip() + "\n")
        atexit.register(_release_lock, lock, me)
        return ""
    return f"ERROR: cannot lock {lock}"


def _process_start(pid: int) -> str:
    """Start time of a process (clock ticks since boot, /proc/<pid>/stat field 22); '' where there is no /proc."""
    try:
        return Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()[19]
    except (OSError, IndexError):
        return ""


def _lock_holder_alive(pid: int, started: str = "") -> bool:
    """Whether the run of a lock still runs: a process of another user counts (os.kill raises PermissionError; it was
    taken for a stale lock); a pid now used by another process (its start time differs from the one in the lock)
    does not (N8, 2026-11-05)."""
    try:
        os.kill(pid, 0)
    except PermissionError:
        return True
    except OSError:
        return False
    now = _process_start(pid)
    return not (started and now and now != started)


def _release_lock(lock: Path, me: str) -> None:
    """atexit: remove the lock when it is still ours (D19: it may be gone or taken over meanwhile)."""
    try:
        if lock.read_text().split()[:1] == [me]:
            lock.unlink(missing_ok=True)
    except OSError:
        pass


def error_exit_code(message: str) -> int:
    """Exit code of a refused run: 1 when a program or the reference data is missing (the run cannot be done here),
    2 when the command line or the input is wrong (MitoM's exit codes: 0 done, 1 failure, 2 usage)."""
    return 1 if "not found" in message else 2
