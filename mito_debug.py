#!/usr/bin/env python3
"""mito_debug.py - what happened in a MitoM run, read from its log (`mitom debug <output folder | log file>`).

The run log (run.log of `mitom run` / `import`, pipeline.log of `annotate` / `reorient`, the log of `dedup`,
`transfer_metadata`, `replace_organism`) holds lines of a fixed form (mitom_log.py); with the output folder the
reports are read too. The report:
  run        command, MitoM version, start, end, duration; status: finished (exit code), stopped on an error (the
             exception, the step it stopped in, the traceback's last lines), or did not finish (killed / still running)
  steps      time and result of each step and of its parts
  records    every record that failed or was left out, with the first error of its own log (reports/failed_logs/,
             work/ with --keep_intermediate yes / --debug)
  commands   every program that failed: exit code, time, folder, the end of its stderr
  programs   programs that were missing, and how to install them
  warnings   grouped by kind, with one example each
  NCBI       records with table2asn ERRORs (ncbi_submission/validation_messages.tsv)
  reproduce  the run again in debug mode, the failed program by hand, a failed record on its own
Exit code: 0 = the run finished without a problem, 1 = problems found (above), 2 = no MitoM log found.

Usage:  mitom debug mitom_out          (or: mitom debug mitom_out/run.log)
"""
from __future__ import annotations

import argparse
import csv
import re
import shlex
import sys
import time
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402  (install hints)

LOGS = ("run.log", "pipeline.log", "dedup.log", "transfer_metadata.log", "replace_organism.log")
STDERR_MARK = "    | "              # mito_common.STDERR_MARK (not imported: mito_common needs Biopython)
STAMP = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)(?:,\d+)? \[?(DEBUG|INFO|WARNING|ERROR|CRITICAL)\]? ?(.*)$")
STEP = re.compile(r"^\[STEP (\d+)/(\d+) (\w+)\] (.*)$")
STEP_END = re.compile(r"^(?:(\w+) )?(done|failed) ([\d.]+)s(?:: (.*))?$", re.S)
CMD = re.compile(r"^\[CMD (?:rc=(-?\d+)|failed) ([\d.]+)s\] (.*)$")
RECORD = re.compile(r"^\[RECORD (.+?)\] (.*)$", re.S)
DONE = re.compile(r"^\[DONE\] exit=(-?\d+)(?: (.*))?$")
STEP_OLD = re.compile(r"^\[Step (\d+)/(\d+)\] (.*)$")                      # logs of MitoM 1.0.0
END_OLD = re.compile(r"^-> (.*?)\s*\((?:(\d+) h )?(?:(\d+) min )?(?:(\d+) s)?\)$", re.S)
MISSING = [re.compile(r"MitoM needs (\S+) for (.+?) \("), re.compile(r"\[core\] (\S+) not found"),
           re.compile(r"\[mito_\w+\] missing ([\w, -]+?) \("), re.compile(r"(MiTFi) \(mitfi_import\.jar")]
ERRORISH = re.compile(r"Error\b|ERROR|\bfailed\b|not found|No such file|[Kk]illed|out of memory|Cannot allocate memory")
FAILED_LISTS = {"annotation_failed.txt": "annotation failed", "final_write_failed.txt": "final GenBank not written",
                "reorient_failed.txt": "re-orientation failed", "metadata_transfer_failed.txt": "header transfer failed"}


@dataclass
class Step:
    n: int
    total: int
    key: str
    title: str = ""
    status: str = "running"            # done / failed / running (no end line: the run stopped in it)
    seconds: float = 0.0
    result: str = ""
    parts: List[Tuple[str, str, float, str]] = field(default_factory=list)   # (name, done / failed, seconds, result)


@dataclass
class Run:
    log: Path
    out: Path
    start: str = ""
    end: str = ""
    version: str = ""
    command: str = ""
    exit: Optional[int] = None
    error: str = ""                    # the exception the run stopped on, with the end of its traceback
    steps: List[Step] = field(default_factory=list)
    records: "OrderedDict[str, List[str]]" = field(default_factory=OrderedDict)
    commands: List[Dict] = field(default_factory=list)                         # failed programs
    warnings: "OrderedDict[str, List[str]]" = field(default_factory=OrderedDict)
    missing: "OrderedDict[str, str]" = field(default_factory=OrderedDict)      # program -> what needs it
    ncbi: "OrderedDict[str, str]" = field(default_factory=OrderedDict)         # record -> its table2asn ERRORs
    notes: List[str] = field(default_factory=list)                             # run-level failures ([FAIL] ...)


def find_log(path: Path) -> Optional[Path]:
    """The MitoM log of an output folder (or the file itself)."""
    if path.is_file():
        return path
    return next((path / n for n in LOGS if (path / n).is_file()), None)


def entries(log: Path) -> List[Tuple[str, str, str]]:
    """(time, level, message) of each log line; lines without a time stamp continue the message before them."""
    out: List[List[str]] = []
    for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
        m = STAMP.match(line)
        if m:
            out.append([m.group(1), m.group(2), m.group(3)])
        elif out:
            out[-1][2] += "\n" + line
    return [tuple(e) for e in out]


def split_cmd(text: str) -> Dict:
    """A [CMD ...] entry: exit code (None: did not start / timed out), seconds, command, folder, exception, the end of
    its stderr (lines after STDERR_MARK), and the programs that failed inside it (its own [CMD] lines, debug mode)."""
    first, _, rest = text.partition("\n")
    m = CMD.match(first)
    cmd, exc = m.group(3), ""
    if m.group(1) is None and ": " in cmd:                 # [CMD failed 0.0s] <command>: <exception>
        cmd, exc = cmd.rsplit(": ", 1)
    cmd, _, cwd = cmd.partition("  (in ")
    stderr = [ln[len(STDERR_MARK):] if ln.startswith(STDERR_MARK) else ln for ln in rest.splitlines()]
    inner = list(dict.fromkeys(mm.group(0) for ln in stderr
                               for mm in [re.search(r"\[CMD rc=(-?[1-9]\d*|-\d+) [\d.]+s\] .*", ln)] if mm))
    return {"rc": int(m.group(1)) if m.group(1) is not None else None, "seconds": float(m.group(2)), "cmd": cmd,
            "cwd": cwd.rstrip(")"), "exception": exc, "stderr": "\n".join(stderr).strip(), "inner": inner}


def parse(log: Path, out: Path) -> Run:
    run = Run(log=log, out=out)
    steps: Dict[Tuple[int, str], Step] = {}
    old_done = False                                       # MitoM 1.0.0: '[DONE]' without the exit code
    rows = entries(log)
    if rows:
        run.start, run.end = rows[0][0], rows[-1][0]
    for _, level, text in rows:
        first = text.split("\n", 1)[0]
        if re.match(r"\[ENV\] MitoM (\d|\?)", first):         # [ENV] MitoM <version>  <folder>, git <commit>
            run.version = first[len("[ENV] MitoM "):]
        elif first.startswith("[ENV] command: "):
            run.command = first[len("[ENV] command: "):]
        if STEP.match(first):
            n, total, key, rest = STEP.match(first).groups()
            st = steps.setdefault((int(n), key), Step(int(n), int(total), key))
            if rest.startswith("start: "):
                st.title = rest[len("start: "):]
                run.steps.append(st)
            elif STEP_END.match(rest):
                part, status, sec, result = STEP_END.match(rest).groups()
                if part:
                    st.parts.append((part, status, float(sec), (result or "").strip()))
                else:
                    st.status, st.seconds, st.result = status, float(sec), (result or "").strip()
            continue
        if STEP_OLD.match(first):                          # MitoM 1.0.0: '[Step n/N] title' ... '-> result (time)'
            n, total, title = STEP_OLD.match(first).groups()
            run.steps.append(Step(int(n), int(total), "-", title))
            continue
        if END_OLD.match(text) and run.steps and run.steps[-1].key == "-" and run.steps[-1].status == "running":
            res, h, mi, sec = END_OLD.match(text).groups()
            st = run.steps[-1]
            st.status, st.result = "done", res
            st.seconds = 3600 * int(h or 0) + 60 * int(mi or 0) + int(sec or 0)
            continue
        if first == "[DONE]":                              # MitoM 1.0.0: the end of a run, without its exit code
            old_done = True
            continue
        if CMD.match(first):
            c = split_cmd(text)
            if c["rc"] != 0:
                run.commands.append(c)
                if c["rc"] is None and "FileNotFoundError" in c["exception"]:
                    run.missing.setdefault(Path(shlex.split(c["cmd"])[0]).name if c["cmd"] else "?", "a step")
            continue
        m = RECORD.match(text)
        if m:
            run.records.setdefault(m.group(1), []).append(f"{level}: {m.group(2).strip()}")
            if m.group(2).startswith("NCBI ERROR"):
                run.ncbi[m.group(1)] = m.group(2).strip()
            continue
        if DONE.match(first):
            run.exit = int(DONE.match(first).group(1))
            continue
        if first.startswith("[ERROR] "):
            run.error = text[len("[ERROR] "):]
            continue
        for rx in MISSING:
            mm = rx.search(text)
            if mm:
                for prog in re.split(r",\s*", mm.group(1)):
                    run.missing.setdefault(prog.strip(), mm.group(2) if mm.re.groups > 1 else "a step")
        if first.startswith("[ENV] program ") and first.endswith(": not found"):
            prog = first[len("[ENV] program "):-len(": not found")]
            if prog in mito_env.REQUIRED:
                run.missing.setdefault(prog, "mitom run")
        if first.startswith("[FAIL] "):
            run.notes.append(first[len("[FAIL] "):])
        elif level in ("WARNING", "ERROR", "CRITICAL"):
            tag = re.match(r"^(\[[A-Z0-9_ ]+\])", first)
            run.warnings.setdefault(tag.group(1) if tag else first[:50], []).append(first)
    if run.exit is None and old_done:                      # 1.0.0 wrote '[FAIL] every record failed' after [DONE]
        run.exit = 1 if any("every record failed" in n for n in run.notes) else 0
    return run


def first_error(path: Path) -> str:
    """The error a log ends on: the exception of its last traceback, else its first line that reads like an error."""
    lines = [ln.rstrip() for ln in path.read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()]
    tb = [i for i, ln in enumerate(lines) if ln.startswith("Traceback (most recent call last)")]
    if tb:
        exc = [ln for ln in lines[tb[-1] + 1:] if not ln.startswith((" ", "\t"))]
        if exc:
            return exc[0]
    for ln in lines:
        msg = STAMP.match(ln).group(3) if STAMP.match(ln) else ln
        if ERRORISH.search(msg):
            return msg
    return ""


def record_logs(out: Path, name: str, text: str) -> List[Path]:
    """The logs of a record: its copies in reports/failed_logs/, the sample log named in the [RECORD] line, the logs of
    its work folder (kept with --keep_intermediate yes / --debug)."""
    logs = sorted((out / "reports" / "failed_logs").glob(f"{name}.*")) if (out / "reports" / "failed_logs").is_dir() else []
    m = re.search(r"sample log: (\S+?/(0[13]\.[\w.]+)/([^/\s]+))/sample\.log", text)
    if m:
        for d in (Path(m.group(1)), out / "work" / m.group(2) / m.group(3)):
            if d.is_dir():
                logs += [d / "logs" / "core.stderr.log", d / "sample.log"]
                break
    return [p for p in logs if p.is_file()]


def record_input(out: Path, name: str, events: List[str]) -> Optional[Path]:
    """A file with just this record, to run it alone: its original GenBank record (GenBank input), else its FASTA of
    the run (kept with --keep_intermediate yes / --debug)."""
    if (out / "original_annotation" / f"{name}.gb").is_file():
        return out / "original_annotation" / f"{name}.gb"
    m = re.search(r"/0[13]\.[\w.]+/([^/\s]+)/sample\.log", "\n".join(events))
    for d in ((out / "00.internal_fasta", out / "work" / "00.internal_fasta") if m else ()):
        if (d / f"{m.group(1)}.fasta").is_file():
            return d / f"{m.group(1)}.fasta"
    return None


def add_output(run: Run) -> None:
    """What the output folder adds to the log: the problem lists of reports/, the first error of each failed record's
    own log, the NCBI ERRORs of ncbi_submission/validation_messages.tsv."""
    rep = run.out / "reports"
    lists = dict(FAILED_LISTS)
    lists.update({p.name: "no orientation gene" for p in rep.glob("no_*_samples.txt")} if rep.is_dir() else {})
    for name, what in lists.items():
        f = rep / name
        if f.is_file():
            for ln in f.read_text(encoding="utf-8", errors="replace").splitlines():
                if ln.strip() and not ln.startswith("#") and ln.strip() not in run.records:
                    run.records[ln.strip()] = [f"WARNING: {what} (reports/{name})"]
    for name, events in run.records.items():
        for p in record_logs(run.out, name, "\n".join(events)):
            err = first_error(p)
            if err:
                events.append(f"first error in {p.relative_to(run.out) if run.out in p.parents else p}: {err}")
                break
    val = run.out / "ncbi_submission" / "validation_messages.tsv"
    if val.is_file():
        with val.open(encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("severity") == "ERROR" and row.get("record") not in run.ncbi:
                    run.ncbi[row["record"]] = f"NCBI ERROR: {row.get('code', '')} {row.get('message', '')}"[:200]


def relocated(text: str, out: Path) -> str:
    """Paths of the run's work folders in text as they are after the run: moved into work/ (--keep_intermediate yes,
    --debug)."""
    try:
        from mitom_output import WORK_DIRS
    except ImportError:                     # without Biopython
        return text
    base = out.resolve()
    for name in WORK_DIRS:
        if not (base / name).exists() and (base / "work" / name).exists():
            text = text.replace(f"{base}/{name}", f"{base}/work/{name}")
    return text


def rerun(command: str, inp: Optional[str] = None, out: Optional[str] = None) -> str:
    """The run's command in debug mode, with another input (-i) and / or output folder (-o)."""
    try:
        words = shlex.split(command)
    except ValueError:
        return command + " --debug"
    res: List[str] = []
    i = 0
    while i < len(words):
        w = words[i]
        if w in ("-i", "--input") and inp and i + 1 < len(words):
            res += [w, inp]
            i += 2
            continue
        if w in ("-o", "--out", "--out_root") and out and i + 1 < len(words):
            res += [w, out]
            i += 2
            continue
        res.append(w)
        i += 1
    if "--debug" not in res:
        res.append("--debug")
    return " ".join(shlex.quote(w) for w in res)


def status(run: Run) -> str:
    running = next((s for s in run.steps if s.status == "running"), None)
    if run.error:
        where = f" in step {running.n}/{running.total} {running.key}" if running else ""
        return f"STOPPED ON AN ERROR{where}: {run.error.splitlines()[0]}"
    if run.exit is None:
        age = time.time() - run.log.stat().st_mtime
        where = f" in step {running.n}/{running.total} {running.key} ({running.title})" if running else ""
        return (f"DID NOT FINISH{where}: the log ends at {run.end}"
                + (" - written in the last minute, the run may still be going" if age < 60 else
                   " (killed? out of memory or disk space? see the last lines of the log)"))
    if run.exit == 0:
        return "finished (exit 0)"
    return f"FAILED (exit {run.exit})" + (": " + "; ".join(run.notes) if run.notes else "")


def report(run: Run) -> Tuple[str, bool]:
    """The report text and whether it found a problem."""
    L: List[str] = [f"MitoM run: {run.out}  (log: {run.log.name})"]
    if run.command:
        L.append(f"  command : {run.command}")
    if run.version:
        L.append(f"  MitoM   : {run.version}")
    L.append(f"  time    : {run.start} - {run.end}")
    L.append(f"  status  : {status(run)}")
    if run.error:
        L += ["", "Error (end of the traceback; the whole one is in the log)"]
        L += [f"  {ln}" for ln in run.error.splitlines()[-8:]]
    if run.steps:
        L += ["", "Steps"]
        for s in run.steps:
            L.append(f"  {s.n:>2}/{s.total} {s.key:12} {s.seconds:8.1f} s  {s.status:7} {s.result[:110]}")
            for name, st, sec, res in s.parts:
                L.append(f"        {name:12} {sec:8.1f} s  {st:7} {res[:110]}")
    if run.records:
        L += ["", f"Records with problems ({len(run.records)})"]
        for name, events in run.records.items():
            L.append(f"  {name}")
            L += [f"      {e[:300]}" for e in events]
    if run.commands:
        L += ["", f"Failed programs ({len(run.commands)})"]
        for c in run.commands[:20]:
            how = f"exit {c['rc']}" if c["rc"] is not None else f"did not run: {c['exception']}"
            L.append(f"  [{how}, {c['seconds']:.1f} s] {c['cmd'][:400]}" + (f"\n      in {c['cwd']}" if c["cwd"] else ""))
            L += [f"      failed inside it: {x[:300]}" for x in c["inner"][:5]]
            L += [f"      | {ln}" for ln in c["stderr"].splitlines()[-6:]]
    if run.missing:
        L += ["", "Missing programs"]
        for prog, what in run.missing.items():
            L.append(f"  {prog} (needed by {what}): {mito_env.install_hint(prog)}")
    if run.ncbi:
        L += ["", f"NCBI ERRORs ({len(run.ncbi)} records; ncbi_submission/validation_messages.tsv)"]
        L += [f"  {rec}: {text}" for rec, text in run.ncbi.items()]
    if run.warnings:
        L += ["", f"Warnings ({sum(len(v) for v in run.warnings.values())}, by kind)"]
        for tag, msgs in sorted(run.warnings.items(), key=lambda kv: -len(kv[1])):
            L.append(f"  {len(msgs):4}x  {msgs[0][:160]}")
    problem = bool(run.error or run.exit != 0 or run.records or run.commands or run.missing or run.ncbi)
    L += ["", "Reproduce"]
    if run.command:
        L.append(f"  the run in debug mode : {rerun(run.command, out=str(run.out) + '_debug')}")
    if run.commands:
        c = {k: relocated(v, run.out) if isinstance(v, str) else v for k, v in run.commands[0].items()}
        gone = c["cwd"] and not Path(c["cwd"]).is_dir()
        L.append("  the first failed program by hand:"
                 + (f"\n      cd {shlex.quote(c['cwd'])} &&" if c["cwd"] else "") + f"\n      {c['cmd']}"
                 + ("\n      (its folder was removed after the run: run with --debug to keep it)" if gone else ""))
    for name in list(run.records)[:3]:
        one = record_input(run.out, name, run.records[name])
        if run.command and one:
            L.append(f"  record {name} alone   : {rerun(run.command, inp=str(one), out=str(run.out) + '_' + name)}")
    if not problem:
        L.append("  (nothing to reproduce: no problem found)")
    return "\n".join(L) + "\n", problem


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("path", help="output folder of a MitoM run, or its log file (run.log, pipeline.log ...)")
    a = p.parse_args(argv)
    path = Path(a.path).expanduser()
    log = find_log(path)
    if log is None:
        print(f"mitom debug: no MitoM log ({', '.join(LOGS)}) in {path}", file=sys.stderr)
        return 2
    run = parse(log, log.parent)
    add_output(run)
    text, problem = report(run)
    groups = sorted((log.parent / "reports" / "groups").glob("*.run.log"))
    if groups:                       # --clade auto over several phyla: one run per group (mitom_groups.py)
        text = text.split("\nSteps", 1)[0].rstrip() + f"\n\n{len(groups)} groups (phylum / genetic code), " \
               "one run each (reports/groups/<group>.run.log):\n"
        for g in groups:
            gr = parse(g, log.parent)
            gtext, gproblem = report(gr)
            problem = problem or gproblem
            name = g.name[:-len(".run.log")]
            text += f"  {name:<28} {status(gr)}\n" if not gproblem else \
                f"\n---- group {name} ----\n" + gtext + "\n"
    sys.stdout.write(text)
    return 1 if problem else 0


if __name__ == "__main__":
    sys.exit(main())
