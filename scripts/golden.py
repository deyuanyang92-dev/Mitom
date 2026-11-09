#!/usr/bin/env python3
"""golden.py - whole-workflow regression check: the same input must give the same output (refactoring guard).

Each case runs one `mitom` command on files in tests/data; its text outputs are normalised (output folder, repository
and home paths replaced by <OUT>, <REPO>, <HOME>) and compared with tests/data/golden/<case>/. Not compared: run logs
(times), plots (binary), .sqn (table2asn date stamp), intermediate folders (01.*, 02.* ...; only the arranged output).

The reference outputs depend on the external programs (BLAST, exonerate, Infernal, HMMER, table2asn, refseq89m):
they were made on the machine named in tests/data/golden/MACHINE.txt. On another machine, or after an intended change
of the output, rebuild them with --update and commit them together with the change that explains the difference.

Usage: python3 scripts/golden.py [--update] [--cases a,b] [--jobs 6] [--keep DIR]
       MITOM_GOLDEN=1 pytest tests/test_golden.py         (the same check as a test; skipped without MITOM_GOLDEN)
"""
import argparse
import concurrent.futures as cf
import difflib
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "tests" / "data"
GOLDEN = DATA / "golden"
COMMON = ["--threads", "2"]
CLASSIC = ["--output_layout", "classic"]          # the layout the first references were made in
# case -> (input builder, mitom_pipeline.py / script arguments); {IN} = input path, {OUT} = output folder
CASES: Dict[str, List[str]] = {
    "run_annelid_gb": ["mitom_pipeline.py", "run", "-i", "{DATA}/NC_000931.gb", "-o", "{OUT}", *COMMON, *CLASSIC],
    "run_coral_intron": ["mitom_pipeline.py", "run", "-i", "{DATA}/NC_006902.gb", "-o", "{OUT}", *COMMON, *CLASSIC],
    "run_two_phyla": ["mitom_pipeline.py", "run", "-i", "{IN}/two", "-o", "{OUT}", *COMMON, *CLASSIC],
    "run_fasta_blind": ["mitom_pipeline.py", "run", "-i", "{IN}/NC_010217.fasta", "-o", "{OUT}", *COMMON, *CLASSIC],
    "import_mitoz": ["mitom_pipeline.py", "import", "-i", "{DATA}/KY753833_mitoz.gbf", "-o", "{OUT}", "--clade",
                     "Annelida-segmented-worms", "--genetic_code", "5", *COMMON, *CLASSIC],
    # issues #7 / #8 / #9 (NCBI MZ052209.1, OZ287881.1, OZ346613.1), the default grouped layout
    "run_issues": ["mitom_pipeline.py", "run", "-i", "{DATA}/issues", "-o", "{OUT}", *COMMON],
    "annotate_gb": ["mitom_pipeline.py", "annotate", "-i", "{DATA}/NC_000931.gb", "-o", "{OUT}", *COMMON],
    "curate_mitoz": ["mito_curate.py", "-i", "{DATA}/OW387151_mitoz.gbf", "-o", "{OUT}"],
}
SKIP_SUFFIX = (".pdf", ".png", ".sqn")


def make_inputs(d: Path) -> None:
    """Inputs built from tests/data: a folder with two records of two phyla, and a bare FASTA (sequence only)."""
    from Bio import SeqIO
    (d / "two").mkdir(parents=True, exist_ok=True)
    for acc in ("NC_001636", "NC_008164"):                      # Mollusca (code 5) + Cnidaria (code 4)
        shutil.copy(DATA / f"{acc}.gb", d / "two" / f"{acc}.gb")
    rec = SeqIO.read(str(DATA / "NC_010217.gb"), "genbank")
    (d / "NC_010217.fasta").write_text(f">NC_010217\n{rec.seq}\n")


def keep(rel: str) -> bool:
    top = rel.split("/")[0]
    return not (re.match(r"\d\d\.", top) or rel.endswith(SKIP_SUFFIX) or Path(rel).name.endswith(("run.log", "pipeline.log")))


DATE = re.compile(r"\b\d\d-(?:JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)-\d{4}\b")


def normalise(text: str, out: Path, inp: Path) -> str:
    for p, tag in ((out, "<OUT>"), (inp, "<IN>"), (REPO, "<REPO>"), (Path.home(), "<HOME>")):
        text = text.replace(str(p.resolve()), tag).replace(str(p), tag)
    return DATE.sub("<DATE>", text)                 # (the LOCUS date of records written today)


def outputs(out: Path, inp: Path) -> Dict[str, str]:
    got = {}
    for p in sorted(out.rglob("*")):
        rel = p.relative_to(out).as_posix()
        if p.is_file() and keep(rel):
            got[rel] = normalise(p.read_text(errors="replace"), out, inp)
    return got


def run_case(name: str, work: Path, inp: Path) -> Dict[str, str]:
    out = work / name
    shutil.rmtree(out, ignore_errors=True)
    args = [a.format(DATA=DATA, IN=inp, OUT=out) for a in CASES[name]]
    env = {k: v for k, v in os.environ.items() if k not in ("MITOM_LAUNCHER", "MITOM_DEBUG")}
    r = subprocess.run([sys.executable, "-W", "ignore", str(REPO / args[0]), *args[1:]], capture_output=True,
                       text=True, env=env, cwd=work)
    if r.returncode != 0:
        raise RuntimeError(f"{name}: exit {r.returncode}\n{r.stderr[-3000:]}")
    return outputs(out, inp)


def compare(name: str, got: Dict[str, str]) -> List[str]:
    """Differences (a few lines each) between the outputs of a case and its reference."""
    ref_dir = GOLDEN / name
    want = {p.relative_to(ref_dir).as_posix(): DATE.sub("<DATE>", p.read_text()) for p in sorted(ref_dir.rglob("*")) if p.is_file()}
    msgs = [f"{name}: missing output {f}" for f in sorted(set(want) - set(got))]
    msgs += [f"{name}: new output {f}" for f in sorted(set(got) - set(want))]
    for f in sorted(set(want) & set(got)):
        if want[f] != got[f]:
            d = list(difflib.unified_diff(want[f].splitlines(), got[f].splitlines(), "golden/" + f, "now/" + f, n=0,
                                          lineterm=""))
            msgs.append(f"{name}: {f} differs\n" + "\n".join(d[:30]) + ("\n..." if len(d) > 30 else ""))
    return msgs


def write(name: str, got: Dict[str, str]) -> None:
    d = GOLDEN / name
    shutil.rmtree(d, ignore_errors=True)
    for rel, text in got.items():
        (d / rel).parent.mkdir(parents=True, exist_ok=True)
        (d / rel).write_text(text)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--update", action="store_true", help="rewrite tests/data/golden from this run")
    p.add_argument("--cases", default=",".join(CASES))
    p.add_argument("--jobs", type=int, default=6)
    p.add_argument("--keep", type=Path, help="run in this folder and keep it (default: a temporary folder)")
    a = p.parse_args(argv)
    names = a.cases.split(",")
    work = a.keep or Path(tempfile.mkdtemp(prefix="mitom_golden_"))
    work.mkdir(parents=True, exist_ok=True)
    inp = work / "inputs"
    make_inputs(inp)
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        res = dict(zip(names, ex.map(lambda n: run_case(n, work, inp), names)))
    if a.update:
        for n, got in res.items():
            write(n, got)
        (GOLDEN / "MACHINE.txt").write_text(f"{platform.node()}\n")
        print(f"[golden] {len(names)} case(s) written to {GOLDEN}")
        msgs = []
    else:
        msgs = [m for n, got in res.items() for m in compare(n, got)]
        print("\n".join(msgs) if msgs else f"[golden] {len(names)} case(s): outputs identical")
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    return 1 if msgs else 0


if __name__ == "__main__":
    sys.exit(main())
