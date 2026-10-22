"""Run annotation tools on the bare benchmark records and lift every result to NCBI coordinates.

Tools:
  ours    batch_mitoz.py run (current defaults) + mito_curate.py
  mitom   mitom run (one command: cores, MitoM steps, curation; final_annotation/)
  mitoz   MitoZ `annotate` unchanged (no post-processing)
  mitos2  MITOS2 runmitos, full mode (proteins, rRNA, tRNA), refseq89m
  mitos1  MITOS 1 (bioconda mitos 1.1.7, Infernal 1.0, MITOS1 reference data refseq39 of Zenodo 2683856):
          runmitos.py of $MITOS1_BIN (default ~/miniforge3/envs/mitos1/bin), -r $MITOS1_REFDIR
  mitofinder  MitoFinder -a (annotation of an assembly) with MiTFi tRNAs; the reference (-r) is the closest
          record of --ref_pool (longest shared NCBI lineage, same genus allowed, same species excluded):
          --mitofinder <MitoFinder script> [--py2 python2.7]
A record whose phylum MitoZ does not support (empty mitoz_clade) is 'unsupported' for ours and mitoz.

Per tool and record: <root>/runs/<tool>/<acc>/{final.gbf, status.json, log.txt}. final.gbf is in the
NCBI record's orientation and origin (lift_to_ncbi). Records with a status.json are skipped (resume).

Usage: python3 mito_bench_run.py --set sets/mito_bench_2026-10.tsv --root projects/mito_benchmark
         [--split dev] [--tools ours,mitoz,mitos2] [--jobs 4] [--threads 4] [--timeout 3600] [--acc A,B]
         [--recollect]
"""
import argparse
import csv
import json
import os
import re
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, ExactPosition, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
PIPE = HERE.parent
sys.path.insert(0, str(PIPE))
import mito_env  # noqa: E402
from mito_rescue import mitos_to_record  # noqa: E402,F401  (used by collect(); tests import it from here)
sys.path.insert(0, str(HERE))
from mito_names import PCG, gene_name  # noqa: E402

MIN_REF_PCG = 12                          # a MitoFinder reference must carry (nearly) all genes with standard names

TOOLS = ("ours", "mitoz", "mitos2", "mitos1", "mitom", "mitom_auto", "mitom_native", "mitom_lo", "mitofinder")
# mitom_lo = mitom_native without the record's own proteins among the core's references (refseq89m holds the proteins of
# 8,012 RefSeq records, MITOS2 uses the same set): the leave-out accuracy a new genome can expect
MITOFINDER: Dict[str, str] = {}          # {"script": ..., "py2": ...} from the command line
REF_POOL: List[dict] = []                # --ref_pool records: {acc, organism, lineage, record}
# MitoM extended clades (batch_mitoz.EXTRA_CLADES) for phyla MitoZ lacks; any other phylum -> all-animal database
EXTRA_CLADE = {"Cnidaria": "Cnidaria", "Platyhelminthes": "Platyhelminthes-flatworms",
               "Chaetognatha": "Chaetognatha-arrow-worms"}
ANCHOR = 50


class LiftError(Exception):
    pass


# ── coordinates ──────────────────────────────────────────────────────────────
def _merge(parts: List[FeatureLocation]) -> List[FeatureLocation]:
    """Join consecutive parts that touch (an origin split undone by rotation)."""
    out: List[FeatureLocation] = []
    for p in parts:
        if out and out[-1].strand == p.strand:
            q = out[-1]
            if p.strand == -1 and int(p.end) == int(q.start):
                out[-1] = FeatureLocation(p.start, q.end, -1)
                continue
            if p.strand != -1 and int(q.end) == int(p.start):
                out[-1] = FeatureLocation(q.start, p.end, p.strand)
                continue
        out.append(p)
    return out


def _shift(part: FeatureLocation, off: int, n: int) -> List[FeatureLocation]:
    """Move one part by off on a circle of length n; split at the origin (biological order kept). Partial markers
    (<start, >end) stay on their position: table2asn accepts a 5'- / 3'-partial CDS without start / stop codon."""
    s = (int(part.start) + off) % n
    e = s + int(part.end) - int(part.start)
    sp = BeforePosition if isinstance(part.start, BeforePosition) else ExactPosition
    ep = AfterPosition if isinstance(part.end, AfterPosition) else ExactPosition
    if e <= n:
        return [FeatureLocation(sp(s), ep(e), part.strand)]
    a, b = FeatureLocation(sp(s), n, part.strand), FeatureLocation(0, ep(e - n), part.strand)
    return [b, a] if part.strand == -1 else [a, b]


def lift_to_ncbi(pred: SeqRecord, ncbi_seq: Seq) -> SeqRecord:
    """Return pred re-expressed on ncbi_seq (same molecule, possibly rotated and/or reverse-complemented)."""
    n = len(ncbi_seq)
    pseq = Seq(str(pred.seq))
    if len(pseq) != n:
        raise LiftError(f"length {len(pseq)} != NCBI {n}")
    ring = str(ncbi_seq).upper() * 2
    for rc in (False, True):
        q = str(pseq.reverse_complement() if rc else pseq).upper()
        hits = [m.start() for m in re.finditer(f"(?={re.escape(q[:ANCHOR])})", ring[:n + ANCHOR - 1])]
        hits = [h for h in hits if ring[h:h + n] == q]
        if len(hits) == 1:
            off = hits[0]
            break
    else:
        raise LiftError("no unique match of the predicted sequence on the NCBI sequence")
    out = SeqRecord(ncbi_seq, id=pred.id, name=pred.name, description=pred.description)
    out.annotations = dict(pred.annotations)
    for f in pred.features:
        parts = list(f.location.parts)
        if rc:                                           # a <start becomes a >end and vice versa
            parts = [FeatureLocation(BeforePosition(n - int(p.end)) if isinstance(p.end, AfterPosition) else n - int(p.end),
                                     AfterPosition(n - int(p.start)) if isinstance(p.start, BeforePosition) else n - int(p.start),
                                     -(p.strand or 1)) for p in parts]
        new = _merge([x for p in parts for x in _shift(p, off, n)])
        loc = new[0] if len(new) == 1 else CompoundLocation(new)
        q = dict(f.qualifiers)
        if "transl_except" in q:                   # positions inside the qualifier move with the feature
            q["transl_except"] = [_lift_except(t, rc, off, n) for t in q["transl_except"]]
        out.features.append(SeqFeature(loc, type=f.type, qualifiers=q))
    return out


def _lift_except(text: str, rc: bool, off: int, n: int) -> str:
    """'(pos:X..Y,aa:TERM)' / '(pos:complement(X),aa:TERM)' re-expressed after reverse-complement and rotation."""
    m = re.match(r"\(pos:(complement\()?(\d+)(?:\.\.(\d+))?\)?,(.*)$", text)
    if not m:
        return text
    def move(x: int) -> int:
        i = n - x if rc else x - 1
        return (i + off) % n + 1
    a, b = move(int(m[2])), move(int(m[3] or m[2]))
    lo, hi = (b, a) if rc else (a, b)              # rc reverses the order; lo > hi only when the span wraps
    if lo <= hi:
        pos = f"{lo}..{hi}" if lo != hi else f"{lo}"
    else:                                          # a 2-nt stop split by the origin
        pos = f"join({lo}{'' if lo == n else f'..{n}'},1{'' if hi == 1 else f'..{hi}'})"
    return f"(pos:complement({pos}),{m[4]}" if bool(m[1]) != rc else f"(pos:{pos},{m[4]}"


# ── MITOS2 output: parser shared with the pipeline (mito_rescue.mitos_to_record) ──


# ── running ──────────────────────────────────────────────────────────────────
def _env_path(exe: Optional[str]) -> Dict[str, str]:
    env = dict(os.environ)
    if exe:
        env["PATH"] = str(Path(exe).parent) + os.pathsep + env.get("PATH", "")
    return env


def build_cmd(tool: str, row: dict, bare: Path, out: Path, threads: int) -> Optional[List[List[str]]]:
    """Command sequence for one tool and record; None when the tool does not support the record's phylum."""
    acc, code = row["accession"], str(row["genetic_code"])
    if tool in ("ours", "mitoz") and not row["mitoz_clade"]:
        return None
    if tool in ("mitom", "mitom_auto", "mitom_native", "mitom_lo"):   # one command: annotation, MitoM steps, curation
        clade = row["mitoz_clade"] or EXTRA_CLADE.get(row.get("phylum", ""), "Animal")
        return [[*(["env", f"MITOM_EXCLUDE_REF={acc}"] if tool == "mitom_lo" else []),
                 sys.executable, str(PIPE / "batch_mitoz.py"), "run", "-i", str(bare), "-o", str(out / "raw"),
                 "--clade", clade, "--genetic_code", code, "--threads", str(threads), "--max_tasks", "1",
                 *(["--if_reannotation_2", "auto"] if tool == "mitom_auto" else []),
                 *(["--core", "native"] if tool in ("mitom_native", "mitom_lo") else []),
                 # the phylum `--clade auto` reads from the taxonomy (an explicit --clade Animal has none): rRNA HMM
                 *(["--phylum", row["phylum"]] if row.get("phylum") else []),
                 *(["--keep_intermediate", "yes"] if tool == "mitom_lo" else []),   # pre-curation files: rule A/B tests
                 *(["--rrna_hmm_lib", os.environ["MITOM_RRNA_LIB"]] if os.environ.get("MITOM_RRNA_LIB") else [])]]
    if tool == "ours":
        clade = row["mitoz_clade"]
        return [[sys.executable, str(PIPE / "batch_mitoz.py"), "run", "-i", str(bare), "-o", str(out / "raw"),
                 "--clade", clade, "--genetic_code", code, "--threads", str(threads), "--max_tasks", "1",
                 "--keep_intermediate", "yes"],
                [sys.executable, str(PIPE / "mito_curate.py"), "-i", str(out / "raw" / "work" / "04.final_gb"),
                 "-o", str(out / "curated"), "--genetic_code", code, "--no_validate"]]
    fa = str(out / "in.fa")
    if tool == "mitoz":
        exe = mito_env.find("mitoz", None) or "mitoz"
        return [[exe, "annotate", "--workdir", str(out / "w"), "--outprefix", acc, "--thread_number", str(threads),
                 "--clade", row["mitoz_clade"], "--genetic_code", code, "--fastafiles", fa]]
    if tool == "mitofinder":
        ncbi = SeqIO.read(str(out.parent.parent.parent / "gb" / f"{acc}.gb"), "genbank")
        ref = pick_reference(ncbi, REF_POOL)
        if ref is None:
            raise ValueError(f"{acc}: no MitoFinder reference in --ref_pool")
        SeqIO.write(reference_for_mitofinder(ref["record"]), str(out / "ref.gb"), "genbank")
        (out / "ref.txt").write_text(f"{ref['acc']}\t{ref['organism']}\n")
        return [[MITOFINDER["py2"], MITOFINDER["script"], "-j", acc, "-a", fa, "-r", str(out / "ref.gb"),
                 "-o", code, "-p", str(threads), "-m", "10", "--max-contig-size", str(len(ncbi.seq) + 1000)]]
    if tool == "mitos2":
        (out / "mitos").mkdir(parents=True, exist_ok=True)          # runmitos refuses a missing -o directory
        exe = mito_env.find("runmitos", None) or "runmitos"
        cmd = [exe, "-i", fa, "-c", code, "-o", str(out / "mitos"), "-r", "refseq89m",
               "--refdir", str(mito_env.mitos_refdir(None, "refseq89m") or ""), "--noplots"]
        return [cmd + (["--linear"] if row.get("topology") == "linear" else [])]
    if tool == "mitos1":                                      # no --linear / --noplots in MITOS 1
        (out / "mitos").mkdir(parents=True, exist_ok=True)
        b = Path(os.environ.get("MITOS1_BIN", Path.home() / "miniforge3/envs/mitos1/bin"))
        return [[str(b / "python"), str(b / "runmitos.py"), "-i", fa, "-c", code, "-o", str(out / "mitos"),
                 "-r", os.environ.get("MITOS1_REFDIR", str(Path.home() / "mitos_refs1/mitos1-refdata"))]]
    raise ValueError(tool)


def reference_for_mitofinder(rec: SeqRecord) -> SeqRecord:
    """Copy of a reference record without the features MitoFinder rejects: it stops on any CDS whose name it does
    not know ('hypothetical protein' ORFs of cnidarian records, NC_020459.1), so only the standard protein genes,
    rRNAs and tRNAs are kept."""
    out = rec[:]
    out.annotations = dict(rec.annotations)
    dropped = {str(f.location) for f in rec.features if f.type == "CDS" and gene_name(f.qualifiers, "CDS") not in PCG}
    out.features = [f for f in rec.features
                    if not (f.type in ("CDS", "gene") and str(f.location) in dropped)]
    return out


def usable_reference(rec: SeqRecord) -> bool:
    """>= MIN_REF_PCG standard protein genes and both rRNAs (MitoFinder annotates only what the reference has;
    UNVERIFIED / unannotated records are useless)."""
    pcg = {gene_name(f.qualifiers, "CDS") for f in rec.features if f.type == "CDS"} & PCG
    rrna = {gene_name(f.qualifiers, "rRNA") for f in rec.features if f.type == "rRNA"}
    return len(pcg) >= MIN_REF_PCG and {"rrnS", "rrnL"} <= rrna


def load_ref_pool(paths: List[Path]) -> List[dict]:
    """GenBank records that may serve as MitoFinder reference (usable_reference only)."""
    pool = []
    for path in paths:
        for rec in SeqIO.parse(str(path), "genbank"):
            if usable_reference(rec):
                pool.append({"acc": rec.id, "organism": str(rec.annotations.get("organism", "")),
                             "lineage": list(rec.annotations.get("taxonomy", [])), "record": rec})
    return pool


def pad_to_molecule(pred: SeqRecord, ncbi_seq: Seq) -> SeqRecord:
    """MitoFinder trims the overlapping end of a circular contig: its sequence is a substring of the (rotated,
    possibly reverse-complemented) NCBI molecule. Return pred with the full molecule in pred's orientation and
    origin, so that lift_to_ncbi applies; unchanged when the lengths already agree."""
    if len(pred.seq) == len(ncbi_seq):
        return pred
    n = len(ncbi_seq)
    for full in (ncbi_seq, ncbi_seq.reverse_complement()):
        i = str(full + full).find(str(pred.seq))
        if 0 <= i < n:
            out = SeqRecord(full[i:] + full[:i], id=pred.id, name=pred.name, description=pred.description)
            out.annotations = dict(pred.annotations)
            out.features = list(pred.features)
            return out
    raise LiftError(f"length {len(pred.seq)} != NCBI {n} and not a substring of the molecule")


def pick_reference(ncbi: SeqRecord, pool: List[dict]) -> Optional[dict]:
    """Closest pool record to `ncbi`: longest shared lineage prefix, then same genus, RefSeq preferred; never the
    record itself, its species, or a record of the same organism name."""
    organism = str(ncbi.annotations.get("organism", ""))
    genus = organism.split(" ")[0]
    lineage = list(ncbi.annotations.get("taxonomy", []))
    base = ncbi.id.split(".")[0]

    def shared(a: List[str], b: List[str]) -> int:
        n = 0
        for x, y in zip(a, b):
            if x != y:
                break
            n += 1
        return n

    best, best_key = None, None
    for r in pool:
        if r["acc"].split(".")[0] == base or r["organism"] == organism or not r["organism"]:
            continue
        key = (shared(lineage, r["lineage"]), r["organism"].split(" ")[0] == genus,
               r["acc"].startswith("NC_"), -len(r["organism"]), r["acc"])
        if best_key is None or key > best_key:
            best, best_key = r, key
    return best


def _killpg(p: subprocess.Popen) -> None:
    """SIGTERM the process group, SIGKILL whatever is left after 5 s."""
    for sig, wait in ((signal.SIGTERM, 5), (signal.SIGKILL, 5)):
        try:
            os.killpg(p.pid, sig)
        except ProcessLookupError:
            return
        try:
            p.wait(wait)
            os.killpg(p.pid, 0)              # leader gone; children may still be in the group
        except (subprocess.TimeoutExpired, ProcessLookupError, PermissionError):
            pass
        time.sleep(0.2)


def _watch_rss(pid: int, stop: threading.Event, st: dict) -> None:
    """Peak resident memory (MB) of the process and all its descendants, sampled every second -> st["peak_rss_mb"]
    (left unset without the psutil package)."""
    try:
        import psutil
    except ImportError:
        return
    while not stop.wait(1.0):
        try:
            root = psutil.Process(pid)
            procs = [root] + root.children(recursive=True)
        except psutil.Error:
            continue
        total = 0
        for q in procs:
            try:
                total += q.memory_info().rss
            except psutil.Error:
                pass
        st["peak_rss_mb"] = max(st.get("peak_rss_mb", 0), round(total / 2 ** 20))


def run_steps(cmds: List[List[str]], out: Path, timeout: int) -> dict:
    """Run commands in order (each with its program's directory first on PATH); status dict."""
    out.mkdir(parents=True, exist_ok=True)
    st: dict = {"status": "ok", "seconds": 0.0, "step_seconds": {}, "error": ""}
    with (out / "log.txt").open("a") as log:
        for i, cmd in enumerate(cmds):
            t0 = time.time()
            log.write("$ " + " ".join(cmd) + "\n")
            log.flush()
            try:
                # own process group, so a timeout also stops the tools it started (MitoZ -> perl -> tblastn)
                p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.PIPE, text=True, errors="replace",
                                     start_new_session=True, cwd=str(out),
                                     env=_env_path(cmd[0] if os.sep in cmd[0] else None))
            except OSError as e:
                st.update(status="failed", error=f"step {i + 1}: {e}")
            else:
                stop = threading.Event()
                mem = threading.Thread(target=_watch_rss, args=(p.pid, stop, st), daemon=True)
                mem.start()
                try:
                    _, err = p.communicate(timeout=timeout)
                except subprocess.TimeoutExpired:
                    _killpg(p)
                    st.update(status="timeout", error=f"step {i + 1} exceeded {timeout} s")
                else:
                    log.write(err)
                    if p.returncode != 0:
                        st.update(status="failed", error=f"step {i + 1} exit {p.returncode}: {err.strip()[-500:]}")
                finally:
                    stop.set()
                    mem.join()
            st["step_seconds"][Path(cmd[1] if cmd[0] == sys.executable else cmd[0]).name] = round(time.time() - t0, 1)
            st["seconds"] = round(st["seconds"] + time.time() - t0, 1)
            if st["status"] != "ok":
                break
    return st


def write_status(out: Path, st: dict) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "status.json").write_text(json.dumps(st, indent=1))


def collect(tool: str, out: Path, ncbi: SeqRecord) -> SeqRecord:
    """The tool's annotation of this record, before lift-over."""
    if tool in ("mitos2", "mitos1"):
        return mitos_to_record(out / "mitos", Seq(str(ncbi.seq)), str(ncbi.id))
    pattern = {"ours": "curated/*.gbf", "mitom": "raw/final_annotation/*.gb",
               "mitom_auto": "raw/final_annotation/*.gb", "mitom_native": "raw/final_annotation/*.gb",
               "mitom_lo": "raw/final_annotation/*.gb", "mitoz": "w/**/*.gbf",
               "mitofinder": "*/*_Final_Results/*_mtDNA_contig.gb"}[tool]
    files = sorted(out.glob(pattern), key=lambda p: p.stat().st_size, reverse=True)
    if not files:
        raise LiftError(f"no GenBank output ({pattern})")
    rec = SeqIO.read(str(files[0]), "genbank")
    return pad_to_molecule(rec, Seq(str(ncbi.seq))) if tool == "mitofinder" else rec


def run_one(tool: str, row: dict, root: Path, threads: int, timeout: int, recollect: bool = False) -> dict:
    """Run one tool on one record (skipped when status.json exists). recollect: rebuild final.gbf from the
    kept tool output of an 'ok' run, or of a run that failed only in collect/lift, without running the tool
    again (after a parser fix)."""
    acc = row["accession"]
    out = root / "runs" / tool / acc
    old = json.loads((out / "status.json").read_text()) if (out / "status.json").exists() else None
    reusable = old and (old["status"] == "ok" or old["error"].startswith("collect/lift:"))
    if old and not (recollect and reusable):
        return old
    ncbi = SeqIO.read(str(root / "gb" / f"{acc}.gb"), "genbank")
    out.mkdir(parents=True, exist_ok=True)
    if old:
        st = {"status": "ok", "seconds": old["seconds"], "step_seconds": old["step_seconds"], "error": "",
              **({"peak_rss_mb": old["peak_rss_mb"]} if old.get("peak_rss_mb") else {})}
    else:
        (out / "in.fa").write_text(f">{acc}\n{ncbi.seq}\n")
        cmds = build_cmd(tool, row, root / "bare" / acc, out, threads)
        if cmds is None:
            st = {"status": "unsupported", "seconds": 0.0, "step_seconds": {}, "error": "phylum not supported by MitoZ"}
        else:
            st = run_steps(cmds, out, timeout)
    if st["status"] == "ok":
        try:
            lifted = lift_to_ncbi(collect(tool, out, ncbi), Seq(str(ncbi.seq)))
            SeqIO.write(lifted, str(out / "final.gbf"), "genbank")
        except Exception as e:                   # any collect / lift problem is this record's failure, not the run's
            st.update(status="failed", error=f"collect/lift: {type(e).__name__}: {e}")
    st = {"tool": tool, "accession": acc, **st}
    write_status(out, st)
    return st


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", required=True)
    p.add_argument("--root", required=True)
    p.add_argument("--split", default="dev", help="dev, test or all")
    p.add_argument("--tools", default=",".join(TOOLS))
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--timeout", type=int, default=3600)
    p.add_argument("--acc", help="only these accessions (comma-separated)")
    p.add_argument("--recollect", action="store_true",
                   help="rebuild final.gbf of finished runs from their kept tool output (no tool is re-run)")
    p.add_argument("--mitofinder", help="MitoFinder script (tool mitofinder)")
    p.add_argument("--py2", help="python2.7 for MitoFinder (default: first found in PATH / conda environments)")
    p.add_argument("--ref_pool", action="append", default=[],
                   help="GenBank file(s) the MitoFinder reference is chosen from (repeatable)")
    a = p.parse_args(argv)
    if "mitofinder" in a.tools.split(","):
        if not (a.mitofinder and a.ref_pool):
            p.error("tool mitofinder needs --mitofinder and --ref_pool")
        py2 = a.py2 or mito_env.find("python2.7")
        if not py2:
            p.error("python2.7 not found (--py2)")
        MITOFINDER.update(script=str(Path(a.mitofinder).resolve()), py2=py2)
        REF_POOL.extend(load_ref_pool([Path(x) for x in a.ref_pool]))
        print(f"[bench] MitoFinder reference pool: {len(REF_POOL)} records", flush=True)
    rows = list(csv.DictReader(open(a.set, encoding="utf-8"), delimiter="\t"))
    rows = [r for r in rows if a.split == "all" or r["split"] == a.split]
    if a.acc:
        rows = [r for r in rows if r["accession"] in a.acc.split(",")]
    tasks = [(t, r) for r in rows for t in a.tools.split(",")]
    root = Path(a.root).resolve()                  # the tools run with cwd = their own run directory
    with ThreadPoolExecutor(a.jobs) as ex:
        for st in ex.map(lambda tr: run_one(tr[0], tr[1], root, a.threads, a.timeout, a.recollect), tasks):
            print(f"{st['tool']:7s} {st['accession']:14s} {st['status']:11s} {st['seconds']:8.1f}s {st['error'][:100]}",
                  flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
