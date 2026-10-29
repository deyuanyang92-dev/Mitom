"""mitom_steps.py - the steps after the annotation, each optional and none allowed to end the run: MITOS2 rescue
(mito_rescue.py), rRNA HMM fallback (mito_rrna_hmm.py), introns (mito_intron.py de novo, then mito_intron_ref.py site
transfer), boundary curation and NCBI names (mito_curate.py), NCBI submission files (ncbi_validate.py, table2asn),
genome maps (mito_plot.py)."""
from __future__ import annotations

import argparse
import functools
import logging
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable, List, Optional, Tuple, TypeVar

from Bio import SeqIO

from mito_common import write_nonempty_table
from mito_taxa import abbreviated_overlap, aggressive_start_genes, no_upstream_atg
from mitom_inputs import find_fq_pair_for_sample, parse_suffix_fq


def run_curate_and_submit(args: argparse.Namespace, final_dir: Path, out_root: Path) -> Path:
    """MitoM final steps: mito_curate.py (boundary rules, NCBI gene names) -> 05.curated/, then the NCBI
    submission inputs per record in 06.ncbi_submission/: <id>.fsa + <id>.tbl (BankIt), <id>.sqn + <id>.val
    (table2asn; --sbt adds the submitter template) and submission_summary.tsv. Returns the folder of the
    final GenBank files (05.curated, or final_dir with --curate no)."""
    cur = run_curate(args, final_dir, out_root)
    if cur != final_dir:
        run_submission(args, cur, out_root)
    return cur


def run_curate(args: argparse.Namespace, final_dir: Path, out_root: Path) -> Path:
    """mito_curate.py on final_dir -> 05.curated/ (returned); final_dir when --curate no or curation failed."""
    if getattr(args, "curate", "yes") == "no":
        return final_dir
    here = Path(__file__).resolve().parent
    cur = out_root / "05.curated"
    p = subprocess.run([sys.executable, "-W", "ignore", str(here / "mito_curate.py"), "-i", str(final_dir), "-o", str(cur),
                        "--genetic_code", str(int(args.genetic_code)), "--no_validate",
                        "--aggressive_start_genes", aggressive_start_genes(getattr(args, "clade", None)),
                        *abbreviated_overlap(getattr(args, "clade", None)),
                        *no_upstream_atg(getattr(args, "clade", None))],
                       capture_output=True, text=True)
    if p.returncode:
        logging.error(f"[CURATE] failed (exit {p.returncode}): {p.stderr[-300:]} - final files left uncurated")
        return final_dir
    logging.info(f"[CURATE] {(p.stdout.strip().splitlines() or [''])[-1]}  -> {cur}")
    return cur


def run_submission(args: argparse.Namespace, cur: Path, out_root: Path) -> None:
    """NCBI table2asn validation of every final record and the submission files -> 06.ncbi_submission/."""
    here = Path(__file__).resolve().parent
    if str(here) not in sys.path:
        sys.path.insert(0, str(here))
    import mito_env
    import ncbi_validate as nv
    t2a = mito_env.find("table2asn")
    if not t2a:
        logging.info("[SUBMISSION] table2asn not found: no .sqn / validation (install.sh installs it)")
        return
    sub = out_root / "06.ncbi_submission"
    sub.mkdir(parents=True, exist_ok=True)
    summary, msgs = [], []
    for gbf in sorted(cur.rglob("*.gbf")):
        for rec in SeqIO.parse(str(gbf), "genbank"):
            rows = nv.validate_record(rec, sub / "_work", t2a, getattr(args, "sbt", None))
            d = sub / "_work" / re.sub(r"[^\w.]+", "_", rec.id)
            for ext in (".fsa", ".tbl", ".sqn", ".val"):
                if (d / f"rec{ext}").exists():
                    shutil.copy2(d / f"rec{ext}", sub / f"{d.name}{ext}")
            n = {sev: sum(r["severity"] == sev for r in rows) for sev in ("ERROR", "WARNING")}
            summary.append(f"{rec.id}\t{n['ERROR']}\t{n['WARNING']}\t"
                           + "; ".join(sorted({r["code"] for r in rows if r["severity"] == "ERROR"})))
            msgs += [f"{r['record']}\t{r['severity']}\t{r['code']}\t{r['message']}" for r in rows]
    shutil.rmtree(sub / "_work", ignore_errors=True)
    write_nonempty_table(sub / "submission_summary.tsv", "record\tERROR\tWARNING\terror_codes", summary)
    write_nonempty_table(sub / "validation_messages.tsv", "record\tseverity\tcode\tmessage", msgs)
    bad = [x.split("\t")[0] for x in summary if x.split("\t")[1] != "0"]
    logging.info(f"[SUBMISSION] {len(summary)} record(s) -> {sub}; NCBI ERROR in {len(bad)}: {bad[:10]}"
                 + (" - see validation_messages.tsv" if bad else ""))


F = TypeVar("F", bound=Callable[..., Any])


def optional_step(tag: str) -> Callable[[F], F]:
    """An optional post-annotation step that must never abort the run: any exception is logged as
    '[TAG] failed: ...' and the workflow goes on (intron fix, curation, submission files, merge)."""
    def deco(fn: F) -> F:
        @functools.wraps(fn)
        def wrapper(*a: Any, **k: Any) -> Any:
            try:
                return fn(*a, **k)
            except Exception as exc:
                logging.error(f"[{tag}] failed: {type(exc).__name__}: {exc} - step skipped, the run continues")
                return None
        return wrapper  # type: ignore[return-value]
    return deco


@optional_step("RESCUE")
def run_rescue(args: argparse.Namespace, final_dir: Path, out_root: Path) -> None:
    """Add genes the core missed from a MITOS2 annotation (mito_rescue.py, --rescue yes; MITOS2 is optional), in
    place; report rescue_report.tsv."""
    if getattr(args, "rescue", "no") == "no":
        return
    if args.rescue == "auto":
        # the native core retries other loci / blastx placement itself; MITOS2 took 5-7 min and ~1 GB per record
        # lacking a gene (mostly atp8 of nematodes / flatworms, which have none) for 2 genes in 24 records
        logging.info("[RESCUE] skipped: native core (--rescue yes runs MITOS2 for records missing a gene)")
        return
    if not (args.mitos_bin and args.mitos_refdir):
        msg = "[RESCUE] skipped: MITOS2 not found (--mitos_bin / --mitos_refdir)"
        logging.error(msg) if args.rescue == "yes" else logging.info(msg)
        return
    import mito_rescue
    rows = mito_rescue.rescue_dir(final_dir, int(args.genetic_code), jobs=int(getattr(args, "max_tasks", 1) or 1),
                                  mitos_bin=args.mitos_bin, refdir=args.mitos_refdir)
    write_nonempty_table(out_root / "rescue_report.tsv", "file\trecord\tgene\tlocation\tstatus\thint",
                         ["\t".join([r[k] for k in ("file", "record", "gene", "location", "status")] + [r.get("hint", "")])
                          for r in rows])
    added = [f"{r['record']}:{r['gene']}" for r in rows if r["status"] == "added"]
    failed = [r["record"] for r in rows if r["status"] != "added"]
    logging.info(f"[RESCUE] {len(added)} gene(s) added from MITOS2 {added}" + (f"; MITOS2 failed: {failed}" if failed else ""))


@optional_step("RRNA_HMM")
def run_rrna_hmm(args: argparse.Namespace, final_dir: Path, out_root: Path) -> None:
    """Last-resort 12S/16S search with per-phylum profile HMMs for rRNAs still missing after the cores and
    the MITOS2 rescue (mito_rrna_hmm.py), in place; report rrna_hmm_report.tsv."""
    mode = getattr(args, "rrna_hmm", "no")
    if mode == "no":
        return
    import mito_env
    import mito_rrna_hmm
    lib = Path(args.rrna_hmm_lib) if getattr(args, "rrna_hmm_lib", None) else mito_rrna_hmm.LIB
    if not (mito_rrna_hmm.read_manifest(lib) and mito_env.find("nhmmer")):
        msg = f"[RRNA_HMM] skipped: no HMM library in {lib} or nhmmer not found"
        logging.error(msg) if mode == "yes" else logging.info(msg)
        return
    rows = mito_rrna_hmm.search_dir(final_dir, lib, getattr(args, "phylum", None))
    write_nonempty_table(out_root / "rrna_hmm_report.tsv",
                         "file\trecord\tphylum\tgene\tlocation\tstatus\tcoverage\tnote",
                         ["\t".join([r[k] for k in ("file", "record", "phylum", "gene", "location", "status")]
                                    + [str(r.get(k, "")) for k in ("coverage", "note")]) for r in rows])
    added = [f"{r['record']}:{r['gene']}" for r in rows if str(r["status"]).startswith("added")]
    logging.info(f"[RRNA_HMM] {len(added)} rRNA(s) added from per-phylum HMMs {added}")
    n_skip = sum(str(r["status"]).startswith(("skipped: no library", "skipped: no model")) for r in rows)
    if n_skip:
        msg = (f"[RRNA_HMM] WARNING: {n_skip} missing rRNA(s) were not searched because no HMM model was chosen for "
               f"their phylum (phylum unknown, or no model for it); give the phylum with --phylum, e.g. --phylum Mollusca "
               f"(see rrna_hmm_report.tsv)")
        logging.warning(msg)
        print(msg, file=sys.stderr)


def run_intron_fix(args: argparse.Namespace, final_dir: Path, out_root: Path) -> str:
    """Add introns to final GenBank files in place with mito_intron.py (de novo; protein library MFannot's, or
    --intron_lib), then mito_intron_ref.py (intron sites of annotated relatives, run_intron_ref); returns the line for
    the step summary."""
    if args.intron_fix == "no":
        return "skipped (--intron_fix no)"
    script = Path(__file__).resolve().with_name("mito_intron.py")
    exe = (str(Path(args.intron_bin) / "exonerate") if args.intron_bin else shutil.which("exonerate")) or ""
    if args.intron_bin and not (Path(args.intron_bin) / "blastp").exists():
        exe = ""
    if not script.exists() or not Path(exe).exists():          # on the console: introns would be missing silently
        msg = "[INTRON_FIX] skipped: exonerate / blastp not found (install.sh installs them; or --intron_bin)"
        logging.error(msg) if args.intron_fix == "yes" else logging.warning(msg)
        return "skipped: exonerate / blastp not found (--intron_bin)"
    cmd = [sys.executable, str(script), "-i", str(final_dir), "-o", str(final_dir),
           "--genetic_code", str(int(args.genetic_code)), "--threads", str(int(args.threads)),
           "--jobs", str(int(getattr(args, "max_tasks", 1) or 1))]
    for lib in args.intron_lib:
        cmd += ["--lib", lib]
    for ref in args.intron_ref_gb:
        cmd += ["--ref_gb", ref]
    if args.intron_bin:
        cmd += ["--bin", args.intron_bin]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    report = final_dir / "intron_report.tsv"
    if report.exists():
        shutil.move(str(report), str(out_root / "intron_report.tsv"))
    last = (proc.stdout.strip().splitlines() or [f"failed (exit {proc.returncode}), see intron_fix.log"])[-1]
    logging.info(f"[INTRON_FIX] {last}")
    ref, ref_log = run_intron_ref(args, final_dir, out_root)
    (out_root / "intron_fix.log").write_text(proc.stdout + proc.stderr + ref_log, encoding="utf-8")
    return (last.replace("[mito_intron] ", "").split(" -> ")[0] + f"; {ref}; per intron: intron_report.tsv, "
            "intron_ref_report.tsv")


def run_intron_ref(args: argparse.Namespace, final_dir: Path, out_root: Path) -> Tuple[str, str]:
    """mito_intron_ref.py on the de novo output, in place: the intron sites of the library models/intron_sites.*
    (193 NCBI records, verified gold standard 2026-10-29) re-place de novo boundaries and add introns de novo misses
    (mini-exons, trans-spliced pieces, rRNA introns). --phylum gives the lineage of records without a taxonomy.
    Returns (summary, log text)."""
    script = Path(__file__).resolve().with_name("mito_intron_ref.py")
    cmd = [sys.executable, str(script), "annotate", "-i", str(final_dir), "-o", str(final_dir),
           "--genetic_code", str(int(args.genetic_code)), "--jobs", str(int(getattr(args, "max_tasks", 1) or 1))]
    if args.intron_bin:
        cmd += ["--bin", args.intron_bin]
    if getattr(args, "phylum", None):
        cmd += ["--lineage", f"Eukaryota;Metazoa;{args.phylum}"]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    report = final_dir / "intron_ref_report.tsv"
    if report.exists():
        shutil.move(str(report), str(out_root / "intron_ref_report.tsv"))
    last = (proc.stderr.strip().splitlines() or [f"failed (exit {proc.returncode}), see intron_fix.log"])[-1]
    if proc.returncode:
        logging.warning(f"[INTRON_REF] {last}")
        return f"site transfer failed (exit {proc.returncode}), see intron_fix.log", proc.stdout + proc.stderr
    logging.info(f"[INTRON_REF] {last}")
    return "site transfer: " + last.replace("[mito_intron_ref] ", "").split(" -> ")[0], proc.stdout + proc.stderr


def count_rows(table: Path, column: Optional[str] = None, value: Any = None) -> int:
    """Data rows of a TSV report (0 when missing); with column, rows where it equals value / startswith value /
    value(cell) is true."""
    if not table.is_file():
        return 0
    lines = [x for x in table.read_text(encoding="utf-8", errors="ignore").splitlines() if x and not x.startswith("#")]
    if len(lines) < 2:
        return 0
    if column is None:
        return len(lines) - 1
    head = lines[0].split("\t")
    if column not in head:
        return 0
    i = head.index(column)
    cells = [r.split("\t")[i] if len(r.split("\t")) > i else "" for r in lines[1:]]
    if callable(value):
        return sum(bool(value(c)) for c in cells)
    return sum(c == value or str(c).startswith(str(value)) for c in cells)


def _plot_file(job: Tuple[Path, Path, List[str], int]) -> Tuple[int, int, List[str]]:
    """run_plots worker, one final GenBank file: (maps written, maps with read depth, messages for the log)."""
    import mito_plot
    gbf, out, reads, threads = job
    n = depth_n = 0
    msgs: List[str] = []
    for rec in SeqIO.parse(str(gbf), "genbank"):
        depth = None
        if reads:
            try:
                depth = mito_plot.depth_from_reads(rec, reads, threads=threads)
                depth_n += 1
            except Exception as exc:                       # noqa: BLE001  (a figure must not end the run)
                msgs.append(f"[PLOT] {rec.id}: no depth ring ({exc})")
        try:
            mito_plot.plot_record(rec, out, depth=depth)
            n += 1
        except Exception as exc:                           # noqa: BLE001
            msgs.append(f"[PLOT] {rec.id}: {exc!r}")
    return n, depth_n, msgs


def run_plots(args: argparse.Namespace, gb_dir: Path, out_root: Path) -> str:
    """mito_plot.py: one circular map per final record -> out_root/plots/<record>.png / .pdf (--max_tasks processes,
    ~1.7 s per record); with --infq the reads of the record's sample add a sequencing-depth ring (minimap2 when installed, else bwa mem; --threads
    each). Returns the line for the step summary; never raises."""
    if getattr(args, "plot", "yes") == "no":
        return "skipped (--plot no)"
    try:
        import mito_plot  # noqa: F401  (matplotlib)
    except Exception as exc:                               # noqa: BLE001  (a figure must not end the run)
        logging.error(f"[PLOT] mito_plot.py not usable: {exc!r}")
        return "skipped: mito_plot.py not usable (see run.log)"
    infq = Path(args.infq).resolve() if getattr(args, "infq", None) else None
    s1, s2 = parse_suffix_fq(args.suffix_fq) if infq else ("", "")
    jobs = [(gbf, out_root / "plots",
             [str(x) for x in find_fq_pair_for_sample(gbf.stem, infq, args.fq_position, s1, s2) if x] if infq else [],
             int(args.threads)) for gbf in sorted(gb_dir.rglob("*.gbf"))]
    n = depth_n = 0
    try:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(max(1, min(len(jobs), int(getattr(args, "max_tasks", 1) or 1)))) as pool:
            for k, d, msgs in pool.imap(_plot_file, jobs):          # matplotlib is not thread-safe: processes
                n, depth_n = n + k, depth_n + d
                for m in msgs:
                    logging.warning(m)
    except Exception as exc:                               # noqa: BLE001
        logging.error(f"[PLOT] failed: {exc!r}")
    return f"{n} map(s) in plots/" + (f", {depth_n} with read depth" if depth_n else "")
