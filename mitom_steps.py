"""mitom_steps.py - the steps after the annotation: MITOS2 rescue (mito_rescue.py), rRNA HMM fallback
(mito_rrna_hmm.py), introns (mito_intron.py de novo, then mito_intron_ref.py site transfer), boundary curation and NCBI
names (mito_curate.py), NCBI submission files (ncbi_validate.py, table2asn), original vs MitoM (mito_compare.py),
genome maps (mito_plot.py).

One interface: run_<step>(args, gb_dir, out_root) works on the GenBank files in gb_dir (rescue, rRNA HMM and introns
edit them in place; curation writes 05.curated/), writes its reports into out_root and returns a StepResult: the line
for the step summary of `run`, the folder of the GenBank files after the step, and whether the step did its job."""
from __future__ import annotations

import argparse
import functools
import logging
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, List, Optional, Tuple, TypeVar

from Bio import SeqIO

from mito_common import reports_dir, run_tool, write_table
from mito_taxa import abbreviated_overlap, aggressive_start_genes, no_upstream_atg
from mitom_inputs import find_fq_pair_for_sample, parse_suffix_fq


def summary_line(text: str, prefix: str, default: str) -> str:
    """The last line of a step script's stderr that starts with its prefix ('[mito_curate] ...'), else default (the
    scripts print their one-line summary there)."""
    return next((ln for ln in reversed((text or "").splitlines()) if ln.startswith(prefix)), default)


@dataclass
class StepResult:
    """What one step after the annotation did. summary: its line in the step summary of `run` (terminal and run.log);
    gb_dir: the folder of the GenBank files after the step (curation: 05.curated/; None = the step's input folder);
    ok: False when the step failed or could not run although asked to (the run goes on; run.log has the reason)."""
    summary: str = ""
    gb_dir: Optional[Path] = None
    ok: bool = True


F = TypeVar("F", bound=Callable[..., Any])


def optional_step(tag: str, line: Optional[Callable[[argparse.Namespace, Path, Path], str]] = None) -> Callable[[F], F]:
    """A step after the annotation that must never abort the run: any exception is logged as '[TAG] failed: ...' and
    the step returns StepResult(line(args, gb_dir, out_root), ok=False) (line: the step's summary line)."""
    def deco(fn: F) -> F:
        @functools.wraps(fn)
        def wrapper(args: argparse.Namespace, gb_dir: Path, out_root: Path, *a: Any, **k: Any) -> StepResult:
            try:
                return fn(args, gb_dir, out_root, *a, **k)
            except Exception as exc:
                logging.error(f"[{tag}] failed: {type(exc).__name__}: {exc} - step skipped, the run continues")
                return StepResult(line(args, gb_dir, out_root) if line else "", ok=False)
        return wrapper  # type: ignore[return-value]
    return deco


def run_curate_and_submit(args: argparse.Namespace, gb_dir: Path, out_root: Path) -> StepResult:
    """Curation (run_curate) and, when the records were curated, the NCBI submission files (run_submission); the
    result is the curation's (gb_dir: 05.curated/, or the input folder with --curate no or when curation failed)."""
    cur = run_curate(args, gb_dir, out_root)
    if cur.gb_dir != gb_dir:
        run_submission(args, cur.gb_dir, out_root)
    return cur


def run_curate(args: argparse.Namespace, gb_dir: Path, out_root: Path) -> StepResult:
    """mito_curate.py on gb_dir -> 05.curated/ (boundary rules, NCBI gene names; curation_report.tsv); the folder
    stays gb_dir with --curate no or when curation failed."""
    if getattr(args, "curate", "yes") == "no":
        return StepResult("skipped (--curate no)", gb_dir=gb_dir)
    here = Path(__file__).resolve().parent
    cur = out_root / "05.curated"
    p = run_tool([sys.executable, "-W", "ignore", str(here / "mito_curate.py"), "-i", str(gb_dir), "-o", str(cur),
                  "--genetic_code", str(int(args.genetic_code)), "--no_validate",
                  "--aggressive_start_genes", aggressive_start_genes(getattr(args, "clade", None)),
                  *abbreviated_overlap(getattr(args, "clade", None)),
                  *no_upstream_atg(getattr(args, "clade", None))])
    if p.returncode:
        logging.error(f"[CURATE] failed (exit {p.returncode}): {p.stderr[-300:]} - final files left uncurated")
        return StepResult("curation failed (see run.log): the annotation before curation is used", gb_dir=gb_dir,
                          ok=False)
    logging.info(f"[CURATE] {summary_line(p.stderr, '[mito_curate]', '')}  -> {cur}")
    report = reports_dir(out_root) / "curation_report.tsv"
    if (cur / "curation_report.tsv").exists():
        shutil.move(str(cur / "curation_report.tsv"), str(report))
    return StepResult(f"{count_rows(report)} changes, each with rule and reason in reports/curation_report.tsv",
                      gb_dir=cur)


def _submission_line(out_root: Path) -> str:
    sub_sum = out_root / "06.ncbi_submission" / "submission_summary.tsv"
    n_rec = count_rows(sub_sum)
    n_err = count_rows(sub_sum, "ERROR", lambda v: v != "0")
    return (f"{n_rec} records checked: {n_rec - n_err} without NCBI ERROR, {n_err} with ERROR "
            f"(ncbi_submission/validation_messages.tsv)" if n_rec else
            ("validation produced no records (see run.log)" if sub_sum.parent.is_dir() else
             "table2asn not found: skipped"))


def run_submission(args: argparse.Namespace, gb_dir: Path, out_root: Path) -> StepResult:
    """NCBI table2asn validation of every final record and the submission files -> 06.ncbi_submission/: per record
    <id>.fsa + <id>.tbl (BankIt), <id>.sqn + <id>.val (table2asn; --sbt adds the submitter template),
    submission_summary.tsv and validation_messages.tsv."""
    here = Path(__file__).resolve().parent
    if str(here) not in sys.path:
        sys.path.insert(0, str(here))
    import mito_env
    import ncbi_validate as nv
    t2a = mito_env.find("table2asn")
    if not t2a:
        logging.info("[SUBMISSION] table2asn not found: no .sqn / validation (install.sh installs it)")
        return StepResult(_submission_line(out_root), ok=False)
    sub = out_root / "06.ncbi_submission"
    sub.mkdir(parents=True, exist_ok=True)
    summary, msgs = [], []
    for gbf in sorted(gb_dir.rglob("*.gbf")):
        for rec in SeqIO.parse(str(gbf), "genbank"):
            rows = nv.validate_record(rec, sub / "_work", t2a, getattr(args, "sbt", None))
            d = nv.record_dir(sub / "_work", rec.id)
            for ext in (".fsa", ".tbl", ".sqn", ".val"):
                if (d / f"rec{ext}").exists():
                    shutil.copy2(d / f"rec{ext}", sub / f"{d.name}{ext}")
            n = {sev: sum(r["severity"] == sev for r in rows) for sev in ("ERROR", "WARNING")}
            summary.append(f"{rec.id}\t{n['ERROR']}\t{n['WARNING']}\t"
                           + "; ".join(sorted({r["code"] for r in rows if r["severity"] == "ERROR"})))
            msgs += [f"{r['record']}\t{r['severity']}\t{r['code']}\t{r['message']}" for r in rows]
    shutil.rmtree(sub / "_work", ignore_errors=True)
    write_table(sub / "submission_summary.tsv", "record\tERROR\tWARNING\terror_codes", summary)
    write_table(sub / "validation_messages.tsv", "record\tseverity\tcode\tmessage", msgs)
    bad = [x.split("\t")[0] for x in summary if x.split("\t")[1] != "0"]
    logging.info(f"[SUBMISSION] {len(summary)} record(s) -> {sub}; NCBI ERROR in {len(bad)}: {bad[:10]}"
                 + (" - see validation_messages.tsv" if bad else ""))
    return StepResult(_submission_line(out_root))


def _rescue_line(args: argparse.Namespace, gb_dir: Path, out_root: Path) -> str:
    return f"genes added by MITOS2 rescue : {count_rows(out_root / 'reports' / 'rescue_report.tsv', 'status', 'added')}"


@optional_step("RESCUE", _rescue_line)
def run_rescue(args: argparse.Namespace, gb_dir: Path, out_root: Path) -> StepResult:
    """Add genes the core missed from a MITOS2 annotation (mito_rescue.py, --rescue yes; MITOS2 is optional), in
    place; report rescue_report.tsv."""
    if getattr(args, "rescue", "no") == "no":
        return StepResult(_rescue_line(args, gb_dir, out_root))
    if args.rescue == "auto":
        # the native core retries other loci / blastx placement itself; MITOS2 took 5-7 min and ~1 GB per record
        # lacking a gene (mostly atp8 of nematodes / flatworms, which have none) for 2 genes in 24 records
        logging.info("[RESCUE] skipped: native core (--rescue yes runs MITOS2 for records missing a gene)")
        return StepResult(_rescue_line(args, gb_dir, out_root))
    if not (args.mitos_bin and args.mitos_refdir):
        msg = "[RESCUE] skipped: MITOS2 not found (--mitos_bin / --mitos_refdir)"
        logging.error(msg) if args.rescue == "yes" else logging.info(msg)
        return StepResult(_rescue_line(args, gb_dir, out_root), ok=args.rescue != "yes")
    import mito_rescue
    rows = mito_rescue.rescue_dir(gb_dir, int(args.genetic_code), jobs=int(getattr(args, "max_tasks", 1) or 1),
                                  mitos_bin=args.mitos_bin, refdir=args.mitos_refdir)
    write_table(reports_dir(out_root) / "rescue_report.tsv", "file\trecord\tgene\tlocation\tstatus\thint",
                         ["\t".join([r[k] for k in ("file", "record", "gene", "location", "status")] + [r.get("hint", "")])
                          for r in rows])
    added = [f"{r['record']}:{r['gene']}" for r in rows if r["status"] == "added"]
    failed = [r["record"] for r in rows if r["status"] != "added"]
    logging.info(f"[RESCUE] {len(added)} gene(s) added from MITOS2 {added}" + (f"; MITOS2 failed: {failed}" if failed else ""))
    return StepResult(_rescue_line(args, gb_dir, out_root), ok=not failed)


def _rrna_hmm_line(args: argparse.Namespace, gb_dir: Path, out_root: Path) -> str:
    return f"rRNAs added by rRNA HMM      : {count_rows(out_root / 'reports' / 'rrna_hmm_report.tsv', 'status', 'added')}"


@optional_step("RRNA_HMM", _rrna_hmm_line)
def run_rrna_hmm(args: argparse.Namespace, gb_dir: Path, out_root: Path) -> StepResult:
    """Last-resort 12S/16S search with per-phylum profile HMMs for rRNAs still missing after the cores and
    the MITOS2 rescue (mito_rrna_hmm.py), in place; report rrna_hmm_report.tsv."""
    mode = getattr(args, "rrna_hmm", "no")
    if mode == "no":
        return StepResult(_rrna_hmm_line(args, gb_dir, out_root))
    import mito_env
    import mito_rrna_hmm
    lib = Path(args.rrna_hmm_lib) if getattr(args, "rrna_hmm_lib", None) else mito_rrna_hmm.LIB
    if not (mito_rrna_hmm.read_manifest(lib) and mito_env.find("nhmmer")):
        msg = f"[RRNA_HMM] skipped: no HMM library in {lib} or nhmmer not found"
        logging.error(msg) if mode == "yes" else logging.info(msg)
        return StepResult(_rrna_hmm_line(args, gb_dir, out_root), ok=mode != "yes")
    rows = mito_rrna_hmm.search_dir(gb_dir, lib, getattr(args, "phylum", None))
    write_table(reports_dir(out_root) / "rrna_hmm_report.tsv",
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
    return StepResult(_rrna_hmm_line(args, gb_dir, out_root),
                      ok=not any(str(r["status"]).startswith("failed") for r in rows))


INTRON_LABEL = "introns                      : "


def run_intron_fix(args: argparse.Namespace, gb_dir: Path, out_root: Path) -> StepResult:
    """Add introns to the GenBank files in place with mito_intron.py (de novo; protein library MFannot's, or
    --intron_lib), then mito_intron_ref.py (intron sites of annotated relatives, run_intron_ref); reports
    intron_report.tsv, intron_ref_report.tsv and the programs' output intron_fix.log."""
    if args.intron_fix == "no":
        return StepResult(INTRON_LABEL + "skipped (--intron_fix no)")
    script = Path(__file__).resolve().with_name("mito_intron.py")
    exe = (str(Path(args.intron_bin) / "exonerate") if args.intron_bin else shutil.which("exonerate")) or ""
    if args.intron_bin and not (Path(args.intron_bin) / "blastp").exists():
        exe = ""
    if not script.exists() or not Path(exe).exists():          # on the console: introns would be missing silently
        msg = "[INTRON_FIX] skipped: exonerate / blastp not found (install.sh installs them; or --intron_bin)"
        logging.error(msg) if args.intron_fix == "yes" else logging.warning(msg)
        return StepResult(INTRON_LABEL + "skipped: exonerate / blastp not found (--intron_bin)",
                          ok=args.intron_fix != "yes")
    cmd = [sys.executable, str(script), "-i", str(gb_dir), "-o", str(gb_dir),
           "--genetic_code", str(int(args.genetic_code)), "--threads", str(int(args.threads)),
           "--jobs", str(int(getattr(args, "max_tasks", 1) or 1))]
    for lib in args.intron_lib:
        cmd += ["--lib", lib]
    for ref in args.intron_ref_gb:
        cmd += ["--ref_gb", ref]
    if args.intron_bin:
        cmd += ["--bin", args.intron_bin]
    proc = run_tool(cmd)
    report = gb_dir / "intron_report.tsv"
    if report.exists():
        shutil.move(str(report), str(reports_dir(out_root) / "intron_report.tsv"))
    last = summary_line(proc.stderr, "[mito_intron]", f"failed (exit {proc.returncode}), see intron_fix.log")
    logging.info(f"[INTRON_FIX] {last}")
    ref, ref_log, ref_ok = run_intron_ref(args, gb_dir, out_root)
    (reports_dir(out_root) / "intron_fix.log").write_text(proc.stdout + proc.stderr + ref_log, encoding="utf-8")
    return StepResult(INTRON_LABEL + last.replace("[mito_intron] ", "").split(" -> ")[0] + f"; {ref}; per intron: "
                      "intron_report.tsv, intron_ref_report.tsv", ok=proc.returncode == 0 and ref_ok)


def run_intron_ref(args: argparse.Namespace, final_dir: Path, out_root: Path) -> Tuple[str, str, bool]:
    """mito_intron_ref.py on the de novo output, in place: the intron sites of the library models/intron_sites.*
    (193 NCBI records, verified gold standard 2026-10-29) re-place de novo boundaries and add introns de novo misses
    (mini-exons, trans-spliced pieces, rRNA introns). --phylum gives the lineage of records without a taxonomy.
    Returns (summary, log text, ok)."""
    script = Path(__file__).resolve().with_name("mito_intron_ref.py")
    cmd = [sys.executable, str(script), "annotate", "-i", str(final_dir), "-o", str(final_dir),
           "--genetic_code", str(int(args.genetic_code)), "--jobs", str(int(getattr(args, "max_tasks", 1) or 1))]
    if args.intron_bin:
        cmd += ["--bin", args.intron_bin]
    if getattr(args, "phylum", None):
        cmd += ["--lineage", f"Eukaryota;Metazoa;{args.phylum}"]
    proc = run_tool(cmd)
    report = final_dir / "intron_ref_report.tsv"
    if report.exists():
        shutil.move(str(report), str(reports_dir(out_root) / "intron_ref_report.tsv"))
    last = summary_line(proc.stderr, "[mito_intron_ref]", f"failed (exit {proc.returncode}), see intron_fix.log")
    if proc.returncode:
        logging.warning(f"[INTRON_REF] {last}")
        return f"site transfer failed (exit {proc.returncode}), see intron_fix.log", proc.stdout + proc.stderr, False
    logging.info(f"[INTRON_REF] {last}")
    return ("site transfer: " + last.replace("[mito_intron_ref] ", "").split(" -> ")[0], proc.stdout + proc.stderr,
            True)


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


def run_compare(args: argparse.Namespace, gb_dir: Path, out_root: Path, original: Optional[Path] = None) -> StepResult:
    """mito_compare.py: the original annotation of the input records (`original`: original_annotation/, GenBank input
    only) vs MitoM's, gene by gene -> original_vs_mitom.tsv. A report: a failure never ends the run."""
    if original is None or not original.is_dir():
        return StepResult("FASTA input has no original annotation: skipped")
    try:
        import mito_compare
        n = mito_compare.compare_dirs(original, gb_dir, out_root / "original_vs_mitom.tsv", int(args.genetic_code))
        return StepResult(mito_compare.summary_text(n) + "; one row per gene: original_vs_mitom.tsv")
    except Exception as exc:                   # noqa: BLE001  (a report: must not end the run after annotation)
        logging.error(f"[COMPARE] failed: {exc!r}")
        return StepResult("comparison failed (see run.log); the annotation itself is complete", ok=False)


def run_plots(args: argparse.Namespace, gb_dir: Path, out_root: Path) -> StepResult:
    """mito_plot.py: one circular map per final record -> out_root/plots/<record>.png / .pdf (--max_tasks processes,
    ~1.7 s per record); with --infq the reads of the record's sample add a sequencing-depth ring (minimap2 when
    installed, else bwa mem; --threads each). Never raises."""
    if getattr(args, "plot", "yes") == "no":
        return StepResult("skipped (--plot no)")
    try:
        import mito_plot  # noqa: F401  (matplotlib)
    except Exception as exc:                               # noqa: BLE001  (a figure must not end the run)
        logging.error(f"[PLOT] mito_plot.py not usable: {exc!r}")
        return StepResult("skipped: mito_plot.py not usable (see run.log)", ok=False)
    infq = Path(args.infq).resolve() if getattr(args, "infq", None) else None
    s1, s2 = parse_suffix_fq(args.suffix_fq) if infq else ("", "")
    jobs = [(gbf, out_root / "plots",
             [str(x) for x in find_fq_pair_for_sample(gbf.stem, infq, args.fq_position, s1, s2) if x] if infq else [],
             int(args.threads)) for gbf in sorted(gb_dir.rglob("*.gbf"))]
    n = depth_n = 0
    ok = True
    try:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(max(1, min(len(jobs), int(getattr(args, "max_tasks", 1) or 1)))) as pool:
            for k, d, msgs in pool.imap(_plot_file, jobs):          # matplotlib is not thread-safe: processes
                n, depth_n = n + k, depth_n + d
                for m in msgs:
                    logging.warning(m)
    except Exception as exc:                               # noqa: BLE001
        logging.error(f"[PLOT] failed: {exc!r}")
        ok = False
    return StepResult(f"{n} map(s) in plots/" + (f", {depth_n} with read depth" if depth_n else ""), ok=ok)
