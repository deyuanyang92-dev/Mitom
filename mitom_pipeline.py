#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""mitom_pipeline.py - the MitoM workflow: `mitom run`, `import`, `annotate`, `reorient`, `transfer_metadata`,
`replace_organism`, `dedup` and `intron` (`mitom <command> -h` lists the options; batch_mitoz.py is the old name).

run      GenBank input: records deduplicated (RefSeq NC_ kept), the original annotation saved for the comparison;
         FASTA input: one sample per record -> round-1 annotation (MitoM's core mito_core.py, or --engine mitoz /
         mitos2 via mito_engines.py) -> re-orientation (plus strand, cox1 first) -> round 2 only for the records
         whose annotation the shift breaks (--if_reannotation_2 auto) -> final GenBank with the original header and
         /organism -> MITOS2 rescue (optional), rRNA HMM fallback, introns (mito_intron.py de novo, then
         mito_intron_ref.py site transfer) -> boundary rules and NCBI names (mito_curate.py) -> NCBI submission files
         (ncbi_validate.py, table2asn) -> original vs MitoM (mito_compare.py) -> genome maps (mito_plot.py) -> the
         output layout of OUTPUT_LAYOUT (README.txt in the output folder).
import   the same steps on existing MitoZ / MITOS2 results or any annotated GenBank, without annotating again.
annotate one round of annotation and the steps after it, in the work layout (04.final_gb/, 05.curated/ ...).

The annotation runs --jobs samples at a time with --threads each, in threads or with Snakemake
(--scheduler_backend snakemake / auto; Snakemake is optional and not installed with MitoM); everything after the
annotation runs in this process.

This file holds the workflows, the command line and the environment preflight; their parts are:
  mitom_inputs.py    input type, --clade auto, GenBank deduplication and conversion, samples and internal IDs
  mitom_stage.py     the annotation stage (engine commands, one job per sample; mitom_snakemake.py: Snakemake)
  mitom_orient.py    re-orientation, the round-2 check, rotation after the intron step
  mitom_metadata.py  original GenBank header and /organism restored
  mitom_output.py    final GenBank per sample, summaries, the output layout
  mitom_steps.py     rescue, rRNA HMM fallback, introns, curation, NCBI submission files, genome maps
  mitom_log.py       run log, numbered steps and progress on the terminal
  mito_taxa.py       clades, phyla, genetic codes; mito_common.py: shared helpers
Every name of these modules stays available as mitom_pipeline.<name> (scripts that import mitom_pipeline or
batch_mitoz keep working).

Biopython: Cock et al. (2009) Bioinformatics doi:10.1093/bioinformatics/btp163
Snakemake: Mölder et al. (2021) F1000Res doi:10.12688/f1000research.29032.2
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

try:
    import Bio  # noqa: F401
except ImportError:                       # `mitom` checks this before it starts the script
    if __name__ == "__main__":
        print("ERROR: Biopython not installed.  pip install biopython", file=sys.stderr)
        raise SystemExit(2)
    raise
from Bio import SeqIO  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
# the parts of the workflow (see above); every name is re-exported, so mitom_pipeline.<name> keeps working
from mito_common import (command_line, debug_mode, GB_SUFFIXES, gb_fasta_record, iter_gb_files, old_spelling,  # noqa: E402,F401
                         reports_dir, run_tool,  # noqa: F401
                         safe_copy, safe_name, write_nonempty_list, write_nonempty_table, write_table)  # noqa: F401
from mito_taxa import (abbreviated_overlap, aggressive_start_genes, ALL_CLADES, BASE_CLADES,  # noqa: E402,F401
                       CLADE_CODES, CLADE_PHYLUM, EXTRA_CLADES, METAZOAN_PHYLA, no_upstream_atg,  # noqa: F401
                       PHYLUM_CLADE, PHYLUM_CODE, PHYLUM_SYNONYMS)  # noqa: F401
from mitom_log import (close_sample_logger, duration, get_sample_logger, _LockedHandler, _LOG_LOCK,  # noqa: E402,F401
                       Progress, record_event, say, setup_root_logger, Steps, tail_text, watch_files)  # noqa: F401
from mitom_metadata import (_acc_variants, build_metadata_index, build_organism_index, _clean_acc,  # noqa: E402,F401
                            diagnose_unmatched, find_metadata_from_index,  # noqa: F401
                            gb_file_has_meaningful_metadata, gb_record_accession,  # noqa: F401
                            _get_record_accession, interactive_resolve, _LOCUS_BPCOUNT_RE,  # noqa: F401
                            lookup_organism, _meta_is_meaningful, _NON_REFSEQ_VERSION_PAT,  # noqa: F401
                            _normalize_acc, _REFSEQ_PREFIX_PAT, replace_metadata_in_text,  # noqa: F401
                            replace_organism_in_gb_file, replace_organism_in_record, _similarity,  # noqa: F401
                            _split_gb_records, _split_meta_body, _sync_organism_qualifier,  # noqa: F401
                            transfer_gb_metadata_file)  # noqa: F401
from mitom_inputs import (auto_clade, check_transl_tables, collect_genbank_files, collect_sequence_files, record_format,  # noqa: E402,F401
                          convert_gb_to_fasta_for_annotation, dedup_gb_records, _FASTA_EXTS,  # noqa: F401
                          find_fq_pair_for_sample, _GB_EXTS, _gb_version_int, GbRecordMeta,  # noqa: F401
                          header_first_token, infer_topology_from_description, _is_fasta_file,  # noqa: F401
                          _is_gb_file, is_refseq_acc, _is_topology_tag, iter_inputs, LAST_DEDUP_STATS,  # noqa: F401
                          make_internal_id, MIN_CIRCULAR_LEN, parse_suffix_fq, REFSEQ_ACC_RE,  # noqa: F401
                          refseq_source_accessions, REFSEQ_SOURCE_RE, SampleInfo, _sniff_input_type,  # noqa: F401
                          split_or_sanitize_fastas, write_dedup_report, write_one_record_fasta)  # noqa: F401
from mitom_orient import (edge_gene_truncated, feature_gene_names, feature_start_pos,  # noqa: E402,F401
                          find_orientation_feature, gb_to_internal_fasta, _GENE_ALIAS_GROUPS,  # noqa: F401
                          gene_aliases, _LINEAR_SKIP_SEEN, merge_abutting_parts, normalize_gene_key,  # noqa: F401
                          ORIENTATION_FALLBACK, orientation_genes, orientation_text, _rc_location,  # noqa: F401
                          REANNOT_EDGE, reorient_after_introns, reorient_check, reorient_one_gb,  # noqa: F401
                          reverse_complement_record, rotate_record, RRNA_MIN, select_reannotation,  # noqa: F401
                          shift_location, STANDARD_GENES, starts_with_gene, TRNA_LEN,  # noqa: F401
                          write_orientation_fallback)  # noqa: F401
from mitom_snakemake import (AnnotTask, int_mem_mb, render_snakemake_snakefile,  # noqa: E402,F401
                             render_snakemake_worker_script, run_snakemake_for_annotation_tasks,  # noqa: F401
                             snakemake_available, which_or_none)  # noqa: F401
from mitom_stage import (annotate_one, engine_cmd, ENGINE_LABEL, MITO_CORE, MITO_ENGINES,  # noqa: E402,F401
                         native_core_cmd, run_annotation_stage, run_cmd_capture, run_phase2,  # noqa: F401
                         sample_lineage, STAGE_LABEL, write_input_annotation)  # noqa: F401
from mitom_output import (cleanup_intermediates, counts_from_genbank, finalize_genbank_output, new_reports,  # noqa: E402,F401
                          INTERMEDIATE_DIRS, _LOCUS_NAME_MAXLEN, merge_gb_files, organize_outputs,  # noqa: F401
                          OUTPUT_LAYOUT, refresh_summary, REPORT_FILES, REPORT_NOTES,  # noqa: F401
                          rewrite_genbank_record_ids, RunSummary, WORK_DIRS, write_summary)  # noqa: F401
from mitom_topology import apply_topology, linear_lineages_first  # noqa: E402
from mitom_steps import (count_rows, intron_denovo_cmd, INTRON_LABEL, intron_ref_cmd, optional_step,  # noqa: E402,F401
                         _plot_file, run_compare, run_curate, run_curate_and_submit, run_intron_fix,  # noqa: F401
                         run_intron_ref, run_plots, run_rescue, run_rrna_hmm, run_submission, StepResult)  # noqa: F401


# ── environment preflight ────────────────────────────────────────────────────────────────────────────────────────────
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


# ── before a run: clade, programs, imported results ──────────────────────────────────────────────────────────────────
def announce(msg: str) -> None:
    """A decision made before the run's log file exists (resolve_tools): log it and show it on the console."""
    logging.info(msg)
    print(msg, file=sys.stderr)


def resolve_tools(args: argparse.Namespace) -> str:
    """Fill unset tool locations from mito_env (config file, PATH, conda environments) and check the clade.
    Returns an error message, '' when fine."""
    if str(Path(__file__).resolve().parent) not in sys.path:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
    import mito_env
    if args.clade == "auto":
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
    SeqIO.write(recs, str(out / "imported.gb"), "genbank")
    for note in notes:
        announce(f"[IMPORT] {note}")
    announce(f"[IMPORT] {len(recs)} annotated record(s) from {args.input}")
    args.input, args.engine = str(out / "imported.gb"), "input"
    return ""


# ── workflows ────────────────────────────────────────────────────────────────────────────────────────────────────────
KNOWN_LIMITATION = ("[KNOWN_LIMITATION] the annotation core does not annotate introns (group I/II); they are added "
                    "afterwards by mito_intron.py (MFannot method) and mito_intron_ref.py (intron sites of "
                    "annotated relatives; --intron_fix) when exonerate is available - check intron_report.tsv "
                    "and intron_ref_report.tsv.")


def workflow_annotate(args: argparse.Namespace) -> int:
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    new_reports(out_root)
    setup_root_logger(out_root / "pipeline.log", debug=getattr(args, "debug", False))
    log_environment(args)
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))
    preflight(args)

    in_path     = Path(args.input).resolve()
    split_dir   = out_root / "00.internal_fasta"
    gb_conv_dir = out_root / "00.gb_converted_fasta"

    in_files, input_type, acc_to_orig_gb = iter_inputs(in_path, gb_to_fasta_dir=gb_conv_dir,
                                                       dedup_report=reports_dir(out_root) / "dedup_report.tsv")
    logging.info(f"[INPUT_TYPE] {input_type}  ({len(in_files)} file(s))")
    samples, id_map = split_or_sanitize_fastas(
        in_files, split_dir, args.if_split_fasta == "yes", int(args.id_maxlen),
        acc_to_orig_gb=acc_to_orig_gb, topology=args.topology)
    samples = linear_lineages_first(samples, args)
    logging.info(f"[ID_MAP] {id_map}  (total: {len(samples)})")

    results = run_annotation_stage(
        samples=samples, phase_dirname="01.anno_round1",
        out_root=out_root, args=args,
        fasta_selector=lambda si: si.internal_fasta)
    samples, results = apply_topology(samples, results, args, out_root, reports_dir(out_root))

    summary    = RunSummary()
    final_dir  = out_root / "04.final_gb"
    final_dir.mkdir(parents=True, exist_ok=True)
    summary_rows: List[str] = []

    for si in samples:
        r       = results.get(si.internal_id, {})
        best_gb = Path(r["best_gb"]) if r.get("best_gb") else None
        if not r.get("ok") or not best_gb or not best_gb.exists():
            summary.annotation_failed.append(si.sample_name)
            continue
        dest = final_dir / f"{si.sample_name}.gbf"
        try:
            pcg, trna, rrna, mode, orient_ok = finalize_genbank_output(
                sample_name=si.sample_name, body_src_gb=best_gb,
                metadata_src_gb=si.orig_gb_path, dest_gb=dest,
                which_gene_first=orientation_genes(args),
                force_forward=args.force_forward,
                preserve_metadata=(args.preserve_phase1_metadata == "yes")
                    if hasattr(args, "preserve_phase1_metadata") else False,
            )
            summary_rows.append(f"{si.sample_name}\t{pcg}\t{trna}\t{rrna}\t{orientation_text(orient_ok, si)}\t{mode}")
        except Exception as exc:
            record_event(si.sample_name, f"final GenBank not written: {exc}", logging.ERROR)
            summary.final_write_failed.append(si.sample_name)

    write_nonempty_list(reports_dir(out_root) / "annotation_failed.txt",
                        summary.annotation_failed, "# annotation core non-zero exit")
    write_nonempty_list(reports_dir(out_root) / "final_write_failed.txt",
                        summary.final_write_failed, "# Failed during final write")
    write_table(
        out_root / "annotation-summary.tsv",
        "sample_name\tPCG\ttRNA\trRNA\tcox1_orientation_ok\tmetadata_mode",
        summary_rows)
    run_rescue(args, final_dir, out_root)
    run_rrna_hmm(args, final_dir, out_root)
    run_intron_fix(args, final_dir, out_root)
    reorient_after_introns(args, final_dir)
    merge_dir = run_curate_and_submit(args, final_dir, out_root).gb_dir
    refresh_summary(out_root, merge_dir)
    all_final_gbs = sorted(                      # rglob: the no_<gene>_gene/ records are merged too
        p for p in merge_dir.rglob("*.gbf")
        if p.is_file() and p.stat().st_size > 0)
    if all_final_gbs:
        merge_gb_files(all_final_gbs, out_root / "all.final.gbf")
    logging.info(f"[PLOT] {run_plots(args, merge_dir, out_root).summary}")
    (reports_dir(out_root) / "run_params.json").write_text(
        json.dumps(vars(args), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    summary.log()
    logging.warning(KNOWN_LIMITATION)
    write_orientation_fallback(out_root, samples)
    if args.keep_intermediate == "no":
        cleanup_intermediates(out_root, samples, summary.all_problem_samples())
    if samples and len(summary.annotation_failed) == len(samples):
        logging.warning(f"[FAIL] every record failed in the annotation core ({len(samples)}): see reports/failed_logs/")
        return 1
    return 0


def workflow_reorient(args: argparse.Namespace) -> int:
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "pipeline.log", debug=getattr(args, "debug", False))
    log_environment(args)
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))

    in_path  = Path(args.input).resolve()
    gb_files = collect_genbank_files(in_path)
    if not gb_files:
        logging.error("[ERROR] No GB files found.")
        return 2

    out_dir     = out_root / "reoriented"
    no_gene_dir = out_root / f"no_{args.which_gene_first}_gene"
    out_dir.mkdir(parents=True, exist_ok=True)
    no_gene_dir.mkdir(parents=True, exist_ok=True)

    summary      = RunSummary()
    summary_rows: List[str] = []

    for gb_in in gb_files:
        sample     = safe_name(gb_in.stem)
        gb_out     = out_dir     / f"{sample}.gbf"
        gb_no_gene = no_gene_dir / f"{sample}.gbf"
        try:
            n_ok, n_no, no_gene_ids = reorient_one_gb(
                gb_in, gb_out, gb_no_gene, orientation_genes(args), args.force_forward)
            if n_no > 0:
                logging.warning(f"[NO_TARGET_GENE] '{sample}': "
                                 f"{n_no} record(s) missing '{args.which_gene_first}'")
                summary.no_target_gene.append(sample)
            if n_ok > 0 and gb_out.exists():
                pcg, trna, rrna = counts_from_genbank(gb_out)
                summary_rows.append(f"{sample}\t{pcg}\t{trna}\t{rrna}\tyes")
            elif n_ok == 0 and n_no > 0:
                safe_copy(gb_in, gb_out)
                pcg, trna, rrna = counts_from_genbank(gb_out)
                summary_rows.append(f"{sample}\t{pcg}\t{trna}\t{rrna}\tno")
                logging.warning(f"[FALLBACK_COPY] '{sample}': original copied unchanged")
        except Exception as exc:
            logging.error(f"[FAIL] {sample}: {exc}")
            summary.reorient_failed.append(sample)

    write_nonempty_list(out_root / "no_target_gene_samples.txt",
                        summary.no_target_gene, f"# samples without {args.which_gene_first}")
    write_nonempty_list(out_root / "failed_samples.txt",
                        summary.all_problem_samples(), "# all problematic samples")
    write_table(out_root / "annotation-summary.tsv", "sample_name\tPCG\ttRNA\trRNA\tcox1_orientation_ok",
                summary_rows)
    summary.log()
    return 0


def fallback_src_tag(mode: str, iid: str, reannotated_ids: set) -> str:
    """annotation-summary body_source of a sample whose final body is the shifted round-1 file: 'shifted_round1'
    when --if_reannotation_2 auto decided that round 2 is not needed, 'phase1b_reoriented' when round 2 was
    meant to run but failed (or is off)."""
    return "shifted_round1" if mode == "auto" and iid not in reannotated_ids else "phase1b_reoriented"


def workflow_run(args: argparse.Namespace) -> int:
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    new_reports(out_root)
    setup_root_logger(out_root / "run.log", console_level=logging.WARNING, mode="w",   # one run, one log
                      debug=getattr(args, "debug", False))
    log_environment(args)
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))
    preflight(args)
    steps = Steps(11)
    say("")
    say("MitoM - mitochondrial genome annotation")
    say(f"    input        : {Path(args.input).resolve()}")
    say(f"    output       : {out_root}")
    say(f"    clade        : {args.clade}    genetic code: {args.genetic_code}")
    say(f"    parallel     : {args.max_tasks} records at a time x {args.threads} threads")
    say(f"    annotation   : {ENGINE_LABEL[getattr(args, 'engine', 'native')]}")
    say(f"    round 2      : {args.if_reannotation_2}  (auto = only records whose annotation is broken by the cox1 shift)")
    if getattr(args, "engine", "native") == "native":
        say(f"    topology     : {args.topology}" + ("  (as labelled, checked at the sequence ends)"
                                                     if args.topology == "auto" else ""))
    say(f"    full log     : {out_root / 'run.log'}")

    in_path     = Path(args.input).resolve()
    split_dir   = out_root / "00.internal_fasta"
    gb_conv_dir = out_root / "00.gb_converted_fasta"

    steps.start("input", "Read input and remove duplicate records")
    in_files, input_type, acc_to_orig_gb = iter_inputs(in_path, gb_to_fasta_dir=gb_conv_dir,
                                                       orig_dir=out_root / "original_annotation",
                                                       dedup_report=reports_dir(out_root) / "dedup_report.tsv")
    logging.info(f"[INPUT_TYPE] {input_type}  ({len(in_files)} file(s))")
    samples, id_map = split_or_sanitize_fastas(
        in_files, split_dir, args.if_split_fasta == "yes", int(args.id_maxlen),
        acc_to_orig_gb=acc_to_orig_gb, topology=args.topology)
    samples = linear_lineages_first(samples, args)
    total = len(samples)
    logging.info(f"[ID_MAP] {id_map}  (total: {total})")
    if input_type == "converted_from_gb":
        st = LAST_DEDUP_STATS
        steps.info(f"GenBank records read: {st.get('total_in', '?')}")
        dropped = [(st.get(k, 0), why) for k, why in (
            ("dup_refseq", "GenBank record that a RefSeq NC_ record was made from (NC_ kept)"),
            ("dup_version", "older version of the same accession"), ("dup_exact", "same accession twice"),
            ("dup_seq", "identical sequence"), ("no_sequence", "no sequence"))]
        for n, why in dropped:
            if n:
                steps.info(f"removed {n}: {why}")
        steps.info("the original annotation of every kept record is saved for comparison (original_annotation/)")
        steps.done(f"{total} records to annotate; reasons per record: reports/dedup_report.tsv")
    else:
        steps.done(f"{total} FASTA sequences to annotate")

    summary = RunSummary()
    check_transl_tables(out_root / "original_annotation", int(args.genetic_code))

    steps.start("annotate", f"Round-1 annotation ({ENGINE_LABEL[getattr(args, 'engine', 'native')]}: protein genes, rRNAs, tRNAs)")
    results1 = run_annotation_stage(
        samples=samples, phase_dirname="01.anno_round1",
        out_root=out_root, args=args,
        fasta_selector=lambda si: si.internal_fasta)
    for si in samples:
        if not results1.get(si.internal_id, {}).get("ok"):
            summary.annotation_failed.append(si.sample_name)
    samples, results1 = apply_topology(samples, results1, args, out_root, reports_dir(out_root), steps.info)
    steps.done(f"{total - len(summary.annotation_failed)} annotated, {len(summary.annotation_failed)} failed"
               + (f" ({', '.join(summary.annotation_failed[:5])}{' ...' if len(summary.annotation_failed) > 5 else ''})"
                  if summary.annotation_failed else ""))

    steps.start("reorient", f"Re-orient: plus strand, sequence starts with {args.which_gene_first} "
                f"(if missing: {args.fallback_genes})")

    logging.info(
        f"[PHASE1b] reorient gene='{args.which_gene_first}' force_forward={args.force_forward}")
    reor_fasta_dir  = out_root / "00.reoriented_fasta_internal"
    reor_fasta_dir.mkdir(parents=True, exist_ok=True)
    reor_ok_ids:    List[str]        = []
    no_gene_iids:   set              = set()
    partial_no_gene_gb: Dict[str, Path] = {}     # contigs without the gene: merged back into the final record

    for si in samples:
        r1      = results1.get(si.internal_id, {})
        best_gb = Path(r1["best_gb"]) if r1.get("best_gb") else None
        if not best_gb or not best_gb.exists():
            continue
        reor_dir   = out_root / "01.anno_round1" / si.internal_id / "reoriented"
        reor_dir.mkdir(parents=True, exist_ok=True)
        gb_reor    = reor_dir / f"{si.internal_id}.reoriented.gbf"
        gb_no_gene = reor_dir / f"{si.internal_id}.no_{args.which_gene_first}.gbf"
        try:
            n_ok, n_no, no_gene_ids = reorient_one_gb(
                best_gb, gb_reor, gb_no_gene, orientation_genes(args), args.force_forward)
            if n_ok == 0 and n_no > 0:
                record_event(si.sample_name, f"no orientation gene ({args.which_gene_first}, {args.fallback_genes}) "
                             f"in its {n_no} sequence(s): kept as annotated, no round 2 ({si.internal_id})")
                summary.no_target_gene.append(si.sample_name)
                no_gene_iids.add(si.internal_id)
            elif n_no > 0:
                logging.info(
                    f"[PARTIAL_NO_GENE] '{si.sample_name}' ({si.internal_id}): "
                    f"{n_no} minor contig(s) lack '{args.which_gene_first}'; "
                    f"primary reoriented OK — proceeding to Phase-2.")
                if gb_no_gene.exists():
                    partial_no_gene_gb[si.internal_id] = gb_no_gene
            if n_ok > 0 and gb_reor.exists():
                reor_ok_ids.append(si.internal_id)
                gb_to_internal_fasta(gb_reor, reor_fasta_dir / f"{si.internal_id}.fasta",
                                     si.internal_id)
        except Exception as exc:
            record_event(si.sample_name, f"re-orientation failed ({si.internal_id}): {exc}", logging.ERROR)
            summary.reorient_failed.append(si.sample_name)      # an error, not a record without the gene
            no_gene_iids.add(si.internal_id)

    fallback_used = sum(1 for si in samples if si.internal_id in reor_ok_ids
                        and any(k.startswith(si.internal_id) for k in ORIENTATION_FALLBACK))
    steps.done(f"{len(reor_ok_ids)} re-oriented"
               + (f" ({fallback_used} by a fallback gene)" if fallback_used else "")
               + (f", {len(no_gene_iids)} without any orientation gene (kept as annotated)" if no_gene_iids else ""))

    steps.start("round2", "Check the annotation after re-orientation; annotate again where it is offset (round 2)")
    do_phase2 = args.if_reannotation_2 in ("yes", "auto")
    results2:  Dict[str, Dict] = {}
    reannotated_ids: set = set()
    if do_phase2:
        phase2_samples = [si for si in samples if si.internal_id in reor_ok_ids]
        n_checked = len(phase2_samples)
        if args.if_reannotation_2 == "auto":
            phase2_samples = select_reannotation(phase2_samples, results1, out_root,
                                                 int(getattr(args, "genetic_code", None) or 5))
            steps.info(f"checked {n_checked}: {len(phase2_samples)} need round 2 (gene cut at the old ends, gene "
                       f"changed by the shift, gene across the new start, standard gene missing); "
                       f"{n_checked - len(phase2_samples)} keep the shifted round-1 annotation")
            steps.info("decision and reason per record: reports/reannotation_decision.tsv")
        else:
            steps.info(f"--if_reannotation_2 yes: all {n_checked} re-oriented records are annotated again")
        results2 = run_phase2(phase2_samples, args, out_root, reor_fasta_dir)
        reannotated_ids = {si.internal_id for si in phase2_samples}
        n_fail2 = 0
        for si in phase2_samples:
            if not results2.get(si.internal_id, {}).get("ok"):
                n_fail2 += 1
                record_event(si.sample_name, "round-2 annotation failed: the shifted round-1 annotation is kept")
        steps.done(f"round 2: {len(phase2_samples) - n_fail2} annotated"
                   + (f", {n_fail2} failed (they keep the shifted round-1 annotation)" if n_fail2 else ""))
    else:
        steps.done("--if_reannotation_2 no: skipped")

    steps.start("write", "Write the annotation; copy the original GenBank header (organism, references ...) back")

    final_dir         = out_root / "04.final_gb"
    no_gene_final_dir = final_dir / f"no_{args.which_gene_first}_gene"
    final_dir.mkdir(parents=True, exist_ok=True)
    no_gene_final_dir.mkdir(parents=True, exist_ok=True)

    summary_rows:    List[str] = []
    meta_log_rows:   List[str] = []

    for si in samples:
        sample      = si.sample_name
        iid         = si.internal_id
        is_no_gene  = iid in no_gene_iids
        r1          = results1.get(iid, {})
        metadata_src = si.orig_gb_path

        chosen_gb: Optional[Path] = None
        src_tag   = "phase1_raw"

        if do_phase2 and not is_no_gene:
            r2 = results2.get(iid, {})
            if r2.get("best_gb") and Path(r2["best_gb"]).exists():
                chosen_gb = Path(r2["best_gb"])
                src_tag   = "phase2"

        if chosen_gb is None:
            gb_reor = (out_root / "01.anno_round1" / iid / "reoriented"
                       / f"{iid}.reoriented.gbf")
            if gb_reor.exists():
                chosen_gb = gb_reor
                src_tag   = fallback_src_tag(args.if_reannotation_2 if do_phase2 else "no", iid, reannotated_ids)

        if chosen_gb is None and r1.get("best_gb") and Path(r1["best_gb"]).exists():
            chosen_gb = Path(r1["best_gb"])

        if chosen_gb is None or not chosen_gb.exists():
            if sample not in summary.annotation_failed:
                summary.final_write_failed.append(sample)
            continue

        body_for_final = chosen_gb
        extra_no_gene_gb = partial_no_gene_gb.get(iid)
        if extra_no_gene_gb and extra_no_gene_gb.exists():
            chosen_is_phase1_raw = False
            if r1.get("best_gb"):
                try:
                    chosen_is_phase1_raw = chosen_gb.resolve() == Path(r1["best_gb"]).resolve()
                except Exception:
                    chosen_is_phase1_raw = False
            if not chosen_is_phase1_raw:
                merged_body_dir = out_root / "04.final_body_merged"
                merged_body_dir.mkdir(parents=True, exist_ok=True)
                merged_body = merged_body_dir / f"{iid}.with_no_gene.gbf"
                merge_gb_files([chosen_gb, extra_no_gene_gb], merged_body)
                body_for_final = merged_body
                src_tag += "+partial_no_gene_merged"
                logging.info(
                    f"[FINAL_PARTIAL_NO_GENE_MERGE] '{sample}': merged "
                    f"{extra_no_gene_gb.name} back into final body.")

        dest = (no_gene_final_dir if is_no_gene else final_dir) / f"{sample}.gbf"

        try:
            pcg, trna, rrna, meta_mode, orient_ok = finalize_genbank_output(
                sample_name=sample,
                body_src_gb=body_for_final,                          # with the contigs without the gene
                metadata_src_gb=metadata_src,
                dest_gb=dest,
                which_gene_first=orientation_genes(args),
                force_forward=args.force_forward,
                preserve_metadata=(args.preserve_phase1_metadata == "yes"),
                smart_metadata=(args.smart_metadata == "yes"),
                auto_accept=(args.auto_accept_metadata == "yes"),
            )
            summary_rows.append(
                f"{sample}\t{pcg}\t{trna}\t{rrna}\t{orientation_text(orient_ok, si)}\t{src_tag}\t{meta_mode}")
            meta_log_rows.append(
                f"{sample}\t{iid}\t{metadata_src or 'NA'}\t{body_for_final}\t{dest}\t{meta_mode}")

            if "transferred" in meta_mode:
                summary.metadata_transferred.append(sample)
                if "compat" in meta_mode:
                    summary.metadata_compat_match.append(sample)
                    summary.metadata_key_mismatch.append(sample)
            elif meta_mode in ("metadata_no_match", "metadata_src_not_usable"):
                summary.metadata_skipped.append(sample)
            elif meta_mode == "metadata_failed":
                summary.metadata_transfer_failed.append(sample)

        except Exception as exc:
            record_event(sample, f"final GenBank not written: {exc}", logging.ERROR)
            summary.final_write_failed.append(sample)

    write_table(
        out_root / "annotation-summary.tsv",
        "sample_name\tPCG\ttRNA\trRNA\tcox1_orientation_ok\tbody_source\tmetadata_mode",
        summary_rows)
    steps.done(f"{len(summary_rows)} records written"
               + (f", original header restored in {len(summary.metadata_transferred)}" if summary.metadata_transferred else "")
               + (f", {len(summary.final_write_failed)} failed" if summary.final_write_failed else ""))
    write_table(
        reports_dir(out_root) / "metadata_transfer_log.tsv",
        "sample_name\tinternal_id\tmetadata_src\tbody_src\tdest\tmode",
        meta_log_rows)
    write_nonempty_list(reports_dir(out_root) / "annotation_failed.txt",
                        summary.annotation_failed, "# annotation core returned non-zero exit")
    write_nonempty_list(reports_dir(out_root) / f"no_{args.which_gene_first}_samples.txt",
                        summary.no_target_gene,
                        f"# samples where '{args.which_gene_first}' was not found")
    write_nonempty_list(reports_dir(out_root) / "final_write_failed.txt",
                        summary.final_write_failed, "# failed during final GB write")
    write_nonempty_list(reports_dir(out_root) / "reorient_failed.txt",
                        summary.reorient_failed, "# failed during reorientation")
    write_nonempty_list(reports_dir(out_root) / "metadata_transfer_failed.txt",
                        summary.metadata_transfer_failed, "# metadata transfer failed")
    write_nonempty_list(reports_dir(out_root) / "metadata_key_mismatch.txt",
                        summary.metadata_key_mismatch,
                        "# accession key mismatch during metadata transfer")
    write_nonempty_list(reports_dir(out_root) / "all_problem_samples.txt",
                        summary.all_problem_samples(), "# all problematic samples")
    (reports_dir(out_root) / "run_params.json").write_text(
        json.dumps(vars(args), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    steps.start("add_missing", "Add what the core missed: MITOS2 rescue (optional), 12S/16S by rRNA HMM, introns")
    for name, step in (("rescue", run_rescue), ("rrna_hmm", run_rrna_hmm), ("introns", run_intron_fix)):
        steps.info(steps.part(name, step, args, final_dir, out_root).summary)
    n_rot = reorient_after_introns(args, final_dir)
    if n_rot:
        steps.info(f"re-oriented to the first exon : {n_rot} record(s) whose {args.which_gene_first} gained an exon upstream")
    steps.done()

    steps.start("curate", "Boundary rules for protein genes and rRNAs, NCBI gene names")
    curated = run_curate(args, final_dir, out_root)
    merge_dir = curated.gb_dir
    steps.done(curated.summary, curated.ok)

    steps.start("ncbi", "NCBI check (table2asn) and submission files")
    submitted = run_submission(args, merge_dir, out_root) if merge_dir != final_dir else curated
    steps.done(submitted.summary, submitted.ok)
    refresh_summary(out_root, merge_dir)
    all_final_gbs = sorted(                      # rglob: the no_<gene>_gene/ records are merged too
        p for p in merge_dir.rglob("*.gbf")
        if p.is_file() and p.stat().st_size > 0
    )
    if all_final_gbs:
        n_merged = merge_gb_files(all_final_gbs, out_root / "all.final.gbf")
        logging.info(f"[MERGED_ALL] {out_root / 'all.final.gbf'}  ({n_merged} records)")
    else:
        logging.warning("[MERGED_ALL] No final GB files to merge.")

    n_ok = len(summary_rows)
    logging.info(f"[FINAL_DIR] {final_dir}")
    logging.info(f"[OUTPUT_COUNT] {n_ok}/{total}")
    logging.info(f"[METADATA_MODE] preserve={args.preserve_phase1_metadata}")
    summary.log()
    logging.info(KNOWN_LIMITATION)
    write_orientation_fallback(out_root, samples)
    if args.keep_intermediate == "no":
        cleanup_intermediates(out_root, samples, summary.all_problem_samples())

    steps.start("compare", "Compare with the original annotation")
    compared = run_compare(args, merge_dir, out_root, original=out_root / "original_annotation"
                           if input_type == "converted_from_gb" else None)
    steps.done(compared.summary, compared.ok)

    steps.start("maps", "Genome maps (circular plots)")
    maps = run_plots(args, merge_dir, out_root)
    steps.done(maps.summary, maps.ok)

    steps.start("output", "Organise the output")
    organize_outputs(out_root, merge_dir, keep=args.keep_intermediate == "yes",
                     reports_named=(f"no_{args.which_gene_first}_samples.txt",))
    for name, what in OUTPUT_LAYOUT:
        if (out_root / name).exists():
            steps.info(f"{name:28s} {what}")
    steps.done(f"explained in {out_root / 'README.txt'}")
    say("")
    say(f"MitoM finished in {duration(time.time() - steps.run_t0)}: {out_root}")
    if total and len(summary.annotation_failed) == total:
        logging.warning(f"[FAIL] every record failed in the annotation core ({total}): see reports/failed_logs/")
        return 1
    return 0


def workflow_transfer_metadata(args: argparse.Namespace) -> int:
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "transfer_metadata.log", debug=getattr(args, "debug", False))
    log_environment(args)
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))

    metadata_files = collect_genbank_files(Path(args.metadata).resolve())
    body_files     = collect_genbank_files(Path(args.new_body).resolve())
    if not metadata_files:
        logging.error(f"[ERROR] No metadata GB files: {args.metadata}")
        return 2
    if not body_files:
        logging.error(f"[ERROR] No body GB files: {args.new_body}")
        return 2

    out_dir = out_root / "transferred_gb"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_rows: List[str] = []

    def _process(body_file: Path, meta_file: Optional[Path], out_file: Path) -> str:
        if args.if_replace_metadata == "no":
            if args.overwrite_new_body == "no":
                safe_copy(body_file, out_file)
            log_rows.append(f"{body_file}\t{meta_file or 'NA'}\t{out_file}\tskipped_by_user\t0\t0")
            return "skip"
        if meta_file is None:
            if args.overwrite_new_body == "no":
                safe_copy(body_file, out_file)
            log_rows.append(f"{body_file}\tNA\t{out_file}\tno_match\t0\t0")
            return "no_match"
        ok, n_rep, n_un, unmatched, reason, compat = transfer_gb_metadata_file(
            meta_file, body_file, out_file,
            smart=(args.smart_metadata == "yes"), auto_accept=(args.auto_accept_metadata == "yes"))
        if not ok and args.overwrite_new_body == "no":
            safe_copy(body_file, out_file)
        log_rows.append(
            f"{body_file}\t{meta_file}\t{out_file}\t"
            f"{'ok' if ok else 'no_usable_metadata'}\t{n_rep}\t{n_un}")
        return "ok" if ok else "skip"

    is_single = (len(metadata_files) == 1 and len(body_files) == 1
                 and Path(args.metadata).is_file() and Path(args.new_body).is_file())
    if is_single:
        out_name = getattr(args, "output_name", None) or body_files[0].name
        out_file = body_files[0] if args.overwrite_new_body == "yes" else (out_dir / out_name)
        _process(body_files[0], metadata_files[0], out_file)
        write_table(out_root / "transfer_metadata.tsv",
                             "body_file\tmetadata_file\toutput_file\tstatus\treplaced\tunmatched",
                             log_rows)
        return 0

    meta_by_stem: Dict[str, Path] = {}           # the same stem in two folders is an error, never a silent choice
    dup_stems: Dict[str, List[Path]] = {}
    for mf in metadata_files:
        if mf.stem in meta_by_stem:
            dup_stems.setdefault(mf.stem, [meta_by_stem[mf.stem]]).append(mf)
        else:
            meta_by_stem[mf.stem] = mf
    if dup_stems:
        for stem, paths in sorted(dup_stems.items()):
            logging.error(
                f"[DUP_METADATA_STEM] '{stem}' appears in multiple metadata files: "
                + ", ".join(str(p) for p in paths))
        return 2
    n_ok = n_skip = n_fail = 0
    for bf in body_files:
        matched  = meta_by_stem.get(bf.stem)
        out_file = bf if args.overwrite_new_body == "yes" else (out_dir / bf.name)
        try:
            status = _process(bf, matched, out_file)
            if status == "ok":
                n_ok += 1
            else:
                n_skip += 1
        except Exception:
            if args.overwrite_new_body == "no":
                safe_copy(bf, out_file)
            n_fail += 1
            log_rows.append(f"{bf}\t{matched or 'NA'}\t{out_file}\tfail\t0\t0")
    write_table(out_root / "transfer_metadata.tsv",
                         "body_file\tmetadata_file\toutput_file\tstatus\treplaced\tunmatched",
                         log_rows)
    logging.info(f"[SUMMARY] transferred={n_ok}  skipped={n_skip}  failed={n_fail}")
    return 0


def workflow_replace_organism(args: argparse.Namespace) -> int:
    """Replace /organism= qualifiers in FEATURES source blocks (`run` does this in step 3.5 of
    finalize_genbank_output)."""
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "replace_organism.log", debug=getattr(args, "debug", False))
    log_environment(args)
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))

    meta_files = collect_genbank_files(Path(args.metadata).resolve())
    body_files = collect_genbank_files(Path(args.new_body).resolve())
    if not meta_files:
        logging.error(f"[ERROR] No metadata GB files: {args.metadata}")
        return 2
    if not body_files:
        logging.error(f"[ERROR] No body GB files: {args.new_body}")
        return 2

    logging.info(f"[ORG_INDEX] Building organism+source index from {len(meta_files)} file(s)…")
    org_index = build_organism_index(meta_files)
    logging.info(f"[ORG_INDEX] {len(org_index)} accession variant(s) indexed.")
    if not org_index:
        logging.error("[ERROR] Organism index is empty — check metadata GB files.")
        return 2

    out_dir = out_root / "replaced_organism"
    out_dir.mkdir(parents=True, exist_ok=True)

    log_rows: List[str] = []
    n_ok = n_miss = n_fail = 0

    for bf in body_files:
        out_file = bf if args.overwrite_new_body == "yes" else (out_dir / bf.name)
        try:
            n_rep, n_nf, n_tot = replace_organism_in_gb_file(bf, out_file, org_index)
            log_rows.append(f"{bf.name}\t{out_file}\t{n_rep}\t{n_nf}\t{n_tot}\t")         # note column empty
            if n_rep > 0:
                n_ok += 1
                logging.info(
                    f"[DONE] {bf.name}: {n_rep}/{n_tot} record(s) updated, "
                    f"{n_nf} unchanged.")
            else:
                n_miss += 1
                logging.warning(
                    f"[NO_MATCH] {bf.name}: 0 records matched — organism/source not replaced.")
        except Exception as exc:
            n_fail += 1
            logging.error(f"[FAIL] {bf.name}: {exc}")
            log_rows.append(f"{bf.name}\t{out_file}\t0\t0\t0\tFAIL: {exc}")
            safe_copy(bf, out_file)

    write_table(
        out_root / "replace_organism.tsv",
        "source_file\toutput_file\treplaced\tnot_found\ttotal\tnote",
        log_rows)
    logging.info(
        f"[SUMMARY] files_updated={n_ok}  no_match={n_miss}  failed={n_fail}")
    return 0


def workflow_dedup(args: argparse.Namespace) -> int:
    """Standalone deduplication of GenBank and / or FASTA input: dedup.gb (the kept GenBank records), dedup.fasta
    (every kept record; GenBank ones as >ACCESSION.V definition topology=..., FASTA ones unchanged),
    dedup_report.tsv."""
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "dedup.log", debug=getattr(args, "debug", False))
    log_environment(args)
    files = collect_sequence_files(Path(args.input).resolve())
    if not files:
        print(f"mitom dedup: no GenBank or FASTA file in {args.input}", file=sys.stderr)
        return 2
    kept, _, stats = dedup_gb_records(files, dedup_report_path=out_root / "dedup_report.tsv")
    gb = [m.record for m in kept if record_format(m.source_file) == "genbank"]
    if gb:
        SeqIO.write(gb, str(out_root / "dedup.gb"), "genbank")
    SeqIO.write([gb_fasta_record(m.record) if record_format(m.source_file) == "genbank" else m.record for m in kept],
                str(out_root / "dedup.fasta"), "fasta")
    logging.info(f"[DEDUP] in={stats['total_in']} kept={stats['kept']} no_seq={stats['no_sequence']} "
                 f"dup_exact={stats['dup_exact']} dup_version={stats['dup_version']} "
                 f"dup_refseq={stats['dup_refseq']} dup_seq={stats['dup_seq']} -> {out_root}")
    return 0


def workflow_gb2fasta(args: argparse.Namespace) -> int:
    """GenBank (any GenBank suffix: NCBI .gb, MitoZ .gbf ...) -> FASTA: one <name>.fasta per input file, or one
    <accession>.fasta per record with --split; headers >ACCESSION.V definition topology=circular|linear, ready for
    `mitom run` (the topology tag keeps genomes circular in a multi-record FASTA)."""
    files = iter_gb_files(Path(args.input))
    if not files:
        print(f"mitom gb2fasta: no GenBank file in {args.input} (suffixes {', '.join(GB_SUFFIXES)})", file=sys.stderr)
        return 2
    out = Path(args.out_root)
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    for f in files:
        recs = [gb_fasta_record(r) for r in SeqIO.parse(str(f), "genbank")]
        n += len(recs)
        if args.split:
            for r in recs:
                SeqIO.write([r], str(out / f"{safe_name(r.id)}.fasta"), "fasta")
        elif recs:
            SeqIO.write(recs, str(out / f"{f.stem}.fasta"), "fasta")
    print(f"[gb2fasta] {n} record(s) from {len(files)} GenBank file(s) -> {out}", file=sys.stderr)
    return 0


def workflow_intron(args: argparse.Namespace) -> int:
    """`mitom intron`: introns as `run` adds them - mito_intron.py (de novo, MFannot method) -> OUT, then
    mito_intron_ref.py (intron sites of annotated relatives) on OUT in place; --method denovo / ref runs one of them.
    Options this command does not know (args.extra) go to mito_intron.py, or to mito_intron_ref.py with --method ref.
    Returns the exit code of the first script that failed (the site transfer does not run after a failed de novo run)."""
    out = Path(args.out_root)
    if args.method != "ref":
        rc = run_tool(intron_denovo_cmd(args.input, out, args.genetic_code, args.threads, args.jobs, bin_dir=args.bin)
                      + args.extra, capture=False).returncode
        if rc or args.method == "denovo":
            return rc
    src, extra = (args.input, args.extra) if args.method == "ref" else (out, [])
    return run_tool(intron_ref_cmd(src, out, args.genetic_code, args.jobs, args.bin, args.lineage) + extra,
                    capture=False).returncode


# ── command line ─────────────────────────────────────────────────────────────────────────────────────────────────────
def positive_int_arg(x: str) -> int:
    try:
        v = int(x)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(f"expected positive integer, got {x!r}") from exc
    if v < 1:
        raise argparse.ArgumentTypeError(f"expected positive integer, got {x!r}")
    return v


def id_maxlen_arg(x: str) -> int:
    v = positive_int_arg(x)
    if v < len("S00001"):
        raise argparse.ArgumentTypeError(
            f"--id_maxlen must be >= {len('S00001')} to keep internal IDs unique")
    return v


_CLADE_DEFAULT = "auto"


def default_threads() -> int:
    """Threads per annotation job when --threads is not given: min(8, CPUs of this machine)."""
    return min(8, os.cpu_count() or 8)


def default_max_tasks() -> int:
    """Concurrent jobs when --jobs is not given: CPUs // default threads (at least 1)."""
    return max(1, (os.cpu_count() or 8) // default_threads())


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="mitom" if os.environ.get("MITOM_LAUNCHER") else "mitom_pipeline.py",   # sub-commands show `mitom run`
        description=(
            "MitoM - annotation workflow for animal mitochondrial genomes (MitoM's own annotation core, mito_core.py)\n\n"
            "  run               - full workflow: annotate -> re-orient -> round 2 where needed -> final GenBank ->\n"
            "                      rRNA HMM fallback, introns, curation, NCBI submission files, comparison, genome maps\n"
            "  import            - the same steps on existing MitoZ / MITOS2 results or any annotated GenBank\n"
            "  annotate          - one round of annotation and the steps after it (work layout: 04.final_gb/ ...)\n"
            "  reorient          - re-orient GenBank files to start at a target gene\n"
            "  transfer_metadata - copy GenBank headers from the original records\n"
            "  replace_organism  - restore /organism= of the FEATURES source from the original records\n"
            "  dedup             - deduplicate GenBank records (RefSeq NC_ kept)\n"
            "  intron            - introns as `run` adds them: de novo (MFannot method), then the intron sites of\n"
            "                      annotated relatives\n\n"
            "The annotation runs in threads or with Snakemake (--scheduler_backend; cluster: --smk_extra '--profile\n"
            "slurm'); everything after it runs in this process. `mitom <command> -h` lists the options.\n\n"
            "Citation:\n"
            "  Snakemake:  Mölder et al. (2021) F1000Res doi:10.12688/f1000research.29032.2\n"
            "  Biopython:  Cock et al. (2009) Bioinformatics doi:10.1093/bioinformatics/btp163\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("-i", "--input",  required=True,
                        help="Input FASTA file/dir or GenBank file/dir.")
        sp.add_argument("-o", "--out", dest="out_root", metavar="OUT", default="mitom_out",
                        help="Output folder (default: mitom_out).")
        old_spelling(sp, "--out_root", "out_root")
        sp.add_argument("--if_split_fasta", choices=["yes", "no"], default="yes")
        sp.add_argument("--id_maxlen", type=id_maxlen_arg, default=15,
                        help="Max length of internal sample IDs (default: 15).")
        sp.add_argument("--which_gene_first", default="cox1",
                        help="Gene used as sequence start (default: cox1).")
        sp.add_argument("--fallback_genes", default="cox2,cox3,cob,nad1",
                        help="Genes tried in order when --which_gene_first is absent "
                             "(default: cox2,cox3,cob,nad1; 'none' = leave such records unrotated). "
                             "Records oriented this way are listed in oriented_by_fallback.txt.")
        sp.add_argument("--force_forward", action="store_true", default=True,
                        help="RC the sequence when target gene is on minus strand (default: on).")
        sp.add_argument("--no_force_forward", action="store_false", dest="force_forward")
        sp.add_argument("--debug", action="store_true",
                        help="Debug mode: DEBUG lines, program versions and the environment in the log; the output of "
                             "every program in the logs; intermediate files kept (--keep_intermediate yes). "
                             "`mitom debug OUT` reads the log.")

    def add_annotation(sp: argparse.ArgumentParser) -> None:
        g = sp.add_argument_group("Annotation")
        g.add_argument("--core", choices=["native"], default="native", help=argparse.SUPPRESS)   # old scripts
        g.add_argument("--mitos_bin", default=None,
                       help="MITOS2 environment bin (runmitos). Default: found automatically (mito_env.py).")
        g.add_argument("--refdir", dest="mitos_refdir", metavar="DIR", default=None,
                       help="Folder that holds the MITOS2 reference data refseq89m. Default: found automatically.")
        old_spelling(g, "--mitos_refdir", "mitos_refdir")
        g.add_argument("--plot", choices=["yes", "no"], default="yes",
                       help="Circular genome map per record (mito_plot.py, matplotlib) -> plots/<record>.png and .pdf: "
                            "genes by strand and type, GC content, GC skew; with --infq also the read depth (minimap2 when installed, else bwa mem). "
                            "Default yes (about 2 s per record, --jobs in parallel).")
        g.add_argument("--plot_style", choices=["mitom", "circos"], default="mitom",
                       help="Map style: mitom (default) = genes on both sides of one circle, GC content and GC skew rings; "
                            "circos (after mitoPlotR) = a framed lane per strand, names along the circle, a position axis, "
                            "GC content, read depth.")
        g.add_argument("--curate", choices=["yes", "no"], default="yes",
                       help="yes (default): mito_curate.py boundary rules and NCBI gene names -> 05.curated/, and "
                            "NCBI submission files (.fsa/.tbl/.sqn/.val) -> 06.ncbi_submission/")
        g.add_argument("--sbt", default=None,
                       help="NCBI submitter template (template.sbt from submit.ncbi.nlm.nih.gov/genbank/template) "
                            "for the .sqn files")
        g.add_argument("--rescue", choices=["auto", "yes", "no"], default="auto",
                       help="Add protein genes / rRNAs the core missed from a MITOS2 annotation (mito_rescue.py; "
                            "nothing annotated is changed; needs the optional MITOS2 installation). auto (default) / "
                            "no = not run (the core retries other loci itself); yes = run for records missing a gene.")
        g.add_argument("--rrna_hmm", choices=["auto", "yes", "no"], default="auto",
                       help="Last-resort 12S/16S search for rRNAs still missing after the cores and the rescue, "
                            "with per-phylum profile HMMs built from NCBI rRNAs (mito_rrna_hmm.py; only missing "
                            "rRNAs are searched). auto = run when the library and nhmmer are found (default).")
        g.add_argument("--rrna_hmm_lib", default=None,
                       help="Per-phylum rRNA HMM library (mito_rrna_hmm.py build). Default: models/rrna_hmm.")
        g.add_argument("--rrna_method", choices=["hmm", "cm"], default="hmm",
                       help="12S/16S search of the core: hmm (default, faster) = per-phylum profile HMMs, "
                            "MITOS2's covariance models where no model or hit; cm = MITOS2's covariance models only. "
                            "Same accuracy (2026-10-28, 450 leave-out records: found 98.5 / 98.5 %%, both ends exact "
                            "47.4 / 47.6 %%); hmm takes seconds and ~0.1 GB, cm ~10 s and ~0.5 GB per gene.")
        g.add_argument("--phylum", default=None,
                       help="Phylum of all input records for --rrna_hmm (default: the phylum found by --clade auto, "
                            "or the phylum of an explicit --clade; with --clade Animal and FASTA input give it here, "
                            "otherwise the record's own taxonomy is used when it has one).")
        g.add_argument("--intron_fix", choices=["auto", "yes", "no"], default="auto",
                       help="Add introns with mito_intron.py (MFannot method: Exonerate "
                            "protein2genome + MFannot boundary rules + Rfam intron type), then with "
                            "mito_intron_ref.py (intron sites of annotated relatives, models/intron_sites). "
                            "auto = run when exonerate is found (default).")
        g.add_argument("--intron_bin", default=None,
                       help="Directory containing exonerate, cmsearch, blastp, tblastn and blastn (e.g. a conda env "
                            "bin).")
        g.add_argument("--intron_lib", action="append", default=[],
                       help="Protein library for intron models in MFannot header format "
                            "('>Species mt cox1; ...'), repeatable. Default: MFannot's curated "
                            "collections bundled in models/mfannot_proteins.")
        g.add_argument("--intron_ref_gb", action="append", default=[],
                       help="GenBank file/dir of close relatives used as extra homologs for intron models "
                            "(repeatable).")
        g.add_argument("--keep_intermediate", choices=["yes", "no"], default="no",
                       help="Keep the intermediate folders (00.*, 01.*, 03.*, 04.*, 05.*, _snakemake_*; run: moved to "
                            "work/). Default no: deleted after the final output is written; core logs of problem "
                            "samples are always kept (reports/failed_logs/).")
        g.add_argument("--threads",      type=positive_int_arg, default=default_threads(),
                       help="Threads per annotation job (default: min(8, number of CPUs) = %(default)s).")
        g.add_argument("--jobs", dest="max_tasks", metavar="N", type=positive_int_arg, default=default_max_tasks(),
                       help="Records annotated at the same time (default: CPUs // threads = %(default)s).")
        old_spelling(g, "--max_tasks", "max_tasks", type=positive_int_arg)
        g.add_argument("--clade",        default=_CLADE_DEFAULT,
                       help="Taxonomic group (sets the genetic code and group rules). auto (default): from the taxonomy of the GenBank "
                            "input (one phylum per run; class names such as Vertebrata, Gastropoda, Insecta count "
                            "for their phylum; a group below or Animal). A taxonomy that names no "
                            "known phylum, several phyla, or CDS transl_tables that differ between records are "
                            "errors; input with NO taxonomy at all (FASTA, GenBank without a lineage) falls back to "
                            "Annelida-segmented-worms with a WARNING - give the clade explicitly for any other "
                            f"taxon: {', '.join(ALL_CLADES)}.")
        g.add_argument("--genetic_code", type=int, default=None,
                       help="Genetic code (default with --clade auto: the CDS transl_table of the input, else the "
                            "phylum's usual code; otherwise the code of the clade - 2 Chordata, 9 Echinodermata, "
                            "4 Porifera, else 5; required for Cnidaria, Platyhelminthes, Chaetognatha, Animal).")
        g.add_argument("--species_name", default=None, help=argparse.SUPPRESS)   # hidden; passed to MitoZ (--engine mitoz)

        fq = sp.add_argument_group("FASTQ (optional - read depth ring of the genome maps, --plot)")
        fq.add_argument("--infq", default=None,
                        help="Directory with the paired FASTQ files of the samples (matched by sample name).")
        fq.add_argument("--fq_position", choices=["flat", "sample_dir"], default="flat")
        fq.add_argument("--suffix_fq", default="_1.clean.fq.gz,_2.clean.fq.gz")

        sch = sp.add_argument_group("Scheduler")
        sch.add_argument("--scheduler_backend",
                         choices=["auto", "snakemake", "threads"], default="auto",
                         help=(
                             "Dispatch backend (default: auto).\n"
                             "  auto      — use Snakemake if available, else threads\n"
                             "  snakemake — always use Snakemake\n"
                             "  threads   — always use ThreadPoolExecutor\n"
                             "Annotation is parallelised; finalization runs\n"
                             "in the main process regardless of backend."
                         ))
        sch.add_argument("--snakemake_bin",    default=None)
        sch.add_argument("--smk_cores",        type=positive_int_arg, default=64)
        sch.add_argument("--smk_mem_mb",       type=positive_int_arg, default=256000)
        sch.add_argument("--smk_latency_wait", type=positive_int_arg, default=60)
        sch.add_argument("--smk_extra",        action="append", default=[],
                         metavar="FLAG",
                         help="Extra Snakemake flags appended verbatim.\n"
                              "Example: --smk_extra '--profile slurm'")

    def add_topology(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--topology", choices=["auto", "label", "circular", "linear"], default="auto",
                        help="Circular or linear. auto (default): as labelled - GenBank LOCUS line; FASTA header "
                             "topology=circular|linear; else one sequence of >= 10 kb circular, several in a file or a "
                             "shorter one linear - and the protein genes at the two ends checked against their "
                             "reference proteins: a gene running from the end across to the start makes the record "
                             "circular, a gene cut at an end of a circular-labelled record makes it linear, and so do "
                             "jellyfish / hydroids / calcareous sponges (linear mtDNA) without a gene across the ends; "
                             "such records are annotated again (reports/topology_decision.tsv). label: the label only. "
                             "circular / linear: every record.")

    def add_engine(sp: argparse.ArgumentParser) -> None:
        g = sp.add_argument_group("Annotation engine")
        g.add_argument("--engine", choices=["native", "mitoz", "mitos2"], default="native",
                       help="Program that finds the genes: native (default) = MitoM's own core (mito_core.py); mitoz = "
                            "MitoZ 3.6, mitos2 = MITOS2 2.1.10 (not installed with MitoM - install them yourself, see "
                            "mito_engines.py; found on PATH or in a conda environment). Every later MitoM step runs on "
                            "their annotation: reorientation, round 2, curation, NCBI files, genome maps.")
        g.add_argument("--mitoz_bin", default=None,
                       help="Folder with the mitoz program (--engine mitoz). Default: found automatically.")

    def add_metadata(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--preserve_phase1_metadata", choices=["yes", "no"], default="yes",
                        help="Restore original GenBank header metadata into final files "
                             "(default: yes; requires GenBank input).")
        sp.add_argument("--smart_metadata",       choices=["yes", "no"], default="no",
                        help="Enable fuzzy/interactive accession matching (default: no).")
        sp.add_argument("--auto_accept_metadata", choices=["yes", "no"], default="no",
                        help="Auto-accept best fuzzy match (requires --smart_metadata yes).")

    # ── annotate ──────────────────────────────────────────────────────────
    a = sub.add_parser("annotate", help="Batch annotation only (one round of the annotation engine).")
    add_common(a)
    add_annotation(a)
    add_engine(a)
    add_topology(a)
    a.add_argument("--preserve_phase1_metadata", choices=["yes", "no"], default="no")

    # ── reorient ──────────────────────────────────────────────────────────
    r = sub.add_parser("reorient", help="Reorient GenBank files to start at a target gene.",
                       description="Turn circular records to the strand of the target gene (default cox1) and start "
                                   "them at it; fallback genes when it is missing; linear records are not rotated.")
    add_common(r)

    # ── run ───────────────────────────────────────────────────────────────
    run_p = sub.add_parser(
        "run",
        help="Full pipeline: annotate → reorient → re-annotate → guaranteed final reorientation.")
    add_common(run_p)
    add_annotation(run_p)
    add_engine(run_p)
    add_topology(run_p)
    run_p.add_argument("--if_reannotation_2", choices=["yes", "auto", "no"], default="auto",
                       help="Phase-2 re-annotation of the reoriented sequences. auto (default): round 1 runs in "
                            "full and only samples that fail the check after the cox1 shift (gene at the original "
                            "sequence ends, standard gene missing, gene changed by the shift, gene across the new "
                            "origin) are re-annotated, decisions in reannotation_decision.tsv; yes: every sample; "
                            "no: never.")
    add_metadata(run_p)

    # ── import ────────────────────────────────────────────────────────────
    im = sub.add_parser(
        "import",
        help="Existing MitoZ / MITOS2 results or any annotated GenBank -> the MitoM steps without annotating again: "
             "reorientation, curation and NCBI names, NCBI submission files, genome maps (-i: a result folder, "
             "a .gbf / .gb file, or a folder of them).")
    add_common(im)
    add_annotation(im)
    add_metadata(im)
    im.add_argument("--fasta", action="append", default=[],
                    help="Genome FASTA of MITOS2 results (sequence named as in result.bed; repeatable). Default: "
                         "FASTA files in the MITOS2 result folder or the folder above it.")
    im.add_argument("--topology", choices=["auto", "circular", "linear"], default="auto",
                    help="auto (default): MITOS2 from result.gff, GenBank from the LOCUS line, MitoZ output circular "
                         "(MitoZ writes linear for every record; give linear for contigs).")

    # ── transfer_metadata ─────────────────────────────────────────────────
    tm = sub.add_parser("transfer_metadata",
                        help="Standalone GenBank header transfer between two files/dirs.")
    tm.add_argument("-m", "--metadata", required=True,
                    help="Original GenBank file/dir (metadata source).")
    tm.add_argument("-n", "--new_body", required=True,
                    help="Re-annotated GenBank file/dir (body source).")
    tm.add_argument("-o", "--out", dest="out_root", metavar="OUT", default="transfer_metadata_out")
    old_spelling(tm, "--out_root", "out_root")
    tm.add_argument("--if_replace_metadata",  choices=["yes", "no"], default="yes")
    tm.add_argument("--overwrite_new_body",   choices=["yes", "no"], default="no")
    tm.add_argument("--output_name", default=None,
                    help="Output filename (single-file mode only).")
    tm.add_argument("--smart_metadata",       choices=["yes", "no"], default="no")
    tm.add_argument("--auto_accept_metadata", choices=["yes", "no"], default="no")

    # ── replace_organism ──────────────────────────────────────────────────
    ro = sub.add_parser(
        "replace_organism",
        help=(
            "Restore /organism= and SOURCE text in FEATURES source blocks.\n"
            "Fixes placeholder organism names (e.g. 'Test sp.') that survive\n"
            "transfer_metadata because they live in the body, not the header.\n"
            "In 'run' mode this is applied automatically as Step 3.5."
        ),
    )
    ro.add_argument("-m", "--metadata", required=True,
                    help="Original GenBank file/dir (organism/source name source).")
    ro.add_argument("-n", "--new_body", required=True,
                    help="Re-annotated GenBank file/dir whose /organism= needs fixing.")
    ro.add_argument("-o", "--out", dest="out_root", metavar="OUT", default="replace_organism_out",
                    help="Output folder (default: replace_organism_out).")
    old_spelling(ro, "--out_root", "out_root")
    ro.add_argument("--overwrite_new_body", choices=["yes", "no"], default="no",
                    help="Overwrite the input files in-place (default: no).")

    # ── dedup ─────────────────────────────────────────────────────────────
    dd = sub.add_parser(
        "dedup",
        help="Deduplicate GenBank records (same accession / older version / RefSeq copy of an "
             "INSDC submission / identical sequence; RefSeq NC_ kept) and write GenBank + FASTA.",
        description="Remove duplicate genomes, as `run` does before annotating. The same genome is often in NCBI "
                    "more than once: the authors' record (e.g. MK123456) and the RefSeq copy NCBI made of it "
                    "(NC_...), two versions of one accession (MK123456.1 / .2), or another accession with an "
                    "identical sequence. One record per genome is kept: the NC_ record, the newest version. "
                    "dedup_report.tsv: every input record, kept or dropped, why, and the record kept instead.")
    dd.add_argument("-i", "--input", required=True,
                    help="GenBank and / or FASTA file or folder. FASTA records: identical sequences, and versions of one "
                         "accession when the header's first word is an accession (MK123456.1); other names are never "
                         "merged.")
    dd.add_argument("-o", "--out", dest="out_root", metavar="OUT", default="dedup_out", help="Output folder (default: dedup_out).")
    old_spelling(dd, "--out_root", "out_root")

    # ── gb2fasta ──────────────────────────────────────────────────────────
    g2 = sub.add_parser(
        "gb2fasta",
        help="GenBank (.gb .gbk .gbf .gbff .genbank) -> FASTA, headers with topology= for mitom run.",
        description="Convert GenBank files (NCBI .gb, MitoZ .gbf ...) to FASTA: one <name>.fasta per input file, or "
                    "one <accession>.fasta per record with --split. Headers: >ACCESSION.V definition "
                    "topology=circular|linear - the topology tag keeps a circular genome circular when `mitom run` "
                    "reads several records from one FASTA file.")
    g2.add_argument("-i", "--input", required=True, help="GenBank file or folder (subfolders included).")
    g2.add_argument("-o", "--out", dest="out_root", metavar="OUT", default="gb2fasta_out",
                    help="Output folder (default: gb2fasta_out).")
    g2.add_argument("--split", action="store_true", help="One FASTA file per record (<accession>.fasta).")

    # ── intron ────────────────────────────────────────────────────────────
    it = sub.add_parser(
        "intron", allow_abbrev=False, formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Introns as `run` adds them: mito_intron.py (de novo, MFannot method), then mito_intron_ref.py (intron "
             "sites of annotated relatives) on its output.",
        description="Introns as `run` adds them: mito_intron.py (de novo, MFannot method) -> OUT, then "
                    "mito_intron_ref.py (intron\nsites of annotated relatives) on OUT in place. Reports: "
                    "OUT/intron_report.tsv, OUT/intron_ref_report.tsv.",
        epilog="Other options go to mito_intron.py (python mito_intron.py -h: --lib, --ref_gb, --maxintron, "
               "--timeout ...),\nwith --method ref to mito_intron_ref.py (python mito_intron_ref.py annotate -h).")
    it.add_argument("-i", "--input", required=True, help="GenBank file or folder.")
    it.add_argument("-o", "--out", dest="out_root", metavar="OUT", required=True,
                    help="Output folder (may be the input folder: updated in place).")
    it.add_argument("--method", choices=["both", "denovo", "ref"], default="both",
                    help="both (default): de novo, then the site transfer; denovo: mito_intron.py only (the "
                         "`mitom intron` of MitoM <= 1.0.0); ref: mito_intron_ref.py only.")
    it.add_argument("--genetic_code", type=int, default=5, help="Default 5.")
    it.add_argument("--threads", type=int, default=1, help="Threads for blastp (de novo). Default 1.")
    it.add_argument("--jobs", type=int, default=1, help="GenBank files processed in parallel. Default 1.")
    it.add_argument("--bin", default=None, help="Folder with exonerate, blastp, tblastn, blastn, cmsearch. "
                                                "Default: found automatically.")
    it.add_argument("--lineage", default="",
                    help="NCBI lineage of records without their own taxonomy, for the site transfer "
                         "(e.g. 'Eukaryota;Metazoa;Cnidaria').")

    return p


def error_exit_code(message: str) -> int:
    """Exit code of a refused run: 1 when a program or the reference data is missing (the run cannot be done here),
    2 when the command line or the input is wrong (MitoM's exit codes: 0 done, 1 failure, 2 usage)."""
    return 1 if "not found" in message else 2


def main() -> int:
    parser = build_parser()
    args, extra = parser.parse_known_args()
    if args.cmd == "intron":
        args.extra = extra                       # options of mito_intron.py / mito_intron_ref.py, passed on
    elif extra:
        parser.error("unrecognized arguments: " + " ".join(extra))

    if args.cmd == "import":
        err = prepare_import(args)
        if err:
            print(err, file=sys.stderr)
            return error_exit_code(err)
    if args.cmd in ("annotate", "run", "import"):
        err = resolve_tools(args)
        if err:
            print(err, file=sys.stderr)
            return error_exit_code(err)

    if getattr(args, "suffix_fq", None):
        parse_suffix_fq(args.suffix_fq)

    if (getattr(args, "auto_accept_metadata", "no") == "yes"
            and getattr(args, "smart_metadata", "no") != "yes"):
        print("ERROR: --auto_accept_metadata yes requires --smart_metadata yes", file=sys.stderr)
        return 2


    dispatch = {
        "annotate":          workflow_annotate,
        "reorient":          workflow_reorient,
        "run":               workflow_run,
        "import":            workflow_run,
        "transfer_metadata": workflow_transfer_metadata,
        "replace_organism":  workflow_replace_organism,
        "dedup":             workflow_dedup,
        "gb2fasta":          workflow_gb2fasta,
        "intron":            workflow_intron,
    }
    fn = dispatch.get(args.cmd)
    if fn is None:
        print(f"Unknown command: {args.cmd}", file=sys.stderr)
        return 2
    if getattr(args, "debug", False) or debug_mode():     # --debug (or `mitom --debug <command>`: MITOM_DEBUG=1)
        args.debug = True
        os.environ["MITOM_DEBUG"] = "1"                   # the programs MitoM starts log what they do too
        if hasattr(args, "keep_intermediate"):
            args.keep_intermediate = "yes"
    t0 = time.time()
    try:
        rc = fn(args)
    except Exception as exc:                 # into the run's log with the traceback (mitom debug reads it), then exit 1
        if Steps.current is not None:
            Steps.current.failed(exc)
        logging.exception(f"[ERROR] {type(exc).__name__}: {exc}")
        logging.info(f"[DONE] exit=1 {duration(time.time() - t0)}")
        out = getattr(args, "out_root", None)
        print(f"MitoM stopped on an error: {type(exc).__name__}: {exc}" + (f"\n  diagnosis: mitom debug {out}" if out else ""),
              file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        logging.error(f"[ERROR] interrupted (Ctrl-C) after {duration(time.time() - t0)}")
        raise
    logging.info(f"[DONE] exit={rc} {duration(time.time() - t0)}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
