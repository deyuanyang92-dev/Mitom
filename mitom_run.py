"""mitom_run.py - the workflows `mitom run` / `import` (workflow_run: the stages of RUN_STAGES, which hand on a
RunState), `annotate` and the run per clade group (--clade auto with several phyla / genetic codes)."""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import mitom_groups
from mito_common import reports_dir, write_nonempty_list, write_table
from mitom_context import RunContext
from mitom_inputs import SampleInfo, check_transl_tables, iter_inputs, split_or_sanitize_fastas
from mitom_log import Steps, duration, record_event, say, setup_root_logger
from mitom_orient import (
    gb_to_internal_fasta,
    orientation_genes,
    orientation_text,
    reorient_after_introns,
    reorient_one_gb,
    select_reannotation,
    write_orientation_fallback,
)
from mitom_output import (
    GROUPED_LAYOUT,
    OUTPUT_LAYOUT,
    RunSummary,
    cleanup_intermediates,
    clear_previous_work,
    finalize_genbank_output,
    merge_gb_files,
    new_reports,
    organize_outputs,
    regroup_outputs,
    where,
    refresh_summary,
)
from mitom_prepare import log_environment, preflight
from mitom_stage import ENGINE_LABEL, run_annotation_stage, run_phase2
from mitom_steps import (
    StepResult,
    run_compare,
    run_curate,
    run_intron_fix,
    run_ncr,
    run_plots,
    run_rescue,
    run_rrna_hmm,
    run_submission,
)
from mitom_topology import apply_topology, linear_lineages_first

KNOWN_LIMITATION = ("[KNOWN_LIMITATION] the annotation core does not annotate introns (group I/II); they are added "
                    "afterwards by mito_intron.py (MFannot method) and mito_intron_ref.py (intron sites of "
                    "annotated relatives; --intron_fix) when exonerate is available - check intron_report.tsv "
                    "and intron_ref_report.tsv.")


def _inputs(args: argparse.Namespace) -> List[str]:
    """-i as given and as read (00.input/ copy): neither may be removed as a previous run's output."""
    return [x for x in (getattr(args, "input", None), getattr(args, "input_given", None)) if x]


def _log_removed(out_root: Path, removed: List[str]) -> None:
    if removed:
        logging.info(f"[CLEAN] removed the previous run's files from {out_root}: {', '.join(removed)}")


def workflow_annotate(args: argparse.Namespace, ctx: Optional[RunContext] = None) -> int:
    """`mitom annotate`: one round of annotation and the same steps after it as `mitom run` (ANNOTATE_STAGES: no
    re-orientation round 2, no comparison, the work layout 04.final_gb/ ... instead of the output layout)."""
    return run_stages(args, ctx, ANNOTATE_STAGES)


def fallback_src_tag(mode: str, iid: str, reannotated_ids: set) -> str:
    """annotation-summary body_source of a sample whose final body is the shifted round-1 file: 'shifted_round1'
    when --if_reannotation_2 auto decided that round 2 is not needed, 'phase1b_reoriented' when round 2 was
    meant to run but failed (or is off)."""
    return "shifted_round1" if mode == "auto" and iid not in reannotated_ids else "phase1b_reoriented"


@dataclass
class RunState:
    """What the stages of one `mitom run` / `import` hand on to the next ones (workflow_run)."""
    args: argparse.Namespace
    ctx: RunContext
    out_root: Path
    steps: Steps
    summary: RunSummary = field(default_factory=RunSummary)
    samples: List[SampleInfo] = field(default_factory=list)
    input_type: str = ""
    total: int = 0
    results1: Dict[str, Dict] = field(default_factory=dict)            # round 1: internal id -> ok, best_gb ...
    reor_fasta_dir: Optional[Path] = None
    reor_ok_ids: List[str] = field(default_factory=list)               # re-oriented: round-2 candidates
    no_gene_iids: Set[str] = field(default_factory=set)                # no orientation gene (or failed): as annotated
    partial_no_gene_gb: Dict[str, Path] = field(default_factory=dict)  # contigs without the gene: merged back later
    do_phase2: bool = False
    results2: Dict[str, Dict] = field(default_factory=dict)
    reannotated_ids: Set[str] = field(default_factory=set)
    final_dir: Optional[Path] = None                                   # 04.final_gb
    summary_rows: List[str] = field(default_factory=list)              # annotation-summary.tsv
    merge_dir: Optional[Path] = None                                   # the final records (curated when --curate yes)
    curated: Optional[StepResult] = None
    all_final_gbs: List[Path] = field(default_factory=list)


def workflow_run(args: argparse.Namespace, ctx: Optional[RunContext] = None) -> int:
    """`mitom run` / `mitom import`: the stages of RUN_STAGES in order (each one numbered step on the terminal); a
    stage that returns an exit code ends the run."""
    return run_stages(args, ctx, RUN_STAGES)


def run_stages(args: argparse.Namespace, ctx: Optional[RunContext], stages) -> int:
    """The stages in order, one numbered step each; a stage that returns an exit code ends the run."""
    r = start_run(args, ctx if ctx is not None else RunContext(), stages)
    for stage in stages:
        rc = stage(r)
        if rc is not None:
            return rc
    return run_exit_code(r)


def start_run(args: argparse.Namespace, ctx: RunContext, stages) -> RunState:
    """Output folder emptied of a previous run's work, run.log, [ENV] / [PARAMS] lines, the summary on the terminal;
    one numbered step per stage of `stages` (run: 11, annotate: 7)."""
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    removed = clear_previous_work(out_root, keep=("00.input", "00.imported"),   # made by this run (input copies, `import`)
                                  protect=_inputs(args))
    new_reports(out_root)
    setup_root_logger(out_root / "run.log", console_level=logging.WARNING, mode="w",   # one run, one log
                      debug=getattr(args, "debug", False))
    _log_removed(out_root, removed)
    log_environment(args)
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))
    preflight(args)
    steps = Steps(len(stages))
    say("")
    say("MitoM - mitochondrial genome annotation")
    say(f"    input        : {Path(args.input).resolve()}")
    say(f"    output       : {out_root}")
    say(f"    clade        : {args.clade}    genetic code: {args.genetic_code}")
    say(f"    parallel     : {args.max_tasks} records at a time x {args.threads} threads")
    say(f"    annotation   : {ENGINE_LABEL[getattr(args, 'engine', 'native')]}")
    if stage_round2 in stages:                            # (`mitom annotate`: one round)
        say(f"    round 2      : {args.if_reannotation_2}  (auto = only records whose annotation is broken by the cox1 shift)")
    if getattr(args, "engine", "native") == "native":
        say(f"    topology     : {args.topology}" + ("  (as labelled, checked at the sequence ends)"
                                                     if args.topology == "auto" else ""))
    say(f"    full log     : {out_root / 'run.log'}")
    return RunState(args=args, ctx=ctx, out_root=out_root, steps=steps)


def stage_input(r: RunState) -> Optional[int]:
    """Input read, GenBank records deduplicated (original annotation saved), one sample per record."""
    args, ctx, out_root, steps = r.args, r.ctx, r.out_root, r.steps
    in_path     = Path(args.input).resolve()
    split_dir   = out_root / "00.internal_fasta"
    gb_conv_dir = out_root / "00.gb_converted_fasta"

    steps.start("input", "Read input and remove duplicate records")
    in_files, r.input_type, acc_to_orig_gb = iter_inputs(in_path, gb_to_fasta_dir=gb_conv_dir,
                                                         orig_dir=out_root / "original_annotation",
                                                         dedup_report=reports_dir(out_root) / "dedup_report.tsv",
                                                         stats_out=ctx.dedup_stats)
    logging.info(f"[INPUT_TYPE] {r.input_type}  ({len(in_files)} file(s))")
    samples, id_map = split_or_sanitize_fastas(
        in_files, split_dir, args.if_split_fasta == "yes", int(args.id_maxlen),
        acc_to_orig_gb=acc_to_orig_gb, topology=args.topology)
    r.samples = linear_lineages_first(samples, args)
    total = r.total = len(r.samples)
    logging.info(f"[ID_MAP] {id_map}  (total: {total})")
    if r.input_type == "converted_from_gb":
        st = ctx.dedup_stats
        steps.info(f"GenBank records read: {st.get('total_in', '?')}")
        dropped = [(st.get(k, 0), why) for k, why in (
            ("dup_refseq", "GenBank record that a RefSeq NC_ record was made from (NC_ kept)"),
            ("dup_version", "older version of the same accession"), ("dup_exact", "same accession twice"),
            ("dup_seq", "identical sequence"), ("no_sequence", "no sequence"))]
        for n, why in dropped:
            if n:
                steps.info(f"removed {n}: {why}")
        steps.info("the original annotation of every kept record is saved for comparison "
                   f"({where(args, 'original_annotation/')})")
        steps.done(f"{total} records to annotate; reasons per record: {where(args, 'reports/dedup_report.tsv')}")
    else:
        steps.done(f"{total} FASTA sequences to annotate")
    if total == 0:                                       # unreadable, empty or all-N input: nothing to annotate
        msg = (f"ERROR: no sequence to annotate in {in_path} (empty, unreadable or no A/C/G/T; the warnings above "
               f"and {out_root / 'run.log'} name the file)")
        logging.error(msg)                               # (on the console too)
        return 2

    if ctx.plan is not None:                            # --clade auto: the phylum and code of every record and why
        mitom_groups.write_decision(ctx.plan, reports_dir(out_root) / "clade_decision.tsv")
    check_transl_tables(out_root / "original_annotation", int(args.genetic_code))
    return None


def stage_round1(r: RunState) -> Optional[int]:
    """Round-1 annotation of every sample (the engine), topology checked at the sequence ends."""
    args, out_root, steps, summary = r.args, r.out_root, r.steps, r.summary
    steps.start("annotate", f"Round-1 annotation ({ENGINE_LABEL[getattr(args, 'engine', 'native')]}: protein genes, rRNAs, tRNAs)")
    results1 = run_annotation_stage(
        samples=r.samples, phase_dirname="01.anno_round1",
        out_root=out_root, args=args,
        fasta_selector=lambda si: si.internal_fasta)
    for si in r.samples:
        if not results1.get(si.internal_id, {}).get("ok"):
            summary.annotation_failed.append(si.sample_name)
    r.samples, r.results1 = apply_topology(r.samples, results1, args, out_root, reports_dir(out_root), steps.info)
    steps.done(f"{r.total - len(summary.annotation_failed)} annotated, {len(summary.annotation_failed)} failed"
               + (f" ({', '.join(summary.annotation_failed[:5])}{' ...' if len(summary.annotation_failed) > 5 else ''})"
                  if summary.annotation_failed else ""))
    return None


def stage_reorient(r: RunState) -> Optional[int]:
    """Re-orientation: plus strand, the sequence starts with --which_gene_first (else a fallback gene)."""
    args, ctx, out_root, steps, summary = r.args, r.ctx, r.out_root, r.steps, r.summary
    steps.start("reorient", f"Re-orient: plus strand, sequence starts with {args.which_gene_first} "
                f"(if missing: {args.fallback_genes})")

    logging.info(
        f"[PHASE1b] reorient gene='{args.which_gene_first}' force_forward={args.force_forward}")
    reor_fasta_dir = r.reor_fasta_dir = out_root / "00.reoriented_fasta_internal"
    reor_fasta_dir.mkdir(parents=True, exist_ok=True)

    for si in r.samples:
        r1      = r.results1.get(si.internal_id, {})
        best_gb = Path(r1["best_gb"]) if r1.get("best_gb") else None
        if not best_gb or not best_gb.exists():
            continue
        reor_dir   = out_root / "01.anno_round1" / si.internal_id / "reoriented"
        reor_dir.mkdir(parents=True, exist_ok=True)
        gb_reor    = reor_dir / f"{si.internal_id}.reoriented.gbf"
        gb_no_gene = reor_dir / f"{si.internal_id}.no_{args.which_gene_first}.gbf"
        try:
            n_ok, n_no, no_gene_ids = reorient_one_gb(
                best_gb, gb_reor, gb_no_gene, orientation_genes(args), args.force_forward, ctx.orient)
            if n_ok == 0 and n_no > 0:
                record_event(si.sample_name, f"no orientation gene ({args.which_gene_first}, {args.fallback_genes}) "
                             f"in its {n_no} sequence(s): kept as annotated, no round 2 ({si.internal_id})")
                summary.no_target_gene.append(si.sample_name)
                r.no_gene_iids.add(si.internal_id)
            elif n_no > 0:
                logging.info(
                    f"[PARTIAL_NO_GENE] '{si.sample_name}' ({si.internal_id}): "
                    f"{n_no} minor contig(s) lack '{args.which_gene_first}'; "
                    f"primary reoriented OK — proceeding to Phase-2.")
                if gb_no_gene.exists():
                    r.partial_no_gene_gb[si.internal_id] = gb_no_gene
            if n_ok > 0 and gb_reor.exists():
                r.reor_ok_ids.append(si.internal_id)
                gb_to_internal_fasta(gb_reor, reor_fasta_dir / f"{si.internal_id}.fasta",
                                     si.internal_id)
        except Exception as exc:
            record_event(si.sample_name, f"re-orientation failed ({si.internal_id}): {exc}", logging.ERROR)
            summary.reorient_failed.append(si.sample_name)      # an error, not a record without the gene
            r.no_gene_iids.add(si.internal_id)

    fallback_used = sum(1 for si in r.samples if si.internal_id in r.reor_ok_ids
                        and any(k.startswith(si.internal_id) for k in ctx.orient.fallback))
    steps.done(f"{len(r.reor_ok_ids)} re-oriented"
               + (f" ({fallback_used} by a fallback gene)" if fallback_used else "")
               + (f", {len(r.no_gene_iids)} without any orientation gene (kept as annotated)" if r.no_gene_iids else ""))
    return None


def stage_round2(r: RunState) -> Optional[int]:
    """Round 2: the re-oriented records annotated again where the shift broke the annotation (--if_reannotation_2)."""
    args, out_root, steps, summary = r.args, r.out_root, r.steps, r.summary
    steps.start("round2", "Check the annotation after re-orientation; annotate again where it is offset (round 2)")
    r.do_phase2 = args.if_reannotation_2 in ("yes", "auto")
    if not r.do_phase2:
        steps.done("--if_reannotation_2 no: skipped")
        return None
    phase2_samples = [si for si in r.samples if si.internal_id in r.reor_ok_ids]
    n_checked = len(phase2_samples)
    if args.if_reannotation_2 == "auto":
        phase2_samples = select_reannotation(phase2_samples, r.results1, out_root,
                                             int(getattr(args, "genetic_code", None) or 5))
        steps.info(f"checked {n_checked}: {len(phase2_samples)} need round 2 (gene cut at the old ends, gene "
                   f"changed by the shift, gene across the new start, standard gene missing); "
                   f"{n_checked - len(phase2_samples)} keep the shifted round-1 annotation")
        steps.info(f"decision and reason per record: {where(args, 'reports/reannotation_decision.tsv')}")
    else:
        steps.info(f"--if_reannotation_2 yes: all {n_checked} re-oriented records are annotated again")
    r.results2 = run_phase2(phase2_samples, args, out_root, r.reor_fasta_dir)
    r.reannotated_ids = {si.internal_id for si in phase2_samples}
    n_fail2 = 0
    for si in phase2_samples:
        if not r.results2.get(si.internal_id, {}).get("ok"):
            n_fail2 += 1
            summary.reannotation_failed.append(si.sample_name)
            record_event(si.sample_name, "round-2 annotation failed: the shifted round-1 annotation is kept")
    steps.done(f"round 2: {len(phase2_samples) - n_fail2} annotated"
               + (f", {n_fail2} failed (they keep the shifted round-1 annotation)" if n_fail2 else ""))
    return None


def final_body(r: RunState, si: SampleInfo) -> Tuple[Optional[Path], str]:
    """The annotation a sample's final record is written from and the body_source of annotation-summary.tsv: round 2,
    else the shifted round 1, else round 1 as annotated (no orientation gene); (None, ...) when there is none."""
    iid = si.internal_id
    r1 = r.results1.get(iid, {})
    chosen_gb: Optional[Path] = None
    src_tag   = "phase1_raw"

    if r.do_phase2 and iid not in r.no_gene_iids:
        r2 = r.results2.get(iid, {})
        if r2.get("best_gb") and Path(r2["best_gb"]).exists():
            chosen_gb = Path(r2["best_gb"])
            src_tag   = "phase2"

    if chosen_gb is None:
        gb_reor = (r.out_root / "01.anno_round1" / iid / "reoriented"
                   / f"{iid}.reoriented.gbf")
        if gb_reor.exists():
            chosen_gb = gb_reor
            src_tag   = fallback_src_tag(r.args.if_reannotation_2 if r.do_phase2 else "no", iid, r.reannotated_ids)

    if chosen_gb is None and r1.get("best_gb") and Path(r1["best_gb"]).exists():
        chosen_gb = Path(r1["best_gb"])
    return chosen_gb, src_tag


def with_partial_contigs(r: RunState, si: SampleInfo, chosen_gb: Path, src_tag: str) -> Tuple[Path, str]:
    """The contigs of a sample that lack the orientation gene merged back into its re-oriented body (not into the
    round-1 file, which still holds them)."""
    iid = si.internal_id
    r1 = r.results1.get(iid, {})
    extra_no_gene_gb = r.partial_no_gene_gb.get(iid)
    if not (extra_no_gene_gb and extra_no_gene_gb.exists()):
        return chosen_gb, src_tag
    chosen_is_phase1_raw = False
    if r1.get("best_gb"):
        try:
            chosen_is_phase1_raw = chosen_gb.resolve() == Path(r1["best_gb"]).resolve()
        except Exception:
            chosen_is_phase1_raw = False
    if chosen_is_phase1_raw:
        return chosen_gb, src_tag
    merged_body_dir = r.out_root / "04.final_body_merged"
    merged_body_dir.mkdir(parents=True, exist_ok=True)
    merged_body = merged_body_dir / f"{iid}.with_no_gene.gbf"
    merge_gb_files([chosen_gb, extra_no_gene_gb], merged_body)
    logging.info(
        f"[FINAL_PARTIAL_NO_GENE_MERGE] '{si.sample_name}': merged "
        f"{extra_no_gene_gb.name} back into final body.")
    return merged_body, src_tag + "+partial_no_gene_merged"


def stage_write_final(r: RunState) -> Optional[int]:
    """Final GenBank per sample (04.final_gb/): re-oriented once more, the original header and /organism restored."""
    args, ctx, out_root, steps, summary = r.args, r.ctx, r.out_root, r.steps, r.summary
    steps.start("write", "Write the annotation; copy the original GenBank header (organism, references ...) back")

    final_dir         = r.final_dir = out_root / "04.final_gb"
    no_gene_final_dir = final_dir / f"no_{args.which_gene_first}_gene"
    final_dir.mkdir(parents=True, exist_ok=True)
    no_gene_final_dir.mkdir(parents=True, exist_ok=True)

    summary_rows:    List[str] = r.summary_rows
    meta_log_rows:   List[str] = []

    for si in r.samples:
        sample      = si.sample_name
        iid         = si.internal_id
        metadata_src = si.orig_gb_path
        chosen_gb, src_tag = final_body(r, si)
        if chosen_gb is None or not chosen_gb.exists():
            if sample not in summary.annotation_failed:
                summary.final_write_failed.append(sample)
            continue
        body_for_final, src_tag = with_partial_contigs(r, si, chosen_gb, src_tag)
        dest = (no_gene_final_dir if iid in r.no_gene_iids else final_dir) / f"{sample}.gbf"

        try:
            pcg, trna, rrna, meta_mode, orient_ok = finalize_genbank_output(
                sample_name=sample,
                body_src_gb=body_for_final,                          # with the contigs without the gene
                metadata_src_gb=metadata_src,
                dest_gb=dest,
                which_gene_first=orientation_genes(args),
                force_forward=args.force_forward,
                preserve_metadata=(args.preserve_phase1_metadata == "yes"),
                smart_metadata=(getattr(args, "smart_metadata", "no") == "yes"),     # (not options of annotate)
                auto_accept=(getattr(args, "auto_accept_metadata", "no") == "yes"),
                orient_notes=ctx.orient,
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
    write_problem_lists(r)
    (reports_dir(out_root) / "run_params.json").write_text(
        json.dumps(vars(args), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return None


def write_problem_lists(r: RunState) -> None:
    """reports/<problem>.txt: the samples of each problem (only problems that occurred)."""
    reports, summary, gene = reports_dir(r.out_root), r.summary, r.args.which_gene_first
    write_nonempty_list(reports / "annotation_failed.txt",
                        summary.annotation_failed, "# annotation core returned non-zero exit")
    write_nonempty_list(reports / f"no_{gene}_samples.txt",
                        summary.no_target_gene,
                        f"# samples where '{gene}' was not found")
    write_nonempty_list(reports / "final_write_failed.txt",
                        summary.final_write_failed, "# failed during final GB write")
    write_nonempty_list(reports / "reorient_failed.txt",
                        summary.reorient_failed, "# failed during reorientation")
    write_nonempty_list(reports / "metadata_transfer_failed.txt",
                        summary.metadata_transfer_failed, "# metadata transfer failed")
    write_nonempty_list(reports / "metadata_key_mismatch.txt",
                        summary.metadata_key_mismatch,
                        "# accession key mismatch during metadata transfer")
    write_nonempty_list(reports / "all_problem_samples.txt",
                        summary.all_problem_samples(), "# all problematic samples")


def stage_add_missing(r: RunState) -> Optional[int]:
    """What the core missed: MITOS2 rescue (optional), 12S/16S by rRNA HMM, introns; re-oriented to the first exon."""
    args, steps = r.args, r.steps
    steps.start("add_missing", "Add what the core missed: MITOS2 rescue (optional), 12S/16S by rRNA HMM, introns")
    for name, step in (("rescue", run_rescue), ("rrna_hmm", run_rrna_hmm), ("introns", run_intron_fix)):
        steps.info(steps.part(name, step, args, r.final_dir, r.out_root).summary)
    n_rot = reorient_after_introns(args, r.final_dir, r.ctx.orient)
    if n_rot:
        steps.info(f"re-oriented to the first exon : {n_rot} record(s) whose {args.which_gene_first} gained an exon upstream")
    steps.done()
    return None


def stage_curate(r: RunState) -> Optional[int]:
    """Boundary rules and NCBI gene names (mito_curate.py -> 05.curated/), non-coding regions."""
    steps = r.steps
    steps.start("curate", "Boundary rules for protein genes and rRNAs, NCBI gene names")
    r.curated = run_curate(r.args, r.final_dir, r.out_root)
    r.merge_dir = r.curated.gb_dir
    steps.info(f"non-coding regions : {run_ncr(r.args, r.merge_dir, r.out_root).summary}")
    steps.done(r.curated.summary, r.curated.ok)
    return None


def stage_ncbi(r: RunState) -> Optional[int]:
    """NCBI check (table2asn) and submission files; summary.tsv; all records in one GenBank file."""
    args, out_root, steps, merge_dir = r.args, r.out_root, r.steps, r.merge_dir
    steps.start("ncbi", "NCBI check (table2asn) and submission files")
    submitted = run_submission(args, merge_dir, out_root) if merge_dir != r.final_dir else r.curated
    steps.done(submitted.summary, submitted.ok)
    refresh_summary(out_root, merge_dir)
    r.all_final_gbs = sorted(                    # rglob: the no_<gene>_gene/ records are merged too
        p for p in merge_dir.rglob("*.gbf")
        if p.is_file() and p.stat().st_size > 0
    )
    if r.all_final_gbs:
        n_merged = merge_gb_files(r.all_final_gbs, out_root / "all.final.gbf")
        logging.info(f"[MERGED_ALL] {out_root / 'all.final.gbf'}  ({n_merged} records)")
    else:
        logging.warning("[MERGED_ALL] No final GB files to merge.")

    logging.info(f"[FINAL_DIR] {r.final_dir}")
    logging.info(f"[OUTPUT_COUNT] {len(r.summary_rows)}/{r.total}")
    logging.info(f"[METADATA_MODE] preserve={args.preserve_phase1_metadata}")
    r.summary.log()
    logging.info(KNOWN_LIMITATION)
    write_orientation_fallback(out_root, r.samples, r.ctx.orient)
    if args.keep_intermediate == "no":
        cleanup_intermediates(out_root, r.samples, r.summary.all_problem_samples())
    return None


def stage_compare(r: RunState) -> Optional[int]:
    """Original annotation vs MitoM, gene by gene (GenBank input)."""
    r.steps.start("compare", "Compare with the original annotation")
    compared = run_compare(r.args, r.merge_dir, r.out_root, original=r.out_root / "original_annotation"
                           if r.input_type == "converted_from_gb" else None)
    r.steps.done(compared.summary, compared.ok)
    return None


def stage_maps(r: RunState) -> Optional[int]:
    """Circular genome maps (--plot)."""
    r.steps.start("maps", "Genome maps (circular plots)")
    maps = run_plots(r.args, r.merge_dir, r.out_root)
    r.steps.done(maps.summary, maps.ok)
    return None


def stage_output(r: RunState) -> Optional[int]:
    """The output layout (OUTPUT_LAYOUT, README.txt)."""
    args, out_root, steps = r.args, r.out_root, r.steps
    steps.start("output", "Organise the output")
    organize_outputs(out_root, r.merge_dir, keep=args.keep_intermediate == "yes",
                     reports_named=(f"no_{args.which_gene_first}_samples.txt",))
    grouped = getattr(args, "output_layout", "grouped") == "grouped"
    if grouped:
        regroup_outputs(out_root)
    for name, what in (GROUPED_LAYOUT if grouped else OUTPUT_LAYOUT):
        if (out_root / name).exists():
            steps.info(f"{name:28s} {what}")
    steps.done(f"explained in {out_root / 'README.txt'}")
    say("")
    say(f"MitoM finished in {duration(time.time() - steps.run_t0)}: {out_root}")
    return None


def run_exit_code(r: RunState) -> int:
    """1 when every record failed in the annotation core or no final GenBank was written, else 0."""
    total = r.total
    if total and len(r.summary.annotation_failed) == total:
        logging.warning(f"[FAIL] every record failed in the annotation core ({total}): see "
                        f"{where(r.args, 'reports/failed_logs/')}")
        return 1
    if total and not (r.summary_rows and r.all_final_gbs):  # (Codex R7, 2026-11-05: exit 0 with no final GenBank)
        logging.warning(f"[FAIL] no final GenBank written for any of the {total} record(s): see {where(r.args, 'reports/')}")
        return 1
    return 0


# `mitom run` / `import` in order; one numbered step each on the terminal (run_stages: Steps(len(stages)))
RUN_STAGES = (stage_input, stage_round1, stage_reorient, stage_round2, stage_write_final, stage_add_missing,
              stage_curate, stage_ncbi, stage_compare, stage_maps, stage_output)
# `mitom annotate`: one round, the steps after it as in `run`, the work layout (no comparison, no output layout)
ANNOTATE_STAGES = (stage_input, stage_round1, stage_write_final, stage_add_missing, stage_curate, stage_ncbi,
                   stage_maps)


def workflow_groups(args: argparse.Namespace, ctx: RunContext) -> int:
    """--clade auto with several phyla / genetic codes: `mitom run` per group, merged (mitom_groups.run_groups)."""
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "run.log", console_level=logging.WARNING, mode="w",
                      debug=getattr(args, "debug", False))
    log_environment(args)
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))
    return mitom_groups.run_groups(args, ctx.plan, sys.argv[1:], say=say)
