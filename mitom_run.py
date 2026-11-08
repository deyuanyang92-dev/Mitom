"""mitom_run.py - the workflows `mitom run` / `import` (workflow_run), `annotate` and the run per clade group
(--clade auto with several phyla / genetic codes)."""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import mitom_groups
from mitom_context import RunContext
from mito_common import reports_dir, write_nonempty_list, write_table
from mitom_inputs import check_transl_tables, iter_inputs, split_or_sanitize_fastas
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
    OUTPUT_LAYOUT,
    RunSummary,
    cleanup_intermediates,
    clear_previous_work,
    finalize_genbank_output,
    merge_gb_files,
    new_reports,
    organize_outputs,
    refresh_summary,
)
from mitom_prepare import log_environment, preflight
from mitom_stage import ENGINE_LABEL, run_annotation_stage, run_phase2
from mitom_steps import (
    run_compare,
    run_curate,
    run_curate_and_submit,
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
    ctx = ctx if ctx is not None else RunContext()
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    removed = clear_previous_work(out_root, keep=("00.input", "00.imported"), protect=_inputs(args))
    new_reports(out_root)
    setup_root_logger(out_root / "pipeline.log", debug=getattr(args, "debug", False))
    _log_removed(out_root, removed)
    log_environment(args)
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))
    preflight(args)

    in_path     = Path(args.input).resolve()
    split_dir   = out_root / "00.internal_fasta"
    gb_conv_dir = out_root / "00.gb_converted_fasta"

    in_files, input_type, acc_to_orig_gb = iter_inputs(in_path, gb_to_fasta_dir=gb_conv_dir,
                                                       dedup_report=reports_dir(out_root) / "dedup_report.tsv",
                                                       stats_out=ctx.dedup_stats)
    logging.info(f"[INPUT_TYPE] {input_type}  ({len(in_files)} file(s))")
    samples, id_map = split_or_sanitize_fastas(
        in_files, split_dir, args.if_split_fasta == "yes", int(args.id_maxlen),
        acc_to_orig_gb=acc_to_orig_gb, topology=args.topology)
    samples = linear_lineages_first(samples, args)
    logging.info(f"[ID_MAP] {id_map}  (total: {len(samples)})")
    if not samples:                                      # as `mitom run` (Codex R7, 2026-11-05: annotate exited 0)
        logging.error(f"ERROR: no sequence to annotate in {in_path} (empty, unreadable or no A/C/G/T; the warnings "
                      f"above and {out_root / 'pipeline.log'} name the file)")
        return 2

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
                orient_notes=ctx.orient,
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
    reorient_after_introns(args, final_dir, ctx.orient)
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
    write_orientation_fallback(out_root, samples, ctx.orient)
    if args.keep_intermediate == "no":
        cleanup_intermediates(out_root, samples, summary.all_problem_samples())
    if samples and len(summary.annotation_failed) == len(samples):
        logging.warning(f"[FAIL] every record failed in the annotation core ({len(samples)}): see reports/failed_logs/")
        return 1
    if samples and not (summary_rows and all_final_gbs):  # (Codex R7, 2026-11-05: exit 0 with no final GenBank)
        logging.warning(f"[FAIL] no final GenBank written for any of the {len(samples)} record(s): see reports/")
        return 1
    return 0


def fallback_src_tag(mode: str, iid: str, reannotated_ids: set) -> str:
    """annotation-summary body_source of a sample whose final body is the shifted round-1 file: 'shifted_round1'
    when --if_reannotation_2 auto decided that round 2 is not needed, 'phase1b_reoriented' when round 2 was
    meant to run but failed (or is off)."""
    return "shifted_round1" if mode == "auto" and iid not in reannotated_ids else "phase1b_reoriented"


def workflow_run(args: argparse.Namespace, ctx: Optional[RunContext] = None) -> int:
    ctx = ctx if ctx is not None else RunContext()
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
                                                       dedup_report=reports_dir(out_root) / "dedup_report.tsv",
                                                       stats_out=ctx.dedup_stats)
    logging.info(f"[INPUT_TYPE] {input_type}  ({len(in_files)} file(s))")
    samples, id_map = split_or_sanitize_fastas(
        in_files, split_dir, args.if_split_fasta == "yes", int(args.id_maxlen),
        acc_to_orig_gb=acc_to_orig_gb, topology=args.topology)
    samples = linear_lineages_first(samples, args)
    total = len(samples)
    logging.info(f"[ID_MAP] {id_map}  (total: {total})")
    if input_type == "converted_from_gb":
        st = ctx.dedup_stats
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
    if total == 0:                                       # unreadable, empty or all-N input: nothing to annotate
        msg = (f"ERROR: no sequence to annotate in {in_path} (empty, unreadable or no A/C/G/T; the warnings above "
               f"and {out_root / 'run.log'} name the file)")
        logging.error(msg)                               # (on the console too)
        return 2

    if ctx.plan is not None:                            # --clade auto: the phylum and code of every record and why
        mitom_groups.write_decision(ctx.plan, reports_dir(out_root) / "clade_decision.tsv")
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
                best_gb, gb_reor, gb_no_gene, orientation_genes(args), args.force_forward, ctx.orient)
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
                        and any(k.startswith(si.internal_id) for k in ctx.orient.fallback))
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
                summary.reannotation_failed.append(si.sample_name)
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
    n_rot = reorient_after_introns(args, final_dir, ctx.orient)
    if n_rot:
        steps.info(f"re-oriented to the first exon : {n_rot} record(s) whose {args.which_gene_first} gained an exon upstream")
    steps.done()

    steps.start("curate", "Boundary rules for protein genes and rRNAs, NCBI gene names")
    curated = run_curate(args, final_dir, out_root)
    merge_dir = curated.gb_dir
    steps.info(f"non-coding regions : {run_ncr(args, merge_dir, out_root).summary}")
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
    write_orientation_fallback(out_root, samples, ctx.orient)
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
    if total and not (n_ok and all_final_gbs):            # (Codex R7, 2026-11-05: exit 0 with no final GenBank)
        logging.warning(f"[FAIL] no final GenBank written for any of the {total} record(s): see reports/")
        return 1
    return 0


def workflow_groups(args: argparse.Namespace, ctx: RunContext) -> int:
    """--clade auto with several phyla / genetic codes: `mitom run` per group, merged (mitom_groups.run_groups)."""
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    setup_root_logger(out_root / "run.log", console_level=logging.WARNING, mode="w",
                      debug=getattr(args, "debug", False))
    log_environment(args)
    logging.info("[PARAMS]\n" + json.dumps(vars(args), indent=2, ensure_ascii=False))
    return mitom_groups.run_groups(args, ctx.plan, sys.argv[1:], say=say)
