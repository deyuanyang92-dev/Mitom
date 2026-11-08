"""mitom_commands.py - the smaller commands of mitom_pipeline.py: reorient, transfer_metadata, replace_organism,
dedup, gb2fasta, intron."""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO

from mito_common import (
    GB_SUFFIXES,
    gb_fasta_record,
    iter_gb_files,
    run_tool,
    safe_copy,
    safe_name,
    unique_path,
    write_gb,
    write_nonempty_list,
    write_table,
)
from mitom_inputs import collect_genbank_files, collect_sequence_files, dedup_gb_records, record_format
from mitom_log import setup_root_logger
from mitom_metadata import build_organism_index, replace_organism_in_gb_file, transfer_gb_metadata_file
from mitom_orient import orientation_genes, reorient_one_gb
from mitom_output import RunSummary, counts_from_genbank
from mitom_prepare import log_environment
from mitom_steps import intron_denovo_cmd, intron_ref_cmd


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

    used_samples: set = set()
    for gb_in in gb_files:
        sample = base = safe_name(gb_in.stem)
        k = 2
        while sample in used_samples:                  # a/rec.gb and b/rec.gb: rec and rec_2, not one overwritten
            sample, k = f"{base}_{k}", k + 1
        used_samples.add(sample)
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
    if summary.reorient_failed:                        # (2026-10-31 review: exit 0 when every input failed)
        print(f"mitom reorient: {len(summary.reorient_failed)} of {len(gb_files)} file(s) failed: "
              f"{', '.join(summary.reorient_failed[:5])} (see {out_root / 'pipeline.log'})", file=sys.stderr)
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
    used_names: set = set()
    for bf in body_files:
        matched  = meta_by_stem.get(bf.stem)
        out_file = bf if args.overwrite_new_body == "yes" else unique_path(out_dir / bf.name, used_names)
        try:
            status = _process(bf, matched, out_file)
            if status == "ok":
                n_ok += 1
            else:
                n_skip += 1
        except Exception as exc:
            logging.error(f"[FAIL] {bf}: {exc!r}")      # the reason, not only 'fail' (GLM 1.7, 2026-11-05)
            if args.overwrite_new_body == "no":
                safe_copy(bf, out_file)
            n_fail += 1
            log_rows.append(f"{bf}\t{matched or 'NA'}\t{out_file}\tfail\t0\t0")
    write_table(out_root / "transfer_metadata.tsv",
                         "body_file\tmetadata_file\toutput_file\tstatus\treplaced\tunmatched",
                         log_rows)
    logging.info(f"[SUMMARY] transferred={n_ok}  skipped={n_skip}  failed={n_fail}")
    return 1 if n_fail else 0                         # (N10, 2026-11-05: failures were exit 0)


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

    used_names: set = set()

    for bf in body_files:
        out_file = bf if args.overwrite_new_body == "yes" else unique_path(out_dir / bf.name, used_names)
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
    return 1 if n_fail else 0


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
        write_gb(gb, out_root / "dedup.gb")
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
    used: set = set()
    for f in files:
        recs = []
        for r in SeqIO.parse(str(f), "genbank"):
            try:
                recs.append(gb_fasta_record(r))
            except Exception:                            # noqa: BLE001  (no sequence: CONTIG / WGS record)
                print(f"[gb2fasta] {r.id} ({f.name}): no sequence (CONTIG / WGS record): skipped", file=sys.stderr)
        n += len(recs)
        if args.split:
            for r in recs:
                SeqIO.write([r], str(unique_path(out / f"{safe_name(r.id)}.fasta", used)), "fasta")
        elif recs:                                   # (a.gb in two subfolders: two files, not one overwritten)
            SeqIO.write(recs, str(unique_path(out / f"{f.stem}.fasta", used)), "fasta")
    print(f"[gb2fasta] {n} record(s) from {len(files)} GenBank file(s) -> {out}", file=sys.stderr)
    return 0 if n else 2                              # (0 GenBank records: FASTA or empty input)


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
