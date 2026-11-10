"""mitom_output.py - the final GenBank of each sample (internal ID -> sample name, orientation, original header and
/organism, gene counts; written as a draft that replaces the final file at the end) and the output folder of `run`:
problem lists, annotation summary, summary.tsv, the layout of OUTPUT_LAYOUT with README.txt, the intermediate folders
kept (work/) or removed."""
from __future__ import annotations

import logging
import os
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from mito_common import REPORTS, reports_dir, safe_copy, write_gb
from mitom_log import record_event
from mitom_inputs import SampleInfo
from mitom_metadata import (_sync_organism_qualifier, build_organism_index, gb_file_has_meaningful_metadata,
                            replace_organism_in_gb_file, transfer_gb_metadata_file)
from mitom_orient import OrientNotes, reorient_one_gb, starts_with_gene


@dataclass
class RunSummary:
    annotation_failed:        List[str] = field(default_factory=list)
    final_write_failed:       List[str] = field(default_factory=list)
    no_target_gene:           List[str] = field(default_factory=list)
    metadata_transferred:     List[str] = field(default_factory=list)
    metadata_compat_match:    List[str] = field(default_factory=list)
    metadata_key_mismatch:    List[str] = field(default_factory=list)
    metadata_skipped:         List[str] = field(default_factory=list)
    metadata_transfer_failed: List[str] = field(default_factory=list)
    reorient_failed:          List[str] = field(default_factory=list)       # reorient workflow (nothing annotated)
    reannotation_failed:      List[str] = field(default_factory=list)       # round 2: logs kept (GLM 4.1, 2026-11-05)

    def all_problem_samples(self) -> List[str]:
        return sorted(set(
            self.annotation_failed + self.final_write_failed
            + self.no_target_gene  + self.metadata_transfer_failed
            + self.metadata_key_mismatch                           # header matched only by a compatible accession
            + self.reorient_failed + self.reannotation_failed
        ))

    def log(self) -> None:
        logging.info(
            "[SUMMARY] annotation_failed=%d  final_write_failed=%d  no_target_gene=%d"
            "  reorient_failed=%d",
            len(self.annotation_failed), len(self.final_write_failed), len(self.no_target_gene),
            len(self.reorient_failed))
        logging.info(
            "[SUMMARY] metadata_transferred=%d  compat_match=%d  key_mismatch=%d  "
            "skipped=%d  transfer_failed=%d",
            len(self.metadata_transferred), len(self.metadata_compat_match),
            len(self.metadata_key_mismatch), len(self.metadata_skipped),
            len(self.metadata_transfer_failed))


def merge_gb_files(gb_files: List[Path], merged_out: Path) -> int:
    """The files' text one after the other (each checked by parsing it): a Biopython round trip would write the
    value-less /pseudo of mito_curate P10 as /pseudo=""."""
    merged_out.parent.mkdir(parents=True, exist_ok=True)
    n, chunks = 0, []
    for gb in gb_files:
        try:
            k = sum(1 for _ in SeqIO.parse(str(gb), "genbank"))
        except Exception as exc:
            logging.warning(f"[MERGE_SKIP] {gb}: {exc}")
            continue
        if k:
            text = gb.read_text()
            chunks.append(text if text.endswith("\n") else text + "\n")
            n += k
    if chunks:                                         # a temporary file, then replace: never a half-written merge
        tmp = merged_out.with_name(f".{merged_out.name}.tmp")     # (Grok 3.12, 2026-11-05)
        tmp.write_text("".join(chunks))
        os.replace(str(tmp), str(merged_out))
    return n


def refresh_summary(out_root: Path, gb_dir: Path) -> None:
    """Recount PCG / tRNA / rRNA of annotation-summary.tsv from the final GenBank files (after rescue, rRNA
    HMM fallback and curation); rows without a final file are kept unchanged."""
    f = out_root / "annotation-summary.tsv"
    if not f.exists():
        return
    finals = {p.stem: p for p in Path(gb_dir).rglob("*.gbf")}
    lines = f.read_text(encoding="utf-8").splitlines()
    out = lines[:1]
    for line in lines[1:]:
        c = line.split("\t")
        if c and c[0] in finals:
            c[1:4] = [str(x) for x in counts_from_genbank(finals[c[0]])]
        out.append("\t".join(c))
    f.write_text("\n".join(out) + "\n", encoding="utf-8")


def counts_from_genbank(gb_path: Path) -> Tuple[int, int, int]:
    pcg = trna = rrna = 0
    for rec in SeqIO.parse(str(gb_path), "genbank"):
        for feat in rec.features:
            if feat.type == "CDS":
                pcg  += 1
            elif feat.type == "tRNA":
                trna += 1
            elif feat.type == "rRNA":
                rrna += 1
    return pcg, trna, rrna


_LOCUS_NAME_MAXLEN = 16


def rewrite_genbank_record_ids(gb_in: Path, gb_out: Path, new_base_id: str) -> Tuple[int, int, int]:
    recs = list(SeqIO.parse(str(gb_in), "genbank"))
    if not recs:
        raise ValueError(f"No records in {gb_in}")
    pcg = trna = rrna = 0
    out_recs: List["SeqRecord"] = []
    for i, r in enumerate(recs, start=1):
        rid = new_base_id if i == 1 else f"{new_base_id}_{i}"
        locus = rid[:_LOCUS_NAME_MAXLEN]                  # a new record, not r[:] (which copies every feature)
        if len(rid) > _LOCUS_NAME_MAXLEN:
            logging.info(
                f"[LOCUS_TRUNCATE] '{rid}' ({len(rid)} chars) truncated to '{locus}' "
                f"in LOCUS line (GenBank limit={_LOCUS_NAME_MAXLEN}). "
                f"ACCESSION/VERSION unaffected.")
        new_annotations = dict(r.annotations)
        new_annotations["accessions"] = [rid]
        new_annotations.pop("sequence_version", None)
        rr = SeqRecord(
            r.seq,
            id=rid,
            name=locus,
            description=r.description or "",
            dbxrefs=list(getattr(r, "dbxrefs", [])),
            annotations=new_annotations,
            features=list(r.features),
        )
        out_recs.append(rr)
        for feat in r.features:
            if feat.type == "CDS":
                pcg  += 1
            elif feat.type == "tRNA":
                trna += 1
            elif feat.type == "rRNA":
                rrna += 1
    gb_out.parent.mkdir(parents=True, exist_ok=True)
    write_gb(out_recs, gb_out)                       # /pseudo stays value-less (GLM 1.9, 2026-11-05)
    return pcg, trna, rrna


def finalize_genbank_output(
    *,
    sample_name:       str,
    body_src_gb:       Path,
    metadata_src_gb:   Optional[Path],
    dest_gb:           Path,
    which_gene_first:  str,
    force_forward:     bool,
    preserve_metadata: bool,
    smart_metadata:    bool = False,
    auto_accept:       bool = False,
    orient_notes:      Optional[OrientNotes] = None,
) -> Tuple[int, int, int, str, bool]:
    """Produce the final, publication-grade GenBank file for one sample.

    Steps
    -----
    1. ID rewrite      — replace internal S0xxxx IDs with sample_name.
    2. Reorientation   — guarantee the record starts at which_gene_first.
    3. Metadata overlay — (when preserve_metadata) replace GB header section.
    3.5 Organism sync  — fix /organism= inside the FEATURES source feature.
    4. Recount         — recount PCG / tRNA / rRNA from the definitive file.

    Every step writes a draft next to dest_gb that replaces it only at the end: a failure never leaves a half-final
    file at the final path.
    """
    pcg = trna = rrna = 0
    orientation_ok = False

    td = Path(tempfile.mkdtemp(prefix="mitom_final_"))
    dest_gb.parent.mkdir(parents=True, exist_ok=True)
    draft_gb = dest_gb.with_name(f".{dest_gb.name}.tmp")
    committed = False
    try:
        if draft_gb.exists():
            draft_gb.unlink()

        # ── Step 1: ID rewrite ────────────────────────────────────────
        id_gb = td / "id_rewrite.gbf"
        try:
            pcg, trna, rrna = rewrite_genbank_record_ids(body_src_gb, id_gb, sample_name)
        except Exception as exc:
            raise RuntimeError(f"ID rewrite failed for {sample_name}: {exc}") from exc

        # ── Step 2: Final guaranteed reorientation ────────────────────
        reor_gb    = td / "reoriented.gbf"
        no_gene_gb = td / "no_gene.gbf"
        try:
            n_ok, n_no, _ = reorient_one_gb(
                id_gb, reor_gb, no_gene_gb, which_gene_first, force_forward, orient_notes)
            if n_ok > 0 and reor_gb.exists():
                orientation_ok = True
                if n_no > 0 and no_gene_gb.exists():
                    merged_reor = td / "merged_reoriented.gbf"
                    merge_gb_files([reor_gb, no_gene_gb], merged_reor)
                    working_gb = merged_reor
                    logging.info(
                        f"[FINAL_MULTICONTIG] {sample_name}: merged {n_ok} reoriented "
                        f"+ {n_no} un-reoriented contig(s) into final GB.")
                else:
                    working_gb = reor_gb
            else:
                working_gb = id_gb
                logging.info(
                    f"[FINAL_NO_GENE] {sample_name}: '{which_gene_first}' not found "
                    f"in chosen_gb; orientation unchanged.")
        except Exception as exc:
            record_event(sample_name, f"final re-orientation failed ({exc}): written as annotated")
            working_gb = id_gb

        # ── Step 3: Optional metadata overlay (header section, text-level) ──
        metadata_mode = "rewritten_no_metadata"

        if preserve_metadata and metadata_src_gb and metadata_src_gb.exists():
            if gb_file_has_meaningful_metadata(metadata_src_gb):
                meta_out = td / "with_metadata.gbf"
                ok, n_replaced, n_unmatched, unmatched, reason, compat = transfer_gb_metadata_file(
                    metadata_src_gb=metadata_src_gb,
                    new_body_gb=working_gb,
                    output_gb=meta_out,
                    smart=smart_metadata,
                    auto_accept=auto_accept,
                )
                if ok and n_replaced > 0:
                    safe_copy(meta_out, draft_gb)
                    metadata_mode = "metadata_transferred"
                    if compat:
                        metadata_mode += "+compat"
                    if unmatched:
                        logging.warning(
                            f"[META_UNMATCHED] {sample_name}: {len(unmatched)} record(s) "
                            f"unmatched: {unmatched}")
                elif ok and n_replaced == 0:
                    safe_copy(working_gb, draft_gb)
                    metadata_mode = "metadata_no_match"
                    logging.info(
                        f"[META_NO_MATCH] {sample_name}: 0 records matched in "
                        f"{metadata_src_gb.name}; using ID-rewritten GB.")
                else:
                    safe_copy(working_gb, draft_gb)
                    metadata_mode = "metadata_failed"
                    logging.warning(f"[META_FAIL] {sample_name}: {reason}")
            else:
                safe_copy(working_gb, draft_gb)
                metadata_mode = "metadata_src_not_usable"
        else:
            safe_copy(working_gb, draft_gb)

        # ── Step 3.5: Sync /organism= in FEATURES source ──
        try:
            org_fixed = td / "org_fixed.gbf"
            applied = False

            if metadata_src_gb and metadata_src_gb.exists():
                org_index = build_organism_index([metadata_src_gb])
                if org_index:
                    n_rep, n_nf, n_tot = replace_organism_in_gb_file(draft_gb, org_fixed, org_index)
                    if n_rep > 0:
                        safe_copy(org_fixed, draft_gb)
                        applied = True
                        logging.info(
                            f"[ORG_SYNC] {sample_name}: /organism= restored in "
                            f"{n_rep}/{n_tot} record(s) from original GB"
                            + (f" ({n_nf} not matched)" if n_nf else ""))

            if not applied:
                n_synced = _sync_organism_qualifier(draft_gb, org_fixed)
                if n_synced > 0:
                    safe_copy(org_fixed, draft_gb)
                    logging.info(
                        f"[ORG_SYNC] {sample_name}: /organism= synced from header "
                        f"in {n_synced} record(s)")

        except Exception as exc:
            logging.warning(
                f"[ORG_SYNC_FAIL] {sample_name}: organism qualifier sync failed: {exc}")

        # ── Step 4: recount from definitive file ──────────────────────
        try:
            pcg, trna, rrna = counts_from_genbank(draft_gb)
            # cox1_first from the file itself: a linear record is never rotated (rotate_record), whatever
            # reorient_one_gb counted
            orientation_ok = starts_with_gene(draft_gb, which_gene_first, force_forward)
        except Exception as exc:             # the counts of the summary table stay those before this step (N4)
            logging.warning(f"[FINAL] {sample_name}: genes not recounted from {dest_gb.name}: {exc}")
        os.replace(str(draft_gb), str(dest_gb))
        committed = True

    finally:
        if not committed:
            try:
                draft_gb.unlink()
            except Exception:
                pass
        shutil.rmtree(td, ignore_errors=True)

    return pcg, trna, rrna, metadata_mode, orientation_ok


INTERMEDIATE_DIRS = (
    "00.internal_fasta", "00.gb_converted_fasta", "00.reoriented_fasta_internal",
    "01.anno_round1", "02.anno_topology", "03.anno_round2_reoriented", "04.final_body_merged", "00.imported",
    "00.input",
)


def cleanup_intermediates(out_root: Path, samples: List[SampleInfo],
                          problem_samples: List[str]) -> None:
    """Delete internal-ID intermediates once final output is written.

    Kept: 04.final_gb/, all.final.gbf, reports/, run.log. The core logs of problem samples are copied to
    reports/failed_logs/ before deletion. `run` then arranges what is kept into its final layout (organize_outputs).
    """
    import mitom_log
    problems = set(problem_samples) | mitom_log.RECORD_NAMES     # 02 / round-2 failures too (Kimi 4.1)
    for si in samples:
        if si.sample_name not in problems:
            continue
        for phase in ("01.anno_round1", "02.anno_topology", "03.anno_round2_reoriented"):
            wd = out_root / phase / si.internal_id
            for log in [wd / "sample.log", *sorted((wd / "logs").glob("*.log"))]:
                if log.is_file():
                    dst = reports_dir(out_root) / "failed_logs" / f"{si.sample_name}.{phase}.{log.name}"
                    dst.parent.mkdir(exist_ok=True)
                    shutil.copy2(log, dst)
    freed = 0
    for name in INTERMEDIATE_DIRS:
        d = out_root / name
        if d.is_dir():
            freed += sum(p.stat().st_size for p in d.rglob("*") if p.is_file())
            shutil.rmtree(d, ignore_errors=True)
    logging.info(f"[CLEANUP] intermediates removed, {freed / 1e6:.1f} MB freed "
                 f"(use --keep_intermediate yes to keep them)")


# the output folder after `run` (organize_outputs); README.txt in the folder repeats this
OUTPUT_LAYOUT = [
    ("final_annotation/", "MitoM annotation, one GenBank file per record (<accession>.gb) - the result"),
    ("final_annotation_all.gb", "the same records in one file"),
    ("original_annotation/", "the input records unchanged (GenBank input only), same file names"),
    ("original_vs_mitom.tsv", "original vs MitoM annotation, one row per gene"),
    ("summary.tsv", "one row per record: gene counts, orientation, round, NCBI errors"),
    ("ncbi_submission/", "NCBI submission files per record (.fsa .tbl .sqn) and table2asn results"),
    ("plots/", "circular genome map per record (.png, .pdf; --plot_style circos or mitom): genes, GC content, read depth with --infq"),
    ("reports/", "one report per step (deduplication, round-2 decisions, curation, ...)"),
    ("run.log", "everything the run did, with all commands"),
    ("work/", "intermediate files (only with --keep_intermediate yes)"),
]
# --output_layout grouped (default since 2026-11-09, issue #9): the same files, results / checks / logs apart
GROUPED_LAYOUT = [
    ("summary.tsv", "one row per record: gene counts, orientation, round, NCBI errors - read this first"),
    ("annotation/", "the result: <accession>.gb per record (cox1 first, plus strand), all_records.gb, genes.fna "
                    "(protein genes and rRNAs), proteins.faa, plots/ (genome maps)"),
    ("submission/", "NCBI submission files, one folder per record (<accession>/: .fsa .tbl .sqn .val); "
                    "submission_summary.tsv, validation_messages.tsv"),
    ("qc/", "checks and decisions: original_vs_mitom.tsv and one report per step (deduplication, curation, ...)"),
    ("input/", "the input records unchanged (GenBank input only), same file names"),
    ("logs/", "run.log (everything the run did, with all commands), run_params.json, intron_fix.log, failed records"),
    ("work/", "intermediate files (only with --keep_intermediate yes)"),
]
GROUPED_DIRS = ("annotation", "submission", "qc", "input", "logs")
LOG_REPORTS = ("run_params.json", "intron_fix.log", "failed_logs", "groups")   # reports/ -> logs/ (the rest: qc/)
_GROUPED_DIR = {"ncbi_submission": "submission", "final_annotation": "annotation", "original_annotation": "input",
                "plots": "annotation/plots", "06.ncbi_submission": "submission"}


def where(args, rel: str) -> str:
    """rel, a path of the classic layout (reports/curation_report.tsv), as the user finds it after the run: the
    grouped layout moves reports to qc/ or logs/ (--output_layout; `annotate` keeps its work layout)."""
    if getattr(args, "output_layout", "grouped") != "grouped" or getattr(args, "cmd", "run") == "annotate":
        return rel
    head, _, rest = rel.partition("/")
    if head == REPORTS:
        first = rest.split("/")[0]
        return ("logs/" if rest and (first in LOG_REPORTS or first.endswith(".log")) else "qc/") + rest
    return _GROUPED_DIR.get(head, head) + (f"/{rest}" if rest or rel.endswith("/") else "")


REPORT_NOTES = {
    "dedup_report.tsv": "deduplication: one row per input record with a sequence, kept / dropped and why",
    "clade_decision.tsv": "--clade auto: phylum, clade and genetic code of every record, from its taxonomy or by blastx",
    "topology_decision.tsv": "circular or linear: each record's label, where it came from, the topology annotated, why",
    "reannotation_decision.tsv": "after re-orientation: round 2 or not, and why",
    "curation_report.tsv": "every boundary / name change: rule, old and new location, reason",
    "rescue_report.tsv": "genes added from MITOS2",
    "pcg_hmm_report.tsv": "missing protein genes searched with a profile HMM of the nearest relatives: found or not, E-value",
    "rrna_hmm_report.tsv": "rRNAs searched / added by the per-phylum rRNA HMMs",
    "intron_report.tsv": "introns found",
    "intron_ref_report.tsv": "introns placed from the intron sites of annotated relatives (models/intron_sites)",
    "orf_report.tsv": "ORFs of unknown function in unannotated sequence: written or to check, with the evidence",
    "metadata_transfer_log.tsv": "where each record's GenBank header came from",
    "oriented_by_fallback.txt": "records oriented by another gene because cox1 was not found",
    "no_cox1_samples.txt": "records without cox1 or a fallback gene: not re-oriented",
    "all_problem_samples.txt": "records with any problem during the run",
    "run_params.json": "all options of this run",
}


REPORT_FILES = tuple(REPORT_NOTES) + (
    "intron_fix.log", "annotation_failed.txt", "final_write_failed.txt", "reorient_failed.txt",
    "metadata_transfer_failed.txt", "metadata_key_mismatch.txt", "no_target_gene_samples.txt", "failed_samples.txt",
    "annotation-summary.tsv", "pipeline.log",
)                                                   # top-level files `run` writes; organize_outputs moves these to reports/
WORK_DIRS = INTERMEDIATE_DIRS + ("04.final_gb", "05.curated", "_snakemake_01.anno_round1",
                                 "_snakemake_02.anno_topology",
                                 "_snakemake_03.anno_round2_reoriented")   # folders of the run: work/ or deleted


def write_summary(out_root: Path) -> None:
    """summary.tsv from annotation-summary.tsv (internal names) + NCBI error counts, with readable columns."""
    src = out_root / "annotation-summary.tsv"
    if not src.is_file():
        return
    errs: Dict[str, str] = {}
    sub = out_root / "ncbi_submission" / "submission_summary.tsv"
    if sub.is_file():
        for line in sub.read_text(encoding="utf-8").splitlines()[1:]:
            c = line.split("\t")
            if len(c) > 1:
                errs[c[0]] = c[1]
    rows = ["record\tprotein_genes\ttRNAs\trRNAs\tcox1_first\tannotation_used\toriginal_header_restored\tncbi_errors"]
    for line in src.read_text(encoding="utf-8").splitlines()[1:]:
        c = line.split("\t") + [""] * 7
        name, pcg, trn, rrn, ori, body, meta = c[:7]
        used = ("round 2" if body.startswith("phase2") else
                "round 1, shifted" if body.startswith(("shifted_round1", "phase1b")) else
                "round 1, not re-oriented" if body.startswith("phase1_raw") else body)
        if "partial_no_gene_merged" in body:
            used += " + extra contigs"
        rows.append("\t".join([name, pcg, trn, rrn, ori, used, "yes" if "transferred" in meta else "no",
                               errs.get(name, "")]))
    (out_root / "summary.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    src.unlink()


def clear_previous_work(out_root: Path, keep: tuple = ("00.input",), protect=()) -> List[str]:
    """Remove the intermediate folders a previous run left in out_root (interrupted, or --keep_intermediate yes):
    they were picked up as this run's results (2026-10-31 review: a stale 04.final_gb/*.gbf was published, a stale
    reoriented body of another genome became a final record). Only MitoM's own folder names; `keep` = made by this run.
    Returns what was removed (for the log: the run's logger does not exist yet)."""
    removed: List[str] = []
    for name in WORK_DIRS + ("06.ncbi_submission", "work"):
        d = Path(out_root) / name
        if name not in keep and d.is_dir() and not d.is_symlink() and not any(_overlaps(d, x) for x in protect):
            shutil.rmtree(d)
            removed.append(name + "/")
    # and the previous run's results: a table / map this run does not write (FASTA input: no original_vs_mitom.tsv)
    # stayed beside the new ones (Codex R9 / N14, 2026-11-05). run.log is this run's; reports/: new_reports().
    # Only when the folder holds a MitoM output (README.txt of `run`, the log of `run` / `annotate`): a user's own
    # README.txt / plots/ in -o are not ours to delete (code review, 2026-11-06)
    names = [n.rstrip("/") for n, _ in OUTPUT_LAYOUT if n not in ("run.log", "reports/", "work/")] + \
        ["annotation-summary.tsv", "all.final.gbf", "README.txt", *GROUPED_DIRS]
    for name in names if _mitom_output(Path(out_root)) else ():
        d = Path(out_root) / name
        if name in keep or d.is_symlink() or any(_overlaps(d, x) for x in protect):
            continue                                 # (-i out/final_annotation -o out: never the input)
        if d.is_dir():
            shutil.rmtree(d)
            removed.append(name + "/")
        elif d.is_file():
            d.unlink()
            removed.append(name)
    return removed


def _mitom_output(out_root: Path) -> bool:
    """out_root holds the output of an earlier MitoM run: its README.txt, or a run.log / pipeline.log MitoM wrote."""
    for name, mark in (("README.txt", "MitoM output"), ("run.log", "[ENV] MitoM "), ("logs/run.log", "[ENV] MitoM "),
                       ("pipeline.log", "[ENV] MitoM ")):
        try:
            with (out_root / name).open(encoding="utf-8", errors="replace") as fh:
                if mark in fh.read(4096):
                    return True
        except OSError:
            continue
    return False


def _overlaps(a: Path, b) -> bool:
    try:
        a, b = Path(a).resolve(), Path(b).resolve()
    except OSError:
        return True
    return a == b or a in b.parents or b in a.parents


def new_reports(out_root: Path) -> Path:
    """An empty reports/ folder for a new run into out_root (a previous run's reports are removed)."""
    d = Path(out_root) / REPORTS
    if d.is_dir() and not d.is_symlink():
        shutil.rmtree(d)
    return reports_dir(out_root)


def organize_outputs(out_root: Path, final_gb_dir: Path, keep: bool = False, reports_named: tuple = ()) -> None:
    """Put the run's output into the layout of OUTPUT_LAYOUT (final_annotation/, reports/ ...) and write
    README.txt. Intermediate folders go to work/ (keep) or are deleted. Only the files and folders MitoM writes
    are touched (the output folder may hold the user's own files); the results of a previous run into the same
    folder are replaced. reports_named: further report files of this run whose names depend on its options
    (no_<gene>_samples.txt)."""
    def remove(p: Path) -> None:
        if p.is_dir() and not p.is_symlink():
            shutil.rmtree(p)
        elif p.exists() or p.is_symlink():
            p.unlink()

    def replace(src: Path, dst: Path) -> None:                # a previous run's copy must not swallow the new one
        remove(dst)
        shutil.move(str(src), str(dst))

    for name in ("final_annotation", "final_annotation_all.gb", "ncbi_submission", "summary.tsv"):
        remove(out_root / name)                                # results of a previous run into the same folder
    fin = out_root / "final_annotation"
    fin.mkdir()
    for f in sorted(final_gb_dir.rglob("*.gbf")):
        shutil.move(str(f), str(fin / f"{f.stem}.gb"))
    reports = reports_dir(out_root)                            # the steps wrote their reports there (new_reports())
    for f in final_gb_dir.glob("*.tsv"):                      # a report a step left with the GenBank files
        shutil.move(str(f), str(reports / f.name))
    if (out_root / "all.final.gbf").exists():
        shutil.move(str(out_root / "all.final.gbf"), str(out_root / "final_annotation_all.gb"))
    if (out_root / "06.ncbi_submission").is_dir():
        shutil.move(str(out_root / "06.ncbi_submission"), str(out_root / "ncbi_submission"))
    write_summary(out_root)
    for name in REPORT_FILES + tuple(reports_named):           # MitoM's own files only: the folder may hold the user's
        if (out_root / name).is_file():
            replace(out_root / name, reports / name)
    if (out_root / "failed_logs").is_dir():
        replace(out_root / "failed_logs", reports / "failed_logs")
    for name in WORK_DIRS:
        d = out_root / name
        if not d.is_dir():
            continue
        if keep:
            (out_root / "work").mkdir(exist_ok=True)
            replace(d, out_root / "work" / name)
        else:
            shutil.rmtree(d, ignore_errors=True)
    lines = ["MitoM output", "=" * 60, ""]
    lines += [f"{n:28s} {w}" for n, w in OUTPUT_LAYOUT if (out_root / n.rstrip("/")).exists()]
    lines += ["", "reports/", "-" * 60]
    lines += [f"  {n:28s} {w}" for n, w in REPORT_NOTES.items() if (reports / n).exists()]
    others = sorted(f.name for f in reports.iterdir() if f.name not in REPORT_NOTES)
    if others:
        lines += [f"  {n}" for n in others]
    lines += ["", "summary.tsv columns", "-" * 60,
              "  annotation_used: round 1, shifted = round-1 annotation moved with the sequence (passed the check);",
              "                   round 2 = annotated again after re-orientation; not re-oriented = no orientation gene",
              "  original_header_restored: the GenBank header of the input record (organism, references) was copied back",
              "  ncbi_errors: NCBI table2asn ERROR count (0 = ready for submission)",
              "", "original_vs_mitom.tsv result column", "-" * 60,
              "  same / ends_differ / name_differs / strand_differs / only_original / only_mitom / same_place",
              "  locations are on the original record's sequence; difference_nt = larger shift of the two ends"]
    (out_root / "README.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_sequences(gb_files: List[Path], genes_out: Path, proteins_out: Path) -> Tuple[int, int]:
    """genes.fna (protein genes and rRNAs, nucleotides) and proteins.faa (/translation of each CDS) of the final
    records; names <accession>|<gene>|<organism> (spaces: _), ready for alignment and trees (issue #9).
    Returns (genes, proteins) written."""
    from Bio import SeqIO
    genes: List[str] = []
    prots: List[str] = []
    for gb in gb_files:
        for rec in SeqIO.parse(str(gb), "genbank"):
            org = rec.annotations.get("organism", "") or next(
                (f.qualifiers.get("organism", [""])[0] for f in rec.features if f.type == "source"), "")
            org = "_".join(org.split())
            for f in rec.features:
                if f.type not in ("CDS", "rRNA"):
                    continue
                gene = f.qualifiers.get("gene", f.qualifiers.get("product", ["?"]))[0]
                name = f"{rec.id}|{'_'.join(gene.split())}|{org}"
                genes.append(f">{name}\n{f.extract(rec.seq)}\n")
                tr = f.qualifiers.get("translation", [""])[0]
                if f.type == "CDS" and tr and "pseudo" not in f.qualifiers:
                    prots.append(f">{name}\n{tr}\n")
    genes_out.write_text("".join(genes), encoding="utf-8")
    proteins_out.write_text("".join(prots), encoding="utf-8")
    return len(genes), len(prots)


def regroup_outputs(out_root: Path) -> None:
    """The classic layout organize_outputs / mitom_groups.merge_outputs wrote, arranged as GROUPED_LAYOUT
    (--output_layout grouped): annotation/ (+ genes.fna, proteins.faa), submission/<record>/, qc/, input/, logs/;
    README.txt rewritten. The open run.log is moved (same file system: the log keeps writing to it)."""
    out_root = Path(out_root)
    for name in GROUPED_DIRS:                                  # a previous grouped run's folders
        d = out_root / name
        if d.is_dir() and not d.is_symlink():
            shutil.rmtree(d)
    ann, sub, qc, inp, logs = (out_root / n for n in GROUPED_DIRS)
    for d in (ann, qc, logs):
        d.mkdir()
    if (out_root / "final_annotation").is_dir():
        shutil.move(str(out_root / "final_annotation"), str(ann / "_tmp"))
        for f in sorted((ann / "_tmp").iterdir()):
            shutil.move(str(f), str(ann / f.name))
        (ann / "_tmp").rmdir()
    if (out_root / "final_annotation_all.gb").is_file():
        shutil.move(str(out_root / "final_annotation_all.gb"), str(ann / "all_records.gb"))
    if (out_root / "plots").is_dir():
        shutil.move(str(out_root / "plots"), str(ann / "plots"))
    n_gene, n_prot = write_sequences(sorted(f for f in ann.glob("*.gb") if f.name != "all_records.gb"),
                                     ann / "genes.fna", ann / "proteins.faa")
    logging.info(f"[OUTPUT] annotation/genes.fna: {n_gene} genes, annotation/proteins.faa: {n_prot} proteins")
    if (out_root / "ncbi_submission").is_dir():
        shutil.move(str(out_root / "ncbi_submission"), str(sub))
        for f in sorted(p for p in sub.iterdir() if p.is_file() and p.suffix != ".tsv"):
            rec = f.name[:-len("".join(f.suffixes[-1:]))]       # <accession>.fsa -> <accession>/
            (sub / rec).mkdir(exist_ok=True)
            shutil.move(str(f), str(sub / rec / f.name))
    if (out_root / "original_vs_mitom.tsv").is_file():
        shutil.move(str(out_root / "original_vs_mitom.tsv"), str(qc / "original_vs_mitom.tsv"))
    rep = out_root / REPORTS
    if rep.is_dir():
        for f in sorted(rep.iterdir()):
            shutil.move(str(f), str((logs if f.name in LOG_REPORTS or f.suffix == ".log" else qc) / f.name))
        rep.rmdir()
    if (out_root / "original_annotation").is_dir():
        shutil.move(str(out_root / "original_annotation"), str(inp))
    if (out_root / "run.log").is_file():
        shutil.move(str(out_root / "run.log"), str(logs / "run.log"))
    for f in logs.glob("*.log"):                                # the report paths a step's log names (issue #6)
        if f.name != "run.log":
            f.write_text(f.read_text(encoding="utf-8", errors="replace").replace(f"{out_root / REPORTS}/", f"{qc}/"),
                         encoding="utf-8")
    old = (out_root / "README.txt").read_text(encoding="utf-8") if (out_root / "README.txt").is_file() else ""
    tail = old[old.index("\nsummary.tsv columns"):] if "\nsummary.tsv columns" in old else ""
    lines = ["MitoM output", "=" * 60, ""]
    lines += [f"{n:14s} {w}" for n, w in GROUPED_LAYOUT if (out_root / n.rstrip("/")).exists()]
    lines += ["", "qc/", "-" * 60]
    notes = {"original_vs_mitom.tsv": "original vs MitoM annotation, one row per gene (GenBank input)", **REPORT_NOTES}
    lines += [f"  {n:28s} {w}" for n, w in notes.items() if (qc / n).exists()]
    others = sorted(f.name for f in qc.iterdir() if f.name not in notes)
    lines += [f"  {n}" for n in others]
    tail = tail.replace("reports/clade_decision.tsv", "qc/clade_decision.tsv").replace("reports/groups", "logs/groups")
    (out_root / "README.txt").write_text("\n".join(lines) + "\n" + tail, encoding="utf-8")
