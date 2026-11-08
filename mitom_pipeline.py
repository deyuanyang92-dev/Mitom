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

This file is the entry point (main: checks before the run, then the workflow of the command); the parts are:
  mitom_cli.py       the command line (build_parser)
  mitom_prepare.py   environment check and [ENV] log lines, tools and --clade auto, `import` results, output lock
  mitom_run.py       the workflows run / import / annotate and the run per clade group
  mitom_commands.py  reorient, transfer_metadata, replace_organism, dedup, gb2fasta, intron
  mitom_inputs.py    input type, --clade auto, GenBank deduplication and conversion, samples and internal IDs
  mitom_stage.py     the annotation stage (engine commands, one job per sample; mitom_snakemake.py: Snakemake)
  mitom_orient.py    re-orientation, the round-2 check, rotation after the intron step
  mitom_metadata.py  original GenBank header and /organism restored
  mitom_output.py    final GenBank per sample, summaries, the output layout
  mitom_steps.py     rescue, rRNA HMM fallback, introns, curation, NCBI submission files, genome maps
  mitom_log.py       run log, numbered steps and progress on the terminal
  mito_taxa.py       clades, phyla, genetic codes; mito_common.py: shared helpers
Import names from the module that defines them (mitom_pipeline has only main()).

Biopython: Cock et al. (2009) Bioinformatics doi:10.1093/bioinformatics/btp163
Snakemake: Mölder et al. (2021) F1000Res doi:10.12688/f1000research.29032.2
"""

from __future__ import annotations

import logging
import os
import shutil
import sys
import time
from pathlib import Path

try:
    import Bio  # noqa: F401
except ImportError:                       # `mitom` checks this before it starts the script
    if __name__ == "__main__":
        print("ERROR: Biopython not installed.  pip install biopython", file=sys.stderr)
        raise SystemExit(2)
    raise

import mito_common  # noqa: E402
import mitom_cli  # noqa: E402
import mitom_commands  # noqa: E402
import mitom_context  # noqa: E402
import mitom_inputs  # noqa: E402
import mitom_log  # noqa: E402
import mitom_prepare  # noqa: E402
import mitom_run  # noqa: E402


def main() -> int:
    ctx = mitom_context.RunContext()                    # what the steps of this run share (no module globals)
    parser = mitom_cli.build_parser()
    args, extra = parser.parse_known_args()
    if args.cmd == "intron":
        args.extra = extra                       # options of mito_intron.py / mito_intron_ref.py, passed on
    elif extra:
        parser.error("unrecognized arguments: " + " ".join(extra))

    if args.cmd in ("dedup", "reorient", "gb2fasta", "transfer_metadata", "replace_organism"):
        def check_inputs() -> int:                       # missing / empty / unreadable input: one line, exit 2
            for attr, what in (("input", "-i"), ("metadata", "-m"), ("new_body", "-n")):
                if getattr(args, attr, None):
                    setattr(args, attr, mito_common.prepare_input(getattr(args, attr), what))
                    mito_common.require_records(mito_common.count_seq_records(getattr(args, attr)), "sequence records", getattr(args, attr))
            return 0
        rc = mito_common.cli_guard(check_inputs)
        if rc:
            return rc
    if args.cmd in ("annotate", "run", "import"):
        out_new = not Path(args.out_root).exists()      # an error before the run leaves no folder it created
        err = mitom_prepare.lock_output(Path(args.out_root))
        if err:
            print(err, file=sys.stderr)
            return 2
        if args.cmd == "import":                         # writes 00.imported/: only into a folder locked for us (R6)
            err = mitom_prepare.prepare_import(args)
            if err:
                print(err, file=sys.stderr)
                if out_new:
                    shutil.rmtree(args.out_root, ignore_errors=True)
                return mitom_prepare.error_exit_code(err)
        if args.cmd in ("annotate", "run") and Path(args.input).exists():   # .gz unpacked, Latin-1 text converted
            new_input, notes = mitom_inputs.normalize_input(Path(args.input), Path(args.out_root))
            args.input_given = args.input              # the user's -i, protected from clear_previous_work
            args.input = str(new_input)
            for n in notes:
                print(f"    [INPUT] {n} ({Path(args.out_root) / mitom_inputs.NORMALIZED_DIR})", file=sys.stderr)
        err = mitom_prepare.resolve_tools(args, ctx)
        if err:
            print(err, file=sys.stderr)
            if out_new:
                shutil.rmtree(args.out_root, ignore_errors=True)
            else:
                shutil.rmtree(Path(args.out_root) / mitom_inputs.NORMALIZED_DIR, ignore_errors=True)
            return mitom_prepare.error_exit_code(err)

    if getattr(args, "suffix_fq", None):
        mitom_inputs.parse_suffix_fq(args.suffix_fq)

    if (getattr(args, "auto_accept_metadata", "no") == "yes"
            and getattr(args, "smart_metadata", "no") != "yes"):
        print("ERROR: --auto_accept_metadata yes requires --smart_metadata yes", file=sys.stderr)
        return 2

    dispatch = {
        "annotate":          lambda a: mitom_run.workflow_annotate(a, ctx),
        "reorient":          mitom_commands.workflow_reorient,
        "run":               lambda a: mitom_run.workflow_run(a, ctx),
        "import":            lambda a: mitom_run.workflow_run(a, ctx),
        "transfer_metadata": mitom_commands.workflow_transfer_metadata,
        "replace_organism":  mitom_commands.workflow_replace_organism,
        "dedup":             mitom_commands.workflow_dedup,
        "gb2fasta":          mitom_commands.workflow_gb2fasta,
        "intron":            mitom_commands.workflow_intron,
    }
    fn = (lambda a: mitom_run.workflow_groups(a, ctx)) if getattr(args, "clade_groups", 0) else dispatch.get(args.cmd)
    if fn is None:
        print(f"Unknown command: {args.cmd}", file=sys.stderr)
        return 2
    if getattr(args, "debug", False) or mito_common.debug_mode():     # --debug (or `mitom --debug <command>`: MITOM_DEBUG=1)
        args.debug = True
        os.environ["MITOM_DEBUG"] = "1"                   # the programs MitoM starts log what they do too
        if hasattr(args, "keep_intermediate"):
            args.keep_intermediate = "yes"
    t0 = time.time()
    try:
        rc = fn(args)
    except mito_common.InputError as exc:                # a problem of the input found during the run: one line, exit 2
        logging.error(f"[ERROR] {exc}")
        logging.info(f"[DONE] exit=2 {mitom_log.duration(time.time() - t0)}")
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:                 # into the run's log with the traceback (mitom debug reads it), then exit 1
        if mitom_log.Steps.current is not None:
            mitom_log.Steps.current.failed(exc)
        logging.exception(f"[ERROR] {type(exc).__name__}: {exc}")
        logging.info(f"[DONE] exit=1 {mitom_log.duration(time.time() - t0)}")
        out = getattr(args, "out_root", None)
        print(f"MitoM stopped on an error: {type(exc).__name__}: {exc}" + (f"\n  diagnosis: mitom debug {out}" if out else ""),
              file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        logging.error(f"[ERROR] interrupted (Ctrl-C) after {mitom_log.duration(time.time() - t0)}")
        raise
    logging.info(f"[DONE] exit={rc} {mitom_log.duration(time.time() - t0)}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
