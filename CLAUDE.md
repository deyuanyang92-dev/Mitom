# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

- Install for development: `pip install -e .[test]` (external programs: `install.sh` / `env/mitotools.yml`; `mitom setup` fetches refseq89m, `mitom doctor` lists found programs).
- Tests: `python3 -m pytest -q tests`; one file/test: `pytest -q tests/test_mito_trna.py -k name`.
- Whole-workflow regression (needs BLAST, exonerate, Infernal, HMMER, table2asn, refseq89m): `MITOM_GOLDEN=1 pytest tests/test_golden.py` or `python3 scripts/golden.py [--cases a,b]`. References in `tests/data/golden/` were made on the machine in `MACHINE.txt`; after an intended output change rebuild with `--update` and commit with the change.
- After changing any CLI option: `python3 scripts/make_options_doc.py` (regenerates `docs/OPTIONS.md`; `tests/test_docs.py` fails otherwise).
- Run: `mitom run -i genomes.gb -o out` (FASTA: add `--clade X --genetic_code N`, or `--clade auto`); `mitom debug OUT` explains a failed run; `--debug` on any command keeps intermediates.

## Architecture

- **Flat layout, no package.** Every module is a top-level `*.py` (`setup.py` globs them) and each script is also runnable as `python3 <script>.py`. `mitom_main.py` is the `mitom` entry point: a `COMMANDS` table mapping subcommand -> (script, fixed args) that launches that script with a Python having Biopython/scikit-learn. A new command = a row there + options doc.
- **Naming:** `mito_*.py` = biology/algorithm modules (`mito_core` own annotation core: blastx vs refseq89m + exonerate, `mito_trna`/`trna_ml` tRNA callers, `mito_intron*` MFannot-method introns + known-site transfer, `mito_curate*` NCBI boundary rules P1-P11, `mito_rna*` RNA-seq evidence, `mito_nuc*` nuclear rDNA/histones, `mito_engines` MitoZ/MITOS2 as engines or imported results). `mitom_*.py` = workflow plumbing (`mitom_pipeline.py` is only the entry point; `mitom_run.py` holds `RUN_STAGES` / `ANNOTATE_STAGES` that hand a `RunState` along; `mitom_steps.py` the post-annotation `run_<step>(args, gb_dir, out_root)` steps; `mitom_groups` per-clade-group runs for `--clade auto`; `mitom_output` output layout; `mitom_snakemake`).
- **Shared state** goes through `RunContext` (`mitom_context.py`), never module globals.
- **Structure rules enforced by `tests/test_structure.py`:** no `global` statements, no `sys.path` edits, GenBank written only via `mito_common.write_gb`, external programs started only via `mito_common.run_tool` (exceptions listed in the test), `mitom_pipeline.py` has only `main`.
- **Data:** `models/` (HMMs, intron-site library, tRNA/start models, reference proteins) ships with the repo and carries mixed licences (see `pyproject.toml`). `benchmark/` holds evaluation scripts (`benchmark/frozen/` = pre-registered frozen evaluation, not yet to be run). `knowledge/` is the cited knowledge base (`INDEX.md`; check with `python3 benchmark/check_knowledge.py knowledge/`).

## Notes for agents

- **Tests: check `knowledge/TEST_REGISTRY.md` before running any benchmark / real-data test** (two machines run tests:
  `deyuan` = /media/deyuan/18T_B/yzz/mitom_work, `cody` = /media/cody/18T_R1). Do not repeat a test that is listed
  unless the code it measured has changed; say which row you rely on. After any test add one row there (date, machine,
  question, data, tools/versions, key numbers, conclusion/commit, where recorded) and keep details in
  `benchmark/results/<name_date>/README.md`. Commit both so the other machine sees them.
- Run errors: cause and fix go to the troubleshooting tables (docs/USAGE_zh.md section 9, docs/MITOS2.md section 6).
- Docs for users: Chinese + English, short and plain. A measured number (accuracy, time, memory, counts from a benchmark) is written once, in a table with bilingual headers. The Chinese paragraph and the English paragraph point at that table and do not repeat the number. Option defaults and rule thresholds stay next to the option or the rule, in one place. Full options: `python3 scripts/make_options_doc.py` after CLI changes
  (tests/test_docs.py checks it).
