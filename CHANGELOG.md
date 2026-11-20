# Changelog

Notable changes of MitoM. Versions follow semantic versioning (major.minor.patch); `mitom --version` prints the
installed one.

## Unreleased (0.0.2)

### Added
- `mito_intron --short_exons` (experimental, off by default): homolog-vote search for short divergent exons Exonerate does
  not align (issue #15). No gain on the 134 NCBI records; benchmark/results/short_exons_2026-11-20/README.md.

### Fixed
- Intron transfer: a strong, continuous model is no longer refused because of a weak conflicting core CDS inside its
  exon (MW592986.1 atp6, issue #15); the weak call is replaced. 1 of 134 NCBI records changed.
  benchmark/results/weak_conflict_2026-11-20/README.md.

## 0.0.1 - 2026-11-20 (beta, first tagged test release)

MitoM is still being tested: check introns, tRNAs and gene boundaries by hand before publishing or submitting to NCBI.

### Fixed (review of 2026-11-20; benchmark/results/review_check_2026-11-20/README.md)
- Exonerate second batch: the better 5'-partial placement was chosen by comparing a window end with a score.
- `complete_claim`: every sample of a multi-record GenBank file took the first record's DEFINITION (4 of 170 annelid
  RefSeq records would have been turned circular wrongly).
- Curation: a stop-looking codon read by an internal `/transl_except` is no internal stop (the CDS was written
  `/pseudo`); an imported `/pseudo` CDS is left as it is and gets no `/translation`; O1 keeps internal exceptions.
- rRNA fill (R1) no longer stretches an rRNA to the end of a linear sequence; `genes.fna` leaves out `/pseudo` CDS.
- `mito_nuc` reads `.gbff`; `mito_check` checks the last coding codon of an abbreviated-stop CDS; the intron-tool
  check of `mitom run` works when exonerate is missing.
- Known-site intron transfer keeps the de novo model when neither the merged nor the transferred model translates.
- `mito_pcg_hmm`: `/translation` written, genes across the origin compared on the circle, `--k` help corrected.
- `.mitom.lock` is created atomically (48 runs started together: exactly one gets the folder).

### Added
- `mitom run --coverage auto|yes|no` (issue #12): with `--infq`, qc/coverage.tsv (mean / median / min / max depth,
  % bases >= 1x / 10x / 100x, mapped reads, topology, GC %) and qc/depth.tsv.gz; summary.tsv gets mean_depth and
  breadth_1x_pct. `mitom coverage` writes the same columns and `--depth_gz`.
- `mitom setup` downloads refseq89m from Zenodo or its copy on GitHub (md5 checked), with curl as a fallback
  (issue #11).

## Earlier changes (before the first tag)

### Fixed (frameshift joins)
- Curation P9: before joining a protein gene over a "frameshift", an uninterrupted ORF in another frame covering
  >= 85 % of the gene is taken instead - the core had put the gene in the wrong frame. Travisia forbesii ATP8
  (OZ479805): 49-aa join -> 53-aa single ORF as in T. sanrikuensis and another Travisia (RNA-seq and HiFi keep the
  'skipped' base). Real frameshifts unchanged (turtle nad3 joined at position 174); 1 of 873 protein genes of 67
  annelid genome-project mitogenomes changed. benchmark/results/annelid_rna_2026-11-10/README.md.

### Added (RNA-seq)
- `mitom rna` (`mito_rna.py`): RNA-seq mapped to the mitogenome only; per gene / ORF / non-coding region the depth
  relative to the 13 genes, the share transcribed, an expressed-like-the-genes call, the library strand (reverse
  stranded recognised) and poly(A) 3' ends; few poly(A) reads or too little depth are reported, not silent. On the
  Polydora RNA-seq of Selifanova et al. 2023 it reproduces all three ORFan conclusions.
  benchmark/results/rnaseq_2026-11-10/.

### Added (parts of one genome)
- `mito_multipart.py`, in step 9: records of one assembly (COMMENT) or BioSample (DBLINK) whose protein genes are
  complete only together are reported as parts of one genome (multipartite mtDNA or pieces) in
  multipart_report.tsv with a WARNING; copies of one sequence are counted once. On 16 multi-record genome-project
  assemblies: 6 groups found (earwig, millipede and nematode two-circle mtDNA; three butterflies in two pieces).
  benchmark/results/ebp_mito_2026-11-10/.
- Benchmark on 221 complete genome-project (DToL / ERGA / VGP) mitogenomes; the unexpected missing genes are a
  three-chromosome sponge, a two-circle nematode and a myxozoan gene loss (MitoM right).

### Fixed (benchmark)
- `benchmark/mito_sample.write_bare` kept no DBLINK (record slicing drops dbxrefs).

### Added (missing protein genes)
- `mito_pcg_hmm.py`, part of step 6 (`--pcg_hmm yes`, default): protein genes still missing after the core are
  searched with profile HMMs (MAFFT + HMMER) of the gene in the record's nearest relatives (models/pcg_refs.faa.gz,
  140,824 proteins of RefSeq + NCBI annelid mitogenomes; `scripts/build_pcg_refs.py`), E <= 0.01 x the prior odds of
  the gene among the relatives. Leave-genus-out: 2,412 of 2,448 genes found at the right place, atp8 187/189; 335
  unannotated NCBI annelids: atp8 found in 29 of 38 records missing it (tube worms 19/26).
  benchmark/results/pcg_hmm_2026-11-10/.

### Fixed (ORFs)
- A tandem repeat is masked only when the array is >= 30 nt (3 copies of CTT, Leu-Leu-Leu, were masked in real genes).
- An ORF whose start codon lies beyond the origin of a circular record got a start past the record end and crashed the
  ORF step (OZ210882.1, 47,540 nt).

### Changed (structure review of fb526c3)
- ORFs and non-coding regions are their own step (`stage_orf_ncr`, step 8: `run` has 12 steps, `annotate` 8);
  `stage_curate` does the curation only.
- The ORF step runs `mito_orf.py` as a program, like the curation and the introns (crash isolation; when it fails the
  records stay as they were); `tests/test_step_argv.py` checks its arguments. METHODS states which steps run as
  programs and which in process (NCR, write-back, maps: light, pure Python, exceptions logged and skipped).
- Golden case `annotate_rotated`: `annotate` on NC_000931 with its start moved 5,000 nt locks the rotation back to
  cox1 (the existing `annotate_gb` input already starts at cox1).

### Changed (ORFs, after a literature / GitHub survey)
- A frame that runs into a tRNA without a stop ends before the tRNA (tRNA punctuation; open 3' end when the tail is
  not T / TA), as Selifanova et al. 2023 define the Polydora ORFans: all three (ORFan-544 / -1060 / -427) are now
  written in P. cf. ciliata OQ078742 and its relatives. Pfam families (models/orf_families.hmm.gz, hmmsearch
  --cut_ga) name homing endonucleases, intron reverse transcriptases, polB and MutS. Conserved = both tests at 0.05
  (instead of Fisher <= 1e-4). RNAcode was tested and adds almost nothing; not used. Known ORFs written 228 -> 233.

### Changed (ORFs of unknown function, after a recall test on 342 NCBI-annotated non-standard CDS)
- Introns are searched (homing endonucleases, intron mutS); substitutions are summed over up to 20 homologs (the
  larger of one per extra base kind per position and half the MST of their differences); the frame test uses the
  full codon pattern (position 3 most, 2 least; chi-square), which also picks the right strand of overlapping ORFs;
  the two tests are Fisher-combined (<= 1e-4, each <= 0.05). Known ORFs written 175 -> 228 of 342; intron-free
  RefSeq records 5 writes (ORF finding by length >= 75 codons: 357). ORF library rebuilt with these rules.

### Added (ORFs of unknown function)
- `mito_orf.py`, step of `run` / `annotate` / `import` after the curation (`--orf yes`, default; `all`, `report`,
  `no`; `--orf_ref` relatives): an ORF in unannotated sequence is written as CDS `orf<codons>` "hypothetical protein"
  when a homolog stayed open through enough substitutions to show selection and the differences fall mostly at codon
  position 3 (tandem repeats masked, shadow frames dropped), or is a known protein family (atp9 named). ORFs only
  longer than chance, pieces of the 13 genes and frames running into the next gene are WARNINGs, rows in
  `orf_report.tsv` and text in the NCR /note. ORF library `models/orf_library.fna.gz` (16,073 RefSeq animal
  mitogenomes, `scripts/build_orf_library.py`). benchmark/results/orf_2026-11-09/.
- Annelid intron evaluation, including the 335 NCBI annelid records without any annotation:
  benchmark/results/intron_annelid_2026-11-09/.

### Changed
- The intron_report.tsv row for a CDS that still has an internal stop after the intron step says the curation may
  still correct it (two Chaetopteridae cox2, corrected by rule P7).

### Changed (introns)
- An intron of no known type (no Rfam group I / II hit) is written when the gene's homologs run on across its
  junction (no residue missing or added; a discontinuous junction is searched +-30 nt for a continuous placement) and
  splicing is needed (reading on in place covers the homolog >= 20 BLOSUM62 points worse), instead of a group I splice
  score >= 14 and >= 80 % homolog coverage (`mitom run --intron_rule splice`, `mito_intron.py --untyped_rule splice`).
  Bryozoan introns found 12 -> 31 of 50; gold standard and 362 intron-free records unchanged
  (benchmark/results/intron_algo_2026-11-09/).
- A gene whose annotated CDS has an internal stop codon is an intron candidate; a CDS that still has one after the
  intron step is listed in intron_report.tsv. intron_report.tsv gains junction_span, junction_missing_aa,
  junction_extra_aa, splice_score, splicing_gain.

### Fixed
- `benchmark/mito_bench_run.py` and `benchmark/eval_introns.py` read the grouped output layout (annotation/,
  qc/intron_report.tsv) as well as the classic one; they found no output after the layout change of #9.

### Added
- Intron test on the NCBI animal mitogenomes with introns outside the gold standard (134 records):
  benchmark/results/intron_ncbi_2026-11-09/.

### Fixed (GitHub issues, 2026-11-09; checked on NCBI MZ052209.1, OZ287881.1, OZ346613.1)
- #6: `logs/intron_fix.log` (and run.log) named `04.final_gb/intron_report.tsv`, a folder removed after the run; the
  paths are now where the reports end up (also in merged group runs).
- #8: GenBank records labelled linear whose DEFINITION says they are whole genomes (`complete genome`, or the DToL /
  ATLASea / MitoHiFi `genome assembly, organelle: mitochondrion`), with the complete gene set, no gene cut at an end
  and not of a linear-mtDNA lineage, are now circular and rotated to cox1 (OZ287881: the 6.1 kb non-coding region
  split over both ends is one NCR1 of 6,121 bp). A record turned to the strand of its orientation gene keeps
  `source 1..L` (was `complement(1..L)`); the log says `[LINEAR_REVCOMP]` when a linear record is reverse-complemented.
- #7 (3 trnP in MZ052209, reported on 8a0debe): not reproduced on the current code (one trnP at 7832..7894, 22 tRNAs);
  MZ052209 is now a golden case (run_issues).

### Changed
- `start_run` takes the stage list: the step count and the round-2 line on the terminal follow it (was a default of 11
  and `hasattr(args, ...)`); stage docstrings no longer carry run step numbers. `mitom annotate --help` says what it
  does (7 steps; records rotated to `--which_gene_first` without the round-2 check). Docs: the METHODS flowchart
  (04.final_gb / annotation-summary era) is redrawn from the 11 run stages; benchmark numbers that were in two places
  (table2asn ERROR count, intron gold standard, nuclear-marker accuracy and time, P3 20 nt limit) are in one.
- The NCBI step submits when the curation says it wrote its folder (`StepResult.curated`), not by comparing folder
  paths. Log tags `[PHASE1b]` / `[PHASE2]` are `[REORIENT]` / `[ROUND2]`. minimap2 speed (2.7 s vs bwa 21 s) and
  the time per record are given once (USAGE 1.5, README time table); other places point there.
- #9: output of `run` / `import` (`--output_layout grouped`, default): summary.tsv, README.txt and annotation/
  (<record>.gb, all_records.gb, genes.fna, proteins.faa named `<accession>|<gene>|<organism>`, plots/),
  submission/<record>/ (.fsa .tbl .sqn .val; submission_summary.tsv), qc/ (original_vs_mitom.tsv and the step
  reports), input/ (the input records), logs/ (run.log, run_params.json, intron_fix.log, failed_logs/, groups/).
  `--output_layout classic` keeps the layout before (final_annotation/ ncbi_submission/ plots/ reports/ run.log).
  Messages name the new paths; `mitom debug` reads both layouts. Wrappers that read reports/intron_report.tsv: read
  qc/intron_report.tsv, or run with `--output_layout classic`.

### Added
- `--mitos_args` (run, annotate): further runmitos options for `--engine mitos2`, written `--mitos_args="--best --finovl 30"`;
  options MitoM sets itself (input, output, -c, -r / --refdir, --linear, --json) are refused before the run; without it
  the MITOS2 command is unchanged. Real run (NC_000931): the option reaches runmitos (run.log, core.stderr.log).
- docs/MITOS2.md (Chinese / English): MitoM as batch post-processing of MITOS2 - GenBank output, orientation, curation,
  NCBI files; `--engine mitos2` and its options, `mitom import`, every output file, Geneious -> `mitom curate`,
  accuracy, from one real run (NC_000931).

- benchmark/frozen/: frozen evaluation for the paper numbers - plan (zh/en), leak audit and test-set selection
  against everything in the repository (frozen_set.py), freeze manifest (freeze.sh), one guarded run (run_frozen.sh),
  family-bootstrap intervals and paired differences (frozen_stats.py); benchmark tool `mitom_blind` (sequence only).
  Measured: 169 of 170 annelid RefSeq records (2025-03 pool) were already seen, so the test set draws on GenBank too.

### Changed (structure; outputs unchanged - scripts/golden.py)
- `mitom annotate` runs the stage functions of `run` (`ANNOTATE_STAGES`: input, round 1, write, add missing, curate,
  NCBI, maps) instead of its own copy of the steps after annotation. Results identical (golden case annotate_gb);
  new: the terminal steps, `original_annotation/`, `reports/clade_decision.tsv` and `metadata_transfer_log.tsv`,
  `annotation-summary.tsv` column body_source; the log is run.log (was pipeline.log). Golden: 7 cases (annotate_gb added), LOCUS dates normalised.
- README: the benchmarks have their own section (Development benchmarks) after the statement; data sizes on one
  data line per benchmark. MITOS2.md: the 20 nt stop rule in a table.
- mito_curate.py (1433 lines) split by rule family, code moved verbatim: `mito_curate_base.py` (limits, notes, span
  helpers), `mito_curate_pcg.py` (P1-P11, D1, /translation), `mito_curate_origin.py` (O1 / O2), `mito_curate_rna.py`
  (R1, T1, C1, N1); `Curator` (479 lines with the CLI) combines them. mito_intron_ref.py (1655 -> 712 lines):
  `mito_intron_junction.py` (query genome, junction calling) and `mito_intron_models.py` (gene models).
- `setup.py` installs every top-level *.py (no hand-kept py-modules list); the packaging test builds a wheel and
  compares. README step list = the 11 terminal steps (details stay in USAGE_zh.md section 4).
- Docs: a measured number is written once, in a table (bilingual headers). Chinese and English prose point at the table and do not repeat it. README benchmark blocks, MITOS2.md timings, METHODS exceptions, and the USAGE accuracy summary follow this. `mito_intron.py` and `mito_nuc.py` stay one file each.
- `tests/test_step_argv.py`: the arguments the pipeline passes to mito_curate.py / mito_intron.py / mito_intron_ref.py
  are parsed with those scripts' own parsers (they stay separate processes: a crash affects one step only).
  Docs: METHODS.md opens with the code layout (`mitom` -> mitom_main -> mitom_run.RUN_STAGES) instead of the stale
  "mitom_pipeline.py is the main program" and a second copy of the step list; README / MITOS2.md accuracy numbers
  are dated as development benchmarks (2026-10), not the frozen evaluation.
- Refactoring guard: `scripts/golden.py` runs 6 whole workflows on tests/data (GenBank, intron, two phyla, bare FASTA,
  MitoZ import, curate) and compares every text output with tests/data/golden (`MITOM_GOLDEN=1 pytest tests/test_golden.py`).
- `pyproject.toml`: MitoM is a Python project (`pip install -e .` gives the `mitom` command; version from
  mitom_main.VERSION). The launcher's code moved from the `mitom` file to `mitom_main.py` (the file stays, 10 lines).
  The modules no longer edit `sys.path` (23 files): a script imports its neighbours from its own folder, as Python does;
  `python -I script.py` / PYTHONSAFEPATH are not supported for direct script runs (`mitom` clears PYTHONSAFEPATH).
- mitom_pipeline.py (1678 lines: workflows, command line, preflight, ~200 re-exported names) split, code unchanged:
  `mitom_cli.py` (build_parser), `mitom_prepare.py` (environment check, tools, --clade auto, import, output lock),
  `mitom_run.py` (run / import / annotate, clade groups), `mitom_commands.py` (reorient, transfer_metadata,
  replace_organism, dedup, gb2fasta, intron); mitom_pipeline.py keeps main() (183 lines). No re-exports: an old
  `mitom_pipeline.<name>` no longer resolves (no script used it): import from the defining module.
- No module globals that carry one run's state: `mitom_context.RunContext` (made by main(), passed to the run
  workflows) holds the --clade auto plan (was mitom_groups.LAST_PLAN), the GenBank dedup counts (mitom_inputs.
  LAST_DEDUP_STATS) and the re-orientation notes (`mitom_orient.OrientNotes`: fallback start genes, linear records
  reported; were ORIENTATION_FALLBACK, _LINEAR_SKIP_SEEN); the annotation progress bar is passed to annotate_one
  (mitom_stage._STAGE_PROGRESS). In the leaf scripts the values set from the command line are passed down instead of
  rewriting module globals: mito_trna circular length (CIRC_L -> circ_len), mito_core --exclude_ref and
  --exonerate_model (tools["exclude_refs"], tools["exonerate_model"]), mito_intron --timeout. mito_intron's
  failure-tolerant program call is named `try_program` (was a second `run_tool`, the name of mito_common's).
- `workflow_run` (350 lines) = 11 stage functions in `mitom_run.RUN_STAGES` (input, round 1, reorient, round 2, final
  GenBank, add missing, curate, NCBI, compare, maps, output; one numbered step each) handing on a `RunState`;
  workflow_run itself is 10 lines. Code of each stage unchanged; the final-body choice is `final_body()` /
  `with_partial_contigs()`. run.log and the terminal lines are the same as before (checked on NC_000931).
- In-process API of the leaf commands: `mito_curate`, `mito_intron`, `mito_intron_ref`, `mito_nuc` have
  `build_parser()`, `options(...)` (the command-line defaults, required options checked; `mito_common.options`) and
  `run(options)`; `main()` only parses. E.g. `mito_curate.run(mito_curate.options(input="gb/", out="cur/"))`
  (tests/test_api.py: same output as the command). The workflow still starts them as programs (timeouts, logs,
  Snakemake). mito_nuc's program helper `run(cmd)` is `program_output(cmd)`.
- External programs go through `mito_common.run_tool` (logged, `mitom debug`) in mito_nuc (blastn, exonerate,
  table2asn) and mitom_groups too. tests/test_structure.py keeps the rules: no `global`, no sys.path edits, GenBank
  written only by write_gb, subprocess only in run_tool and the listed probes / pipes, mitom_pipeline.py = main().
- The largest files split along their sections, code unchanged, names still importable from the old module:
  mito_nuc.py (2320 lines) -> `mito_nuc_lib.py` (settings, programs, reference records and library),
  `mito_nuc_rdna.py`, `mito_nuc_histone.py`, mito_nuc.py (input, output, NCBI files, annotate, command line; 1046);
  mito_intron.py (1668) -> `mito_intron_rules.py` (gene model, MFannot boundary rules; GPL-3.0 like mito_intron.py)
  + mito_intron.py (1242); mito_intron_ref.py (2057) -> `mito_intron_sites.py` (intron-site library) +
  mito_intron_ref.py (1655). mito_curate.py stays one file (one Curator class whose rules share its state).
  Checked: golden outputs identical; `mitom nuc annotate` on 3 rDNA contigs identical before / after.
  pyproject license = CC0-1.0 AND GPL-3.0-only AND MIT (was the CC0 file only).

### Fixed
- MITOS2 installed as bioconda mitos 2.1.9 (only `runmitos.py`, no `runmitos`) was not found (`mitom doctor`:
  runmitos not found): both names are accepted, by every lookup (mito_env.find: doctor, debug, mito_trna).
- `--mitos_args` let MitoM's own options through in the forms argparse accepts: attached (`-c2`, `-o/tmp`) and
  abbreviated (`--cod 2`, `--outd x`); they overrode the genetic code / output. All forms are refused now (Codex
  release-readiness review 2026-11-08; tests/test_review_2026_11_08.py).
- `mitom gb2tbl` read no `.gb.gz` (a folder of them: "no GenBank / EMBL file"); now unpacked, output named without
  `.gb`. The other commands already read them (checked one by one).

## 0.0.1 - 2026-11-07

First numbered test release. MitoM is still being tested, so the version starts again at 0.0.x; the
internal numbers 0.9.0 and 1.0.0 below were never released (logs of those runs still say them).

### Changed
- `--clade auto` (default) per record, no clade needed: a GenBank record by its taxonomy and CDS /transl_table (no CDS:
  the code of its lineage in the reference genomes, e.g. tunicates 13); FASTA and records without a taxonomy by
  blastx to 8 conserved proteins of refseq89m (mito_classify.py, models/ref_taxa.tsv.gz built from local RefSeq by
  scripts/build_ref_taxa.py): phylum of the closest genomes, their genetic code checked on the genome's own AGA AGG
  ATA AAA TGA TAA codons (codes reading AGR as stop left out when AGR occurs inside the genes; codons seen no more
  often than TAG are noise); no close relative (phylum support < 0.7): clade Animal, code from the codons. An input
  of several phyla / codes is no longer an error: `mitom run` annotates each group with its own clade and code (in
  parallel, --jobs shared) and merges the results into one output folder (mitom_groups.py; reports/clade_decision.tsv,
  reports/groups/). FASTA input no longer falls back to Annelida. `annotate` / `import` still need one phylum.
- Topology without a label (FASTA without topology=): decided from the sequence for every record of >= 10 kb, also
  when the file holds several (they were all linear: 81 of 83 complete genomes started at a gene boundary stayed
  linear). In order: a gene across the ends -> circular; a gene cut at an end -> linear; a lineage with linear mtDNA
  -> linear; else the gene set: >= 12 of the 13 protein genes and both rRNAs -> circular, otherwise linear (partial
  genome, contig). Before annotation: an assembler's terminal overlap (the end repeats the start, >= 50 nt, not a
  tandem repeat, as MITObim circules.py / MitoHiFi / NOVOPlasty) -> circular, the repeat removed; header words
  ('(circular)' of GetOrganelle, NCBI 'complete genome' / 'partial genome' / 'nearly complete') are a label.
  Benchmark (wgs_skim_test/mitom_topology, 364 unlabelled sequences of 91 NCBI genomes, 7 phyla): 64.4 -> 84.3 %
  correct; complete genomes 2 -> 81 of 83 (NCBI start), 48 -> 81 (random start); fragments of 3-12 kb 85 -> 89 of 91.
  Limits: a contig missing only non-coding sequence or tRNAs reads as a complete circle (56 of 91 contigs missing
  0.5-2 kb, 88 before); NCBI 'partial' genomes holding every gene too; myxozoans (8-11 protein genes) linear.

### Fixed
- `mitom import` reads gzip-compressed GenBank (`.gb.gz`, `.gbf.gz`): a single `.gb.gz` was `no GenBank`, and
  in a folder with `.gb` and `.gb.gz` files the compressed ones were left out without a note.
- Install and self-check: the `mitom` command picks a Python with Biopython, numpy and scikit-learn (the
  mitotools environment first; an old symlinked command ran a Python without scikit-learn and the MitoM-ML tRNAs
  failed); `mitom doctor` shows the version, folder, git commit and Python in use, checks scikit-learn / numpy, and
  flags a `mitom` on PATH that runs another copy (with the fix); README: how to update (`git pull && bash install.sh`).
- Cursor review of the fixes (Mitom_review/FIX_REVIEW_2026-11-05.md): explicit `--clade Chordata` with
  records of different codes (vertebrates 2 + tunicates 13) is an error, not code 2 for all; `--clade auto` reads
  .gb.gz; a.fa.gz.gz unpacks to a.fa; a failed group's tables stay out of the merged reports; a tRNA call between
  two separate same-aa tRNAs no longer merges them; a site at the origin of a circular record is join(n..n,1..1).
  Bugbot on PR #3: tRNA conflicts are checked on every member of a group, a call without anticodon no longer
  chains trnL1 and trnL2; `mitom curate` keeps the inferred /transl_table when a record fails.
  Code review: previous results are cleared only from a MitoM output folder (README.txt / run log), the given -i
  is protected and every removal is logged; x.embl.gz becomes GenBank; `import` keeps GenBank files beside a MitoZ
  *.result/; curate never writes the fallback code 5; dropped features are [RECORD] events; a failed subtype-model
  download is retried after a day. Codex review: `import` protects its -i; tRNA calls past the end of a linear
  record are left out; a table2asn timeout makes the submission step fail like a crash.
- Model reviews round 2 (Grok 4.7, GLM 5.3, Kimi K3, 11-03 static review; each item checked and reproduced,
  Mitom_review/FIX_STATUS_2026-11-05.md). Annotation: input tRNAs across the origin are no longer the whole record;
  tRNA calls between two others join their group; intron-step neighbours read by part (a join across the origin or
  a gene inside its intron was skipped); atp8 / nad4L fallback ORFs across the origin; ARWEN anticodons kept;
  `mitom curate` without --genetic_code uses the taxonomy's code and writes /transl_table (was 5); CDS with
  /codon_start 2 / 3 left out of the boundary rules; a CDS named `atp` is no longer renamed ATP6 and then deleted as a
  second ATP6; "cytochrome c oxidase subunit 1" products recognised; explicit `--clade Chordata` and Urochordata
  give tunicates 13, pterobranchs 24 / 33; atp9 kept from MITOS2. Output: every GenBank writer keeps a value-less
  /pseudo (the final records of `mitom run` had /pseudo=""), atomic writes, n^1 sites, plastids. Logs / debug: the
  logs of topology-phase and round-2 failures kept, failed programs' messages, step failures as [RECORD] lines.
  Fetch: "complete mitochondrial genome" titles. Python 3.8.0-3.11.3 tar extraction. Tests: an environment leak
  made two tRNA tests fail only in the full suite. Not changed (checked, not bugs): GLM 2.1 (its proposed change
  would drop the gene), Grok 1.1 / 2.4, Kimi 2.6 / 3.5. 49 NCBI records old vs new: identical coordinates.
- Model reviews 2026-11-04 (Mitom_review/VERIFICATION_REPORT.md, NEW_BUGS_2026-11-04.md; checked and fixed
  2026-11-05, tests/test_review_2026_11_05.py):
  - NCBI check: a table2asn that crashed or could not run (no .val) was counted as 0 ERROR (run, curate, gb2tbl,
    nuc); now MITOM.Table2asnFailed, the submission step fails. Messages of severity `Critical` (table2asn 1.28) were
    not read at all; REJECT / Critical / FATAL now count in the ERROR column. table2asn's exit code is no sign of a
    crash (it exits 2 on records with validation errors). A gb2tbl timeout was counted to SeqID '' (no row).
  - .tbl: value-less qualifiers (/pseudo, /trans_splicing, /ribosomal_slippage) as the key alone (they were written
    `pseudo<TAB>None` from records in memory, e.g. `mitom curate`'s check after the rules); '(MitoM)' left out of
    the NCR /note.
  - `--topology linear` no longer cuts an end that repeats the start; a GenBank record without topology on its LOCUS
    line is linear below 10 kb (as FASTA). `annotate` with nothing to annotate exits 2, and `run` / `annotate` exit 1
    when no final GenBank was written; transfer_metadata / replace_organism exit 1 when a file failed.
  - `import --clade auto` of a GenBank without CDS /transl_table: the lineage's code as `run` (tunicates 13, was 2).
  - Intron step (mito_intron.py): exonerate gets the 64-letter table for codes it lacks (24, 33), as the core does.
  - install.sh: the macOS Miniforge download (MacOSX, not Darwin); mitotools always built (`--envs mitos`).
  - Output folder: `import` writes only after the folder is locked; a lock of another user's run is not taken over;
    reusing -o removes the previous run's results (a table / map this run does not write stayed), never the input;
    `same.fa` + `same.fa.gz` is an error (one replaced the other), gzip of gzip unpacked.
  - Groups: a failed group's annotation is not merged (its reports are); tables with different headers merged by
    column name. Rescue / round-2 check know genes named cox1 / rrnS (COX1 was added twice); the MITOS2 runs of the
    rescue share the CPUs. A feature that cannot be moved to a new origin is left out (it kept its old coordinates).
    Import: GenBank files beside MitoZ folders are imported, renamed duplicates never take an existing name.
    Metadata: the accession from the full id, the LOCUS length of smart matches. GenBank results written atomically.

- tRNA Leu / Ser numbers follow MITOS2 2.1.10 for every taxon, vertebrates included (L1 = CUN, L2 = UUR, S1 = AGN,
  S2 = UCN; NCBI's vertebrate records number them the other way): unchanged behaviour, now documented and tested.
- Robustness round 2 (2026-10-31: every command fuzzed with 24 good / malformed inputs, four module reviews):
  - Reverse complement reversed the parts of every join (exons, trans-spliced and origin-crossing features came out
    scrambled after a record with cox1 on the minus strand was turned; 6 of 208 features in the test records), and
    rotation wrote a minus-strand feature across the new origin in the wrong order (50 of 272 random rotations).
    Now every feature's sequence survives rotations, reverse complements and both (tested on 5 records x 12 shifts).
  - A previous run's 04.final_gb/ or reoriented body in the same -o was published as this run's result (another
    genome under the record's header): the run clears its own intermediate folders first.
  - Curation: O2 cut complete genes across the origin in half (NCBI / MITOS2 / rotated input; ~25 genes in a rotation
    test), C1 / P4 judged a join across the origin by its hull 1..L and removed partial cox1 / nad5, `mitom curate`
    used code 5 for every record (vertebrate AGG stops read through; now each record's /transl_table), every
    /transl_except was deleted (an internal Trp codon gave X in the translation; now kept and applied), rRNA names
    l-rRNA / small subunit ribosomal RNA / RNR1 / RNR2 were not read, /gene was added to tRNAs / rRNAs without a
    gene feature (22 table2asn GeneXrefWithoutGene on NC_001636.1).
  - A tRNA across the origin got a truncated duplicate from MiTFi (vertebrate trnE cut to 32 of 69 nt; written
    twice with --trna mitfi); MitoM-ML ignored cmsearch / cmalign failures (0 tRNAs, empty anticodons: now an error,
    the core falls back to MiTFi with a warning); intron steps ignored failing exonerate / tblastn / blastn.
  - table2asn hung for ever on NC_012920.1 + NC_012920.1_2 in one file: gb2tbl renames duplicates ID-2 (never an ID
    already used) and every table2asn / core program call has a timeout.
  - Every command: a missing, empty, unreadable (truncated, not GenBank / FASTA), gzip, BOM, Latin-1 or EMBL input
    now gives one ERROR line and exit 2 (or is unpacked / converted) instead of a traceback or 'finished, 0 records'
    with exit 0 (dedup, reorient, gb2fasta, curate, check, compare, plot, validate, rrna-hmm, trna, intron, nuc,
    transfer_metadata, replace_organism); protein FASTA is not annotated as DNA; cut-off GenBank records (a killed
    download) are skipped and `mitom fetch` downloads them again; CONTIG records without sequence are skipped.
  - Same-name inputs no longer overwrite each other (dedup of LOCUS-only names, gb2fasta, reorient, plot,
    transfer_metadata / replace_organism, grouped runs); `reorient` exits 1 when files failed; group runs killed by a
    signal fail the run; one bad record no longer turns off --clade auto for all; /pseudo stays value-less and the
    NCR step writes atomically; `mitom compare` matches NCBI's complement(join(1,16306..16371)) and two files of any
    name.
- Every record of a multi-record GenBank file (the output of `mitom fetch`) took the lineage of the FIRST record
  (mitom_stage.sample_lineage): Metridium senile (Anthozoa, circular) after Aurelia aurita was made linear as a
  jellyfish; the MitoM-ML tRNA clade models used the wrong lineage too. Now the record's own lineage (read once per file).
- Stress test (2026-10-31): alignment gaps ('-', '.') and '*' in a sequence made exonerate fail on every protein gene
  (they are removed, odd letters read as N, with a warning); empty and all-N sequences were annotated (an all-N
  sequence got NCBI submission files) - now skipped; a run with nothing to annotate (empty, unreadable or all-N input)
  ended 'finished' with exit 0 - now an error (exit 2); a missing input file under --clade auto gave a traceback;
  .gz input and Latin-1 headers (Windows programs) ended with 0 records - now unpacked / converted (out/00.input);
  a second FASTA sequence with the same name was renamed silently - now a warning.
- rRNA fallback (mito_rrna_hmm.py) for a record without any phylum (unlabelled FASTA of uncertain phylum): the
  library whose 16S / 12S models fit best is used instead of no search (two ctenophores had lost both rRNAs).
- End check: a spliced join across the junction (66e0e09) is accepted only when the part it adds aligns by itself on
  the far side (>= half its residues, >= half the score per residue): contigs were joined to chance matches (atp6
  with an 800-nt 'intron', annelid nad4 with an 11-kb one). A gene that stops short of the junction with the residues
  towards it missing is checked spliced too (an intron between the exon and the junction: placozoan nad5). Intron
  benchmark: the four rotations of 187 records give the unrotated introns (300 of 300 rows each).
- Threads per record lowered so that the records annotated at once fit the cores (16 x 48 threads on 64 cores ran one
  record at a time: 148 records in 57 min instead of about 4).

### Added
- Documentation kept up with the code: `docs/OPTIONS.md` lists every option of every command as its `-h` prints it
  (`scripts/make_options_doc.py`; `tests/test_docs.py` fails while an option or a command is missing from the docs).
  The manual, methods and README now cover the `joined_intron` end check, the vertebrate stop rule, the learned start
  codon (P11), non-coding regions, the new map drawing (hatched introns, NCR shading) and gb2tbl's stress test; the
  example maps in `docs/images/` are redrawn by the current code, with a new one for introns and NCRs
  (`map_introns.png`, *Metridium senile*).
- Non-coding regions (mito_ncr.py; `run` / `annotate` / `import` `--ncr yes|no`, default yes, `--ncr_min_len`,
  default 100): every stretch of >= 100 nt no feature covers is written as misc_feature
  /note="non-coding region NCR<k>, <length> bp: no gene annotated (MitoM)" (the unnamed control region, sequence that
  may hold unidentified genes - it showed in neither the GenBank record nor the map); a joined gene covers its
  introns. table2asn: no new message (22 MJ20260820041 records: 0 ERROR, the same 18 WARNING). 20 NCBI Hesionidae
  mitogenomes: 33 NCRs, 24,280 nt; where NCBI annotated the record, the same stretches (Sirsoe methanicola 2,641 nt,
  one with a 288-codon ORF), plus NCBI's control-region features.
- Genome maps: introns hatched and labelled with gene and length (intron features, or the gaps of >= 142 nt of a
  join); non-coding regions shaded and labelled NCR<k> as in the GenBank (`mitom plot --min_noncoding`, 0: none);
  the name of a gene with an intron on its largest exon. On the 270 records of the intron benchmark the maps find
  the 311 introns eval_introns.py counts.

### Fixed
- Genome maps: a joined feature spans its first to its last part in transcript order. An intron NCBI gives by its
  two ends that holds other genes (hexacoral nad5, 10 kb of 17 kb) was taken the short way round the circle, and
  the line joining the gene's exons with it.
- Runs of several records were 4-100x slow: NumPy / OpenBLAS / OpenMP (tRNA classifier) started one thread per
  machine core in every process. Capped to 1 thread unless the user set them (mito_env): 22 annelid mitogenomes in
  154 s (7 s each), 4 records 47 min -> 34 s on a 128-core machine.
- A gene with an intron across the origin of a circular record (record started in exon 2: annelid cox1 group II
  introns, MJ20260820041 CHY / *Travisia*, CWSC / Nephtyidae): `--topology auto` saw exon 2 stop at the junction
  (protein2dna has no introns) and made the record linear, so the intron was never annotated (cox1 = exon 2 only:
  1,287 or 848 nt). The end check now aligns the gene again with exonerate protein2genome when it looks cut
  (verdict `joined_intron` -> circular); the de novo intron's description is kept on the CDS while the intron
  crosses the origin and the intron feature is restored after the record is rotated to the gene's first exon
  (mitom_orient.add_exon_gap_introns). 22 MJ mitogenomes: 13 CHY + 5 CWSC records gain their cox1 intron (checked
  with the spliced protein, read depth and junction reads, and LC677172 / Nephtys cox1 introns), the 4 others are
  unchanged; table2asn 0 errors.
- The same check missed introns of 6 kb or more lying wholly on one side of the junction (hexacoral nad5 intron,
  10-12 kb, holding other genes): the spliced re-check now aligns on the whole rest of the circle before the gene,
  then after it (no intron length limit). Benchmark (193 NCBI records with 309 introns: Cnidaria, Annelida, Placozoa,
  Porifera; 76 without; intron sites of the record's genus left out): records rotated so the origin lies in an exon,
  in the middle of an intron or 15 nt inside either end of it give the same introns as the unrotated records (300 of
  300 rows); before the fix 15 nt inside the right / left end: 232 / 190 of 295 introns found.

### Added
- `mitom gb2fasta`: GenBank (.gb .gbk .gbf .gbff .genbank) to FASTA, one file per input or per record (`--split`);
  headers `>ACCESSION.V definition topology=circular|linear`, ready for `mitom run`.
- `mitom dedup` takes FASTA too: identical sequences, and versions of one accession when the header's first word is
  an accession (MK123456.1); local names such as NODE_1 are never merged on their name.
- `mitom fetch --fasta`: a FASTA of the downloaded records next to the GenBank file.
- `mitom gb2tbl` (`mito_tbl.py`): any annotated GenBank (files and folders, any mix) -> NCBI submission files, the
  5-column feature table `.tbl` and the `.fsa` with source modifiers: every record in one pair (default), one pair
  per input file (`--split file`) or per record (`--split record`). A record repeated across inputs is written once,
  clashing IDs are renamed ID_2 ..., IDs cleaned to NCBI's characters; location (mitochondrion / chloroplast / nucleus),
  genetic code, taxid and source qualifiers (isolate, voucher, country -> geo_loc_name, collection date ...) in the
  FASTA headers; `--validate` runs table2asn (`.sqn`, `.val`, errors per record in gb2tbl_report.tsv; `--sbt`).
  Cleaned for NCBI and reported: qualifiers outside INSDC (Geneious / SnapGene / ApE) left out, non-INSDC feature keys
  as misc_feature, text made ASCII, features past the end and 100^101 sites flagged; unreadable files and records
  without sequence skipped with a message; `--genetic_code` (else the phylum's code when no /transl_table). Stress
  test against NCBI's gbf2tbl (14 cases, table2asn): 0 ERROR in every valid case vs 16 for gbf2tbl
  (benchmark/results/gb2tbl_stress_2026-10-30/).
  Taken over from the other converters: EMBL / ENA input and `-i -` for a pipe (GB2sequin, gbf2tbl); `--source_table`
  (one row per record, NCBI .src or CSV; columns by modifier name, others reported) and `--modifiers` for all records
  (GB2sequin's source modifier table, table2asn -src-file / -j); missing /product of mitochondrial genes filled from
  the gene name or /anticodon, `--product_table` (GB2sequin's gene product table); the DEFINITION as the FASTA title
  and 70-column sequence lines (gbf2fsa); Windows-1252 / Latin-1 files read. Stress test now 17 cases: 0 ERROR in all 15 valid ones
  vs 51 for gbf2tbl + gbf2fsa, which read no EMBL.
- `--topology auto|label|circular|linear` (`mitom run`, `mitom annotate`; default auto): records are annotated as
  labelled (GenBank LOCUS line; FASTA `topology=`; else one sequence >= 10 kb circular, several or shorter linear)
  and the protein genes at the two ends are aligned to their reference protein across the junction
  (`mito_core.py --check_ends`): a gene running from the end across to the start makes a record circular, a gene cut
  at an end of a circular-labelled record makes it linear, and circular-labelled jellyfish / hydroids / calcareous
  sponges (linear mtDNA) are annotated as linear unless a gene runs across their ends. Changed records are annotated
  again; `reports/topology_decision.tsv` gives each record's label, topology and reason.

- Genome maps in a new style, the default (`mitom plot --style circos`, `mitom run --plot_style circos`), after
  mitoPlotR (K. Onthank, MIT licence): a framed lane per strand, blocks coloured by type, names written along the
  circle, a position axis, GC content, read depth in 25 nt bins scaled to the 99th percentile with the mean as a
  dashed line, the organism in the centre; tRNA codes added and colour-blind-safe hues. The previous maps:
  `--style mitom` / `--plot_style mitom`. Examples: docs/images/.

### Changed
- `dedup.fasta` headers carry `topology=circular|linear` (without it `mitom run` takes the records of a multi-record
  FASTA for linear contigs).
- Help and docs in plain words: every command with a usage section (docs/USAGE_zh.md sections 5-6), what dedup
  does, a note that MitoM is still being tested and used inside the group.
- Annelids (`--clade` Annelida): a protein gene whose first complete stop lies 11-30 nt inside the next same-strand
  protein gene ends with T / TA before it, as at the sequence origin (`mito_curate.py --abbreviated_pcg_overlap`;
  NCBI annelid records 142 vs 43, mostly nad2 before cox1 and nad6 before cob; other phyla keep the complete stop:
  Cnidaria 109 vs 5, Mollusca 66 vs 12). Curation of 1,718 leave-out benchmark records: annelid protein genes with
  both ends as NCBI 7,357 -> 7,414, no other phylum changed. The same nad2 no longer depends on whether the record
  was rotated.

- Start codons (`mito_curate.py` P11, `mito_start.py`): after P1, the in-frame start codon NCBI records would choose
  is taken - gradient-boosted trees learned from 1,683 leave-out benchmark records (codon, distance to the preceding
  feature, the core's alignment, gene, phylum), stored as plain JSON trees (models/start_model/, no pickle). Exact
  starts, cross-validated by family: 89.1 -> 90.7 % (arthropods 77.8 -> 82.3 %); `--start_model no` turns it off;
  retrain with benchmark/train_start_model.py.
- Vertebrates: a complete stop running 1-2 nt into the next protein gene is kept (atp6's TAA sharing its A with
  cox3's ATG); `--abbreviated_overlap` (T / TA) applies before tRNAs only. 1,719 leave-out records: protein genes with
  both ends as NCBI 17,727 -> 17,951.

### Fixed
- Genome maps: no position tick right before the origin (16 kb of a 16,063 bp record covered "0").
- The GenBank header copied from the input record keeps the topology MitoM annotated in its LOCUS line (it took the
  input's back, which every later step reads).

## 1.0.0 - 2026-10-30

### Added
- `--engine mitoz|mitos2` (`mitom run`, `mitom annotate`): MitoZ 3.6 or MITOS2 2.1.10 find the genes when they are
  installed (they are not installed with MitoM); every later MitoM step runs on their annotation (reorientation,
  round 2, curation and NCBI names, NCBI submission files, genome maps). `mito_engines.py`.
- `mitom import`: existing MitoZ results (result folders or .gbf files), MITOS2 result folders (result.bed and the
  genome FASTA) or any annotated GenBank -> the same MitoM steps, without annotating again.
- Circular genome maps (`mito_plot.py`, `--plot`, default yes): genes by strand and class, GC content, GC skew; read
  depth with `--infq` (bwa; minimap2 when it is installed, about 8x faster).
- `mitom setup`: downloads the MITOS2 reference data refseq89m (md5 checked) and writes the config file.
- `mitom --version`, `CHANGELOG.md`, a bioconda recipe
  (`packaging/bioconda/`).
- Introns from the intron sites of annotated relatives (`mito_intron_ref.py`, run after `mito_intron.py`; report
  `intron_ref_report.tsv`): the site library `models/intron_sites.*` (79 KB) is built from a gold standard of all 193
  local NCBI records with introns, each NCBI intron checked (spliced product, junction against homologs, group I / II
  splice rules, related records): 309 introns, 97 placed differently from NCBI (`benchmark/results/
  intron_gold_2026-10-29/`). `mitom run` against it: 309 / 309 exact with the whole library, 299 / 309 when the
  record is left out (a new genome), no false intron on 76 intron-free records, no NCBI ERROR.
- `benchmark/eval_introns.py` (introns vs NCBI or the gold standard, nested introns, per gold class) and
  `benchmark/bench_intron_ref.py` (own / leave-one-out / leave-species / leave-genus).
- `mitom debug <output folder | log file>` (`mito_debug.py`): what happened in a run, from its log - status (finished
  with its exit code, stopped on an error with the step and the end of the traceback, or did not finish), the time
  of each step and of its parts, every record that failed or was left out with the first error of its own log,
  every program that failed (exit code, time, folder, end of its stderr; with `--debug` also the programs that
  failed inside the annotation core), missing programs and how to install them, NCBI ERRORs, warnings grouped by
  kind, and commands to reproduce (the run in debug mode, the failed program by hand, a failed record alone).
  Exit code 0 (no problem), 1 (problems found), 2 (no MitoM log). Reads the logs of MitoM 1.0.0 too.
- `--debug` (every `mitom` command; also an option of `run`, `import`, `annotate`, `reorient`): DEBUG lines in the
  log; the environment (MitoM version and git commit, the command as typed, Python, the machine, every program with
  its version, the Python packages, the config file, reference and data folders, environment variables, memory,
  free disk space); the output of every program each step runs; intermediate files kept (`work/`).
- `mitom plot`: the circular genome maps of `run` as a command of their own (`mito_plot.py`); a folder input now
  includes its subfolders and every GenBank suffix (.gb .gbk .gbf .genbank .gbff), as the other step commands.
- `mitom intron --method both | denovo | ref`: `mitom intron` now adds introns as `run` does, de novo
  (`mito_intron.py`) and then at the intron sites of annotated relatives (`mito_intron_ref.py`, report
  `intron_ref_report.tsv`); `--method denovo` is the former `mitom intron`; options of `mito_intron.py` (`--lib`,
  `--ref_gb`, `--maxintron` ...) are passed on.

### Changed
- One conda environment (`env/mitotools.yml`) holds everything, the Python packages included (Biopython,
  scikit-learn 1.7.2, matplotlib); `mitom` runs on that environment's Python and uses its programs first; nothing is
  installed into the system Python any more.
- The main program is `mitom_pipeline.py` (was `batch_mitoz.py`; the old name still works).
- Work folders `01.anno_round1/`, `03.anno_round2_reoriented/`; default output folder `mitom_out`.
- The Docker image is built from that one environment, so it includes scikit-learn (MitoM-ML tRNAs instead of the
  MiTFi fallback).
- Java (MiTFi; its classes are Java 5, so any Java >= 8 runs it): the computer's own Java is used; install.sh adds
  openjdk 8 to the environment only when there is none (openjdk 25 alone was 0.33 GB). Environment 1.6 GB, 157 packages.
- Step 6 of `run` prints what the intron step did (introns written, low-confidence candidates not written, or skipped
  because exonerate / blastp were not found - now a warning on the console, before only in run.log).

- Second copies of a protein gene are kept when they are >= 95 % identical to it (length within 3 nt): the core keeps
  a single such copy (`COPY_IDENTITY`; before only copies in a block of >= 2 genes), and curation rule D1 no longer
  removes them (before it needed both copies to end with a stop, which the core's genes lack at that point).
  NCBI annotates such copies (heron / cormorant / turaco nad6, frog cob, squid cox2): on the records concerned
  MitoM gains 9 genes NCBI has, 1 it does not (NC_049135.1 atp8, 98 %); imported MitoZ / MITOS2 results gain 9 / 7,
  all in NCBI; nothing else changed.
- `mito_intron.py`: a gene across the origin of a circular record is modelled on a rotated copy (the cox1-first
  rotation put the origin inside annelid cox1 introns, half of them were missed); exons chained in protein order only;
  a model whose intron holds a copy of its own exons is refused (duplicated genes, Hexamermis atp6); trans-spliced
  genes written as partial CDS pieces; partial core genes searched with the core's reference protein; intron ends
  refined independently by the median score over up to five homologs; introns next to an inserted mini-exon placed
  again. Against the gold standard 285 -> 296 of 309 exact on its own.
- tRNA step: the sensitive pass for MiTFi models without a fast-pass hit runs `cmsearch --nohmm --F6 3e-4 --notrunc`
  outside the interiors of the core's protein genes / rRNAs instead of `--max` on the whole record: NC_014885.1
  297 -> 18 s (2 threads), 542 benchmark records 11.9 -> 1.0 h CPU; recall and ends unchanged, 14 weak calls inside
  genes / rRNAs (none in NCBI) gone.
- `mitom nuc annotate` on whole-genome assemblies: 6 s - 2.2 min at 8 threads on 30 SPAdes assemblies up to 1.7 Gb
  (before up to 22 min and 16 GB): low-complexity masking in the prefilter (microsatellites in the reference spacers
  had selected ~10,000 junk contigs), per-gene candidate regions, parallel BLAST database, indexed FASTA, one tblastn
  per distinct histone protein, fixed work limits for assemblies; partial H3 / H4 models below 80 % identity are
  no longer reported (degenerate H3 retrocopies: SRR2006110 3,688 -> 107 H3). The NCBI benchmark is unchanged.
- After the intron step a circular record is rotated again so that the orientation gene's first exon starts it;
  rotation keeps partial ends (< / >; dropped before, which gave table2asn ERRORs on trans-spliced pieces).
- The run log (`run.log`, `pipeline.log` of `annotate`) has lines of a fixed form that `mitom debug` reads
  (`mitom_log.py`): `[STEP 6/11 add_missing] start | done | failed <seconds>s`, `[CMD rc=<exit> <seconds>s]
  <command>` for every program the workflow runs (always, with the end of its stderr when it failed),
  `[RECORD <name>] ...` for every record that failed or was left out, `[ERROR]` with the whole traceback when a run
  stops on an error (before, the traceback went to the terminal only; now the terminal names `mitom debug OUT`),
  `[DONE] exit=<code> <time>`, `[ENV]` lines with the MitoM version and the command. The programs of the de novo
  intron step and of the Snakemake backend are logged the same way.
- One spelling per option, the old one still accepted: `-o / --out` (was `--out_root`, `--out_dir`), `--jobs` (was
  `--max_tasks`), `--refdir` (was `--mitos_refdir`), `--threads` (was `--cpu`), `--manifest_only` (was
  `--manifest-only`).
- Exit codes: 0 done, 1 failure (a program, the reference data or a Python package missing; every record failed),
  2 wrong command line or input. `mitom annotate` exits 1 when every record failed, as `run` does.
- The step commands (`curate`, `check`, `compare`, `trna`, `validate`, `rescue`, `rrna-hmm search`) print their
  one-line summary on stderr; stdout holds data only.
- Report tables have LF line ends (the module reports had CRLF) and always a header line. `run` / `import` /
  `annotate` write every report into `reports/` when the step runs; `reports/rrna_hmm_report.tsv` and
  `ncbi_submission/validation_messages.tsv` exist (header only) also when nothing was added or reported; `mitom
  annotate` now keeps `intron_report.tsv`, `intron_ref_report.tsv`, `intron_fix.log`, `run_params.json`,
  `curation_report.tsv` and `failed_logs/` in `reports/`; `no_<gene>_samples.txt` is in `reports/` (it was left at
  the top of the output folder).
- `mitom doctor` lists every program each step needs (new rows cmalign, makeblastdb, hmmbuild; tblastn / blastn for
  the intron-site transfer); the check before a run covers the steps the run will do and names the step.
- MitoM-ML unpacks its clade models and trains its classifiers in a user data folder (`$MITOM_DATA`, else
  `$XDG_DATA_HOME/mitom`, else `~/.local/share/mitom`), prepared by `mitom setup`, instead of the install folder.
- GenBank files rewritten in place (rescue, introns, rotation after the intron step) are written atomically.
- run.log no longer holds the INFO lines of fontTools / matplotlib (about 40 per genome map).
- The hint of `rescue_report.tsv` no longer compares MITOS2 with MitoZ.
- Every external program is run through one function that logs it at DEBUG level (command, folder, exit code,
  seconds, end of stderr on failure); the foundation of a debug mode.
- Code layout: `mitom_pipeline.py` is split into `mitom_inputs.py`, `mitom_stage.py`, `mitom_snakemake.py`,
  `mitom_orient.py`, `mitom_metadata.py`, `mitom_output.py`, `mitom_steps.py`, `mitom_log.py`, `mito_taxa.py`
  (every name still importable from `mitom_pipeline`); gene names in one module (`mito_names.py`); shared helpers in
  `mito_common.py`; one MITOS2 command builder and result.bed reader (`mito_engines.py`).

### Fixed
- `mitom compare` and `mitom check` crashed in the bioconda package and the Docker image (they imported from
  `benchmark/`, which is not shipped), and step 9 of `run` (original vs MitoM) then failed silently.
- A record turned to the plus strand (orientation gene on the minus strand) lost the partial ends (< / >) of its
  genes (table2asn StartCodon / NoStop errors).
- `/transl_except` and `/anticodon` positions did not follow a rotated or reverse-complemented record (intron
  CDS, NCBI records given to `mitom import`).
- MitoM-ML wrote into the install folder at its first run (a read-only install failed; parallel first runs wrote
  the same files); standalone `trna_ml.py` did not find refseq89m given the folder that holds it.
- `ncbi_validate` (`validation.tsv` column `ours`, `--only_ours`) did not count messages on the intron steps'
  features as MitoM's.
- `mitom doctor` said "All required programs found" without cmalign (MitoM-ML tRNAs) or makeblastdb (introns).
- `mito_intron.py` on a record whose introns were annotated already (an NCBI record, MitoM output) wrote a second
  intron feature for each intron.
- A folder whose name ends in .gb was taken for a GenBank file by compare / curate / validate / the intron-site
  transfer.

### Removed
- MitoZ as a part of MitoM: `--core mitoz|hybrid`, `--trna`, `--circos`, `--conda_env_mitoz`, `--mitoz_path`,
  `--mitoz_extra`, `env/mitoz.yml`, `scripts/mitoz_fast.py` (MitoZ can still be the engine, see above).

## 0.9.0 - 2026-10-27
- NCBI-submittable output (no table2asn ERROR on 2,497 benchmark records), MiTFi run as MITOS2 runs it, tRNA second
  opinion (MitoM-ML with MiTFi), `--rrna_method hmm|cm`; MitoM's own annotation core `mito_core.py` as the default.
