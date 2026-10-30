# Changelog

Notable changes of MitoM. Versions follow semantic versioning (major.minor.patch); `mitom --version` prints the
installed one.

## Unreleased

### Added
- `mitom gb2fasta`: GenBank (.gb .gbk .gbf .gbff .genbank) to FASTA, one file per input or per record (`--split`);
  headers `>ACCESSION.V definition topology=circular|linear`, ready for `mitom run`.
- `mitom dedup` takes FASTA too: identical sequences, and versions of one accession when the header's first word is
  an accession (MK123456.1); local names such as NODE_1 are never merged on their name.
- `mitom fetch --fasta`: a FASTA of the downloaded records next to the GenBank file.
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
