# MitoM 全部命令与参数 / All commands and options

由 `python3 scripts/make_options_doc.py` 从每个命令的 `-h` 生成（命令行改动后重新生成；`tests/test_docs.py` 检查这里有没有漏掉参数）。各命令的用法、示例和输出见 [USAGE_zh.md](USAGE_zh.md)。

Generated from each command's `-h` by `scripts/make_options_doc.py`; usage, examples and outputs: [USAGE_zh.md](USAGE_zh.md).

| 命令 / command | 程序 / script |
|---|---|
| [`mitom run`](#mitom-run) | `mitom_pipeline.py run` |
| [`mitom import`](#mitom-import) | `mitom_pipeline.py import` |
| [`mitom annotate`](#mitom-annotate) | `mitom_pipeline.py annotate` |
| [`mitom reorient`](#mitom-reorient) | `mitom_pipeline.py reorient` |
| [`mitom transfer_metadata`](#mitom-transfer_metadata) | `mitom_pipeline.py transfer_metadata` |
| [`mitom replace_organism`](#mitom-replace_organism) | `mitom_pipeline.py replace_organism` |
| [`mitom dedup`](#mitom-dedup) | `mitom_pipeline.py dedup` |
| [`mitom gb2fasta`](#mitom-gb2fasta) | `mitom_pipeline.py gb2fasta` |
| [`mitom setup`](#mitom-setup) | `mito_env.py setup` |
| [`mitom doctor`](#mitom-doctor) | `mito_env.py doctor` |
| [`mitom curate`](#mitom-curate) | `mito_curate.py` |
| [`mitom intron`](#mitom-intron) | `mitom_pipeline.py intron` |
| [`mitom trna`](#mitom-trna) | `mito_trna.py` |
| [`mitom rescue`](#mitom-rescue) | `mito_rescue.py` |
| [`mitom rrna-hmm`](#mitom-rrna-hmm) | `mito_rrna_hmm.py` |
| [`mitom validate`](#mitom-validate) | `ncbi_validate.py` |
| [`mitom gb2tbl`](#mitom-gb2tbl) | `mito_tbl.py` |
| [`mitom fetch`](#mitom-fetch) | `mito_fetch.py` |
| [`mitom check`](#mitom-check) | `mito_check.py` |
| [`mitom compare`](#mitom-compare) | `mito_compare.py` |
| [`mitom plot`](#mitom-plot) | `mito_plot.py` |
| [`mitom nuc`](#mitom-nuc) | `mito_nuc.py` |
| [`mitom debug`](#mitom-debug) | `mito_debug.py` |

## mitom run

```text
usage: mitom run [-h] -i INPUT [-o OUT] [--if_split_fasta {yes,no}] [--id_maxlen ID_MAXLEN]
                 [--which_gene_first WHICH_GENE_FIRST] [--fallback_genes FALLBACK_GENES] [--force_forward]
                 [--no_force_forward] [--debug] [--mitos_bin MITOS_BIN] [--refdir DIR] [--plot {yes,no}]
                 [--plot_style {circos,mitom}] [--ncr {yes,no}] [--ncr_min_len NCR_MIN_LEN]
                 [--curate {yes,no}] [--sbt SBT] [--rescue {auto,yes,no}] [--rrna_hmm {auto,yes,no}]
                 [--rrna_hmm_lib RRNA_HMM_LIB] [--rrna_method {hmm,cm}] [--phylum PHYLUM]
                 [--intron_fix {auto,yes,no}] [--intron_bin INTRON_BIN] [--intron_lib INTRON_LIB]
                 [--intron_ref_gb INTRON_REF_GB] [--keep_intermediate {yes,no}] [--threads THREADS]
                 [--jobs N] [--clade CLADE] [--genetic_code GENETIC_CODE] [--infq INFQ]
                 [--fq_position {flat,sample_dir}] [--suffix_fq SUFFIX_FQ]
                 [--scheduler_backend {auto,snakemake,threads}] [--snakemake_bin SNAKEMAKE_BIN]
                 [--smk_cores SMK_CORES] [--smk_mem_mb SMK_MEM_MB] [--smk_latency_wait SMK_LATENCY_WAIT]
                 [--smk_extra FLAG] [--engine {native,mitoz,mitos2}] [--mitoz_bin MITOZ_BIN]
                 [--topology {auto,label,circular,linear}] [--if_reannotation_2 {yes,auto,no}]
                 [--preserve_phase1_metadata {yes,no}] [--smart_metadata {yes,no}]
                 [--auto_accept_metadata {yes,no}]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTA file/dir or GenBank file/dir.
  -o OUT, --out OUT     Output folder (default: mitom_out).
  --if_split_fasta {yes,no}
  --id_maxlen ID_MAXLEN
                        Max length of internal sample IDs (default: 15).
  --which_gene_first WHICH_GENE_FIRST
                        Gene used as sequence start (default: cox1).
  --fallback_genes FALLBACK_GENES
                        Genes tried in order when --which_gene_first is absent (default: cox2,cox3,cob,nad1;
                        'none' = leave such records unrotated). Records oriented this way are listed in
                        oriented_by_fallback.txt.
  --force_forward       RC the sequence when target gene is on minus strand (default: on).
  --no_force_forward
  --debug               Debug mode: DEBUG lines, program versions and the environment in the log; the output
                        of every program in the logs; intermediate files kept (--keep_intermediate yes).
                        `mitom debug OUT` reads the log.
  --topology {auto,label,circular,linear}
                        Circular or linear. auto (default): as labelled - GenBank LOCUS line; FASTA header
                        topology=circular|linear; header words ('circular', 'complete genome' / 'partial
                        genome'); an assembler's terminal overlap (removed): circular; else a sequence of >=
                        10 kb circular (one or several in a file), a shorter one linear - and the protein
                        genes at the two ends checked against their reference proteins: a gene running from
                        the end across to the start makes the record circular, a gene cut at an end of a
                        circular-labelled record makes it linear, and so do jellyfish / hydroids /
                        calcareous sponges (linear mtDNA) without a gene across the ends; an unlabelled
                        record no gene at the ends decides is circular with >= 12 of the 13 protein genes
                        and both rRNAs, else linear; changed records are annotated again
                        (reports/topology_decision.tsv). label: the label only. circular / linear: every
                        record.
  --if_reannotation_2 {yes,auto,no}
                        Phase-2 re-annotation of the reoriented sequences. auto (default): round 1 runs in
                        full and only samples that fail the check after the cox1 shift (gene at the original
                        sequence ends, standard gene missing, gene changed by the shift, gene across the new
                        origin) are re-annotated, decisions in reannotation_decision.tsv; yes: every sample;
                        no: never.
  --preserve_phase1_metadata {yes,no}
                        Restore original GenBank header metadata into final files (default: yes; requires
                        GenBank input).
  --smart_metadata {yes,no}
                        Enable fuzzy/interactive accession matching (default: no).
  --auto_accept_metadata {yes,no}
                        Auto-accept best fuzzy match (requires --smart_metadata yes).

Annotation:
  --mitos_bin MITOS_BIN
                        MITOS2 environment bin (runmitos). Default: found automatically (mito_env.py).
  --refdir DIR          Folder that holds the MITOS2 reference data refseq89m. Default: found automatically.
  --plot {yes,no}       Circular genome map per record (mito_plot.py, matplotlib) -> plots/<record>.png and
                        .pdf: genes by strand and type, GC content; with --infq also the read depth
                        (minimap2 when installed, else bwa mem). Default yes (about 2 s per record, --jobs
                        in parallel).
  --plot_style {circos,mitom}
                        Map style: circos (default, after mitoPlotR) = a framed lane per strand, names along
                        the circle, a position axis, GC content, read depth; mitom = genes on both sides of
                        one circle, GC content and GC skew rings.
  --ncr {yes,no}        yes (default): every stretch of >= --ncr_min_len nt no feature covers becomes a
                        misc_feature /note="non-coding region NCR<k>, <length> bp: no gene annotated
                        (MitoM)" in the GenBank and NCBI files, shaded and labelled on the maps (unnamed
                        control region, sequence that may hold unidentified genes); no: none (mito_ncr.py)
  --ncr_min_len NCR_MIN_LEN
                        shortest non-coding region written (nt, default 100; intergenic spacers are shorter)
  --curate {yes,no}     yes (default): mito_curate.py boundary rules and NCBI gene names -> 05.curated/, and
                        NCBI submission files (.fsa/.tbl/.sqn/.val) -> 06.ncbi_submission/
  --sbt SBT             NCBI submitter template (template.sbt from submit.ncbi.nlm.nih.gov/genbank/template)
                        for the .sqn files
  --rescue {auto,yes,no}
                        Add protein genes / rRNAs the core missed from a MITOS2 annotation (mito_rescue.py;
                        nothing annotated is changed; needs the optional MITOS2 installation). auto
                        (default) / no = not run (the core retries other loci itself); yes = run for records
                        missing a gene.
  --rrna_hmm {auto,yes,no}
                        Last-resort 12S/16S search for rRNAs still missing after the cores and the rescue,
                        with per-phylum profile HMMs built from NCBI rRNAs (mito_rrna_hmm.py; only missing
                        rRNAs are searched). auto = run when the library and nhmmer are found (default).
  --rrna_hmm_lib RRNA_HMM_LIB
                        Per-phylum rRNA HMM library (mito_rrna_hmm.py build). Default: models/rrna_hmm.
  --rrna_method {hmm,cm}
                        12S/16S search of the core: hmm (default, faster) = per-phylum profile HMMs,
                        MITOS2's covariance models where no model or hit; cm = MITOS2's covariance models
                        only. Same accuracy (2026-10-28, 450 leave-out records: found 98.5 / 98.5 %, both
                        ends exact 47.4 / 47.6 %); hmm takes seconds and ~0.1 GB, cm ~10 s and ~0.5 GB per
                        gene.
  --phylum PHYLUM       Phylum of all input records for --rrna_hmm (default: the phylum found by --clade
                        auto, or the phylum of an explicit --clade; with --clade Animal and FASTA input give
                        it here, otherwise the record's own taxonomy is used when it has one).
  --intron_fix {auto,yes,no}
                        Add introns with mito_intron.py (MFannot method: Exonerate protein2genome + MFannot
                        boundary rules + Rfam intron type), then with mito_intron_ref.py (intron sites of
                        annotated relatives, models/intron_sites). auto = run when exonerate is found
                        (default).
  --intron_bin INTRON_BIN
                        Directory containing exonerate, cmsearch, blastp, tblastn and blastn (e.g. a conda
                        env bin).
  --intron_lib INTRON_LIB
                        Protein library for intron models in MFannot header format ('>Species mt cox1;
                        ...'), repeatable. Default: MFannot's curated collections bundled in
                        models/mfannot_proteins.
  --intron_ref_gb INTRON_REF_GB
                        GenBank file/dir of close relatives used as extra homologs for intron models
                        (repeatable).
  --keep_intermediate {yes,no}
                        Keep the intermediate folders (00.*, 01.*, 03.*, 04.*, 05.*, _snakemake_*; run:
                        moved to work/). Default no: deleted after the final output is written; core logs of
                        problem samples are always kept (reports/failed_logs/).
  --threads THREADS     Threads per annotation job (default: min(8, number of CPUs) = 8).
  --jobs N              Records annotated at the same time (default: CPUs // threads = 9).
  --clade CLADE         Taxonomic group (sets the genetic code and group rules). auto (default): per record
                        - the taxonomy and CDS transl_table of a GenBank record; FASTA and records without a
                        taxonomy: blastx to the reference proteins (mito_classify.py: phylum of the closest
                        genomes, their genetic code checked on the codons; no close relative: Animal and the
                        code of the codons). Several phyla / codes: `run` annotates each group with its own
                        clade and code and merges the results (reports/clade_decision.tsv); `import` and
                        `annotate` need one. Or give one of: Chordata, Arthropoda, Echinodermata, Annelida-
                        segmented-worms, Bryozoa, Mollusca, Nematoda, Nemertea-ribbon-worms, Porifera-
                        sponges, Cnidaria, Platyhelminthes-flatworms, Chaetognatha-arrow-worms, Animal.
  --genetic_code GENETIC_CODE
                        Genetic code (default with --clade auto: the CDS transl_table of each record, else
                        the code of its lineage or of the closest reference genomes; otherwise the code of
                        the clade - 2 Chordata, 9 Echinodermata, 4 Porifera, else 5; required for Cnidaria,
                        Platyhelminthes, Chaetognatha, Animal).

FASTQ (optional - read depth ring of the genome maps, --plot):
  --infq INFQ           Directory with the paired FASTQ files of the samples (matched by sample name).
  --fq_position {flat,sample_dir}
  --suffix_fq SUFFIX_FQ

Scheduler:
  --scheduler_backend {auto,snakemake,threads}
                        Dispatch backend (default: auto). auto — use Snakemake if available, else threads
                        snakemake — always use Snakemake threads — always use ThreadPoolExecutor Annotation
                        is parallelised; finalization runs in the main process regardless of backend.
  --snakemake_bin SNAKEMAKE_BIN
  --smk_cores SMK_CORES
  --smk_mem_mb SMK_MEM_MB
  --smk_latency_wait SMK_LATENCY_WAIT
  --smk_extra FLAG      Extra Snakemake flags appended verbatim. Example: --smk_extra '--profile slurm'

Annotation engine:
  --engine {native,mitoz,mitos2}
                        Program that finds the genes: native (default) = MitoM's own core (mito_core.py);
                        mitoz = MitoZ 3.6, mitos2 = MITOS2 2.1.10 (not installed with MitoM - install them
                        yourself, see mito_engines.py; found on PATH or in a conda environment). Every later
                        MitoM step runs on their annotation: reorientation, round 2, curation, NCBI files,
                        genome maps.
  --mitoz_bin MITOZ_BIN
                        Folder with the mitoz program (--engine mitoz). Default: found automatically.
```

## mitom import

```text
usage: mitom import [-h] -i INPUT [-o OUT] [--if_split_fasta {yes,no}] [--id_maxlen ID_MAXLEN]
                    [--which_gene_first WHICH_GENE_FIRST] [--fallback_genes FALLBACK_GENES]
                    [--force_forward] [--no_force_forward] [--debug] [--mitos_bin MITOS_BIN] [--refdir DIR]
                    [--plot {yes,no}] [--plot_style {circos,mitom}] [--ncr {yes,no}]
                    [--ncr_min_len NCR_MIN_LEN] [--curate {yes,no}] [--sbt SBT] [--rescue {auto,yes,no}]
                    [--rrna_hmm {auto,yes,no}] [--rrna_hmm_lib RRNA_HMM_LIB] [--rrna_method {hmm,cm}]
                    [--phylum PHYLUM] [--intron_fix {auto,yes,no}] [--intron_bin INTRON_BIN]
                    [--intron_lib INTRON_LIB] [--intron_ref_gb INTRON_REF_GB] [--keep_intermediate {yes,no}]
                    [--threads THREADS] [--jobs N] [--clade CLADE] [--genetic_code GENETIC_CODE]
                    [--infq INFQ] [--fq_position {flat,sample_dir}] [--suffix_fq SUFFIX_FQ]
                    [--scheduler_backend {auto,snakemake,threads}] [--snakemake_bin SNAKEMAKE_BIN]
                    [--smk_cores SMK_CORES] [--smk_mem_mb SMK_MEM_MB] [--smk_latency_wait SMK_LATENCY_WAIT]
                    [--smk_extra FLAG] [--preserve_phase1_metadata {yes,no}] [--smart_metadata {yes,no}]
                    [--auto_accept_metadata {yes,no}] [--fasta FASTA] [--topology {auto,circular,linear}]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTA file/dir or GenBank file/dir.
  -o OUT, --out OUT     Output folder (default: mitom_out).
  --if_split_fasta {yes,no}
  --id_maxlen ID_MAXLEN
                        Max length of internal sample IDs (default: 15).
  --which_gene_first WHICH_GENE_FIRST
                        Gene used as sequence start (default: cox1).
  --fallback_genes FALLBACK_GENES
                        Genes tried in order when --which_gene_first is absent (default: cox2,cox3,cob,nad1;
                        'none' = leave such records unrotated). Records oriented this way are listed in
                        oriented_by_fallback.txt.
  --force_forward       RC the sequence when target gene is on minus strand (default: on).
  --no_force_forward
  --debug               Debug mode: DEBUG lines, program versions and the environment in the log; the output
                        of every program in the logs; intermediate files kept (--keep_intermediate yes).
                        `mitom debug OUT` reads the log.
  --preserve_phase1_metadata {yes,no}
                        Restore original GenBank header metadata into final files (default: yes; requires
                        GenBank input).
  --smart_metadata {yes,no}
                        Enable fuzzy/interactive accession matching (default: no).
  --auto_accept_metadata {yes,no}
                        Auto-accept best fuzzy match (requires --smart_metadata yes).
  --fasta FASTA         Genome FASTA of MITOS2 results (sequence named as in result.bed; repeatable).
                        Default: FASTA files in the MITOS2 result folder or the folder above it.
  --topology {auto,circular,linear}
                        auto (default): MITOS2 from result.gff, GenBank from the LOCUS line, MitoZ output
                        circular (MitoZ writes linear for every record; give linear for contigs).

Annotation:
  --mitos_bin MITOS_BIN
                        MITOS2 environment bin (runmitos). Default: found automatically (mito_env.py).
  --refdir DIR          Folder that holds the MITOS2 reference data refseq89m. Default: found automatically.
  --plot {yes,no}       Circular genome map per record (mito_plot.py, matplotlib) -> plots/<record>.png and
                        .pdf: genes by strand and type, GC content; with --infq also the read depth
                        (minimap2 when installed, else bwa mem). Default yes (about 2 s per record, --jobs
                        in parallel).
  --plot_style {circos,mitom}
                        Map style: circos (default, after mitoPlotR) = a framed lane per strand, names along
                        the circle, a position axis, GC content, read depth; mitom = genes on both sides of
                        one circle, GC content and GC skew rings.
  --ncr {yes,no}        yes (default): every stretch of >= --ncr_min_len nt no feature covers becomes a
                        misc_feature /note="non-coding region NCR<k>, <length> bp: no gene annotated
                        (MitoM)" in the GenBank and NCBI files, shaded and labelled on the maps (unnamed
                        control region, sequence that may hold unidentified genes); no: none (mito_ncr.py)
  --ncr_min_len NCR_MIN_LEN
                        shortest non-coding region written (nt, default 100; intergenic spacers are shorter)
  --curate {yes,no}     yes (default): mito_curate.py boundary rules and NCBI gene names -> 05.curated/, and
                        NCBI submission files (.fsa/.tbl/.sqn/.val) -> 06.ncbi_submission/
  --sbt SBT             NCBI submitter template (template.sbt from submit.ncbi.nlm.nih.gov/genbank/template)
                        for the .sqn files
  --rescue {auto,yes,no}
                        Add protein genes / rRNAs the core missed from a MITOS2 annotation (mito_rescue.py;
                        nothing annotated is changed; needs the optional MITOS2 installation). auto
                        (default) / no = not run (the core retries other loci itself); yes = run for records
                        missing a gene.
  --rrna_hmm {auto,yes,no}
                        Last-resort 12S/16S search for rRNAs still missing after the cores and the rescue,
                        with per-phylum profile HMMs built from NCBI rRNAs (mito_rrna_hmm.py; only missing
                        rRNAs are searched). auto = run when the library and nhmmer are found (default).
  --rrna_hmm_lib RRNA_HMM_LIB
                        Per-phylum rRNA HMM library (mito_rrna_hmm.py build). Default: models/rrna_hmm.
  --rrna_method {hmm,cm}
                        12S/16S search of the core: hmm (default, faster) = per-phylum profile HMMs,
                        MITOS2's covariance models where no model or hit; cm = MITOS2's covariance models
                        only. Same accuracy (2026-10-28, 450 leave-out records: found 98.5 / 98.5 %, both
                        ends exact 47.4 / 47.6 %); hmm takes seconds and ~0.1 GB, cm ~10 s and ~0.5 GB per
                        gene.
  --phylum PHYLUM       Phylum of all input records for --rrna_hmm (default: the phylum found by --clade
                        auto, or the phylum of an explicit --clade; with --clade Animal and FASTA input give
                        it here, otherwise the record's own taxonomy is used when it has one).
  --intron_fix {auto,yes,no}
                        Add introns with mito_intron.py (MFannot method: Exonerate protein2genome + MFannot
                        boundary rules + Rfam intron type), then with mito_intron_ref.py (intron sites of
                        annotated relatives, models/intron_sites). auto = run when exonerate is found
                        (default).
  --intron_bin INTRON_BIN
                        Directory containing exonerate, cmsearch, blastp, tblastn and blastn (e.g. a conda
                        env bin).
  --intron_lib INTRON_LIB
                        Protein library for intron models in MFannot header format ('>Species mt cox1;
                        ...'), repeatable. Default: MFannot's curated collections bundled in
                        models/mfannot_proteins.
  --intron_ref_gb INTRON_REF_GB
                        GenBank file/dir of close relatives used as extra homologs for intron models
                        (repeatable).
  --keep_intermediate {yes,no}
                        Keep the intermediate folders (00.*, 01.*, 03.*, 04.*, 05.*, _snakemake_*; run:
                        moved to work/). Default no: deleted after the final output is written; core logs of
                        problem samples are always kept (reports/failed_logs/).
  --threads THREADS     Threads per annotation job (default: min(8, number of CPUs) = 8).
  --jobs N              Records annotated at the same time (default: CPUs // threads = 9).
  --clade CLADE         Taxonomic group (sets the genetic code and group rules). auto (default): per record
                        - the taxonomy and CDS transl_table of a GenBank record; FASTA and records without a
                        taxonomy: blastx to the reference proteins (mito_classify.py: phylum of the closest
                        genomes, their genetic code checked on the codons; no close relative: Animal and the
                        code of the codons). Several phyla / codes: `run` annotates each group with its own
                        clade and code and merges the results (reports/clade_decision.tsv); `import` and
                        `annotate` need one. Or give one of: Chordata, Arthropoda, Echinodermata, Annelida-
                        segmented-worms, Bryozoa, Mollusca, Nematoda, Nemertea-ribbon-worms, Porifera-
                        sponges, Cnidaria, Platyhelminthes-flatworms, Chaetognatha-arrow-worms, Animal.
  --genetic_code GENETIC_CODE
                        Genetic code (default with --clade auto: the CDS transl_table of each record, else
                        the code of its lineage or of the closest reference genomes; otherwise the code of
                        the clade - 2 Chordata, 9 Echinodermata, 4 Porifera, else 5; required for Cnidaria,
                        Platyhelminthes, Chaetognatha, Animal).

FASTQ (optional - read depth ring of the genome maps, --plot):
  --infq INFQ           Directory with the paired FASTQ files of the samples (matched by sample name).
  --fq_position {flat,sample_dir}
  --suffix_fq SUFFIX_FQ

Scheduler:
  --scheduler_backend {auto,snakemake,threads}
                        Dispatch backend (default: auto). auto — use Snakemake if available, else threads
                        snakemake — always use Snakemake threads — always use ThreadPoolExecutor Annotation
                        is parallelised; finalization runs in the main process regardless of backend.
  --snakemake_bin SNAKEMAKE_BIN
  --smk_cores SMK_CORES
  --smk_mem_mb SMK_MEM_MB
  --smk_latency_wait SMK_LATENCY_WAIT
  --smk_extra FLAG      Extra Snakemake flags appended verbatim. Example: --smk_extra '--profile slurm'
```

## mitom annotate

```text
usage: mitom annotate [-h] -i INPUT [-o OUT] [--if_split_fasta {yes,no}] [--id_maxlen ID_MAXLEN]
                      [--which_gene_first WHICH_GENE_FIRST] [--fallback_genes FALLBACK_GENES]
                      [--force_forward] [--no_force_forward] [--debug] [--mitos_bin MITOS_BIN]
                      [--refdir DIR] [--plot {yes,no}] [--plot_style {circos,mitom}] [--ncr {yes,no}]
                      [--ncr_min_len NCR_MIN_LEN] [--curate {yes,no}] [--sbt SBT] [--rescue {auto,yes,no}]
                      [--rrna_hmm {auto,yes,no}] [--rrna_hmm_lib RRNA_HMM_LIB] [--rrna_method {hmm,cm}]
                      [--phylum PHYLUM] [--intron_fix {auto,yes,no}] [--intron_bin INTRON_BIN]
                      [--intron_lib INTRON_LIB] [--intron_ref_gb INTRON_REF_GB]
                      [--keep_intermediate {yes,no}] [--threads THREADS] [--jobs N] [--clade CLADE]
                      [--genetic_code GENETIC_CODE] [--infq INFQ] [--fq_position {flat,sample_dir}]
                      [--suffix_fq SUFFIX_FQ] [--scheduler_backend {auto,snakemake,threads}]
                      [--snakemake_bin SNAKEMAKE_BIN] [--smk_cores SMK_CORES] [--smk_mem_mb SMK_MEM_MB]
                      [--smk_latency_wait SMK_LATENCY_WAIT] [--smk_extra FLAG]
                      [--engine {native,mitoz,mitos2}] [--mitoz_bin MITOZ_BIN]
                      [--topology {auto,label,circular,linear}] [--preserve_phase1_metadata {yes,no}]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTA file/dir or GenBank file/dir.
  -o OUT, --out OUT     Output folder (default: mitom_out).
  --if_split_fasta {yes,no}
  --id_maxlen ID_MAXLEN
                        Max length of internal sample IDs (default: 15).
  --which_gene_first WHICH_GENE_FIRST
                        Gene used as sequence start (default: cox1).
  --fallback_genes FALLBACK_GENES
                        Genes tried in order when --which_gene_first is absent (default: cox2,cox3,cob,nad1;
                        'none' = leave such records unrotated). Records oriented this way are listed in
                        oriented_by_fallback.txt.
  --force_forward       RC the sequence when target gene is on minus strand (default: on).
  --no_force_forward
  --debug               Debug mode: DEBUG lines, program versions and the environment in the log; the output
                        of every program in the logs; intermediate files kept (--keep_intermediate yes).
                        `mitom debug OUT` reads the log.
  --topology {auto,label,circular,linear}
                        Circular or linear. auto (default): as labelled - GenBank LOCUS line; FASTA header
                        topology=circular|linear; header words ('circular', 'complete genome' / 'partial
                        genome'); an assembler's terminal overlap (removed): circular; else a sequence of >=
                        10 kb circular (one or several in a file), a shorter one linear - and the protein
                        genes at the two ends checked against their reference proteins: a gene running from
                        the end across to the start makes the record circular, a gene cut at an end of a
                        circular-labelled record makes it linear, and so do jellyfish / hydroids /
                        calcareous sponges (linear mtDNA) without a gene across the ends; an unlabelled
                        record no gene at the ends decides is circular with >= 12 of the 13 protein genes
                        and both rRNAs, else linear; changed records are annotated again
                        (reports/topology_decision.tsv). label: the label only. circular / linear: every
                        record.
  --preserve_phase1_metadata {yes,no}

Annotation:
  --mitos_bin MITOS_BIN
                        MITOS2 environment bin (runmitos). Default: found automatically (mito_env.py).
  --refdir DIR          Folder that holds the MITOS2 reference data refseq89m. Default: found automatically.
  --plot {yes,no}       Circular genome map per record (mito_plot.py, matplotlib) -> plots/<record>.png and
                        .pdf: genes by strand and type, GC content; with --infq also the read depth
                        (minimap2 when installed, else bwa mem). Default yes (about 2 s per record, --jobs
                        in parallel).
  --plot_style {circos,mitom}
                        Map style: circos (default, after mitoPlotR) = a framed lane per strand, names along
                        the circle, a position axis, GC content, read depth; mitom = genes on both sides of
                        one circle, GC content and GC skew rings.
  --ncr {yes,no}        yes (default): every stretch of >= --ncr_min_len nt no feature covers becomes a
                        misc_feature /note="non-coding region NCR<k>, <length> bp: no gene annotated
                        (MitoM)" in the GenBank and NCBI files, shaded and labelled on the maps (unnamed
                        control region, sequence that may hold unidentified genes); no: none (mito_ncr.py)
  --ncr_min_len NCR_MIN_LEN
                        shortest non-coding region written (nt, default 100; intergenic spacers are shorter)
  --curate {yes,no}     yes (default): mito_curate.py boundary rules and NCBI gene names -> 05.curated/, and
                        NCBI submission files (.fsa/.tbl/.sqn/.val) -> 06.ncbi_submission/
  --sbt SBT             NCBI submitter template (template.sbt from submit.ncbi.nlm.nih.gov/genbank/template)
                        for the .sqn files
  --rescue {auto,yes,no}
                        Add protein genes / rRNAs the core missed from a MITOS2 annotation (mito_rescue.py;
                        nothing annotated is changed; needs the optional MITOS2 installation). auto
                        (default) / no = not run (the core retries other loci itself); yes = run for records
                        missing a gene.
  --rrna_hmm {auto,yes,no}
                        Last-resort 12S/16S search for rRNAs still missing after the cores and the rescue,
                        with per-phylum profile HMMs built from NCBI rRNAs (mito_rrna_hmm.py; only missing
                        rRNAs are searched). auto = run when the library and nhmmer are found (default).
  --rrna_hmm_lib RRNA_HMM_LIB
                        Per-phylum rRNA HMM library (mito_rrna_hmm.py build). Default: models/rrna_hmm.
  --rrna_method {hmm,cm}
                        12S/16S search of the core: hmm (default, faster) = per-phylum profile HMMs,
                        MITOS2's covariance models where no model or hit; cm = MITOS2's covariance models
                        only. Same accuracy (2026-10-28, 450 leave-out records: found 98.5 / 98.5 %, both
                        ends exact 47.4 / 47.6 %); hmm takes seconds and ~0.1 GB, cm ~10 s and ~0.5 GB per
                        gene.
  --phylum PHYLUM       Phylum of all input records for --rrna_hmm (default: the phylum found by --clade
                        auto, or the phylum of an explicit --clade; with --clade Animal and FASTA input give
                        it here, otherwise the record's own taxonomy is used when it has one).
  --intron_fix {auto,yes,no}
                        Add introns with mito_intron.py (MFannot method: Exonerate protein2genome + MFannot
                        boundary rules + Rfam intron type), then with mito_intron_ref.py (intron sites of
                        annotated relatives, models/intron_sites). auto = run when exonerate is found
                        (default).
  --intron_bin INTRON_BIN
                        Directory containing exonerate, cmsearch, blastp, tblastn and blastn (e.g. a conda
                        env bin).
  --intron_lib INTRON_LIB
                        Protein library for intron models in MFannot header format ('>Species mt cox1;
                        ...'), repeatable. Default: MFannot's curated collections bundled in
                        models/mfannot_proteins.
  --intron_ref_gb INTRON_REF_GB
                        GenBank file/dir of close relatives used as extra homologs for intron models
                        (repeatable).
  --keep_intermediate {yes,no}
                        Keep the intermediate folders (00.*, 01.*, 03.*, 04.*, 05.*, _snakemake_*; run:
                        moved to work/). Default no: deleted after the final output is written; core logs of
                        problem samples are always kept (reports/failed_logs/).
  --threads THREADS     Threads per annotation job (default: min(8, number of CPUs) = 8).
  --jobs N              Records annotated at the same time (default: CPUs // threads = 9).
  --clade CLADE         Taxonomic group (sets the genetic code and group rules). auto (default): per record
                        - the taxonomy and CDS transl_table of a GenBank record; FASTA and records without a
                        taxonomy: blastx to the reference proteins (mito_classify.py: phylum of the closest
                        genomes, their genetic code checked on the codons; no close relative: Animal and the
                        code of the codons). Several phyla / codes: `run` annotates each group with its own
                        clade and code and merges the results (reports/clade_decision.tsv); `import` and
                        `annotate` need one. Or give one of: Chordata, Arthropoda, Echinodermata, Annelida-
                        segmented-worms, Bryozoa, Mollusca, Nematoda, Nemertea-ribbon-worms, Porifera-
                        sponges, Cnidaria, Platyhelminthes-flatworms, Chaetognatha-arrow-worms, Animal.
  --genetic_code GENETIC_CODE
                        Genetic code (default with --clade auto: the CDS transl_table of each record, else
                        the code of its lineage or of the closest reference genomes; otherwise the code of
                        the clade - 2 Chordata, 9 Echinodermata, 4 Porifera, else 5; required for Cnidaria,
                        Platyhelminthes, Chaetognatha, Animal).

FASTQ (optional - read depth ring of the genome maps, --plot):
  --infq INFQ           Directory with the paired FASTQ files of the samples (matched by sample name).
  --fq_position {flat,sample_dir}
  --suffix_fq SUFFIX_FQ

Scheduler:
  --scheduler_backend {auto,snakemake,threads}
                        Dispatch backend (default: auto). auto — use Snakemake if available, else threads
                        snakemake — always use Snakemake threads — always use ThreadPoolExecutor Annotation
                        is parallelised; finalization runs in the main process regardless of backend.
  --snakemake_bin SNAKEMAKE_BIN
  --smk_cores SMK_CORES
  --smk_mem_mb SMK_MEM_MB
  --smk_latency_wait SMK_LATENCY_WAIT
  --smk_extra FLAG      Extra Snakemake flags appended verbatim. Example: --smk_extra '--profile slurm'

Annotation engine:
  --engine {native,mitoz,mitos2}
                        Program that finds the genes: native (default) = MitoM's own core (mito_core.py);
                        mitoz = MitoZ 3.6, mitos2 = MITOS2 2.1.10 (not installed with MitoM - install them
                        yourself, see mito_engines.py; found on PATH or in a conda environment). Every later
                        MitoM step runs on their annotation: reorientation, round 2, curation, NCBI files,
                        genome maps.
  --mitoz_bin MITOZ_BIN
                        Folder with the mitoz program (--engine mitoz). Default: found automatically.
```

## mitom reorient

```text
usage: mitom reorient [-h] -i INPUT [-o OUT] [--if_split_fasta {yes,no}] [--id_maxlen ID_MAXLEN]
                      [--which_gene_first WHICH_GENE_FIRST] [--fallback_genes FALLBACK_GENES]
                      [--force_forward] [--no_force_forward] [--debug]

Turn circular records to the strand of the target gene (default cox1) and start them at it; fallback genes
when it is missing; linear records are not rotated.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTA file/dir or GenBank file/dir.
  -o OUT, --out OUT     Output folder (default: mitom_out).
  --if_split_fasta {yes,no}
  --id_maxlen ID_MAXLEN
                        Max length of internal sample IDs (default: 15).
  --which_gene_first WHICH_GENE_FIRST
                        Gene used as sequence start (default: cox1).
  --fallback_genes FALLBACK_GENES
                        Genes tried in order when --which_gene_first is absent (default: cox2,cox3,cob,nad1;
                        'none' = leave such records unrotated). Records oriented this way are listed in
                        oriented_by_fallback.txt.
  --force_forward       RC the sequence when target gene is on minus strand (default: on).
  --no_force_forward
  --debug               Debug mode: DEBUG lines, program versions and the environment in the log; the output
                        of every program in the logs; intermediate files kept (--keep_intermediate yes).
                        `mitom debug OUT` reads the log.
```

## mitom transfer_metadata

```text
usage: mitom transfer_metadata [-h] -m METADATA -n NEW_BODY [-o OUT] [--if_replace_metadata {yes,no}]
                               [--overwrite_new_body {yes,no}] [--output_name OUTPUT_NAME]
                               [--smart_metadata {yes,no}] [--auto_accept_metadata {yes,no}]

options:
  -h, --help            show this help message and exit
  -m METADATA, --metadata METADATA
                        Original GenBank file/dir (metadata source).
  -n NEW_BODY, --new_body NEW_BODY
                        Re-annotated GenBank file/dir (body source).
  -o OUT, --out OUT
  --if_replace_metadata {yes,no}
  --overwrite_new_body {yes,no}
  --output_name OUTPUT_NAME
                        Output filename (single-file mode only).
  --smart_metadata {yes,no}
  --auto_accept_metadata {yes,no}
```

## mitom replace_organism

```text
usage: mitom replace_organism [-h] -m METADATA -n NEW_BODY [-o OUT] [--overwrite_new_body {yes,no}]

options:
  -h, --help            show this help message and exit
  -m METADATA, --metadata METADATA
                        Original GenBank file/dir (organism/source name source).
  -n NEW_BODY, --new_body NEW_BODY
                        Re-annotated GenBank file/dir whose /organism= needs fixing.
  -o OUT, --out OUT     Output folder (default: replace_organism_out).
  --overwrite_new_body {yes,no}
                        Overwrite the input files in-place (default: no).
```

## mitom dedup

```text
usage: mitom dedup [-h] -i INPUT [-o OUT]

Remove duplicate genomes, as `run` does before annotating. The same genome is often in NCBI more than once:
the authors' record (e.g. MK123456) and the RefSeq copy NCBI made of it (NC_...), two versions of one
accession (MK123456.1 / .2), or another accession with an identical sequence. One record per genome is kept:
the NC_ record, the newest version. dedup_report.tsv: every input record, kept or dropped, why, and the
record kept instead.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        GenBank and / or FASTA file or folder. FASTA records: identical sequences, and
                        versions of one accession when the header's first word is an accession (MK123456.1);
                        other names are never merged.
  -o OUT, --out OUT     Output folder (default: dedup_out).
```

## mitom gb2fasta

```text
usage: mitom gb2fasta [-h] -i INPUT [-o OUT] [--split]

Convert GenBank files (NCBI .gb, MitoZ .gbf ...) to FASTA: one <name>.fasta per input file, or one
<accession>.fasta per record with --split. Headers: >ACCESSION.V definition topology=circular|linear - the
topology tag keeps a circular genome circular when `mitom run` reads several records from one FASTA file.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        GenBank file or folder (subfolders included).
  -o OUT, --out OUT     Output folder (default: gb2fasta_out).
  --split               One FASTA file per record (<accession>.fasta).
```

## mitom setup

```text
usage: mitom setup [-h] [--refdir REFDIR]

Download the MITOS2 reference data refseq89m (Zenodo 4284483, md5 checked) into refdir (default
$MITOS_REFDIR or ~/mitos_refs) unless it is found, then save the found locations (doctor --write-config).

options:
  -h, --help       show this help message and exit
  --refdir REFDIR  folder for refseq89m (default: $MITOS_REFDIR or ~/mitos_refs)
```

## mitom doctor

```text
usage: mitom doctor [--write-config]

Which external programs MitoM finds (config, PATH, conda environments): OK / -- (optional, missing) /
!! (required, missing), with what each one is used for.

options:
  --write-config  save the locations found to mito_config.json (mitom setup does it too)
```

## mitom curate

```text
usage: mito_curate.py [-h] -i INPUT -o OUT [--genetic_code GENETIC_CODE] [--rrna_bounds {trim,fill}]
                      [--max_fill MAX_FILL] [--max_start_shift MAX_START_SHIFT] [--aggressive_start]
                      [--no_upstream_atg] [--abbreviated_overlap] [--abbreviated_pcg_overlap]
                      [--aggressive_start_genes AGGRESSIVE_START_GENES] [--start_model START_MODEL]
                      [--dump_starts DUMP_STARTS] [--phylum PHYLUM] [--table2asn TABLE2ASN] [--no_validate]

mito_curate.py - rule-based correction of mitochondrial genome annotation (PCGs, 12S/16S rRNA).

Applies the annotation criteria of Fourdrilis et al. 2018 (Sci. Rep. 8:17964) as summarised by
Ghiselli et al. 2021 (Phil. Trans. R. Soc. B 376:20200159, doi:10.1098/rstb.2020.0159, section 3)
to automatic annotation (MitoM's core, or another annotator's GenBank such as MITOS output):

  P1  PCG start: the first eligible in-frame start codon at the 5' end, i.e. the start codon
      nearest the preceding gene WITHOUT overlapping it (criterion i). Only applied when the
      annotated start overlaps the preceding same-strand feature (move downstream), or when an
      in-frame start lies upstream before the preceding feature and the open frame reaches it
      (move upstream, at most --max_start_shift codons; only when the annotated start codon is
      not a valid start, or for --aggressive_start_genes; annelids: nad5 and nad4L; or, unless
      --no_upstream_atg, to the most upstream in-frame ATG when the start codon is not ATG). A start that
      overlaps the preceding gene is moved downstream only when a start codon lies within 15 nt. A start placed
      upstream of a whole same-strand tRNA (in the gene's first third) moves to the first in-frame start codon after
      it (table2asn BadCDScomponentOverlapTRNA: unionid nad1 across trnL2).
  P11 PCG start, after P1: of the in-frame start codons within 30 codons of P1's start (reading frame open in
      between), the one NCBI records would choose, as scored by gradient-boosted trees learned from 1,683 NCBI
      records (mito_start.py, models/start_model/start_model.json; --start_model no turns it off): codon, distance to
      the preceding feature, the core's alignment (residues missing at the N-terminus), gene, phylum. Cross-validated
      by family: exact starts 89.1 -> 90.7 % (arthropods 77.8 -> 82.3 %, cnidarians 79.8 -> 83.2 %).
  P2  Two PCGs on the same strand and in the same frame cannot overlap (criterion ii) - reported.
  D1  A protein gene annotated twice at different places (e.g. MitoZ on a duplicated region): only one
      copy is kept - the complete one (start + stop, no internal stop), the longer, then the first in
      the record; the other CDS and its gene feature are removed (NCBI records annotate one copy), unless
      the other copy is complete too and its sequence is >= 95 % identical with the same length +-3 nt, for
      >= 2 genes (a duplicated block, annotated twice by NCBI: NC_009093.1 mite, 9 genes; one repeated
      atp8 is annotated once: NC_049135.1).
      An overlapping copy on the other strand that is incomplete or has internal stops is removed as well.
  P3  PCG end: the first in-frame complete stop codon, or an abbreviated stop (T / TA) ending
      immediately before the downstream same-strand tRNA, completed by polyadenylation
      (criterion iii; tRNA punctuation model, Ojala et al. 1981). /transl_except + note added.
      When no tRNA follows and the first complete stop lies > 30 nt inside the next same-strand
      protein gene (not the bicistronic pairs atp8-atp6, nad4L-nad4; rRNA ends are trimmed by R1
      instead; > 10 nt for any other feature; an rRNA start does not bound an end without any stop signal, R1 moves
      the rRNA), the gene is not extended into it: the annotated end is kept when it is an abbreviated T / TA stop (the end of
      the protein homology, e.g. cox1 before cox2), else the T / TA right before the next gene
      (e.g. nad6 before cob) is used. An end that stops short of the first complete stop with
      nothing annotated in between (core homology lost at the C-terminus, 100-500 nt in divergent
      phyla) is extended to that stop (at most 600 nt). A gene that would contain a whole same-strand tRNA and has
      no T / TA before it ends right before the tRNA (3'-partial).
  P7  start and end in different frames (internal stops): replaced by the longer ORF - annotated start
      to its first stop, or annotated stop with the in-frame start nearest the annotated start - when
      it keeps >= half the length (metazoa benchmark 2026-10: 14 CDS with internal stops, NCBI 0).
  P8  start and end in one frame with internal stops (a start placed upstream across a stop, an end past the gene's
      stop): the longer of start .. first stop and first start codon after the last internal stop .. end (>= half).
  P9  one untranslated extra nucleotide (bird / turtle nad3 at position 174, Mindell et al. 1998; rarely other genes):
      join(start..p-1,p+1..end) as NCBI writes it, with its note, when no single-frame repair keeps 85 % of the gene
      and a start / stop codon within 30 nt of the annotated ends make an open frame over one skipped nucleotide.
  P10 a CDS that still has internal stops is written /pseudo (gene and CDS, not translated; table2asn InternalStop).
  O2  a join across the origin of a circular record written by the core with <= 300 nt past the origin is cut back to
      the sequence end; P1 / P3 / O1 place its end again (annelid nad2 before cox1 ended without a stop: NoStop).
  T1  a weak tRNA call (MiTFi E > 1e-3) lying wholly inside a same-strand protein gene is removed when it is a second
      copy or lies in the gene's middle third; in its first / last third P1 / P3 move the gene end instead.
  C1  a tRNA still lying wholly inside a same-strand protein gene: the gene is removed when it has no start or stop
      codon or is /pseudo (a false call), else the conflict is reported.
  O1  circular record: a CDS that ends without a stop <= 2 nt before the transcript end, or further when its
      open frame runs to the end with nothing annotated in between (the core ends at the last aligned codon), is
      continued across the origin to its first in-frame stop (join; genes cut by the cox1-first rotation); at most
      60 nt where another feature starts at the origin (genuine overlaps are <= ~50 nt), else 300 nt.
  P5/P6  illegal start / no stop codon: CDS marked 5'- / 3'-partial with 'start / stop codon not
      determined' (NCBI practice; table2asn SEQ_FEAT.StartCodon / NoStop otherwise). A stop-less CDS ending 1-2 nt
      before the sequence end is first extended to it (the cut codon). A single-interval CDS the core already marked
      partial (< / >) goes through P1 / P3 / O1 like any other; only joins are left alone (P4).
  R1  rRNA (12S rrnS, 16S rrnL) must not overlap flanking genes; boundaries are trimmed to the
      neighbours (criterion vi) and, by default, extended to the flanking genes (gaps <= --max_fill
      nt, the GenBank convention; --rrna_bounds fill). --rrna_bounds trim only removes overlaps.
  R3  rRNA names standardised (rrnS / rrnL; 12S / 16S ribosomal RNA); rRNA with a large CDS
      overlap is reported.
  N1  NCBI gene names and products (COX1 ... ND4L, rrnS / rrnL, trnX with L1 / L2, S1 / S2 by anticodon); the
      notes the annotation core writes for itself (reference protein and score, model columns, MitoM-ML / MiTFi
      scores) are removed, so the GenBank and the NCBI submission files carry only biological notes.
tRNA boundaries are kept as annotated (criterion v). Every change is listed in
curation_report.tsv with the rule, old and new coordinates and the reason.

Each run validates every record with NCBI table2asn before and after curation (when table2asn is
available) and writes validation_compare.tsv: every validator message type per feature marked
fixed / remaining / new, so the effect of the curation is measured, not assumed.

Usage: python3 mito_curate.py -i gbf_file_or_dir -o out_dir --genetic_code 5 [--rrna_bounds trim|fill]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
  -o OUT, --out OUT
  --genetic_code GENETIC_CODE
  --rrna_bounds {trim,fill}
                        fill (default): remove overlaps with flanking genes and extend rRNA ends to them
                        (gaps <= --max_fill), the GenBank convention; trim: only remove overlaps (Fourdrilis
                        criterion vi). 2026-10-21 dev benchmark (90 rRNAs vs NCBI): both ends exact trim 22,
                        fill(300) 33
  --max_fill MAX_FILL   Largest gap filled with --rrna_bounds fill (nt; default 300: exact rRNAs 30 / 33 /
                        33 / 31 at 100 / 300 / 500 / 1000 on the dev set)
  --max_start_shift MAX_START_SHIFT
                        Largest start/stop move (codons)
  --aggressive_start    Always move PCG starts to the in-frame start codon nearest the preceding gene
                        (criterion i literally); default only when the annotated start codon is invalid
  --no_upstream_atg     P1: keep a non-ATG start codon instead of moving to the most upstream in-frame ATG
                        (mitom run sets it for nematodes and arthropods)
  --abbreviated_overlap
                        P3: a complete stop overlapping the next same-strand tRNA by 1-2 nt becomes the T /
                        TA before it (vertebrate RefSeq convention; before a protein gene the complete stop
                        is kept: atp6 TAA / cox3 ATG; mitom run sets it for the Chordata clade)
  --abbreviated_pcg_overlap
                        P3: a complete stop more than 10 nt inside the next same-strand protein gene becomes
                        the T / TA before it, as O1 at the origin (NCBI annelid records; others allow 30 nt;
                        mitom run sets it for the Annelida clade)
  --aggressive_start_genes AGGRESSIVE_START_GENES
                        Genes whose start is always moved to the in-frame start codon nearest the preceding
                        gene (default none; `mitom run` passes ND5,ND4L for the Annelida clade). Benchmark
                        2026-10: right for nad5 / nad4L in every differing annelid record (40), wrong in
                        most other phyla (173 metazoa records: 17 right, 52 wrong).
  --start_model START_MODEL
                        P11: auto (default) = models/start_model/start_model.json when present; no = off; or
                        a model file (benchmark/train_start_model.py)
  --dump_starts DUMP_STARTS
                        write every P11 start candidate with its features to this TSV (training data)
  --phylum PHYLUM       phylum of records without a taxonomy (a feature of the start model)
  --table2asn TABLE2ASN
                        NCBI table2asn for validation before/after curation (default: found by mito_env.py)
  --no_validate         Skip the NCBI validation before/after curation
```

## mitom intron

```text
usage: mitom intron [-h] -i INPUT -o OUT [--method {both,denovo,ref}] [--genetic_code GENETIC_CODE]
                    [--threads THREADS] [--jobs JOBS] [--bin BIN] [--lineage LINEAGE]

Introns as `run` adds them: mito_intron.py (de novo, MFannot method) -> OUT, then mito_intron_ref.py (intron
sites of annotated relatives) on OUT in place. Reports: OUT/intron_report.tsv, OUT/intron_ref_report.tsv.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        GenBank file or folder.
  -o OUT, --out OUT     Output folder (may be the input folder: updated in place).
  --method {both,denovo,ref}
                        both (default): de novo, then the site transfer; denovo: mito_intron.py only (the
                        `mitom intron` of MitoM <= 1.0.0); ref: mito_intron_ref.py only.
  --genetic_code GENETIC_CODE
                        Default 5.
  --threads THREADS     Threads for blastp (de novo). Default 1.
  --jobs JOBS           GenBank files processed in parallel. Default 1.
  --bin BIN             Folder with exonerate, blastp, tblastn, blastn, cmsearch. Default: found
                        automatically.
  --lineage LINEAGE     NCBI lineage of records without their own taxonomy, for the site transfer (e.g.
                        'Eukaryota;Metazoa;Cnidaria').

Other options go to mito_intron.py (python mito_intron.py -h: --lib, --ref_gb, --maxintron, --timeout ...),
with --method ref to mito_intron_ref.py (python mito_intron_ref.py annotate -h).
```

## mitom trna

```text
usage: mito_trna.py [-h] -i INPUT -o OUT [--genetic_code GENETIC_CODE] [--tools TOOLS] [--priority PRIORITY]
                    [--select {auto,mitom,mitos2_or_support,support,mitos2}] [--min_support MIN_SUPPORT]
                    [--bin BIN] [--mitos_bin MITOS_BIN] [--refdir DIR] [--mitos_ref MITOS_REF]
                    [--lineage LINEAGE] [--threads THREADS] [--mitos_mode {fast,sensitive,full}]

mito_trna.py - multi-tool tRNA annotation of mitochondrial genomes, GeSeq style.

Several tRNA predictors are run on each record and ALL their predictions are written into one
GenBank file for review (e.g. in Geneious), as GeSeq does with its third-party tRNA tools
(Tillich et al. 2017, Nucleic Acids Res 45:W6, doi:10.1093/nar/gkx391). A second GenBank keeps one
tRNA per locus for submission. The number of tRNAs is never imposed: it varies widely in animal
mtDNA (e.g. 1-2 in most anthozoans, 23-24 with duplications in some annelids), so every
prediction is judged on its own evidence.

Sources (any subset, --tools):
  mitos2    MITOS2 runmitos, tRNA module only (metazoan reference data, e.g. refseq89m)
  arwen     ARWEN (metazoan mitochondrial tRNA; recommended by GeSeq for metazoan mtDNA)
  aragorn   ARAGORN (-mt)
  trnascan  tRNAscan-SE 2 (-O organellar mode)
  input     tRNA features already present in the input GenBank
  mitom     MitoM-ML (trna_ml.py): MITOS2 tRNA CMs + clade CMs on the record's NCBI lineage, a learned caller
            using ARWEN / ARAGORN as extra evidence, learned ends. Unseen-family benchmark (1,069 RefSeq invertebrate
            mitogenomes, 22,099 tRNAs): recall 98.5 % at 0.17 extra tRNAs per genome vs 94.9 % / 0.19 for the
            MiTFi CMs alone; 0 calls on 150 shuffled genomes (projects/invert_trna/REPORT.md). Lineage: GenBank
            taxonomy of the record, else --lineage (without one only the pan-metazoan CMs are used)

Support is counted in independent method families (MITOS2; ARWEN+ARAGORN as one, same author and
algorithm; tRNAscan-SE; input annotation), not in tools.

Locus groups: predictions on the same strand overlapping >= 50 % of the shorter one, with the
same amino acid AND anticodon (so trnL1/trnL2 and trnS1/trnS2 stay separate). Overlapping
predictions with different anticodons form separate groups flagged as a conflict.

Selection --select mitom (default when MitoM-ML is installed: models/trna_ml, cmsearch, cmalign, refseq89m): the
selected tRNAs are the MitoM-ML calls (its boundaries and anticodons); the other sources stay in the review GenBank.
Selection --select mitos2_or_support (default otherwise): a group is selected when MITOS2 found it, or when
>= --min_support (2) independent method families agree. Overlapping groups with different anticodons
(conflict): the one with stronger evidence (MITOS2, then more families) is kept, ties drop both.
--select mitos2: the input tRNAs are replaced by the MITOS2 tRNAs (tools mitos2,input;
the input tRNAs stay visible in the review GenBank). Boundaries otherwise come from the first source
in --priority (MITOS2 first: most exact boundaries vs NCBI).

MITOS2 runs in fast mode by default (scripts/runmitos_fast.py): runmitos re-searches every tRNA model it did not
find with a sensitive cmsearch (--max, E 0.1, single CPU): 640 s vs 8 s on one anthozoan (Metridium), plus spurious
tRNAs; fast mode skips that pass, --mitos_mode sensitive runs it in parallel (fast mode missed 28/28 of the weak tRNAs
found only by that pass, 2026-10-20 dev benchmark), --mitos_mode full is the unmodified runmitos. `mitom run` does
not run this script: its tRNAs come from the annotation core (mito_core.py, MitoM-ML with MiTFi).

The annotation core's tRNA step is here too (core_trnas, called by mito_core.py for each record): MitoM-ML calls
kept when sure enough or confirmed by MiTFi, MiTFi's ends for vertebrates, MiTFi's calls of the families MitoM-ML does
not call; MiTFi alone (cmsearch with the MiTFi models of refseq89m, a sensitive pass outside the protein gene / rRNA
interiors, mitfi_import.jar) when MitoM-ML is not installed.

Outputs per input file: review/<name>.gbf (all predictions; /note = source, score, anticodon,
group, support), selected/<name>.gbf (selected tRNAs only; input for mito_curate.py -i out/selected),
and trna_report.tsv (one row per group).

Usage:
  python3 mito_trna.py -i gb_or_fasta -o out_dir --genetic_code 5       --mitos_refdir DIR --mitos_ref refseq89m --mitos_bin ENV/bin --bin ENV/bin

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        GenBank or FASTA file / directory
  -o OUT, --out OUT
  --genetic_code GENETIC_CODE
  --tools TOOLS         default: mitom,arwen,aragorn,input with --select mitom;
                        mitos2,arwen,aragorn,trnascan,input with mitos2_or_support / support; mitos2,input
                        with --select mitos2
  --priority PRIORITY   source whose boundaries are used for the selected tRNA of a group
  --select {auto,mitom,mitos2_or_support,support,mitos2}
                        selected tRNAs: auto (default) = mitom when MitoM-ML is installed, else
                        mitos2_or_support; mitom = the MitoM-ML calls (a record whose run fails keeps its
                        input tRNAs); mitos2 = replace the input tRNAs with the MITOS2 tRNAs (MITOS2
                        boundaries and anticodons; a record whose MITOS2 run fails is kept unchanged);
                        mitos2_or_support (default) = found by MITOS2, or supported by >= --min_support
                        independent method families; support = families only. Benchmark (39 annelid + 71
                        anthozoan/sponge/placozoan NCBI mitogenomes): MITOS2 alone had the fewest false
                        positives; ARWEN/ARAGORN alone gave hundreds in tRNA-poor genomes
  --min_support MIN_SUPPORT
                        independent method families needed for a group to be selected (ARWEN+ARAGORN count
                        once)
  --bin BIN             directory with arwen, aragorn, tRNAscan-SE
  --mitos_bin MITOS_BIN
                        directory with runmitos (MITOS2 environment bin)
  --refdir DIR          folder that holds the MITOS2 reference data refseq89m
  --mitos_ref MITOS_REF
  --lineage LINEAGE     NCBI lineage (';'-separated) for MitoM-ML clade models when the input has no GenBank
                        taxonomy
  --threads THREADS     CPUs for MitoM-ML cmsearch / cmalign
  --mitos_mode {fast,sensitive,full}
                        fast (default): MITOS2 without its sensitive-mode fallback
                        (scripts/runmitos_fast.py); full: unmodified runmitos (much slower when tRNAs are
                        missing, e.g. anthozoans)
```

## mitom rescue

```text
usage: mito_rescue.py [-h] -i INPUT --genetic_code GENETIC_CODE [--mitos_bin MITOS_BIN] [--refdir DIR]

mito_rescue.py - add genes the annotation core missed, taken from a MITOS2 annotation of the same sequence
(optional: --rescue yes; needs the MITOS2 installation).

The primary annotation can miss short or divergent genes that MITOS2 finds. 2026-10-20 dev benchmark, with MitoZ as
the primary annotator then (NCBI NC_ annotation as reference, offline simulation):
  MitoZ-supported phyla (45 records): protein genes found 579 -> 580 (atp8), rRNA 81 -> 87, no new FP.
  Other phyla (48 records, MitoZ with a borrowed protein database): protein genes 525 -> 548, rRNA 75 -> 82.
MITOS2 boundaries were less exact (52 % vs MitoZ's 80 % both ends exact), so nothing annotated is changed: only
genes absent from the record are added (each listed in rescue_report.tsv with a hint to check it). `mitom run`
runs this step only with --rescue yes (MITOS2 is not installed with MitoM).

Rules (each from the benchmark and from the 2026-10-21 review):
  - protein gene: added when no CDS of that gene exists and the MITOS2 gene does not lie mostly
    (> 50 %) on an existing CDS, rRNA or tRNA, or on a gene added earlier in the same call (longest first;
    that is a naming disagreement, not a missed gene).
    Containment release: an existing CDS/rRNA stops blocking only when it looks like an intron-spanning or
    merged call: it is longer than 1.5 x MAX_LEN of its own gene (MAX_LEN = longest expected length in nt)
    AND wholly contains >= 2 MITOS2 genes of other names that are absent from the record (NC_006902.1:
    MitoZ nad5 over its 11 kb group I intron, 10 genes inside). A normal-length gene never releases, so
    homologue fragments (nad4 / nad2 inside nad5) are not added as bogus CDS.
  - rRNA: MITOS2 often splits a divergent rRNA into short pieces (29-200 nt); the pieces (-a/-b) are not
    joined, only the longest piece of each rRNA, and only when >= 300 nt, is added (rescued-rRNA FP 13 -> 3
    on the other phyla). CDS fragments (-a/-b) are joined only on the same strand with a gap of 0-3000 nt.
  - /translation of a rescued CDS: trimmed to whole codons, first codon M if it is a start codon of the
    table, trailing stop removed; omitted when an internal stop remains.
  - tRNAs are not handled here (the annotation core's tRNA step).
  - The feature /note is neutral (NOTE); the advice to check the boundaries is the "hint" column of
    rescue_report.tsv.

MITOS2: Donath et al. 2019, Nucleic Acids Res 47:10543 (doi:10.1093/nar/gkz833). It runs through
scripts/runmitos_fast.py (Infernal searches in parallel, same results as runmitos).

Usage: python3 mito_rescue.py -i 04.final_gb/ --genetic_code 5 [--mitos_bin DIR] [--mitos_refdir DIR]
Writes the GenBank files in place and rescue_report.tsv next to them.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        directory of final GenBank files (*.gbf), edited in place
  --genetic_code GENETIC_CODE
  --mitos_bin MITOS_BIN
  --refdir DIR          folder that holds the MITOS2 reference data refseq89m
```

## mitom rrna-hmm

```text
usage: mito_rrna_hmm.py [-h] {build,search} ...

mito_rrna_hmm.py - 12S/16S with per-phylum profile HMMs (MitoM): the annotation core's rRNA step and a last-resort
search.

The core's rRNA step (annotate_rrnas, called by mito_core.py for each record): the per-phylum profile HMM of the
library (nhmmer) when it has the record's phylum, else MITOS2's covariance model (cmsearch); hits joined by the rules
below (best_spans), a span over less than MIN_COVER of the model not taken.

Most genomes are annotated well by the core (per-phylum HMMs / MITOS2 covariance models); this step runs only for
an rRNA that is still missing after them and the MITOS2 rescue. It searches the genome with an nhmmer
profile HMM built from the NCBI rRNAs of the record's own phylum, as EMMA does with its (vertebrate)
HMMs (knowledge/tools/emma.md).
Evidence (2026-10-21 dev benchmark, 186 NCBI rRNAs of 93 records, leave-genus-out, knowledge/topics/
rrna_annotation_methods.md section 6): at search level MITOS2 CM (fragments merged) found 166, CM then
per-phylum HMM 175 (knowledge/experiments/rrna_ref/rr_detect_summary_dev.txt). That experiment merged
all same-strand hits within 300 nt and applied none of the filters below, so it is not a measurement of
this code. In the full pipeline (HMM step before mito_curate, library without any benchmark genus;
pipeline_ab_dev.txt, run with an earlier version of this step) missed rRNAs went 17 -> 14, false
positives 9 -> 10, exactly matching ones 59 -> 59 (present but partial rRNAs count as present). The
coverage / partial-end / window rules below were added after that run and are covered by unit tests only.

Rules:
  - only rRNAs absent from the record are searched (names as in mito_names.RRNA_ALIAS and free text);
    the phylum is taken from the record's taxonomy (or --phylum); phyla with fewer than MIN_REFS
    reference rRNAs have no model and are skipped;
  - hits of one model are joined only when they are pieces of one gene: same strand, disjoint and
    consecutive model columns in the same order as in the genome, <= MAX_GAP nt apart. Hits over the
    same model columns are copies (mitochondrial genes can be duplicated) and are never joined; all such
    clusters are built first and the one with the highest summed score wins. A second copy is only
    reported ("second copy at X..Y not added"), never added;
  - circular records are searched as seq + first 2000 nt, so a gene across the origin is one hit and
    becomes a join location; clusters that start in the appended part are duplicates and dropped;
  - the span must cover >= MIN_COVER (0.4) of the model columns (LENG of the .hmm), unless a linear record's
    sequence end cuts it on the side of its missing columns (the gene runs on beyond the sequenced part; the
    other end complete), and be >= MIN_RRNA
    (200) nt; <= LEN_HI (1.3) x the 95th percentile length of the references (manifest p95_len) and, when
    complete at both ends, >= LEN_LO (0.7) x their 5th percentile; it must not lie mostly (> 50 %) on an
    existing CDS / rRNA. A span whose model columns start after column 30 / end before LENG - 30 is added
    with < / > on the 5' / 3' end ("added (partial 5')"); report columns coverage and note;
  - other ends are not refined here (mito_curate.py R1); the feature /note is a plain statement of the
    method, the "check boundaries manually" hint is only in the report.

Usage:
  python3 scripts/fetch_rrna_refs.py -o refs --all-metazoa [--set benchmark/sets/mito_bench_2026-10.tsv]
  python3 mito_rrna_hmm.py build --gb refs/*.gb -o models/rrna_hmm [--exclude accessions.txt] [--manifest_only]
  python3 mito_rrna_hmm.py search -i 04.final_gb/ [--lib models/rrna_hmm] [--phylum Ctenophora]
build: one GenBank file per phylum (file name = phylum), made by scripts/fetch_rrna_refs.py (reusable form
of knowledge/experiments/rrna_ref/fetch_refs.py, which also writes refs_<phylum>.txt accession lists);
rRNAs of 300-4000 nt with <= 1 % N; --manifest_only recomputes manifest.tsv (lengths, percentiles) without
re-aligning. search writes the GenBank files in place (all *.gbf / *.gb below the directory; failures
become "failed: ..." rows) and rrna_hmm_report.tsv next to them.

HMMER: Eddy 2011, PLoS Comput Biol 7:e1002195 (doi:10.1371/journal.pcbi.1002195); nhmmer: Wheeler &
Eddy 2013, Bioinformatics 29:2487 (doi:10.1093/bioinformatics/btt403); MAFFT: Katoh & Standley 2013,
Mol Biol Evol 30:772 (doi:10.1093/molbev/mst010).

positional arguments:
  {build,search}
    build         build per-phylum rRNA HMMs from reference GenBank files
    search        search missing rRNAs in final GenBank files (in place)

options:
  -h, --help      show this help message and exit
```

## mitom validate

```text
usage: ncbi_validate.py [-h] -i INPUT -o OUT [--table2asn TABLE2ASN] [--focus FOCUS] [--only_ours]
                        [--nuclear]

ncbi_validate.py - check GenBank annotation against NCBI's own submission validator.

Each record is converted to the NCBI submission inputs (.fsa with source modifiers and
5-column feature table .tbl) and validated with table2asn -V vb, the program NCBI uses
when it processes a submission. Messages are summarised per severity and code.
Submitter/publication messages (no template) are ignored: only sequence and feature
problems (SEQ_FEAT, SEQ_INST, SEQ_DESCR on the record) count.

Usage:
  python3 ncbi_validate.py -i gb_file_or_dir -o out_dir [--table2asn path] [--focus intron]
Outputs: out_dir/<record>/{.fsa,.tbl,.val}, out_dir/validation.tsv (one row per message),
and a summary on stdout.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
  -o OUT, --out OUT
  --table2asn TABLE2ASN
                        table2asn executable (default: found by mito_env.py)
  --focus FOCUS         Only count messages whose text matches this regex (e.g. intron|CDS)
  --only_ours           Only count messages on features written by MitoM's intron steps (mito_intron.py,
                        mito_intron_ref.py)
  --nuclear             nuclear sequences (mitom nuc): no location, [gcode=]
```

## mitom gb2tbl

```text
usage: mitom gb2tbl [-h] -i INPUT [INPUT ...] -o OUT [--split {all,file,record}] [--name NAME]
                    [--source_modifiers {yes,no}] [--source_table SOURCE_TABLE] [--modifiers MODIFIERS]
                    [--product_table PRODUCT_TABLE] [--genetic_code GENETIC_CODE] [--validate] [--sbt SBT]
                    [--table2asn TABLE2ASN]

mito_tbl.py - GenBank / EMBL flat files -> NCBI submission files: 5-column feature table (.tbl) and FASTA with
source modifiers (.fsa), the pair BankIt and table2asn take (`mitom gb2tbl`).

Input: one or more GenBank (.gb .gbk .gbf .gbff .genbank) or EMBL / ENA (.embl .ebl .ena) files or folders (every such
file below a folder), any mix, or '-' for standard input (efetch -format gb | mitom gb2tbl -i - ...); the format is read
from the file (LOCUS / ID line). MitoM's final_annotation, MitoZ / MITOS2 / Geneious output, NCBI and ENA downloads.
Output (-o OUT):
  --split all (default)  OUT/<name>.tbl + OUT/<name>.fsa with every record (one batch submission); <name>: --name,
                         else the stem of a single input file, else 'submission'
  --split file           one pair per input file: OUT/<file stem>.tbl / .fsa
  --split record         one pair per record:     OUT/<SeqID>.tbl / .fsa
  gb2tbl_report.tsv      every record: its SeqID, source file, output pair, organism, length, topology, location,
                         genetic code, gene counts, warnings (and the NCBI ERROR / WARNING counts with --validate)
Records:
  - the same record in two inputs (same ID and sequence) is written once; different records with one ID are renamed
    ID_2, ID_3 ... (NCBI needs unique SeqIDs in a submission); an ID is cleaned to NCBI's characters
    (A-Z a-z 0-9 _ - . : * #) and cut to 50
  - the FASTA header carries [organism=] (the source /organism), [location=mitochondrion] / [location=chloroplast]
    (from /organelle; none for a nuclear record), [topology=], [completeness=complete] for a circular record,
    [mgcode=] / [gcode=] (the CDS transl_table), [taxid=] (/db_xref="taxon:..."), the source qualifiers NCBI takes as
    modifiers: isolate, strain, specimen_voucher, country / geo_loc_name, collection_date, lat_lon, collected_by,
    identified_by, isolation_source, haplotype, host, sex, tissue_type, dev_stage, note ... (--source_modifiers no: only
    the first ones), then the DEFINITION as the title (as NCBI's gbf2fsa; not a placeholder such as MitoZ's
    'topology=linear'); the sequence in lines of 70
  - the feature table: every feature but source, partial ends < >, joins, qualifiers (not /translation, /db_xref,
    /protein_id, MitoM's internal notes) - the table `mitom run` writes to ncbi_submission/
Your source information and products (the tables of GB2sequin, https://chlorobox.mpimp-golm.mpg.de/GenBank2Sequin.html):
  --source_table FILE    one row per record, tab-separated (NCBI's .src) or CSV: the first column the record (its
                         SeqID, ID with or without version, LOCUS name, or the stem of its file), the others FASTA
                         modifiers by name (isolate, specimen_voucher, country, collection_date, lat_lon, organism,
                         location, topology, genetic_code, definition, bioproject, biosample ...; case, spaces, '_'
                         and '-' do not matter); other columns are left out and reported
  --modifiers TEXT       the same for every record, as table2asn -j: "[location=mitochondrion] [topology=circular]";
                         a --source_table value wins over --modifiers, both over the record's own
  --product_table FILE   two columns, gene and product: the /product of the CDS, tRNA and rRNA with that /gene (any
                         case). Products still missing are filled from the gene name for the mitochondrial genes
                         (cox1 -> cytochrome c oxidase subunit I, trnL1 -> tRNA-Leu, rrnL -> 16S ribosomal RNA): a CDS
                         without one is an NCBI ERROR (SEQ_FEAT.MissingProteinName)
Cleaning (what other programs leave in a GenBank file and NCBI does not take; every change in the report):
  - qualifiers outside the INSDC feature table definition (Geneious /created_by /modified_by, SnapGene / ApE /label
    /ApEinfo_*, ...) are left out; a feature key outside it (Geneious 'ORF' ...) becomes misc_feature with the old key
    as a note
  - text is made ASCII (NCBI rejects anything else, SEQ_INST.NonAsciiAsn): accents dropped (Muller, Koln), mu -> u,
    typographic quotes and dashes plain; other characters removed. A file is read as UTF-8, else Windows-1252 /
    Latin-1 (older Windows programs)
  - a feature running past the end of the sequence is reported; a site between two bases (100^101) is written as the
    2-nt interval around it (as NCBI's gbf2tbl)
  - a file that is not readable GenBank / EMBL, and a record without sequence (CONTIG / WGS master records) are
    reported and skipped; the other inputs are written
  - genetic code: the CDS /transl_table; --genetic_code N for records without one (else, for a mitochondrion, the
    usual code of its phylum, reported)
--validate: table2asn (found by mito_env.py) on every pair -> <pair>.sqn (the file to upload; --sbt adds the submitter
template from https://submit.ncbi.nlm.nih.gov/genbank/template/submission/) and <pair>.val; NCBI ERROR / WARNING per
record in the report.

Usage:
  mitom gb2tbl -i final_annotation/ -o ncbi                      # all records -> ncbi/final_annotation.tbl + .fsa
  mitom gb2tbl -i a.gbf b.gb c.embl -o ncbi --name batch1         # three files -> ncbi/batch1.tbl + .fsa
  mitom gb2tbl -i genomes/ -o ncbi --split record --validate --sbt template.sbt
  mitom gb2tbl -i genomes/ -o ncbi --source_table samples.tsv --modifiers "[topology=circular]"

options:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        GenBank / EMBL files and / or folders, '-' for standard input (repeatable: -i a.gbf
                        b.embl -i more/)
  -o OUT, --out OUT     output folder
  --split {all,file,record}
                        all (default): one .tbl + .fsa with every record; file: one pair per input file;
                        record: one pair per record
  --name NAME           base name of the pair with --split all (default: the stem of a single input file,
                        else 'submission')
  --source_modifiers {yes,no}
                        copy the source qualifiers NCBI takes as FASTA modifiers (isolate, specimen_voucher,
                        country, collection_date ...) into the .fsa headers (default yes)
  --source_table SOURCE_TABLE
                        source information per record: tab-separated (NCBI .src) or CSV; first column the
                        record (SeqID, ID, LOCUS name or file stem), the others FASTA modifiers by name
                        (isolate, specimen_voucher, country, collection_date, lat_lon, organism, topology,
                        definition ...)
  --modifiers MODIFIERS
                        FASTA modifiers for every record, as table2asn -j: "[location=mitochondrion]
                        [topology=circular]" (a --source_table value wins)
  --product_table PRODUCT_TABLE
                        gene and product (two columns): the /product of the CDS / tRNA / rRNA with that
                        /gene (missing products of mitochondrial genes are filled from the gene name anyway)
  --genetic_code GENETIC_CODE
                        genetic code of every record (default: the CDS /transl_table; records without one:
                        the usual mitochondrial code of their phylum, reported)
  --validate            run NCBI table2asn on every pair: .sqn (the file to upload) and .val, errors in the
                        report
  --sbt SBT             NCBI submitter template (.sbt) for the .sqn files (with --validate)
  --table2asn TABLE2ASN
                        table2asn executable (default: found by mito_env.py)
```

## mitom fetch

```text
usage: mito_fetch.py [-h] (--taxon TAXON | --query QUERY) [-o OUT] [--name NAME] [--refseq_only]
                     [--min_len MIN_LEN] [--max_len MAX_LEN] [--max MAX] [--email EMAIL] [--api_key API_KEY]
                     [--dry_run] [--fasta]

mito_fetch.py - download complete mitochondrial genomes (GenBank flat files) from NCBI nuccore.

Records are selected by taxon (scientific name or taxid; the whole subtree) or by a raw Entrez query.
By default both RefSeq (NC_) and INSDC (GenBank) records are taken; `mitom dedup` later keeps one copy
of each genome (the RefSeq one when it exists). Output: one GenBank file with every record, a report
table, and the query that was used. A second run in the same directory only downloads records that
are not in the file yet (resume / update).

Usage:
  mitom fetch --taxon Annelida -o fetch_out                 # name or taxid, RefSeq + GenBank
  mitom fetch --taxon 6340 --refseq_only -o fetch_out       # RefSeq (NC_) only
  mitom fetch --query 'Nereididae[Organism] AND mitochondrion[Filter] AND "complete genome"[Title]' -o out
  mitom fetch --taxon Annelida -o fetch_out --fasta      # also fetch_out/Annelida.fasta
  then:  mitom dedup -i fetch_out/Annelida.gb -o dedup_out     or     mitom run -i fetch_out/Annelida.gb -o out

NCBI asks for an e-mail with every E-utilities request (--email, $NCBI_EMAIL, or "ncbi_email" in the
mito_env config file) and allows 10 requests/s with an API key (--api_key, $NCBI_API_KEY, "ncbi_api_key")
instead of 3/s without.

Outputs in -o:  <name>.gb  fetch_report.tsv  query.txt  fetch.log  (+ <name>.fasta with --fasta)

options:
  -h, --help         show this help message and exit
  --taxon TAXON      Scientific name or taxid; the whole subtree is taken (e.g. Annelida, 6340).
  --query QUERY      Raw Entrez query for nuccore (used instead of the taxon + mitogenome filters).
  -o OUT, --out OUT  Output folder (default: fetch_out).
  --name NAME        Base name of the GenBank file (default: taxon name, or 'ncbi' for --query).
  --refseq_only      RefSeq (NC_) records only.
  --min_len MIN_LEN  Minimum sequence length (bp).
  --max_len MAX_LEN  Maximum sequence length (bp).
  --max MAX          Stop after this many new records (test runs).
  --email EMAIL      E-mail for NCBI E-utilities (or $NCBI_EMAIL / config ncbi_email).
  --api_key API_KEY  NCBI API key (or $NCBI_API_KEY / config ncbi_api_key).
  --dry_run          Only print the query and the record count.
  --fasta            Also write <name>.fasta (every record of <name>.gb; headers >ACCESSION.V definition
                     topology=circular|linear, ready for mitom run). GenBank stays the main output: its
                     taxonomy sets the clade and genetic code, its header is copied into the results.
```

## mitom check

```text
usage: mito_check.py [-h] -i INPUT [-o OUT] [--genetic_code GENETIC_CODE] [--min_trna MIN_TRNA]

mito_check.py - audit the existing annotation of mitochondrial genome GenBank records.

A report only: it shows what is wrong with an existing annotation (e.g. records downloaded from NCBI, or a
final MitoM result); it does not change any file. `mitom run` re-annotates every record regardless.
An ERROR means the annotation is unreliable; a WARNING is a smaller problem or legitimate variation.

ERROR   no_annotation      no CDS feature at all
        missing_pcg        one of the 13 protein genes absent (cox1-3, cob, atp6, atp8, nad1-6, nad4l)
        missing_rrna       rrnS (12S) or rrnL (16S) absent
        duplicate_pcg      the same protein gene annotated twice (split genes COX1a/COX1b are a warning)
        internal_stop      stop codon inside a CDS (translated with its transl_table)
        frame              CDS length not a multiple of 3 without partial marker or transl_except
        pcg_length         PCG shorter than 60 % or longer than 160 % of the usual length
        pcg_overlap        two protein/rRNA genes overlap by more than 100 nt
        trna_count         fewer tRNAs than --min_trna (22 in most bilaterians)
WARNING no_transl_table    CDS without /transl_table
        start_codon        first codon is not a start codon of the table
        no_stop            last codon is not a stop and no partial marker / transl_except
        split_gene         gene annotated in two parts (COX1a / COX1b)
        trna_set           amino acids without a tRNA, or extra copies
        orientation        cox1 not on the plus strand, or not the first gene (the run re-orients)
        unknown_cds        a CDS that is not one of the 13 genes (ORF ...)

Usage:
  mitom check -i genomes.gb -o check_out [--genetic_code 5] [--min_trna 22]
Outputs in -o:  check_report.tsv (one row per record)  check_details.tsv (one row per issue)
                check_flagged.txt (records with an ERROR)

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        GenBank file or directory.
  -o OUT, --out OUT     Output folder (default: check_out).
  --genetic_code GENETIC_CODE
                        Translation table for CDS without /transl_table (default: the record's own tables).
  --min_trna MIN_TRNA   Fewer tRNAs than this is an ERROR (default 22; 0 = off).
```

## mitom compare

```text
usage: mito_compare.py [-h] --original ORIGINAL --mitom MITOM [-o OUT] [--genetic_code GENETIC_CODE]

mito_compare.py - compare the original annotation of each record (e.g. as downloaded from NCBI) with the
MitoM annotation of the same sequence, gene by gene.

MitoM rotates / reverse-complements a record so that cox1 comes first; the MitoM genes are mapped back onto the
original sequence first, so both annotations are compared in the original coordinates.

One row per gene (protein genes, rRNAs, tRNAs, introns):
  result  same            same name, strand and both ends
          ends_differ     same gene, at least one end differs (difference_nt = larger end shift)
          name_differs    same place, different gene name
          strand_differs  same place, other strand
          only_original   annotated only in the original record
          only_mitom      annotated only by MitoM
          same_place      tRNA-Leu / tRNA-Ser at the same place, which of the two genes cannot be told (no anticodon)

Usage:  mitom compare -o original_vs_mitom.tsv --original original_annotation/ --mitom final_annotation/ [--genetic_code 5]
        (files are paired by name: original_annotation/<acc>.gb <-> final_annotation/<acc>.gb)

options:
  -h, --help            show this help message and exit
  --original ORIGINAL   Original GenBank file or directory.
  --mitom MITOM         MitoM GenBank file or directory.
  -o OUT, --out OUT
  --genetic_code GENETIC_CODE
```

## mitom plot

```text
usage: mito_plot.py [-h] -i INPUT -o OUT [--fq1 FQ1] [--fq2 FQ2] [--threads THREADS] [--formats FORMATS]
                    [--style {circos,mitom}] [--min_noncoding MIN_NONCODING]

mito_plot.py - circular map of an annotated mitochondrial genome (MitoM), drawn with matplotlib.

Two styles. circos (default, --style circos), after mitoPlotR (K. Onthank, github.com/KirtOnthank/mitoPlotR, MIT
licence): a framed lane per strand (+ outside, - inside) with blocks coloured by type (protein gene, tRNA, rRNA,
other), gene names along the circle outside, tRNA codes, a position axis, GC content, read depth in 25 nt bins
scaled to the 99th percentile with the mean as a dashed line, the organism in the centre.

mitom (--style mitom; the default until 2026-10-30), from the outside in:
  genes       + strand outside the backbone circle, - strand inside, coloured by class (COX, NADH, ATP, CYTB,
              rRNA, tRNA; other = D-loop, misc_feature, rep_origin, repeat_region and CDS of other names, e.g.
              intron ORFs). Protein genes and rRNAs are labelled by /gene, tRNAs by their one-letter code (L1 / L2,
              S1 / S2 kept); crowded labels are fanned out with leader lines. Every part of a join is drawn (a gene
              across the origin, a frameshift join, the exons of a gene with an intron, linked by a line).
  GC content  deviation of each window from the genome mean (WINDOW nt windows every STEP nt)
  GC skew     (G - C) / (G + C) of the same windows
  depth       only with reads: minimap2 -x sr -c when it is installed (about 8x faster), else bwa mem; the
              reference span of every primary alignment; a circular record is mapped as sequence + its first
              ORIGIN_EXT nt and positions past the end are folded back onto the start
  centre      organism, length, topology, accession; a kb ruler runs inside the innermost ring
Position 1 is at the top and positions run clockwise; a linear record leaves a gap between its end and start.

Usage:
  mitom plot -i <GenBank file or folder> -o <out_dir> [--fq1 R1.fq.gz [--fq2 R2.fq.gz]] [--threads N]
             [--formats png,pdf] [--style circos|mitom]         (= python3 mito_plot.py ...)
Every record of the GenBank file, or of every GenBank file in the folder and its subfolders (.gb .gbk .gbf .genbank
.gbff), becomes <out_dir>/<record id>.png (300 dpi) and .pdf; reads are used only when the input holds one record.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        GenBank file, or a folder of GenBank files (subfolders included)
  -o OUT, --out OUT     folder for the maps
  --fq1 FQ1             reads for the depth ring (FASTQ, .gz allowed); only with a single record
  --fq2 FQ2             mates of --fq1
  --threads THREADS     aligner threads (minimap2, else bwa mem; default 1)
  --formats FORMATS     comma-separated (default png,pdf)
  --style {circos,mitom}
                        circos (default, after mitoPlotR): a framed lane per strand, names along the circle,
                        a position axis, GC content, read depth; mitom: genes on both sides of one circle,
                        GC content and GC skew rings
  --min_noncoding MIN_NONCODING
                        shortest unannotated stretch shaded and labelled 'NCR <length>' (default 100 nt; 0:
                        none)
```

## mitom nuc

```text
usage: mito_nuc.py [-h] {fetch,build,annotate} ...

mito_nuc.py - nuclear marker genes in contigs: rDNA (18S, ITS1, 5.8S, ITS2, 28S) and histones (H1, H2A, H2B, H3, H4).

Reference-based, like MitoFinder for mitogenomes: the contigs are searched with complete NCBI annotations and
each gene boundary is projected from the reference alignments onto the contig.

rDNA: blastn (megablast, word 16) of the candidate regions against the reference rDNA units; per reference the best
collinear chain of HSPs is kept (several chains = several copies on one contig, e.g. a tandem array); every
annotated reference boundary (18S/5.8S/28S start and end) is mapped through the alignment column by column
(or projected <= EXTRAP nt past the end of an HSP) and the TOPN best references vote, weighted by chain bit
score; the position with the most weight wins. A boundary that no reference reaches is projected without
alignment from the chain end; when that lands outside the contig the gene is partial (< / >). ITS1 / ITS2 are
the spacers between the rRNAs.
Histones: tblastn of the reference proteins, HSPs grouped into loci, exonerate protein2genome and protein2dna
(intronless) with the TOPK best proteins of the locus, best score wins (spliced models x INTRON_PEN per
intron); the CDS is extended in frame to the next stop codon and, when the model does not start with ATG at
protein position 1, to the farthest upstream in-frame ATG before a stop. Nuclear code (transl_table 1). Each copy
is classed replication-dependent (tandem cluster) / replacement (H3.3, spliced) / single-copy from H3 residues
31/87/89/90, introns and its coverage relative to the single-copy coverage (histone_class); the per-gene consensus
takes a cluster copy first (pick_histone).
Massive assemblies: both searches first run with one reference per genus (rDNA; low-complexity stretches do not
seed) / per gene and phylum (histones) against a BLAST database of all contigs; only the regions hit (+-PAD nt) go
to the steps above, each gene set in its own regions, and at most MAX_RDNA_REGIONS rDNA regions and
MAX_HIST_REGIONS histone regions per gene of them (best first; the time is bounded whatever the assembly).
2026-10-30: 30 SPAdes assemblies of 0.005-1.7 Gb took 6 s - 2.2 min on 8 threads, <= 1.6 GB (before: up to 20 min
and 16 GB on 16 threads; barrnap 0.9 on 8 threads 2 s - 4.6 min).

Per-gene consensus (multi-copy rDNA split over many contigs of one assembly): see consensus(); output
consensus/<gene>.fasta and consensus.tsv (which contigs were used, which were dropped and why).

Evidence (2026-10-23, benchmark/nuc_bench.py, knowledge/topics/nuclear_markers.md): NCBI metazoan records with
complete 18S-ITS1-5.8S-ITS2-28S annotation, one per species (378), references of the same genus left out:
boundaries exact 85.2 %, within 5 nt 92.6 % (barrnap 0.9 --kingdom euk: 1.8 % / 42.1 %); 677 complete histone
CDS: whole model (start, stop, exons, name) exact 97.5 %.

Usage:
  mitom nuc fetch -o nuc_refs                           # NCBI GenBank records: rdna.gb, histone.gb
  mitom nuc build --gb nuc_refs/rdna.gb nuc_refs/histone.gb [my_geneious.gb] -o models/nuc
  mitom nuc annotate -i contigs.fasta -o out [--genes rdna,histone] [--threads 16] [--lib models/nuc]
Outputs of annotate in -o: nuc_features.tsv (one row per feature), nuc_annotation.gb (contigs with features),
genes/<gene>.fasta (every copy; + histone .faa), consensus/<gene>.fasta + consensus.tsv (one sequence per
gene), nuc.log.

positional arguments:
  {fetch,build,annotate}
    fetch               download the NCBI reference records (rdna.gb, histone.gb)
    build               build the reference library from annotated GenBank files
    annotate            annotate rDNA and histone genes in contigs (FASTA or GenBank)

options:
  -h, --help            show this help message and exit
```

## mitom debug

```text
usage: mito_debug.py [-h] path

mito_debug.py - what happened in a MitoM run, read from its log (`mitom debug <output folder | log file>`).

The run log (run.log of `mitom run` / `import`, pipeline.log of `annotate` / `reorient`, the log of `dedup`,
`transfer_metadata`, `replace_organism`) holds lines of a fixed form (mitom_log.py); with the output folder the
reports are read too. The report:
  run        command, MitoM version, start, end, duration; status: finished (exit code), stopped on an error (the
             exception, the step it stopped in, the traceback's last lines), or did not finish (killed / still running)
  steps      time and result of each step and of its parts
  records    every record that failed or was left out, with the first error of its own log (reports/failed_logs/,
             work/ with --keep_intermediate yes / --debug)
  commands   every program that failed: exit code, time, folder, the end of its stderr
  programs   programs that were missing, and how to install them
  warnings   grouped by kind, with one example each
  NCBI       records with table2asn ERRORs (ncbi_submission/validation_messages.tsv)
  reproduce  the run again in debug mode, the failed program by hand, a failed record on its own
Exit code: 0 = the run finished without a problem, 1 = problems found (above), 2 = no MitoM log found.

Usage:  mitom debug mitom_out          (or: mitom debug mitom_out/run.log)

positional arguments:
  path        output folder of a MitoM run, or its log file (run.log, pipeline.log ...)

options:
  -h, --help  show this help message and exit
```
