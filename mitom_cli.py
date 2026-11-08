"""mitom_cli.py - the command line of mitom_pipeline.py (`mitom run`, `import`, `annotate` ...): build_parser()."""
from __future__ import annotations

import argparse
import os

from mito_common import old_spelling
from mito_taxa import ALL_CLADES


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
                            "genes by strand and type, GC content; with --infq also the read depth (minimap2 when installed, else bwa mem). "
                            "Default yes (about 2 s per record, --jobs in parallel).")
        g.add_argument("--plot_style", choices=["circos", "mitom"], default="circos",
                       help="Map style: circos (default, after mitoPlotR) = a framed lane per strand, names along the "
                            "circle, a position axis, GC content, read depth; mitom = genes on both sides of one circle, "
                            "GC content and GC skew rings.")
        g.add_argument("--ncr", choices=["yes", "no"], default="yes",
                       help="yes (default): every stretch of >= --ncr_min_len nt no feature covers becomes a "
                            "misc_feature /note=\"non-coding region NCR<k>, <length> bp: no gene annotated (MitoM)\" "
                            "in the GenBank and NCBI files, shaded and labelled on the maps (unnamed control region, "
                            "sequence that may hold unidentified genes); no: none (mito_ncr.py)")
        g.add_argument("--ncr_min_len", type=int, default=100,
                       help="shortest non-coding region written (nt, default 100; intergenic spacers are shorter)")
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
                       help="Taxonomic group (sets the genetic code and group rules). auto (default): per record - "
                            "the taxonomy and CDS transl_table of a GenBank record; FASTA and records without a "
                            "taxonomy: blastx to the reference proteins (mito_classify.py: phylum of the closest "
                            "genomes, their genetic code checked on the codons; no close relative: Animal and the "
                            "code of the codons). Several phyla / codes: `run` annotates each group with its own "
                            "clade and code and merges the results (reports/clade_decision.tsv); `import` and "
                            "`annotate` need one. Or give one of: "
                            f"{', '.join(ALL_CLADES)}.")
        g.add_argument("--genetic_code", type=int, default=None,
                       help="Genetic code (default with --clade auto: the CDS transl_table of each record, else the "
                            "code of its lineage or of the closest reference genomes; otherwise the code of the clade - 2 Chordata, 9 Echinodermata, "
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
                             "topology=circular|linear; header words ('circular', 'complete genome' / 'partial "
                             "genome'); an assembler's terminal overlap (removed): circular; else a sequence of >= 10 kb "
                             "circular (one or several in a file), a shorter one linear - and the protein genes at the "
                             "two ends checked against their reference proteins: a gene running from the end across to "
                             "the start makes the record circular, a gene cut at an end of a circular-labelled record "
                             "makes it linear, and so do jellyfish / hydroids / calcareous sponges (linear mtDNA) without "
                             "a gene across the ends; an unlabelled record no gene at the ends decides is circular with "
                             ">= 12 of the 13 protein genes and both rRNAs, else linear; changed records are annotated "
                             "again (reports/topology_decision.tsv). label: the label only. "
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
        g.add_argument("--mitos_args", default="",
                       help='Further runmitos options for --engine mitos2, as one quoted string written with "=", '
                            'e.g. --mitos_args="--best --finovl 30" (runmitos -h lists them). MitoM sets the input, output, '
                            '-c (--genetic_code), -r/--refdir (refseq89m, --refdir) and --linear (topology) itself; '
                            'default: none (runmitos defaults, --noplots).')

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
