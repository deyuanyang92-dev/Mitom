MitoM output
============================================================

final_annotation/            MitoM annotation, one GenBank file per record (<accession>.gb) - the result
final_annotation_all.gb      the same records in one file
original_annotation/         the input records unchanged (GenBank input only), same file names
original_vs_mitom.tsv        original vs MitoM annotation, one row per gene
summary.tsv                  one row per record: gene counts, orientation, round, NCBI errors
ncbi_submission/             NCBI submission files per record (.fsa .tbl .sqn) and table2asn results
plots/                       circular genome map per record (.png, .pdf; --plot_style circos or mitom): genes, GC content, read depth with --infq
reports/                     one report per step (deduplication, round-2 decisions, curation, ...)
run.log                      everything the run did, with all commands

reports/
------------------------------------------------------------
  dedup_report.tsv             deduplication: one row per input record with a sequence, kept / dropped and why
  topology_decision.tsv        circular or linear: each record's label, where it came from, the topology annotated, why
  reannotation_decision.tsv    after re-orientation: round 2 or not, and why
  curation_report.tsv          every boundary / name change: rule, old and new location, reason
  multipart_report.tsv         records of one assembly / BioSample: parts of one genome (multipartite mtDNA) or copies
  pcg_hmm_report.tsv           missing protein genes searched with a profile HMM of the nearest relatives: found or not, E-value
  rrna_hmm_report.tsv          rRNAs searched / added by the per-phylum rRNA HMMs
  intron_report.tsv            introns found
  intron_ref_report.tsv        introns placed from the intron sites of annotated relatives (models/intron_sites)
  orf_report.tsv               ORFs of unknown function in unannotated sequence: written or to check, with the evidence
  metadata_transfer_log.tsv    where each record's GenBank header came from
  run_params.json              all options of this run
  intron_fix.log

summary.tsv columns
------------------------------------------------------------
  annotation_used: round 1, shifted = round-1 annotation moved with the sequence (passed the check);
                   round 2 = annotated again after re-orientation; not re-oriented = no orientation gene
  original_header_restored: the GenBank header of the input record (organism, references) was copied back
  ncbi_errors: NCBI table2asn ERROR count (0 = ready for submission)

original_vs_mitom.tsv result column
------------------------------------------------------------
  same / ends_differ / name_differs / strand_differs / only_original / only_mitom / same_place
  locations are on the original record's sequence; difference_nt = larger shift of the two ends

This run had several phyla / genetic codes (--clade auto): each group was annotated with its own clade and code (reports/clade_decision.tsv) and the results merged here; the log of each group is reports/groups/<group>.run.log (mitom debug reads them all).
