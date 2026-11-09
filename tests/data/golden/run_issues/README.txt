MitoM output
============================================================

summary.tsv    one row per record: gene counts, orientation, round, NCBI errors - read this first
annotation/    the result: <accession>.gb per record (cox1 first, plus strand), all_records.gb, genes.fna (protein genes and rRNAs), proteins.faa, plots/ (genome maps)
submission/    NCBI submission files, one folder per record (<accession>/: .fsa .tbl .sqn .val); submission_summary.tsv, validation_messages.tsv
qc/            checks and decisions: original_vs_mitom.tsv and one report per step (deduplication, curation, ...)
input/         the input records unchanged (GenBank input only), same file names
logs/          run.log (everything the run did, with all commands), run_params.json, intron_fix.log, failed records

qc/
------------------------------------------------------------
  original_vs_mitom.tsv        original vs MitoM annotation, one row per gene (GenBank input)
  dedup_report.tsv             deduplication: one row per input record with a sequence, kept / dropped and why
  clade_decision.tsv           --clade auto: phylum, clade and genetic code of every record, from its taxonomy or by blastx
  topology_decision.tsv        circular or linear: each record's label, where it came from, the topology annotated, why
  reannotation_decision.tsv    after re-orientation: round 2 or not, and why
  curation_report.tsv          every boundary / name change: rule, old and new location, reason
  rrna_hmm_report.tsv          rRNAs searched / added by the per-phylum rRNA HMMs
  intron_report.tsv            introns found
  intron_ref_report.tsv        introns placed from the intron sites of annotated relatives (models/intron_sites)
  orf_report.tsv               ORFs of unknown function in unannotated sequence: written or to check, with the evidence
  metadata_transfer_log.tsv    where each record's GenBank header came from

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
