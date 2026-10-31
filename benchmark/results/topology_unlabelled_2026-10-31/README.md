# Circular or linear without a label (2026-10-31)

Question: can MitoM annotate correctly when the input does not say circular or linear (one FASTA, several sequences,
no topology=)?

Set (scripts/make_set.py, seed 1, --max 40 per GenBank file): 91 NCBI genomes of 7 phyla - the 76 intron-free
records and intron records of intron_gold_2026-10-29, 20 Hesionidae - each as 4 unlabelled sequences (truth.tsv):

| kind | sequence | truth |
|---|---|---|
| circ_orig | the genome as in NCBI | NCBI topology |
| circ_rot | rotated to a random start | NCBI topology |
| frag | a random piece of 3-12 kb | linear |
| gap | the genome without a random piece of 500-2000 nt, opened at the gap | linear |

Run per phylum (phylum.tsv) with its genetic code and clade: `scripts/run_set.sh <Mitom dir> <set dir> <label>`
(2 threads x 24 records); score: `scripts/score.py <set dir> <label>` (score_<commit>.tsv: per sequence).

| kind (truth) | n | a50402d | 2dae5b3 | cc8e82f |
|---|---|---|---|---|
| complete, NCBI start (circular) | 83 | 2 | 80 | 81 |
| complete, random start (circular) | 83 | 48 | 80 | 81 |
| NCBI 'partial / nearly complete' genome (linear; 7 + 7 rotations, 1 ambiguous 'complete genome' labelled linear left out) | 14 | 9 | 0 | 0 |
| fragment (linear) | 90 | 85 | 88 | 89 |
| gap (linear) | 90 | 88 | 61 | 56 |
| all | 360 | 64.4 % | 85.8 % | 84.3 % |

cc8e82f = 2dae5b3 without the rule ">= 21 tRNAs in phyla that keep all 22" (recomputed from the decision reasons; the
two complete NCBI genomes the rule made linear - Hexamermis agrotis, 17 tRNAs annotated, and Antechinus arktos, 20 -
rerun: circular). score_2dae5b3.tsv = run of 243f03a (the same code before the rebase).

Rules of cc8e82f (mitom_inputs, mitom_topology): an assembler's terminal overlap (end repeats start, >= 50 nt, not a
tandem repeat) -> circular, repeat removed; header words ('(circular)', 'complete genome' / 'partial genome', 'nearly
complete') are a label; otherwise >= 10 kb annotated circular first, then a gene across the ends -> circular, a gene
cut at an end -> linear, a lineage with linear mtDNA -> linear, else >= 12 of 13 protein genes and both rRNAs ->
circular, otherwise linear.

Limits: a contig missing only non-coding sequence or tRNAs holds every protein gene and rRNA and reads as a closed
circle (the 'gap' misses, NCBI 'partial' genomes missing the control region); lineages with reduced gene sets
(myxozoans, 8-11 protein genes) come out linear without a taxonomy. Reads spanning the junction would decide these.
