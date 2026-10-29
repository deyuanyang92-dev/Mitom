# MitoM, MitoZ, MITOS2: raw and after MitoM's curation (2026-10-28)

24 RefSeq (NC_) mitogenomes, 8 phyla x 3 (`benchmark/sets/mito_bench_mini_2026-10.tsv`), bare sequences, compared
gene by gene with the NCBI annotation (`mito_bench_eval.py --table2asn`). MitoZ has no clade for Cnidaria and
Platyhelminthes, so the side-by-side table uses the 18 records all three tools annotate.

| | PCG both ends exact | PCG within 3 nt | rRNA within 3 nt | tRNA found | tRNA both ends exact | records with NCBI ERROR | gene ends > 30 nt off | genes missing |
|---|---|---|---|---|---|---|---|---|
| MitoM `--curate no` | 58.7 % | 70.4 % | 29.4 % | 99.5 % | 70.5 % | 18 / 18 | 24 | 1 |
| MitoZ 3.6 raw | 75.7 % | 83.0 % | 26.5 % | 92.9 % | 69.2 % | 13 / 18 | 23 | 4 |
| MITOS2 2.1.10 raw | 50.0 % | 56.5 % | 23.5 % | 95.9 % | 68.2 % | 18 / 18 | 56 | 3 |
| MitoM (default) | 83.0 % | 89.6 % | 61.8 % | 99.5 % | 70.5 % | 0 / 18 | 10 | 0 |
| MitoZ + MitoM curation | 78.7 % | 84.8 % | 61.8 % | 92.9 % | 68.7 % | 0 / 18 | 15 | 4 |
| MITOS2 + MitoM curation | 56.5 % | 62.2 % | 64.7 % | 95.9 % | 67.7 % | 0 / 18 | 48 | 3 |

(230 protein genes, 34 rRNAs, 393 tRNAs in the 18 records; all 24 records: MitoM 83.0 / 52.2 / 99.1 %, MITOS2 +
curation 58.8 / 52.2 / 95.9 % for PCG exact / rRNA within 3 nt / tRNA found.)

- Raw MitoZ / MITOS2: the runs of `mitom_work/bench_mini` (unchanged tool versions). MitoM: current code
  (`--curate no` and default). "+ MitoM curation" = `mitom import` of that raw output with curation only
  (`--rrna_hmm no --intron_fix no --rescue no`): reorientation, end rules and NCBI names, NCBI files (`run.sh`).
- Curation is what makes every tool's output NCBI-valid (0 records with a table2asn ERROR) and moves rRNA ends
  (24-29 % -> 62-65 %). It cannot add genes a tool missed: tRNAs found and missing genes stay as the tool left them.
- MitoM's core ends a protein gene at its last aligned codon (the stop codon only when it follows directly) and
  leaves the stop search and several start rules to the curation step (mito_core.place_cds), so its raw PCG ends are
  less exact than MitoZ's; after curation MitoM is ahead (83.0 vs 78.7 %), finds more tRNAs (99.5 vs 92.9 %) and
  misses no gene (MitoZ 4).
- MITOS2's protein ends stay far off after curation (56.5 % exact; 48 ends more than 30 nt away).
