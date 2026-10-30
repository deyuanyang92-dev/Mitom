# MitoM benchmarks

Scripts that build the benchmark record sets, run MitoM and other tools on them, and score the results against the
NCBI annotation (or a verified gold standard). Results of past runs are in `results/`; the method descriptions and
numbers are in [../docs/METHODS.md](../docs/METHODS.md). Nothing here is needed to run MitoM.

## Typical run

```bash
# 1. a record set: NC_ records stratified by order from a local RefSeq pool (or mito_sample.py for the Metazoa set)
python3 make_large_set.py --pool Metazoa.gb --out bench --tsv sets/my_set.tsv --per_order 3
# 2. annotate the bare sequences with one or more tools; results are lifted to NCBI coordinates
python3 mito_bench_run.py --set sets/my_set.tsv --root bench --split all --tools mitom_native,mitos2 --jobs 8 --threads 2
#    mitom_lo = MitoM without the record's own proteins among the core's references (a new genome)
# 3. score them against NCBI
python3 mito_bench_eval.py --set sets/my_set.tsv --root bench --split all --tools mitom_native,mitos2   # -> bench/eval/all/
python3 bench_report.py --eval bench/eval/all --tools mitom_native,mitos2 --out report.md             # Markdown summary
```

## Scripts

| script | purpose |
|---|---|
| `mito_sample.py`, `make_large_set.py` | build record sets (`sets/*.tsv`, `gb/`, `bare/`) |
| `mito_bench_run.py` | run MitoM / MitoZ / MITOS2 / MITOS1 / MitoFinder on bare records; lift results to NCBI coordinates |
| `mito_bench_eval.py`, `bench_report.py` | gene-by-gene agreement with NCBI (PCG / rRNA ends, tRNAs, NCBI ERRORs); report |
| `eval_features.py`, `eval_cds.py` | per-gene boundaries; CDS-level agreement of intron models |
| `eval_trna.py` | tRNA methods of `mito_trna.py` against NCBI tRNAs |
| `eval_introns.py` | introns against NCBI or the gold standard (`--gb gold/gb`): exact / within 6 nt / more / missed / false, per gold class |
| `bench_intron_ref.py` | intron-site transfer: own / leave-one-out / leave-species / leave-genus |
| `train_start_model.py` | start-codon model of `mito_curate.py` P11 (`models/start_model/start_model.json`): candidates from `mito_curate.py --dump_starts` on the leave-out runs, labelled by NCBI, cross-validated by family |
| `compare.py`, `bidirectional.py`, `make_mitoz_like_input.py` | early intron evaluation vs real MFannot (2026-10-09/10, see below) |
| `nuc_bench.py`, `nuc_gold_bench.py`, `nuc_cm_gold.py`, `nuc_cm_speed.py`, `nuc_family_bench.py`, `nuc_fp_bench.py` | `mitom nuc`: rDNA boundaries, histone models, precision, Rfam end placement |
| `check_knowledge.py` | checks the knowledge-base Markdown files |
| `mito_names.py` | kept for old imports; the gene names live in the top-level `mito_names.py` |
| `trna_ml/` | MitoM-ML tRNA training / evaluation scripts and notes |

## Results (`results/`)

| folder | content |
|---|---|
| `six_tools_2026-10-22/` | 24 records x 6 tools (MitoM, MitoZ, MITOS2, MITOS1, MitoFinder, ...) |
| `metazoa_mini_2026-10/`, `metazoa_173_2026-10/`, `annelida_2026-10/` | MitoM on the 24-record, 173-record and annelid sets |
| `full_v30_2026-10-27/` | all benchmark sets (2,497 records): reports per set, tRNA tables, NCBI ERRORs |
| `curation_cmp_2026-10-28/` | MitoM, MitoZ, MITOS2 raw and after MitoM's curation |
| `intron_gold_2026-10-29/` | verified gold standard of the introns of 193 NCBI records; MitoM end to end against it |

## Early intron evaluation (mito_intron.py, 2026-10-09)

1. Positive set: NCBI nuccore `Metazoa[Organism] AND mitochondrion[Filter] AND complete genome[Title] AND intron[Feature key] AND refseq[Filter]` (85 records); negative set: 20 RefSeq annelid mitogenomes.
2. Split into `bench/fasta/<acc>.gb|.fasta` and `bench/records.tsv` (acc, set, genetic code, organism).
3. `python3 make_mitoz_like_input.py bench` -> `bench/sim/<acc>/` (introns removed, split genes cut to exon 1).
4. MitoM: `mito_intron.py -i bench/sim/<acc> -o bench/ours_full/<acc> --genetic_code <code> --bin <env>/bin` (and `--exclude_own_species` -> `ours_loo`).
5. MFannot: `docker run --rm -v bench:/data -w /data/mfannot nbeck/mfannot mfannot -g <code> -o /data/mfannot/<acc>.mf /data/fasta/<acc>.fasta`.
6. `python3 compare.py bench` -> summary + `bench/compare.tsv`.

## Early tRNA evaluation (eval_trna.py)

Records: `trna_sets.tsv` (85 intron-containing RefSeq mitogenomes, genetic code 4; 40 annelid RefSeq mitogenomes with 22
NCBI tRNAs). Strip all features but `source`, run `mito_trna.py -i <acc_dir> -o out/<acc> --genetic_code <code> --tools
mitos2,arwen,aragorn,trnascan`, then `python3 eval_trna.py <ncbi_gb_dir> out <label>`.
