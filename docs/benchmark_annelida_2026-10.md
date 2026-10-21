# Annelida benchmark, 2026-10-21: MitoM vs MitoZ vs MITOS2 vs MitoFinder

Full tables: `benchmark/results/annelida_2026-10/report.md` (summary.tsv, errors.tsv, timing.tsv, features.tsv.gz
next to it). Reproduce: see "How it was run" below.

## Data

40 annelid RefSeq mitogenomes (`benchmark/sets/annelida_2026-10.tsv`): drawn from the 167 Annelida NC_ records at NCBI
(2026-10-21), at most 4 per order, 15 orders, 29 families, genetic code 5, 13.7-21.9 kb. Truth = the NCBI RefSeq
annotation. Every tool received the bare sequence (source feature only). The RefSeq annotations are themselves
imperfect (rRNA ends in particular), so "exact" means "agrees with the NCBI curators", not "biologically correct".

## Tools

| tool | version / mode |
|---|---|
| MitoM | this repository, `mitom run --clade Annelida-segmented-worms --genetic_code 5` (MitoZ 3.6 + MITOS2 2.1.10 cores, all MitoM steps, curation; run with `--if_reannotation_2 yes`, the default at the time - `auto`, the default since, gave identical annotations in the A/B of knowledge/experiments/reannotation) |
| MitoZ | 3.6 `annotate`, unchanged |
| MITOS2 | 2.1.10 `runmitos`, full mode, refseq89m |
| MitoFinder | 1.4.2 (GitHub source, MiTFi tRNAs), `-a` on the bare sequence; reference = the closest annotated Annelida record at NCBI (longest shared lineage, same species excluded) |

## Results (all 40 records)

Records with an *obvious error* (missing gene, internal stop, duplicate gene, strand or name error, boundary > 30 nt, run failure):

| tool | records ok | with obvious error | median s/record |
|---|---|---|---|
| MitoM | 40 | 16 | 100 |
| MitoZ | 40 | 36 | 244 |
| MITOS2 | 40 | 39 | 353 |
| MitoFinder | 39 | 34 | 206 |

Protein genes (517 in the truth):

| tool | recall | both ends exact | within 3 nt | > 3 nt | missed | extra |
|---|---|---|---|---|---|---|
| MitoM | 99.8% | 80.1% | 87.4% | 64 | 1 | 2 |
| MitoZ | 99.6% | 63.2% | 70.6% | 150 | 2 | 5 |
| MITOS2 | 99.4% | 53.6% | 60.0% | 204 | 3 | 17 |
| MitoFinder | 95.0% | 48.4% | 61.7% | 172 | 26 | 2 |

rRNAs (80):

| tool | recall | both ends exact | within 3 nt |
|---|---|---|---|
| MitoM | 100% | 30.0% | 56.2% |
| MitoZ | 100% | 17.5% | 18.8% |
| MITOS2 | 98.8% | 20.0% | 21.2% |
| MitoFinder | 75.0% | 12.5% | 27.5% |

tRNAs (886): all three MiTFi-based tools are equal (MitoM 75.8% exact / 88.7% within 3 nt, MitoZ 77.2 / 88.8,
MITOS2 74.9 / 87.8); MitoFinder 69.2 / 80.5 with 22 Leu/Ser name errors.

Obvious errors by type:

| error | MitoM | MitoZ | MITOS2 | MitoFinder |
|---|---|---|---|---|
| boundary > 30 nt | 28 | 60 | 82 | 95 |
| missing gene | 1 | 2 | 4 | 31 |
| name error (tRNA) | 5 | 4 | 5 | 22 |
| duplicate gene | 0 | 2 | 8 | 0 |
| internal stop | 0 | 2 | 4 | 1 |
| failed | 0 | 0 | 0 | 1 |

## What was fixed in MitoM during this benchmark

The first pass (same cores) had MitoM at 74% exact protein ends. Three curation rules came out of the error clusters
(`mito_curate.py`), all regression-tested:

1. **P3, stop codon**: the first complete stop was taken even when it lay inside the next protein gene (cox1 → cox2,
   nad6 → cob: 22-85 nt overlaps). NCBI records end such genes with an abbreviated T/TA stop. Now, when no tRNA
   follows and the first complete stop is > 10 nt inside the next protein gene (bicistronic atp8-atp6 and nad4L-nad4
   excepted), the annotated T/TA end is kept, else the T/TA right before the next gene is used.
2. **D1, duplicate gene**: MitoZ annotated atp8 twice on a repeated region (NC_049135.1); one copy is kept.
3. **P1, start codon, per gene**: moving the start to the most upstream in-frame start codon (criterion i taken
   literally, `--aggressive_start`) is wrong for cob / cox2 / cox3 / nad4 (NCBI keeps the downstream start) but
   right for nad5 and nad4L in every differing record (7 of 7). `mitom run` passes `--aggressive_start_genes ND5,ND4L`
   for the Annelida clade only: on 173 records of 30 phyla the same rule was wrong 52 times and right 17 times
   (`docs/benchmark_metazoa_2026-10.md`), so the mito_curate.py default is empty.

Also fixed in the benchmark itself: NCBI writes Leu/Ser tRNA genes either with the anticodon (`trnL(tag)`) or the
codon (`trnL1(cta)`); the evaluator now compares the two Leu / Ser genes by codon family (`mito_bench_eval.ls_class`),
which removed 7 false name errors per tool. MitoFinder references are now restricted to records that carry all genes
(an UNVERIFIED record without CDS had been chosen), and a contig that MitoFinder trims at the circular junction is
padded back before lift-over.

## Remaining MitoM differences from the RefSeq annotation (16 records)

- rRNA ends (13 of 28 boundary differences): the NCBI records themselves disagree with each other on where 12S / 16S
  end (gaps of 30-80 nt, or 12S/16S overlapping each other by 12 nt in NC_049135/NC_049136); no rule can settle these.
- Start codons where NCBI chose a later in-frame start than criterion i gives (nad4l NC_037947, cox3 NC_049136,
  cob NC_069203, nad1 NC_013569 end): ambiguous without protein evidence.
- NC_061743.1 (Lumbriculida, 21.9 kb): wrong for every tool (MitoZ 1 of 13 genes exact, MITOS2 5, MitoM 4).
- Ser/Arg tRNA swap in NC_037236 / NC_060759 (NCBI `trnS(CGA)` where the anticodon UCG reads Arg): MITOS2 and MitoZ
  call the same genes; the NCBI naming looks inconsistent.
- Eight tandem tRNA-Pro copies in NC_060759 are all annotated by MITOS2 (NCBI annotates one).

## How it was run

```bash
cd benchmark
python3 mito_sample.py --out <root> --phyla Annelida --quota 40 --order_cap 4 --tsv sets/annelida_2026-10.tsv
python3 mito_bench_run.py --set sets/annelida_2026-10.tsv --root <root> --split all \
    --tools mitom,mitoz,mitos2,mitofinder --jobs 20 --threads 6 \
    --mitofinder <MitoFinder>/mitofinder --ref_pool <Annelida.gb from mitom fetch>
python3 mito_bench_eval.py --set sets/annelida_2026-10.tsv --root <root> --split all --tools mitom,mitoz,mitos2,mitofinder
python3 bench_report.py --eval <root>/eval/all --out report.md
```

Machine: 72 cores; wall time about 50 minutes for 4 tools x 40 records at 20 parallel jobs.
