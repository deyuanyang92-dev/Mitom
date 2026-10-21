# Metazoa benchmark, 2026-10-22: MitoM vs MitoZ vs MITOS2 vs MitoFinder

动物线粒体基因组注释软件的系统对比：时间、内存、准确度（与 NCBI RefSeq 注释逐基因比较）、能否直接通过 NCBI 校验、功能、安装难易、便捷性。
A systematic comparison of animal mitogenome annotation tools: run time, memory, accuracy against the NCBI RefSeq annotation,
NCBI submission readiness, features, installation and convenience.

Full tables: `benchmark/results/metazoa_mini_2026-10/report.md` (summary.tsv, errors.tsv, timing.tsv, pcg_checks.tsv,
features.tsv.gz next to it). Reproduce: "How it was run" below.

## Data

24 RefSeq mitogenomes (`benchmark/sets/mito_bench_mini_2026-10.tsv`): 8 phyla x 3 records (Annelida, Arthropoda, Chordata,
Cnidaria, Echinodermata, Mollusca, Nematoda, Platyhelminthes), different classes / orders within a phylum where available,
genetic codes 2, 4, 5 and 9, 12.6-20.1 kb. Truth = the NCBI RefSeq annotation. Every tool received the bare sequence
(source feature only). RefSeq annotations are themselves imperfect (rRNA ends and start codons in particular), so "exact"
means "agrees with the NCBI curators".

The 24 were drawn from a 173-record set of 30 phyla (`benchmark/sets/mito_bench_2026-10.tsv`) on which MitoM alone was run
first; the error clusters of that run drove the MitoM changes listed under "What was fixed". The 4-tool comparison uses the
small set so that every tool could be run under the same conditions in about an hour.

## Tools

| tool | version | how it was run |
|---|---|---|
| MitoM | this repository (commit of this document) | `mitom run -i <bare.gb> -o out --clade <MitoZ clade or Animal> --genetic_code <code>`; all steps, curation, table2asn |
| MitoZ | 3.6 (bioconda) | `mitoz annotate --fastafiles in.fa --clade <clade> --genetic_code <code> --thread_number 2` |
| MITOS2 | 2.1.10 (bioconda `mitos`), reference data refseq89m | `runmitos.py -i in.fa -c <code> -o out -r refseq89m --refdir <dir> --noplots` |
| MitoFinder | 1.4.2 (GitHub source; MiTFi tRNAs) | `mitofinder -j <id> -a in.fa -r ref.gb -o <code> -p 2 -m 10 --max-contig-size <len+1000>`; ref.gb = the closest RefSeq record (same genus or family, never the species itself), one per record, listed in `benchmark/sets/mito_bench_mini_2026-10.mitofinder_refs.tsv`; CDS MitoFinder does not recognise ("hypothetical protein") were removed from the reference because MitoFinder refuses such a file |

Machine: 72 cores, 503 GB RAM. 8 records at a time per tool, 2 threads each. Per record: wall-clock time and the peak
resident memory of the whole process tree (sampled every second).

## Results

### Records

| tool | records | ok | failed | timeout | unsupported | with obvious error | median s/record | median peak RAM (MB) | max peak RAM (MB) |
|---|---|---|---|---|---|---|---|---|---|
| mitom | 24 | 24 | 0 | 0 | 0 | 14 | 175 | 732 | 1113 |
| mitoz | 24 | 18 | 0 | 0 | 6 | 17 | 249 | 492 | 719 |
| mitos2 | 24 | 24 | 0 | 0 | 0 | 24 | 385 | 953 | 1929 |
| mitofinder | 24 | 22 | 2 | 0 | 0 | 22 | 204 | 238 | 283 |

### Protein genes (all records a tool finished)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom | 306 | 305 | 99.7% | 240 (78.4%) | 18 (84.3% cum.) | 47 | 0 | 0 | 1 | 2 |
| mitoz | 230 | 230 | 100.0% | 174 (75.7%) | 17 (83.0% cum.) | 39 | 0 | 0 | 0 | 1 |
| mitos2 | 306 | 301 | 98.4% | 163 (53.3%) | 16 (58.5% cum.) | 122 | 0 | 0 | 5 | 5 |
| mitofinder | 306 | 271 | 88.6% | 175 (57.2%) | 29 (66.7% cum.) | 67 | 0 | 0 | 35 | 2 |

### rRNAs (all records a tool finished)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom | 46 | 46 | 100.0% | 13 (28.3%) | 9 (47.8% cum.) | 24 | 0 | 0 | 0 | 2 |
| mitoz | 34 | 30 | 88.2% | 8 (23.5%) | 1 (26.5% cum.) | 21 | 0 | 0 | 4 | 2 |
| mitos2 | 46 | 45 | 97.8% | 8 (17.4%) | 3 (23.9% cum.) | 34 | 0 | 0 | 1 | 6 |
| mitofinder | 46 | 39 | 84.8% | 3 (6.5%) | 8 (23.9% cum.) | 28 | 0 | 0 | 7 | 3 |

### tRNAs (all records a tool finished)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom | 464 | 453 | 97.6% | 319 (68.8%) | 58 (81.2% cum.) | 16 | 0 | 0 | 11 | 5 |
| mitoz | 393 | 365 | 92.9% | 272 (69.2%) | 40 (79.4% cum.) | 6 | 0 | 0 | 28 | 0 |
| mitos2 | 464 | 445 | 95.9% | 313 (67.5%) | 58 (80.0% cum.) | 16 | 0 | 0 | 19 | 4 |
| mitofinder | 464 | 406 | 87.5% | 282 (60.8%) | 50 (71.6% cum.) | 10 | 0 | 7 | 51 | 35 |

### Same records for every tool (the 17 records all four finished)

| tool | kind | truth | recall | both ends exact | within 3 nt |
|---|---|---|---|---|---|
| mitom | PCG | 218 | 100.0% | 82.1% | 89.9% |
| mitoz | PCG | 218 | 100.0% | 78.4% | 86.2% |
| mitos2 | PCG | 218 | 99.1% | 50.5% | 57.8% |
| mitofinder | PCG | 218 | 97.7% | 59.2% | 72.0% |
| mitom | rRNA | 32 | 100.0% | 37.5% | 65.6% |
| mitoz | rRNA | 32 | 87.5% | 25.0% | 28.1% |
| mitos2 | rRNA | 32 | 96.9% | 15.6% | 25.0% |
| mitofinder | rRNA | 32 | 90.6% | 9.4% | 28.1% |
| mitom | tRNA | 371 | 97.8% | 69.5% | 82.2% |
| mitoz | tRNA | 371 | 92.5% | 68.5% | 79.2% |
| mitos2 | tRNA | 371 | 96.0% | 68.2% | 80.9% |
| mitofinder | tRNA | 371 | 96.5% | 67.1% | 79.0% |

### Protein genes, both ends exact, per phylum

| phylum | records | mitom | mitoz | mitos2 | mitofinder |
|---|---|---|---|---|---|
| Annelida | 3 | 71.8% | 61.5% | 43.6% | 33.3% |
| Arthropoda | 3 | 97.4% | 89.7% | 53.8% | 69.2% |
| Chordata | 3 | 89.7% | 82.1% | 79.5% | 69.2% |
| Cnidaria | 3 | 84.6% | - | 69.2% | 87.2% |
| Echinodermata | 3 | 87.2% | 92.3% | 48.7% | 74.4% |
| Mollusca | 3 | 63.2% | 47.4% | 42.1% | 39.5% |
| Nematoda | 3 | 75.0% | 80.6% | 27.8% | 50.0% |
| Platyhelminthes | 3 | 56.8% | - | 59.5% | 32.4% |

### Obvious errors by type

| error | mitom | mitoz | mitos2 | mitofinder |
|---|---|---|---|---|
| boundary_gt30 | 27 | 23 | 79 | 39 |
| table2asn_error | 7 | 14 | 24 | 22 |
| missing_gene | 1 | 4 | 6 | 13 |
| internal_stop | 1 | 3 | 0 | 3 |
| name_error | 0 | 0 | 0 | 7 |
| duplicate_gene | 0 | 0 | 5 | 0 |
| failed | 0 | 0 | 0 | 2 |
| records with any error | 14 | 17 | 24 | 24 |

### Obvious errors by type, records every tool finished (17 records)

| error | mitom | mitoz | mitos2 | mitofinder |
|---|---|---|---|---|
| boundary_gt30 | 6 | 14 | 48 | 29 |
| table2asn_error | 2 | 13 | 17 | 17 |
| missing_gene | 0 | 4 | 3 | 8 |
| name_error | 0 | 0 | 0 | 7 |
| internal_stop | 0 | 3 | 0 | 3 |
| duplicate_gene | 0 | 0 | 3 | 0 |
| records with any error | 7 | 16 | 17 | 17 |


### Reading the tables

- *both ends exact*: start and end agree with RefSeq; *within 3 nt*: both ends within 3 nt.
- *obvious error*: a missing gene, an internal stop, a duplicate gene, a strand / name error, a boundary more than 30 nt off,
  or an NCBI table2asn ERROR; a record "with obvious error" has at least one.
- *table2asn_error*: the GenBank file the tool produced (MITOS2: converted from its feature table) was validated with NCBI
  table2asn; the count is records with at least one ERROR, i.e. not acceptable for submission as they are.
- *memory*: MitoM's peak (0.7 GB median) is the MITOS2 tRNA search it runs inside MitoZ: one `cmsearch` (Infernal 1.1,
  glocal) holds ~0.4 GB; MitoZ's own MiTFi (Infernal 1.0, Java) needs less, MITOS2 standalone runs several `cmsearch`
  at once (1.9 GB max). MitoM's own processes add ~0.25 GB (Python, Snakemake, `conda run`). This is the one measure
  where MitoZ and MitoFinder are lower.
- MitoZ has no protein profiles for Cnidaria and Platyhelminthes (6 records: *unsupported*). MitoFinder failed on 2
  records for which the closest RefSeq reference is a different order (no megablast hit). The "same records for every
  tool" table compares the 17 records all four tools finished.

## Features

| | MitoM | MitoZ 3.6 | MITOS2 2.1.10 | MitoFinder 1.4.2 |
|---|---|---|---|---|
| input | GenBank (clade and genetic code taken from the record's taxonomy) or FASTA | FASTA | FASTA | FASTA + a GenBank reference of a close relative (user supplied) |
| taxa | all Metazoa: MitoZ profiles for its 9 clades, MITOS2 for the rest and for genes MitoZ misses | 9 clades (Chordata, Arthropoda, Echinodermata, Annelida, Bryozoa, Mollusca, Nematoda, Nemertea, Porifera) | all Metazoa | whatever the reference has (a gene absent from the reference is not annotated: "ATP8 is not in the reference file") |
| genes | protein genes, rRNAs, tRNAs, group I / II introns, genes across the sequence start | protein genes, rRNAs, tRNAs | protein genes, rRNAs, tRNAs, introns, OL / OH | protein genes, rRNAs, tRNAs |
| boundary curation (start / stop codon, abbreviated stops, rRNA ends, duplicates, frame) | yes, every change logged with its rule | no | no | no |
| NCBI names and qualifiers (`/transl_except`, partial markers, standard gene names) | yes | partly (names) | no | partly |
| NCBI validation | table2asn on every record; ERROR count per record; `.sqn` / `.tbl` / `.fsa` submission files | no report | no | no |
| output | one GenBank per record + merged file; `summary.tsv`; original-vs-MitoM comparison; reports of every step; README in the folder | `.gbf`, `.tbl`, gene / CDS / rRNA / tRNA FASTA, circos figure, summary | bed, gff, FASTA, protein FASTA, NCBI feature table, gene order, figure (no GenBank) | `.gb`, `.tbl`, `.gff`, FASTA |
| batch | many records in one run (GenBank file or folder), deduplication (RefSeq copy kept), parallel (threads or Snakemake) | one FASTA per run | one sequence per run | one assembly per run |
| orientation | rotated / reverse-complemented to start with cox1 (fallback cox2, cox3, cob, nad1) | no | no | no |
| NCBI download / audit of existing annotation | `mitom fetch`, `mitom check`, `mitom compare` | no | no | no |
| transparency | 10 numbered steps printed with counts, times and result files; `run.log` | log | log | log |

## Installation

Verified on this machine (Linux x86_64, miniforge). Sizes are the installed conda environments.

| | MitoM | MitoZ | MITOS2 | MitoFinder |
|---|---|---|---|---|
| command | `bash install.sh` (creates the three environments below, downloads refseq89m with md5 check, installs the `mitom` command, writes the tool paths); `mitom doctor` checks every program; Dockerfile for other systems | `conda create -n mitoz -c bioconda mitoz=3.6` | `conda create -n mitos -c bioconda mitos=2.1.10`, then download the reference data (refseq89m, Zenodo 4284483) by hand and pass `--refdir` | bioconda package exists; here it did not run (missing `install.sh.ok` / config), so the GitHub source was used: Python 2.7 and Java, BLAST / ARWEN / tRNAscan-SE from conda, MiTFi's infernal compiled, two config files edited (absolute paths, shebang) |
| disk | 4.8 GB (mitoz) + 2.8 GB (mitos) + 0.4 GB (mitotools) + 93 MB reference data | 4.8 GB | 2.8 GB + 93 MB | 0.5 GB environment + 0.65 GB source |
| other requirements | Python 3 with Biopython | - | - | a GenBank reference per species, chosen by the user |
| platforms | Linux, Windows WSL2; macOS and others through Docker | Linux (bioconda) | Linux (bioconda) | Linux |

## Convenience

- What the user has to know: MitoM with GenBank input needs nothing (clade and code come from the record); with FASTA
  input the clade and genetic code. MitoZ: clade and code. MITOS2: code and the reference data path. MitoFinder: code and a
  close reference genome - the choice of the reference decides the result (two of 24 records failed because the closest
  RefSeq genome is a different order; a reference carrying a "hypothetical protein" CDS stops the program).
- What comes out: MitoM gives a GenBank file that is already in NCBI form (names, `/transl_except`, partial markers), the
  table2asn verdict per record and the submission files; the other tools give a draft that needs boundary editing and
  validation before submission (see *table2asn_error* above).
- Checking the result: MitoM prints each step with counts and writes `summary.tsv` (genes per record, orientation,
  NCBI ERROR count) and, for GenBank input, `original_vs_mitom.tsv` (every gene, original vs new coordinates), so a
  re-annotation of published records can be audited gene by gene.
- Run time per record (median): MitoM 175 s, MitoFinder 204 s, MitoZ 249 s, MITOS2 385 s (2 threads each). MitoM runs
  MitoZ and MITOS2 inside, so its memory is theirs (median 0.7 GB); MitoFinder needs the least (0.24 GB).

## What was fixed in MitoM during this benchmark

From the 173-record run (`benchmark/sets/mito_bench_2026-10.tsv`), each with a regression test and measured on the same
173 records before / after (A/B by re-curating the same core output):

1. **Protein gene ends stopping 100-500 nt short** (52 genes: nad5, nad4, nad2, cox1 in Acanthocephala, Brachiopoda,
   Ctenophora, Placozoa ...): the cores take a T inside the gene as an abbreviated stop once the C-terminal homology is lost;
   NCBI ends the gene at the first complete stop (50 of 52 exactly). P3 now extends an end to the first in-frame stop when
   nothing is annotated in between (at most 600 nt). 3' ends exact 1881 -> 1937 of 2178.
2. **Internal stop codons** (14 CDS): start and end in different frames (core start 1-2 nt off) -> P7 keeps the longer open
   reading frame; an antisense copy of a real gene -> D1 removes it; in-frame internal stops (frameshift / nonsense) stay
   reported for the user. 14 -> 8, none hidden.
3. **Stop codons near the next gene**: a complete stop 11-30 nt inside the next same-strand protein gene is what
   NCBI annotates (17 vs 6 abbreviated ends in the 173 records), so P3 now retracts to an abbreviated stop only beyond
   30 nt; an end that already is an abbreviated stop is no longer moved deeper than 10 nt into a following rRNA or
   antisense tRNA (NC_067900.1 nad1 had been extended 50 nt into the 16S rRNA: no NCBI record does that).
4. **Start codons**: the annelid rule "nad5 / nad4L start at the most upstream in-frame start" is wrong in other phyla
   (173 records: 52 wrong, 17 right) and is now applied for the Annelida clade only; a start overlapping the preceding tRNA is
   moved to the next start codon only when that is within 15 nt (moves <= 15 nt: 31 right / 12 wrong; longer: 14 / 28 -
   NCBI keeps the overlap, nematode genomes especially).
5. **Benchmark harness**: `/transl_except` positions were not lifted with the feature (236 false NCBI "NoStop" errors);
   Leu / Ser tRNA genes are now classified by anticodon before the submitter's number; MitoFinder references are stripped of
   CDS it does not recognise; peak memory is recorded.
6. A review of the new code (four independent reviewers) found 30 issues, all fixed: output folder handling (only MitoM's own
   files are moved or deleted, a second run into the same folder replaces the first), `mitom fetch` resume and incomplete
   downloads, `mitom check` without a genetic code, `mitom compare` import on a plain installation, P7 / P3 edge cases.

## Remaining MitoM differences from RefSeq (24-record set)

- rRNA ends: 24 of 46 rRNAs differ by more than 3 nt; the RefSeq records disagree among themselves on where 12S / 16S end
  (the other tools are further off).
- Start codons where NCBI chose a downstream in-frame start without protein evidence (Mollusca, Platyhelminthes).
- NC_027103.1 (*Rhodactis indosinensis*, Cnidaria, 20 kb): nad5 spans 19.8 kb with cox1 and nad1 inside its introns; wrong
  for every tool.
- tRNAs: all tools that use MiTFi (MitoM, MitoZ, MITOS2) find the same genes (~69 % both ends exact, ~81 % within 3 nt);
  nematode and acanthocephalan tRNAs with reduced structures are missed by every tool.

## How it was run

```bash
cd benchmark
# 24 records: 3 per phylum from the 173-record set (gb/ and bare/ come from mito_sample.py)
python3 mito_bench_run.py --set sets/mito_bench_mini_2026-10.tsv --root <root> --split all \
    --tools mitom,mitoz,mitos2,mitofinder --jobs 8 --threads 2 --timeout 7200 \
    --mitofinder <MitoFinder>/mitofinder --ref_pool <refs.gb>     # one closest RefSeq record per species
python3 mito_bench_eval.py --set sets/mito_bench_mini_2026-10.tsv --root <root> --split all \
    --tools mitom,mitoz,mitos2,mitofinder --table2asn
python3 bench_report.py --eval <root>/eval/all --tools mitom,mitoz,mitos2,mitofinder \
    --set sets/mito_bench_mini_2026-10.tsv --out report.md
```
