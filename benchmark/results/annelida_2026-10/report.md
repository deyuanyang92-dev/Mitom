# Benchmark: mitom, mitoz, mitos2, mitofinder (/media/deyuan/18T_B/yzz/mitom_work/bench/eval/all)

## Records

| tool | records | ok | failed | timeout | with obvious error | median s/record |
|---|---|---|---|---|---|---|
| mitom | 40 | 40 | 0 | 0 | 16 | 100 |
| mitoz | 40 | 40 | 0 | 0 | 36 | 244 |
| mitos2 | 40 | 40 | 0 | 0 | 39 | 353 |
| mitofinder | 40 | 39 | 1 | 0 | 36 | 206 |

## PCG (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom | 517 | 516 | 99.8% | 414 (80.1%) | 38 (87.4% cum.) | 64 | 0 | 0 | 1 | 2 |
| mitoz | 517 | 515 | 99.6% | 327 (63.2%) | 38 (70.6% cum.) | 150 | 0 | 0 | 2 | 5 |
| mitos2 | 517 | 514 | 99.4% | 277 (53.6%) | 33 (60.0% cum.) | 204 | 0 | 0 | 3 | 17 |
| mitofinder | 517 | 491 | 95.0% | 250 (48.4%) | 69 (61.7% cum.) | 172 | 0 | 0 | 26 | 2 |

## rRNA (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom | 80 | 80 | 100.0% | 24 (30.0%) | 21 (56.2% cum.) | 35 | 0 | 0 | 0 | 0 |
| mitoz | 80 | 80 | 100.0% | 14 (17.5%) | 1 (18.8% cum.) | 65 | 0 | 0 | 0 | 0 |
| mitos2 | 80 | 79 | 98.8% | 16 (20.0%) | 1 (21.2% cum.) | 62 | 0 | 0 | 1 | 0 |
| mitofinder | 80 | 60 | 75.0% | 10 (12.5%) | 12 (27.5% cum.) | 38 | 0 | 0 | 20 | 8 |

## tRNA (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom | 886 | 880 | 99.3% | 672 (75.8%) | 114 (88.7% cum.) | 16 | 0 | 5 | 1 | 10 |
| mitoz | 886 | 873 | 98.5% | 684 (77.2%) | 103 (88.8% cum.) | 10 | 0 | 4 | 9 | 13 |
| mitos2 | 886 | 871 | 98.3% | 664 (74.9%) | 114 (87.8% cum.) | 16 | 0 | 5 | 10 | 10 |
| mitofinder | 886 | 841 | 94.9% | 613 (69.2%) | 100 (80.5% cum.) | 11 | 0 | 22 | 23 | 14 |

## Same records for every tool (common_ok)

| tool | kind | truth | recall | both ends exact | within 3 nt |
|---|---|---|---|---|---|
| mitom | PCG | 504 | 99.8% | 80.2% | 87.5% |
| mitoz | PCG | 504 | 99.6% | 63.9% | 71.2% |
| mitos2 | PCG | 504 | 99.6% | 54.0% | 60.3% |
| mitofinder | PCG | 504 | 97.4% | 49.6% | 63.3% |
| mitom | rRNA | 78 | 100.0% | 30.8% | 56.4% |
| mitoz | rRNA | 78 | 100.0% | 17.9% | 19.2% |
| mitos2 | rRNA | 78 | 98.7% | 20.5% | 21.8% |
| mitofinder | rRNA | 78 | 76.9% | 12.8% | 28.2% |
| mitom | tRNA | 864 | 99.3% | 75.7% | 88.4% |
| mitoz | tRNA | 864 | 98.5% | 76.6% | 88.5% |
| mitos2 | tRNA | 864 | 98.4% | 74.9% | 87.6% |
| mitofinder | tRNA | 864 | 97.3% | 70.9% | 82.5% |

## Obvious errors by type

| error | mitom | mitoz | mitos2 | mitofinder |
|---|---|---|---|---|
| boundary_gt30 | 28 | 60 | 82 | 95 |
| missing_gene | 1 | 2 | 4 | 31 |
| name_error | 5 | 4 | 5 | 22 |
| duplicate_gene | 0 | 2 | 8 | 0 |
| internal_stop | 0 | 2 | 4 | 1 |
| failed | 0 | 0 | 0 | 1 |

## Error records per tool

**mitom** (16 records)

- NC_013569.1: boundary_gt30:nad1
- NC_027854.1: boundary_gt30:rrnL
- NC_037236.1: name_error:S, name_error:R
- NC_037947.1: boundary_gt30:nad4l, boundary_gt30:rrnL
- NC_049135.1: boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_049136.1: boundary_gt30:cox3, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_052721.1: boundary_gt30:nad5, boundary_gt30:rrnL
- NC_060759.1: name_error:S, name_error:R
- NC_061323.1: boundary_gt30:atp8, boundary_gt30:rrnL
- NC_061743.1: boundary_gt30:cox1, boundary_gt30:nad5, boundary_gt30:nad4, missing_gene:atp6, boundary_gt30:nad2, boundary_gt30:cox2, boundary_gt30:rrnL, boundary_gt30:W
- NC_064503.1: boundary_gt30:rrnS
- NC_065392.1: boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_069203.1: boundary_gt30:nad2, boundary_gt30:cob
- NC_069297.1: boundary_gt30:rrnS, name_error:S
- NC_072218.1: boundary_gt30:rrnL
- NC_072223.1: boundary_gt30:rrnL

**mitoz** (36 records)

- NC_006321.1: boundary_gt30:rrnL
- NC_006379.1: boundary_gt30:rrnS
- NC_007933.1: boundary_gt30:rrnL
- NC_011826.1: boundary_gt30:rrnL
- NC_013569.1: boundary_gt30:nad6, boundary_gt30:nad1, boundary_gt30:rrnL
- NC_026862.1: boundary_gt30:rrnL
- NC_026995.1: boundary_gt30:rrnL
- NC_027257.1: boundary_gt30:rrnL
- NC_029867.1: boundary_gt30:rrnL
- NC_034783.1: boundary_gt30:rrnL
- NC_037236.1: boundary_gt30:rrnL, name_error:S, name_error:R
- NC_037947.1: boundary_gt30:nad2, boundary_gt30:nad1, boundary_gt30:nad5, boundary_gt30:nad4l, boundary_gt30:rrnL
- NC_046790.1: boundary_gt30:rrnL
- NC_049135.1: boundary_gt30:rrnL, duplicate_gene:atp8
- NC_049136.1: boundary_gt30:cox3, boundary_gt30:rrnL
- NC_052721.1: boundary_gt30:nad5, boundary_gt30:rrnL
- NC_058285.1: boundary_gt30:rrnL
- NC_058588.1: boundary_gt30:rrnL
- NC_060759.1: boundary_gt30:rrnL, name_error:R
- NC_061323.1: missing_gene:atp8, boundary_gt30:nad6, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_061377.1: boundary_gt30:atp6, boundary_gt30:rrnL
- NC_061743.1: boundary_gt30:cox1, boundary_gt30:nad5, boundary_gt30:nad4, missing_gene:atp6, boundary_gt30:nad2, boundary_gt30:nad4l, boundary_gt30:cox3, boundary_gt30:nad1, boundary_gt30:cox2, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_063475.1: boundary_gt30:nad5, boundary_gt30:rrnL, internal_stop:atp6, internal_stop:atp6, duplicate_gene:atp6
- NC_065095.1: boundary_gt30:rrnL
- NC_065392.1: boundary_gt30:rrnL
- NC_069058.1: boundary_gt30:cox2
- NC_069059.1: boundary_gt30:cox2
- NC_069203.1: boundary_gt30:nad2, boundary_gt30:cob, boundary_gt30:rrnL
- NC_069204.1: boundary_gt30:nad4l, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_069297.1: name_error:S
- NC_071182.1: boundary_gt30:rrnL
- NC_072218.1: boundary_gt30:rrnL
- NC_072223.1: boundary_gt30:rrnL
- NC_072606.1: boundary_gt30:rrnL
- NC_082191.1: boundary_gt30:rrnL
- NC_082300.1: boundary_gt30:rrnL

**mitos2** (39 records)

- NC_006321.1: boundary_gt30:nad5, boundary_gt30:nad4, boundary_gt30:nad1
- NC_006379.1: boundary_gt30:atp6, boundary_gt30:nad5, boundary_gt30:rrnS
- NC_007933.1: boundary_gt30:nad6, boundary_gt30:rrnL
- NC_011826.1: boundary_gt30:cox1, boundary_gt30:nad3, boundary_gt30:rrnL
- NC_013569.1: boundary_gt30:cox1, boundary_gt30:atp8, boundary_gt30:cox3, missing_gene:nad6, boundary_gt30:atp6, boundary_gt30:nad1, boundary_gt30:rrnL
- NC_026862.1: boundary_gt30:nad5
- NC_026995.1: boundary_gt30:nad2
- NC_027257.1: boundary_gt30:nad4, boundary_gt30:nad2
- NC_027854.1: boundary_gt30:nad5, boundary_gt30:rrnL
- NC_029867.1: boundary_gt30:nad4
- NC_037236.1: boundary_gt30:nad6, name_error:S, name_error:R
- NC_037947.1: boundary_gt30:atp8, boundary_gt30:nad2, boundary_gt30:nad5, boundary_gt30:nad4l, boundary_gt30:rrnL
- NC_046790.1: boundary_gt30:cox1, boundary_gt30:rrnL
- NC_049135.1: boundary_gt30:cox3, boundary_gt30:nad3, internal_stop:nad5, duplicate_gene:atp8, duplicate_gene:nad5
- NC_049136.1: boundary_gt30:rrnL, duplicate_gene:nad5
- NC_052721.1: missing_gene:nad2, boundary_gt30:nad5, boundary_gt30:nad4, internal_stop:nad5, duplicate_gene:nad5
- NC_058285.1: boundary_gt30:atp8, boundary_gt30:rrnL
- NC_058588.1: boundary_gt30:nad5, boundary_gt30:nad2
- NC_059817.1: boundary_gt30:cox1, boundary_gt30:nad5
- NC_060759.1: boundary_gt30:nad5, name_error:S, name_error:R
- NC_061323.1: boundary_gt30:cox1, boundary_gt30:atp8, boundary_gt30:nad6, boundary_gt30:nad5, boundary_gt30:nad4, boundary_gt30:rrnL
- NC_061377.1: boundary_gt30:cox3, boundary_gt30:atp6, boundary_gt30:nad5, boundary_gt30:cox2, boundary_gt30:nad1, boundary_gt30:rrnL
- NC_061743.1: boundary_gt30:cox1, boundary_gt30:nad4l, boundary_gt30:cox2, boundary_gt30:rrnS, boundary_gt30:rrnL, boundary_gt30:W
- NC_063475.1: boundary_gt30:rrnL, internal_stop:atp6, internal_stop:atp6, duplicate_gene:atp6
- NC_064503.1: boundary_gt30:cox1, boundary_gt30:rrnL
- NC_065095.1: boundary_gt30:nad5
- NC_065392.1: boundary_gt30:nad5, boundary_gt30:rrnS, missing_gene:rrnL
- NC_066720.1: boundary_gt30:cox1
- NC_069058.1: boundary_gt30:cox2, duplicate_gene:nad5
- NC_069059.1: boundary_gt30:cox2, duplicate_gene:nad5
- NC_069203.1: missing_gene:nad2, boundary_gt30:cob
- NC_069204.1: boundary_gt30:nad5, boundary_gt30:rrnL
- NC_069297.1: name_error:S
- NC_071182.1: boundary_gt30:nad6
- NC_072218.1: boundary_gt30:nad5, boundary_gt30:atp6, boundary_gt30:rrnL
- NC_072223.1: boundary_gt30:rrnL
- NC_072606.1: boundary_gt30:nad5, boundary_gt30:cox1
- NC_082191.1: boundary_gt30:nad5
- NC_082300.1: duplicate_gene:nad5

**mitofinder** (37 records)

- NC_006321.1: boundary_gt30:cob, boundary_gt30:nad4, missing_gene:rrnS, missing_gene:rrnL
- NC_007933.1: boundary_gt30:cox1, missing_gene:atp8, missing_gene:nad6, missing_gene:nad4l, boundary_gt30:nad4, boundary_gt30:nad1, missing_gene:nad2, missing_gene:rrnS, missing_gene:rrnL, internal_stop:nad4
- NC_011826.1: boundary_gt30:nad6, boundary_gt30:nad5, boundary_gt30:rrnS, missing_gene:rrnL
- NC_013569.1: boundary_gt30:atp6, boundary_gt30:nad5, boundary_gt30:nad2
- NC_026862.1: boundary_gt30:nad2, boundary_gt30:rrnL
- NC_026995.1: missing_gene:atp8, missing_gene:nad6, boundary_gt30:nad5, boundary_gt30:nad4, boundary_gt30:nad1, missing_gene:rrnS, missing_gene:rrnL, name_error:S
- NC_027257.1: boundary_gt30:rrnL
- NC_027854.1: boundary_gt30:rrnL, name_error:S
- NC_029867.1: boundary_gt30:nad5
- NC_034783.1: missing_gene:atp8, boundary_gt30:nad2, boundary_gt30:rrnS, missing_gene:rrnL
- NC_037236.1: boundary_gt30:nad5, boundary_gt30:rrnL, name_error:S, name_error:R
- NC_037947.1: boundary_gt30:nad3, boundary_gt30:nad5, boundary_gt30:nad4, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_046790.1: missing_gene:atp8, missing_gene:nad6, boundary_gt30:nad4, boundary_gt30:nad2, missing_gene:rrnS, missing_gene:rrnL, name_error:S
- NC_049135.1: boundary_gt30:atp8, boundary_gt30:nad6, boundary_gt30:nad2, boundary_gt30:rrnS, missing_gene:rrnL
- NC_049136.1: missing_gene:atp8, boundary_gt30:cox3, boundary_gt30:nad6, boundary_gt30:nad3, boundary_gt30:nad4, boundary_gt30:nad1, boundary_gt30:nad5, boundary_gt30:nad2, missing_gene:rrnS, missing_gene:rrnL, name_error:S
- NC_052721.1: failed
- NC_058588.1: boundary_gt30:rrnL, name_error:S
- NC_059817.1: boundary_gt30:nad4, boundary_gt30:nad1, boundary_gt30:rrnS, missing_gene:rrnL, name_error:S
- NC_060759.1: boundary_gt30:nad5, boundary_gt30:rrnL, name_error:S, name_error:R
- NC_061323.1: missing_gene:atp8, boundary_gt30:cox3, missing_gene:atp6, boundary_gt30:nad5, boundary_gt30:nad2, boundary_gt30:rrnS, boundary_gt30:rrnL, name_error:S
- NC_061377.1: boundary_gt30:nad4, boundary_gt30:cox3, boundary_gt30:nad4l, boundary_gt30:rrnL
- NC_061743.1: boundary_gt30:cox1, boundary_gt30:nad4, boundary_gt30:cob, boundary_gt30:atp6, boundary_gt30:nad3, boundary_gt30:nad4l, boundary_gt30:cox2, missing_gene:rrnS, missing_gene:rrnL
- NC_063475.1: boundary_gt30:rrnL, name_error:S
- NC_064503.1: boundary_gt30:rrnL
- NC_065095.1: boundary_gt30:nad2, boundary_gt30:rrnL, name_error:S
- NC_065392.1: boundary_gt30:nad5, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_066720.1: boundary_gt30:rrnL, name_error:S
- NC_069058.1: name_error:L, name_error:S
- NC_069059.1: name_error:L, name_error:S
- NC_069203.1: boundary_gt30:nad1, boundary_gt30:nad2, boundary_gt30:cob, boundary_gt30:rrnL
- NC_069204.1: boundary_gt30:nad6, boundary_gt30:nad5, boundary_gt30:nad2, missing_gene:atp8, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_069297.1: boundary_gt30:nad4l, boundary_gt30:nad4, boundary_gt30:nad1, boundary_gt30:nad2, missing_gene:rrnL, missing_gene:rrnS
- NC_071182.1: boundary_gt30:nad5, boundary_gt30:nad1, boundary_gt30:rrnL, name_error:S
- NC_072218.1: boundary_gt30:nad1, boundary_gt30:nad2, boundary_gt30:rrnL, name_error:S
- NC_072223.1: boundary_gt30:nad6, name_error:S
- NC_072606.1: boundary_gt30:nad5, boundary_gt30:nad1, boundary_gt30:rrnL
- NC_082300.1: boundary_gt30:rrnL, name_error:S
