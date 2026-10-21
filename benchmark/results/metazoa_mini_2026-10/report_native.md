# Benchmark: mitom, mitom_native, mitoz, mitos2, mitofinder (/media/deyuan/18T_B/yzz/mitom_work/bench_mini/eval/all)

## Records

| tool | records | ok | failed | timeout | unsupported | with obvious error | median s/record | median peak RAM (MB) | max peak RAM (MB) |
|---|---|---|---|---|---|---|---|---|---|
| mitom | 24 | 24 | 0 | 0 | 0 | 13 | 175 | 732 | 1113 |
| mitom_native | 24 | 24 | 0 | 0 | 0 | 14 | 109 | 320 | 1126 |
| mitoz | 24 | 18 | 0 | 0 | 6 | 17 | 249 | 492 | 719 |
| mitos2 | 24 | 24 | 0 | 0 | 0 | 24 | 385 | 953 | 1929 |
| mitofinder | 24 | 22 | 2 | 0 | 0 | 22 | 204 | 238 | 283 |

## PCG (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom | 306 | 305 | 99.7% | 240 (78.4%) | 18 (84.3% cum.) | 47 | 0 | 0 | 1 | 2 |
| mitom_native | 306 | 304 | 99.3% | 247 (80.7%) | 17 (86.3% cum.) | 40 | 0 | 0 | 2 | 2 |
| mitoz | 230 | 230 | 100.0% | 174 (75.7%) | 17 (83.0% cum.) | 39 | 0 | 0 | 0 | 1 |
| mitos2 | 306 | 301 | 98.4% | 163 (53.3%) | 16 (58.5% cum.) | 122 | 0 | 0 | 5 | 5 |
| mitofinder | 306 | 271 | 88.6% | 175 (57.2%) | 29 (66.7% cum.) | 67 | 0 | 0 | 35 | 2 |

## rRNA (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom | 46 | 46 | 100.0% | 13 (28.3%) | 9 (47.8% cum.) | 24 | 0 | 0 | 0 | 2 |
| mitom_native | 46 | 45 | 97.8% | 11 (23.9%) | 12 (50.0% cum.) | 22 | 0 | 0 | 1 | 3 |
| mitoz | 34 | 30 | 88.2% | 8 (23.5%) | 1 (26.5% cum.) | 21 | 0 | 0 | 4 | 2 |
| mitos2 | 46 | 45 | 97.8% | 8 (17.4%) | 3 (23.9% cum.) | 34 | 0 | 0 | 1 | 6 |
| mitofinder | 46 | 39 | 84.8% | 3 (6.5%) | 8 (23.9% cum.) | 28 | 0 | 0 | 7 | 3 |

## tRNA (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom | 464 | 453 | 97.6% | 319 (68.8%) | 58 (81.2% cum.) | 16 | 0 | 0 | 11 | 5 |
| mitom_native | 464 | 448 | 96.6% | 316 (68.1%) | 58 (80.6% cum.) | 16 | 0 | 1 | 15 | 6 |
| mitoz | 393 | 365 | 92.9% | 272 (69.2%) | 40 (79.4% cum.) | 6 | 0 | 0 | 28 | 0 |
| mitos2 | 464 | 445 | 95.9% | 313 (67.5%) | 58 (80.0% cum.) | 16 | 0 | 0 | 19 | 4 |
| mitofinder | 464 | 406 | 87.5% | 282 (60.8%) | 50 (71.6% cum.) | 10 | 0 | 7 | 51 | 35 |

## Same records for every tool (common_ok)

| tool | kind | truth | recall | both ends exact | within 3 nt |
|---|---|---|---|---|---|
| mitom | PCG | 218 | 100.0% | 82.1% | 89.9% |
| mitom_native | PCG | 218 | 99.5% | 82.6% | 89.9% |
| mitoz | PCG | 218 | 100.0% | 78.4% | 86.2% |
| mitos2 | PCG | 218 | 99.1% | 50.5% | 57.8% |
| mitofinder | PCG | 218 | 97.7% | 59.2% | 72.0% |
| mitom | rRNA | 32 | 100.0% | 37.5% | 65.6% |
| mitom_native | rRNA | 32 | 100.0% | 28.1% | 62.5% |
| mitoz | rRNA | 32 | 87.5% | 25.0% | 28.1% |
| mitos2 | rRNA | 32 | 96.9% | 15.6% | 25.0% |
| mitofinder | rRNA | 32 | 90.6% | 9.4% | 28.1% |
| mitom | tRNA | 371 | 97.8% | 69.5% | 82.2% |
| mitom_native | tRNA | 371 | 96.8% | 69.0% | 81.7% |
| mitoz | tRNA | 371 | 92.5% | 68.5% | 79.2% |
| mitos2 | tRNA | 371 | 96.0% | 68.2% | 80.9% |
| mitofinder | tRNA | 371 | 96.5% | 67.1% | 79.0% |

## Protein genes, both ends exact, per phylum

| phylum | records | mitom | mitom_native | mitoz | mitos2 | mitofinder |
|---|---|---|---|---|---|---|
| Annelida | 3 | 71.8% | 76.9% | 61.5% | 43.6% | 33.3% |
| Arthropoda | 3 | 97.4% | 94.9% | 89.7% | 53.8% | 69.2% |
| Chordata | 3 | 89.7% | 87.2% | 82.1% | 79.5% | 69.2% |
| Cnidaria | 3 | 84.6% | 92.3% | - | 69.2% | 87.2% |
| Echinodermata | 3 | 87.2% | 87.2% | 92.3% | 48.7% | 74.4% |
| Mollusca | 3 | 63.2% | 52.6% | 47.4% | 42.1% | 39.5% |
| Nematoda | 3 | 75.0% | 80.6% | 80.6% | 27.8% | 50.0% |
| Platyhelminthes | 3 | 56.8% | 73.0% | - | 59.5% | 32.4% |

## Obvious errors by type

| error | mitom | mitom_native | mitoz | mitos2 | mitofinder |
|---|---|---|---|---|---|
| boundary_gt30 | 27 | 23 | 23 | 79 | 39 |
| table2asn_error | 2 | 1 | 13 | 24 | 22 |
| missing_gene | 1 | 3 | 4 | 6 | 13 |
| name_error | 0 | 1 | 0 | 0 | 7 |
| internal_stop | 1 | 0 | 3 | 0 | 3 |
| duplicate_gene | 0 | 0 | 0 | 5 | 0 |
| failed | 0 | 0 | 0 | 0 | 2 |
| records with any error | 13 | 14 | 17 | 24 | 24 |

## Obvious errors by type, records every tool finished (17 records)

| error | mitom | mitom_native | mitoz | mitos2 | mitofinder |
|---|---|---|---|---|---|
| boundary_gt30 | 6 | 7 | 14 | 48 | 29 |
| table2asn_error | 0 | 0 | 12 | 17 | 17 |
| missing_gene | 0 | 1 | 4 | 3 | 8 |
| name_error | 0 | 0 | 0 | 0 | 7 |
| internal_stop | 0 | 0 | 3 | 0 | 3 |
| duplicate_gene | 0 | 0 | 0 | 3 | 0 |
| records with any error | 6 | 7 | 16 | 17 | 17 |

## Error records per tool

**mitom** (13 records)

- NC_001861.1: boundary_gt30:rrnS
- NC_013831.1: boundary_gt30:rrnL
- NC_027103.1: missing_gene:cox1, boundary_gt30:nad2, boundary_gt30:nad6, boundary_gt30:rrnL, internal_stop:nad5, table2asn_error
- NC_031213.1: boundary_gt30:nad5, boundary_gt30:rrnL
- NC_034715.1: boundary_gt30:rrnL
- NC_035255.1: boundary_gt30:cob, boundary_gt30:cox3, boundary_gt30:nad6, boundary_gt30:nad4, boundary_gt30:rrnS
- NC_046775.1: boundary_gt30:rrnL, table2asn_error
- NC_062047.1: boundary_gt30:cox1, boundary_gt30:cox2, boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_063787.1: boundary_gt30:rrnL
- NC_067900.1: boundary_gt30:rrnS
- NC_067930.1: boundary_gt30:rrnS
- NC_068797.1: boundary_gt30:rrnL
- NC_083122.1: boundary_gt30:cox1, boundary_gt30:nad6, boundary_gt30:atp6, boundary_gt30:nad5, boundary_gt30:rrnS

**mitom_native** (14 records)

- NC_001861.1: boundary_gt30:rrnS
- NC_005928.1: boundary_gt30:rrnL
- NC_027103.1: missing_gene:nad5, boundary_gt30:rrnL
- NC_031213.1: boundary_gt30:nad5, boundary_gt30:rrnL
- NC_034715.1: boundary_gt30:rrnL
- NC_035255.1: boundary_gt30:cob, boundary_gt30:nad3, boundary_gt30:cox3, boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_046775.1: boundary_gt30:rrnL, table2asn_error
- NC_062047.1: boundary_gt30:nad6, boundary_gt30:cox1, name_error:Q
- NC_067900.1: boundary_gt30:rrnS
- NC_067930.1: boundary_gt30:rrnS
- NC_068797.1: boundary_gt30:rrnL
- NC_069205.1: boundary_gt30:nad2
- NC_083122.1: boundary_gt30:cox1, boundary_gt30:cox2, boundary_gt30:nad4, boundary_gt30:rrnL, missing_gene:rrnS
- NC_086642.1: boundary_gt30:nad4l, missing_gene:atp8

**mitoz** (17 records)

- NC_001861.1: missing_gene:rrnS, missing_gene:rrnL
- NC_005928.1: boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_007933.1: boundary_gt30:rrnL
- NC_009732.1: table2asn_error
- NC_013831.1: missing_gene:rrnL, missing_gene:rrnS, internal_stop:nad4, table2asn_error
- NC_023772.1: boundary_gt30:rrnL, table2asn_error
- NC_045231.1: table2asn_error
- NC_053107.1: table2asn_error
- NC_063787.1: table2asn_error
- NC_065213.1: boundary_gt30:rrnL, table2asn_error
- NC_065650.1: table2asn_error
- NC_067900.1: boundary_gt30:rrnL, table2asn_error
- NC_067930.1: boundary_gt30:rrnL, boundary_gt30:rrnS, table2asn_error
- NC_068797.1: boundary_gt30:rrnL
- NC_069205.1: boundary_gt30:rrnL, table2asn_error
- NC_083122.1: boundary_gt30:cox1, boundary_gt30:nad6, boundary_gt30:atp6, boundary_gt30:cox2, boundary_gt30:cox3, boundary_gt30:nad5, boundary_gt30:cob, boundary_gt30:nad2, boundary_gt30:rrnS, table2asn_error
- NC_086642.1: boundary_gt30:cox2, boundary_gt30:atp6, boundary_gt30:rrnL, boundary_gt30:rrnS, internal_stop:nad5, internal_stop:atp6, table2asn_error

**mitos2** (24 records)

- NC_001861.1: boundary_gt30:nad2, boundary_gt30:nad4, boundary_gt30:cox1, boundary_gt30:nad6, boundary_gt30:cob, boundary_gt30:atp6, boundary_gt30:rrnS, missing_gene:rrnL, duplicate_gene:rrnL, duplicate_gene:rrnS, table2asn_error
- NC_005928.1: boundary_gt30:cox2, missing_gene:nad6, boundary_gt30:cox3, boundary_gt30:rrnL, table2asn_error
- NC_007933.1: boundary_gt30:nad6, boundary_gt30:rrnL, table2asn_error
- NC_009732.1: boundary_gt30:nad2, duplicate_gene:nad5, table2asn_error
- NC_013831.1: boundary_gt30:nad5, boundary_gt30:atp6, boundary_gt30:cob, boundary_gt30:cox3, boundary_gt30:rrnL, boundary_gt30:rrnS, table2asn_error
- NC_023772.1: boundary_gt30:cox2, boundary_gt30:atp6, boundary_gt30:cox3, boundary_gt30:nad5, table2asn_error
- NC_027103.1: missing_gene:nad5, missing_gene:cox1, boundary_gt30:nad2, boundary_gt30:nad6, boundary_gt30:rrnL, boundary_gt30:rrnS, duplicate_gene:nad5, table2asn_error
- NC_031213.1: boundary_gt30:nad2, boundary_gt30:nad5, boundary_gt30:cox1, boundary_gt30:rrnL, table2asn_error
- NC_034715.1: boundary_gt30:cox3, boundary_gt30:rrnL, boundary_gt30:rrnS, table2asn_error
- NC_035255.1: missing_gene:cob, boundary_gt30:nad3, boundary_gt30:cox3, boundary_gt30:nad6, boundary_gt30:cox2, boundary_gt30:rrnL, boundary_gt30:rrnS, table2asn_error
- NC_045231.1: boundary_gt30:nad2, table2asn_error
- NC_046775.1: boundary_gt30:nad2, boundary_gt30:nad5, table2asn_error
- NC_053107.1: table2asn_error
- NC_062047.1: boundary_gt30:nad6, boundary_gt30:nad5, boundary_gt30:cox1, boundary_gt30:rrnL, duplicate_gene:rrnL, table2asn_error
- NC_063787.1: boundary_gt30:nad1, boundary_gt30:cox2, boundary_gt30:atp6, boundary_gt30:cox3, boundary_gt30:cob, table2asn_error
- NC_065213.1: boundary_gt30:cox1, boundary_gt30:nad6, table2asn_error
- NC_065650.1: boundary_gt30:nad2, table2asn_error
- NC_067900.1: boundary_gt30:cox2, boundary_gt30:cox3, boundary_gt30:nad5, boundary_gt30:nad1, boundary_gt30:rrnS, table2asn_error
- NC_067930.1: boundary_gt30:nad2, boundary_gt30:nad5, table2asn_error
- NC_068797.1: boundary_gt30:nad5, boundary_gt30:cox2, boundary_gt30:atp6, table2asn_error
- NC_069205.1: boundary_gt30:cox1, table2asn_error
- NC_083035.1: table2asn_error
- NC_083122.1: boundary_gt30:cox1, boundary_gt30:cox2, boundary_gt30:cox3, boundary_gt30:nad5, boundary_gt30:cob, boundary_gt30:nad2, boundary_gt30:rrnL, boundary_gt30:rrnS, table2asn_error
- NC_086642.1: boundary_gt30:nad5, missing_gene:nad4l, boundary_gt30:cob, boundary_gt30:cox2, boundary_gt30:nad2, boundary_gt30:rrnL, table2asn_error

**mitofinder** (24 records)

- NC_001861.1: boundary_gt30:nad2, boundary_gt30:nad4, boundary_gt30:nad4l, boundary_gt30:nad3, boundary_gt30:nad5, table2asn_error
- NC_005928.1: name_error:L, table2asn_error
- NC_007933.1: boundary_gt30:cox1, missing_gene:atp8, missing_gene:nad6, missing_gene:nad4l, boundary_gt30:nad4, boundary_gt30:nad1, missing_gene:nad2, missing_gene:rrnS, missing_gene:rrnL, internal_stop:nad4, table2asn_error
- NC_009732.1: boundary_gt30:rrnS, name_error:S, name_error:L, table2asn_error
- NC_013831.1: boundary_gt30:nad3, table2asn_error
- NC_023772.1: table2asn_error
- NC_027103.1: missing_gene:nad5, missing_gene:cox1, boundary_gt30:nad4, boundary_gt30:rrnS, table2asn_error
- NC_031213.1: boundary_gt30:cox1, boundary_gt30:rrnL, table2asn_error
- NC_034715.1: boundary_gt30:rrnL, table2asn_error
- NC_035255.1: failed
- NC_045231.1: table2asn_error
- NC_046775.1: boundary_gt30:rrnL, table2asn_error
- NC_053107.1: boundary_gt30:rrnL, name_error:S, name_error:L, table2asn_error
- NC_062047.1: missing_gene:nad6, boundary_gt30:nad5, boundary_gt30:cob, missing_gene:atp6, missing_gene:nad2, boundary_gt30:rrnL, boundary_gt30:rrnS, table2asn_error
- NC_063787.1: boundary_gt30:nad5, boundary_gt30:rrnL, table2asn_error
- NC_065213.1: table2asn_error
- NC_065650.1: boundary_gt30:nad5, table2asn_error
- NC_067900.1: boundary_gt30:nad5, boundary_gt30:nad1, table2asn_error
- NC_067930.1: boundary_gt30:atp8, boundary_gt30:rrnL, table2asn_error
- NC_068797.1: boundary_gt30:nad6, boundary_gt30:nad2, boundary_gt30:nad4l, missing_gene:atp8, boundary_gt30:rrnS, missing_gene:rrnL, table2asn_error
- NC_069205.1: boundary_gt30:nad5, boundary_gt30:nad1, boundary_gt30:rrnL, table2asn_error
- NC_083035.1: boundary_gt30:nad5, name_error:S, name_error:L, table2asn_error
- NC_083122.1: failed
- NC_086642.1: boundary_gt30:nad1, boundary_gt30:nad4l, boundary_gt30:rrnL, internal_stop:atp8, internal_stop:nad4, table2asn_error
