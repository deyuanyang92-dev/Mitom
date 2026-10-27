# Benchmark: mitom_native, mitom_native_v27, mitos2 (/media/deyuan/18T_B/yzz/mitom_work/eval_v30/anchors/eval/all)

## Records

| tool | records | ok | failed | timeout | unsupported | with obvious error | median s/record | median peak RAM (MB) | max peak RAM (MB) |
|---|---|---|---|---|---|---|---|---|---|
| mitom_native | 13 | 13 | 0 | 0 | 0 | 3 | 175 | 291 | 473 |
| mitom_native_v27 | 13 | 13 | 0 | 0 | 0 | 5 | 224 | 195 | 474 |
| mitos2 | 13 | 13 | 0 | 0 | 0 | 13 | 922 | 1000 | 2142 |

## PCG (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom_native | 165 | 165 | 100.0% | 144 (87.3%) | 13 (95.2% cum.) | 8 | 0 | 0 | 0 | 2 |
| mitom_native_v27 | 165 | 165 | 100.0% | 142 (86.1%) | 13 (93.9% cum.) | 10 | 0 | 0 | 0 | 2 |
| mitos2 | 165 | 165 | 100.0% | 79 (47.9%) | 9 (53.3% cum.) | 77 | 0 | 0 | 0 | 5 |

## rRNA (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom_native | 26 | 26 | 100.0% | 15 (57.7%) | 2 (65.4% cum.) | 9 | 0 | 0 | 0 | 0 |
| mitom_native_v27 | 26 | 26 | 100.0% | 14 (53.8%) | 3 (65.4% cum.) | 9 | 0 | 0 | 0 | 0 |
| mitos2 | 26 | 24 | 92.3% | 1 (3.8%) | 5 (23.1% cum.) | 18 | 0 | 0 | 2 | 1 |

## tRNA (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom_native | 288 | 284 | 98.6% | 212 (73.6%) | 27 (83.0% cum.) | 5 | 0 | 4 | 0 | 1 |
| mitom_native_v27 | 288 | 281 | 97.6% | 205 (71.2%) | 30 (81.6% cum.) | 7 | 0 | 4 | 3 | 3 |
| mitos2 | 288 | 279 | 96.9% | 203 (70.5%) | 29 (80.6% cum.) | 8 | 0 | 4 | 5 | 1 |

## Same records for every tool (common_ok)

| tool | kind | truth | recall | both ends exact | within 3 nt |
|---|---|---|---|---|---|
| mitom_native | PCG | 165 | 100.0% | 87.3% | 95.2% |
| mitom_native_v27 | PCG | 165 | 100.0% | 86.1% | 93.9% |
| mitos2 | PCG | 165 | 100.0% | 47.9% | 53.3% |
| mitom_native | rRNA | 26 | 100.0% | 57.7% | 65.4% |
| mitom_native_v27 | rRNA | 26 | 100.0% | 53.8% | 65.4% |
| mitos2 | rRNA | 26 | 92.3% | 3.8% | 23.1% |
| mitom_native | tRNA | 288 | 98.6% | 73.6% | 83.0% |
| mitom_native_v27 | tRNA | 288 | 97.6% | 71.2% | 81.6% |
| mitos2 | tRNA | 288 | 96.9% | 70.5% | 80.6% |

## Protein genes, both ends exact, per phylum

| phylum | records | mitom_native | mitom_native_v27 | mitos2 |
|---|---|---|---|---|
| Arthropoda | 5 | 90.8% | 90.8% | 53.8% |
| Chordata | 3 | 73.7% | 73.7% | 63.2% |
| Echinodermata | 2 | 100.0% | 100.0% | 38.5% |
| Nematoda | 2 | 91.7% | 83.3% | 29.2% |
| Mollusca | 1 | 75.0% | 75.0% | 25.0% |

## Obvious errors by type

| error | mitom_native | mitom_native_v27 | mitos2 |
|---|---|---|---|
| boundary_gt30 | 3 | 4 | 39 |
| table2asn_error | 0 | 1 | 13 |
| name_error | 4 | 4 | 4 |
| missing_gene | 0 | 0 | 2 |
| internal_stop | 0 | 0 | 1 |
| duplicate_gene | 0 | 0 | 1 |
| records with any error | 3 | 5 | 13 |

## Obvious errors by type, records every tool finished (13 records)

| error | mitom_native | mitom_native_v27 | mitos2 |
|---|---|---|---|
| boundary_gt30 | 3 | 4 | 39 |
| table2asn_error | 0 | 1 | 13 |
| name_error | 4 | 4 | 4 |
| missing_gene | 0 | 0 | 2 |
| internal_stop | 0 | 0 | 1 |
| duplicate_gene | 0 | 0 | 1 |
| records with any error | 3 | 5 | 13 |

## Error records per tool

**mitom_native** (3 records)

- JN248424.2: boundary_gt30:rrnS
- NC_006886.2: boundary_gt30:cox3, boundary_gt30:rrnL
- NC_024511.2: name_error:L, name_error:S, name_error:S, name_error:L

**mitom_native_v27** (5 records)

- JN248424.2: boundary_gt30:rrnS
- NC_001328.1: boundary_gt30:cox1
- NC_002177.1: table2asn_error
- NC_006886.2: boundary_gt30:cox3, boundary_gt30:rrnL
- NC_024511.2: name_error:L, name_error:S, name_error:S, name_error:L

**mitos2** (13 records)

- JN248424.2: boundary_gt30:nad6, boundary_gt30:nad1, table2asn_error
- MK374364.1: boundary_gt30:nad3, boundary_gt30:nad6, table2asn_error
- NC_001322.1: boundary_gt30:cox2, boundary_gt30:nad5, boundary_gt30:nad1, table2asn_error
- NC_001327.1: boundary_gt30:nad6, boundary_gt30:nad1, boundary_gt30:nad2, boundary_gt30:cox3, boundary_gt30:cox1, boundary_gt30:nad3, boundary_gt30:nad5, boundary_gt30:rrnS, missing_gene:rrnL, internal_stop:nad3, table2asn_error
- NC_001328.1: boundary_gt30:cob, boundary_gt30:cox1, boundary_gt30:rrnS, boundary_gt30:rrnL, duplicate_gene:rrnL, table2asn_error
- NC_001453.1: boundary_gt30:cox2, boundary_gt30:atp6, boundary_gt30:cox3, boundary_gt30:nad5, table2asn_error
- NC_001572.1: boundary_gt30:cox2, boundary_gt30:atp6, boundary_gt30:cox3, boundary_gt30:nad4, boundary_gt30:nad5, table2asn_error
- NC_002081.1: table2asn_error
- NC_002177.1: boundary_gt30:nad4l, boundary_gt30:rrnS, boundary_gt30:rrnL, table2asn_error
- NC_005089.1: table2asn_error
- NC_006817.1: boundary_gt30:nad5, boundary_gt30:rrnL, table2asn_error
- NC_006886.2: boundary_gt30:cob, boundary_gt30:cox2, boundary_gt30:cox3, boundary_gt30:cox1, boundary_gt30:atp6, missing_gene:rrnL, table2asn_error
- NC_024511.2: boundary_gt30:rrnL, name_error:L, name_error:S, name_error:S, name_error:L, table2asn_error
