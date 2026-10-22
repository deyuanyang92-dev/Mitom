# Benchmark: mitom_native, mitom_lo, mitoz, mitos2, mitofinder (/media/deyuan/18T_B/yzz/mitom_work/bench_new/eval/all)

## Records

| tool | records | ok | failed | timeout | unsupported | with obvious error | median s/record | median peak RAM (MB) | max peak RAM (MB) |
|---|---|---|---|---|---|---|---|---|---|
| mitom_native | 6 | 6 | 0 | 0 | 0 | 0 | 93 | 194 | 236 |
| mitom_lo | 6 | 6 | 0 | 0 | 0 | 0 | 93 | 193 | 237 |
| mitoz | 6 | 6 | 0 | 0 | 0 | 5 | 261 | 556 | 720 |
| mitos2 | 6 | 6 | 0 | 0 | 0 | 6 | 455 | 1029 | 2091 |
| mitofinder | 6 | 6 | 0 | 0 | 0 | 6 | 261 | 250 | 262 |

## PCG (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom_native | 78 | 78 | 100.0% | 74 (94.9%) | 0 (94.9% cum.) | 4 | 0 | 0 | 0 | 0 |
| mitom_lo | 78 | 78 | 100.0% | 74 (94.9%) | 0 (94.9% cum.) | 4 | 0 | 0 | 0 | 0 |
| mitoz | 78 | 77 | 98.7% | 71 (91.0%) | 0 (91.0% cum.) | 6 | 0 | 0 | 1 | 0 |
| mitos2 | 78 | 78 | 100.0% | 50 (64.1%) | 0 (64.1% cum.) | 28 | 0 | 0 | 0 | 0 |
| mitofinder | 78 | 73 | 93.6% | 50 (64.1%) | 1 (65.4% cum.) | 22 | 0 | 0 | 5 | 0 |

## rRNA (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom_native | 12 | 12 | 100.0% | 4 (33.3%) | 3 (58.3% cum.) | 5 | 0 | 0 | 0 | 0 |
| mitom_lo | 12 | 12 | 100.0% | 4 (33.3%) | 3 (58.3% cum.) | 5 | 0 | 0 | 0 | 0 |
| mitoz | 12 | 12 | 100.0% | 3 (25.0%) | 3 (50.0% cum.) | 6 | 0 | 0 | 0 | 0 |
| mitos2 | 12 | 12 | 100.0% | 4 (33.3%) | 3 (58.3% cum.) | 5 | 0 | 0 | 0 | 0 |
| mitofinder | 12 | 10 | 83.3% | 1 (8.3%) | 0 (8.3% cum.) | 9 | 0 | 0 | 2 | 1 |

## tRNA (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom_native | 132 | 130 | 98.5% | 107 (81.1%) | 16 (93.2% cum.) | 0 | 0 | 0 | 2 | 4 |
| mitom_lo | 132 | 130 | 98.5% | 107 (81.1%) | 16 (93.2% cum.) | 0 | 0 | 0 | 2 | 4 |
| mitoz | 132 | 129 | 97.7% | 103 (78.0%) | 19 (92.4% cum.) | 0 | 0 | 0 | 3 | 4 |
| mitos2 | 132 | 131 | 99.2% | 107 (81.1%) | 17 (93.9% cum.) | 0 | 0 | 0 | 1 | 4 |
| mitofinder | 132 | 128 | 97.0% | 94 (71.2%) | 19 (85.6% cum.) | 0 | 0 | 2 | 2 | 4 |

## Same records for every tool (common_ok)

| tool | kind | truth | recall | both ends exact | within 3 nt |
|---|---|---|---|---|---|
| mitom_native | PCG | 78 | 100.0% | 94.9% | 94.9% |
| mitom_lo | PCG | 78 | 100.0% | 94.9% | 94.9% |
| mitoz | PCG | 78 | 98.7% | 91.0% | 91.0% |
| mitos2 | PCG | 78 | 100.0% | 64.1% | 64.1% |
| mitofinder | PCG | 78 | 93.6% | 64.1% | 65.4% |
| mitom_native | rRNA | 12 | 100.0% | 33.3% | 58.3% |
| mitom_lo | rRNA | 12 | 100.0% | 33.3% | 58.3% |
| mitoz | rRNA | 12 | 100.0% | 25.0% | 50.0% |
| mitos2 | rRNA | 12 | 100.0% | 33.3% | 58.3% |
| mitofinder | rRNA | 12 | 83.3% | 8.3% | 8.3% |
| mitom_native | tRNA | 132 | 98.5% | 81.1% | 93.2% |
| mitom_lo | tRNA | 132 | 98.5% | 81.1% | 93.2% |
| mitoz | tRNA | 132 | 97.7% | 78.0% | 92.4% |
| mitos2 | tRNA | 132 | 99.2% | 81.1% | 93.9% |
| mitofinder | tRNA | 132 | 97.0% | 71.2% | 85.6% |

## Protein genes, both ends exact, per phylum

| phylum | records | mitom_native | mitom_lo | mitoz | mitos2 | mitofinder |
|---|---|---|---|---|---|---|
| Echinodermata | 6 | 94.9% | 94.9% | 91.0% | 64.1% | 64.1% |

## Obvious errors by type

| error | mitom_native | mitom_lo | mitoz | mitos2 | mitofinder |
|---|---|---|---|---|---|
| boundary_gt30 | 0 | 0 | 2 | 19 | 21 |
| table2asn_error | 0 | 0 | 3 | 6 | 6 |
| missing_gene | 0 | 0 | 1 | 0 | 7 |
| name_error | 0 | 0 | 0 | 0 | 2 |
| records with any error | 0 | 0 | 5 | 6 | 6 |

## Obvious errors by type, records every tool finished (6 records)

| error | mitom_native | mitom_lo | mitoz | mitos2 | mitofinder |
|---|---|---|---|---|---|
| boundary_gt30 | 0 | 0 | 2 | 19 | 21 |
| table2asn_error | 0 | 0 | 3 | 6 | 6 |
| missing_gene | 0 | 0 | 1 | 0 | 7 |
| name_error | 0 | 0 | 0 | 0 | 2 |
| records with any error | 0 | 0 | 5 | 6 | 6 |

## Error records per tool

**mitom_native** (0 records)


**mitom_lo** (0 records)


**mitoz** (5 records)

- NC_007689.1: table2asn_error
- NC_045937.1: boundary_gt30:rrnL
- NC_066490.1: table2asn_error
- NC_068254.1: missing_gene:atp8, boundary_gt30:rrnL
- NC_069630.1: table2asn_error

**mitos2** (6 records)

- NC_007689.1: boundary_gt30:cox2, boundary_gt30:atp6, boundary_gt30:cox3, boundary_gt30:nad5, boundary_gt30:rrnL, table2asn_error
- NC_010692.1: boundary_gt30:cox2, table2asn_error
- NC_045937.1: boundary_gt30:cox2, boundary_gt30:atp6, boundary_gt30:cox3, boundary_gt30:nad5, boundary_gt30:nad2, boundary_gt30:nad1, table2asn_error
- NC_066490.1: boundary_gt30:nad1, boundary_gt30:cox2, boundary_gt30:atp6, table2asn_error
- NC_068254.1: boundary_gt30:nad4, boundary_gt30:cox2, table2asn_error
- NC_069630.1: boundary_gt30:cox2, boundary_gt30:atp6, table2asn_error

**mitofinder** (6 records)

- NC_007689.1: boundary_gt30:nad2, boundary_gt30:nad1, missing_gene:atp8, boundary_gt30:nad5, boundary_gt30:rrnS, boundary_gt30:rrnL, name_error:S, table2asn_error
- NC_010692.1: boundary_gt30:nad4, boundary_gt30:nad5, boundary_gt30:nad2, boundary_gt30:nad1, boundary_gt30:rrnS, boundary_gt30:rrnL, table2asn_error
- NC_045937.1: missing_gene:atp8, boundary_gt30:nad3, boundary_gt30:nad4, missing_gene:nad6, missing_gene:nad2, boundary_gt30:nad1, missing_gene:rrnL, missing_gene:rrnS, table2asn_error
- NC_066490.1: boundary_gt30:nad2, boundary_gt30:nad5, name_error:S, table2asn_error
- NC_068254.1: boundary_gt30:nad1, boundary_gt30:nad4l, table2asn_error
- NC_069630.1: missing_gene:atp8, boundary_gt30:nad5, boundary_gt30:nad1, boundary_gt30:rrnL, table2asn_error

