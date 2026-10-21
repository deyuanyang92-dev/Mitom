# Benchmark: mitom_v1, mitom_native (/media/deyuan/18T_B/yzz/mitom_work/bench_metazoa/eval/all)

## Records

| tool | records | ok | failed | timeout | unsupported | with obvious error | median s/record | median peak RAM (MB) | max peak RAM (MB) |
|---|---|---|---|---|---|---|---|---|---|
| mitom_v1 | 173 | 173 | 0 | 0 | 0 | 161 | 325 | - | - |
| mitom_native | 173 | 173 | 0 | 0 | 0 | 121 | 146 | 383 | 3197 |

## PCG (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom_v1 | 2178 | 2138 | 98.2% | 1453 (66.7%) | 138 (73.0% cum.) | 547 | 0 | 1 | 39 | 31 |
| mitom_native | 2178 | 2154 | 98.9% | 1720 (79.0%) | 134 (85.1% cum.) | 300 | 0 | 3 | 21 | 29 |

## rRNA (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom_v1 | 347 | 336 | 96.8% | 113 (32.6%) | 43 (45.0% cum.) | 180 | 0 | 0 | 11 | 8 |
| mitom_native | 347 | 327 | 94.2% | 112 (32.3%) | 52 (47.3% cum.) | 163 | 0 | 0 | 20 | 10 |

## tRNA (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom_v1 | 3335 | 3050 | 91.5% | 2179 (65.3%) | 367 (76.3% cum.) | 175 | 24 | 27 | 234 | 121 |
| mitom_native | 3335 | 3023 | 90.6% | 2166 (64.9%) | 369 (76.0% cum.) | 165 | 24 | 32 | 256 | 119 |

## Same records for every tool (common_ok)

| tool | kind | truth | recall | both ends exact | within 3 nt |
|---|---|---|---|---|---|
| mitom_v1 | PCG | 2178 | 98.2% | 66.7% | 73.0% |
| mitom_native | PCG | 2178 | 98.9% | 79.0% | 85.1% |
| mitom_v1 | rRNA | 347 | 96.8% | 32.6% | 45.0% |
| mitom_native | rRNA | 347 | 94.2% | 32.3% | 47.3% |
| mitom_v1 | tRNA | 3335 | 91.5% | 65.3% | 76.3% |
| mitom_native | tRNA | 3335 | 90.6% | 64.9% | 76.0% |

## Protein genes, both ends exact, per phylum

| phylum | records | mitom_v1 | mitom_native |
|---|---|---|---|
| Annelida | 12 | 82.7% | 78.8% |
| Arthropoda | 12 | 78.8% | 82.7% |
| Chordata | 12 | 84.0% | 83.3% |
| Cnidaria | 12 | 76.8% | 87.1% |
| Echinodermata | 12 | 87.2% | 86.5% |
| Mollusca | 12 | 76.5% | 76.5% |
| Nematoda | 12 | 66.4% | 79.9% |
| Platyhelminthes | 12 | 73.8% | 80.7% |
| Acanthocephala | 6 | 19.4% | 68.1% |
| Brachiopoda | 6 | 35.1% | 67.5% |
| Bryozoa | 6 | 68.8% | 77.9% |
| Hemichordata | 6 | 73.1% | 71.8% |
| Nemertea | 6 | 79.5% | 83.3% |
| Porifera | 6 | 93.6% | 96.2% |
| Xenacoelomorpha | 6 | 38.7% | 81.3% |
| Chaetognatha | 5 | 63.6% | 81.8% |
| Ctenophora | 4 | 16.7% | 45.2% |
| Onychophora | 4 | 30.8% | 69.2% |
| Kinorhyncha | 3 | 39.5% | 86.8% |
| Placozoa | 3 | 48.6% | 78.4% |
| Rotifera | 3 | 37.5% | 79.2% |
| Tardigrada | 3 | 56.4% | 76.9% |
| Entoprocta | 2 | 61.5% | 84.6% |
| Gnathostomulida | 2 | 26.9% | 80.8% |
| Nematomorpha | 2 | 38.5% | 23.1% |
| Gastrotricha | 1 | 46.2% | 100.0% |
| Orthonectida | 1 | 23.1% | 38.5% |
| Phoronida | 1 | 61.5% | 84.6% |
| Priapulida | 1 | 38.5% | 61.5% |

## Obvious errors by type

| error | mitom_v1 | mitom_native |
|---|---|---|
| boundary_gt30 | 384 | 208 |
| table2asn_error | 113 | 11 |
| missing_gene | 50 | 41 |
| name_error | 28 | 35 |
| strand_error | 24 | 25 |
| internal_stop | 14 | 6 |
| duplicate_gene | 1 | 0 |
| records with any error | 161 | 121 |

## Obvious errors by type, records every tool finished (173 records)

| error | mitom_v1 | mitom_native |
|---|---|---|
| boundary_gt30 | 384 | 208 |
| table2asn_error | 113 | 11 |
| missing_gene | 50 | 41 |
| name_error | 28 | 35 |
| strand_error | 24 | 25 |
| internal_stop | 14 | 6 |
| duplicate_gene | 1 | 0 |
| records with any error | 161 | 121 |

## Error records per tool

**mitom_v1** (161 records)

- NC_001627.1: strand_error:F, strand_error:N, table2asn_error
- NC_001861.1: boundary_gt30:nad4l, boundary_gt30:rrnS, table2asn_error
- NC_001887.1: boundary_gt30:cox1, internal_stop:cox1, table2asn_error
- NC_002082.1: table2asn_error
- NC_002322.1: boundary_gt30:nad5, boundary_gt30:nad3, boundary_gt30:nad2, boundary_gt30:R, table2asn_error
- NC_002545.1: table2asn_error
- NC_002681.1: boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_003086.1: boundary_gt30:cox3, boundary_gt30:nad1, boundary_gt30:nad4l, boundary_gt30:nad2
- NC_005928.1: boundary_gt30:nad5, table2asn_error
- NC_006083.1: boundary_gt30:cox2, boundary_gt30:nad3, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_006386.1: boundary_gt30:nad5, boundary_gt30:cox3, boundary_gt30:rrnL, boundary_gt30:rrnS, table2asn_error
- NC_006886.2: boundary_gt30:cox3, table2asn_error
- NC_006902.1: boundary_gt30:nad2, boundary_gt30:nad6, boundary_gt30:nad4l, boundary_gt30:rrnL
- NC_007438.1: name_error:K, table2asn_error
- NC_007933.1: table2asn_error
- NC_008047.1: missing_gene:rrnL, missing_gene:rrnS, name_error:Y
- NC_008151.3: boundary_gt30:nad6, missing_gene:cox1, boundary_gt30:nad2, missing_gene:cox1, missing_gene:rrnL, boundary_gt30:rrnL, boundary_gt30:rrnS, name_error:L, name_error:S, internal_stop:cox1, table2asn_error
- NC_008192.1: boundary_gt30:nad5, boundary_gt30:rrnL, boundary_gt30:rrnS, name_error:Y, table2asn_error
- NC_008556.1: boundary_gt30:nad6, strand_error:E
- NC_008557.1: boundary_gt30:nad6, boundary_gt30:nad4, boundary_gt30:nad5, table2asn_error
- NC_008640.1: missing_gene:nad3, missing_gene:nad3, boundary_gt30:rrnS, boundary_gt30:rrnL, table2asn_error
- NC_009082.1: boundary_gt30:nad2, boundary_gt30:rrnS, boundary_gt30:rrnL, name_error:M, name_error:C, name_error:E, internal_stop:nad3, table2asn_error
- NC_009732.1: table2asn_error
- NC_009938.1: table2asn_error
- NC_010209.1: boundary_gt30:rrnL, name_error:I
- NC_010432.1: boundary_gt30:nad2, boundary_gt30:cox3, internal_stop:cox3, table2asn_error
- NC_010472.1: boundary_gt30:nad2, boundary_gt30:atp6, boundary_gt30:rrnL
- NC_010484.1: boundary_gt30:nad4, table2asn_error
- NC_010961.1: boundary_gt30:rrnL, boundary_gt30:rrnS, name_error:V, name_error:T, name_error:S, name_error:H, table2asn_error
- NC_011638.1: table2asn_error
- NC_012821.1: boundary_gt30:nad4l, table2asn_error
- NC_013568.1: boundary_gt30:atp6, boundary_gt30:nad4, boundary_gt30:nad2, boundary_gt30:rrnL, name_error:L
- NC_013810.1: boundary_gt30:cox2, boundary_gt30:nad5, boundary_gt30:nad4l, boundary_gt30:rrnL, table2asn_error
- NC_013811.1: boundary_gt30:nad5, boundary_gt30:nad4l, boundary_gt30:cox1, boundary_gt30:cox2, boundary_gt30:rrnL, table2asn_error
- NC_013814.1: boundary_gt30:cox2, boundary_gt30:nad4l, boundary_gt30:rrnS, boundary_gt30:rrnL, table2asn_error
- NC_013831.1: boundary_gt30:nad5, boundary_gt30:rrnL, internal_stop:nad4, table2asn_error
- NC_013877.1: boundary_gt30:nad5, boundary_gt30:cob, strand_error:S, table2asn_error
- NC_014273.1: boundary_gt30:rrnL, boundary_gt30:rrnS, name_error:K, boundary_gt30:F, table2asn_error
- NC_014578.1: boundary_gt30:cox2, boundary_gt30:atp6, boundary_gt30:nad5, boundary_gt30:nad4, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_014848.1: boundary_gt30:nad4, boundary_gt30:rrnL, internal_stop:nad1, table2asn_error
- NC_014868.1: boundary_gt30:rrnL, boundary_gt30:rrnS, name_error:I
- NC_014869.1: boundary_gt30:nad4l, boundary_gt30:cox1, table2asn_error
- NC_015309.1: boundary_gt30:nad2, boundary_gt30:nad6, boundary_gt30:cox2, missing_gene:cox1, boundary_gt30:rrnS, boundary_gt30:rrnL, boundary_gt30:rrnL, name_error:W, internal_stop:cox1, duplicate_gene:rrnL, table2asn_error
- NC_015646.1: boundary_gt30:nad5, boundary_gt30:cox2, boundary_gt30:nad4l, boundary_gt30:cox3, boundary_gt30:rrnS, name_error:N
- NC_015829.1: missing_gene:rrnL, table2asn_error
- NC_015890.1: boundary_gt30:nad3, boundary_gt30:nad6, boundary_gt30:nad4l, boundary_gt30:rrnS, boundary_gt30:rrnL, boundary_gt30:V, boundary_gt30:Y, internal_stop:nad6, table2asn_error
- NC_016117.1: boundary_gt30:cox1, boundary_gt30:nad5, boundary_gt30:rrnL
- NC_016722.1: table2asn_error
- NC_016754.1: boundary_gt30:atp6, boundary_gt30:nad3, boundary_gt30:nad4, boundary_gt30:nad2, boundary_gt30:rrnS, name_error:T, name_error:N
- NC_017613.1: table2asn_error
- NC_017750.1: boundary_gt30:atp8, boundary_gt30:rrnS, boundary_gt30:rrnL, boundary_gt30:E, name_error:A
- NC_018355.1: name_error:V
- NC_018362.1: boundary_gt30:nad4, table2asn_error
- NC_018761.1: boundary_gt30:nad6, table2asn_error
- NC_020774.1: boundary_gt30:nad4l, boundary_gt30:rrnS, strand_error:H, table2asn_error
- NC_021422.1: boundary_gt30:nad4l, boundary_gt30:cox2, boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_022829.1: boundary_gt30:nad2, boundary_gt30:nad6, boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_022921.1: boundary_gt30:nad6, boundary_gt30:atp6, missing_gene:nad4l, boundary_gt30:nad4, boundary_gt30:nad5, boundary_gt30:cox3, boundary_gt30:nad2, boundary_gt30:rrnL, missing_gene:rrnS, internal_stop:nad3, table2asn_error
- NC_023523.1: table2asn_error
- NC_023772.1: table2asn_error
- NC_023929.1: table2asn_error
- NC_023931.1: table2asn_error
- NC_025335.1: boundary_gt30:rrnL
- NC_026285.1: table2asn_error
- NC_026516.1: boundary_gt30:nad6, missing_gene:atp6, boundary_gt30:nad3, boundary_gt30:nad4, boundary_gt30:nad5, boundary_gt30:nad1, boundary_gt30:nad2, boundary_gt30:rrnL, table2asn_error
- NC_026983.1: boundary_gt30:nad5, missing_gene:nad6, boundary_gt30:nad4l, boundary_gt30:nad4, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_026984.1: boundary_gt30:nad4, boundary_gt30:nad5, boundary_gt30:atp6, boundary_gt30:nad2, boundary_gt30:rrnS, missing_gene:rrnL, table2asn_error
- NC_026985.1: boundary_gt30:nad5, boundary_gt30:atp6, boundary_gt30:nad2, boundary_gt30:rrnS, boundary_gt30:rrnL, table2asn_error
- NC_027063.1: table2asn_error
- NC_027103.1: missing_gene:cox1, boundary_gt30:nad2, boundary_gt30:nad6, boundary_gt30:rrnL, internal_stop:nad5, table2asn_error
- NC_027109.1: table2asn_error
- NC_027697.1: boundary_gt30:nad6, boundary_gt30:rrnL, strand_error:L
- NC_028229.1: table2asn_error
- NC_028416.1: table2asn_error
- NC_029217.1: boundary_gt30:nad5, missing_gene:nad6
- NC_029219.1: missing_gene:nad6
- NC_029766.1: missing_gene:nad6, missing_gene:atp6, boundary_gt30:nad3, boundary_gt30:nad4l, boundary_gt30:nad4, boundary_gt30:nad5, boundary_gt30:nad2, boundary_gt30:rrnL, boundary_gt30:rrnS, table2asn_error
- NC_031213.1: boundary_gt30:nad5, boundary_gt30:nad4l, boundary_gt30:rrnL
- NC_031407.1: boundary_gt30:nad6, boundary_gt30:cob, boundary_gt30:rrnL, table2asn_error
- NC_031872.1: missing_gene:atp8, boundary_gt30:nad2, boundary_gt30:atp6, boundary_gt30:nad5, boundary_gt30:rrnL, table2asn_error
- NC_031873.1: boundary_gt30:nad5, boundary_gt30:rrnL
- NC_033868.1: boundary_gt30:nad1, boundary_gt30:atp6, boundary_gt30:rrnS, boundary_gt30:Q, boundary_gt30:D, table2asn_error
- NC_034715.1: boundary_gt30:rrnL, table2asn_error
- NC_034946.1: missing_gene:nad6, boundary_gt30:atp6, name_error:atp8, missing_gene:cox2, boundary_gt30:nad1, boundary_gt30:nad4, boundary_gt30:nad5, boundary_gt30:rrnL, table2asn_error
- NC_034948.1: boundary_gt30:nad1, boundary_gt30:nad5, boundary_gt30:nad2, boundary_gt30:nad6, boundary_gt30:cox2, missing_gene:atp6, boundary_gt30:nad4, boundary_gt30:cob, boundary_gt30:rrnL, boundary_gt30:rrnS, strand_error:N, strand_error:H
- NC_035061.1: boundary_gt30:rrnS, strand_error:Q, strand_error:M, strand_error:V, table2asn_error
- NC_035255.1: boundary_gt30:cob, boundary_gt30:nad1, boundary_gt30:cox3, boundary_gt30:nad6, boundary_gt30:nad2, boundary_gt30:nad4, boundary_gt30:rrnS, table2asn_error
- NC_035983.1: boundary_gt30:nad5, boundary_gt30:rrnL, boundary_gt30:rrnS, name_error:I
- NC_036051.1: boundary_gt30:nad2, boundary_gt30:cox2, table2asn_error
- NC_036308.1: missing_gene:cox3, boundary_gt30:atp6, boundary_gt30:nad5, boundary_gt30:nad6, boundary_gt30:rrnS, table2asn_error
- NC_036679.1: boundary_gt30:cox1, boundary_gt30:cox2, boundary_gt30:nad2, boundary_gt30:cob, boundary_gt30:atp6, boundary_gt30:nad6, boundary_gt30:nad3, boundary_gt30:nad5, boundary_gt30:cox3
- NC_037188.1: boundary_gt30:nad4l
- NC_038065.1: missing_gene:nad6, boundary_gt30:nad4, missing_gene:nad4l, boundary_gt30:nad5, boundary_gt30:cob, boundary_gt30:nad1, boundary_gt30:rrnS, boundary_gt30:rrnL, table2asn_error
- NC_038192.1: boundary_gt30:nad5, boundary_gt30:nad4l, table2asn_error
- NC_039685.1: boundary_gt30:cob, boundary_gt30:atp6, boundary_gt30:nad3, boundary_gt30:cox3, boundary_gt30:nad4, boundary_gt30:rrnL, boundary_gt30:R, table2asn_error
- NC_040107.1: boundary_gt30:cox1, boundary_gt30:cox2, missing_gene:atp8, strand_error:Y, table2asn_error
- NC_040108.1: boundary_gt30:nad4l, boundary_gt30:cox1, boundary_gt30:cox2, table2asn_error
- NC_040109.1: table2asn_error
- NC_040150.1: boundary_gt30:rrnL, name_error:L, table2asn_error
- NC_041549.1: boundary_gt30:nad6, missing_gene:nad1, boundary_gt30:nad2, boundary_gt30:nad4l, missing_gene:rrnL, boundary_gt30:rrnS, boundary_gt30:rrnL, internal_stop:cox1, table2asn_error
- NC_044092.1: boundary_gt30:cox3, boundary_gt30:cox1, boundary_gt30:rrnL
- NC_044095.1: missing_gene:atp8, boundary_gt30:atp6, boundary_gt30:nad4, boundary_gt30:nad5, missing_gene:nad6, table2asn_error
- NC_044669.1: boundary_gt30:nad4l, table2asn_error
- NC_045069.1: boundary_gt30:cox3, missing_gene:nad6, boundary_gt30:nad4, missing_gene:rrnS, missing_gene:rrnL
- NC_045231.1: table2asn_error
- NC_045305.1: boundary_gt30:cox3, missing_gene:nad4, missing_gene:nad4l, boundary_gt30:nad5, boundary_gt30:rrnS, boundary_gt30:rrnL, table2asn_error
- NC_045864.1: boundary_gt30:cox1, boundary_gt30:cox2, boundary_gt30:cox3, missing_gene:nad2, boundary_gt30:nad4, boundary_gt30:cob, boundary_gt30:nad5, missing_gene:nad6, missing_gene:nad4l, boundary_gt30:rrnL
- NC_046405.1: boundary_gt30:cox1, boundary_gt30:rrnS, boundary_gt30:rrnL, table2asn_error
- NC_046775.1: boundary_gt30:rrnL, table2asn_error
- NC_046851.1: boundary_gt30:nad4l, table2asn_error
- NC_048467.1: boundary_gt30:nad5, boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_049998.1: table2asn_error
- NC_053107.1: table2asn_error
- NC_054230.1: boundary_gt30:rrnL, table2asn_error
- NC_057981.1: table2asn_error
- NC_060482.1: boundary_gt30:atp6, missing_gene:nad4l, boundary_gt30:nad4, boundary_gt30:nad5, boundary_gt30:cox3, boundary_gt30:nad2, boundary_gt30:rrnL, boundary_gt30:rrnS, table2asn_error
- NC_060484.1: boundary_gt30:atp6, boundary_gt30:nad3, missing_gene:nad4l, boundary_gt30:nad4, boundary_gt30:nad5, boundary_gt30:nad1, boundary_gt30:cox2, boundary_gt30:cox3, boundary_gt30:rrnL, boundary_gt30:rrnS, table2asn_error
- NC_060553.1: boundary_gt30:nad4, missing_gene:nad6, missing_gene:atp8, missing_gene:nad4l, boundary_gt30:nad2, boundary_gt30:nad5, missing_gene:rrnL, boundary_gt30:rrnS, strand_error:F, table2asn_error
- NC_061205.1: table2asn_error
- NC_061323.1: boundary_gt30:atp8, boundary_gt30:rrnL, name_error:S
- NC_061414.1: table2asn_error
- NC_061910.1: boundary_gt30:nad4l, boundary_gt30:rrnS, strand_error:M, strand_error:Y, strand_error:C, strand_error:W, strand_error:Q, strand_error:G, strand_error:E, strand_error:T
- NC_062009.1: boundary_gt30:nad2, boundary_gt30:cox1, boundary_gt30:rrnL, boundary_gt30:rrnS, internal_stop:cox1, table2asn_error
- NC_062023.1: boundary_gt30:atp8, boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_062047.1: boundary_gt30:cox1, boundary_gt30:cox2, boundary_gt30:rrnL, boundary_gt30:rrnS, table2asn_error
- NC_062277.1: table2asn_error
- NC_063745.1: table2asn_error
- NC_063753.1: boundary_gt30:nad3, boundary_gt30:nad4, boundary_gt30:nad2, boundary_gt30:rrnS, table2asn_error
- NC_063769.1: boundary_gt30:nad4l, boundary_gt30:rrnS
- NC_063787.1: boundary_gt30:rrnL, table2asn_error
- NC_065088.1: boundary_gt30:rrnS
- NC_065650.1: table2asn_error
- NC_066979.1: boundary_gt30:nad4l, boundary_gt30:rrnS, table2asn_error
- NC_067797.1: table2asn_error
- NC_067900.1: boundary_gt30:nad1, boundary_gt30:rrnS, table2asn_error
- NC_067930.1: boundary_gt30:rrnS, table2asn_error
- NC_068082.1: strand_error:T, strand_error:P, table2asn_error
- NC_068797.1: boundary_gt30:rrnL, table2asn_error
- NC_069055.1: boundary_gt30:nad5, boundary_gt30:rrnL, boundary_gt30:rrnS, internal_stop:nad5, table2asn_error
- NC_069205.1: table2asn_error
- NC_069206.1: boundary_gt30:cob, table2asn_error
- NC_069241.1: boundary_gt30:cox1, boundary_gt30:nad4l, table2asn_error
- NC_069297.1: boundary_gt30:rrnS, name_error:S
- NC_072223.1: boundary_gt30:rrnL
- NC_072228.1: boundary_gt30:nad4l
- NC_073494.1: boundary_gt30:nad2, boundary_gt30:rrnS, boundary_gt30:rrnL, strand_error:M
- NC_080261.1: boundary_gt30:nad1, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_080948.1: boundary_gt30:rrnL
- NC_080956.1: boundary_gt30:nad4l, boundary_gt30:rrnL
- NC_082192.1: table2asn_error
- NC_083035.1: table2asn_error
- NC_083122.1: boundary_gt30:cox1, boundary_gt30:nad6, boundary_gt30:atp6, boundary_gt30:cox2, boundary_gt30:cox3, boundary_gt30:nad2, boundary_gt30:rrnS, table2asn_error
- NC_083273.1: boundary_gt30:rrnL, table2asn_error
- NC_084457.1: boundary_gt30:nad6, boundary_gt30:nad5
- NC_085462.1: boundary_gt30:cox3, boundary_gt30:nad1, boundary_gt30:cox2, boundary_gt30:nad2, boundary_gt30:rrnL
- NC_085563.1: boundary_gt30:cox1, boundary_gt30:cox3, boundary_gt30:nad2, boundary_gt30:atp6, missing_gene:nad6, boundary_gt30:nad5, boundary_gt30:rrnS, table2asn_error
- NC_085569.1: missing_gene:atp8, missing_gene:nad4l, boundary_gt30:nad4, boundary_gt30:nad5, missing_gene:nad6, missing_gene:rrnL, name_error:A, table2asn_error
- NC_085811.1: boundary_gt30:rrnL, table2asn_error
- NC_085825.1: boundary_gt30:nad2, boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_086642.1: boundary_gt30:nad5, internal_stop:nad5, table2asn_error
- NC_088493.1: boundary_gt30:nad6, table2asn_error

**mitom_native** (121 records)

- NC_001627.1: boundary_gt30:rrnL, strand_error:F, strand_error:N
- NC_001861.1: boundary_gt30:rrnS
- NC_001887.1: boundary_gt30:cox1, internal_stop:cox1, table2asn_error
- NC_005928.1: boundary_gt30:rrnL
- NC_006083.1: boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_006386.1: boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_006886.2: boundary_gt30:cox3, boundary_gt30:rrnL, table2asn_error
- NC_006902.1: boundary_gt30:rrnL
- NC_007438.1: name_error:K
- NC_008047.1: missing_gene:rrnL, missing_gene:rrnS, name_error:Y
- NC_008151.3: missing_gene:cox1, missing_gene:cox1, missing_gene:rrnL, boundary_gt30:rrnL, boundary_gt30:rrnS, name_error:L, name_error:S
- NC_008192.1: boundary_gt30:rrnL, boundary_gt30:rrnS, name_error:Y
- NC_008556.1: strand_error:E
- NC_008557.1: boundary_gt30:nad4, boundary_gt30:nad5
- NC_008640.1: missing_gene:nad3, missing_gene:nad3, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_009082.1: boundary_gt30:nad3, boundary_gt30:rrnS, boundary_gt30:rrnL, name_error:M, name_error:C, name_error:E
- NC_010209.1: boundary_gt30:rrnL, name_error:I
- NC_010472.1: boundary_gt30:rrnL
- NC_010961.1: boundary_gt30:rrnL, boundary_gt30:rrnS, name_error:V, name_error:T, name_error:H, table2asn_error
- NC_013568.1: missing_gene:cox3, boundary_gt30:rrnL, name_error:L
- NC_013810.1: boundary_gt30:rrnL
- NC_013811.1: boundary_gt30:cox1, boundary_gt30:rrnL
- NC_013814.1: boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_013877.1: boundary_gt30:cob, strand_error:S
- NC_014273.1: boundary_gt30:rrnL, boundary_gt30:rrnS, name_error:K, boundary_gt30:F, name_error:T
- NC_014578.1: boundary_gt30:rrnS, name_error:Q
- NC_014848.1: boundary_gt30:rrnL
- NC_014868.1: boundary_gt30:rrnL, boundary_gt30:rrnS, name_error:I
- NC_015309.1: missing_gene:cox1, boundary_gt30:rrnS, boundary_gt30:rrnL, missing_gene:rrnL, name_error:W
- NC_015646.1: boundary_gt30:nad5, boundary_gt30:cox2, boundary_gt30:rrnS, name_error:N
- NC_015829.1: missing_gene:rrnL
- NC_015890.1: missing_gene:cox3, boundary_gt30:rrnL, boundary_gt30:Y
- NC_016117.1: missing_gene:rrnS, missing_gene:rrnL
- NC_016722.1: boundary_gt30:cox2
- NC_016754.1: boundary_gt30:rrnL, boundary_gt30:rrnS, name_error:T, name_error:N
- NC_017750.1: boundary_gt30:E, name_error:A
- NC_018355.1: name_error:V
- NC_018362.1: boundary_gt30:nad4
- NC_020774.1: strand_error:H
- NC_021422.1: boundary_gt30:cox2, boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_022829.1: boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_022921.1: boundary_gt30:nad4, boundary_gt30:nad1, boundary_gt30:cox3, boundary_gt30:rrnL, missing_gene:rrnS
- NC_023523.1: table2asn_error
- NC_025335.1: boundary_gt30:rrnL
- NC_026516.1: boundary_gt30:rrnL
- NC_026983.1: missing_gene:nad5, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_026984.1: boundary_gt30:nad2, boundary_gt30:rrnS, missing_gene:rrnL
- NC_026985.1: boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_027103.1: missing_gene:nad5, boundary_gt30:rrnL
- NC_027697.1: boundary_gt30:cox1, strand_error:L
- NC_028229.1: table2asn_error
- NC_028416.1: boundary_gt30:rrnL
- NC_029766.1: boundary_gt30:rrnL
- NC_031213.1: boundary_gt30:nad5, boundary_gt30:rrnL
- NC_031407.1: boundary_gt30:rrnL
- NC_031872.1: boundary_gt30:rrnL
- NC_031873.1: boundary_gt30:rrnL
- NC_033868.1: boundary_gt30:nad1, boundary_gt30:atp6, boundary_gt30:rrnS, boundary_gt30:Q, boundary_gt30:D
- NC_034715.1: boundary_gt30:rrnL
- NC_034946.1: boundary_gt30:nad6, boundary_gt30:nad1, boundary_gt30:rrnL
- NC_034948.1: boundary_gt30:cob, boundary_gt30:rrnL, strand_error:N, strand_error:H
- NC_035061.1: strand_error:Q, strand_error:M, name_error:Y, strand_error:V
- NC_035255.1: boundary_gt30:cob, boundary_gt30:nad3, boundary_gt30:cox3, boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_035983.1: boundary_gt30:rrnS, name_error:I
- NC_036308.1: boundary_gt30:cob, name_error:cox3
- NC_036679.1: name_error:Q
- NC_038065.1: missing_gene:nad6, missing_gene:nad3, boundary_gt30:nad4l, boundary_gt30:nad5, boundary_gt30:cob, boundary_gt30:nad1, missing_gene:rrnS, boundary_gt30:rrnL, internal_stop:atp6, internal_stop:nad6, internal_stop:cob, table2asn_error
- NC_038192.1: boundary_gt30:nad5, boundary_gt30:rrnL
- NC_039685.1: boundary_gt30:atp6, boundary_gt30:rrnL
- NC_040107.1: boundary_gt30:cox1, boundary_gt30:cox2, strand_error:Y
- NC_040108.1: boundary_gt30:nad5, boundary_gt30:cox1, boundary_gt30:cox2, name_error:Y
- NC_040109.1: boundary_gt30:nad5
- NC_040150.1: boundary_gt30:rrnL
- NC_041549.1: boundary_gt30:nad6, boundary_gt30:nad5, boundary_gt30:nad2, boundary_gt30:intron, strand_error:intron, missing_gene:rrnL, boundary_gt30:rrnS, missing_gene:rrnL
- NC_044092.1: boundary_gt30:cox1, boundary_gt30:rrnL
- NC_044095.1: missing_gene:atp8, boundary_gt30:nad5, missing_gene:nad2
- NC_044109.1: boundary_gt30:atp8
- NC_045069.1: boundary_gt30:nad5, boundary_gt30:rrnS, missing_gene:rrnL, internal_stop:nad5, table2asn_error
- NC_045305.1: name_error:nad4l, boundary_gt30:nad5, missing_gene:rrnS, missing_gene:rrnL
- NC_045864.1: name_error:nad2, boundary_gt30:cob, boundary_gt30:nad5, missing_gene:nad6, missing_gene:nad4l, missing_gene:rrnL, missing_gene:rrnS
- NC_046405.1: boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_046775.1: boundary_gt30:rrnL, table2asn_error
- NC_046790.1: boundary_gt30:nad2
- NC_048467.1: boundary_gt30:nad4l, boundary_gt30:cox2, boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_054230.1: boundary_gt30:rrnL
- NC_060482.1: boundary_gt30:rrnL, boundary_gt30:rrnS
- NC_060484.1: boundary_gt30:nad6, boundary_gt30:atp6, boundary_gt30:nad3, boundary_gt30:nad4, boundary_gt30:rrnL, boundary_gt30:rrnS, table2asn_error
- NC_060553.1: boundary_gt30:nad4, missing_gene:nad6, missing_gene:atp8, boundary_gt30:nad2, missing_gene:rrnL, boundary_gt30:rrnS, strand_error:F
- NC_061323.1: boundary_gt30:rrnL, name_error:S
- NC_061910.1: boundary_gt30:nad4l, boundary_gt30:rrnS, strand_error:M, strand_error:Y, strand_error:C, strand_error:W, strand_error:Q, strand_error:G, strand_error:E, strand_error:T
- NC_062009.1: boundary_gt30:nad2, boundary_gt30:nad3, boundary_gt30:rrnS
- NC_062023.1: boundary_gt30:atp8, boundary_gt30:rrnS
- NC_062047.1: boundary_gt30:nad6, boundary_gt30:cox1, name_error:Q
- NC_063753.1: boundary_gt30:nad3, boundary_gt30:nad4, boundary_gt30:rrnS
- NC_063769.1: boundary_gt30:rrnS
- NC_065088.1: boundary_gt30:rrnS
- NC_066979.1: boundary_gt30:nad5
- NC_067900.1: boundary_gt30:rrnS
- NC_067930.1: boundary_gt30:rrnS
- NC_068082.1: strand_error:T, strand_error:P
- NC_068797.1: boundary_gt30:rrnL
- NC_069055.1: boundary_gt30:rrnL, boundary_gt30:rrnS, internal_stop:nad4l, table2asn_error
- NC_069205.1: boundary_gt30:nad2
- NC_069206.1: boundary_gt30:nad1, boundary_gt30:cob
- NC_069241.1: boundary_gt30:nad4l
- NC_069297.1: boundary_gt30:rrnS, name_error:S
- NC_072223.1: boundary_gt30:nad2, boundary_gt30:rrnL
- NC_073494.1: boundary_gt30:nad2, boundary_gt30:rrnS, boundary_gt30:rrnL, strand_error:M
- NC_080261.1: boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_080948.1: boundary_gt30:rrnL
- NC_080956.1: boundary_gt30:rrnL
- NC_083122.1: boundary_gt30:cox1, boundary_gt30:cox2, boundary_gt30:nad4, boundary_gt30:rrnL, missing_gene:rrnS
- NC_083273.1: boundary_gt30:rrnL
- NC_084457.1: boundary_gt30:nad6, missing_gene:atp8
- NC_085462.1: boundary_gt30:cox3, boundary_gt30:nad1, boundary_gt30:cox2, missing_gene:nad4, boundary_gt30:rrnL
- NC_085563.1: boundary_gt30:nad4, boundary_gt30:nad2, boundary_gt30:nad4l, boundary_gt30:nad6, boundary_gt30:nad3, boundary_gt30:rrnS
- NC_085569.1: missing_gene:atp8, boundary_gt30:atp6, boundary_gt30:nad3, boundary_gt30:nad4, boundary_gt30:nad5, missing_gene:rrnL, name_error:A
- NC_085811.1: boundary_gt30:nad5, table2asn_error
- NC_085825.1: boundary_gt30:nad2, boundary_gt30:nad5, boundary_gt30:rrnL, boundary_gt30:rrnS, name_error:Q
- NC_086642.1: boundary_gt30:nad4l, missing_gene:atp8
- NC_088493.1: boundary_gt30:nad6
