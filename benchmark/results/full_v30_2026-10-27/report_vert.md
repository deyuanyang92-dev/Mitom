# Benchmark: mitom_native, mitom_native_v28, mitos2 (/media/deyuan/18T_B/yzz/mitom_work/eval_v30/vert/eval/all)

## Records

| tool | records | ok | failed | timeout | unsupported | with obvious error | median s/record | median peak RAM (MB) | max peak RAM (MB) |
|---|---|---|---|---|---|---|---|---|---|
| mitom_native | 188 | 188 | 0 | 0 | 0 | 46 | 77 | 241 | 326 |
| mitom_native_v28 | 188 | 188 | 0 | 0 | 0 | 68 | 103 | 240 | 327 |
| mitos2 | 188 | 188 | 0 | 0 | 0 | 188 | 447 | 831 | 2094 |

## PCG (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom_native | 2444 | 2443 | 100.0% | 2172 (88.9%) | 206 (97.3% cum.) | 65 | 0 | 0 | 1 | 1 |
| mitom_native_v28 | 2444 | 2443 | 100.0% | 2146 (87.8%) | 219 (96.8% cum.) | 78 | 0 | 0 | 1 | 1 |
| mitos2 | 2444 | 2444 | 100.0% | 1778 (72.7%) | 433 (90.5% cum.) | 233 | 0 | 0 | 0 | 55 |

## rRNA (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom_native | 376 | 376 | 100.0% | 254 (67.6%) | 85 (90.2% cum.) | 37 | 0 | 0 | 0 | 0 |
| mitom_native_v28 | 376 | 376 | 100.0% | 254 (67.6%) | 85 (90.2% cum.) | 37 | 0 | 0 | 0 | 0 |
| mitos2 | 376 | 376 | 100.0% | 139 (37.0%) | 92 (61.4% cum.) | 145 | 0 | 0 | 0 | 3 |

## tRNA (records where every tool finished: stratum `common_ok` below; this table: all)

| tool | truth | found | recall | both ends exact | within 3 nt | >3 nt | strand err | name err | missed | extra |
|---|---|---|---|---|---|---|---|---|---|---|
| mitom_native | 4140 | 4073 | 98.4% | 3453 (83.4%) | 181 (87.8% cum.) | 23 | 31 | 20 | 16 | 9 |
| mitom_native_v28 | 4140 | 4075 | 98.4% | 3455 (83.5%) | 181 (87.8% cum.) | 23 | 31 | 20 | 14 | 9 |
| mitos2 | 4140 | 4070 | 98.3% | 3452 (83.4%) | 180 (87.7% cum.) | 23 | 31 | 27 | 12 | 15 |

## Same records for every tool (common_ok)

| tool | kind | truth | recall | both ends exact | within 3 nt |
|---|---|---|---|---|---|
| mitom_native | PCG | 2444 | 100.0% | 88.9% | 97.3% |
| mitom_native_v28 | PCG | 2444 | 100.0% | 87.8% | 96.8% |
| mitos2 | PCG | 2444 | 100.0% | 72.7% | 90.5% |
| mitom_native | rRNA | 376 | 100.0% | 67.6% | 90.2% |
| mitom_native_v28 | rRNA | 376 | 100.0% | 67.6% | 90.2% |
| mitos2 | rRNA | 376 | 100.0% | 37.0% | 61.4% |
| mitom_native | tRNA | 4140 | 98.4% | 83.4% | 87.8% |
| mitom_native_v28 | tRNA | 4140 | 98.4% | 83.5% | 87.8% |
| mitos2 | tRNA | 4140 | 98.3% | 83.4% | 87.7% |

## Protein genes, both ends exact, per phylum

| phylum | records | mitom_native | mitom_native_v28 | mitos2 |
|---|---|---|---|---|
| Chordata | 188 | 88.9% | 87.8% | 72.7% |

## Obvious errors by type

| error | mitom_native | mitom_native_v28 | mitos2 |
|---|---|---|---|
| table2asn_error | 0 | 31 | 188 |
| boundary_gt30 | 27 | 22 | 89 |
| strand_error | 31 | 31 | 31 |
| name_error | 20 | 20 | 27 |
| internal_stop | 0 | 21 | 7 |
| duplicate_gene | 0 | 0 | 20 |
| missing_gene | 1 | 1 | 0 |
| records with any error | 46 | 68 | 188 |

## Obvious errors by type, records every tool finished (188 records)

| error | mitom_native | mitom_native_v28 | mitos2 |
|---|---|---|---|
| table2asn_error | 0 | 31 | 188 |
| boundary_gt30 | 27 | 22 | 89 |
| strand_error | 31 | 31 | 31 |
| name_error | 20 | 20 | 27 |
| internal_stop | 0 | 21 | 7 |
| duplicate_gene | 0 | 0 | 20 |
| missing_gene | 1 | 1 | 0 |
| records with any error | 46 | 68 | 188 |

## Error records per tool

**mitom_native** (46 records)

- NC_039356.1: boundary_gt30:rrnL
- NC_039438.1: strand_error:Q
- NC_039669.1: strand_error:E, strand_error:P
- NC_039950.1: boundary_gt30:rrnL
- NC_041673.1: boundary_gt30:atp8
- NC_042750.1: strand_error:S
- NC_043950.1: strand_error:E
- NC_044877.1: strand_error:N
- NC_044884.1: strand_error:N
- NC_045217.1: strand_error:I, strand_error:Q
- NC_045916.1: strand_error:S
- NC_047241.1: strand_error:Q
- NC_050988.1: strand_error:I, strand_error:Q
- NC_051473.1: boundary_gt30:rrnS
- NC_052782.1: boundary_gt30:nad3
- NC_053822.1: strand_error:Y, strand_error:S
- NC_053898.1: boundary_gt30:rrnS
- NC_054268.1: boundary_gt30:cox2
- NC_054291.1: name_error:L, name_error:S, name_error:S, name_error:L
- NC_056403.1: name_error:L, name_error:S, name_error:S, name_error:L
- NC_056807.1: boundary_gt30:rrnL
- NC_057282.1: boundary_gt30:atp8
- NC_057518.1: boundary_gt30:cox2, strand_error:Y
- NC_057532.1: strand_error:Q
- NC_058991.1: boundary_gt30:rrnL
- NC_060993.1: boundary_gt30:rrnL
- NC_062636.1: boundary_gt30:rrnL
- NC_063494.1: boundary_gt30:atp8
- NC_065055.1: strand_error:Q, strand_error:A, strand_error:N, strand_error:C, strand_error:Y, strand_error:S
- NC_065150.1: boundary_gt30:rrnL
- NC_066464.1: strand_error:Y
- NC_069308.1: boundary_gt30:atp6
- NC_071234.1: missing_gene:nad6
- NC_071873.1: boundary_gt30:atp8, boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_071908.1: boundary_gt30:rrnL
- NC_072211.1: name_error:L, name_error:S, name_error:S, name_error:L
- NC_077264.1: strand_error:S, strand_error:Y, strand_error:C, strand_error:N, strand_error:A, strand_error:Q
- NC_079638.1: name_error:E, name_error:D, name_error:N, name_error:Q
- NC_083010.1: boundary_gt30:atp8
- NC_083073.1: boundary_gt30:cob
- NC_083983.1: boundary_gt30:rrnL, strand_error:E
- NC_084075.1: name_error:L, name_error:S, name_error:S, name_error:L
- NC_086643.1: boundary_gt30:rrnL
- NC_086960.1: boundary_gt30:rrnL
- NC_087724.1: boundary_gt30:nad3
- NC_087991.1: boundary_gt30:rrnL

**mitom_native_v28** (68 records)

- NC_038220.1: internal_stop:nad3, table2asn_error
- NC_039356.1: boundary_gt30:rrnL
- NC_039438.1: strand_error:Q
- NC_039669.1: strand_error:E, strand_error:P
- NC_039950.1: boundary_gt30:rrnL, internal_stop:nad3, table2asn_error
- NC_040004.1: internal_stop:nad3, table2asn_error
- NC_040005.1: internal_stop:nad3, table2asn_error
- NC_042750.1: strand_error:S
- NC_043950.1: strand_error:E, table2asn_error
- NC_044672.1: internal_stop:nad3, table2asn_error
- NC_044673.1: internal_stop:nad3, table2asn_error
- NC_044877.1: strand_error:N
- NC_044884.1: strand_error:N
- NC_045217.1: strand_error:I, strand_error:Q
- NC_045377.1: table2asn_error
- NC_045916.1: strand_error:S
- NC_047241.1: strand_error:Q
- NC_048989.1: internal_stop:nad3, table2asn_error
- NC_050835.1: table2asn_error
- NC_050988.1: strand_error:I, strand_error:Q
- NC_051473.1: boundary_gt30:rrnS, internal_stop:nad3, table2asn_error
- NC_052782.1: boundary_gt30:nad3
- NC_052784.1: table2asn_error
- NC_052798.1: internal_stop:nad3, table2asn_error
- NC_052799.1: internal_stop:nad3, table2asn_error
- NC_052805.1: internal_stop:nad3, table2asn_error
- NC_052809.1: internal_stop:nad3, table2asn_error
- NC_052811.1: internal_stop:nad3, table2asn_error
- NC_052818.1: internal_stop:nad3, table2asn_error
- NC_052822.1: internal_stop:nad3, table2asn_error
- NC_053822.1: strand_error:Y, strand_error:S
- NC_053843.1: internal_stop:nad3, table2asn_error
- NC_053898.1: boundary_gt30:rrnS
- NC_054268.1: boundary_gt30:cox2
- NC_054291.1: name_error:L, name_error:S, name_error:S, name_error:L
- NC_056091.1: table2asn_error
- NC_056403.1: name_error:L, name_error:S, name_error:S, name_error:L
- NC_056807.1: boundary_gt30:rrnL
- NC_057518.1: boundary_gt30:cox2, strand_error:Y
- NC_057532.1: strand_error:Q
- NC_058991.1: boundary_gt30:rrnL
- NC_060993.1: boundary_gt30:rrnL
- NC_062636.1: boundary_gt30:rrnL
- NC_062960.1: internal_stop:nad3, table2asn_error
- NC_065055.1: strand_error:Q, strand_error:A, strand_error:N, strand_error:C, strand_error:Y, strand_error:S
- NC_065150.1: boundary_gt30:rrnL
- NC_066464.1: strand_error:Y
- NC_066800.1: table2asn_error
- NC_068663.1: table2asn_error
- NC_068903.1: internal_stop:nad3, table2asn_error
- NC_069308.1: boundary_gt30:atp6
- NC_071234.1: missing_gene:nad6, table2asn_error
- NC_071873.1: boundary_gt30:rrnS, boundary_gt30:rrnL
- NC_071908.1: boundary_gt30:rrnL
- NC_072211.1: name_error:L, name_error:S, name_error:S, name_error:L
- NC_072569.1: internal_stop:nad3, table2asn_error
- NC_077264.1: strand_error:S, strand_error:Y, strand_error:C, strand_error:N, strand_error:A, strand_error:Q, table2asn_error
- NC_079638.1: name_error:E, name_error:D, name_error:N, name_error:Q
- NC_082168.1: internal_stop:nad3, table2asn_error
- NC_083073.1: boundary_gt30:cob
- NC_083983.1: boundary_gt30:rrnL, strand_error:E
- NC_084075.1: name_error:L, name_error:S, name_error:S, name_error:L
- NC_086643.1: boundary_gt30:rrnL
- NC_086816.1: table2asn_error
- NC_086960.1: boundary_gt30:rrnL
- NC_087724.1: boundary_gt30:nad3
- NC_087991.1: boundary_gt30:rrnL
- NC_088448.1: internal_stop:nad3, table2asn_error

**mitos2** (188 records)

- NC_038175.1: boundary_gt30:cox1, boundary_gt30:cox2, table2asn_error
- NC_038220.1: boundary_gt30:nad3, table2asn_error
- NC_039135.1: table2asn_error
- NC_039356.1: boundary_gt30:rrnL, table2asn_error
- NC_039438.1: boundary_gt30:rrnL, strand_error:Q, table2asn_error
- NC_039669.1: boundary_gt30:nad4, strand_error:E, strand_error:P, table2asn_error
- NC_039717.1: name_error:D, table2asn_error
- NC_039822.1: table2asn_error
- NC_039950.1: boundary_gt30:rrnL, table2asn_error
- NC_039987.1: boundary_gt30:nad2, table2asn_error
- NC_040004.1: boundary_gt30:nad3, table2asn_error
- NC_040005.1: boundary_gt30:nad3, table2asn_error
- NC_041284.1: boundary_gt30:nad4l, boundary_gt30:cob, boundary_gt30:rrnL, table2asn_error
- NC_041673.1: table2asn_error
- NC_042393.1: boundary_gt30:rrnL, table2asn_error
- NC_042750.1: strand_error:S, table2asn_error
- NC_042753.1: duplicate_gene:nad5, table2asn_error
- NC_042754.1: table2asn_error
- NC_043830.1: internal_stop:nad4, duplicate_gene:nad4, table2asn_error
- NC_043853.1: table2asn_error
- NC_043875.1: boundary_gt30:cox2, table2asn_error
- NC_043950.1: boundary_gt30:nad3, strand_error:E, table2asn_error
- NC_044427.1: name_error:L, table2asn_error
- NC_044459.1: table2asn_error
- NC_044672.1: boundary_gt30:nad3, duplicate_gene:rrnS, table2asn_error
- NC_044673.1: boundary_gt30:nad3, duplicate_gene:rrnS, table2asn_error
- NC_044716.1: table2asn_error
- NC_044745.1: table2asn_error
- NC_044769.1: table2asn_error
- NC_044858.1: table2asn_error
- NC_044877.1: strand_error:N, table2asn_error
- NC_044884.1: strand_error:N, table2asn_error
- NC_045213.1: boundary_gt30:rrnL, table2asn_error
- NC_045217.1: strand_error:I, strand_error:Q, name_error:D, table2asn_error
- NC_045377.1: boundary_gt30:nad3, table2asn_error
- NC_045916.1: strand_error:S, table2asn_error
- NC_046684.1: table2asn_error
- NC_047231.1: table2asn_error
- NC_047241.1: boundary_gt30:rrnL, strand_error:Q, table2asn_error
- NC_047480.1: table2asn_error
- NC_048989.1: boundary_gt30:nad3, duplicate_gene:nad3, table2asn_error
- NC_049119.1: table2asn_error
- NC_049866.1: duplicate_gene:nad4, table2asn_error
- NC_050278.1: table2asn_error
- NC_050835.1: boundary_gt30:nad3, table2asn_error
- NC_050988.1: strand_error:I, strand_error:Q, name_error:D, table2asn_error
- NC_050994.1: table2asn_error
- NC_051473.1: boundary_gt30:nad3, table2asn_error
- NC_051875.1: table2asn_error
- NC_052782.1: boundary_gt30:nad3, internal_stop:nad6, duplicate_gene:nad6, table2asn_error
- NC_052784.1: boundary_gt30:nad3, boundary_gt30:rrnL, table2asn_error
- NC_052796.1: boundary_gt30:rrnL, table2asn_error
- NC_052798.1: boundary_gt30:nad3, duplicate_gene:nad3, table2asn_error
- NC_052799.1: boundary_gt30:rrnL, table2asn_error
- NC_052805.1: boundary_gt30:nad3, table2asn_error
- NC_052809.1: boundary_gt30:nad3, internal_stop:atp6, table2asn_error
- NC_052811.1: boundary_gt30:nad3, table2asn_error
- NC_052818.1: boundary_gt30:nad3, table2asn_error
- NC_052821.1: internal_stop:atp6, table2asn_error
- NC_052822.1: boundary_gt30:nad3, table2asn_error
- NC_053556.1: table2asn_error
- NC_053637.1: boundary_gt30:rrnL, table2asn_error
- NC_053752.1: table2asn_error
- NC_053758.1: table2asn_error
- NC_053822.1: strand_error:Y, strand_error:S, table2asn_error
- NC_053843.1: boundary_gt30:nad3, table2asn_error
- NC_053898.1: boundary_gt30:cox1, boundary_gt30:rrnS, boundary_gt30:rrnL, name_error:D, duplicate_gene:rrnL, table2asn_error
- NC_053978.1: boundary_gt30:cob, table2asn_error
- NC_054244.1: table2asn_error
- NC_054268.1: name_error:D, table2asn_error
- NC_054291.1: name_error:L, name_error:S, name_error:S, name_error:L, table2asn_error
- NC_056091.1: boundary_gt30:nad3, table2asn_error
- NC_056170.1: table2asn_error
- NC_056291.1: table2asn_error
- NC_056308.1: table2asn_error
- NC_056334.1: table2asn_error
- NC_056403.1: boundary_gt30:rrnL, name_error:L, name_error:S, name_error:S, name_error:L, table2asn_error
- NC_056807.1: table2asn_error
- NC_057224.1: table2asn_error
- NC_057232.1: boundary_gt30:rrnL, table2asn_error
- NC_057260.1: boundary_gt30:atp6, boundary_gt30:rrnL, table2asn_error
- NC_057282.1: boundary_gt30:rrnL, table2asn_error
- NC_057518.1: strand_error:Y, name_error:D, table2asn_error
- NC_057532.1: strand_error:Q, table2asn_error
- NC_057650.1: table2asn_error
- NC_058991.1: table2asn_error
- NC_060993.1: table2asn_error
- NC_062636.1: table2asn_error
- NC_062945.1: table2asn_error
- NC_062960.1: boundary_gt30:nad3, table2asn_error
- NC_063128.1: table2asn_error
- NC_063494.1: table2asn_error
- NC_063691.1: boundary_gt30:rrnL, table2asn_error
- NC_063704.1: table2asn_error
- NC_063714.1: table2asn_error
- NC_063867.1: table2asn_error
- NC_063943.1: table2asn_error
- NC_064206.1: boundary_gt30:nad6, table2asn_error
- NC_064339.1: boundary_gt30:rrnL, duplicate_gene:nad5, table2asn_error
- NC_065055.1: strand_error:Q, strand_error:A, strand_error:N, strand_error:C, strand_error:Y, strand_error:S, duplicate_gene:nad6, table2asn_error
- NC_065150.1: table2asn_error
- NC_065160.1: table2asn_error
- NC_065197.1: table2asn_error
- NC_066464.1: boundary_gt30:rrnL, strand_error:Y, table2asn_error
- NC_066800.1: table2asn_error
- NC_066903.1: boundary_gt30:atp6, boundary_gt30:cob, duplicate_gene:nad4, table2asn_error
- NC_066907.1: table2asn_error
- NC_067821.1: table2asn_error
- NC_067864.1: table2asn_error
- NC_068516.1: table2asn_error
- NC_068663.1: boundary_gt30:nad3, table2asn_error
- NC_068903.1: boundary_gt30:nad3, boundary_gt30:rrnL, internal_stop:atp6, duplicate_gene:nad5, table2asn_error
- NC_069308.1: boundary_gt30:cob, boundary_gt30:atp6, table2asn_error
- NC_069654.1: table2asn_error
- NC_069655.1: table2asn_error
- NC_070016.1: table2asn_error
- NC_070063.1: table2asn_error
- NC_070068.1: table2asn_error
- NC_071234.1: boundary_gt30:nad3, duplicate_gene:cob, duplicate_gene:nad6, table2asn_error
- NC_071873.1: table2asn_error
- NC_071908.1: table2asn_error
- NC_071939.1: table2asn_error
- NC_071975.1: boundary_gt30:nad5, table2asn_error
- NC_072211.1: boundary_gt30:rrnL, name_error:L, name_error:S, name_error:S, name_error:L, table2asn_error
- NC_072260.1: table2asn_error
- NC_072535.1: table2asn_error
- NC_072569.1: boundary_gt30:nad3, table2asn_error
- NC_072726.1: boundary_gt30:nad3, table2asn_error
- NC_073503.1: table2asn_error
- NC_073504.1: boundary_gt30:rrnL, table2asn_error
- NC_073529.1: boundary_gt30:nad1, table2asn_error
- NC_077264.1: boundary_gt30:nad3, boundary_gt30:rrnL, strand_error:S, strand_error:Y, strand_error:C, strand_error:N, strand_error:A, strand_error:Q, table2asn_error
- NC_077553.1: table2asn_error
- NC_077648.1: duplicate_gene:cox3, table2asn_error
- NC_079638.1: name_error:E, name_error:D, name_error:N, name_error:Q, table2asn_error
- NC_080894.1: table2asn_error
- NC_081024.1: internal_stop:atp6, duplicate_gene:cob, table2asn_error
- NC_081027.1: table2asn_error
- NC_081030.1: duplicate_gene:nad5, table2asn_error
- NC_081040.1: boundary_gt30:rrnL, table2asn_error
- NC_081043.1: table2asn_error
- NC_081187.1: boundary_gt30:rrnL, table2asn_error
- NC_082168.1: boundary_gt30:nad3, table2asn_error
- NC_082528.1: boundary_gt30:nad2, table2asn_error
- NC_082550.1: table2asn_error
- NC_082554.1: table2asn_error
- NC_082563.1: boundary_gt30:rrnL, table2asn_error
- NC_082565.1: table2asn_error
- NC_082578.1: boundary_gt30:cox1, table2asn_error
- NC_082848.1: table2asn_error
- NC_082993.1: table2asn_error
- NC_082995.1: table2asn_error
- NC_082996.1: table2asn_error
- NC_082997.1: table2asn_error
- NC_082998.1: boundary_gt30:rrnL, table2asn_error
- NC_082999.1: table2asn_error
- NC_083000.1: boundary_gt30:rrnL, table2asn_error
- NC_083001.1: boundary_gt30:rrnL, table2asn_error
- NC_083006.1: table2asn_error
- NC_083007.1: table2asn_error
- NC_083010.1: boundary_gt30:rrnL, table2asn_error
- NC_083011.1: boundary_gt30:nad2, table2asn_error
- NC_083027.1: table2asn_error
- NC_083030.1: table2asn_error
- NC_083035.1: table2asn_error
- NC_083064.1: table2asn_error
- NC_083073.1: table2asn_error
- NC_083076.1: table2asn_error
- NC_083162.1: boundary_gt30:rrnL, table2asn_error
- NC_083188.1: table2asn_error
- NC_083983.1: boundary_gt30:rrnL, strand_error:E, table2asn_error
- NC_083998.1: table2asn_error
- NC_084075.1: name_error:L, name_error:S, name_error:S, name_error:L, table2asn_error
- NC_084105.1: table2asn_error
- NC_085258.1: table2asn_error
- NC_085756.1: boundary_gt30:nad2, boundary_gt30:nad4, internal_stop:nad4, duplicate_gene:nad2, table2asn_error
- NC_086643.1: table2asn_error
- NC_086816.1: boundary_gt30:rrnL, table2asn_error
- NC_086960.1: boundary_gt30:rrnL, duplicate_gene:nad5, table2asn_error
- NC_086961.1: table2asn_error
- NC_087724.1: boundary_gt30:nad3, table2asn_error
- NC_087898.1: table2asn_error
- NC_087984.1: table2asn_error
- NC_087991.1: table2asn_error
- NC_088014.1: boundary_gt30:rrnL, table2asn_error
- NC_088027.1: boundary_gt30:rrnL, table2asn_error
- NC_088448.1: boundary_gt30:nad3, table2asn_error
- NC_133448.1: table2asn_error
