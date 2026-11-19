# 覆盖度计算比较 / Mitogenome read-depth methods compared (2026-11-19, cody)

数据 / Data: *Katharina tunicata* NC_001636.1（15,532 bp，环状 / circular）; wgsim 2×150 bp, error 0.2 %, simulated from a copy
rotated by half the length, so reads across the origin exist and the true depth of every base is known. Divergence =
wgsim SNP rate between the simulated individual and the reference. Wall time = index + map + count, best of 2, 128-core
machine. Script: `benchmark/coverage_bench.py`; raw numbers: `results.tsv`.

## 方法 / Methods (what each tool does, from its code)

| | 比对 / Mapper | 深度 / Depth | 环状 / Circular | 图 / Plot |
|---|---|---|---|---|
| using_bwa_draw_coveragev1.2.py | bwa mem, `view -F12`, sort | `samtools depth` (no `-a`) | no | flat (matplotlib) |
| MitoZ visualize (3.6) | bwa mem + samtools | `samtools depth` file | no | circos (Perl) depth ring |
| GetOrganelle `evaluate_assembly_using_mapping.py` | bowtie2 `--sensitive --dovetail --no-discordant -X` | SAM parsed in Python | **yes**: reference padded with the maximum insert size | flat, dots, log-like density, per-contig labels |
| `pyasm` (PyPI 1.1.5) | – | – | – | a Python **bytecode** library, nothing to do with sequences; not compared. 您指的可能是 pysam / mosdepth / pyfastx？ |
| `mitom coverage` (new) | bwa / bwa-mem2 / minimap2 / bowtie2 | SAM on a pipe → difference array | **yes**: first 1,000 nt appended, folded back | flat, circos, mitom |

## 结果 / Results

| 场景 / Scenario | 方法 / Method | 时间 / Wall s | Pearson r | 相对误差 / Rel. abs. error | 原点深度比 / Origin depth ratio |
|---|---|---|---|---|---|
| 200k pairs, 0 %, 8 thr | samtools-bwa (using_bwa / MitoZ) | 4.4 | 0.9986 | 0.0004 | **0.94** |
| | samtools-bwa-circ | 4.4 | 1.0000 | 0.0000 | 1.00 |
| | stream-bwa | 3.1 | 1.0000 | 0.0000 | 1.00 |
| | stream-bwa-mem2 | 4.7 | 1.0000 | 0.0000 | 1.00 |
| | stream-minimap2 | **2.4** | 1.0000 | 0.0000 | 1.00 |
| | stream-bowtie2 (GetOrganelle setting) | 6.5 | 0.9999 | 0.0027 | 1.00 |
| 200k pairs, 5 % SNP, 8 thr | samtools-bwa | 8.7 | 0.9948 | 0.0022 | 0.90 |
| | stream-bwa | 6.2 | **0.9985** | **0.0014** | 1.00 |
| | stream-bwa-mem2 | 11.7 | 0.9985 | 0.0014 | 1.00 |
| | stream-minimap2 | **2.5** | 0.9954 | 0.0038 | 1.00 |
| | stream-bowtie2 | 7.6 | 0.9938 | 0.0145 | 0.99 |
| 1M pairs, 0 %, 8 thr | samtools-bwa | 13.1 | 0.9986 | 0.0004 | 0.94 |
| | stream-bwa | 9.4 | 1.0000 | 0.0000 | 1.00 |
| | stream-bwa-mem2 | 9.6 | 1.0000 | 0.0000 | 1.00 |
| | stream-minimap2 | 8.7 | 1.0000 | 0.0000 | 1.00 |
| | stream-bowtie2 | 29.8 | 1.0000 | 0.0027 | 1.00 |
| 200k pairs, 0 %, 1 thr | samtools-bwa | 15.1 | 0.9986 | 0.0004 | 0.94 |
| | stream-bwa / bwa-mem2 / minimap2 | 12.9 / 14.0 / 13.6 | 1.0000 | 0.0000 | 1.00 |
| | stream-bowtie2 | 46.7 | 0.9999 | 0.0027 | 1.00 |

## 结论 / Conclusions

- **圆环处理决定准确性 / Circular handling matters most**: without it (using_bwa, MitoZ) the depth at the origin is 6–10 %
  low in the simulation; the same pipeline on an extended reference removes the error (samtools-bwa-circ).
- **流式计数更快 / Streaming count is faster** than BAM + sort + index + `samtools depth` (about 25–30 % less time, no BAM
  on disk) with identical numbers, so it is the default counter.
- **bwa mem 是默认比对器 / bwa mem is the default mapper**: most accurate when the reads differ from the reference
  (5 %); minimap2 `-x sr` is the fastest (2.5 s) and exact at 0 % but loses accuracy at 5 % (r 0.9954 vs 0.9985): use
  `--mapper minimap2` for same-individual reads or quick looks; bwa-mem2 gives the same numbers as bwa and is slower on
  this 16 kb reference (index time); bowtie2 is the slowest (3–5x) and the least accurate in relative error (end-to-end
  alignment loses read ends: 0.3–1.5 %).
- On 1M pairs the Python SAM parser, not the mapper, is the limit (minimap2 8.7 s vs bwa 9.4 s). If this matters, count in C.
- Simulated reads only: no real-data run yet (adapters, duplicates, NUMTs, low-complexity regions).
