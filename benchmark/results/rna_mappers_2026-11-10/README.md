# Which mapper for RNA-seq on a mitochondrial genome (2026-11-10, cody)

问题：`mitom rna` 把转录组 reads 比对到线粒体基因组上，看各基因的表达和 poly(A) 3′ 端。比对软件很多，哪个最合适？
Question: `mitom rna` maps RNA-seq reads to the mitogenome (expression per gene, poly(A) 3′ ends). Which mapper?

线粒体转录本的特点决定了评价标准 / what matters here:

- 读到 3′ 端的 reads 末尾带 poly(A)，比对软件必须软剪切（soft clip）它，才能定 3′ 端；
  reads at a 3′ end carry a poly(A) tail that must be soft-clipped;
- 动物线粒体基因没有剪接内含子，剪接比对软件报出的 "N" 都是假的；no spliceosomal introns: any N is false;
- 转录组样本常常不是组装那一个体（例如 *Erpobdella* 差 3.5 %），比对要耐受分化；the RNA-seq individual often differs
  from the assembled one;
- 环状基因组，跨原点的 reads；reads crossing the origin of a circular genome.

## Method / 方法

1. 模拟（`simulate.py`，有真值）：3 个环节动物线粒体基因组（*Polydora* OQ078742、*Erpobdella* OZ022248、多鳞虫
   OZ414225，均为 MitoM 注释），每个蛋白基因和 rRNA 的成熟转录本带 20–60 nt poly(A)，表达量对数正态；15 % 为 1–4 kb
   前体片段（含跨原点）；个体与参考差 0 / 1 / 3 %（SNP + 1/10 的 1–3 nt 插入缺失），测序错误 0.3 %；双端 2×61 和
   2×150，每组 20 万对，dUTP 链特异。共 18 组 × 13 种设置。
   Simulation with truth: 3 genomes × divergence 0 / 1 / 3 % × 2×61 / 2×150, 200,000 pairs each, 13 settings.
2. 评价（`evaluate.py`）：位置正确（最左比对碱基与真值差 ≤ 5 nt）；MAPQ ≥ 20 中错位的比例；剪接 reads；带 ≥ 6 nt 尾的
   reads 中尾被看见的比例；用 `mito_rna` 自己的规则找 poly(A) 位点，真 3′ 端找回率和位点正确率（±3 nt）；基因深度与
   真值的中位 |log2 比|；跨原点 reads 比上的比例；时间、内存（8 线程）。
3. 真实数据（`evaluate_real.py`，无真值，前 400 万对）：*Polydora* cf. *ciliata* SRR26830859（61 bp，随机引物）对
   OQ078742；*Erpobdella octoculata* SRR10997416（cDNA）对同种长读长组装 OZ022248。
4. 软件 / versions（conda env `mapbench`）：minimap2 2.31，bwa 0.7.19，bwa-mem2，bowtie2 2.5.5，HISAT2 2.2.3，STAR
   2.7.10b，BBMap 40.05，Subread，segemehl。命令见 `run_mappers.py`。

## Results / 结果

### Simulated / 模拟（18 组平均；"3 %, 61 bp" 列为最难的一组）

| 设置 / setting | 位置正确 % / correct | 3 %, 61 bp | MAPQ≥20 中错位 % / wrong at MAPQ≥20 | 假剪接 reads / false spliced | 尾被看见 % / tails seen | poly(A) 3′ 端找回 / recall | 位点正确率 / precision | 深度误差 / depth \|log2\| | 跨原点比上 % / origin | 秒 / s | MB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **bwa mem** | **98.2** | **96.7** | 0.41 | 0 | 81 | **0.96** | **0.99** | **0.003** | 100 | 5.8 | 281 |
| bwa-mem2 | 98.2 | 96.7 | 0.41 | 0 | 81 | 0.96 | 0.99 | 0.003 | 100 | 5.5 | 374 |
| BBMap local=t | 97.8 | 96.5 | 0.26 | 0 | 57 | 0.88 | 0.72 | 0.012 | 89 | 42.5 | 1802 |
| bowtie2 --very-sensitive-local | 97.4 | 94.8 | 0.84 | 0 | 82 | 0.91 | 0.90 | 0.009 | 100 | 13.4 | 132 |
| segemehl | 95.3 | 94.3 | 1.40 | 0 | 0 | 0.00 | — | 0.019 | 19 | 27.5 | 133 |
| minimap2 -ax splice:sr | 93.7 | 80.9 | 1.61 | 2750 | 54 | 0.94 | 0.96 | 0.076 | 81 | 2.6 | 206 |
| bowtie2（端到端 / end-to-end） | 93.6 | 83.6 | 0.00 | 0 | 0 | 0.00 | — | 0.135 | 19 | 5.6 | 96 |
| minimap2 -ax sr（旧默认 / old default） | 93.5 | 76.3 | 0.36 | 0 | 66 | 0.93 | 0.99 | 0.163 | 89 | 4.8 | 213 |
| STAR | 89.4 | 91.5 | 1.99 | 2865 | 65 | 0.90 | 0.61 | 0.161 | 61 | 51.0 | 938 |
| STAR --alignIntronMax 1 | 89.2 | 91.0 | 3.61 | 0 | 66 | 0.88 | 0.74 | 0.130 | 60 | 8.0 | 940 |
| Subread | 81.9 | 84.9 | 2.92 | 0 | 70 | 0.71 | 0.92 | 1.536 | 95 | 12.1 | 1166 |
| HISAT2 | 77.2 | 53.7 | 0.22 | 1627 | 5 | 0.54 | 0.91 | 0.386 | 12 | 2.5 | 27 |
| HISAT2 --no-spliced-alignment | 77.0 | 53.8 | 0.18 | 0 | 5 | 0.57 | 1.00 | 0.384 | 12 | 2.5 | 22 |

每组的数字：`simulated.tsv`（一行一个软件 × 一组条件）。Per condition: `simulated.tsv`.

### Real / 真实数据（*Erpobdella*，前 400 万对；转录组个体与组装差 3.5 % / RNA-seq individual 3.5 % from the assembly）

| 设置 / setting | MAPQ≥20 reads / 1000 对 / pairs | 蛋白基因深度中位数 / PCG median depth | poly(A) reads | 位点在基因 3′ 端 / sites at gene ends | 假剪接 / spliced | 秒 / s |
|---|---|---|---|---|---|---|
| bowtie2 --very-sensitive-local | 75.5 | 251 | 59343 | 10 | 0 | 29.6 |
| **bwa mem** / bwa-mem2 | 69.2 | **289** | 10175 | 5 | 0 | 66.4 / 29.6 |
| minimap2 -ax splice:sr | 41.1 | 158 | 900 | 3 | 419 | 21.5 |
| BBMap local=t | 40.0 | 265 | 311 | 9 | 0 | 25.4 |
| segemehl | 33.6 | 219 | 0 | 0 | 0 | 708.7 |
| minimap2 -ax sr（旧默认 / old default） | 32.8 | 80 | 464 | 1 | 0 | 20.2 |
| Subread | 20.4 | 11 | 920 | 3 | 0 | 167.3 |
| STAR --alignIntronMax 1 | 17.2 | 82 | 65 | 4 | 0 | 40.2 |
| HISAT2 | 17.6 | 4 | 18 | 0 | 28 | 64.0 |
| bowtie2（端到端 / end-to-end） | 16.8 | 15 | 0 | 0 | 0 | 24.0 |
| STAR | 16.2 | 68 | 49 | 3 | 704 | 52.0 |

*Polydora* 前 400 万对中线粒体 reads 太少（蛋白基因深度中位数 4），各软件结果相同，不能区分（`real.tsv`）。
Polydora's first 4 M pairs hold too few mitochondrial reads to tell the mappers apart (`real.tsv`).

## Conclusion / 结论

- 默认改为 **bwa mem**（`mito_rna.py --mapper bwa`）：模拟中位置、poly(A)、深度都最好，个体分化 3 % 时仍稳定；真实
  分化样本上比旧默认 minimap2 多比上一倍多的 reads。minimap2（`--mapper minimap2`）和 bowtie2 局部模式
  （`--mapper bowtie2`）保留为选项。bwa mem is the new default; minimap2 and bowtie2 local stay as options.
- 不要用：端到端比对（bowtie2 默认、segemehl 不软剪切，看不到 poly(A)）、剪接比对（STAR、HISAT2、minimap2 splice
  报假内含子）、HISAT2 / Subread（分化时丢 reads）。Avoid end-to-end mappers (no tails), spliced mappers (false
  introns), HISAT2 / Subread (lose reads when the individual differs).
- 真实数据另发现 oligo-dT 的内部引发（internal priming）：*Erpobdella* reads 最多的 "3′ 端"（7,080 条）停在基因组
  自身 9 个 A 处，不是真末端。`mito_rna` 现在把位点前 10 nt 含 ≥ 7 个连续 A 或共 ≥ 9 个 A 的位点标为
  internal_priming，不计入 3′ 端证据（67 个环节动物基因组 1,047 个真 3′ 端中只有 1 个满足此条件）。
  Real data also showed oligo-dT internal priming (the 7,080-read "end" sits on a genomic A9 run); such sites are now
  flagged and not counted.

## Files / 文件

`simulate.py`, `run_mappers.py`, `evaluate.py`, `evaluate_real.py`, `run.sh`, `run_real.sh`; `simulated.tsv`,
`real.tsv`. Work data on cody: `/media/cody/18T_R1/yzz/_realtest/rna_mappers_2026-11-10/`. libmamba's TLS fails
through `ALL_PROXY=socks://...`; create the env with `env -u ALL_PROXY -u all_proxy mamba create ...`.
