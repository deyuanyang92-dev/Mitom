# 主要功能的评测摘要 / Benchmark summary of the main features

本页只放**用户最常问的数字**。方法、留一设置、逐条结果在 [METHODS.md](METHODS.md) 和 `benchmark/results/`。开发期基准（2026-10 / 11），不是论文冻结评测（见 `benchmark/frozen/PLAN.md`）。

This page holds only the numbers users ask about most. Methods, the leave-one-out setup and per-record results are in [METHODS.md](METHODS.md) and `benchmark/results/`. These are development benchmarks (2026-10 / 11), not the frozen evaluation for the paper (see `benchmark/frozen/PLAN.md`).

使用与示例见 [USAGE_zh.md](USAGE_zh.md)、[EXAMPLES.md](EXAMPLES.md)。
Usage and examples: [USAGE_zh.md](USAGE_zh.md), [EXAMPLES.md](EXAMPLES.md).

---

## 1. 端到端注释（相对 NCBI） / End-to-end annotation (vs NCBI)

| 指标 / Metric | 结果 / Result | 说明 / Notes |
|---|---|---|
| 蛋白基因两端与 NCBI 完全一致 / Protein-gene ends identical to NCBI | 1,050 条多类群 RefSeq：24 个门全部高于 MITOS2 / multi-phylum RefSeq records: all 24 phyla above MITOS2 | METHODS 校正 / tRNA 相关节 / METHODS curation and tRNA sections |
| NCBI table2asn ERROR | 2,497 条基准记录：**0**（此前 199）/ benchmark records: **0** (199 before) | 含多类群 RefSeq、环节动物、未参与训练的无脊椎、脊椎动物等 / multi-phylum RefSeq, annelids, untrained invertebrates, vertebrates, etc. |
| tRNA（未参与训练的科）/ tRNA (families not in training) | 召回高于 MITOS2，精确率持平或更高 / recall above MITOS2, precision equal or higher | METHODS tRNA 节 / METHODS tRNA section |

---

## 2. 类群与遗传密码（`--clade auto`）
*Clade and genetic code (`--clade auto`)*

数据：2026-10-31，本机分层抽样（每门 ≤ 15 条），去掉同科参考。
Data: 2026-10-31, stratified sample on this machine (at most 15 records per phylum), same-family references removed.

| 测试 / Test | 记录 / Records | 遗传密码正确 / Genetic code correct | 门：正确 / 不确定 / 错误 / Phylum: correct / uncertain / wrong |
|---|---|---|---|
| 完整基因组 / Complete genomes | 238 | 93.3% | 63.0 / 32.4 / 4.6% |
| 8 kb / 5 kb / 3 kb 片段 / fragments | 238×3 | 91.6 / 85.3 / 74.4% | 正确比例随片段变短下降 / correct share falls as fragments get shorter |
| 全流程 55 条无标签 FASTA（29 门，排除自身参考）/ full run, unlabelled FASTA (29 phyla, own reference excluded) | 55 | 52 / 55 | 29 / 23 / 3；NCBI ERROR 0 |

门判错的多为参考库里只有 1–2 个科的小门；遗传密码仍常正确。知道类群时请写 `--clade` 与 `--genetic_code`。
Wrong phyla are mostly small phyla with very few families in the reference set; the genetic code is often still right. If you know the clade, give `--clade` and `--genetic_code`.

细节：[METHODS.md § 类群与遗传密码](METHODS.md#类群与遗传密码自动判断----clade-auto-mito_classifypy-mitom_groupspy)
Details: [METHODS.md § clade and genetic code](METHODS.md#类群与遗传密码自动判断----clade-auto-mito_classifypy-mitom_groupspy)

---

## 3. 注释后校正（相对 NCBI） / Post-annotation curation (vs NCBI)

31 条鳞沙蚕类真实 MitoZ 注释，对照 NCBI：
Real MitoZ annotations of scale worms, compared with NCBI:

| | MitoZ 原始 / MitoZ raw | 校正后 / After curation |
|---|---|---|
| PCG 5′ / 3′ 一致（284 个）/ PCG 5′ / 3′ identical | 253 / 211 | **265 / 254** |
| 12S 5′ / 3′（22 个） | 5 / 2 | **17 / 15** |
| 16S 5′ / 3′（22 个） | 0 / 0 | **17 / 11** |
| rRNA 与相邻基因重叠 / rRNA overlapping a neighbour gene | 44 | **0** |
| table2asn ERROR | 128 | **4**（报告中已标需人工检查 / flagged in the report for manual check） |

细节：[METHODS.md § 注释后校正](METHODS.md#注释后校正--post-annotation-curation-mito_curatepy)
Details: [METHODS.md § post-annotation curation](METHODS.md#注释后校正--post-annotation-curation-mito_curatepy)

---

## 4. 内含子
*Introns*

金标准：193 条核实的 NCBI 带内含子记录 → 309 个内含子。`mitom run`（裸序列）：
Gold standard: verified NCBI records with introns (see table). `mitom run` on bare sequences:

| 口径 / Setting | 完全一致 / Exact match | 漏检 / Missed | 76 条无内含子记录的假阳性 / False positives on 76 intron-free records | NCBI ERROR |
|---|---|---|---|---|
| 默认（完整位点库），194 条 / Default (full site library), 194 records | **309 / 309** | 0 | 0 | 0 |
| 留一（本记录不进库），270 条 / Leave-one-out (own record not in library), 270 records | **299 / 309** | 6 | 0 | 0 |

细节：[METHODS.md § intron](METHODS.md#intron-注释--introns-mito_intronpy-mito_intron_refpy)，目录 `benchmark/results/intron_*`
Details: [METHODS.md § introns](METHODS.md#intron-注释--introns-mito_intronpy-mito_intron_refpy); folders `benchmark/results/intron_*`

---

## 5. 其它功能（入口） / Other features (entry points)

| 功能 / Feature | 结果目录 / 说明 / Results folder / notes |
|---|---|
| 漏注蛋白基因（近缘谱 HMM）/ Missed protein genes (relative-profile HMM) | `benchmark/results/pcg_hmm_2026-11-10/` |
| 未知 ORF / Unknown ORFs | `benchmark/results/orf_2026-11-09/` |
| RNA-seq 证据 / RNA-seq evidence | `benchmark/results/rnaseq_2026-11-10/` |
| 多分体线粒体 / Multipartite mitogenomes | `benchmark/results/ebp_mito_2026-11-10/` |
| 六工具对比 / Six-tool comparison | [benchmark_six_tools_2026-10-22.md](benchmark_six_tools_2026-10-22.md) |
| 环节 / 多门大表 / Annelid / multi-phylum tables | [benchmark_annelida_2026-10.md](benchmark_annelida_2026-10.md)、[benchmark_metazoa_2026-10.md](benchmark_metazoa_2026-10.md) |
| 核 rDNA / 组蛋白 / Nuclear rDNA / histones | [METHODS.md § 核基因](METHODS.md#核基因标记--nuclear-markers-rdna-and-histones-mito_nucpy) |

---

## 6. 安装体积与用时（参考）
*Install size and run time (reference)*

全新安装（2026-10-28，Linux x86_64：conda + refseq89m）：
Fresh install (2026-10-28, Linux x86_64: conda + refseq89m):

| 软件 / Tool | 下载 / Download | 安装后 / Installed |
|---|---|---|
| MitoM | 0.48 GB | 1.8 GB（首次再解压约 0.19 GB tRNA 模型 / tRNA models unpacked on first run） |
| MITOS2 2.1.10 | 0.75 GB | 2.9 GB |
| MitoZ 3.6 | 1.37 GB | 4.5 GB |

单条与批量用时见仓库 README「用时」表与 METHODS 各节。
Times for single and batch runs: see the tables in the repository README and the METHODS sections.
