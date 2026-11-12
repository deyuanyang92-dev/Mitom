# 主要功能的评测摘要

本页只放**用户最常问的数字**。方法、留一设置、逐条结果在 [METHODS.md](METHODS.md) 和 `benchmark/results/`。开发期基准（2026-10 / 11），不是论文冻结评测（见 `benchmark/frozen/PLAN.md`）。

使用与示例见 [USAGE_zh.md](USAGE_zh.md)、[EXAMPLES.md](EXAMPLES.md)。

---

## 1. 端到端注释（相对 NCBI）

| 指标 | 结果 | 说明 |
|---|---|---|
| 蛋白基因两端与 NCBI 完全一致 | 1,050 条多类群 RefSeq：24 个门全部高于 MITOS2 | METHODS 校正 / tRNA 相关节 |
| NCBI table2asn ERROR | 2,497 条基准记录：**0**（此前 199） | 含多类群 RefSeq、环节动物、未参与训练的无脊椎、脊椎动物等 |
| tRNA（未参与训练的科） | 召回高于 MITOS2，精确率持平或更高 | METHODS tRNA 节 |

---

## 2. 类群与遗传密码（`--clade auto`）

数据：2026-10-31，本机分层抽样（每门 ≤ 15 条），去掉同科参考。

| 测试 | 记录 | 遗传密码正确 | 门：正确 / 不确定 / 错误 |
|---|---|---|---|
| 完整基因组 | 238 | 93.3% | 63.0 / 32.4 / 4.6% |
| 8 kb / 5 kb / 3 kb 片段 | 238×3 | 91.6 / 85.3 / 74.4% | 正确比例随片段变短下降 |
| 全流程 55 条无标签 FASTA（29 门，排除自身参考） | 55 | 52 / 55 | 29 / 23 / 3；NCBI ERROR 0 |

门判错的多为参考库里只有 1–2 个科的小门；遗传密码仍常正确。知道类群时请写 `--clade` 与 `--genetic_code`。

细节：[METHODS.md § 类群与遗传密码](METHODS.md#类群与遗传密码自动判断----clade-auto-mito_classifypy-mitom_groupspy)

---

## 3. 注释后校正（相对 NCBI）

31 条鳞沙蚕类真实 MitoZ 注释，对照 NCBI：

| | MitoZ 原始 | 校正后 |
|---|---|---|
| PCG 5′ / 3′ 一致（284 个） | 253 / 211 | **265 / 254** |
| 12S 5′ / 3′（22 个） | 5 / 2 | **17 / 15** |
| 16S 5′ / 3′（22 个） | 0 / 0 | **17 / 11** |
| rRNA 与相邻基因重叠 | 44 | **0** |
| table2asn ERROR | 128 | **4**（报告中已标需人工检查） |

细节：[METHODS.md § 注释后校正](METHODS.md#注释后校正--post-annotation-curation-mito_curatepy)

---

## 4. 内含子

金标准：193 条核实的 NCBI 带内含子记录 → 309 个内含子。`mitom run`（裸序列）：

| 口径 | 完全一致 | 漏检 | 76 条无内含子记录的假阳性 | NCBI ERROR |
|---|---|---|---|---|
| 默认（完整位点库），194 条 | **309 / 309** | 0 | 0 | 0 |
| 留一（本记录不进库），270 条 | **299 / 309** | 6 | 0 | 0 |

细节：[METHODS.md § intron](METHODS.md#intron-注释--introns-mito_intronpy-mito_intron_refpy)，目录 `benchmark/results/intron_*`

---

## 5. 其它功能（入口）

| 功能 | 结果目录 / 说明 |
|---|---|
| 漏注蛋白基因（近缘谱 HMM） | `benchmark/results/pcg_hmm_2026-11-10/` |
| 未知 ORF | `benchmark/results/orf_2026-11-09/` |
| RNA-seq 证据 | `benchmark/results/rnaseq_2026-11-10/` |
| 多分体线粒体 | `benchmark/results/ebp_mito_2026-11-10/` |
| 六工具对比 | [benchmark_six_tools_2026-10-22.md](benchmark_six_tools_2026-10-22.md) |
| 环节 / 多门大表 | [benchmark_annelida_2026-10.md](benchmark_annelida_2026-10.md)、[benchmark_metazoa_2026-10.md](benchmark_metazoa_2026-10.md) |
| 核 rDNA / 组蛋白 | [METHODS.md § 核基因](METHODS.md#核基因标记--nuclear-markers-rdna-and-histones-mito_nucpy) |

---

## 6. 安装体积与用时（参考）

全新安装（2026-10-28，Linux x86_64：conda + refseq89m）：

| 软件 | 下载 | 安装后 |
|---|---|---|
| MitoM | 0.48 GB | 1.8 GB（首次再解压约 0.19 GB tRNA 模型） |
| MITOS2 2.1.10 | 0.75 GB | 2.9 GB |
| MitoZ 3.6 | 1.37 GB | 4.5 GB |

单条与批量用时见仓库 README「用时」表与 METHODS 各节。
