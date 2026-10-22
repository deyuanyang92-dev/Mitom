# 动物线粒体注释软件对比（2026-10-22）：MitoM vs MitoZ vs MITOS2 vs MITOS1 vs MitoFinder

本轮在上一版（`docs/benchmark_metazoa_2026-10.md`）之后，按四个并行审查（蛋白基因逻辑对照 MitoZ、RNA 逻辑对照 MITOS2、流程集成、极端输入压力测试）修复了 MitoM 原生核心与流程中的问题，又在 255 条 RefSeq 记录上以 leave-out 方式调整了起止密码子规则，然后重新对比。

## 数据与方法

- **24 条主测试集**（`benchmark/sets/mito_bench_mini_2026-10.tsv`）：8 个门 × 3 条 NC_ 记录（环节、节肢、脊索、刺胞、棘皮、软体、线虫、扁形）。输入为去掉注释的裸序列；输出换算回 NCBI 坐标后逐基因与 RefSeq 比较，并用 table2asn 检查能否直接提交 NCBI。
- **补充**：58 条新 NC_ 记录（13 个门，`benchmark/sets/mito_bench_new_2026-10.tsv`）中的 6 条棘皮动物；MITOS1 只跑 8 条（每门 1 条），因为每条要约 40 分钟。
- **mitom_lo**：MitoM 在参考蛋白中排除该记录自身蛋白（`MITOM_EXCLUDE_REF`）。MITOS2 的 refseq89m 参考库包含 8,012 条 RefSeq 记录的蛋白，MITOS2 和 `mitom_native` 都用这套库；`mitom_lo` 给出新测序基因组可预期的精度。rRNA 一律使用不含测试记录的 HMM 库。
- 运行环境：同一台机器，同时并行运行多条，每条 2 线程。耗时和内存是每条记录整个进程树的数值（中位数 / 最大值）。

## 结果（24 条）

### 速度与内存

| 工具 | 完成 | 中位耗时 (s) | 中位峰值内存 (MB) | 最大峰值内存 (MB) |
|---|---|---|---|---|
| **MitoM** | 24/24 | **83** | **184** | **283** |
| MitoZ 3.6 | 18/24（无刺胞、扁形 clade） | 249 | 492 | 719 |
| MITOS2 2.1.10 | 24/24 | 385 | 953 | 1929 |
| MitoFinder 1.4.2 | 22/24（2 条缺近缘参考） | 204 | 238 | 283 |
| MITOS1 1.1.7 | 4/8（4 条程序自身崩溃） | 2465 | 775 | 788 |

### 准确度：五个工具都完成的 17 条记录

| 工具 | 蛋白基因两端完全一致 | 蛋白基因 ≤3 nt | rRNA 找回 | rRNA ≤3 nt | tRNA 找回 | tRNA ≤3 nt |
|---|---|---|---|---|---|---|
| **MitoM** | **83.9%** | **91.3%** | **100%** | **59.4%** | **96.8%** | **81.9%** |
| MitoM（leave-out） | 81.7% | 89.0% | 100% | 59.4% | 96.8% | 81.9% |
| MitoZ | 78.4% | 86.2% | 87.5% | 28.1% | 92.5% | 79.2% |
| MITOS2 | 50.5% | 57.8% | 96.9% | 25.0% | 96.0% | 80.9% |
| MitoFinder | 59.2% | 72.0% | 90.6% | 28.1% | 96.5% | 79.0% |

### 明显错误（17 条共同记录）

| 错误 | MitoM | MitoZ | MITOS2 | MitoFinder |
|---|---|---|---|---|
| NCBI table2asn ERROR | **0** | 12 | 17 | 17 |
| 缺基因 | **0** | 4 | 3 | 8 |
| 边界差 >30 nt | **8** | 14 | 48 | 29 |
| 内部终止密码子 | **0** | 3 | 0 | 3 |
| 有任一错误的记录 | **7/17** | 16/17 | 17/17 | 17/17 |

### 按门：蛋白基因两端完全一致

| 门 | MitoM | MitoM leave-out | MitoZ | MITOS2 | MitoFinder |
|---|---|---|---|---|---|
| 环节 | **76.9%** | 74.4% | 61.5% | 43.6% | 33.3% |
| 节肢 | **94.9%** | 94.9% | 89.7% | 53.8% | 69.2% |
| 脊索 | **87.2%** | 87.2% | 82.1% | 79.5% | 69.2% |
| 刺胞 | **94.9%** | 92.3% | 不支持 | 69.2% | 87.2% |
| 棘皮（3 条） | 89.7% | 89.7% | **92.3%** | 48.7% | 74.4% |
| 软体 | **60.5%** | 60.5% | 47.4% | 42.1% | 39.5% |
| 线虫 | **86.1%** | 75.0% | 80.6% | 27.8% | 50.0% |
| 扁形 | **70.3%** | 62.2% | 不支持 | 59.5% | 32.4% |

**棘皮动物**：24 条集中的 3 条只有 39 个蛋白基因，MitoM 比 MitoZ 少对 1 个。这个基因是 NC_023772.1 的 cob。NCBI 序列的起点正好落在它的终止密码子 TAA 里，所以 NCBI 写成末端的 TA 加 `/transl_except`。MitoM 先把序列旋转为 cox1 开头，因此注释的是完整的 TAA，换算回 NCBI 坐标后差 1 nt。这是坐标惯例造成的差异，两种写法都正确。

另加 58 条集中的 6 条棘皮动物，共 9 条、117 个基因：

| 工具 | 两端完全一致 | 有明显错误的记录（6 条） | 中位耗时 | 中位内存 |
|---|---|---|---|---|
| **MitoM** | **109/117（93.2%）** | **0/6** | **93 s** | **194 MB** |
| MitoZ | 107/117（91.5%） | 5/6 | 261 s | 556 MB |
| MITOS2 | 69/117（59.0%） | 6/6 | 455 s | 1029 MB |
| MitoFinder | 79/117（67.5%） | 6/6 | 261 s | 250 MB |

### MITOS1（4 条成功记录）

MITOS1（bioconda mitos 1.1.7、Infernal 1.0、Zenodo 2683856 参考数据）在 8 条中有 4 条崩溃：rRNA 合并步骤里 `cmrealign` 的输出文件缺失。它的 RNAplot 画图步骤也会崩溃，评测中只给这一步加了异常保护，注释逻辑没有改动。

下表是 MITOS1 成功的 4 条记录上各工具的结果：

| 工具 | 中位耗时 | 中位内存 | 蛋白基因找回 | 蛋白两端一致 | 蛋白 ≤3 nt | tRNA ≤3 nt |
|---|---|---|---|---|---|---|
| **MitoM** | **84 s** | **187 MB** | 98% | **80%** | **84%** | 84% |
| MitoZ | 292 s | 495 MB | 100% | 68% | 68% | 88% |
| MITOS2 | 385 s | 1108 MB | 96% | 47% | 47% | 84% |
| MITOS1 | 2483 s | 785 MB | 88% | 0% | 14% | 88% |
| MitoFinder | 211 s | 261 MB | 95% | 74% | 79% | 83% |

MITOS1 不确定精确的起止密码子，所以蛋白基因没有一个两端完全一致。

## 本轮修复与改进（均有回归测试）

逻辑审查发现的问题：

- **跨原点基因**：环状序列改为按"序列 + 前 3 kb"搜索，跨原点的蛋白基因、rRNA、tRNA 能完整找到，并写成 join。以前会被截断或丢失；MITOS2 对此是搜索两个旋转。
- **旁系同源误注**：如 atp8 落在 atp6 内、nad2 落在 nad5 内。现在丢弃主要落在得分更高的另一基因上的 CDS，再到其他位点重试，并把已注释的基因屏蔽掉。
- **exonerate**：跳过反链比对（以前会崩溃）。遗传密码 24 和 33 以密码表字符串传入（以前全部基因失败）。只有在比对框内有 ≥2 个终止密码子时才改读码框。
- **MiTFi**：无反密码子时，MiTFi 在氨基酸列写的是模型名，以前被误读为 Thr。tRNA copy 的阈值改为 E ≤ 0.001。敏感搜索改为按 tRNA 名称判断（与 MITOS2 一致）。
- **rRNA**：按模型列拼接命中；覆盖模型不足 40% 的跨度交给 rRNA HMM 步骤处理。
- **curation 与内含子**：带模糊端点的单区间 CDS 也做修正。终点停在 rRNA 前的基因可以延伸。O1 规则用于距原点更远的末端。被切断的 CDS 补到序列末端。内含子末端只在已有 CDS 末端带终止密码子时才沿用。不再在线性记录两端之间或超过 13 kb 的跨度上拼出"内含子"。
- **流程**：核心失败会在控制台显示，并返回非零退出码。原生核心不再打印 MitoZ/MITOS2 相关信息。安装脚本会更新已有环境。线性记录的 `cox1_first` 字段如实填写。run.log 每次运行重写。遗传密码与记录的 transl_table 不一致时给出警告。多 contig 和 <10 kb 的 FASTA 默认按线性处理。

精度调整（在 255 条 RefSeq 记录上以 leave-out 核心输出做 A/B 测试；PCG 两端一致 2186 → 2261）：

- 比对缺 ≥2 个 N 端残基时，在期望位置附近找起始密码子。
- 非 ATG 起始改为上游最远的同框 ATG。线虫和节肢动物例外：它们常用 TTG/ATT/ATA 起始。
- 脊索动物：完整终止与下游基因重叠 1–2 nt 时写成 T/TA，这是脊椎动物 RefSeq 的惯例。
- 下游没有 T/TA 时，若完整终止伸入下一特征 ≤20 nt 则采用它。核心比对已到达的完整终止若只伸入预测的 tRNA ≤20 nt，则保留。
- 终止只需跨原点 1–2 nt 时，在序列末端写 T/TA。
- exonerate 只比对出片段、而 blastx 覆盖全长时，用 blastx 命中定位基因（与 MITOS2 做法相同）。只有微弱命中的基因再做一次敏感搜索。前 3 个参考失败时，改试接下来 3 个。
- 刺胞动物不做 tRNA `--max` 敏感搜索：它们只编码 1–6 个 tRNA，这一步主要产生假阳性（17 → 4），中位耗时 344 → 50 s。

效率：

- 原生核心默认不再调用 MITOS2 补救（`--rescue yes` 仍可开启）。以前它在每条缺基因的记录上花 5–7 分钟、约 1 GB，24 条里只补回 2 个基因；现在核心自己重试就能找回这两个。
- 样本数不超过 `--max_tasks` 时不启动 Snakemake（省约 120 MB）。
- blastx 每个进程 2 线程、并发进程数减半。
- Infernal 比对矩阵上限设为 128 MB；内含子分型的内存 293 → 183 MB，结果相同。

单条记录的峰值内存 318 → 190 MB，中位耗时 109 → 83 s。

## 安装与功能

本机实测（Linux x86_64，miniforge）。功能对比的细项见 `docs/benchmark_metazoa_2026-10.md` 的 Features 一节；本轮只有 MitoM 的安装方式变了，MITOS1 是新加入的。

| | MitoM | MitoZ 3.6 | MITOS2 2.1.10 | MITOS1 1.1.7 | MitoFinder 1.4.2 |
|---|---|---|---|---|---|
| 安装 | `bash install.sh`：一个 conda 环境（mitotools）+ refseq89m 参考数据（自动下载、md5 校验） | `conda create -n mitoz -c bioconda mitoz=3.6` | bioconda + 手动下载 refseq89m | bioconda + 手动下载两份参考数据（Zenodo 2683856），画图步骤需修补才能跑完 | GitHub 源码 + Python 2.7 + Java + 手改配置；每个物种要自备近缘参考 |
| 磁盘 | 1.1 GB 环境 + 93 MB 参考数据 + 13 MB 模型 | 4.8 GB | 2.8 GB + 93 MB | 2.3 GB + 54 MB | 0.5 GB + 0.65 GB 源码 |
| 类群 | 全部后生动物 | 9 个 clade（无刺胞、扁形） | 全部后生动物 | 全部后生动物 | 取决于参考 |
| 输入 | GenBank（类群、遗传密码从记录读出）或 FASTA，可批量 | 一次一个 FASTA | 一次一条序列 | 一次一条序列 | 一次一个组装 + 参考 |
| 产出 | NCBI 格式 GenBank、table2asn 校验、提交文件、原注释对比、每步报告 | gbf/tbl/FASTA/图 | bed/gff/FASTA/特征表（无 GenBank） | bed/gff/FASTA | gb/tbl/gff/FASTA |
| 内含子 / 跨原点基因 / 边界校正 / 去重 | 有 | 无 | 内含子片段，无边界 | 无 | 无 |

## 复现

```bash
cd benchmark
export MITOM_RRNA_LIB=<不含测试记录的 rRNA HMM 库>      # 公平测试 rRNA
python3 mito_bench_run.py --set sets/mito_bench_mini_2026-10.tsv --root <root> --split all \
    --tools mitom_native,mitom_lo,mitoz,mitos2,mitofinder --jobs 24 --threads 2
python3 mito_bench_run.py --set sets/mito_bench_mini_2026-10.tsv --root <root> --split all --tools mitos1 \
    --acc NC_069205.1,NC_067900.1,NC_009732.1,NC_027103.1,NC_063787.1,NC_083122.1,NC_001861.1,NC_035255.1 \
    --jobs 8 --threads 1 --timeout 10800
python3 mito_bench_eval.py --set sets/mito_bench_mini_2026-10.tsv --root <root> --split all \
    --tools mitom_native,mitom_lo,mitoz,mitos2,mitofinder --table2asn
python3 bench_report.py --eval <root>/eval/all --tools mitom_native,mitom_lo,mitoz,mitos2,mitofinder \
    --set sets/mito_bench_mini_2026-10.tsv
```

结果文件：`benchmark/results/six_tools_2026-10-22/`（`report_24.md`、`report_echinodermata_6.md`、逐基因 `24_features.tsv.gz`、耗时 `24_timing.tsv`）。
