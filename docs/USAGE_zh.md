# MitoM 使用手册

> **说明**：MitoM 目前仍在测试中，供本团队内部使用，尚未正式发布。结果（尤其是内含子、tRNA 和基因边界）请人工核对后再用于发表或提交 NCBI。

MitoM（命令 `mitom`）注释**动物线粒体基因组**：找出 13 个蛋白基因、2 个 rRNA、22 个 tRNA 和内含子，统一方向与起点，按 NCBI 规范修整边界和基因名，生成可提交的文件与环状基因组图。另有 `mitom nuc`：在全基因组组装中注释核 rDNA（18S–ITS1–5.8S–ITS2–28S）和组蛋白。

---

## 本手册怎么读

| 需求 | 去哪 |
|---|---|
| 安装与各命令的输入 / 输出 | **本页** |
| 用仓库测试数据跟做 | [EXAMPLES.md](EXAMPLES.md) |
| 主要功能的数字 | [BENCHMARK.md](BENCHMARK.md) |
| 方法与完整评测表 | [METHODS.md](METHODS.md) |
| 内部参考库怎么制备 | [DATASETS.md](DATASETS.md) |
| 全部参数（由 `-h` 生成） | [OPTIONS.md](OPTIONS.md) |
| 与 MITOS2 联用 | [MITOS2.md](MITOS2.md) |
| 文档总目录 | [README.md](README.md) |

测出来的数字只在 [BENCHMARK.md](BENCHMARK.md) 与 [METHODS.md](METHODS.md)。本页不复述。

---

## 1. 特色

- **一条命令到可提交文件**：`mitom run` 从序列做到定向、补漏、内含子、未知 ORF、边界校正、table2asn 与圈图。
- **多门混合输入**：默认 `--clade auto`，按门 / 遗传密码分组注释，再合并到同一输出目录。
- **校正对着 NCBI 规则**：蛋白基因、rRNA、tRNA、原点规则可复查；报告在 `qc/curation_report.tsv`。
- **内含子两步**：从头识别 + 近缘位点迁移；无类型内含子可用同源连续性规则。
- **未知 ORF**：同源证明受选择才写入；其余进 WARNING 与 NCR 备注。
- **核标记独立命令**：`mitom nuc` 做 rDNA 与组蛋白，不混进线粒体流程。
- **可复现**：`logs/run_params.json`、`mitom doctor`、`mitom debug`。

数字见 [BENCHMARK.md](BENCHMARK.md)。参考库制备见 [DATASETS.md](DATASETS.md)。

---

## 2. 安装

### 2.1 一键安装（推荐）

需要 Linux 或 Windows WSL2，以及 conda / mamba（推荐 [Miniforge](https://github.com/conda-forge/miniforge)）。

```bash
git clone https://github.com/deyuanyang92-dev/Mitom.git MitoM
cd MitoM
bash install.sh        # 建 conda 环境 mitotools，装 mitom，并运行 mitom setup
mitom doctor           # 必需程序都应显示 OK
mitom --version
```

`install.sh` 会：按 `env/mitotools.yml` 建环境；安装 `mitom` 入口；运行 `mitom setup`（下载 MITOS2 的 refseq89m、准备 tRNA 模型、写配置）。不往系统 Python 装包。

### 2.2 Docker

```bash
docker build -t mitom .
docker run --rm -v "$PWD":/data mitom mitom run -i /data/genomes.gb -o /data/out
```

### 2.3 数据与配置在哪

| 内容 | 默认位置 |
|---|---|
| MITOS2 参考数据 refseq89m | `$MITOS_REFDIR` 或 `~/mitos_refs`（`mitom setup --refdir`） |
| MitoM-ML tRNA 模型 | `$MITOM_DATA` 或 `~/.local/share/mitom` |
| 配置文件 | 仓库内 `mito_config.json`（`mitom doctor --write-config`） |

### 2.4 可选软件

| 软件 | 用途 |
|---|---|
| MITOS2 2.1.10 | `--engine mitos2`、`--rescue yes`；`mamba env create -f env/mitos.yml` |
| MitoZ 3.6 | `--engine mitoz`；自行安装 |
| minimap2 | 圈图测序深度环（可选，装上后代替 bwa） |

体积对比见 [BENCHMARK.md § 安装体积](BENCHMARK.md#6-安装体积与用时参考)。

---

## 3. 五分钟上手

```bash
# GenBank（推荐：自带分类）
mitom run -i my_genomes.gb -o out

# 自己组装的 FASTA
mitom run -i contigs.fasta -o out --clade Mollusca --genetic_code 5

# 已有 MitoZ / MITOS2 结果：不重注，只整理
mitom import -i mitoz_or_mitos2_results/ -o out --clade Mollusca --genetic_code 5

# 核 rDNA / 组蛋白
mitom nuc annotate -i scaffolds.fasta -o nuc_out
```

**先看这些输出**

| 路径 | 内容 |
|---|---|
| `out/summary.tsv` | 每条一行总览（基因数、cox1 起点、NCBI ERROR） |
| `out/annotation/` | 最终 GenBank、序列、圈图 |
| `out/submission/` | NCBI `.fsa` / `.tbl` / `.sqn` / `.val` |
| `out/qc/` | 各步报告 |
| `out/logs/run.log` | 完整日志 |

用仓库数据跟做：[EXAMPLES.md](EXAMPLES.md)。

---

## 4. 怎么选命令

| 手上有 | 想要 | 用 |
|---|---|---|
| 还没有数据 | 某类群的 NCBI 线粒体全基因组 | `mitom fetch --taxon Annelida -o fetch_out` |
| GenBank（`.gb`） | 统一重注 | `mitom run -i … -o out` |
| FASTA | 注释 | `mitom run -i … -o out`（可写 `--clade`） |
| MitoZ / MITOS2 结果 | **保留**原基因，只整理 | `mitom import -i … -o out --clade …` |
| MitoZ `.gbf` | **丢掉**原基因、全部重注 | `mitom run -i sample.gbf -o out` |
| 已注释 GenBank | 只要边界 / 命名 | `mitom curate` |
| 组装 contig | 核 rDNA / 组蛋白 | `mitom nuc annotate` |
| 出问题 | 定位哪一步失败 | `mitom debug out` |

**`run` 与 `import`**：`run` 只用序列并全部重注，原注释留在 `input/` 供对比；`import` 保留已有基因，只做定向、补漏、校正、提交与作图。

**为何推荐 GenBank**：自带分类与文献 / 标本信息；FASTA 靠 blastx 推断门与密码，小门请手写 `--clade` 与 `--genetic_code`。

---

## 5. 每个命令：输入 · 命令 · 输出

参数名在各命令里一致：`-i` 输入、`-o` 输出、`--threads`、`--jobs`、`--refdir`。全部参数见 [OPTIONS.md](OPTIONS.md)。每个命令可加 `-h`。

### 5.1 `mitom run`（完整流程，最常用）

| | |
|---|---|
| **做什么** | 读入 → 去重 → 注释 → 定向 → 补漏 → 内含子 → 校正 → ORF / NCR → NCBI 文件 → 对比 → 圈图 |
| **输入** | GenBank 或 FASTA；单文件或目录；可多门混合 |
| **输出** | 见 [§ 6 输出目录](#6-输出目录) |

```bash
mitom run -i fetch_out/Annelida.gb -o out
mitom run -i contigs.fasta -o out --clade Mollusca --genetic_code 5
mitom run -i genomes/ -o out --threads 4 --jobs 8
mitom run -i contigs.fasta -o out --clade Mollusca --genetic_code 5 --engine mitoz
```

常用参数：`--clade`、`--genetic_code`、`--engine native|mitoz|mitos2`、`--sbt`、`--infq`、`--keep_intermediate yes`、`--debug`。12 步与报告文件见 [§ 5.12](#512-mitom-run-的-12-步)。

### 5.2 `mitom import`（整理已有注释）

| | |
|---|---|
| **做什么** | 保留输入里的基因；只做定向、补漏、校正、提交、作图 |
| **输入** | MitoZ 目录 / `.gbf`；MITOS2 目录（需 `--fasta`）；已注释 GenBank |
| **输出** | 与 `run` 相同布局 |

```bash
mitom import -i mitoz_results/ -o out --clade Mollusca --genetic_code 5
mitom import -i mitos2_result/ --fasta genome.fasta -o out --clade Mollusca --genetic_code 5
```

MITOS2 专页：[MITOS2.md](MITOS2.md)。跟做：[EXAMPLES.md 例 5](EXAMPLES.md#例-5-只整理已有注释import)。

### 5.3 `mitom fetch`（从 NCBI 下载）

| | |
|---|---|
| **做什么** | 按类群或 Entrez 检索式下载完整线粒体基因组 |
| **输入** | `--taxon Annelida`（或 taxid）、或 `--query '…'`；需 `--email` 或 `NCBI_EMAIL` |
| **输出** | `<名称>.gb`、`fetch_report.tsv`、`query.txt`、`fetch.log`；`--fasta` 时另有 FASTA |

```bash
mitom fetch --taxon Annelida -o fetch_out --email you@example.com
mitom fetch --taxon Annelida -o fetch_out --refseq_only
mitom fetch --taxon Annelida -o fetch_out --dry_run
```

同一目录可重复运行，只补新记录。下载后建议 `mitom run`（会自动去重）。

### 5.4 `mitom dedup`

| | |
|---|---|
| **做什么** | 同一基因组只留一条（`NC_`、新版本、相同序列留一条） |
| **输入** | GenBank 和 / 或 FASTA |
| **输出** | `dedup.gb`、`dedup.fasta`、`dedup_report.tsv` |

```bash
mitom dedup -i fetch_out/Annelida.gb -o dedup_out
```

### 5.5 `mitom gb2fasta` / `mitom gb2tbl`

**gb2fasta**

| | |
|---|---|
| **输入** | GenBank 文件或目录 |
| **输出** | FASTA；`--split` 时每条一个文件；表头带 `topology=` |

```bash
mitom gb2fasta -i sample.gbf -o fasta_out --split
```

**gb2tbl**（已注释 → NCBI `.tbl` + `.fsa`）

| | |
|---|---|
| **输入** | GenBank / EMBL 文件或目录，或管道 `-i -` |
| **输出** | `.tbl` + `.fsa`；`--validate` 时另有 `.sqn` / `.val` |

```bash
mitom gb2tbl -i annotation/ -o ncbi --validate
```

### 5.6 `mitom annotate` / `mitom curate` / `mitom intron` / `mitom trna`

| 命令 | 输入 | 输出要点 | 示例 |
|---|---|---|---|
| `annotate` | GenBank / FASTA | 工作目录式布局（`04.final_gb/` 等），无 round-2 整理 | `mitom annotate -i NC_000931.gb -o out_ann` |
| `curate` | 已注释 GenBank | 修整边界与命名；`curation_report.tsv` | `mitom curate -i ann/ -o curated --genetic_code 5` |
| `intron` | 已注释目录 | 加内含子；`intron_report.tsv` | `mitom intron -i ann/ -o intron_out --genetic_code 4` |
| `trna` | FASTA / GenBank | 多方法 tRNA 比较 | `mitom trna -i sample.gb -o trna_out --genetic_code 5` |

### 5.7 `mitom rna`（转录组证据）

| | |
|---|---|
| **输入** | 已注释线粒体 GenBank（`-g`）+ FASTQ，或已有 SAM |
| **输出** | `rna_features.tsv`、`rna_polya_sites.tsv`、`rna_depth.tsv.gz`；`--refine`：`rna_refined.gb`、`rna_end_corrections.tsv` |

```bash
mitom rna -g annotation/OQ078742.1.gb -1 R1.fq.gz -2 R2.fq.gz -o rna --threads 8
mitom rna -g annotation/OQ078742.1.gb -o rna --refine --run SRR26830859      # 证据写回注释 / evidence -> GenBank
```

比对默认 bwa mem（`--mapper bwa|minimap2|bowtie2`）。清洗时不要用 fastp `-x` / `-y` / `-D`（会毁掉 polyA 证据）。
转录组必须来自同一基因组项目（最好同一标本）；reads 与基因组差异 > 2 % 会警告。数字见
`benchmark/results/rnaseq_2026-11-10/`、`rna_mappers_2026-11-10/`。
Mapper bwa mem by default; RNA-seq must come from the same genome project (ideally the same specimen).

**`mitom rna` 与 `mito_rnaseq_annot.py` 的区别 / the two RNA-seq tools**

| | `mitom rna` | `python3 mito_rnaseq_annot.py` |
|---|---|---|
| 做什么 / what | 量证据：表达、poly(A) 3′ 端、文库链向；`--refine` 只改 3′ 端 / measures evidence; `--refine` changes 3′ ends only | 用转录组注释：质控、内含子（分段 reads）、重叠终止、3′ 端、RNA 与基因组不一致处、ORF 表达，每个判断标 confirmed / corrected / flagged / not_testable / annotates with RNA-seq, every decision labelled |
| 输出 / output | 证据表 / evidence tables | `rnaseq_evidence.tsv`、`rnaseq_annotated.gb` |
| 是否 `mitom run` 的一步 / part of `mitom run` | 否（单独命令）/ no | 否（单独脚本）/ no |

两者都不改序列；RNA 与基因组不一致处只标记（同一标本 ≥ 80 % reads 不一致多为组装碱基错误，用 HiFi 核实）。
规则与评测：`benchmark/results/annelid_rna_2026-11-10/README.md`。Neither changes the sequence; RNA-genome
mismatches are flagged only.

### 5.8 `mitom nuc`（核标记）

| | |
|---|---|
| **输入** | 组装 FASTA 或已注释 GenBank |
| **输出** | `nuc_features.tsv`、`consensus/`、`submission/` |

```bash
mitom nuc annotate -i scaffolds.fasta -o nuc_out
mitom nuc annotate -i scaffolds.fasta -o nuc_out --genes rdna
```

### 5.9 `mitom plot` / `mitom validate` / `mitom rescue` / `mitom rrna-hmm`

| 命令 | 输入 | 输出 |
|---|---|---|
| `plot` | 已注释 GenBank；可选 `--fq1/--fq2` | `.png` / `.pdf`（`docs/images/` 有样式样例） |
| `validate` | 已注释目录 | 提交文件与 `.val` |
| `rescue` | GenBank **文件夹**（原地改） | 补回漏基因；需 MITOS2 |
| `rrna-hmm search` | GenBank **文件夹** | 补搜 12S / 16S |

```bash
mitom plot -i out/annotation -o maps
mitom validate -i out/annotation -o ncbi_check
```

### 5.10 `mitom debug` / `doctor` / `setup`

```bash
mitom doctor              # 程序是否齐
mitom debug out           # 哪一步、哪条记录、哪个程序失败
mitom run -i in.gb -o out --debug   # 保留中间文件与完整命令输出
mitom setup               # 安装时已跑；可单独重下参考数据
```

### 5.11 其它

| 命令 | 用途 |
|---|---|
| `mitom reorient` | 只统一方向与 cox1 起点 |
| `mitom check` | 审核已有注释，只出报告 |
| `mitom compare` | 原注释 vs MitoM，逐基因 |
| `mitom transfer_metadata` / `mitom replace_organism` | 迁回表头 / 改物种名（`run` 通常已做） |

### 5.12 `mitom run` 的 12 步

| 步 | 做什么 | 报告 |
|---|---|---|
| 1 | 读入与去重 | `qc/dedup_report.tsv` |
| 2 | 第一轮注释；核对环状 / 线状 | `qc/topology_decision.tsv` |
| 3 | 环状记录转到正链、cox1 起始 | `qc/oriented_by_fallback.txt` |
| 4 | 需要时再注释 | `qc/reannotation_decision.tsv` |
| 5 | 写 GenBank，复制原表头 | `qc/metadata_transfer_log.tsv` |
| 6 | 补漏、rRNA HMM、内含子 | `qc/rescue_report.tsv`、`intron_*.tsv` |
| 7 | 边界与 NCBI 命名 | `qc/curation_report.tsv` |
| 8 | 未知 ORF 与 NCR | `qc/orf_report.tsv` |
| 9 | table2asn 与提交文件 | `submission/` |
| 10 | 与原注释对比 | `qc/original_vs_mitom.tsv` |
| 11 | 基因组图 | `annotation/plots/` |
| 12 | 整理目录 | 见下节 |

环状 / 线状、类群自动判断的规则与评测：[METHODS.md](METHODS.md) 对应节；本手册不展开。

---

## 6. 输出目录

默认 `--output_layout grouped`。`--output_layout classic` 为旧平铺布局。

```
out/
├── README.txt
├── summary.tsv              # 先看这个
├── annotation/              # *.gb、all_records.gb、genes.fna、proteins.faa、plots/
├── submission/<登录号>/     # .fsa .tbl .sqn .val
├── qc/                      # 各步报告
├── input/                   # 输入原样（GenBank 时）
├── logs/                    # run.log、run_params.json
└── work/                    # 仅 --keep_intermediate yes
```

**提交 NCBI**：`submission/` 里 ERROR 为 0 的记录可提交。正式提交请加 `--sbt` 提交者模板。NCBI 收 `.sqn`，或 BankIt 上传 `.fsa` + `.tbl`。

---

## 7. 输入要点（摘要）

### 7.1 格式

- GenBank：`.gb .gbk .gbf .gbff .genbank`（单文件或目录）
- FASTA：单文件或目录；表头可写 `topology=circular|linear`

### 7.2 类群与遗传密码

默认 `--clade auto`。GenBank 用记录分类；FASTA 用 blastx + 密码子核验。多门时 `run` 自动分组。`import` / `annotate` 一次只能一个门。

常用 `--clade`：`Chordata`（密码 2）、`Arthropoda` / `Mollusca` / `Annelida-segmented-worms`（5）、`Echinodermata`（9）、`Porifera-sponges`（4）、`Cnidaria` / `Platyhelminthes-flatworms`（需给密码）、`Animal`（其它门）。

评测表：[BENCHMARK.md § 类群](BENCHMARK.md#2-类群与遗传密码---clade-auto)。

### 7.3 去重

`NC_` 优先、新版本优先、相同序列只留一条。报告：`qc/dedup_report.tsv`。

---

## 8. 排错

```bash
mitom doctor
mitom debug out
```

| 现象 | 处理 |
|---|---|
| `annotate` / `import` 报多个门 | 改用 `mitom run`，或显式 `--clade` |
| 门 uncertain | 看 `qc/clade_decision.tsv`；知道类群就手写 |
| 输出目录被锁 | 等当前运行结束，或换 `-o` |
| 输入错误一行 `ERROR:`、退出码 2 | 空文件、截断下载、蛋白序列等；`--debug` 看细节 |
| `ncbi_errors` ≠ 0 | 看 `submission/<id>/<id>.val` |
| 很慢 | `--threads` × `--jobs` ≤ CPU；内含子约 10–20 s / 条 |

更多行见历史手册中的完整 FAQ（CHANGELOG 与 issue 跟踪）。报告问题时附上 `mitom doctor` 与 `mitom debug` 输出。

---

## 9. 引用

请同时引用 MitoM 所调用的软件（MITOS2 参考数据、MiTFi、MFannot、Infernal、BLAST+、Exonerate、HMMER、table2asn、Biopython 等）。见 [CITATION.cff](../CITATION.cff) 与仓库 README「引用」一节。

---

## 附录

- 全部参数：[OPTIONS.md](OPTIONS.md)
- 重建 `models/`：[DATASETS.md](DATASETS.md)
- 评测摘要：[BENCHMARK.md](BENCHMARK.md)
