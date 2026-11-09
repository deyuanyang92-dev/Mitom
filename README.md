# MitoM：动物线粒体基因组注释 / Annotation of animal mitochondrial genomes

> **版本 / Version**：v0.0.1（测试版 / test release）。
>
> **说明**：MitoM 目前仍在测试中，是本团队内部使用的软件，尚未正式发布。结果（尤其是内含子、tRNA 和基因边界）请人工核对后再用于发表或提交 NCBI。
> MitoM is still being tested and is used inside our group only; check the results before publishing or submitting them.

MitoM（命令 `mitom`）找出动物线粒体基因组的 13 个蛋白基因、2 个 rRNA、22 个 tRNA 和内含子，统一方向与起点，
按 NCBI 规范修整边界和基因名，生成可直接提交 NCBI 的文件和环状基因组图；`mitom nuc` 在全基因组组装中注释核 rDNA
（18S–ITS1–5.8S–ITS2–28S）和组蛋白基因。完整使用手册见 **[docs/USAGE_zh.md](docs/USAGE_zh.md)**。

MitoM (command `mitom`) annotates the 13 protein genes, 2 rRNAs, 22 tRNAs and the introns of animal mitochondrial
genomes, orients them, curates gene ends and names to NCBI rules, and writes NCBI submission files and genome maps;
`mitom nuc` annotates nuclear rDNA and histone genes in whole-genome assemblies.

| 精度 / Accuracy（开发期基准，2026-10-27 / 10-29，与 NCBI 比较；不是论文的冻结评测 [benchmark/frozen/PLAN.md](benchmark/frozen/PLAN.md)；细节见 [docs/METHODS.md](docs/METHODS.md)） / development benchmarks, not the frozen evaluation | |
|---|---|
| 蛋白基因两端与 NCBI 完全一致 / protein genes, both ends exact | 1,050 条多类群 RefSeq：24 个门全部高于 MITOS2 |
| NCBI table2asn ERROR | 2,497 条基准记录 / benchmark records（2026-10-27；1,050 多类群 RefSeq、669 环节动物、577 科未参与训练的无脊椎动物、188 脊椎动物、13 锚定基因组）：0（此前 / before：199） |
| tRNA（未参与训练的类群）/ tRNAs (families unseen in training) | 召回率高于 MITOS2，精确率持平或更高 |
| 内含子 / introns | 逐个核实的 NCBI 内含子金标准：默认全部一致，被测记录不进参考库时少数不一致，无内含子记录无假阳性（数字见 [METHODS 内含子表](docs/METHODS.md#intron-注释--introns-mito_intronpy-mito_intron_refpy)） / gold-standard introns: see the METHODS intron table |
| 核 rDNA / 组蛋白 / nuclear markers | 全基因组组装几秒到几分钟；准确度与用时见 [METHODS 核基因节](docs/METHODS.md#核基因标记--nuclear-markers-rdna-and-histones-mito_nucpy) / seconds to minutes per assembly; accuracy and time in METHODS |

## 快速上手 / Quick start

**1. 安装（一条命令）/ Install (one command)**：Linux 或 Windows WSL2；macOS 见 `install.sh` 开头说明，或用 Docker：

```bash
git clone https://github.com/deyuanyang92-dev/Mitom.git MitoM && cd MitoM
bash install.sh            # 一个 conda 环境 mitotools（Python 包和全部程序）+ 命令 mitom + mitom setup（参考数据、配置）；
                           # 没有 conda 时自动装 Miniforge；~/.local/bin 自动写入 PATH（~/.bashrc），之后开一个新终端
mitom doctor               # 检查：第一行是正在用的版本和文件夹；所有必需项都应显示 OK
```

**更新 / Update**：`cd MitoM && git pull && bash install.sh`（已有的环境只更新，可重复运行）。只 `git pull` 而不重跑
`install.sh` 时，命令 `mitom` 可能仍指向旧的副本或别的 Python；`mitom doctor` 会用 `!!` 指出并给出修复命令。
Update: `git pull && bash install.sh`; `mitom doctor` flags a `mitom` command that runs an old copy or another Python.

全部依赖都在一个 conda 环境 `mitotools` 里（`env/mitotools.yml`，Python 包也在内），不往系统 Python 装任何东西；`mitom` 用这个环境的 Python 和程序。
All dependencies, the Python packages included, live in one conda environment (`env/mitotools.yml`); nothing is installed into the system Python.

Docker：`docker build -t mitom . && docker run --rm -v "$PWD":/data mitom mitom run -i /data/genomes.gb -o /data/out`
（bioconda 配方在 `packaging/bioconda/`；发布后可 `mamba create -n mitom -c conda-forge -c bioconda mitom && conda activate mitom && mitom setup`）

**安装大小 / Install size**（2026-10-28，Linux x86_64 全新安装实测：conda 环境 + 参考数据 refseq89m）：

| 软件 / tool | 下载 / download | 安装后 / installed | conda 包 / packages |
|---|---|---|---|
| MitoM | 0.48 GB | 1.8 GB（首次运行再解压 0.19 GB tRNA 模型） | 157 |
| MITOS2 2.1.10 | 0.75 GB | 2.9 GB | 185 |
| MitoZ 3.6 | 1.37 GB | 4.5 GB | 380 |

MitoM 用电脑上已有的 Java（MiTFi 用，任何 ≥ 8 的版本）；没有时 install.sh 另装 openjdk 8（下载约 0.09 GB）。MITOS2 一行含它的参考数据 refseq89m，MitoZ 一行只是它的 conda 环境。conda 下载的安装包缓存可用 `mamba clean -t` 删除。

MitoZ、MITOS2 不随 MitoM 安装；装了它们时，MitoM 可以调用它们注释或导入它们的结果（见下）。

**2. 运行 / Run**

```bash
# 从 NCBI 下载某类群的全部完整线粒体基因组（RefSeq + GenBank；--refseq_only 只要 NC_；--query 自定义检索式）
mitom fetch --taxon Annelida -o fetch_out            # -> fetch_out/Annelida.gb + fetch_report.tsv（可重复运行，只补新记录）

# GenBank 输入（推荐）：每条记录都重新注释（不管原来有没有注释），结果在 out/；run 做了哪些事见下面第 3 点
mitom run -i fetch_out/Annelida.gb -o out

# FASTA 输入（没有分类信息）：不写 --clade 时用 blastx 比对参考蛋白，自动判断每条序列的门和遗传密码
mitom run -i contigs.fasta -o out
# 知道类群时也可以写明（跳过自动判断）
mitom run -i contigs.fasta -o out --clade Mollusca --genetic_code 5

# 用 MitoZ 或 MITOS2 注释（需自己另装；之后照常方向调整、第二轮、校正、NCBI 文件、圈图）
mitom run -i contigs.fasta -o out --clade Mollusca --genetic_code 5 --engine mitoz      # 或 --engine mitos2

# 已有 MitoZ / MITOS2 结果或任何注释好的 GenBank：不再注释，直接方向调整、端点校正与 NCBI 命名、NCBI 提交文件、圈图
mitom import -i mitoz_results/ -o out --clade Mollusca --genetic_code 5               # MitoZ 结果目录或 .gbf
mitom import -i mitos2_result/ --fasta genome.fasta -o out --clade Mollusca --genetic_code 5   # MITOS2 结果目录
```

MITOS2 用户 / MITOS2 users：MitoM 把 MITOS2 结果批量转成 GenBank、统一方向、修整端点、生成 NCBI 文件，见 / see **[docs/MITOS2.md](docs/MITOS2.md)**。

`--clade` 可选：`Chordata Arthropoda Echinodermata Annelida-segmented-worms Bryozoa Mollusca Nematoda
Nemertea-ribbon-worms Porifera-sponges Cnidaria Platyhelminthes-flatworms Chaetognatha-arrow-worms Animal`
（`Animal` = 全动物蛋白库，用于其他门）。常用：`--threads 8 --jobs 4`（每个任务线程数 / 同时任务数）、
`--sbt template.sbt`（NCBI 提交者模板）。`--threads` 默认 `min(8, CPU 数)`，`--jobs` 默认 `CPU 数 // threads`。`--clade auto`（默认）逐条判断：GenBank 记录用自己的分类信息和 CDS 的 `transl_table`；FASTA 和没有分类信息的记录用 blastx 比对 refseq89m 的 8 个保守蛋白，取最近参考基因组的门和遗传密码，再用 AGR / ATA / AAA / TGA 等密码子核验（没有近亲时按 `Animal` 和密码子判断的遗传密码）。一次输入含多个门或多个遗传密码时，`run` 按组分别注释再合并到同一个输出目录；每条记录的判断和依据见 `qc/clade_decision.tsv`。

**3. 运行过程与结果 / What you see and what you get**

**`mitom run` 的 11 步 / The 11 steps of `mitom run`**（与终端的 `[Step k/11]` 一致；每步的细节和报告文件见 [使用手册第 4 节](docs/USAGE_zh.md#4-mitom-run-的流程)）

| 步 / step | 终端显示 / terminal | 做什么 |
|---|---|---|
| 1 | Read input and remove duplicate records | 读入；同一个基因组只留一条（优先 `NC_`、新版本） |
| 2 | Round-1 annotation | 13 个蛋白基因、2 个 rRNA、22 个 tRNA；核对环状 / 线状 |
| 3 | Re-orient | 环状记录转到正链，从 cox1 开始 |
| 4 | Check the annotation after re-orientation (round 2) | 转动后受影响的记录再注释一次 |
| 5 | Write the annotation; copy the original GenBank header back | 写 GenBank，写回物种、文献、标本信息 |
| 6 | Add what the core missed | 补 12S / 16S、内含子；`--rescue yes` 时用 MITOS2 补基因 |
| 7 | Boundary rules for protein genes and rRNAs, NCBI gene names | 按 NCBI 规范修整端点和基因名；标出非编码区（NCR，`--ncr no` 不标） |
| 8 | NCBI check (table2asn) and submission files | 生成 `submission/` |
| 9 | Compare with the original annotation | 逐基因对比原注释（GenBank 输入） |
| 10 | Genome maps | 每条一张圈图（`annotation/plots/`） |
| 11 | Organise the output | 整理成下面的目录结构 |

运行时终端按步骤显示（每步：做什么、多少条、结果、用时、结果文件；注释阶段每 30 秒报告进度），完整细节在 `logs/run.log`：

```text
[Step 1/11] Read input and remove duplicate records
    GenBank records read: 802
    removed 166: GenBank record that a RefSeq NC_ record was made from (NC_ kept)
    -> 636 records to annotate; reasons per record: qc/dedup_report.tsv  (5 s)
[Step 2/11] Round-1 annotation (protein genes, rRNAs, tRNAs)
    round-1 annotation: 120/636 (18%), 6 min 02 s so far, about 26 min left
...
[Step 4/11] Check the annotation after re-orientation; annotate again where it is offset (round 2)
    checked 636: 41 need round 2 (...); 595 keep the shifted round-1 annotation
...
[Step 9/11] Compare with the original annotation
...
[Step 10/11] Genome maps (circular plots)
```

输出目录（`out/README.txt` 逐项说明；`--output_layout classic` 用 2026-11-09 以前的布局）：

| 文件 / 目录 / file / folder | 内容 / contents |
|---|---|
| `summary.tsv` | 每条记录一行：基因数、cox1 起始、用了第几轮注释、是否迁回原 header、NCBI ERROR 数（先看这个） |
| `annotation/<登录号>.gb` | MitoM 最终注释，每条一个文件（cox1 起始、正链、校正、NCBI 命名） |
| `annotation/all_records.gb` | 同上，合并成一个文件 |
| `annotation/genes.fna`、`proteins.faa` | 蛋白基因和 rRNA 的核苷酸序列、蛋白序列；序列名 `登录号\|基因\|物种`，可直接比对、建树 |
| `annotation/plots/<登录号>.png`、`.pdf` | 基因组圈图（`mito_plot.py`，matplotlib）：每条链一个环道（正链在外、负链在内），基因按类别着色、名字沿圆弧书写，向内依次为坐标轴、GC 含量；给出 `--infq` 时再加一圈测序深度（装有 minimap2 时用它，否则用 bwa mem）；`--plot_style mitom` 用原来的样式（另有 GC skew 环），`--plot no` 不画 |
| `submission/<登录号>/` | 每条的 NCBI 提交文件 `.fsa/.tbl/.sqn` 与 table2asn 结果 `.val`；`submission/submission_summary.tsv` ERROR 为 0 的可直接提交 |
| `qc/` | 核对与判断：`original_vs_mitom.tsv`（原注释 vs MitoM，每个基因一行）、`dedup_report.tsv`（去重依据）、`reannotation_decision.tsv`（第二轮与否及原因）、`curation_report.tsv`（每处校正）等 |
| `input/<登录号>.gb` | 输入记录原样（含 NCBI 原注释），文件名与 `annotation/` 对应，便于前后对比 |
| `logs/` | `run.log`（全部过程与调用的命令）、`run_params.json`（全部参数）、`intron_fix.log`、出错记录的日志 |

**4. 全部子命令 / Commands**（每个命令加 `-h` 看参数，如 `mitom run -h`，全部参数汇总在 [docs/OPTIONS.md](docs/OPTIONS.md)；每个命令的用法和示例见 [使用手册第 6 节](docs/USAGE_zh.md#6-每个命令的用法)，怎么选命令见 [第 5 节](docs/USAGE_zh.md#5-怎么选命令)）

| 命令 / command | 作用 / purpose |
|---|---|
| `mitom run` | 完整流程，最常用：从序列到最终注释、NCBI 提交文件和基因组图 |
| `mitom import` | 已有 MitoZ / MITOS2 结果或已注释的 GenBank：不重新注释，只做后面的整理（统一起点、修整、NCBI 文件、作图） |
| `mitom fetch` | 从 NCBI 下载某个类群的线粒体全基因组（GenBank 格式；加 `--fasta` 同时生成 FASTA） |
| `mitom dedup` | 只做“去重”：同一个基因组的多条记录只留一条（说明见第 3 点）；GenBank、FASTA 都可以 |
| `mitom gb2fasta` | GenBank（`.gb`、`.gbf` 都行）转 FASTA，表头带 `topology=`，可直接给 `mitom run` 用 |
| `mitom gb2tbl` | GenBank / EMBL 转 NCBI 提交用的特征表 `.tbl` + `.fsa`：多个文件 / 目录，全部合并、每个文件一对或每条记录一对；`--source_table` 批量填来源信息，缺的产物名自动补；`--validate` 用 table2asn 生成 `.sqn` 并检查 |
| `mitom reorient` | 只统一方向和起点（转到同一条链，从 cox1 开始） |
| `mitom annotate` | 与 `run` 相同，但不整理输出目录（给脚本或集群流程调用） |
| `mitom check` | 检查已有注释有没有问题，只出报告、不改文件 |
| `mitom compare` | 逐个基因比较原注释和 MitoM 的注释 |
| `mitom curate` | 只修整基因起止和基因名（NCBI 规范） |
| `mitom intron` | 只加内含子 |
| `mitom trna` | 用多种方法注释 tRNA 并相互比较 |
| `mitom rrna-hmm` | 补搜缺失的 12S / 16S rRNA |
| `mitom rescue` | 用 MITOS2 补回漏掉的基因（需另装 MITOS2） |
| `mitom validate` | 用 NCBI 的 table2asn 检查，生成提交文件 |
| `mitom plot` | 画环状基因组图（`--style circos` 默认，或 `--style mitom`） |
| `mitom nuc` | 在全基因组组装中注释核 rDNA（18S–28S）和组蛋白基因 |
| `mitom debug` | 出问题时读日志，告诉你哪一步、哪条记录、哪个程序出错 |
| `mitom transfer_metadata` / `replace_organism` | 把原记录的物种、文献等信息复制到新注释 / 修正物种名 |
| `mitom setup` / `mitom doctor` | 下载参考数据（安装时已自动运行）/ 检查需要的程序是否都装好 |

**注释审核报告 / Annotation audit report**：`mitom check -i genomes.gb -o check_out` 审核已有注释（只出报告，不改文件；
可用于查看 NCBI 下载记录的原注释问题，或检查最终结果），`check_report.tsv` 每条一行（ok / flagged），`check_details.tsv` 列出每个问题。ERROR（触发重注）：无注释、缺蛋白基因或
rRNA、基因重复、内部终止子、长度非 3 倍数且无 transl_except、基因长度异常、基因重叠 > 100 nt、tRNA 少于 22（`--min_trna`）；
WARNING：起始/终止密码子、cox1 方向、tRNA 种类、非标准 CDS。

**5. 出问题时 / When something goes wrong**

```bash
mitom doctor                         # 先查安装：版本和文件夹、Python、每个程序；有问题的行标 !!，并给出修复命令
mitom debug out                      # 读 out/logs/run.log：哪一步、哪条记录、哪个程序失败，错误信息，缺什么程序及安装命令，复现命令
mitom run -i in.gb -o out --debug    # 调试模式：记录环境（各程序版本、配置、数据目录）、每条命令的完整输出，保留中间文件
```
求助时请附上 `mitom doctor` 和 `mitom debug out` 的输出。/ When asking for help, include the output of both.

**文档 / Documentation**：完整使用手册 [docs/USAGE_zh.md](docs/USAGE_zh.md)（安装、输入、流程、每个命令、输出、调试、常见问题）；
全部命令与参数 [docs/OPTIONS.md](docs/OPTIONS.md)；方法与评测 [docs/METHODS.md](docs/METHODS.md)；版本变化 [CHANGELOG.md](CHANGELOG.md)；评测数据 `benchmark/results/`。

## MitoM：声明、致谢与引用 / Statement, acknowledgements and citation

**MitoM**（命令名 `mitom`；M = Merge / Metazoa）是一条动物线粒体基因组注释流程。流程本身（下表"MitoM 自有步骤"）是作者多年研究中逐步优化形成的；基因识别由 MitoM 自己的注释核心 `mito_core.py` 完成，各环节用什么方法以 benchmark 结果为依据（`benchmark/`）。MitoZ、MITOS2 不随 MitoM 安装：装了它们时可以用 `--engine mitoz|mitos2` 让它们注释，或用 `mitom import` 读入它们已有的结果，之后的 MitoM 步骤相同；MITOS2 还可用于漏注救回（`--rescue yes`）。核心用到的 MITOS2 参考数据 refseq89m 由 `mitom setup` 下载（install.sh 会运行）。

| 环节 / part | 来源 / origin |
|---|---|
| NCBI 下载（`mito_fetch.py`，按类群 / 检索式）、GenBank 输入去重（保留 RefSeq）、注释审核报告（`mito_check.py`）、内部短 ID、批量与并行 | MitoM 自有步骤 |
| 方向调整（cox1 起始、统一正链，cox1 缺失时依次用 cox2/cox3/cob/nad1） | MitoM 自有步骤 |
| 方向调整后第二轮注释 | MitoM 自有步骤 |
| 元数据迁回（原 GenBank header）、`/organism` 修复 | MitoM 自有步骤 |
| 内含子（`mito_intron.py` 按 MFannot 方法从头识别 + `mito_intron_ref.py` 近缘位点迁移） | MitoM 自有步骤（从头识别方法来自 MFannot；位点库由 193 条核实的 NCBI 记录构建） |
| 漏注基因从 MITOS2 救回（`mito_rescue.py`，可选：`--rescue yes`，需另装 MITOS2） | MitoM 自有步骤 |
| 仍缺的 12S/16S 用按门建立的 rRNA HMM 补搜（`mito_rrna_hmm.py`，库 `models/rrna_hmm/`，nhmmer；只在常规方法找不到时启用） | MitoM 自有步骤（思路参考 EMMA 的 nhmmer rRNA 模型） |
| 蛋白基因与 rRNA 端点校正（`mito_curate.py`；起始密码子另由从 NCBI 记录学习的模型挑选，规则 P11）、NCBI table2asn 校验 | MitoM 自有步骤（规则依据 Fourdrilis 2018、Ghiselli 2021） |
| 环状 / 线状核对（`mitom_topology.py`、`mito_core.py --check_ends`：两端的基因跨接点再比对；跨接点的内含子用 protein2genome） | MitoM 自有步骤 |
| 非编码区（`mito_ncr.py`：≥ 100 nt 无注释的区段写成 NCR misc_feature） | MitoM 自有步骤 |
| 多方法 tRNA 比对（`mito_trna.py`） | MitoM 自有步骤 |
| 基因组圈图（`mito_plot.py`：基因、GC 含量，可选测序深度，内含子画斜线、非编码区画浅灰；默认样式借鉴 mitoPlotR） | MitoM 自有步骤（取代 MitoZ 的 circos 圈图） |
| GenBank / EMBL 转 NCBI 提交文件（`mito_tbl.py`，`mitom gb2tbl`） | MitoM 自有步骤（来源信息表、产物名表的做法借鉴 GB2sequin；管道输入、定义行借鉴 NCBI gbf2tbl / gbf2fsa） |
| 用 MitoZ / MITOS2 注释（`--engine mitoz\|mitos2`）或导入它们已有的结果（`mitom import`；`mito_engines.py`） | 可选：MitoZ、MITOS2 需自己安装 |
| 蛋白基因、rRNA、tRNA 识别 | 注释核心 `mito_core.py`（blastx + exonerate 蛋白基因，按门 rRNA profile HMM，MitoM-ML tRNA：分支协方差模型 + 学习的判别器与端点，未调用的 tRNA 家族由 MiTFi 补，脊椎动物取 MiTFi 端点） |

`mito_core.py` = MitoM 自己的核心（blastx 定位 + exonerate 定蛋白基因边界，参考蛋白来自 MITOS2 的 refseq89m，起始密码子按 MitoZ 的规则取比对上游最近的一个（比对缺少较多残基时取最接近预期位置的）；rRNA 用按门的 profile HMM（`--rrna_method cm` 改用 MITOS2 的协方差模型）；tRNA 用 MitoM-ML（MITOS2 tRNA 模型 + 220 个分支的协方差模型，ARWEN/ARAGORN 作证据，梯度提升判别器与端点模型；它未调用的家族用 MiTFi 补，脊椎动物取 MiTFi 端点；无 scikit-learn 时整个改用 MiTFi）），只需一个 conda 环境，每条记录约 1–2 分钟。

**MitoM** (command name `mitom`; M = Merge / Metazoa) is an annotation workflow for animal mitochondrial genomes. The workflow itself (the "MitoM steps" above: input deduplication, reorientation, round-2 annotation, metadata transfer, introns, rescue of missed genes, per-phylum rRNA HMM fallback, boundary curation, NCBI validation, genome maps) was developed by the author over years of research. Gene finding is done by MitoM's own annotation core `mito_core.py`, its methods chosen by benchmark results (`benchmark/`): blastx + exonerate protein genes, per-phylum rRNA profile HMMs, MitoM-ML tRNAs (clade covariance models with a learned caller and ends, MiTFi for the families it does not call and MiTFi ends for vertebrates); one conda environment. MitoZ and MITOS2 are not installed with MitoM: when they are installed, `--engine mitoz|mitos2` lets them annotate, or `mitom import` reads their existing results, and the MitoM steps that follow are the same; MITOS2 can also rescue missed genes (`--rescue yes`). The MITOS2 reference data refseq89m the core uses are downloaded by `mitom setup` (install.sh runs it).

**声明 / Statement**

- 本软件基于作者多年研究中的一些心得，整合了大家的优秀软件，并使用 Claude Code（Anthropic）开发。
- 作者不保留任何版权：本仓库代码以 CC0 1.0 放入公有领域（见 `LICENSE`），任何人都可以使用、修改、再发布，欢迎持续优化升级。
- 使用本软件时，请引用本 GitHub 仓库，其他没有任何要求。
- This software builds on the author's experience from years of research. It integrates excellent tools written by others and was developed with Claude Code (Anthropic).
- The author claims no copyright: the code in this repository is dedicated to the public domain under CC0 1.0 (see `LICENSE`). Anyone may use, modify and redistribute it, and further improvements are welcome.
- If you use it, please cite this GitHub repository. Nothing else is required.

**第三方软件 / Third-party tools**

注释核心调用的程序（BLAST+、exonerate、Infernal、HMMER、Java、table2asn 等）和测序深度环用的 bwa（装有 minimap2 时自动改用它，快约 8 倍）、Python 包（Biopython、scikit-learn、matplotlib）由 `install.sh` 装入同一个 conda 环境，各自保留原有许可证。随仓库分发的第三方文件：MITOS2 的 MiTFi（`models/mitfi/`，MIT，见其中的 LICENSE）；`mitom setup` 下载 MITOS2 的参考数据 refseq89m（蛋白集与 tRNA/rRNA 模型，MIT）。MitoZ（GPL-3.0）、MITOS2 可选，不随本仓库分发。

The programs the core runs (BLAST+, exonerate, Infernal, HMMER, Java, table2asn ...) and bwa (read depth ring; minimap2 is used instead when installed, about 8x faster) and the Python packages (Biopython, scikit-learn, matplotlib) are installed into one conda environment by `install.sh`; they keep their own licences. Third-party files shipped here: MiTFi from MITOS2 (`models/mitfi/`, MIT); `mitom setup` downloads the MITOS2 reference data refseq89m (protein sets, tRNA / rRNA models, MIT). MitoZ (GPL-3.0) and MITOS2 are optional and not distributed here.

**致谢 / Acknowledgements**

感谢 MitoZ、MITOS2/MITOS、MiTFi、MFannot、Infernal、BLAST+、GeneWise、Exonerate、NCBI table2asn、Biopython、mitoPlotR（圈图样式）、GB2sequin（提交文件转换的做法）等软件的作者。没有他们的工作，就没有本软件。

We thank the authors of MitoZ, MITOS2/MITOS, MiTFi, MFannot, Infernal, BLAST+, GeneWise, Exonerate, NCBI table2asn, Biopython, mitoPlotR (map layout), GB2sequin (submission file conversion) and the other tools listed below. This software would not exist without their work.

**引用 / Citation**

引用本软件的同时，请一并引用实际用到的软件 / When citing this software, please also cite the tools it ran:

- 本软件 / This software: https://github.com/deyuanyang92-dev/Mitom (`CITATION.cff`)
- MitoZ（`--engine mitoz`、方法参考 / methods reference）: Meng G, et al. 2019. Nucleic Acids Res. doi:10.1093/nar/gkz173
- MITOS2: Donath A, et al. 2019. Nucleic Acids Res. doi:10.1093/nar/gkz833
- MITOS: Bernt M, et al. 2013. Mol Phylogenet Evol. doi:10.1016/j.ympev.2012.08.023
- MiTFi (tRNA models; MitoM's core and MITOS2): Jühling F, et al. 2012. Nucleic Acids Res. doi:10.1093/nar/gkr1131
- Infernal: Nawrocki EP, Eddy SR. 2013. Bioinformatics. doi:10.1093/bioinformatics/btt509
- BLAST+: Camacho C, et al. 2009. BMC Bioinformatics. doi:10.1186/1471-2105-10-421
- 内含子 / introns（`mito_intron.py`）— MFannot: Lang BF, et al. 2023. Front Plant Sci. doi:10.3389/fpls.2023.1222186; Exonerate: Slater GSC, Birney E. 2005. BMC Bioinformatics. doi:10.1186/1471-2105-6-31
- 端点校正依据 / boundary curation rules（`mito_curate.py`）: Fourdrilis S, et al. 2018. Sci Rep 8:17964. doi:10.1038/s41598-018-36428-7; Ghiselli F, et al. 2021. Phil Trans R Soc B. doi:10.1098/rstb.2020.0159
- 可选 tRNA 工具 / optional tRNA tools（`mito_trna.py`）: ARWEN doi:10.1093/bioinformatics/btm573; ARAGORN doi:10.1093/nar/gkh152; tRNAscan-SE 2.0 doi:10.1093/nar/gkab688
- 基因组圈图 / genome maps（`mito_plot.py`）— matplotlib: Hunter JD. 2007. Comput Sci Eng 9:90-95. doi:10.1109/MCSE.2007.55; 测序深度 / read depth — BWA-MEM: Li H. 2013. arXiv:1303.3997; minimap2: Li H. 2018. Bioinformatics 34:3094-3100. doi:10.1093/bioinformatics/bty191
- Biopython: Cock PJA, et al. 2009. Bioinformatics. doi:10.1093/bioinformatics/btp163
- NCBI table2asn: https://www.ncbi.nlm.nih.gov/genbank/table2asn/

## 开发期评测 / Development benchmarks（2026-10）

开发期间的基准，不是论文用的冻结评测（[benchmark/frozen/PLAN.md](benchmark/frozen/PLAN.md)）；方法与更多结果见 [docs/METHODS.md](docs/METHODS.md)。测出来的数字只在表里。
Development benchmarks, not the frozen evaluation of the paper; methods and more results in docs/METHODS.md. Measured numbers are only in the tables.

### 六个工具（2026-10-22）/ Six tools

数据 / data：24 条 RefSeq 记录 / records（8 门 / phyla × 3；MITOS1：8）；`benchmark/results/six_tools_2026-10-22/`，`docs/benchmark_six_tools_2026-10-22.md`

裸序列输入，与 NCBI 注释逐基因比较，并用 table2asn 校验。
Bare sequences, compared gene by gene with the NCBI annotation and checked with table2asn.

四个工具都完成的 17 条 / the 17 records all four tools finished:

| 指标 / metric | MitoM | MitoZ | MITOS2 | MitoFinder |
|---|---|---|---|---|
| 蛋白基因两端完全一致 / PCG both ends exact | 84.9% | 78.4% | 50.9% | 59.2% |
| rRNA 两端在 3 nt 内 / rRNA within 3 nt | 59.4% | 28.1% | 25.0% | 28.1% |
| tRNA 找回 / tRNA found | 96.8% | 92.5% | 96.0% | 96.5% |
| NCBI ERROR 条数 / NCBI ERROR count | 0 | 12 | 17 | 17 |

全部记录的用时和内存 / time and memory, all records:

| 工具 / tool | 每条用时中位数 / median time | 峰值内存中位数 / median peak RAM | 峰值内存最大 / max peak RAM |
|---|---|---|---|
| MitoM | 70 s | 182 MB | 270 MB |
| MitoZ | 249 s | 492 MB | 719 MB |
| MITOS2 | 385 s | 953 MB | 1929 MB |
| MitoFinder | 204 s | 238 MB | 283 MB |
| MITOS1 | 约 41 min | 775 MB | |

留一与失败 / leave-out and failures:

| 项目 / item | 结果 / result |
|---|---|
| 去掉该记录自己的参考蛋白后，MitoM 蛋白基因两端一致 / MitoM PCG both ends, reference of that record left out | 82.6% |
| MitoM 不低于 MitoZ 的门 / phyla where MitoM ≥ MitoZ | 8 / 8 |
| 其中最接近的：棘皮动物 MitoM 对 MitoZ / closest: Echinodermata, MitoM vs MitoZ | 94.9% vs 92.3% |
| MITOS1 自身崩溃 / MITOS1 crashed | 4 / 8 |
| MITOS1 成功的记录上蛋白基因两端完全一致 / PCG both ends exact on MITOS1 successes | 0 |

### 多类群（2026-10）/ Metazoa

数据 / data：24 条 RefSeq 记录 / records（环节 Annelida、节肢 Arthropoda、脊索 Chordata、刺胞 Cnidaria、棘皮 Echinodermata、软体 Mollusca、线虫 Nematoda、扁形 Platyhelminthes × 3）；MitoZ：18（无刺胞、扁形 clade / no clade）；MitoFinder：2 条失败，无近缘参考 / failed, no close reference；`benchmark/results/metazoa_mini_2026-10/`，`docs/benchmark_metazoa_2026-10.md`（含安装、功能、便捷性对比 / with installation and features）

裸序列，与 NCBI 注释逐基因比较，table2asn 检查能否直接提交。MitoM 的内存来自其内部的 MITOS2 tRNA 搜索。
Bare sequences, compared with the NCBI annotation; table2asn shows whether a record can be submitted. MitoM's memory comes from its internal MITOS2 tRNA search.

四个工具都完成的 17 条（最后三行：全部记录）/ the 17 records all four tools finished (last three rows: all records):

| 指标 / metric | MitoM | MitoZ | MITOS2 | MitoFinder |
|---|---|---|---|---|
| 蛋白基因两端完全一致 / PCG both ends exact | 82% | 78% | 51% | 59% |
| rRNA 两端一致 / rRNA both ends exact | 38% | 25% | 16% | 9% |
| tRNA 找回 / tRNA found | 98% | 93% | 96% | 97% |
| 无明显错误的记录 / records without a clear error | 10/17 | 1/17 | 0/17 | 0/17 |
| 0 个 NCBI ERROR 的记录 / records with 0 NCBI ERROR | 17/24 | 4/18 | 0 | 0 |
| 每条用时中位数 / median time | 175 s | 249 s | 385 s | 204 s |
| 峰值内存 / peak RAM | 0.7 GB | 0.5 GB | 1.0 GB | 0.24 GB |

原生核心（2026-10-22 起为默认，不再调用 MitoZ），同一批记录 / native core (the default since 2026-10-22, MitoZ not called), the same records. 表 / report：`benchmark/results/metazoa_mini_2026-10/report_native.md`。

| 指标 / metric | 原生核心 / native | MitoZ 路线 / MitoZ route | MitoZ |
|---|---|---|---|
| 蛋白基因两端一致 / PCG both ends exact | 247/306 | 240 | 174/230 |
| 17 条共同记录 / 17 shared records | 82.6% | 82.1% | 78.4% |
| rRNA 3 nt 内 / rRNA within 3 nt | 50% | 48% | 27% |
| tRNA 3 nt 内 / tRNA within 3 nt | 81% | 81% | 79% |
| 每条用时中位 / median time | 109 s | 175 s | 249 s |
| 峰值内存 / peak RAM | 0.32 GB | 0.73 GB | 0.49 GB |
| NCBI ERROR 记录 / records with an NCBI ERROR | 1 | 2 | 13 |

30 个门的 RefSeq（`benchmark/results/metazoa_173_2026-10/`）/ RefSeq across 30 phyla:

| 指标 / metric | 原生核心 / native | MitoZ 路线 / MitoZ route |
|---|---|---|
| 蛋白基因两端一致 / PCG both ends exact | 79.0% | 66.7% |
| 高于另一列的门 / phyla higher than the other column | 25 / 30 | |
| NCBI ERROR 记录 / records with an NCBI ERROR | 11 | 113 |
| 每条用时中位 / median time | 146 s | 325 s |

### 校正前后（2026-10-28）/ Raw and curated

数据 / data：上面记录中 MitoZ 支持的 18 条 / records above that MitoZ supports；`benchmark/results/curation_cmp_2026-10-28/`

箭头左边是校正前，右边是校正后；校正不改变 tRNA。最后一行是 MitoM 关掉校正时的结果（核心把终止子等交给校正规则）。
The arrow is before → after curation; curation does not change tRNAs. The last row is MitoM with curation off (the core leaves stops to the curation rules).

| 指标 / metric | MitoM | MitoZ 校正前 → 后 / MitoZ before → after | MITOS2 校正前 → 后 / MITOS2 before → after |
|---|---|---|---|
| 蛋白基因两端一致 / PCG both ends exact | 83.0% | 75.7% → 78.7% | 50.0% → 56.5% |
| rRNA 3 nt 内 / rRNA within 3 nt | 61.8% | 26.5% → 61.8% | 23.5% → 64.7% |
| tRNA 找回 / tRNA found | 99.5% | 92.9% | 95.9% |
| NCBI ERROR 记录 / records with an NCBI ERROR | 0 | 13/18 → 0 | 18/18 → 0 |
| MitoM 不校正 / MitoM, no curation | 58.7% | | |

### 环节动物（2026-10）/ Annelida

数据 / data：40 条 RefSeq 环节动物 / 40 RefSeq annelids；`benchmark/results/annelida_2026-10/`，`docs/benchmark_annelida_2026-10.md`

裸序列，与 NCBI 注释比较。明显错误 = 缺基因、内部终止、重复或边界差超过 30 nt。
Bare sequences, compared with the NCBI annotation. A clear error = missing gene, internal stop, duplicate or an end more than 30 nt off.

| 指标 / metric | MitoM | MitoZ | MITOS2 | MitoFinder |
|---|---|---|---|---|
| 蛋白基因两端完全一致 / PCG both ends exact | 80% | 63% | 54% | 48% |
| 有明显错误的记录 / records with a clear error | 16/40 | 36/40 | 39/40 | 34/39 |
| 每条用时中位数 / median time | 100 s | 244 s | 353 s | 206 s |

## 依赖 / Requirements

一条命令安装全部依赖：`bash install.sh`（见"快速上手"）；`mitom doctor` 检查。全部在一个 conda 环境 `mitotools`（`env/mitotools.yml`）：

- Python 3.11 及 Biopython、numpy、scikit-learn 1.7.2（MitoM-ML tRNA 判别器按此版本训练）、matplotlib（圈图）
- 注释核心与后续步骤：exonerate、BLAST、Infernal、HMMER、table2asn、ARAGORN、tRNAscan-SE、MAFFT + MITOS2 参考数据 refseq89m（`mitom setup`）
- Java（MiTFi 用，任何 ≥ 8 的版本）：用电脑上已有的；没有时 install.sh 往环境里加 openjdk 8
- 圈图测序深度：bwa；装有 minimap2 时自动改用它，快约 8 倍：`mamba install -n mitotools -c bioconda minimap2`（Linux 上环境多 0.5 GB）
- 可选，自己安装：MITOS2 2.1.10（`--engine mitos2`、`--rescue yes`；`bash install.sh --envs "mitotools mitos"`）、MitoZ 3.6（`--engine mitoz`；`mamba create -n mitoz -c bioconda -c conda-forge mitoz=3.6`）、Snakemake（多任务调度）

Everything is in one conda environment, `mitotools` (`env/mitotools.yml`): Python 3.11 with Biopython, numpy, scikit-learn 1.7.2 and matplotlib; the programs of the core and the later steps; the MITOS2 reference data via `mitom setup`. Java for MiTFi (any version >= 8): the computer's own; install.sh adds openjdk 8 to the environment only when there is none. Optional, installed by you: MITOS2 (`--engine mitos2`, `--rescue yes`), MitoZ (`--engine mitoz`), Snakemake, minimap2 (faster depth ring).

检查帮助：

```bash
mitom -h
mitom run -h
mitom import -h
mitom annotate -h
mitom reorient -h
mitom transfer_metadata -h
mitom replace_organism -h
```

