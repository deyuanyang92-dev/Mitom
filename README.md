# MitoM 使用说明 / User Guide

本文档说明 MitoM（命令 `mitom`，主程序 `batch_mitoz.py`）的安装、输入格式、各子命令、常用参数、输出文件和工作流程。
This document explains installation, inputs, subcommands, key options, outputs, and workflow of MitoM (command `mitom`, main program `batch_mitoz.py`).

## 快速上手 / Quick start

**1. 安装（一条命令）/ Install (one command)** — Linux 或 Windows WSL2；macOS 见 `install.sh` 开头说明，或用 Docker：

```bash
git clone https://github.com/deyuanyang92-dev/Mitom.git MitoM && cd MitoM
bash install.sh            # 一个 conda 环境 mitotools（BLAST、exonerate、Infernal、HMMER、Java、table2asn）+ 参考数据 + 命令 mitom
mitom doctor               # 检查：所有必需程序都应显示 OK
```

Docker：`docker build -t mitom . && docker run --rm -v "$PWD":/data mitom mitom run -i /data/genomes.gb -o /data/out`

**2. 运行 / Run**

```bash
# 从 NCBI 下载某类群的全部完整线粒体基因组（RefSeq + GenBank；--refseq_only 只要 NC_；--query 自定义检索式）
mitom fetch --taxon Annelida -o fetch_out            # -> fetch_out/Annelida.gb + fetch_report.tsv（可重复运行，只补新记录）

# GenBank 输入（推荐）：去重（同一基因组保留 NC_）→ 全部用 MitoM 注释 → cox1 定向 → 定向后有问题的才第二轮注释 → 校正 → NCBI 校验
mitom run -i fetch_out/Annelida.gb -o out            # 不论 NCBI 原注释好坏，每条都按同一流程重新注释

# FASTA 输入（没有分类信息）：不写 --clade 时沿用 Annelida 蛋白库并警告；其他类群请写明类群库和遗传密码
mitom run -i contigs.fasta -o out --clade Mollusca --genetic_code 5
```

`--clade` 可选：`Chordata Arthropoda Echinodermata Annelida-segmented-worms Bryozoa Mollusca Nematoda
Nemertea-ribbon-worms Porifera-sponges Cnidaria Platyhelminthes-flatworms Chaetognatha-arrow-worms Animal`
（`Animal` = 全动物蛋白库，用于其他门）。常用：`--threads 8 --max_tasks 4`（每个任务线程数 / 同时任务数）、
`--sbt template.sbt`（NCBI 提交者模板）。`--threads` 默认 `min(8, CPU 数)`，`--max_tasks` 默认 `CPU 数 // threads`；`--clade auto` 遇到分类信息里没有已知门名（只有纲级名如 `Vertebrata`、`Gastropoda` 会按所属门处理）、含多个门、或各条记录 CDS 的 `transl_table` 不一致时报错（给 `--clade` / `--genetic_code` 即可）。

**3. 运行过程与结果 / What you see and what you get**

运行时终端按步骤显示（每步：做什么、多少条、结果、用时、结果文件；注释阶段每 30 秒报告进度），完整细节在 `run.log`：

```text
[Step 1/10] Read input and remove duplicate records
    GenBank records read: 802
    removed 166: GenBank record that a RefSeq NC_ record was made from (NC_ kept)
    -> 636 records to annotate; reasons per record: reports/dedup_report.tsv  (5 s)
[Step 2/10] Round-1 annotation (core: native; protein genes, rRNAs, tRNAs)
    round-1 annotation: 120/636 (18%), 6 min 02 s so far, about 26 min left
...
[Step 4/10] Check the annotation after re-orientation; annotate again where it is offset (round 2)
    checked 636: 41 need round 2 (...); 595 keep the shifted round-1 annotation
...
[Step 9/10] Compare with the original annotation
```

输出目录（`out/README.txt` 逐项说明）：

| 文件 / 目录 | 内容 |
|---|---|
| `final_annotation/<登录号>.gb` | MitoM 最终注释，每条一个文件（已定向、校正、NCBI 命名） |
| `final_annotation_all.gb` | 同上，合并成一个文件 |
| `original_annotation/<登录号>.gb` | 输入记录原样（含 NCBI 原注释），文件名与上面对应，便于前后对比 |
| `original_vs_mitom.tsv` | 原注释 vs MitoM，每个基因一行：same / ends_differ / only_original / only_mitom … |
| `summary.tsv` | 每条记录一行：基因数、cox1 起始、用了第几轮注释、是否迁回原 header、NCBI ERROR 数 |
| `ncbi_submission/` | 每条的 NCBI 提交文件 `.fsa/.tbl/.sqn` 与 table2asn 结果；`submission_summary.tsv` ERROR 为 0 的可直接提交 |
| `reports/` | 各步骤报告：`dedup_report.tsv`（去重依据）、`reannotation_decision.tsv`（第二轮与否及原因）、`curation_report.tsv`（每处校正）等 |
| `run.log` | 全部过程与调用的命令 |

**4. 单独步骤 / Single steps**：`mitom -h` 列出全部命令（`fetch`、`dedup`、`check`、`curate`、`intron`、`trna`、`rescue`、`rrna-hmm`、`validate` 等），
`mitom <命令> -h` 查看参数。

**注释审核报告 / Annotation audit report**：`mitom check -i genomes.gb -o check_out` 审核已有注释（只出报告，不改文件；
可用于查看 NCBI 下载记录的原注释问题，或检查最终结果），`check_report.tsv` 每条一行（ok / flagged），`check_details.tsv` 列出每个问题。ERROR（触发重注）：无注释、缺蛋白基因或
rRNA、基因重复、内部终止子、长度非 3 倍数且无 transl_except、基因长度异常、基因重叠 > 100 nt、tRNA 少于 22（`--min_trna`）；
WARNING：起始/终止密码子、cox1 方向、tRNA 种类、非标准 CDS。`mitom run` 不使用原注释：去重后每条记录都按同一流程用 MitoM
重新注释；cox1 定向后由 `--if_reannotation_2`（默认 `auto`）检查，注释因定向出现偏移或问题的记录才做第二轮注释。
去重报告 `dedup_report.tsv` 每条一行，
`reason` 列给出丢弃依据（同 accession 取新版本 / RefSeq 记录 COMMENT 注明由该 GenBank 记录生成 / 序列完全相同），`related` 列给出保留的记录。

**默认流程 / Default workflow**：（`fetch` 下载 →）GenBank 去重 → 原生核心注释（`mito_core.py`：蛋白基因、rRNA、tRNA；`--core mitoz|hybrid` 时改用 MitoZ / MITOS2）→ cox1 定向 → 定向后检查，有问题的第二轮注释 →
元数据迁回 →（装有 MITOS2 时）补回漏注基因 → 仍缺的 12S/16S 用按门 HMM 补搜 → 内含子 → 端点校正与 NCBI 命名 → NCBI 校验与提交文件。

## MitoM：声明、致谢与引用 / Statement, acknowledgements and citation

**MitoM**（命令名 `mitom`；M = Merge / Metazoa）是一条动物线粒体基因组注释流程。流程本身（下表"MitoM 自有步骤"）是作者多年研究中逐步优化形成的；基因识别由"注释核心"完成：默认是 MitoM 自己的 `mito_core.py`（`--core native`），也可改用 MitoZ、MITOS2（`--core mitoz|hybrid`），每个环节用哪个核心，以 benchmark 结果为依据（`benchmark/`）。

| 环节 | 来源 |
|---|---|
| NCBI 下载（`mito_fetch.py`，按类群 / 检索式）、GenBank 输入去重（保留 RefSeq）、注释审核报告（`mito_check.py`）、内部短 ID、批量与并行 | MitoM 自有步骤 |
| 方向调整（cox1 起始、统一正链，cox1 缺失时依次用 cox2/cox3/cob/nad1） | MitoM 自有步骤 |
| 方向调整后第二轮注释 | MitoM 自有步骤 |
| 元数据迁回（原 GenBank header）、`/organism` 修复 | MitoM 自有步骤 |
| 内含子（`mito_intron.py`，按 MFannot 方法实现） | MitoM 自有步骤（方法来自 MFannot） |
| 跨序列起点 tRNA 修复、漏注基因从 MITOS2 救回（`mito_rescue.py`，可选：装有 MITOS2 时） | MitoM 自有步骤 |
| 仍缺的 12S/16S 用按门建立的 rRNA HMM 补搜（`mito_rrna_hmm.py`，库 `models/rrna_hmm/`，nhmmer；只在常规方法找不到时启用） | MitoM 自有步骤（思路参考 EMMA 的 nhmmer rRNA 模型） |
| 蛋白基因与 rRNA 端点校正（`mito_curate.py`）、NCBI table2asn 校验 | MitoM 自有步骤（规则依据 Fourdrilis 2018、Ghiselli 2021） |
| 多方法 tRNA 比对（`mito_trna.py`） | MitoM 自有步骤 |
| 蛋白基因、rRNA、tRNA 识别 | 注释核心：默认 `mito_core.py`（`--core native`：blastx + exonerate 蛋白基因，按门 rRNA profile HMM，cmsearch + MiTFi tRNA） |
| 蛋白基因、rRNA 识别（可选 `--core mitoz\|hybrid`） | 注释核心：MitoZ（含 MitoZ 自带但未开放的刺胞、扁形、毛颚、全动物蛋白库） |
| tRNA 识别（可选 `--core hybrid`） | 注释核心：MITOS2（灵敏补搜，并行运行） |

`--core native`（默认）= MitoM 自己的核心（`mito_core.py`：blastx 定位 + exonerate 定蛋白基因边界，参考蛋白来自 MITOS2 的 refseq89m；rRNA 用按门的 profile HMM；tRNA 用 MiTFi 协方差模型 + MiTFi），只需一个 conda 环境，每条记录约 1–2 分钟；`--core hybrid` = MitoZ 蛋白基因/rRNA + MITOS2 tRNA + MITOS2 漏注救回（需另装 mitoz、mitos 环境）；`--core mitoz` = 只用 MitoZ 作为核心。三者都运行全部 MitoM 自有步骤。

**Benchmark（多类群，2026-10）/ Metazoa benchmark**：8 个门 × 3 条 RefSeq 线粒体基因组（环节、节肢、脊索、刺胞、棘皮、软体、线虫、扁形），裸序列输入，与 NCBI 注释逐基因比较，并用 table2asn 检查能否直接提交（`docs/benchmark_metazoa_2026-10.md`，数据 `benchmark/results/metazoa_mini_2026-10/`）。四个工具都完成的 17 条：蛋白基因两端完全一致 MitoM 82%、MitoZ 78%、MITOS2 51%、MitoFinder 59%；rRNA 两端一致 38% / 25% / 16% / 9%；tRNA 找回 98% / 93% / 96% / 97%；记录无明显错误 10/17、1/17、0/17、0/17；0 个 NCBI ERROR 的记录 MitoM 17/24，MitoZ 4/18，MITOS2 与 MitoFinder 0。每条用时中位数 175 s / 249 s / 385 s / 204 s；峰值内存 0.7 / 0.5 / 1.0 / 0.24 GB（MitoM 的内存来自其内部的 MITOS2 tRNA 搜索）。MitoZ 不支持刺胞和扁形动物（无 clade）；MitoFinder 需要用户提供近缘参考基因组，2 条因无近缘参考失败。安装、功能、便捷性的对比见同一文档。 **原生核心（2026-10-22，默认 `--core native`，不再调用 MitoZ）在同一 24 条上**：蛋白基因两端一致 247/306（MitoZ 路线 240，MitoZ 174/230；17 条共同记录 82.6% / 82.1% / 78.4%），rRNA 3 nt 内 50%（48% / 27%），tRNA 3 nt 内 81%（81% / 79%），每条用时中位 109 s（175 / 249），峰值内存 0.32 GB（0.73 / 0.49），NCBI ERROR 记录 1（2 / 13）；只需一个 conda 环境。表见 `benchmark/results/metazoa_mini_2026-10/report_native.md`。 在 30 个门 173 条 RefSeq 记录上：蛋白基因两端一致 79.0%（MitoZ 路线 66.7%），30 个门中 25 个提升，NCBI ERROR 记录 11 vs 113，每条用时中位 146 s vs 325 s（`benchmark/results/metazoa_173_2026-10/`）。
On 17 records all four tools finished (8 phyla), protein genes with both ends exactly as RefSeq: MitoM 82 %, MitoZ 78 %, MITOS2 51 %, MitoFinder 59 %; records with zero NCBI table2asn errors: MitoM 17 of 24, MitoZ 4 of 18, MITOS2 and MitoFinder none. Details, installation and feature comparison: `docs/benchmark_metazoa_2026-10.md`.

**Benchmark（环节动物，2026-10）/ Annelida benchmark**：40 条 RefSeq 环节动物线粒体基因组，裸序列输入，与 NCBI 注释比较（`docs/benchmark_annelida_2026-10.md`，数据 `benchmark/results/annelida_2026-10/`）。蛋白基因两端完全一致：MitoM 80%、MitoZ 63%、MITOS2 54%、MitoFinder 48%；有明显错误（缺基因、内部终止、重复、边界差 >30 nt）的记录：MitoM 16/40、MitoZ 36/40、MITOS2 39/40、MitoFinder 34/39；每条用时中位数 MitoM 100 s、MitoZ 244 s、MITOS2 353 s、MitoFinder 206 s。

**MitoM** (command name `mitom`; M = Merge / Metazoa) is an annotation workflow for animal mitochondrial genomes. The workflow itself (the "MitoM steps" above: input deduplication, reorientation, round-2 annotation, metadata transfer, introns, origin-spanning tRNA fix, rescue of missed genes, per-phylum rRNA HMM fallback, boundary curation, NCBI validation) was developed by the author over years of research. Gene finding is done by an "annotation core", chosen by benchmark results (`benchmark/`): `--core native` (default) = MitoM's own `mito_core.py` (blastx + exonerate protein genes, per-phylum rRNA profile HMMs, cmsearch + MiTFi tRNAs; one conda environment); `--core hybrid` = MitoZ protein genes / rRNAs + MITOS2 tRNAs + MITOS2 rescue; `--core mitoz` = MitoZ alone. All three run every MitoM step.

**声明 / Statement**

- 本软件基于作者多年研究中的一些心得，整合了大家的优秀软件，并使用 Claude Code（Anthropic）开发。
- 作者不保留任何版权：本仓库代码以 CC0 1.0 放入公有领域（见 `LICENSE`），任何人都可以使用、修改、再发布，欢迎持续优化升级。
- 使用本软件时，请引用本 GitHub 仓库，其他没有任何要求。
- This software builds on the author's experience from years of research. It integrates excellent tools written by others and was developed with Claude Code (Anthropic).
- The author claims no copyright: the code in this repository is dedicated to the public domain under CC0 1.0 (see `LICENSE`). Anyone may use, modify and redistribute it, and further improvements are welcome.
- If you use it, please cite this GitHub repository. Nothing else is required.

**第三方软件 / Third-party tools**

默认核心（`--core native`）调用的程序（BLAST+、exonerate、Infernal、HMMER、Java、table2asn 等）由 `install.sh` 装入 conda 环境，各自保留原有许可证。随仓库分发的第三方文件：MITOS2 的 MiTFi（`models/mitfi/`，MIT，见其中的 LICENSE）；`install.sh` 另外下载 MITOS2 的参考数据 refseq89m（蛋白集与 tRNA/rRNA 模型，MIT）。MitoZ（GPL-3.0）只在 `--core hybrid` / `--core mitoz` 时需要，并不随本仓库分发。

The programs the default core runs (BLAST+, exonerate, Infernal, HMMER, Java, table2asn ...) are installed into a conda environment by `install.sh` and keep their own licences. Third-party files shipped here: MiTFi from MITOS2 (`models/mitfi/`, MIT); `install.sh` downloads the MITOS2 reference data refseq89m (protein sets, tRNA / rRNA models, MIT). MitoZ (GPL-3.0) is needed only for `--core hybrid` / `--core mitoz` and is not distributed here.

**致谢 / Acknowledgements**

感谢 MitoZ、MITOS2/MITOS、MiTFi、MFannot、Infernal、BLAST+、GeneWise、Exonerate、NCBI table2asn、Biopython 等软件的作者。没有他们的工作，就没有本软件。

We thank the authors of MitoZ, MITOS2/MITOS, MiTFi, MFannot, Infernal, BLAST+, GeneWise, Exonerate, NCBI table2asn, Biopython and the other tools listed below. This software would not exist without their work.

**引用 / Citation**

引用本软件的同时，请一并引用实际用到的软件 / When citing this software, please also cite the tools it ran:

- 本软件 / This software: https://github.com/deyuanyang92-dev/Mitom (`CITATION.cff`)
- MitoZ: Meng G, et al. 2019. Nucleic Acids Res. doi:10.1093/nar/gkz173
- MITOS2: Donath A, et al. 2019. Nucleic Acids Res. doi:10.1093/nar/gkz833
- MITOS: Bernt M, et al. 2013. Mol Phylogenet Evol. doi:10.1016/j.ympev.2012.08.023
- MiTFi (tRNA models used by MITOS2 and MitoZ): Jühling F, et al. 2012. Nucleic Acids Res. doi:10.1093/nar/gkr1131
- Infernal: Nawrocki EP, Eddy SR. 2013. Bioinformatics. doi:10.1093/bioinformatics/btt509
- BLAST+: Camacho C, et al. 2009. BMC Bioinformatics. doi:10.1186/1471-2105-10-421
- GeneWise (used by MitoZ): Birney E, et al. 2004. Genome Res. doi:10.1101/gr.1865504
- 内含子 / introns（`mito_intron.py`）— MFannot: Lang BF, et al. 2023. Front Plant Sci. doi:10.3389/fpls.2023.1222186; Exonerate: Slater GSC, Birney E. 2005. BMC Bioinformatics. doi:10.1186/1471-2105-6-31
- 端点校正依据 / boundary curation rules（`mito_curate.py`）: Fourdrilis S, et al. 2018. Sci Rep 8:17964. doi:10.1038/s41598-018-36428-7; Ghiselli F, et al. 2021. Phil Trans R Soc B. doi:10.1098/rstb.2020.0159
- 可选 tRNA 工具 / optional tRNA tools（`mito_trna.py`）: ARWEN doi:10.1093/bioinformatics/btm573; ARAGORN doi:10.1093/nar/gkh152; tRNAscan-SE 2.0 doi:10.1093/nar/gkab688
- Biopython: Cock PJA, et al. 2009. Bioinformatics. doi:10.1093/bioinformatics/btp163
- NCBI table2asn: https://www.ncbi.nlm.nih.gov/genbank/table2asn/

## 脚本定位 / What The Script Does

`batch_mitoz.py` 是 MitoM 流程的主程序：批量线粒体注释（默认用原生核心 `mito_core.py`，`--core mitoz|hybrid` 时用 MitoZ / MITOS2）及作者设计的方向调整、元数据迁回、内含子等步骤。

`batch_mitoz.py` is the main program of the MitoM workflow: batch mitochondrial annotation (the native core `mito_core.py` by default; MitoZ / MITOS2 with `--core mitoz|hybrid`) plus the author's own steps - reorientation, metadata transfer, introns and more.

它主要完成：

- 识别 FASTA 或 GenBank 输入。
- 把 GenBank 输入去重并转换为内部 FASTA。
- 为每条记录生成短内部 ID，例如 `S00001`，避免注释核心或 GenBank LOCUS 名称过长。
- 批量调用注释核心（默认 `mito_core.py`；可选 MitoZ）。
- 按目标基因重定向序列，默认 `cox1`。
- 可选对重定向后的序列做第二轮注释。
- 把原始 GenBank header metadata 转移回最终结果。
- 修复核心输出中 `FEATURES source` 的 `/organism=`。
- 输出最终 `.gbf`、统计表、日志、失败列表和合并文件。

It can:

- Detect FASTA or GenBank input.
- Deduplicate GenBank records and convert them to internal FASTA.
- Generate short internal IDs such as `S00001`.
- Run the annotation core in batch (native `mito_core.py` by default; MitoZ optional).
- Reorient records by a target gene, defaulting to `cox1`.
- Optionally run a second annotation pass after reorientation.
- Restore original GenBank header metadata.
- Restore `/organism=` inside `FEATURES source`.
- Write final `.gbf` files, summary tables, logs, failure lists, and merged GenBank archives.

## 依赖 / Requirements

一条命令安装全部依赖：`bash install.sh`（见"快速上手"）；`mitom doctor` 检查。下面是各组件：

- 默认核心（`--core native`）：exonerate、BLAST、Infernal、HMMER、Java、table2asn、MAFFT（env/mitotools.yml）+ MITOS2 参考数据 refseq89m
- 可选：MitoZ（`--core mitoz|hybrid`）、MITOS2 2.1.10（`--core hybrid` 的 tRNA、漏注救回）

基础依赖：

- Python 3.8+
- Biopython
- MitoZ，可选，只在 `--core mitoz|hybrid` 时需要
- Snakemake，可选，用于调度多任务
- Conda，可选，用于 `--conda_env_mitoz`

Basic requirements:

- Python 3.8+
- Biopython
- MitoZ, optional: only for `--core mitoz|hybrid`
- Snakemake, optional scheduler
- Conda, optional when using `--conda_env_mitoz`

检查帮助：

```bash
python3 batch_mitoz.py --help
python3 batch_mitoz.py run --help
python3 batch_mitoz.py annotate --help
python3 batch_mitoz.py reorient --help
python3 batch_mitoz.py transfer_metadata --help
python3 batch_mitoz.py replace_organism --help
```

## 输入格式 / Input Formats

### FASTA 输入 / FASTA Input

支持扩展名：

```text
.fa .fasta .fna .fas
```

FASTA header 中可以包含 topology：

```text
>B9  topology=linear
ATGC...
```

如果没有写 `topology=linear` 或 `topology=circular`，脚本默认按 `circular` 处理。

Supported FASTA extensions are `.fa`, `.fasta`, `.fna`, and `.fas`. If no topology tag is present, the script assumes `circular`.

### GenBank 输入 / GenBank Input

支持扩展名：

```text
.gb .gbk .genbank .gbf .gbff
```

GenBank 输入会先执行：

1. 解析所有 record。
2. 跳过无序列 record，但保留其 metadata 来源索引。
3. 去除完全重复、版本重复或序列重复的 record。
4. 转换为内部 FASTA 后再进入 MitoZ。

GenBank input is parsed, deduplicated, converted to internal FASTA, and then passed to MitoZ.

注意：当前版本已避免把本地样本名 `B64_2` 误认为 accession version，因此 `B64` 和 `B64_2` 会保持为不同样本。

Note: local IDs such as `B64_2` are kept distinct from `B64`.

### 目录扫描与多序列 FASTA / Directory Scanning And Multi-record FASTA

输入目录会递归扫描。若同一目录同时包含 FASTA 和 GenBank，`annotate` / `run` 会优先使用 GenBank 文件，并忽略该目录中的 FASTA 文件。

Input directories are scanned recursively. If a directory contains both FASTA and GenBank files, `annotate` / `run` use the GenBank files and ignore the FASTA files in that mixed directory.

多记录 FASTA 默认会按 record 拆分：

```bash
--if_split_fasta yes
```

如果设置：

```bash
--if_split_fasta no
```

但输入 FASTA 含多个 record，脚本会报错退出。

By default, multi-record FASTA files are split into individual records. With `--if_split_fasta no`, a multi-record FASTA causes the script to exit with an error.

### FASTQ 输入 / Optional FASTQ Input

FASTQ 仅用于 `annotate` / `run`，并且只在同时找到 `fq1` 和 `fq2` 时传给 MitoZ。脚本按最终 `sample_name` 查找双端 FASTQ。

FASTQ options are only used by `annotate` / `run`. Files are passed to MitoZ only when both `fq1` and `fq2` are found. Matching is based on the final `sample_name`.

如果 MitoZ 需要 paired FASTQ，可使用：

```bash
--infq /path/to/fastq_dir
--fq_position flat
--suffix_fq _1.clean.fq.gz,_2.clean.fq.gz
```

`--fq_position flat` 表示 FASTQ 文件直接在 `--infq` 目录下：

```text
<sample_name>_1.clean.fq.gz
<sample_name>_2.clean.fq.gz
```

`--fq_position sample_dir` 表示每个样本有单独子目录：

```text
fastq_dir/<sample_name>/<sample_name>_1.clean.fq.gz
fastq_dir/<sample_name>/<sample_name>_2.clean.fq.gz
```

在 `sample_dir` 模式下，脚本会优先找精确文件名；如果没有精确文件名，会取该样本子目录中第一个匹配 `*suffix1` 和 `*suffix2` 的文件。`--suffix_fq` 必须是两个逗号分隔的后缀。

In `sample_dir` mode, exact filenames are preferred. If exact names are absent, the first files matching `*suffix1` and `*suffix2` in the sample subdirectory are used. `--suffix_fq` must contain exactly two comma-separated suffixes.

## 子命令总览 / Subcommands

| 子命令 | 用途 | English |
|---|---|---|
| `run` | 完整流程：注释、重定向、可选二轮注释、最终 GenBank 生成 | Full pipeline |
| `annotate` | 只做一轮批量注释和最终整理 | Annotation only |
| `reorient` | 只对已有 GenBank 按目标基因重定向 | Reorient existing GenBank |
| `transfer_metadata` | 单独转移 GenBank header metadata | Standalone metadata transfer |
| `replace_organism` | 单独修复 `/organism=` 和 SOURCE | Standalone organism/source restoration |

## 推荐主流程：run / Recommended Main Workflow

`run` 是最完整、最推荐的子命令。

`run` is the recommended full pipeline.

典型流程：

```text
输入 FASTA / GenBank
→ 内部 FASTA 和 id_map.tsv
→ 第一轮注释（注释核心；默认 mito_core.py）
→ 按 cox1 或指定基因重定向
→ 可选第二轮注释
→ 最终 ID 重写
→ metadata 恢复
→ /organism= 同步
→ 统计和合并输出
```

### 用 MitoZ 核心运行（可选）/ Run With The MitoZ Core (optional, --core hybrid)

```bash
python3 batch_mitoz.py run \
  -i /path/to/input_fasta_or_genbank \
  -o mitoz_run_out \
  --core hybrid \
  --conda_env_mitoz mitoz3.6 \
  --scheduler_backend threads \
  --threads 8 \
  --max_tasks 4 \
  --clade Arthropoda \
  --genetic_code 5 \
  --which_gene_first cox1 \
  --if_reannotation_2 yes \
  --preserve_phase1_metadata yes
```

### 用指定 MitoZ 路径运行 / Run With A MitoZ Executable Path

```bash
python3 batch_mitoz.py run \
  -i /path/to/input_fasta_or_genbank \
  -o mitoz_run_out \
  --mitoz_path /path/to/mitoz \
  --scheduler_backend threads \
  --threads 8 \
  --max_tasks 4 \
  --clade Arthropoda \
  --genetic_code 5
```

`--mitoz_path` 可以是 MitoZ 可执行文件，也可以是包含 `mitoz` 的目录。

`--mitoz_path` can point to the MitoZ executable or to a directory containing it.

### 当前目录 FASTA 示例 / Example For The Current Directory

当前目录已有：

```text
summary_all.fasta
rename_summary_all.fasta
```

如果要用已改名后的 FASTA 跑完整流程，可用：

```bash
python3 batch_mitoz.py run \
  -i rename_summary_all.fasta \
  -o mitoz_run_renamed_out \
  --conda_env_mitoz mitoz3.6 \
  --scheduler_backend threads \
  --threads 8 \
  --max_tasks 4 \
  --clade Arthropoda \
  --genetic_code 5 \
  --which_gene_first cox1
```

请把 `mitoz3.6`、`Arthropoda`、`genetic_code` 改成你真实数据需要的值。

Replace `mitoz3.6`, `Arthropoda`, and `genetic_code` with values appropriate for your dataset.

### 带 FASTQ 的 run 示例 / run With Optional FASTQ

```bash
python3 batch_mitoz.py run \
  -i input.fasta \
  -o mitoz_run_out \
  --conda_env_mitoz mitoz3.6 \
  --infq /path/to/fastq_dir \
  --fq_position sample_dir \
  --suffix_fq _1.clean.fq.gz,_2.clean.fq.gz \
  --scheduler_backend threads \
  --threads 8 \
  --max_tasks 4 \
  --clade Arthropoda \
  --genetic_code 5
```

### Snakemake 调度示例 / Snakemake Scheduler Example

```bash
python3 batch_mitoz.py run \
  -i input.fasta \
  -o mitoz_run_out \
  --conda_env_mitoz mitoz3.6 \
  --scheduler_backend snakemake \
  --snakemake_bin /path/to/snakemake \
  --smk_cores 64 \
  --smk_mem_mb 256000 \
  --smk_latency_wait 60 \
  --smk_extra '--profile slurm' \
  --threads 8 \
  --max_tasks 8 \
  --clade Arthropoda \
  --genetic_code 5
```

## 一轮注释：annotate / Annotation Only

`annotate` 只运行第一轮注释（`--core` 选定的核心，默认原生核心），并把结果整理为最终 `.gbf`。

`annotate` runs one annotation pass (the core chosen by `--core`, native by default) and finalizes the `.gbf` output.

注意：`annotate` 不做第二轮注释，但最终 `.gbf` 写出前仍会执行 ID 重写，并按 `--which_gene_first` 做一次保证性重定向；默认目标基因为 `cox1`，负链目标基因默认会反向互补。`--preserve_phase1_metadata` 在 `annotate` 中默认是 `no`，只有 GenBank 输入且显式设为 `yes` 时才恢复原始 header metadata。

Note: `annotate` does not run phase-2 annotation, but final `.gbf` files are still ID-rewritten and reoriented by `--which_gene_first`. Metadata preservation defaults to `no` for `annotate`.

```bash
python3 batch_mitoz.py annotate \
  -i /path/to/fasta_or_genbank \
  -o annotate_out \
  --conda_env_mitoz mitoz3.6 \
  --scheduler_backend threads \
  --threads 8 \
  --max_tasks 4 \
  --clade Arthropoda \
  --genetic_code 5
```

适合：

- 不需要第二轮注释。
- 只想批量运行 MitoZ 并标准化输出文件名。

Use it when one annotation pass is enough.

追加 MitoZ 参数示例：

```bash
python3 batch_mitoz.py annotate \
  -i input.fasta \
  -o annotate_out \
  --conda_env_mitoz mitoz3.6 \
  --mitoz_extra '--your-extra-mitoz-flag value'
```

## 只重定向：reorient / Reorient Only

`reorient` 不运行 MitoZ，只处理已有 GenBank 文件。

`reorient` does not run MitoZ. It only reorients existing GenBank files.

`reorient` 的 `-i` 必须是 GenBank 文件或 GenBank 目录；FASTA 输入不会被该 workflow 处理。由于该子命令复用了通用参数，帮助中会显示 `--if_split_fasta` 和 `--id_maxlen`，但它们对 `reorient` 无实际作用。

For `reorient`, `-i` must be a GenBank file or directory. FASTA input is not processed by this workflow. `--if_split_fasta` and `--id_maxlen` appear in the shared help output but have no practical effect for `reorient`.

```bash
python3 batch_mitoz.py reorient \
  -i /path/to/genbank_dir \
  -o reorient_out \
  --which_gene_first cox1
```

如果目标基因在负链上，默认会反向互补，使目标基因转到正向。如果不想这样做：

```bash
python3 batch_mitoz.py reorient \
  -i /path/to/genbank_dir \
  -o reorient_out \
  --which_gene_first cox1 \
  --no_force_forward
```

By default, minus-strand target genes are reverse-complemented. Use `--no_force_forward` to disable that behavior.

## 单独转移 metadata / Standalone Metadata Transfer

`transfer_metadata` 用原始 GenBank 的 header 区域替换重新注释 GenBank 的 header 区域，保留新文件的 `FEATURES` 和 `ORIGIN` 主体。

`transfer_metadata` replaces the header section of re-annotated GenBank files with metadata from the original GenBank files while keeping the new `FEATURES` and `ORIGIN`.

```bash
python3 batch_mitoz.py transfer_metadata \
  -m /path/to/original_gb_dir \
  -n /path/to/new_body_gb_dir \
  -o transfer_metadata_out \
  --smart_metadata no
```

单文件模式可指定输出文件名：

```bash
python3 batch_mitoz.py transfer_metadata \
  -m original.gbf \
  -n reannotated.gbf \
  -o transfer_metadata_out \
  --output_name final_with_metadata.gbf
```

如果需要自动接受模糊匹配：

```bash
python3 batch_mitoz.py transfer_metadata \
  -m original.gbf \
  -n reannotated.gbf \
  -o transfer_metadata_out \
  --smart_metadata yes \
  --auto_accept_metadata yes
```

注意：目录模式下如果 metadata 目录中存在多个相同 stem 的文件，脚本会报错退出，避免错配 header。

In directory mode, duplicate metadata stems are rejected to avoid header mismatches.

输出与行为：

- 输出目录默认是 `transfer_metadata_out/`。
- 新文件默认写入 `transfer_metadata_out/transferred_gb/`。
- 审计表是 `transfer_metadata.tsv`。
- 日志是 `transfer_metadata.log`。
- 目录模式按文件 stem 配对 `metadata` 与 `new_body`；metadata 目录中重复 stem 会报错退出。
- `--output_name` 仅单文件模式有效。
- `--overwrite_new_body yes` 会直接覆盖 `-n/--new_body` 指定的 body 文件。
- `--if_replace_metadata no` 时不替换 header，只复制或保留 new body。

Outputs and behavior:

- Default output root: `transfer_metadata_out/`.
- New files are written to `transfer_metadata_out/transferred_gb/` by default.
- Audit table: `transfer_metadata.tsv`.
- Log file: `transfer_metadata.log`.
- Directory mode pairs metadata and body files by stem; duplicate metadata stems are rejected.
- `--output_name` is valid only in single-file mode.
- `--overwrite_new_body yes` overwrites the files supplied by `-n/--new_body`.
- `--if_replace_metadata no` keeps or copies the new body without replacing the header.

## 单独修复 organism / Standalone Organism Replacement

`replace_organism` 修复 `FEATURES source` 中的 `/organism=`，也可以同步 SOURCE 文本。

`replace_organism` restores `/organism=` inside `FEATURES source` and can also sync the SOURCE text.

```bash
python3 batch_mitoz.py replace_organism \
  -m /path/to/original_gb_dir \
  -n /path/to/new_body_gb_dir \
  -o replace_organism_out
```

如需覆盖输入文件：

```bash
python3 batch_mitoz.py replace_organism \
  -m /path/to/original_gb_dir \
  -n /path/to/new_body_gb_dir \
  --overwrite_new_body yes
```

Use overwrite mode carefully because it modifies the input body files in place.

输出与匹配规则：

- 输出目录默认是 `replace_organism_out/`。
- 新文件默认写入 `replace_organism_out/replaced_organism/`。
- 审计表是 `replace_organism.tsv`，列为 `source_file`, `output_file`, `replaced`, `not_found`, `total`, `note`。
- 日志是 `replace_organism.log`。
- 匹配依据是 GenBank accession 及其兼容变体。
- 匹配成功时更新 SOURCE 文本和 source feature 的 `/organism=`。
- 未匹配 record 保持不变。

Outputs and matching:

- Default output root: `replace_organism_out/`.
- New files are written to `replace_organism_out/replaced_organism/` by default.
- Audit table: `replace_organism.tsv`.
- Log file: `replace_organism.log`.
- Matching uses GenBank accessions and compatible accession variants.
- On match, SOURCE text and source-feature `/organism=` are updated.
- Unmatched records are left unchanged.

## 关键参数 / Key Options

| 参数 | 适用子命令 | 中文说明 | English |
|---|---|---|---|
| `-i`, `--input` | `run`, `annotate`, `reorient` | FASTA 或 GenBank 文件/目录；`reorient` 只处理 GenBank | Input FASTA or GenBank; `reorient` accepts GenBank only |
| `-o`, `--out_root` | all | 输出根目录；`run`/`annotate`/`reorient` 默认 `mitoz_batch_out`，其他子命令各有默认值 | Output root; defaults vary by subcommand |
| `--which_gene_first` | `run`, `annotate`, `reorient` | 重定向起始基因，默认 `cox1` | Target gene for reorientation, default `cox1` |
| `--fallback_genes` | `run`, `annotate` | 找不到 `--which_gene_first`（默认 cox1）时依次尝试的基因，默认 `cox2,cox3,cob,nad1`；`none` = 不旋转（旧行为）。用了备选基因的样本列在 `oriented_by_fallback.txt`。基因被注释成多段（如 intron 两侧外显子）时取转录方向最 5′ 的一段；线性序列从不旋转 | Fallback start genes when the first gene is absent; 5'-most piece of split genes; linear records never rotated |
| `--force_forward` / `--no_force_forward` | `run`, `annotate`, `reorient` | 目标基因在负链时是否反向互补，默认开启 | Reverse-complement minus-strand target genes by default |
| `--if_split_fasta` | `run`, `annotate` | 多序列 FASTA 是否拆分，默认 `yes`；对 `reorient` 无实际作用 | Split multi-record FASTA, default `yes`; no practical effect for `reorient` |
| `--id_maxlen` | `run`, `annotate` | 内部 ID 最大长度，必须至少 6；对 `reorient` 无实际作用 | Maximum internal ID length, minimum 6; no practical effect for `reorient` |
| `--conda_exe` | `run`, `annotate` | Conda 可执行文件名，默认 `conda` | Conda executable, default `conda` |
| `--conda_env_mitoz` | `run`, `annotate` | 含 MitoZ 的 Conda 环境名 | Conda environment containing MitoZ |
| `--mitoz_path` | `run`, `annotate` | MitoZ 可执行文件或目录 | Path to MitoZ executable or directory |
| `--threads` | `run`, `annotate` | 每个 MitoZ 任务线程数，必须为正整数；默认 `min(8, CPU 数)` | Threads per MitoZ job, positive integer; default `min(8, CPUs)` |
| `--max_tasks` | `run`, `annotate` | 同时运行的任务数，必须为正整数；默认 `max(1, CPU 数 // threads)` | Maximum concurrent jobs, positive integer; default `max(1, CPUs // threads)` |
| `--clade` | `run`, `annotate` | 类群蛋白库，默认 `auto`：从 GenBank 分类信息判断（一次一个门；纲级名按所属门处理；分类信息里没有已知门名、含多个门、或各条记录 CDS `transl_table` 不一致时报错；**完全没有分类信息**（FASTA、无 lineage 的 GenBank）时沿用 Annelida 并警告）；可选值见"快速上手" | Taxon protein database, default `auto` (from the GenBank taxonomy, one phylum per run; class names count for their phylum; a taxonomy without a known phylum, several phyla, or differing CDS transl_tables are errors; input with NO taxonomy at all - FASTA, GenBank without lineage - falls back to Annelida with a warning) |
| `--genetic_code` | `run`, `annotate` | 遗传密码表，无默认固定值：`--clade auto` 时取输入各记录 CDS 的 transl_table（记录间不一致则报错），无则该门常用码；显式 `--clade` 时取 MitoZ 对该类群的密码（Chordata 2、Echinodermata 9、Porifera 4，其余 5；Cnidaria、Platyhelminthes、Chaetognatha、Animal 必须显式给出）；显式给出时优先 | Genetic code (no fixed default): with `--clade auto` from the per-record CDS transl_table (differing records are an error), else the phylum's usual code; with an explicit `--clade` MitoZ's code of that clade (2 Chordata, 9 Echinodermata, 4 Porifera, else 5; required for Cnidaria, Platyhelminthes, Chaetognatha, Animal); an explicit value wins |
| `--core` | `run`, `annotate` | `native`（默认：MitoM 自己的核心 `mito_core.py`——blastx + exonerate 定蛋白基因，按门 rRNA profile HMM，Infernal + MiTFi tRNA 模型；只需 mitotools 环境 + refseq89m）/ `hybrid`（MitoZ 蛋白/rRNA + MITOS2 tRNA + MITOS2 补漏；需 mitoz、mitos 环境）/ `mitoz`（只用 MitoZ 核心，等于 `--trna mitoz --rescue no --rrna_hmm no`） | Annotation core: native (default; MitoM's own core, one environment), hybrid (MitoZ + MITOS2), mitoz (MitoZ alone) |
| `--trna` | `run`, `annotate` | `mitos2`（默认，灵敏补搜并行）/ `mitoz`（MiTFi） | tRNA step: `mitos2` (default) or `mitoz` |
| `--rescue` | `run`, `annotate` | 从 MITOS2 补回核心漏注的蛋白基因/rRNA：`auto`（默认：找到 MITOS2 且核心为 mitoz/hybrid 时运行；native 核心自己重试其他位点）/ `yes`（native 也运行）/ `no` | Rescue genes the core missed from MITOS2: `auto` (default; mitoz / hybrid core only) / `yes` / `no` |
| `--rrna_hmm`, `--rrna_hmm_lib`, `--phylum` | `run`, `annotate` | 仍缺的 12S/16S 用按门 HMM 补搜：`auto`（默认）/ `yes` / `no`；库默认 `models/rrna_hmm`；`--phylum` 默认取 `--clade auto` 判出的门、显式 `--clade` 对应的门（`--clade Animal` 时不确定，FASTA 输入请自己给 `--phylum`），都没有时才用记录自己的分类信息；没有可用模型的记录会在日志和屏幕上提示 | Per-phylum rRNA HMM fallback for rRNAs still missing: `auto` (default) / `yes` / `no`; `--phylum` defaults to the phylum found by `--clade auto` or that of an explicit `--clade` (unknown for `--clade Animal`: give it for FASTA input), then to the record's own taxonomy; records without a model are reported on screen and in the log |
| `--intron_fix` | `run`, `annotate` | 内含子（MFannot 方法）：`auto`（默认，找到 exonerate 即运行）/ `yes` / `no` | Introns (MFannot method): `auto` (default) / `yes` / `no` |
| `--curate`, `--sbt` | `run`, `annotate` | 端点校正 + NCBI 命名 + 提交文件（默认 `yes`）；`--sbt` 为 NCBI 提交者模板 | Curation + NCBI names + submission files (default `yes`); `--sbt` submitter template |
| `--species_name` | `run`, `annotate` | 传给 MitoZ annotate 的物种名 | Species name passed to MitoZ annotate |
| `--mitoz_extra ARG` | `run`, `annotate` | 追加传给 MitoZ annotate 的参数，可重复 | Extra MitoZ annotate arguments, repeatable |
| `--infq`, `--fq_position`, `--suffix_fq` | `run`, `annotate` | 可选双端 FASTQ；`--suffix_fq` 必须有两个逗号分隔后缀 | Optional paired FASTQ; suffix list must contain two comma-separated values |
| `--scheduler_backend` | `run`, `annotate` | `auto` / `threads` / `snakemake`，默认 `auto` | Scheduler backend, default `auto` |
| `--snakemake_bin`, `--smk_cores`, `--smk_mem_mb`, `--smk_latency_wait`, `--smk_extra` | `run`, `annotate` | Snakemake 后端参数；数值项必须为正整数 | Snakemake backend options; numeric values must be positive integers |
| `--if_reannotation_2` | `run` | 第二轮注释：`auto`（默认；第一轮用灵敏 tRNA 模式，平移后检查——原序列两端 50 nt 内有**看起来被截断**的基因（CDS 长度非 3 的倍数且无 T/TA 终止尾、靠序列端的一侧缺起始/终止密码子，按所用遗传密码判断；tRNA 不在 55–95 nt；s-rRNA < 600 nt、l-rRNA < 900 nt）、缺标准基因（COX1-3、CYTB、ND1-6、ND4L、两个 rRNA；ATP6/ATP8 在整个门中可能缺失，不算）、平移前后基因序列不一致、基因跨新起点——只重注未通过的；决定写 `reannotation_decision.tsv`；6 条 A/B 结果完全一致）/ `yes`（每条都重注）/ `no` | Phase-2 re-annotation: `auto` (default: only records failing the post-shift check: truncated-looking genes at the ends of the round-1 sequence, missing standard genes (not ATP6/ATP8), changed gene sequences, genes across the new origin; decisions in `reannotation_decision.tsv`) / `yes` (every record) / `no` |
| `--keep_intermediate` | `run`, `annotate` | 是否保留中间目录（`00.*`、`01.*`、`03.*`、`04.*`、`05.*`、`_snakemake_*`），默认 `no`：最终输出写完后删除以节约空间；`yes` 时 `run` 把它们放进 `work/`。问题样本的 MitoZ 日志总是保留（`run`：`reports/failed_logs/`，`annotate`：`failed_logs/`） | Keep the intermediate folders (run: moved to `work/`); default `no` deletes them after the final output; MitoZ logs of problem samples are always kept (reports/failed_logs/ or failed_logs/) |
| `--preserve_phase1_metadata` | `run`, `annotate` | `run` 默认 `yes`；`annotate` 默认 `no` | Defaults to `yes` for `run`, `no` for `annotate` |
| `--smart_metadata` | `run`, `transfer_metadata` | 开启模糊 metadata 匹配 | Enable fuzzy metadata matching |
| `--auto_accept_metadata` | `run`, `transfer_metadata` | 自动接受最佳模糊匹配；必须配合 `--smart_metadata yes` | Auto-accept best fuzzy match; requires `--smart_metadata yes` |
| `-m`, `--metadata` | `transfer_metadata`, `replace_organism` | 原始 GenBank metadata 来源 | Original GenBank metadata source |
| `-n`, `--new_body` | `transfer_metadata`, `replace_organism` | 新 GenBank body 来源 | New GenBank body source |
| `--if_replace_metadata` | `transfer_metadata` | 是否替换 header metadata，默认 `yes` | Replace header metadata, default `yes` |
| `--overwrite_new_body` | `transfer_metadata`, `replace_organism` | 是否原地覆盖 `new_body` | Overwrite `new_body` in place |
| `--output_name` | `transfer_metadata` | 单文件模式输出文件名 | Output filename in single-file mode |

## 调度建议 / Scheduler Recommendations

### threads 后端 / threads Backend

本地机器推荐先用：

```bash
--scheduler_backend threads --threads 8 --max_tasks 4
```

实际总线程消耗约为：

```text
threads × max_tasks
```

例如 `8 × 4 = 32`。请不要超过机器可承受的 CPU 和内存。

The approximate CPU demand is `threads × max_tasks`.

### Snakemake 后端 / Snakemake Backend

```bash
--scheduler_backend snakemake \
--smk_cores 64 \
--smk_mem_mb 256000 \
--threads 8 \
--max_tasks 8
```

当前版本生成的 Snakefile 会把每个任务声明为真实 `threads`，并声明每任务内存资源，避免调度器误判。

The generated Snakefile declares real per-job `threads` and memory resources.

## 输出目录 / Output Layout

`run` 的输出目录（终端第 10 步会逐项打印，目录里的 `README.txt` 有同样的说明）：

```text
out/
├── final_annotation/           # 最终 GenBank，每条记录一个文件 <登录号>.gb（方向调整、补漏、端点校正、NCBI 命名之后）
├── final_annotation_all.gb     # 上面所有记录合并成一个文件
├── original_annotation/        # 输入 GenBank 记录原样（去重后保留的每一条），用于前后对比；FASTA 输入没有此目录
├── original_vs_mitom.tsv       # 原注释 vs MitoM，每个基因一行：same / ends_differ / name_differs / strand_differs /
│                               #   only_original / only_mitom / same_place（坐标统一换算到原序列）
├── summary.tsv                 # 每条记录一行：蛋白/tRNA/rRNA 数、cox1 是否在首位、用的是第几轮注释、是否恢复原 GenBank 头、NCBI ERROR 数
├── ncbi_submission/            # NCBI 提交文件（需 table2asn）：<登录号>.fsa/.tbl（BankIt 上传）、.sqn/.val（table2asn）、
│                               #   submission_summary.tsv（每条记录 ERROR/WARNING 数）、validation_messages.tsv
├── reports/                    # 各步骤的报告：dedup_report.tsv、reannotation_decision.tsv、curation_report.tsv、validation_compare.tsv、
│                               #   rescue_report.tsv、rrna_hmm_report.tsv、intron_report.tsv、metadata_transfer_log.tsv、
│                               #   oriented_by_fallback.txt、all_problem_samples.txt、*_failed.txt、run_params.json、
│                               #   failed_logs/（问题样本的 MitoZ 日志）
├── run.log                     # 完整日志（终端只显示步骤和警告）
├── work/                       # 只在 --keep_intermediate yes 时存在：00.*、01.*、03.*、04.*、05.*、_snakemake_* 中间目录
└── README.txt                  # 本目录各项说明
```

主要文件：

- `final_annotation/<登录号>.gb`：**最终 GenBank**：端点校正（mito_curate.py）并统一为 NCBI 命名——CDS `COX1…CYTB`、`ND1…ND6`、`ND4L`；rRNA `rrnS`/`rrnL`（12S/16S ribosomal RNA）；tRNA `/gene` 只写符号（`trnL1`/`trnL2`、`trnS1`/`trnS2` 按反密码子，沿用 MITOS 约定 L1 = CUN、L2 = UUR、S1 = AGN、S2 = UCN；脊椎动物 HGNC/RefSeq 的 L1/L2 编号相反，请按投稿类群核对），反密码子写在 `/note`；非法起始密码子按 NCBI 惯例标 5′ partial（"start codon not determined"）。`--curate no` 时为校正前的注释。
- `ncbi_submission/`：每条记录的 NCBI 提交文件和校验报告。`submission_summary.tsv` 中 ERROR 为 0 的记录可直接提交；GenBank 平面文件（.gb）本身不能直接提交给 NCBI。
- `summary.tsv` 的 `annotation_used`：`round 2` = 定向后重注；`round 1, shifted` = 第一轮注释随序列平移（通过检查，不必重注）；`round 1, not re-oriented` = 未找到定向基因。
- 输出目录里用户自己的文件和目录不会被移动或删除；同一目录再次运行会替换上一次的结果。

The `run` output: final_annotation/ (one GenBank per record), final_annotation_all.gb, original_annotation/ (the input records, compared in original_vs_mitom.tsv), summary.tsv, ncbi_submission/, reports/, run.log and, with `--keep_intermediate yes`, work/. README.txt in the folder explains every item. Files the user put into the folder are left alone; a second run into the same folder replaces the previous results.

`annotate`（只注释：不去重、不定向、不比较）仍是平铺布局：

```text
annotate_out/
├── pipeline.log
├── run_params.json
├── 00.internal_fasta/ (id_map.tsv), 00.gb_converted_fasta/ (dedup_report.tsv), 01.mitoz_anno1/   # --keep_intermediate yes 时保留
├── 04.final_gb/                # 端点校正前
├── 05.curated/                 # --curate yes（默认）：端点校正 + NCBI 命名；curation_report.tsv
├── 06.ncbi_submission/         # NCBI 提交文件（需 table2asn）
├── rescue_report.tsv, rrna_hmm_report.tsv, intron_report.tsv, annotation-summary.tsv
├── all.final.gbf
└── annotation_failed.txt, final_write_failed.txt, failed_logs/
```

`annotation-summary.tsv`（`annotate`）中的 `cox1_orientation_ok` 是历史列名。当 `--which_gene_first` 不是 `cox1` 时，该列实际表示指定目标基因是否重定向成功。

The `cox1_orientation_ok` column name is historical. If `--which_gene_first` is not `cox1`, the column still means whether the selected target gene was successfully used for reorientation.

列表和 TSV 在内容为空时可能不会生成。

Some report files are written only when non-empty.

## 注释后校正 / Post-annotation curation (mito_curate.py)

对注释好的 GenBank（MitoZ、MITOS 等输出）单独运行，按 Fourdrilis et al. 2018 的软体动物（腹足类）注释准则修正 PCG 与 12S/16S，准则引自 Ghiselli et al. 2021, *Phil. Trans. R. Soc. B* 376:20200159（doi:10.1098/rstb.2020.0159，第 3 节）：

| 规则 | 准则 | 做法 |
|---|---|---|
| P1 | (i) 起始密码子取离上游基因最近且不重叠者 | 起点与上游 tRNA/rRNA 或同框 PCG 重叠 → 后移到第一个起始密码子，但只在 15 nt 内有起始密码子时才移动，否则保留重叠（173 条多类群真值：≤15 nt 的移动对 31 错 12，更远的移动对 14 错 28）；首密码子非起始密码子时才向上游找；`--aggressive_start_genes` 列出的基因总是按准则找最上游——默认为空，`mitom run` 只在环节动物 clade 传 `ND5,ND4L`（环节动物 benchmark 7/7 正确；其他类群 173 条里错 52 对 17）；`--aggressive_start` 对所有基因如此 |
| P2 | (ii) 同链同框 PCG 不可重叠（不同框可以，如 atp8-atp6、nad4L-nad4） | 只报告 |
| P3 | (iii) 第一个完整终止密码子，或紧接下游 tRNA 前的 T/TA 不完整终止（tRNA 间断模型） | 完整终止落在下游同链 tRNA 内 ≥3 nt → 截到 tRNA 前并加 `/transl_except` 与 polyA note（真值：伸入 1–2 nt 时 NCBI 保留完整终止 130:57，≥5 nt 时多用 T/TA）；完整终止仅与 tRNA 重叠 1–2 nt 时保留；下游无 tRNA 而完整终止伸入下一个同链蛋白基因 > 30 nt（atp8-atp6、nad4L-nad4 除外；真值：11–30 nt 内保留完整终止 17:6，更深则 T/TA 为主）→ 不延伸：保留原 T/TA 末端（如 cox1→cox2），否则取下一基因前的 T/TA（如 nad6→cob）；核心工具已注释的、伸入下一基因的完整终止不回缩；已是 T/TA 结尾的基因不会为了 90 nt 内的完整终止伸入其他特征（rRNA、反链 tRNA）> 10 nt；注释末端停在第一个完整终止之前、且中间没有任何已注释基因 → 延伸到该终止（≤ 600 nt；核心工具在分化大的类群丢了 C 端同源区：nad5/nad4/nad2/cox1 短 100–500 nt；多类群 benchmark 52 例中 50 例延伸后与 NCBI 完全一致） |
| D1 | — | 同一蛋白基因在不同位置注释两次（MitoZ 遇重复区）→ 只保留完整、较长、靠前的一份；另一条链上与之重叠、不完整或含内部终止的同名拷贝（反义链假阳性）也删除 |
| P7 | — | 起点与终点不在同一读码框（内部终止，NCBI ERROR）→ 取两者中较长的开放读码框：保留起点到第一个终止，或保留终止、在其读码框内取离原起点最近的起始密码子；新 ORF 不足原长一半时不改，只报告（多类群 benchmark：内部终止 14 → 5） |
| P4 | — | 内部终止密码子、非起始密码子首位、多段/部分 CDS → 只报告，不改 |
| R1 | (vi) rRNA 边界不扩展到两侧基因 | 与两侧 tRNA/CDS/rRNA 的重叠切掉（默认 `--rrna_bounds trim`）；`fill` 再补到两侧基因（空隙 ≤ `--max_fill` 100 nt，GenBank 常见写法） |
| R3 | 命名 | rrnS / rrnL，12S / 16S ribosomal RNA；rRNA 与 CDS 重叠 > 50 nt 报告 |

tRNA 边界保持注释软件结果（准则 v）。每处修改写入 `curation_report.tsv`（规则、旧/新坐标、理由）。

**内置 NCBI 评估**：有 `table2asn` 时（PATH、`--table2asn` 或 mitointron 环境），每条记录在校正前后各用 NCBI 校验一次，`validation_compare.tsv` 按"记录 + 错误代码 + 特征"列出每类问题 fixed / remaining / new，终端汇总一行（如 `ERROR 128 -> 4 … fixed 124, remaining 4, new 0`）；"new" 即校正引入的新问题，应为 0。`--no_validate` 关闭。

```bash
python3 mito_curate.py -i annotated_gb_dir -o curated_dir --genetic_code 5 [--rrna_bounds trim|fill]
```

**验证**（31 条鳞沙蚕类真实 MitoZ 注释，/mnt/n/DLC_mt_2023_1016/Mitoz_annotate，对照 NCBI 原记录；坐标换算到 NCBI 序列）：

| | MitoZ 原始 | 校正后 |
|---|---|---|
| PCG 5′ / 3′ 与 NCBI 相同（284 个） | 253 / 211 | **265 / 254** |
| 12S 5′ / 3′（22 个） | 5 / 2 | **17 / 15** |
| 16S 5′ / 3′（22 个） | 0 / 0 | **17 / 11** |
| rRNA 与相邻基因重叠 | 44 | **0** |
| table2asn ERROR | 128（126 个 NoStop） | **4**，全部对应报告中标出需人工检查的 3 个基因（跨原点 cox1、含 24 个内部终止的 COX1、tRNA 前无终止的 ATP6） |

剩余差异主要来自 tRNA 边界与 NCBI 不同（NCBI 的 T/TA 终止落在其 tRNA 前），以及 rRNA 写法惯例（NCBI 多为两侧填满；trim 模式按准则不扩展）。准则来自腹足类，此处在环节动物上验证，软体动物数据尚未测试。

## tRNA 多方法注释 / tRNA (mito_trna.py)

仿 GeSeq（Tillich et al. 2017, NAR 45:W6, doi:10.1093/nar/gkx391）：多个 tRNA 软件结果**全部**写进一个 GenBank（`review/<样本>.gbf`，Geneious 里逐个比较、自己选），另出 `selected/<样本>.gbf`（每位点一个 tRNA）和 `trna_report.tsv`。不预设 tRNA 数目（多数珊瑚 1–2 个，部分环节动物 23–24 个）。

- 软件：MITOS2（只取 tRNA，refseq89m 内置后生动物 tRNA 协方差模型）、ARWEN、ARAGORN、tRNAscan-SE 2（-O）、输入 GB 已有 tRNA。
- 证据按**独立方法家族**计：ARWEN 与 ARAGORN 同作者同算法，算一个。
- 默认选择：MITOS2 找到的采用；MITOS2 没找到的需 ≥2 个独立家族支持。同位点不同反密码子（冲突）保留证据强者，平局都不选。边界取 MITOS2（与 NCBI 完全一致最多）。
- MITOS2 慢的原因：runmitos 对没找到的 tRNA 模型再做灵敏搜索（cmsearch --max、E 0.1、单线程），tRNA 少的类群每条 10 分钟以上并带来假 tRNA。默认 fast 模式（`scripts/runmitos_fast.py`）跳过这一步（一条海葵 640 s → 8 s）；`--mitos_mode full` 恢复。
- 跨原点 tRNA 写成 join。

**用 MITOS2 tRNA 替换 MitoZ tRNA**（`--select mitos2`，可选，不在默认流程中）：

```bash
# 默认流程：1 batch_mitoz.py run → 3 mito_curate.py（含 4 NCBI 校验）
# 需要 MITOS2 tRNA 时在 1 与 3 之间插入：
python3 mito_trna.py -i mitoz_run_out/04.final_gb -o trna_mitos2 --select mitos2 --genetic_code 5 \
  --mitos_bin MITOS_ENV/bin --mitos_refdir DIR_WITH_refseq89m
python3 mito_curate.py -i trna_mitos2/selected -o curated --genetic_code 5
```

- 输入必须是**调整方向后的最终 GenBank**（`04.final_gb/`）：MITOS2 在该 GB 自身序列上运行，坐标与其它特征一致；FASTA、`00.*`、`01.mitoz_anno1/`（方向调整前）一律拒绝。
- 只替换 tRNA 及其 gene 特征（命名同 MitoZ，如 `trnN(guu)`），PCG/rRNA 不动；review GB 中保留原 MitoZ tRNA 供对照。
- 某条记录 MITOS2 运行失败 → 该记录原样保留，不会删掉 MitoZ tRNA。

```bash
python3 mito_trna.py -i gb_or_fasta -o out --genetic_code 5 --tools mitos2,arwen,aragorn,trnascan \
  --bin ENV/bin --mitos_bin MITOS_ENV/bin --mitos_refdir DIR_WITH_refseq89m
```

**验证**（裸序列，对照 NCBI tRNA，氨基酸 + 反密码子 + 重叠 ≥50% 算找到；FP = 该处无 NCBI tRNA；`benchmark/eval_trna.py`）：

| | 刺胞/海绵/扁盘 85 条，183 tRNA | 环节动物 RefSeq 40 条，880 tRNA |
|---|---|---|
| MITOS2 | 96.7%，FP 21 | 97.3%，FP 3 |
| ARWEN / ARAGORN | 91.8%，FP 348 / 462 | 95.7% / 96.5%，FP 79 / 117 |
| tRNAscan-SE | 61.2%，FP 1 | 46.1%，FP 0 |
| **默认选择** | **96.7%，FP 21** | **97.5%，FP 3** |

记录清单 `benchmark/trna_sets.tsv`（2026-10-20 重跑）。旧版把 MITOS2 的 `trnM_0(cat)`（同一模型多拷贝）和 `trnQ(---)`（反密码子读不出）丢掉，刺胞 MITOS2 只有 89.6%、默认 93.4%；同一批 85 条修复后 96.7%（另三种工具数字不变，差异全部来自此修复）。

单用 ARWEN/ARAGORN 在 tRNA 少的类群假阳性极多，故只放 review 文件。刺胞的 FP 多来自 MITOS2，是否为 NCBI 漏注未确证。

## intron 注释 / Introns (mito_intron.py)

**MitoZ 缺陷：MitoZ 无法注释 intron（group I / II）。** `run` / `annotate` 写完最终 GenBank 后自动调用 `mito_intron.py` 补上 intron，并整合进 MitoZ 结果（`--intron_fix auto`：找到 exonerate+blastp 时运行；`yes` 强制，`no` 关闭；`--intron_bin <env>/bin`）。

`mito_intron.py` = MFannot（Beck & Lang，https://github.com/BFL-lab/Mfannot ；数据 https://github.com/BFL-lab/MFannot_data ，GPL-3.0；Lang et al. 2023, doi:10.3389/fpls.2023.1222186）intron 流程的 Python 移植，本文件按 GPL-3.0 发布：

| 步骤 | MFannot 对应 | 实现 |
|---|---|---|
| 蛋白库 | protein_collections/identified | 随附 `models/mfannot_proteins/`（313 个基因家族，基因名取自 `mt cox1;` 表头） |
| 预筛 | flip + blast + SelectProteinForExo | ORF ≥40 aa（首个起始密码子起）→ blastp（BLOSUM62，E ≤1e-10）→ 按基因名分组 → 同链相隔 ≥142 nt 的 ≥2 个 ORF 才进入下一步 |
| 基因模型 | FindExonsInHypProtArray / SortAndGroupC4Report | Exonerate protein2genome（minintron 142、maxintron 5000、intronpenalty -35、`models/neutral.pssm`），分数 ≥200 的比对片段按 MFannot 规则串联，可跨越 >5 kb、内含其他基因的 intron |
| 类型 | RNAweasel | Infernal + Rfam RF00028 / RF00029（`models/rfam_group_I_II_introns.cm`） |
| 边界 | Adjust_all_intronic_junctions | ±6 nt 同移，group I（外显子末端 U、intron 末端 ωG）≥19，group II（GUGYG…AY）≥40 |
| 整合 MitoZ | — | 基因起止优先采用 MitoZ 的同框 CDS 端点；无则按同源蛋白延伸到起始/终止（含不完整终止 T/TA） |

与 MFannot 的差异（均写在脚本开头）：类型判定用 Rfam；边界移位以同源蛋白判断（MFannot 以起始模型判断，Exonerate 起始模型错时会拒绝正确边界）；跨环状原点的基因、外显子在两条链上的反式剪接基因不建模（报告中注明）；未移植 MFannot 的小外显子 HMM 步骤。

**验证**（`benchmark/`，2026-10-09）：NCBI RefSeq 中全部带 intron 标注的动物线粒体基因组 85 条（刺胞 82、海绵 2、扁盘 1；真值 138 个 intron = 126 个 intron 特征 + 12 个由 CDS join 推出但缺 intron 特征的），阴性对照环节动物 20 条（含围沙蚕属 4 条）。输入模拟 MitoZ 输出（去掉 intron、拆分基因只留外显子 1）。

| | 完全一致 | ≤6 nt | >6 nt | 漏检 | 多报 | 阴性假阳性 |
|---|---|---|---|---|---|---|
| 真 MFannot（Docker nbeck/mfannot） | 65 | 13 | 7 | 53 | 0 | 0 |
| mito_intron.py | **97** | 23 | 14 | **4** | 1 | **0** |
| mito_intron.py，排除本物种蛋白 | 97 | 23 | 14 | 4 | 1 | 0 |

MFannot 找到的 85 个 intron 中 83 个与本脚本坐标完全相同；MFannot 漏检主要是 >5 kb 的 nad5 大 intron（38 个中漏 38 个，本脚本找到 37 个）。与 NCBI 不一致的 41 处中 19 处 MFannot 与本脚本完全相同（NCBI 边界/终止密码子写法差异）。已知错误：扁盘动物 *Polyplacotoma* 反式剪接 COX1 产生一个 17.7 kb 假 intron（注释带 "UNCONFIRMED … check manually"）；海绵 *Plakinastrella* 两个 intron 间 9 nt 小外显子被合并。**intron 模型仍需人工核对。**

**环节动物**（2026-10-09，NCBI 全部环节动物线粒体全基因组 900 条，去重后 733 条，遗传密码表 5）：环节动物线粒体 intron 主要是 cox1 的 **group II** intron（Rfam RF00029），与刺胞动物的 group I 不同。
- NCBI → 本脚本：NCBI 有 18 个 intron（3 个 intron 特征 + 15 个只体现在 CDS join 中），找到 15 个；漏 3 个（*Decemunciger* ND4 265 nt 断开处，MFannot 亦未检出）。
- 本脚本 → NCBI：21 个预测；3 个对应 NCBI intron 特征，12 个对应 NCBI CDS join，6 个在 NCBI 完全无注释的记录中（*Glycera fallax* ×2、臭海蛹 *Travisia forbesii*、*Chaetopterus dewysee*、*Notomastus* sp. B、*Podarkeopsis*）。
- 与真 MFannot：13 个坐标完全一致（含上述无注释记录中的 6 个里的 5 个，以及臭海蛹 *T. sanrikuensis* / *T. forbesii*）；*Decemunciger* ND1 MFannot 判为 2 个 intron（中间 19 nt 小外显子，本脚本合为 1 个）；*Chaetopterus* 5 个约 300 nt、类型未定的预测 MFannot 未报，**未确证**。
- 臭海蛹属 NCBI 线粒体全基因组 2 条：*T. sanrikuensis* LC677172（cox1 intron 9153..10034，NCBI 9152..10033）、*T. forbesii* OZ479805（NCBI 无注释；cox1 intron 7346..8597，与 MFannot 一致）。

### v5（2026-10-10）：NCBI 提交合规 + 证据门槛

- **NCBI 校验**：`ncbi_validate.py` 把 GenBank 转为 NCBI 提交格式（.fsa + 5 列 .tbl），用 NCBI `table2asn -V vb` 校验（`--only_ours` 只统计本脚本写入的特征）。
- **合规**：新 CDS 带 `/product`（沿用 MitoZ，否则标准名）；T/TA 不完整终止带 `/transl_except=(pos:…,aa:TERM)` 与标准 note；3′ 端延伸到阅读框内第一个终止密码子（允许与相邻基因重叠，`--stop_rule first`，默认）。
- **证据门槛**（不依赖 NCBI 原注释）：拼接后 CDS 有内部终止、起始密码子非法、未定型 intron 剪接分数 < 14、或与已有的另一标准线粒体蛋白基因 CDS 冲突（内含子编码 ORF 如 heg 除外）→ 不写入 GenBank，只在 `intron_report.tsv`（`confidence` / `written` / `evidence` 列）中列出。
- **小外显子**：多个同源蛋白在外显子交界处共同显示缺失时，在 intron 内搜索插入后对所有同源蛋白比对都改善的短外显子（MFannot AnnotateMiniExons 思路）。

| 指标 | 刺胞/海绵/扁盘 85 条（138 intron） | 环节动物 733 条（NCBI 18 intron） |
|---|---|---|
| 本脚本写入特征的 NCBI 校验 ERROR | **0**（v3：202） | **0** |
| intron 与 NCBI 完全一致 / ≤6 nt / >6 nt / 漏检 | 98 / 23 / 6 / 11 | 找到 11，漏检 7 |
| CDS 与 NCBI 完全一致 | 99/126（3′ 端 126/126） | — |
| 阴性记录假阳性 | 0 | 写入 16 个：11 个有 NCBI 支持，5 个在 NCBI 完全无注释的记录中（均与真 MFannot 坐标一致） |

与 NCBI 不一致的 28 处 intron 边界，用独立证据（同源蛋白比对得分、group I ωG 规则）判定：13 处两种边界蛋白相同（写法差异），14 处本脚本更好（多为石珊瑚 COX1，NCBI 边界不符合 ωG 规则），1 处 NCBI 更好。NCBI 原记录自身在 table2asn 下也有错误（刺胞 13 个 `MultiIntervalIntron`，环节 1 个 `AbuttingIntervals`），本脚本不照抄 NCBI 坐标。环节动物漏检 7 个：*Decemunciger* ND4 ×3（MFannot 亦漏）、*Chaetopterus* ×4（约 300 nt、未定型、剪接位点弱、MFannot 未报，NCBI 仅以 CDS join 表示；是否为 intron **未确证**，按门槛不写入）。

**近缘参考（借鉴 MitoFinder）**：`--ref_gb <近缘物种 GenBank>`（`batch_mitoz.py --intron_ref_gb`）把参考 CDS 翻译加入同源蛋白库。MFannot 自带库以真菌/原生生物为主，动物极少、无环节动物，近缘参考可提高环节动物等类群的检出。MitoFinder 本身的 intron 注释只靠 BLASTX 相似性合并片段、不搜边界（其 README 明言），故不替换现有方法。

### v6（2026-10-10）：结构模型扩充 + 裸序列验证

- **类型模型**：Rfam RF00028/RF00029 + 7 个 group II D1–D4 模型（RF01998, RF01999, RF02001, RF02003–05, RF02012；随附）+ Sellés Vidal et al. 2025 的 14 个 group I 亚型模型（IA1–IE3，基于 42,459 个 group I intron 训练，doi:10.1093/bioadv/vbaf020；无许可声明故不随附，首次运行下载到 `~/.cache/mito_intron/`，`--no_subtype_cms` 关闭）。注释写出亚型（如 IB1）。NCBI 真值 intron 判型：刺胞 138 个由 74 → 109 个，环节 27 个不变（18 个 group II）；273 条编码序列 0 误判。
- **起始/不完整 CDS**：首密码子非起始密码子时上游（再下游）搜索起始密码子；阅读框延伸出序列末端时写为 `<`/`>` 不完整 CDS（NCBI 允许），不再因"无起始密码子"丢弃。
- **借鉴 MitoFinder**：`--ref_gb` 近缘参考；外显子总长 ≤ 3×蛋白长+33。
- **其他工具实测/调研**：MITOS2（mamba，refseq89m）在 Metridium 上把 cox1/nad5 注成 -a/-b 片段、外显子边界偏差最多 174 nt、无 intron 特征；在臭海蛹 cox1 上只注出第二外显子。MitoFinder `--allow-intron` 仅按 BLASTX 相似性合并、不搜边界（其 README）。MitoZ、MITOS2、GeSeq、PGA、CPGAVAS2 等均无 intron 模型；目前只有 MFannot/RNAweasel 建模 intron（调研见 benchmark 记录）。

| 验证 | 结果 |
|---|---|
| 刺胞/海绵/扁盘 85 条（138 intron，模拟 MitoZ 输入） | 完全一致 98 / ≤6 nt 23 / >6 nt 6 / 漏检 11（同 v5）；高可信 126、中 3、低 4（v5：105 / 22 / 6）；本脚本写入特征 table2asn ERROR 0 |
| 环节动物 741 条（**裸序列**，去掉全部注释，不用 MitoZ） | NCBI→本脚本：27 个 intron 中完全一致 6 / ≤6 nt 6 / >6 nt 4 / 漏检 11（4 个 intron 伸到片段末端无第二外显子、3 个 *Decemunciger* ND4 与 4 个 *Chaetopterus* 约 300 nt 未定型，MFannot 均未检出）；本脚本→NCBI：写入 21 个，**0 个与 NCBI 矛盾**（8 intron 特征、8 CDS join、5 个在 NCBI 无注释记录中且与 MFannot 一致），全部高可信 |

### v7（2026-10-10）：提速与稳健性

- **瓶颈实测**（cProfile）：v6 的时间几乎全在 cmsearch，且只集中在命中的少数模型（IB3 12 s、RF00028 10 s）；与 intron 长度无关（只扫两端 1.5 kb 不更快，判型也不变）。
- **两步判型**：先用 Infernal `--rfam` 快速过滤，仍未定型的再全灵敏度搜索；同一记录所有 intron 一次调用；只搜正链（`--toponly`）、`--threads` 多线程。165 个 NCBI intron 判型数与大类与全灵敏度完全一致（147/165）。
- **小外显子两步搜索**：先单同源蛋白、交界窗口粗筛，再对前 5 个候选完整验证。
- **稳健性**：外部程序超时（`--timeout`，默认 900 s）不再卡死；单条记录出错不影响整批（原样保留并在报告中注明）；`--jobs N` 多文件并行（batch_mitoz 传 `--max_tasks`）。
- **效果**：单记录 *Porites* 42 → 6 s、*Plakinastrella* 103 → 8.5 s、*Travisia forbesii* 4.2 → 1.2 s；全部 846 条验证记录 12 分钟；坐标与 v6 完全相同，判型 +1（NC_040136 nad5，IB1）；测试 464 → 107 s。

依赖：`conda create -n mitointron -c conda-forge -c bioconda exonerate blast infernal biopython`。

```bash
python3 mito_intron.py -i final_gb_dir -o out_dir --genetic_code 4 --bin /root/miniconda3/envs/mitointron/bin
python3 -m pytest -q tests          # Metridium / Nematostella / Platynereis，排除本物种蛋白
```

## 注意事项 / Practical Notes

1. `--clade` 默认 `auto`：从 GenBank 分类信息判断（一次运行只能一个门；分类信息里没有已知门名、含多个门或各条记录 `transl_table` 不一致时报错）。完全没有分类信息的输入（FASTA、无 lineage 的 GenBank）沿用 `Annelida-segmented-worms` 并给出警告；其他类群请显式给 `--clade` 和 `--genetic_code`。
2. `--genetic_code` 没有固定默认值：`--clade auto` 时取输入 CDS 的 `transl_table`，否则按 `--clade`（Chordata 2、Echinodermata 9、Porifera 4，其余 5）；请按物种确认。
3. `--id_maxlen` 不能小于 6，否则内部 ID 会冲突；脚本现在会拒绝这种参数。
4. 线性序列重定向后会保留 `linear` topology，不会被强制改成 `circular`。
5. 多 contig 样本中缺少目标基因的 contig 会合并回最终结果，不会静默丢弃。
6. 如果 MitoZ 不在 PATH 中，使用 `--conda_env_mitoz` 或 `--mitoz_path`。
7. 如果只需要改名 FASTA header，不需要运行本脚本；可以直接处理 FASTA 文本。

Notes:

1. The default `--clade` is `auto` (from the GenBank taxonomy, one phylum per run; a taxonomy without a known phylum, several phyla, or differing CDS `transl_table` values are errors). Input with no taxonomy at all (FASTA, GenBank without a lineage) falls back to `Annelida-segmented-worms` with a warning; for any other taxon give `--clade` and `--genetic_code`.
2. `--genetic_code` has no fixed default: with `--clade auto` it is the CDS `transl_table` of the input, else MitoZ's code of the `--clade` (2 Chordata, 9 Echinodermata, 4 Porifera, else 5); verify it for your dataset.
3. `--id_maxlen` must be at least 6.
4. Linear records keep linear topology after reorientation.
5. Partial no-target-gene contigs are merged back into final outputs.
6. Use `--conda_env_mitoz` or `--mitoz_path` if MitoZ is not in PATH.
7. Simple FASTA header renaming does not require this MitoZ pipeline.

## 流程图 / Flowcharts

已生成以下流程图文件：

- 中文 SVG: [docs/flowcharts/batch_mitoz_workflow.zh.svg](docs/flowcharts/batch_mitoz_workflow.zh.svg)
- English SVG: [docs/flowcharts/batch_mitoz_workflow.en.svg](docs/flowcharts/batch_mitoz_workflow.en.svg)
- 中文 Mermaid: [docs/flowcharts/batch_mitoz_workflow.zh.mmd](docs/flowcharts/batch_mitoz_workflow.zh.mmd)
- English Mermaid: [docs/flowcharts/batch_mitoz_workflow.en.mmd](docs/flowcharts/batch_mitoz_workflow.en.mmd)
- SVG 生成脚本: [scripts/render_workflow_diagrams.py](scripts/render_workflow_diagrams.py)

重新生成 SVG：

```bash
python3 scripts/render_workflow_diagrams.py
```

### 中文流程图 / Chinese Workflow

```mermaid
flowchart TD
    A["输入数据<br/>FASTA / GenBank；FASTQ 双端可选"]
    B["输入识别与依赖检查<br/>Biopython、MitoZ、Snakemake"]
    C["GenBank 去重与转换<br/>保留 B64_2 与 B64 为不同编号"]
    D["生成内部 FASTA<br/>S00001...；记录 id_map.tsv"]
    E["第一轮 MitoZ 注释<br/>输出 01.mitoz_anno1/ 和 sample.log"]
    F["按目标基因重定向<br/>默认 cox1；保留 linear / circular"]
    G["可选第二轮重注释<br/>if_reannotation_2=yes 时每条运行；auto 时只对未通过平移检查的运行"]
    H["最终 GenBank 生成<br/>ID 重写、metadata 恢复、/organism= 同步"]
    I["统计与问题清单<br/>annotation-summary.tsv、failed/mismatch 列表"]
    J["合并最终结果<br/>04.final_gb/*.gbf 与 all.final.gbf"]
    A --> B --> C --> D --> E --> F --> G --> H --> I --> J
```

### English Workflow

```mermaid
flowchart TD
    A["Input data<br/>FASTA / GenBank; optional paired FASTQ"]
    B["Detect input and check dependencies<br/>Biopython, MitoZ, Snakemake"]
    C["Deduplicate and convert GenBank<br/>Keep local IDs such as B64_2 distinct from B64"]
    D["Create internal FASTA files<br/>S00001... plus id_map.tsv"]
    E["MitoZ annotation, phase 1<br/>Writes 01.mitoz_anno1/ and sample.log files"]
    F["Reorient by target gene<br/>Default cox1; preserve linear / circular topology"]
    G["Optional phase-2 re-annotation<br/>if_reannotation_2=yes: every record; auto: only records failing the shift check"]
    H["Build final GenBank files<br/>Rewrite IDs, restore metadata, synchronize /organism="]
    I["Write statistics and problem lists<br/>annotation-summary.tsv and failed/mismatch files"]
    J["Merge final outputs<br/>04.final_gb/*.gbf and all.final.gbf"]
    A --> B --> C --> D --> E --> F --> G --> H --> I --> J
```
