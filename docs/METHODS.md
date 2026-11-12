# MitoM 方法与评测 / Methods and benchmarks

MitoM 各模块的方法、设计依据和历次评测（2026-10）。使用见 [USAGE_zh.md](USAGE_zh.md)，示例见 [EXAMPLES.md](EXAMPLES.md)，参考库制备见 [DATASETS.md](DATASETS.md)。评测数据在 `../benchmark/results/`。

Methods, design rationale and benchmarks of each MitoM module (2026-10). Usage: [USAGE_zh.md](USAGE_zh.md); examples: [EXAMPLES.md](EXAMPLES.md); preparing reference libraries: [DATASETS.md](DATASETS.md). Benchmark data: `../benchmark/results/`.

Note: a measured figure (benchmark result) is written once, in a table or, where there is no table, in the Chinese text; the English text states methods, rules, conditions and thresholds and does not repeat the figures.

## 程序结构 / Code layout

`mitom run` 做哪些步骤见 [README](../README.md) 第 3 节；本节只说代码在哪里。
What `mitom run` does step by step: [README](../README.md) section 3. This section says where the code is.

命令 `mitom` → `mitom_main.py`（检查环境，启动子命令）→ `mitom_pipeline.main`（只分发命令）→ `mitom_run.RUN_STAGES`（`run` / `import` 的 12 步）。
`mitom` → `mitom_main.py` (environment check, starts the subcommand) → `mitom_pipeline.main` (dispatch only) → `mitom_run.RUN_STAGES` (the 12 steps of `run` / `import`).

| 文件 / file | 内容 / contents |
|---|---|
| `mitom_cli.py` | 命令行参数 / command-line options |
| `mitom_prepare.py` | 环境检查、工具路径、`--clade auto`、输出目录锁 / environment, tools, clade plan, output lock |
| `mitom_run.py` | `run` / `import` 的 12 步，共享 `RunState` / the 12 stages |
| `mitom_steps.py` | 校正、内含子、未知 ORF、NCR、NCBI 文件等步骤 / curation, introns, ORFs, NCR, submission steps |
| `mito_orf.py` | 未知功能的 ORF（同源证明受选择才写）/ ORFs of unknown function |
| `mitom_groups.py` | 多个门或遗传密码时分组运行 / one run per phylum / genetic code group |
| `mitom_commands.py` | 其他子命令（reorient、dedup、intron …）/ the other subcommands |
| `mitom_context.py` | 一次运行的共享状态 `RunContext` / per-run state |
| `mito_core.py`、`mito_trna.py` | 注释核心、tRNA / annotation core, tRNAs |
| `mito_curate*.py` | 注释后校正：`_pcg` 蛋白基因、`_origin` 跨起点、`_rna` rRNA / tRNA / 命名、`_base` 常量 / post-annotation curation by rule family |
| `mito_intron*.py` | 内含子：从头（`mito_intron`、`_rules`）、位点迁移（`_ref`、`_junction`、`_models`、`_sites`）/ introns |
| `mito_nuc*.py` | 核 rDNA、组蛋白 / nuclear rDNA, histones |

I / II 型内含子剪接位点的判定有三处，用途不同、互不调用：`mito_intron_rules.py`（从头识别，MFannot 打分的移植，GPL）、`mito_intron_junction.py` 的 `rule_score`（位点迁移时给候选接头打分）、`mito_intron_sites.py` 的 `canonical_kinds`（金标准的是 / 否判定）。改剪接规则时三处都要看；合并它们会改变结果，需先跑内含子金标准（TEST_REGISTRY 第 68 行）。
Group I / II splice sites are judged in three places with different jobs, none calling another: `mito_intron_rules.py` (de novo recognition; a port of MFannot's scoring, GPL), `rule_score` in `mito_intron_junction.py` (scores candidate junctions during site transfer), and `canonical_kinds` in `mito_intron_sites.py` (the yes / no judgement for the gold standard). When changing the splice rules, check all three; merging them changes the results, so run the intron gold standard first (TEST_REGISTRY line 68).

`mitom annotate` 用同一套阶段函数（`mitom_run.ANNOTATE_STAGES`：不做第二轮、对比和输出整理），后处理只有一份代码。
`mitom annotate` runs the same stage functions (`ANNOTATE_STAGES`: no round 2, comparison or output layout), so the steps after annotation exist once.

校正、内含子与未知 ORF 以独立进程运行（`mito_curate.py`、`mito_intron.py`、`mito_intron_ref.py`、`mito_orf.py`），崩溃只影响该步。NCR、写回、圈图等在主进程内运行；`tests/test_step_argv.py` 检查流水线给子进程的参数。旧名 `batch_mitoz.py` 仍可用。
Curation, introns and ORF run as separate processes; NCR and plots stay in-process. `tests/test_step_argv.py` checks subprocess arguments. The old name `batch_mitoz.py` still works.


## 类群与遗传密码自动判断 / --clade auto (mito_classify.py, mitom_groups.py)

`--clade auto`（默认）逐条记录决定门、类群和遗传密码，结果与依据写入 `qc/clade_decision.tsv`：

`--clade auto` (default) decides the phylum, clade and genetic code of each record; the result and its basis are written to `qc/clade_decision.tsv`:

1. **GenBank 记录**：分类信息里的门（或纲级同义名）；遗传密码取 CDS 的 `/transl_table`（多数票），没有 CDS 时取参考基因组中该谱系最深一级类群的遗传密码（`code_of_lineage`，例如 Tunicata → 13）。分类信息不属于 Metazoa 时报错。

   **GenBank records**: the phylum (or a class-level synonym) from the taxonomy; the genetic code from the CDS `/transl_table` (majority vote), or, if there is no CDS, the genetic code of the deepest clade of the lineage in the reference genomes (`code_of_lineage`, e.g. Tunicata -> 13). A taxonomy outside Metazoa is an error.
2. **FASTA / 无分类信息的记录**：blastx（`-query_gencode 5`，任何动物线粒体密码下都不会提前终止）比对 refseq89m 的 8 个保守蛋白（cox1 cox2 cox3 cob atp6 nad1 nad4 nad5，约 8,000 个基因组）。参考基因组的谱系和遗传密码来自 `models/ref_taxa.tsv.gz`（`scripts/build_ref_taxa.py` 用本机 RefSeq 生成，8,022 个中 7,672 个有谱系）。

   **FASTA / records without taxonomy**: blastx (`-query_gencode 5`, which never terminates early under any animal mitochondrial code) against 8 conserved proteins of refseq89m (cox1 cox2 cox3 cob atp6 nad1 nad4 nad5; about 8,000 genomes). The lineage and genetic code of the reference genomes come from `models/ref_taxa.tsv.gz` (made from the local RefSeq by `scripts/build_ref_taxa.py`).
   - 门：每个基因取前 5 个命中，按（比分 / 该基因最高比分）加权投票；支持度 = 获胜门所占权重。

     Phylum: for each gene the top 5 hits vote, weighted by (score / the gene's highest score); support = the share of the weight held by the winning phylum.
   - 遗传密码：获胜门的命中的 `transl_table` 加权投票。

     Genetic code: weighted vote of the `transl_table` of the hits of the winning phylum.
   - 密码子核验：在这些比对中，基因组自己的 AGA AGG ATA AAA TGA TAA 处对应的参考氨基酸；每个动物线粒体密码（2 4 5 9 13 14 21 24 33）按吻合比例打分。AGR 在比对区内出现（≥ 1/1000 个密码子且多于 3 倍 TAG）时排除把 AGR 当终止的密码 2（文昌鱼 3.0/1000，脊椎动物 0–0.1/1000）；出现次数不多于 3 倍 TAG（所有密码中的终止密码子，比对区内的 TAG 即测序错误 / 移码的噪声）的密码子不计分。打分明显更好（差 ≥ 0.15）的密码取代近亲的密码，并给出警告。

     Codon check: in these alignments, the reference amino acid opposite the genome's own AGA AGG ATA AAA TGA TAA; each animal mitochondrial code (2 4 5 9 13 14 21 24 33) is scored by the fraction of agreement. When AGR occurs in the aligned region (>= 1/1000 codons and more than 3 times TAG), code 2, which treats AGR as a stop, is excluded (amphioxus has AGR in coding regions, vertebrates almost none). A codon that occurs no more than 3 times as often as TAG (TAG is a stop in all codes, so TAG in the aligned region is noise from sequencing errors / frameshifts) is not scored. A code that scores clearly better (difference >= 0.15) replaces the code of the close relatives, with a warning.
   - 支持度 < 0.7（参考库里没有近亲，小门常见）：门判为不确定，类群用 `Animal`，遗传密码用密码子核验的结果。

     Support < 0.7 (no close relative in the reference library, common in small phyla): the phylum is called uncertain, the clade is `Animal`, and the genetic code is the result of the codon check.
   - 命中少于 2 个基因（< 约 4 kb 的片段）：不判断，并入输入中最常见的组。

     Fewer than 2 genes hit (a fragment shorter than about 4 kb): no decision; the record joins the most common group of the input.
3. **多个门 / 遗传密码**：`mitom run` 按（类群，遗传密码，门）分组，每组以显式 `--clade --genetic_code --phylum` 运行一次完整流程（并行；`--jobs` 按记录数分给各组），结果合并到同一输出目录（`annotation/`、`summary.tsv`、`submission/`、`annotation/plots/`、`qc/` 表格合并；各组参数与日志在 `logs/groups/`）。

   **Several phyla / genetic codes**: `mitom run` groups records by (clade, genetic code, phylum) and runs the full pipeline once per group with explicit `--clade --genetic_code --phylum` (in parallel; `--jobs` is divided among the groups by number of records). The results are merged into one output directory (`annotation/`, `summary.tsv`, `submission/`, `annotation/plots/`; the tables in `qc/` are merged; each group's parameters and logs are in `logs/groups/`).
4. 门未知的记录（FASTA、门不确定）缺 rRNA 时，rRNA HMM 后备步骤使用所有门的模型中得分最高的那一个门。

   For records of unknown phylum (FASTA, uncertain phylum) that lack an rRNA, the rRNA HMM fallback step uses the phylum whose models score highest among all phyla.

**评测**（2026-10-31，均为本机数据，未新下载）：

**Benchmark** (2026-10-31, all on local data, nothing newly downloaded):

| 测试 / test | 记录 / records | 遗传密码正确 / genetic code correct | 门 正确 / 不确定 / 错误 / phylum correct / uncertain / wrong |
|---|---|---|---|
| RefSeq 分层抽样（每门 ≤ 15 条 + 各门少数密码 ≤ 3 条，29 门），去掉同科参考 / RefSeq stratified sample (<= 15 per phylum + <= 3 of each phylum's minority code; 29 phyla), same-family references removed | 238 | 93.3 % | 63.0 / 32.4 / 4.6 % |
| 同上，8 kb / 5 kb / 3 kb 随机片段 / as above, random fragments of 8 kb / 5 kb / 3 kb | 238 × 3 | 91.6 / 85.3 / 74.4 % | 58.4 / 50.0 / 37.8 % 正确 / correct |
| 同上，10 % / 20 % 随机替换 / as above, 10 % / 20 % random substitutions | 238 × 2 | 89.1 / 84.0 % | 63.0 / 60.9 % 正确 / correct |
| 同上，20 % 的 N（200 nt 块）、反向互补 + 随机起点 / as above, 20 % N (200-nt blocks); reverse complement + random start | 238 × 2 | 93.3 / 93.3 % | 与完整基因组相同 / same as the full genome |
| 全流程：55 条 29 门无标签 FASTA，自身参考在分类和注释中都排除（`MITOM_EXCLUDE_REF`） / full pipeline: 55 unlabelled FASTA records of 29 phyla, the record's own reference excluded in both classification and annotation (`MITOM_EXCLUDE_REF`) | 55 | 52 / 55 | 29 / 23 / 3；NCBI ERROR 0 |

门判错的都是参考库里只有 1–2 个科的小门（有爪动物、鳃曳动物、线形动物判成节肢动物）；它们的遗传密码仍正确。遗传密码的残余错误：吸虫 14 → 9（两者只差 TAA）、柱头虫 9 / 24、六放海绵 5 / 4。两个门的嵌合序列判为其中之一或不确定（不报警）。知道类群时请显式给 `--clade` 和 `--genetic_code`。

Wrong phylum calls are all in small phyla with only 1–2 families in the reference library (Onychophora, Priapulida and Nematomorpha called as arthropods; their genetic codes are still correct). Residual genetic-code errors: trematodes 14 -> 9 (the two differ only at TAA), acorn worms 9 / 24, hexactinellid sponges 5 / 4. A chimeric sequence of two phyla is called as one of them or as uncertain (no warning). When the clade is known, give `--clade` and `--genetic_code` explicitly.

## 稳健性测试 / Robustness tests (2026-10-31)

- **模糊测试**：19 个命令 × 24 种输入（正常的 GenBank / FASTA / EMBL / 目录，以及空文件、垃圾文本、截断、无 ORIGIN、无特征、gzip、Latin-1、CRLF、蛋白序列、重复登录号、跨起点 / 模糊端点特征、短记录、不存在的路径）：所有命令对坏输入给出一行 `ERROR` 和退出码 2，没有 traceback，也没有"0 条记录"的静默成功（`mitom_work/fuzz_2026-10-31/scripts/fuzz.py`）。

  **Fuzz tests**: 19 commands x 24 kinds of input (normal GenBank / FASTA / EMBL / directories, and empty files, garbage text, truncation, no ORIGIN, no features, gzip, Latin-1, CRLF, protein sequences, duplicate accessions, features across the origin / with fuzzy ends, short records, non-existent paths): every command answers bad input with a one-line `ERROR` and exit code 2, with no traceback and no silent success on "0 records" (`mitom_work/fuzz_2026-10-31/scripts/fuzz.py`).
- **坐标变换**：5 条真实记录（内含子、跨起点、反式剪接）× 12 个随机旋转 × 反向互补：每个特征提取出的序列不变（`tests/test_orientation.py`）。

  **Coordinate transforms**: 5 real records (intron, across the origin, trans-splicing) x 12 random rotations x reverse complement: the sequence extracted for each feature is unchanged (`tests/test_orientation.py`).
- **校正**：把 55 条 29 门的 MitoM 结果分别旋转到起点落在每个基因内部再校正：基因集合不再变化（修复前 70 例被删、约 25 例被截半）。NCBI NC_001636.1 校正后 13 个蛋白的翻译与 NCBI 完全相同。

  **Curation**: 55 MitoM results of 29 phyla were each rotated so that the origin falls inside every gene in turn, then curated: the gene set no longer changes (before the fix, genes were deleted or cut in half). After curation of NCBI NC_001636.1, the translations of its 13 proteins are identical to NCBI's.
- **运行操作**：同一目录重跑、Ctrl-C / kill -9 后重跑、并发运行、threads 与 Snakemake 调度：结果相同；同一输出目录的第二个运行被拒绝（`.mitom.lock`）。

  **Run operations**: re-running in the same directory, re-running after Ctrl-C / kill -9, concurrent runs, threads and Snakemake scheduling: the results are identical; a second run in the same output directory is refused (`.mitom.lock`).
- 四个模块的代码审查报告 32 个问题（另有 2 个附带问题）：29 个复现、修复并加了回归测试（CHANGELOG）；另 3 个：脊椎动物 tRNA 的 L1/L2、S1/S2 按 MITOS2 2.1.10 的规则命名（所有类群相同，与 NCBI 脊椎动物记录的编号相反，有意如此），对已校正记录再次校正时起点会变（2 个，P11 依赖核心的比对注释；待定）。

  A code review of four modules reported issues; those reproduced were fixed and given regression tests (CHANGELOG). Three were not: the L1/L2 and S1/S2 names of vertebrate tRNAs follow the MITOS2 2.1.10 rule (the same for all clades, opposite to the numbering of NCBI vertebrate records; intentional); and, when an already curated record is curated again, the start changes (2 cases; P11 depends on the core's alignment annotation; pending).

## 注释后校正 / Post-annotation curation (mito_curate.py)

对注释好的 GenBank（MitoM 核心、MitoZ、MITOS 等的输出）单独运行，按 Fourdrilis et al. 2018 的软体动物（腹足类）注释准则修正 PCG 与 12S/16S，准则引自 Ghiselli et al. 2021, *Phil. Trans. R. Soc. B* 376:20200159（doi:10.1098/rstb.2020.0159，第 3 节）：

Run on its own on annotated GenBank files (output of the MitoM core, MitoZ, MITOS, etc.). It corrects the PCGs and 12S/16S according to the molluscan (gastropod) annotation criteria of Fourdrilis et al. 2018, as cited in Ghiselli et al. 2021, *Phil. Trans. R. Soc. B* 376:20200159 (doi:10.1098/rstb.2020.0159, section 3):

| 规则 / rule | 准则 / criterion | 做法 / method |
|---|---|---|
| P1 | (i) 起始密码子取离上游基因最近且不重叠者 / (i) take the start codon nearest to the upstream gene that does not overlap it | 起点与上游 tRNA/rRNA 或同框 PCG 重叠 → 后移到第一个起始密码子，但只在 15 nt 内有起始密码子时才移动，否则保留重叠（173 条多类群真值：≤15 nt 的移动对 31 错 12，更远的移动对 14 错 28）；起点放在一个完整同链 tRNA 的上游（基因前三分之一包住整个 tRNA，如蚌类 nad1 包住 trnL2）→ 移到该 tRNA 之后第一个同框起始密码子；首密码子非起始密码子时才向上游找；`--aggressive_start_genes` 列出的基因总是按准则找最上游——默认为空，`mitom run` 只在环节动物 clade 传 `ND5,ND4L`（环节动物 benchmark 7/7 正确；其他类群 173 条里错 52 对 17）；`--aggressive_start` 对所有基因如此 / if the start overlaps an upstream tRNA/rRNA or an in-frame PCG, move it back to the first start codon, but only when a start codon lies within 15 nt, otherwise keep the overlap; if the start is placed upstream of a whole same-strand tRNA (the first third of the gene encloses the whole tRNA, e.g. bivalve nad1 enclosing trnL2), move it to the first in-frame start codon after that tRNA; look upstream only when the first codon is not a start codon; the genes listed in `--aggressive_start_genes` are always moved to the most upstream start by the criterion - the default is empty, and `mitom run` passes `ND5,ND4L` only for the Annelida clade (right on the annelid benchmark; on the other clades it did worse on the multi-clade truth set); `--aggressive_start` does this for all genes |
| P11 | (i) 的补充：起始密码子的选择 / supplement to (i): choice of the start codon | P1 之后，在 P1 的起点前后 30 个密码子内、读码框连续的所有起始密码子中，取 NCBI 记录会选的那个：梯度提升树（`mito_start.py`，模型 `models/start_model/start_model.json`，纯 JSON 树，不用 pickle）根据密码子、与上游特征的距离及其类型和链、核心工具的比对（N 端缺多少残基、覆盖度）、基因、门打分。训练：1,683 条留一基准记录的候选（`benchmark/train_start_model.py`）。按科交叉验证：起点与 NCBI 完全一致 89.1 → 90.7 %（节肢 77.8 → 82.3 %，刺胞 79.8 → 83.2 %，软体 79.5 → 81.3 %，环节 91.5 → 92.9 %，脊索 97.3 → 97.9 %）；`--start_model no` 关闭 / after P1, among all start codons within 30 codons either side of P1's start that keep the reading frame continuous, take the one an NCBI record would choose: a gradient-boosted tree model (`mito_start.py`, model `models/start_model/start_model.json`, plain JSON trees, no pickle) scores by codon, distance to the upstream feature and its type and strand, the core tool's alignment (how many N-terminal residues are missing, coverage), gene and phylum. Training: candidates from 1,683 leave-one-out benchmark records (`benchmark/train_start_model.py`); evaluated by family-level cross-validation. `--start_model no` turns it off |
| P2 | (ii) 同链同框 PCG 不可重叠（不同框可以，如 atp8-atp6、nad4L-nad4） / (ii) same-strand in-frame PCGs must not overlap (different frames may, e.g. atp8-atp6, nad4L-nad4) | 只报告 / report only |
| P3 | (iii) 第一个完整终止密码子，或紧接下游 tRNA 前的 T/TA 不完整终止（tRNA 间断模型） / (iii) the first complete stop codon, or an incomplete T/TA stop immediately before the downstream tRNA (tRNA punctuation model) | 完整终止落在下游同链 tRNA 内 ≥3 nt → 截到 tRNA 前并加 `/transl_except` 与 polyA note（真值：伸入 1–2 nt 时 NCBI 保留完整终止 130:57，≥5 nt 时多用 T/TA）；完整终止仅与 tRNA 重叠 1–2 nt 时保留；下游无 tRNA 而完整终止伸入下一个同链蛋白基因 > 30 nt（环节动物 > 10 nt，`--abbreviated_pcg_overlap`，与原点处 O1 一致；atp8-atp6、nad4L-nad4 除外；真值：11–30 nt 内保留完整终止 17:6，更深则 T/TA 为主；本地全部 NCBI 记录中 11–30 nt 的情形：环节动物 T/TA 142:43，刺胞动物 5:109，软体动物 12:66）→ 不延伸：保留原 T/TA 末端（如 cox1→cox2），否则取下一基因前的 T/TA（如 nad6→cob）；下游特征前没有 T/TA 时，完整终止伸入它 ≤ 20 nt（`OVERLAP_FALLBACK`；预测的 tRNA 端点不精确）则采用；核心比对已到达、只伸入预测 tRNA ≤ 20 nt 的完整终止保留；核心工具已注释的、伸入下一基因的完整终止不回缩；已是 T/TA 结尾的基因不会为了 90 nt 内的完整终止伸入其他特征（rRNA、反链 tRNA）> 10 nt；注释末端停在第一个完整终止之前、且中间没有任何已注释基因 → 延伸到该终止（≤ 600 nt；核心工具在分化大的类群丢了 C 端同源区：nad5/nad4/nad2/cox1 短 100–500 nt；多类群 benchmark 52 例中 50 例延伸后与 NCBI 完全一致）；会把整个同链 tRNA 包进去、且 tRNA 前没有 T/TA 的基因 → 结束在该 tRNA 之前（3′ partial）；脊椎动物（`--abbreviated_overlap`，`mitom run` 对 Chordata 自动加）：完整终止伸入下游同链 tRNA 1–2 nt 时也改为 T/TA，但伸入下一个蛋白基因时保留完整终止（NCBI 脊椎动物记录 atp6 的 TAA 与 cox3 的 ATG 共用 A；1,719 条留一基准记录只重跑校正：两端一致 17,727 → 17,951，80.0 → 81.0 %；脊索动物多对 294、少对 70，辐鳍鱼 atp6 两种写法都有 104 : 60，鸟类 86 : 2，哺乳类 60 : 2；其他门不变；完全不用该规则 17,373） / a complete stop that falls >= 3 nt inside a downstream same-strand tRNA -> cut back to before the tRNA, with `/transl_except` and the polyA note (in the truth set NCBI keeps the complete stop when it extends 1–2 nt, and mostly uses T/TA at >= 5 nt); a complete stop that overlaps the tRNA by only 1–2 nt is kept. With no downstream tRNA, and the complete stop running > 30 nt into the next same-strand protein gene (> 10 nt in annelids, `--abbreviated_pcg_overlap`, consistent with O1 at the origin; atp8-atp6 and nad4L-nad4 excepted; in the truth set, within 11–30 nt the complete stop is kept more often and deeper than that T/TA is the majority; across all local NCBI records the preference for 11–30 nt differs strongly by clade, annelids, cnidarians and molluscs) -> do not extend: keep the original T/TA end (e.g. cox1 -> cox2), otherwise take the T/TA before the next gene (e.g. nad6 -> cob); if there is no T/TA before the downstream feature, a complete stop extending into it by <= 20 nt (`OVERLAP_FALLBACK`; predicted tRNA ends are imprecise) is adopted; a complete stop that the core alignment already reached and that extends <= 20 nt into a predicted tRNA is kept; a complete stop that the core tool annotated and that extends into the next gene is not pulled back; a gene that already ends in T/TA is not extended into other features (rRNA, opposite-strand tRNA) by > 10 nt for a complete stop within 90 nt; an annotated end that stops before the first complete stop, with no annotated gene in between -> extend to that stop (<= 600 nt; the core tool lost the C-terminal homologous region in divergent clades: nad5/nad4/nad2/cox1 short by 100–500 nt); a gene that would enclose a whole same-strand tRNA with no T/TA before the tRNA -> ends before that tRNA (3′ partial). Vertebrates (`--abbreviated_overlap`, added automatically by `mitom run` for Chordata): a complete stop extending 1–2 nt into a downstream same-strand tRNA is also changed to T/TA, but a complete stop extending into the next protein gene is kept (in NCBI vertebrate records the TAA of atp6 and the ATG of cox3 share an A); other phyla unchanged |
| D1 | — | 同一蛋白基因在不同位置注释两次（MitoZ 遇重复区）→ 只保留完整、较长、靠前的一份；另一条链上与之重叠、不完整或含内部终止的同名拷贝（反义链假阳性）也删除 / the same protein gene annotated twice at different positions (MitoZ on repeated regions) -> keep only the complete, longer, more upstream copy; a same-named copy on the other strand that overlaps it and is incomplete or has an internal stop (an antisense-strand false positive) is also deleted |
| P7 | — | 起点与终点不在同一读码框（内部终止，NCBI ERROR）→ 取两者中较长的开放读码框：保留起点到第一个终止，或保留终止、在其读码框内取离原起点最近的起始密码子；新 ORF 不足原长一半时不改，只报告（多类群 benchmark：内部终止 14 → 5） / start and end not in the same reading frame (internal stop, NCBI ERROR) -> take the longer open reading frame of the two: keep the start up to the first stop, or keep the stop and take, within its reading frame, the start codon nearest to the original start; if the new ORF is shorter than half the original length, do not change it, only report |
| P8 | — | 起点与终点同框但有内部终止（起点越过终止子放在上游，或终点越过基因的终止子）→ 取"原起点到第一个终止"与"最后一个内部终止之后的第一个起始密码子到原终点"中较长者（≥ 原长一半） / start and end in the same frame but with an internal stop (the start was placed upstream past a stop, or the end past the gene's stop) -> take the longer of 'the original start to the first stop' and 'the first start codon after the last internal stop to the original end' (>= half the original length) |
| P9 | — | 一个多余碱基不翻译（鸟类、龟类 ND3 第 174 位，Mindell et al. 1998；少数其他基因）→ 按 NCBI 写成 `join(起点..p-1,p+1..终点)` 并加 NCBI 的 note；只在单读码框修正都保不住 85% 基因、且原起止点 30 nt 内有起始/终止密码子能跨过一个碱基构成开放读码框时使用（2026-10 运行：107 个 ND3 中 97 个两端与 NCBI 完全一致）。**拼接之前**先在其他读码框找一个不间断的 ORF（起始密码子在原起点上游 `--max_start_shift` 个密码子内、不进入上游基因，终止在原终点 30 nt 内），覆盖原基因 ≥ 85% 就用它——那是核心放错了读码框，不是移码（臭海蛹 *Travisia forbesii* ATP8，OZ479805：49 aa 拼接 → 53 aa 单框，与 *T. sanrikuensis* 及另一臭海蛹一致，转录组与 HiFi 都保留被“跳过”的碱基；67 条环节动物基因组项目线粒体中只此 1 个基因改变）；真正多一个碱基时（鸟、龟 ND3）各读码框都只覆盖基因的一部分，不受影响 / one extra base is not translated (birds and turtles, ND3 position 174, Mindell et al. 1998; a few other genes) -> written as NCBI does, `join(start..p-1,p+1..end)`, with NCBI's note; used only when no single-reading-frame correction keeps 85 % of the gene and a start/stop codon within 30 nt of the original start/end can complete an open reading frame across one base. **Before the join**, an uninterrupted ORF in another reading frame is looked for (start codon within `--max_start_shift` codons upstream of the original start, not entering the upstream gene, stop within 30 nt of the original end); if it covers >= 85 % of the original gene it is used - the core had put the gene in the wrong reading frame, it is not a frameshift (example: *Travisia forbesii* ATP8, OZ479805: a 49-aa join -> a 53-aa single frame, agreeing with *T. sanrikuensis* and another *Travisia*; the transcriptome and the HiFi reads both keep the 'skipped' base). With a real extra base (bird, turtle ND3), each reading frame covers only part of the gene, so this step does not apply |
| P10 | — | 以上都修不好、仍有内部终止的 CDS → 基因与 CDS 标 `/pseudo`（不翻译），note 写明，需人工核对（table2asn 否则报 InternalStop） / a CDS that none of the steps above can fix and that still has an internal stop -> gene and CDS marked `/pseudo` (not translated), with a note; needs manual checking (otherwise table2asn reports InternalStop) |
| O2 | — | 核心工具写成跨原点 join、原点后只有一小段（≤ 300 nt）的 CDS（cox1 在首位时环节动物 nad2 伸进 cox1 却没有终止子，NoStop）→ 截回到序列末端，再由 P1/P3/O1 重新定终止；O1 在原点后 60 nt 内找不到终止、而原点处紧接下一个基因时，读码框内落在序列末端的 T/TA 作不完整终止（NCBI 写法） / a CDS that the core tool wrote as a join across the origin with only a short piece (<= 300 nt) after the origin (when cox1 is first, annelid nad2 runs into cox1 without a stop, NoStop) -> cut back to the end of the sequence, then P1/P3/O1 set the stop again; if O1 finds no stop within 60 nt after the origin and the next gene follows immediately at the origin, a T/TA in frame at the end of the sequence is taken as an incomplete stop (the NCBI convention) |
| T1 | — | 完全落在同链蛋白基因内部的弱 tRNA 预测（MiTFi E > 1e-3）：是第二拷贝或位于基因中间三分之一 → 删除（tRNA 是转录本的间断点，同链基因中部不可能有完整 tRNA；table2asn BadCDScomponentOverlapTRNA）；位于基因首/尾三分之一 → 保留 tRNA，由 P1/P3 把基因端点移到它之外（棘头动物 nad6 被延伸到 trnD 内的终止子，NCBI 的 nad6 止于 trnD 之前） / a weak tRNA prediction (MiTFi E > 1e-3) lying entirely inside a same-strand protein gene: if it is a second copy or lies in the middle third of the gene -> delete (a tRNA is a punctuation point of the transcript, so a complete tRNA cannot lie in the middle of a same-strand gene; table2asn BadCDScomponentOverlapTRNA); if it lies in the first or last third of the gene -> keep the tRNA, and P1/P3 move the gene end outside it (acanthocephalan nad6 was extended to a stop inside trnD; NCBI's nad6 ends before trnD) |
| C1 | — | 修正后仍完全包住同链 tRNA 的 CDS：无起始/终止密码子或已标 pseudo → 视为假阳性删除（扁形动物"atp8"、毛颚动物"atp6"：这两个门本无此基因）；完整基因只报告 / a CDS that, after correction, still encloses a same-strand tRNA entirely: with no start/stop codon or already marked pseudo -> deleted as a false positive (flatworm 'atp8', chaetognath 'atp6': these two phyla have no such gene); a complete gene is only reported |
| P4 | — | 内部终止密码子、非起始密码子首位、多段/部分 CDS → 只报告，不改 / internal stop codons, a first codon that is not a start codon, multi-part / partial CDSs -> report only, no change |
| R1 | (vi) rRNA 边界不扩展到两侧基因 / (vi) rRNA boundaries are not extended into the genes on either side | 与两侧 tRNA/CDS/rRNA 的重叠切掉（默认 `--rrna_bounds trim`）；`fill` 再补到两侧基因（空隙 ≤ `--max_fill` 100 nt，GenBank 常见写法） / overlaps with the flanking tRNA/CDS/rRNA are cut off (default `--rrna_bounds trim`); `fill` extends the rRNA to the flanking genes (gap <= `--max_fill` 100 nt, the common GenBank convention) |
| R3 | 命名 / naming | rrnS / rrnL，12S / 16S ribosomal RNA；rRNA 与 CDS 重叠 > 50 nt 报告 / rrnS / rrnL, 12S / 16S ribosomal RNA; an rRNA overlapping a CDS by > 50 nt is reported |

tRNA 边界保持注释软件结果（准则 v）。每处修改写入 `curation_report.tsv`（规则、旧/新坐标、理由）。

tRNA boundaries stay as the annotation software gave them (criterion v). Every change is written to `curation_report.tsv` (rule, old / new coordinates, reason).

**内置 NCBI 评估**（单独运行 `mitom curate` 时；`mitom run` 不在这里校验，而在第 8 步统一校验并生成提交文件）：有 `table2asn` 时（PATH、`--table2asn` 或 mitointron 环境），每条记录在校正前后各用 NCBI 校验一次，`validation_compare.tsv` 按"记录 + 错误代码 + 特征"列出每类问题 fixed / remaining / new，终端汇总一行（如 `ERROR 128 -> 4 … fixed 124, remaining 4, new 0`）；"new" 即校正引入的新问题，应为 0。`--no_validate` 关闭。

**Built-in NCBI evaluation** (when `mitom curate` is run alone; `mitom run` does not validate here but validates once in step 8 and makes the submission files): when `table2asn` is available (in PATH, `--table2asn`, or the mitointron environment), each record is validated by NCBI before and after curation; `validation_compare.tsv` lists, per 'record + error code + feature', each class of problem as fixed / remaining / new, and a one-line summary is printed on the terminal (e.g. `ERROR 128 -> 4 … fixed 124, remaining 4, new 0`); 'new' means a new problem introduced by curation and should be 0. `--no_validate` turns it off.

```bash
mitom curate -i annotated_gb_dir -o curated_dir --genetic_code 5 [--rrna_bounds trim|fill]
```

**验证**（31 条鳞沙蚕类真实 MitoZ 注释，/mnt/n/DLC_mt_2023_1016/Mitoz_annotate，对照 NCBI 原记录；坐标换算到 NCBI 序列）：

**Validation** (31 real MitoZ annotations of scale-worm polychaetes (Aphroditiformia), /mnt/n/DLC_mt_2023_1016/Mitoz_annotate, against the original NCBI records; coordinates converted to the NCBI sequence):

| | MitoZ 原始 / MitoZ raw | 校正后 / curated |
|---|---|---|
| PCG 5′ / 3′ 与 NCBI 相同（284 个） / PCG 5′ / 3′ identical to NCBI (284 genes) | 253 / 211 | **265 / 254** |
| 12S 5′ / 3′（22 个） / 12S 5′ / 3′ (22 genes) | 5 / 2 | **17 / 15** |
| 16S 5′ / 3′（22 个） / 16S 5′ / 3′ (22 genes) | 0 / 0 | **17 / 11** |
| rRNA 与相邻基因重叠 / rRNA overlapping an adjacent gene | 44 | **0** |
| table2asn ERROR | 128（126 个 NoStop） / 128 (of which NoStop) | **4**，全部对应报告中标出需人工检查的 3 个基因（跨原点 cox1、含 24 个内部终止的 COX1、tRNA 前无终止的 ATP6） / **4**, all belonging to the 3 genes that the report flags for manual checking (cox1 across the origin, COX1 with 24 internal stops, ATP6 with no stop before the tRNA) |

剩余差异主要来自 tRNA 边界与 NCBI 不同（NCBI 的 T/TA 终止落在其 tRNA 前），以及 rRNA 写法惯例（NCBI 多为两侧填满；trim 模式按准则不扩展）。准则来自腹足类，此处在环节动物上验证，软体动物数据尚未测试。

The remaining differences come mainly from tRNA boundaries that differ from NCBI's (NCBI's T/TA stop lies before its tRNA) and from rRNA conventions (NCBI mostly fills up to both neighbours; the trim mode does not extend, following the criterion). The criteria come from gastropods and are validated here on annelids; molluscan data have not been tested.

## 漏注的蛋白基因：近缘谱 HMM / Missing protein genes by profile HMM (mito_pcg_hmm.py)

注释核心（blastx + exonerate，参考 refseq89m）对高度分化的短基因（atp8、atp6、nad4l、nad6）缺近缘参考时会漏注，如龙介虫科
的 atp8。第 6 步补漏中（MITOS2 补漏之后，`--pcg_hmm yes` 默认），**只对仍缺的蛋白基因**：
1. 从 `models/pcg_refs.faa.gz`（`scripts/build_pcg_refs.py`：RefSeq 与 NCBI 环节动物的 13 个基因蛋白，带分类谱系）选最近缘
   的参考——谱系共享层级最长者优先，层级相同按记录自身 cox1（无则 cob ……）的 3-mer 相似度；
2. 最近 10、30、100 条各建一个谱 HMM（MAFFT、hmmbuild），在没有外显子覆盖的 ORF 中 `hmmsearch --max`，E 值乘 3；
3. 接受 E ≤ 0.01 × min(1, p/(1−p))，p 为最近 30 条已注释近缘记录中有该基因的比例（近缘普遍丢失的基因需要更强证据），
   模型覆盖 ≥ 40 %；从最靠近比对起点的起始密码子到终止密码子（与下一基因重叠 ≤ 30 nt），起点再由校正调整。
每次搜索在 `pcg_hmm_report.tsv` 一行（找到与否、E 值、所用参考与先验），找到的在日志中给 WARNING 提示核查。
测试：`benchmark/results/pcg_hmm_2026-11-10/README.md`（留一属找回率、扁虫 / 线虫的误报、335 条环节动物）。

The core annotation (blastx + exonerate, reference refseq89m) misses highly divergent short genes (atp8, atp6, nad4l, nad6) when no close reference exists, e.g. atp8 in Serpulidae. In step 6, after the MITOS2 gap-filling (`--pcg_hmm yes`, default), **only for protein genes that are still missing**:

1. Choose the closest references from `models/pcg_refs.faa.gz` (`scripts/build_pcg_refs.py`: the 13 gene proteins of RefSeq and NCBI annelids, with taxonomic lineage): the reference sharing the longest lineage wins; ties are broken by the 3-mer similarity of the record's own cox1 (cob if absent, ...).
2. Build one profile HMM (MAFFT, hmmbuild) each from the nearest 10, 30 and 100 references; run `hmmsearch --max` on the ORFs with no exon coverage; E-values are multiplied by 3.
3. Accept E <= 0.01 x min(1, p/(1-p)), where p is the fraction of the 30 nearest annotated relatives that have the gene (a gene commonly lost in relatives needs stronger evidence), and model coverage >= 40 %; the gene runs from the start codon nearest to the alignment start to the stop codon (overlap with the next gene <= 30 nt), and the start is then adjusted by curation.

Each search gets one row in `pcg_hmm_report.tsv` (found or not, E-value, the references used and the prior); a gene that is found gives a WARNING in the log for checking. Tests: `benchmark/results/pcg_hmm_2026-11-10/README.md` (leave-one-genus-out recall, false positives in flatworms / nematodes, 335 annelid records).

## 转录组证据 / Transcript evidence (mito_rna.py, `mitom rna`)

RNA-seq reads 只比对到该线粒体基因组（minimap2 `-ax sr`，否则 bwa mem；MAPQ ≥ 20）。动物线粒体多顺反子转录、在 tRNA 处切开，
覆盖连续，边界不能由覆盖定；给出的证据是：每个基因 / ORF / NCR 的深度相对 13 个蛋白基因中位数之比，深度 ≥ 中位数 10 % 的位置
比例（部分转录），“像基因一样表达”（不低于表达最低蛋白基因的一半；中位深度 < 10 时不判）；文库方向（含反向链特异自动翻转）；
polyA 3′ 端（SAM 中按参考正向书写：比对后截去的 A 串为正链转录本 3′ 端，比对前截去的 T 串为负链；tRNA 3′ 端 CCA 尾不计）。
验证：*Polydora* 公开转录组复现 Selifanova 等 2023 对三个 ORFan 与 ND3 下游区段的结论（`benchmark/results/rnaseq_2026-11-10/`）。

RNA-seq reads are mapped only to this mitogenome (minimap2 `-ax sr`, otherwise bwa mem; MAPQ >= 20). Animal mitochondrial genomes are transcribed polycistronically and cut at the tRNAs, so coverage is continuous and boundaries cannot be set by coverage. The evidence given is: for each gene / ORF / NCR, its depth relative to the median of the 13 protein genes; the share of positions with depth >= 10 % of the median (partial transcription); an "expressed like a gene" call (not lower than half of the lowest-expressed protein gene; not called when the median depth is < 10); library strandedness (reverse-strand-specific libraries are flipped automatically); and polyA 3′ ends (in the SAM, written in the reference forward orientation: an A stretch trimmed after alignment marks the 3′ end of a plus-strand transcript, a T stretch trimmed before alignment marks the minus strand; the CCA tail at the 3′ end of tRNAs is not counted). Validation: a public *Polydora* transcriptome reproduces the conclusions of Selifanova et al. 2023 on the three ORFans and on the segment downstream of ND3 (`benchmark/results/rnaseq_2026-11-10/`).

## 多分体线粒体基因组 / Parts of one genome (mito_multipart.py)

有的线粒体基因组由几个分子组成（海绵 *Lycopodina hypogea* 3 条染色体，nad4l 单独在 291 bp 的小环上；线虫 *Allodiplogaster
pararmata* 一个 HiFi 组装两个环，基因互补）。逐条看会“缺”基因。第 9 步把同一组装（COMMENT “The assembly …”）或同一
BioSample（DBLINK）的记录放在一起：序列相同的只算一次（21-mer ≥ 99 % 相同）；各自缺蛋白基因、合起来 ≥ 12 个的报告为
“parts of one genome”（多分体，或一个基因组的几段），`multipart_report.tsv` 一组一行，日志给 WARNING；长度不足常见长度一半的
基因记为片段。测试：`benchmark/results/ebp_mito_2026-11-10/README.md`。

Some mitochondrial genomes consist of several molecules (the sponge *Lycopodina hypogea* has 3 chromosomes, with nad4l alone on a 291-bp minicircle; the nematode *Allodiplogaster pararmata* has two circles from one HiFi assembly, with complementary genes). Looked at one record at a time, genes seem to be missing. Step 9 groups the records of the same assembly (COMMENT "The assembly …") or the same BioSample (DBLINK): identical sequences are counted once (21-mer >= 99 % identical); records that each lack protein genes but together have >= 12 are reported as "parts of one genome" (multipartite, or several pieces of one genome), one row per group in `multipart_report.tsv`, with a WARNING in the log; a gene shorter than half the usual length is counted as a fragment. Test: `benchmark/results/ebp_mito_2026-11-10/README.md`.

## 非编码区 / Non-coding regions (mito_ncr.py)

每条记录中 ≥ `--ncr_min_len`（默认 100）nt、没有任何特征覆盖的区段写成 misc_feature `/note="non-coding region NCR<k>, <长度> bp: no gene annotated (MitoM)"`：按在记录中的顺序编号，跨原点的写成 join；带内含子的基因覆盖它的内含子；再次运行时替换旧的 NCR。这一步在校正之后（`run` / `annotate` / `import`；`--ncr no` 不写，圈图也不画）。未命名的控制区、可能含未识别基因的序列因此在 GenBank、NCBI 提交文件和圈图里都看得到。

Every stretch of >= `--ncr_min_len` (default 100) nt in a record that is not covered by any feature is written as a misc_feature `/note="non-coding region NCR<k>, <length> bp: no gene annotated (MitoM)"`: numbered in the order of the record; one across the origin is written as a join; a gene with an intron covers its intron; on a re-run the old NCRs are replaced. This step comes after curation (`run` / `annotate` / `import`; `--ncr no` does not write them and the circular plot does not draw them). An unnamed control region, or a sequence that may hold an unrecognised gene, is therefore visible in the GenBank file, the NCBI submission files and the plot.

验证（2026-10-30）：22 条 MJ20260820041 记录写出 26 个 NCR，table2asn 0 ERROR，WARNING 与不写 NCR 时相同（18 个）；20 条 NCBI 的 Hesionidae 线粒体基因组 33 个 NCR（24,280 nt），在 NCBI 有注释的记录上，与 NCBI 未注释的区段（如 *Sirsoe methanicola* 2,641 nt，其中一段含 288 个密码子的 ORF）及其控制区特征一致。

Validation (2026-10-30): on the MJ20260820041 records, table2asn gave no ERROR and the same WARNINGs as without NCRs; on the NCBI Hesionidae mitogenomes, the NCRs written agree with the stretches NCBI left unannotated (e.g. in *Sirsoe methanicola*, one of which holds an ORF) and with their control-region features, in the records that have NCBI annotation.

## 未知功能的 ORF / ORFs of unknown function (mito_orf.py)

13 个蛋白基因以外的编码区（海绵 atp9、八放珊瑚 mutS、淡水蚌 F-ORF / M-ORF、内含子中的归巢内切酶，才女虫科 Polydorini
在 trnQ 与 trnA 之间约 355 个密码子的 ORF 等）：校正之后、标 NCR 之前，在没有外显子覆盖的序列（内含子也算）里找 ≥ 40 个
密码子的 ORF（`--orf yes` 默认）。只看长度不行（ORFfinder / Geneious 式：≥ 75 个密码子在 362 条无内含子 RefSeq 中报出
357 个），判定依据是“受选择而保持开放”：
1. **长度超出随机**：按本记录未注释序列的碱基组成，随机出现至少这么长 ORF 的期望数 `E_len = 2·F·p_start·(1−p_stop)^L`；
2. **保持开放**：同源（blastp E ≤ 1e-5）本身也是 ORF——同一次运行的其他记录、ORF 库 `models/orf_library.fna.gz`
   （16,073 条 RefSeq 动物线粒体基因组，`scripts/build_orf_library.py`）、`--orf_ref` 近缘种。最多 20 个同源放到 ORF 上，
   替换数取两个严格下限中的较大者：每个位置多出的碱基种类数之和；两两差异最小生成树长度的一半（任何带祖先的历史至少这么长）。
   非编码序列中每次替换产生终止的概率为 f（由 ORF 自身的密码子算出），全部碰巧保持开放的概率 `p_open = (1−f)^替换数`；
   近缘种单独不构成证据，多个近缘种可以累加；
3. **编码读框**：在与同源的核苷酸比对中（不受读框影响；蛋白比对会偏向第三位差异），编码读框第三位差异最多、第二位最少
   （第二位的改变一定改变氨基酸）：对 1/3 均匀的卡方检验，只在这个顺序成立时计；与真基因第三位重合的反义读框，第一位
   差异最少，因而不成立；
4. 先屏蔽串联重复（周期 3–150 nt、相邻拷贝 ≥ 85 % 一致）：无终止的重复单元能拼出任意长的 ORF，也能与别的重复比对上；
5. 重叠的 ORF（其他读框、反义链）只保留第二位差异最少、第三位比例最高、证据最强的一个。

`p_open` 与第 3 条的 p 都 ≤ 0.05（两者都是必要条件；非编码序列同时碰巧通过约 0.0025）的（`conserved`），或属已知蛋白家族的
（`known_family`：MFannot 家族，atp9 用其基因名和产物；Pfam 家族，`models/orf_families.hmm.gz`，`hmmsearch --cut_ga` 用 Pfam
自己审定的阈值，命名 LAGLIDADG / GIY-YIG / HNH 归巢内切酶、II 类内含子逆转录酶 / 成熟酶、B 族 DNA 聚合酶、MutS），写成 CDS
`/gene="orf<密码子数>"`，`/product` 为 “hypothetical protein” 或家族产物，note 写明证据。读框读进下一个基因前没有终止的：
读进 **tRNA** 时按线粒体 tRNA 间断加工止于 tRNA 之前——末端是 T / TA 时为不完整终止，否则止于最后一个完整密码子、3′ 端标为
不完整（`>`；table2asn 只给 WARNING）。这也是 Selifanova 等 2023（*Genome Biol Evol*，PMC10721130）定 *Polydora* ORFan 的
方式：无终止区段在 tRNA 边界切开。读进蛋白基因或 rRNA 的不写，给 WARNING（同读框时写明“可能是该基因起点更靠上游”）。只满足 1 的（`long_unconfirmed`）、与 13 个基因之一同源的（`standard_fragment`）也不写，
在日志中给 WARNING、在 `orf_report.tsv` 列出、写进所在 NCR 的 note，不静默。`--orf all` 连 `long_unconfirmed` 也写；
`report` 只报告；`no` 跳过。测过但未采用：密码子使用偏好（GeneMark 式，与本记录 13 个基因比较）区分力太弱；RNAcode
（Washietl 等 2011，比较基因组编码潜力，需把线粒体密码子换成标准同义密码子）与上述两项检验相比几乎没有增益。
测试结果（含与 MITOS2、ORFfinder 式长度判定的比较）：`benchmark/results/orf_2026-11-09/README.md`。

Coding regions other than the 13 protein genes (sponge atp9, octocoral mutS, freshwater mussel F-ORF / M-ORF, homing endonucleases in introns, the ORF of about 355 codons between trnQ and trnA in Polydorini (Spionidae), etc.): after curation and before NCRs are marked, ORFs of >= 40 codons are searched for in sequence not covered by exons (introns included) (`--orf yes`, default). Length alone does not discriminate (an ORFfinder / Geneious-style rule of >= 75 codons reports ORFs in almost all intron-free RefSeq genomes); the criterion is "kept open by selection":

1. **Length beyond chance**: from the base composition of this record's unannotated sequence, the expected number of ORFs at least this long arising at random, `E_len = 2·F·p_start·(1−p_stop)^L`.
2. **Kept open**: a homolog (blastp E <= 1e-5) is itself an ORF - in other records of the same run, in the ORF library `models/orf_library.fna.gz` (16,073 RefSeq animal mitogenomes, `scripts/build_orf_library.py`), or in close relatives given by `--orf_ref`. Up to 20 homologs are placed on the ORF; the number of substitutions is the larger of two strict lower bounds: the sum over positions of the number of extra base kinds, and half the length of the minimum spanning tree of the pairwise differences (any history with an ancestor is at least this long). In non-coding sequence each substitution creates a stop with probability f (computed from the ORF's own codons); the probability that all of them by chance keep the ORF open is `p_open = (1−f)^substitutions`. A close relative alone is not evidence; several close relatives add up.
3. **Coding frame**: in the nucleotide alignment with the homologs (independent of reading frame; a protein alignment is biased toward third-position differences), the coding frame has most differences at position 3 and fewest at position 2 (a change at position 2 always changes the amino acid): a chi-square test against 1/3 uniform, counted only if this order holds; an antisense frame overlapping the third position of a real gene has fewest differences at position 1, so it fails.
4. Tandem repeats are masked first (period 3–150 nt, adjacent copies >= 85 % identical): a repeat unit without a stop can build an ORF of any length, and can also be aligned to other repeats.
5. Overlapping ORFs (other frames, antisense strand) keep only the one with the fewest position-2 differences, the highest position-3 share and the strongest evidence.

An ORF whose `p_open` and the p of item 3 are both <= 0.05 (both are necessary; non-coding sequence passes both by chance with about 0.0025) (`conserved`), or that belongs to a known protein family (`known_family`: MFannot families, with atp9 by its gene name and product; Pfam families, `models/orf_families.hmm.gz`, `hmmsearch --cut_ga` with Pfam's own curated thresholds, naming LAGLIDADG / GIY-YIG / HNH homing endonucleases, group II intron reverse transcriptase / maturase, family B DNA polymerase, MutS), is written as a CDS with `/gene="orf<number of codons>"` and `/product` "hypothetical protein" or the family product, with the evidence in the note. A frame that reads on into the next gene without a stop: into a **tRNA**, it stops before the tRNA by mitochondrial tRNA punctuation - an incomplete stop if the end is T / TA, otherwise it ends at the last complete codon with the 3′ end marked incomplete (`>`; table2asn gives only a WARNING). This is also how Selifanova et al. 2023 (*Genome Biol Evol*, PMC10721130) defined the *Polydora* ORFans: stop-free stretches are cut at tRNA borders. A frame reading into a protein gene or an rRNA is not written and gives a WARNING (for the same frame, saying that it "may be an upstream start of that gene"). Those satisfying only item 1 (`long_unconfirmed`) and those homologous to one of the 13 genes (`standard_fragment`) are not written either; they get a WARNING in the log, a row in `orf_report.tsv` and a mention in the note of the NCR they lie in, never silently. `--orf all` also writes `long_unconfirmed`; `report` only reports; `no` skips. Tried but not adopted: codon usage bias (GeneMark style, compared with this record's 13 genes) discriminates too weakly; RNAcode (Washietl et al. 2011, comparative coding potential, which needs mitochondrial codons converted to standard synonymous codons) adds almost nothing over the two tests above. Test results (including the comparison with MITOS2 and an ORFfinder-style length rule): `benchmark/results/orf_2026-11-09/README.md`.

## 环状还是线状 / Topology (mitom_topology.py, mito_core.py --check_ends)

**问题**：标签不一定对。本地 NCBI 记录中，669 条环节动物“complete genome”有 64 条标为 linear，1,050 条 RefSeq 完整基因组有 19 条标为 linear，
11 条水母类（线状线粒体基因组）标为 circular；自己组装的序列多条放在一个 FASTA 里又没写 `topology=` 时全按线状处理，起点落在基因中间时这个基因会被拆成两段。

**Problem**: the labels are not always right. In the local NCBI records, a number of annelid "complete genome" records are labelled linear, a number of RefSeq complete genomes are labelled linear, and some medusozoan records (whose mitogenomes are linear) are labelled circular; when several self-assembled sequences are in one FASTA without `topology=`, all are treated as linear, and a gene is split in two when the start falls inside it.

**先试过“全部接起来”（2026-10-30，93 条 NCBI 记录，每条按环状、按线状各注释一次，与 NCBI 比较）**：两种方式蛋白基因两端完全一致的比例相差 ≤ 2 个百分点，
NCBI ERROR 都是 0；结果不同的 20 条中，12 条是 nad2 紧接 cox1 时两条终止规则（序列中间的 P3、原点处的 O1）不一致，与拼接本身无关（已统一，见上表 P3）；
真正与两端有关的 8 条里，拼接更准 3 条、更差 5 条（部分基因组末端被截断的基因被错误地“接”到开头）；全部拼接总用时多 22 %（多搜 3,000 nt、定向后坐标偏移的记录再注释一轮）。
所以不全部拼接，只核对两端。

**First tried "join everything" (2026-10-30, 93 NCBI records, each annotated once as circular and once as linear and compared with NCBI)**: the share of protein genes with both ends identical hardly differed between the two ways, and NCBI ERROR was 0 for both; of the records with different results, most were due to the inconsistency of two stop rules (P3 inside the sequence, O1 at the origin) when nad2 directly follows cox1, which has nothing to do with joining (since unified, see P3 above); of those truly related to the ends, joining was more accurate in a few and worse in more (a gene truncated at the end of a partial genome was wrongly "joined" to the start); joining everything also took longer (3,000 nt more to search, and records whose coordinates shifted after reorientation are annotated one more round). So not everything is joined; only the ends are checked.

**两端核对**（`--topology auto`，默认）：先按标签注释（GenBank LOCUS 行；FASTA `topology=`；没有标签时：两端重叠的组装结果去掉重复后为环状，表头写 circular / complete genome 为环状、partial genome 为线状，其余 ≥ 10 kb 的序列（一个文件一条或多条都一样）先按环状注释，再由两端和基因组成决定——没有基因跨两端也没有基因被截断时，蛋白基因 ≥ 12 个且两个 rRNA 齐全为环状，否则为线状；< 10 kb 为线状）。
注释核心对两端的每个蛋白基因（跨原点的，或离序列端点比它缺的残基所需的长度还近的）把序列首尾相接，用该基因的参考蛋白（注释时用的那条）跨接点再比对一次（exonerate）：

**End check** (`--topology auto`, default): the record is first annotated according to its label (GenBank LOCUS line; FASTA `topology=`; with no label: an assembly whose ends overlap, after removing the repeat, is circular; a header saying circular / complete genome is circular, partial genome is linear; other sequences >= 10 kb (one or many per file, no difference) are first annotated as circular and then decided by the ends and gene content - if no gene spans the ends and no gene is truncated, >= 12 protein genes and both rRNAs complete means circular, otherwise linear; < 10 kb is linear).

For every protein gene at the ends of the sequence (one that spans the origin, or lies closer to a sequence end than the length needed for the residues it lacks), the annotation core joins the sequence end to start and aligns the gene's reference protein (the one used in annotation) across the junction once more (exonerate):

| 判定 / call | 条件 / condition | 含义 / meaning |
|---|---|---|
| joined | 参考蛋白跨接点比对，两侧各 ≥ 30 nt；并且两侧各自单独也能比对上、到达接点、残基首尾相接（相差 ≤ 10），或两侧各 ≥ 90 nt、读码框不断、两侧单独的比对不互相矛盾 / the reference protein aligns across the junction with >= 30 nt on each side; and either each side also aligns alone, reaches the junction, and the residues join end to end (differing by <= 10), or each side has >= 90 nt, the reading frame is unbroken, and the separate alignments of the two sides do not contradict each other | 两端是连着的：环状 / the ends are joined: circular |
| cut | 比对停在接点 60 nt 内，那一侧缺 ≥ 30 个残基，读码框一直开放到接点 / the alignment stops within 60 nt of the junction, that side lacks >= 30 residues, and the reading frame stays open up to the junction | 基因被端点截断：部分序列，线状 / the gene is truncated by the sequence end: partial sequence, linear |
| joined_intron | 看似 cut，但用 exonerate protein2genome（允许内含子）在接点前、后剩下的整圈序列上各再比对一次，缺的残基在另一端比对上（中间隔着跨接点的内含子，长度不限） / looks like cut, but when aligned again with exonerate protein2genome (introns allowed) on the remaining whole-circle sequence before and after the junction, the missing residues align at the other end (separated by an intron across the junction, of any length) | 两端是连着的：环状；内含子照常注释 / the ends are joined: circular; the intron is annotated as usual |
| complete | 其他：完整基因、N/C 端分化未比对上，或跨接点的比对只是偶然延伸（两侧单独比对接不上又不够长） / anything else: a complete gene, N/C-terminal divergence that did not align, or an alignment across the junction that is only a chance extension (the two sides do not join when aligned alone and are not long enough) | 不作判断，按标签 / no call; follow the label |

有 joined → 环状；标为环状而有 cut → 线状；标为环状的水母类（Hydrozoa、Scyphozoa、Cubozoa、Staurozoa）和钙质海绵（线状线粒体基因组）没有 joined → 线状，
且从第一轮就按线状注释（不必注释两次）；标为 linear 的 GenBank 记录，DEFINITION 写明是完整基因组（`complete genome`，或 DToL / ATLASea / MitoHiFi 细胞器组装的 `genome assembly, organelle: mitochondrion`），基因齐全（≥ 12 个蛋白基因和两个 rRNA）、两端没有 cut、不属于线状类群 → 环状（2026-11-09，issue #8：OZ287881、OZ346613 是完整环状分子，约 6 kb 非编码区被拆在两端）；其余按标签。与第一轮拓扑不同的记录用新拓扑再注释一遍；每条记录写入 `qc/topology_decision.tsv`。
MitoZ / MITOS2 / `mitom import` 的注释没有核心的比对证据，按标签。

Any joined -> circular; labelled circular but with a cut -> linear; labelled-circular medusozoans (Hydrozoa, Scyphozoa, Cubozoa, Staurozoa) and calcareous sponges (linear mitogenomes) with no joined -> linear, and annotated as linear from the first round (no need to annotate twice); a GenBank record labelled linear whose DEFINITION states a complete genome (`complete genome`, or `genome assembly, organelle: mitochondrion` for DToL / ATLASea / MitoHiFi organelle assemblies), with a full gene set (>= 12 protein genes and two rRNAs), no cut at the ends, and not in a linear-genome group -> circular (2026-11-09, issue #8: OZ287881 and OZ346613 are complete circular molecules whose non-coding region of about 6 kb was split across the two ends); otherwise follow the label. A record whose topology differs from that of the first round is annotated again with the new topology; each record is written to `qc/topology_decision.tsv`.

Annotations from MitoZ / MITOS2 / `mitom import` have no alignment evidence from the core, so they follow the label.

**校准**（上述 93 条记录各按环状、线状注释一次；28 条记录两端有蛋白基因，共 72 个基因次）：真实环状记录跨两端的基因（NC_035740.1 cox1，参考蛋白在接点前 1,548 nt、后 42 nt）两次都是 joined；
部分基因组末端被截断的 nad2（PQ221480.1，残基 3–114/331）两次都是 cut；其余 68 个（位于两端的完整基因、N/C 端分化未比对上的基因）全部 complete，无误判。
在这 93 条按自己的标签运行（同一版代码，`--topology auto` 与 `label` 同时运行）：82 条按标签，10 条标为 circular 的水母类改为线状，NC_035740.1 因 cox1 跨两端保持环状；
蛋白基因两端一致 82.2 % vs 82.2 %，NCBI ERROR 0 vs 0，总用时 164.2 vs 164.6 分钟（每条中位 95 vs 93 秒）。

**Calibration** (the 93 records above, each annotated as circular and as linear; protein genes at the ends): the gene spanning the ends of a truly circular record (NC_035740.1 cox1; the reference protein aligns 1,548 nt before and 42 nt after the junction) was joined in both runs; nad2 truncated at the end of a partial genome (PQ221480.1, residues 3–114 of 331) was cut in both runs; all the others (complete genes at the ends, genes with N/C-terminal divergence that did not align) were complete, with no wrong call. The same 93 records were then run with their own labels (the same code, `--topology auto` and `label` run side by side): most followed the label, the medusozoans labelled circular became linear, and NC_035740.1 stayed circular because its cox1 spans the ends; protein-gene agreement of both ends, NCBI ERROR and total time were the same for the two.

**合成检验**（20 条标为 circular 的完整基因组，各门类，2026-10-30）：

**Synthetic tests** (20 complete genomes labelled circular, from various phyla, 2026-10-30):

| 检验 / test | `--topology auto` | `--topology label` |
|---|---|---|
| 从一个蛋白基因中间切开、标为 linear（组装从基因中间开始） / cut in the middle of a protein gene and labelled linear (the assembly starts in the middle of a gene) | 20/20 判为环状 / called circular | — |
| ↳ 蛋白基因两端与 NCBI 完全一致 / protein genes with both ends identical to NCBI | 78.9 %（与正确标签的原记录相同） / (same as the original record with the correct label) | 72.9 % |
| ↳ 漏注蛋白基因 / 有 NCBI ERROR 的记录 / protein genes missed / records with NCBI ERROR | 2 / 0 | 12 / 2 |
| 去掉 1,500 nt（一个蛋白基因的后半段及其后）、标为 circular（部分序列） / 1,500 nt removed (the second half of a protein gene and what follows) and labelled circular (a partial sequence) | 16/20 判为线状，4 条证据不足按标签 / called linear; 4 with insufficient evidence follow the label | — |
| ↳ 被缺口截断的基因：写成不完整 / 被错误地跨两端接起来 / genes truncated by the gap: written as incomplete / wrongly joined across the ends | 16 / 0 | 6 / 2 |
| ↳ 缺口以外的蛋白基因两端一致 / protein genes outside the gap with both ends identical | 160 | 159 |

**跨接点的内含子**（2026-10-30）：组装从带内含子基因的第 2 个外显子开始时（环节动物 cox1 的 II 类内含子：MJ20260820041 的 CHY / *Travisia*、CWSC / Nephtyidae），protein2dna 比对不允许内含子，外显子 2 停在接点，原先判为 cut、按线状注释，内含子从未被注释（cox1 只有外显子 2：1,287 或 848 nt）。现在看似 cut 的基因再用 protein2genome 比对（判定 joined_intron → 环状）；从头识别的内含子在跨原点时把描述留在 CDS 上，记录旋转到该基因第一个外显子后恢复内含子特征（`mitom_orient.add_exon_gap_introns`）。22 条 MJ 线粒体基因组：13 条 CHY + 5 条 CWSC 得到 cox1 内含子（用剪接后的蛋白、测序深度、跨接点 reads 以及 LC677172 / *Nephtys* 的 cox1 内含子核实），其余 4 条不变，table2asn 0 ERROR。整段位于接点一侧的 ≥ 6 kb 内含子（六放珊瑚 nad5 内含子 10–12 kb，内含其他基因）原先漏掉，重新比对改为在基因前、后剩下的整圈序列上进行（不限内含子长度）。内含子基准（193 条 NCBI 记录、309 个内含子，刺胞、环节、扁盘、海绵动物；另 76 条无内含子；去掉本属的位点）：把记录旋转到原点落在外显子中、内含子中间或距内含子两端 15 nt 处，得到的内含子与未旋转时相同（300 / 300 行）；修正前原点距内含子右 / 左端 15 nt 时只找到 295 个中的 232 / 190 个。

**Introns across the junction** (2026-10-30): when the assembly starts at the second exon of an intron-containing gene (the group II intron of annelid cox1: CHY / *Travisia* and CWSC / Nephtyidae of MJ20260820041), the protein2dna alignment, which does not allow introns, stops at the junction with exon 2, so the call used to be cut and the record was annotated as linear, and the intron was never annotated (cox1 had only exon 2). Now a gene that looks like cut is aligned again with protein2genome (call joined_intron -> circular); for a de novo intron across the origin, the description is left on the CDS, and once the record is rotated to the first exon of that gene the intron feature is restored (`mitom_orient.add_exon_gap_introns`). The CHY and CWSC records got a cox1 intron (checked with the spliced protein, sequencing depth, junction-spanning reads, and the cox1 intron of LC677172 / *Nephtys*), the other MJ mitogenomes were unchanged, and table2asn gave 0 ERROR. An intron of >= 6 kb lying entirely on one side of the junction (the 10–12 kb nad5 intron of hexacorals, containing other genes) used to be missed; the re-alignment is now done on the whole remaining circle before and after the gene (no limit on intron length). Intron benchmark (193 NCBI records, 309 introns, cnidarians, annelids, placozoans, sponges; another 76 records without introns; same-genus sites removed): rotating the record so that the origin falls in an exon, in the middle of an intron, or 15 nt from either end of an intron gives the same introns as the unrotated record; before the fix, an origin 15 nt from the right / left end of an intron found only part of the introns.

早期版本的判定（只要求接点两侧各比对上 30 nt）在部分序列上有 5 条误判为 joined（比对偶然越过假接点 35–338 nt），加上“两侧单独比对首尾相接，或两侧各 ≥ 90 nt 且读码框不断”后为 0；
比对在离端点 20–21 nt 处被 exonerate 截短的基因，按“缺的残基放不下”判断是否位于端点，因此仍被检查。

In an earlier version of the call (requiring only 30 nt aligned on each side of the junction), some partial sequences were wrongly called joined (the alignment by chance crossed a false junction); adding "the two sides align alone and join end to end, or each side has >= 90 nt with an unbroken reading frame" removed these. A gene whose alignment was cut short by exonerate 20–21 nt from the end is still checked, because whether a gene lies at an end is judged by "the missing residues do not fit".

## tRNA 多方法注释 / tRNA (mito_trna.py)

仿 GeSeq（Tillich et al. 2017, NAR 45:W6, doi:10.1093/nar/gkx391）：多个 tRNA 软件结果**全部**写进一个 GenBank（`review/<样本>.gbf`，Geneious 里逐个比较、自己选），另出 `selected/<样本>.gbf`（每位点一个 tRNA）和 `trna_report.tsv`。不预设 tRNA 数目（多数珊瑚 1–2 个，部分环节动物 23–24 个）。

Modelled on GeSeq (Tillich et al. 2017, NAR 45:W6, doi:10.1093/nar/gkx391): the results of several tRNA programs are **all** written into one GenBank file (`review/<sample>.gbf`, to compare one by one in Geneious and choose by hand), and `selected/<sample>.gbf` (one tRNA per locus) and `trna_report.tsv` are also written. The number of tRNAs is not preset (most corals 1–2, some annelids 23–24).

- 软件：MITOS2（只取 tRNA，refseq89m 内置后生动物 tRNA 协方差模型）、ARWEN、ARAGORN、tRNAscan-SE 2（-O）、输入 GB 已有 tRNA。

  Software: MITOS2 (tRNA only; the built-in covariance models of post-metazoan tRNAs in refseq89m), ARWEN, ARAGORN, tRNAscan-SE 2 (-O), and the tRNAs already in the input GB.
- 证据按**独立方法家族**计：ARWEN 与 ARAGORN 同作者同算法，算一个。

  Evidence is counted by **independent method family**: ARWEN and ARAGORN have the same authors and algorithm and count as one.
- 默认选择：MITOS2 找到的采用；MITOS2 没找到的需 ≥2 个独立家族支持。同位点不同反密码子（冲突）保留证据强者，平局都不选。边界取 MITOS2（与 NCBI 完全一致最多）。

  Default selection: a tRNA found by MITOS2 is taken; one that MITOS2 did not find needs support from >= 2 independent families. For different anticodons at the same locus (a conflict) the one with stronger evidence is kept; on a tie neither is chosen. Boundaries are taken from MITOS2 (the most exact matches with NCBI).
- MITOS2 慢的原因：runmitos 对没找到的 tRNA 模型再做灵敏搜索（cmsearch --max、E 0.1、单线程），tRNA 少的类群每条 10 分钟以上并带来假 tRNA。默认 fast 模式（`scripts/runmitos_fast.py`）跳过这一步（一条海葵 640 s → 8 s）；`--mitos_mode full` 恢复。

  Why MITOS2 is slow: runmitos searches again, sensitively (cmsearch --max, E 0.1, single-threaded), for the tRNA models it did not find; in clades with few tRNAs this takes more than 10 minutes per record and brings false tRNAs. The default fast mode (`scripts/runmitos_fast.py`) skips this step (a large speed-up, e.g. on a sea anemone); `--mitos_mode full` restores it.
- 跨原点 tRNA 写成 join。

  A tRNA across the origin is written as a join.

**用 MITOS2 tRNA 替换现有 tRNA**（`--select mitos2`，可选，不在默认流程中）：

**Replacing the existing tRNAs with MITOS2 tRNAs** (`--select mitos2`, optional, not in the default workflow):

```bash
# 默认流程：mitom run（含端点校正 mito_curate.py 与 NCBI 校验）
# 需要 MITOS2 tRNA 时在 1 与 3 之间插入：
mitom trna -i mitom_out/work/04.final_gb -o trna_mitos2 --select mitos2 --genetic_code 5 \
  --mitos_bin MITOS_ENV/bin --refdir DIR_WITH_refseq89m
mitom curate -i trna_mitos2/selected -o curated --genetic_code 5
```

- 输入必须是**调整方向后的最终 GenBank**（`04.final_gb/`）：MITOS2 在该 GB 自身序列上运行，坐标与其它特征一致；FASTA、`00.*`、`01.anno_round1/`（方向调整前）一律拒绝。

  The input must be the **final GenBank after reorientation** (`04.final_gb/`): MITOS2 runs on that GB's own sequence so that the coordinates agree with the other features; FASTA, `00.*` and `01.anno_round1/` (before reorientation) are all rejected.
- 只替换 tRNA 及其 gene 特征（命名同核心，如 `trnN(guu)`），PCG/rRNA 不动；review GB 中保留原 tRNA 供对照。

  Only tRNAs and their gene features are replaced (named as in the core, e.g. `trnN(guu)`); PCGs and rRNAs are untouched; the original tRNAs are kept in the review GB for comparison.
- 某条记录 MITOS2 运行失败 → 该记录原样保留，不会删掉原有 tRNA。

  If MITOS2 fails on a record, that record is left as it was; the existing tRNAs are not deleted.

```bash
mitom trna -i gb_or_fasta -o out --genetic_code 5 --tools mitos2,arwen,aragorn,trnascan \
  --bin ENV/bin --mitos_bin MITOS_ENV/bin --refdir DIR_WITH_refseq89m
```

**验证**（裸序列，对照 NCBI tRNA，氨基酸 + 反密码子 + 重叠 ≥50% 算找到；FP = 该处无 NCBI tRNA；`benchmark/eval_trna.py`）：

**Validation** (bare sequences against NCBI tRNAs; a tRNA counts as found when amino acid + anticodon + overlap >= 50 %; FP = no NCBI tRNA at that place; `benchmark/eval_trna.py`):

| | 刺胞/海绵/扁盘 85 条，183 tRNA / 85 Cnidaria / Porifera / Placozoa records, 183 tRNAs | 环节动物 RefSeq 40 条，880 tRNA / 40 annelid RefSeq records, 880 tRNAs |
|---|---|---|
| MITOS2 | 96.7%，FP 21 | 97.3%，FP 3 |
| ARWEN / ARAGORN | 91.8%，FP 348 / 462 | 95.7% / 96.5%，FP 79 / 117 |
| tRNAscan-SE | 61.2%，FP 1 | 46.1%，FP 0 |
| **默认选择** / **default selection** | **96.7%，FP 21** | **97.5%，FP 3** |

记录清单 `benchmark/trna_sets.tsv`（2026-10-20 重跑）。旧版把 MITOS2 的 `trnM_0(cat)`（同一模型多拷贝）和 `trnQ(---)`（反密码子读不出）丢掉，刺胞 MITOS2 只有 89.6%、默认 93.4%；同一批 85 条修复后 96.7%（另三种工具数字不变，差异全部来自此修复）。

Record list: `benchmark/trna_sets.tsv` (re-run 2026-10-20). The older version dropped MITOS2's `trnM_0(cat)` (several copies of one model) and `trnQ(---)` (anticodon unreadable); after the fix, on the same 85 records, the cnidarian MITOS2 and default figures rose (the other three tools' figures are unchanged; all of the difference comes from this fix).

单用 ARWEN/ARAGORN 在 tRNA 少的类群假阳性极多，故只放 review 文件。刺胞的 FP 多来自 MITOS2，是否为 NCBI 漏注未确证。

ARWEN / ARAGORN alone give very many false positives in clades with few tRNAs, so they are only put in the review file. In cnidarians most FPs come from MITOS2; whether they are tRNAs that NCBI left unannotated is unconfirmed.

**原生核心的 tRNA 与 MITOS2 对比（2026-10-27）**：MitoM-ML（MITOS2 tRNA 模型 + 220 个分支协方差模型，ARWEN/ARAGORN 作证据，学习的判别器与端点）＋ MiTFi 二次确认（p < 0.9 的调用须 MiTFi 也调用；同家族较弱的第二个调用须 p ≥ 0.9）＋ 它没调用的 tRNA 家族由 MiTFi 补 ＋ 脊椎动物取 MiTFi 端点。MiTFi 按 MITOS2 的方式运行（E 值按记录长度计算、tRNA 重叠上限 50 nt，与 MITOS2 自己的 MiTFi 输出逐行一致）；唯一不同：快速搜索未找到的 tRNA 模型的灵敏补搜，MITOS2 用 cmsearch `--max` 搜整条序列（每条几分钟），MitoM 用 `--nohmm --F6 3e-4` 只搜蛋白基因 / rRNA 内部以外的区段（基因两端各 100 nt 仍搜）：542 条需补搜的基准记录 CPU 11.9 → 1.0 小时，同骨海绵 NC_014885.1 297 → 18 s（2 线程），召回率和端点不变，只去掉 14 个落在基因 / rRNA 内的弱 MiTFi 调用（NCBI 均无）。所用记录都没有参与训练（无脊椎动物：所在科未参与训练且不在 refseq89m；脊椎动物：模型从未见过；环节动物：只在 GenBank、不在 RefSeq 的记录）；规则在调参的一半上按事先定的标准选定（精确率比 MITOS2 高 ≥ 0.1 个百分点，再取召回最高），另一半只用来检验。召回率 / 精确率 / 两端与 NCBI 完全一致：

**tRNAs of the native core vs MITOS2 (2026-10-27)**: MitoM-ML (the MITOS2 tRNA models + 220 branch covariance models, ARWEN/ARAGORN as evidence, a learned discriminator and endpoints) + MiTFi second confirmation (a call with p < 0.9 must also be called by MiTFi; a weaker second call of the same family needs p >= 0.9) + tRNA families it did not call are filled in by MiTFi + vertebrates take MiTFi endpoints. MiTFi is run the way MITOS2 runs it (E-values calculated from the record length, tRNA overlap limit 50 nt, identical row by row with MITOS2's own MiTFi output). The only difference is the sensitive back-up search for tRNA models that the fast search did not find: MITOS2 uses cmsearch `--max` over the whole sequence (minutes per record), MitoM uses `--nohmm --F6 3e-4` only on the segments outside the interior of protein genes / rRNAs (100 nt at each gene end are still searched); recall and endpoints are unchanged, and only weak MiTFi calls lying inside genes / rRNAs (none in NCBI) are removed. None of the records used took part in training (invertebrates: the family did not take part in training and is not in refseq89m; vertebrates: never seen by the model; annelids: only in GenBank, not in RefSeq); the rule was chosen on one half (the tuning half) by a criterion set beforehand (precision at least 0.1 percentage points higher than MITOS2, then the highest recall), and the other half was used only for testing. Recall / precision / both ends identical to NCBI:

| 记录 / records | MitoM | MITOS2 |
|---|---|---|
| 无脊椎动物 1,043 条 / 1,043 invertebrate records | 98.55 / 99.04 / 86.48 % | 96.59 / 98.79 / 84.63 % |
| 其中从未用来调参的一半 565 条 / of which the 565 records never used for tuning (the untouched half) | 98.75 / 99.06 / 88.70 % | 96.95 / 98.82 / 87.20 % |
| 脊椎动物 188 条（未用来调参） / 188 vertebrate records (not used for tuning) | 98.38 / 98.55 / 92.66 % | 98.31 / 98.24 / 92.61 % |
| 环节动物（纯 GenBank）501 条 / 501 annelid records (GenBank only) | 99.14 / 99.08 / 83.36 % | 97.61 / 99.01 / 81.47 % |
| 其中从未用来调参的 184 条 / of which the 184 never used for tuning | 98.90 / 98.72 / 80.20 % | 97.87 / 98.71 / 79.69 % |

按门统计（无脊椎动物 1,043 条）的例外只在下表；其余各门（节肢、软体、刺胞、线虫、扁形、棘皮、环节等）三项都高于 MITOS2。报告在各表注明的目录。
Exceptions by phylum (1,043 invertebrate records) are only in the table. On the other invertebrate phyla, all three tRNA figures are higher than MITOS2. Reports are in the directories named on the tables.

| 类群 / group | 与 MITOS2 的差别 / vs MITOS2 |
|---|---|
| 棘头动物 4 条 / Acanthocephala, 4 records | 召回持平 45.5%；精确率 87.0% vs 88.9%（少 1 个 tRNA） / recall equal; precision lower (one tRNA fewer) |
| 动吻动物 1 条 / Kinorhyncha, 1 record | 精确率少 1 个 tRNA / precision lower by one tRNA |

**完整流程与 MITOS2 的对比（2026-10-27，v30：上述全部规则 + NCBI 提交规则，裸序列输入，与 NCBI 注释逐基因比较；两者都完成的记录）**。蛋白基因（PCG）与 rRNA 为两端与 NCBI 完全一致的比例，tRNA 为召回率 / 精确率 / 两端完全一致：

**Full pipeline vs MITOS2 (2026-10-27, v30: all the rules above + the NCBI submission rules, bare-sequence input, compared gene by gene with the NCBI annotation; records that both completed)**. For protein genes (PCG) and rRNA the figure is the share of both ends identical to NCBI; for tRNA it is recall / precision / both ends identical:

| 记录 / records | PCG MitoM / MITOS2 | rRNA MitoM / MITOS2 | tRNA MitoM | tRNA MITOS2 | NCBI ERROR（MitoM） |
|---|---|---|---|---|---|
| 多类群 RefSeq 1,050 条（24 个门） / 1,050 multi-clade RefSeq records (24 phyla) | 84.0 / 59.4 % | 48.9 / 16.9 % | 96.84 / 97.90 / 84.26 % | 95.76 / 97.44 / 83.07 % | 0 条 / 0 records |
| 环节动物 669 条 / 669 annelid records | 81.2 / 51.4 % | 32.9 / 12.7 % | 99.12 / 99.13 / 83.45 % | 97.68 / 99.02 / 81.61 % | 0 条 / 0 records |
| 科未参与训练的无脊椎动物 577 条（不在 refseq89m） / 577 invertebrate records whose family took part in no training (not in refseq89m) | 80.9 / 53.4 % | 33.9 / 19.2 % | 98.76 / 99.06 / 88.70 % | 96.95 / 98.82 / 87.20 % | 0 条 / 0 records |
| 脊椎动物 188 条 / 188 vertebrate records | 88.9 / 72.7 % | 67.6 / 37.0 % | 98.38 / 98.55 / 92.66 % | 98.31 / 98.24 / 92.61 % | 0 条 / 0 records |
| 锚定基因组 13 条 / 13 genomes with rRNA ends mapped experimentally in the literature (anchors) | 87.3 / 47.9 % | 57.7 / 3.8 % | 98.61 / 98.27 / 86.11 % | 96.88 / 98.24 / 82.64 % | 0 条 / 0 records |

不高于 MITOS2 的项目只在下表。按门统计（1,050 条，24 个门），PCG 两端一致、PCG 召回和 rRNA 召回每个门都高于 MITOS2。报告：`benchmark/results/full_v30_2026-10-27/`。
The figures that are not higher than MITOS2 are only in the table. By phylum (1,050 records, 24 phyla), PCG both-ends agreement, PCG recall and rRNA recall are higher than MITOS2 in every phylum. Report: `benchmark/results/full_v30_2026-10-27/`.

| 项目 / item | 结果 / result |
|---|---|
| 纯 GenBank 环节动物、从未调参的 184 条，tRNA 精确率 / 184 untuned GenBank annelids, tRNA precision | 98.70% vs 98.71% |
| 同上，多出来的 NCBI 没有的调用 / extra calls NCBI does not have | 1 / 4000 |
| 其中与 MITOS2 相同的调用 / of those, same call as MITOS2 | 43 |
| MitoM 多找到的真 tRNA / true tRNAs MitoM adds | 41 |
| 动吻动物 1 条 / Kinorhyncha, 1 record | 多 1 个调用 / one extra call |
| 棘头动物 4 条 / Acanthocephala, 4 records | 召回、精确率持平 / recall and precision equal |

ARWEN、ARAGORN、tRNAscan-SE、MiTFi 参数（半全局、重叠、E 值、只按阈值、加分支模型）、多软件投票、重训判别器与端点均测过，见 `benchmark/trna_ml/README.md`。

The parameters of ARWEN, ARAGORN, tRNAscan-SE and MiTFi (semi-global, overlap, E-value, threshold only, adding branch models), multi-program voting, and retraining the discriminator and endpoints were all tested; see `benchmark/trna_ml/README.md`.

## intron 注释 / Introns (mito_intron.py, mito_intron_ref.py)

**注释核心不注释 intron（group I / II；MitoZ 也不能）。** `run` / `annotate` 写完最终 GenBank 后分两步补上 intron，并整合进核心结果（`--intron_fix auto`：找到 exonerate+blastp 时运行；`yes` 强制，`no` 关闭；`--intron_bin <env>/bin`）：

**The annotation core does not annotate introns (group I / II; MitoZ cannot either).** After `run` / `annotate` write the final GenBank, introns are added in two steps and integrated into the core result (`--intron_fix auto`: runs when exonerate + blastp are found; `yes` forces it, `no` turns it off; `--intron_bin <env>/bin`):
1. `mito_intron.py`（从头识别，MFannot 方法，见下）：同源蛋白比对 + 剪接规则定边界；

   `mito_intron.py` (de novo recognition, the MFannot method, see below): homologous-protein alignment + splice rules set the boundaries;
2. `mito_intron_ref.py`（位点迁移）：内含子位点库 `models/intron_sites.*`（193 条 NCBI 带内含子记录、逐个核实后的 309 个内含子，79 KB）中近缘记录的位点映射到本序列：近缘位点重新定位从头识别的边界，补上从头识别漏掉的（小外显子、反式剪接片段、rRNA 内含子）；库中没有的位点保留从头识别结果。报告 `qc/intron_report.tsv`、`qc/intron_ref_report.tsv`。

   `mito_intron_ref.py` (site transfer): sites of close relatives in the intron site library `models/intron_sites.*` (309 introns, each verified, from 193 NCBI records with introns; 79 KB) are mapped onto this sequence: a close relative's site re-places the boundaries recognised de novo and adds what de novo recognition missed (small exons, trans-spliced pieces, rRNA introns); sites not in the library keep the de novo result. Reports: `qc/intron_report.tsv`, `qc/intron_ref_report.tsv`.

### 当前版本（2026-10-29）：核实的金标准 + 位点迁移 / Current version: verified gold standard + site transfer

NCBI 的内含子边界写法因课题组而异（同一位点不同记录可差 1–60 nt，海葵、皮群海葵符合 group I 规则，石珊瑚、角海葵、不少环节动物 group II 不符合），所以先逐个核实本地全部 193 条带内含子的 NCBI 记录（311 个内含子）：剪接后产物（内部终止、起止密码子）、与 12 个他种同源蛋白在拼接点 ±5 aa 内有无插入缺失、group I/II 剪接规则与可滑动窗口、近缘记录同一位点的写法。结果：202 个原样正确，39 个产物相同但边界挪位（改到符合规则的位置），37 个 NCBI 有错且有唯一纠正，16 个证据不足保留 NCBI 写法，2 个 NCBI 内含子特征不是内含子（删除），15 个扁盘动物反式剪接 / rRNA 内含子等只能部分核实；金标准 309 个，其中 97 个位置与 NCBI 不同。详见 `benchmark/results/intron_gold_2026-10-29/`。

NCBI's intron boundaries are written differently by different groups (the same site can differ by 1–60 nt between records; sea anemones and zoanthids fit the group I rules, while stony corals, ceriantharians and many annelids with group II introns do not fit), so all 193 local NCBI records with introns (311 introns) were first verified one by one: the spliced product (internal stops, start and stop codons), insertions or deletions within ±5 aa of the splice junction against 12 homologous proteins of other species, the group I/II splice rules and the window in which the site can slide, and how close relatives write the same site. Result classes: A, correct as written; B, same product but boundary shifted (moved to the position that fits the rules); corrected, NCBI wrong with a unique correction; insufficient evidence, NCBI's version kept; NCBI intron features that are not introns (deleted); trans-spliced / rRNA introns of placozoans, which could be verified only in part. The gold standard has 309 introns, a number of which differ in position from NCBI. Details: `benchmark/results/intron_gold_2026-10-29/`.

真实 `mitom run`（裸序列）对金标准：

Real `mitom run` (bare sequences) against the gold standard:

| 口径 / setting | 完全一致 / exact | ≤6 nt | >6 nt | 漏检 / missed | 76 条无内含子记录的假阳性 / false positives in 76 intron-free records | NCBI ERROR |
|---|---|---|---|---|---|---|
| 默认（完整位点库），194 条 / default (full site library), 194 records | **309 / 309** | 0 | 0 | 0 | 0 | 0 |
| 留一（本记录不进核心参考和位点库，相当于新物种），270 条 / leave-one-out (this record is not in the core reference or the site library, equivalent to a new species), 270 records | **299 / 309** | 4 | 0 | 6 | 0 | 0 |
| 留种 / 留属（同一流程输入，只跑内含子步骤） / leave-one-species / leave-one-genus out (the same pipeline input, only the intron step run) | 296 / 309 | 4 | 3 | 6 | — | — |

留一口径下 A / B / 纠正类全部完全一致（202 / 39 / 37）；差的 10 个是 3 个金标准保留 NCBI 偏离共识写法的珊瑚 nad5 和 7 个库中别无其他记录的扁盘动物反式剪接 / rRNA 内含子。另报出 2 个 NCBI 未注释、很可能真实的内含子（KT989324.1 *Glycera* cox1；丝盘虫 rrnL）。与 NCBI 原注释直接比（311 个）：默认 218 完全一致 / 76 ≤6 nt / 15 >6 nt / 2 漏检（即 NCBI 写错、金标准删除的 2 个）。耗时：从头识别每条 10–20 s，位点迁移 2–5 s。

In the leave-one-out setting, classes A / B / corrected are all exactly identical; the remaining differences are coral nad5 introns for which the gold standard keeps an NCBI boundary that deviates from the consensus, and placozoan trans-spliced / rRNA introns that have no other record in the library. Two further introns that NCBI did not annotate and that are probably real are reported (KT989324.1 *Glycera* cox1; *Trichoplax* rrnL). Compared directly with the original NCBI annotation (311 introns), the introns that are missed are the two that were wrongly written by NCBI and deleted from the gold standard. Run times are given in the Chinese text.

**金标准以外的 NCBI 记录（2026-11-09）**：NCBI 中全部带 intron 特征、≥ 10 kb 的动物线粒体记录里，不在金标准（也不在位点库）中的 134 条（刺胞、苔藓虫、海绵、鸟类）用 `mitom run` 注释，与 NCBI 的内含子比较，位置不同时再比较两种拼接产物与同源蛋白的相似度。刺胞动物和海绵的内含子几乎全部注出，位置不同时多为 NCBI 写错；苔藓虫 Cupuladriidae 的短内含子（约 250 nt，没有 I/II 类特征）旧规则下大多漏检，由此改了无类型内含子的判定（下一段）。数字与逐条结果见 [benchmark/results/intron_ncbi_2026-11-09/](../benchmark/results/intron_ncbi_2026-11-09/README.md)。
**NCBI records outside the gold standard (2026-11-09)**: of all NCBI animal mitochondrial records of >= 10 kb with intron features, those that are not in the gold standard (nor in the site library) - cnidarians, bryozoans, sponges and birds - were annotated with `mitom run` and the introns compared with NCBI's; where the positions differ, the similarity of the two spliced products to homologous proteins was compared. Cnidarian and sponge introns are almost all found, and where the positions differ NCBI's is mostly the wrong one; the short introns (about 250 nt, no group I / II signature) of cupuladriid bryozoans were mostly missed under the old rule, which led to the next paragraph (the judgement of introns of no known type). Numbers and per-record results: [benchmark/results/intron_ncbi_2026-11-09/](../benchmark/results/intron_ncbi_2026-11-09/README.md).

**无类型内含子的判定（2026-11-09，`--intron_rule continuity`，默认）**：Rfam 认不出 I / II 类的内含子（退化的 II 类：环节动物、苔藓虫；Jenkins et al. 2022 *Sci Rep* 12:10889 的 Cupuladriidae 内含子没有 V/VI 结构域、没有 GUGYG…AY，作者也只能按同源蛋白定位）以前要求 I 类剪接分数 ≥ 14 且模型比上同源蛋白 ≥ 80 %——这两个阈值是按数据挑的，I 类剪接规则对 II 类内含子本来就不适用。现在按与类型无关的两条原理判定：
1. **接头连续**（`junction_continuity`）：内含子只去掉 DNA，所以多个同源蛋白在拼接后蛋白的接头处连续——过半同源蛋白跨过接头，同源蛋白缺失或多出的残基中位数 ≤ 1（接头切开的那个密码子）。把分化区或移码区当内含子剪掉时，同源蛋白在接头处缺残基。接头不连续时先在两端各 ±30 nt 内找使其连续的位置。
2. **剪接必要**（`splicing_gain`）：内含子打断编码序列，所以不剪接时基因不能在原位覆盖同源蛋白。两种不剪接的替代（去掉上游外显子、下游外显子在原位按读框向上游读到终止；或反过来）比剪接后的蛋白少 ≥ 20 BLOSUM62 分（约 9 bit；随机序列局部比对的期望分为负）才写入。
有 Rfam 类型的内含子不受这两条限制（结构就是证据）。`--intron_rule splice` 恢复旧规则。另外：注释的 CDS 有内部终止时也作为候选送入内含子步骤（导入的注释、其他引擎）；内含子步骤后 CDS 仍有内部终止的基因在 `intron_report.tsv` 中标出待核查。报告新增 `junction_span`、`junction_missing_aa`、`junction_extra_aa`、`splice_score`、`splicing_gain` 列。验证（金标准 270 条、金标准以外 134 条、按门分层抽取的 362 条无内含子 RefSeq 记录）见 [benchmark/results/intron_algo_2026-11-09/](../benchmark/results/intron_algo_2026-11-09/README.md)。

**Introns of no known type (2026-11-09, `--intron_rule continuity`, default)**: introns that Rfam cannot assign to group I / II (degenerate group II: annelids, bryozoans; the Cupuladriidae introns of Jenkins et al. 2022 *Sci Rep* 12:10889 have no domains V/VI and no GUGYG…AY, and the authors could only place them by homologous protein) used to require a group I splice score >= 14 and a model-to-homologous-protein ratio >= 80 % - thresholds picked to fit the data, and the group I splice rule does not apply to group II introns anyway. They are now decided by two type-independent principles:

1. **Continuous junction** (`junction_continuity`): an intron only removes DNA, so several homologous proteins are continuous across the junction of the spliced protein - more than half of the homologs span the junction, and the median number of residues that the homologs lack or have in excess is <= 1 (the codon cut by the junction). When a divergent or frameshifted region is cut out as an intron, the homologs lack residues at the junction. If the junction is not continuous, a position that makes it continuous is first sought within ±30 nt at each end.
2. **Splicing is necessary** (`splicing_gain`): an intron interrupts the coding sequence, so without splicing the gene cannot cover the homologous protein in place. Each of the two non-spliced alternatives (drop the upstream exon and read the downstream exon in place in frame upstream to a stop; or the reverse) must score >= 20 BLOSUM62 points (about 9 bit; the expected score of a local alignment of random sequences is negative) below the spliced protein before the intron is written.

Introns with an Rfam type are not subject to these two conditions (the structure is the evidence). `--intron_rule splice` restores the old rule. In addition: an annotated CDS with an internal stop is also sent to the intron step as a candidate (imported annotations, other engines); a gene whose CDS still has an internal stop after the intron step is flagged for checking in `intron_report.tsv`. The report has the new columns `junction_span`, `junction_missing_aa`, `junction_extra_aa`, `splice_score`, `splicing_gain`. Validation (the 270-record gold standard, the 134 records outside it, and 362 intron-free RefSeq records stratified by phylum): [benchmark/results/intron_algo_2026-11-09/](../benchmark/results/intron_algo_2026-11-09/README.md).

从头识别（`mito_intron.py`）本轮的改动：跨环状原点的基因在旋转后的副本上建模（原来按 cox1 定向时原点切进 cox1 内含子，环节动物 cox1 内含子因此漏掉一半）；外显子按蛋白顺序串联（刺胞 nad5 的反序片段不再拼成假内含子）；内含子里含有本基因外显子拷贝时拒绝（*Hexamermis* 三个相同 atp6 被拼成 3.7 kb 假内含子）；外显子在两条链上 / 相隔很远的反式剪接基因写成带 `<`/`>` 的部分 CDS；核心的部分基因用核心参考蛋白补搜（*Decemunciger* nad4）；内含子两端可独立微调（按至多 5 个同源蛋白的中位得分）；小外显子插入后两侧内含子重新定边界。`mitom_pipeline` 在内含子步骤后再次旋转，使定向基因的第一个外显子位于序列起点；旋转保留 `<`/`>` 不完整端（此前丢失，反式剪接片段因此出 table2asn ERROR）。

Changes to de novo recognition (`mito_intron.py`) in this round: a gene across the circular origin is modelled on a rotated copy (before, orienting by cox1 let the origin cut into the cox1 intron, so half of the annelid cox1 introns were missed); exons are joined in protein order (reversed-order pieces of cnidarian nad5 are no longer joined into false introns); an intron that contains a copy of the gene's own exon is rejected (three identical atp6 copies in *Hexamermis* had been joined into a false 3.7 kb intron); trans-spliced genes whose exons lie on both strands or far apart are written as partial CDSs with `<`/`>`; partial genes of the core are searched again with the core's reference protein (*Decemunciger* nad4); the two ends of an intron can be fine-tuned independently (by the median score over up to 5 homologous proteins); the boundaries of the introns on both sides are re-set after a small exon is inserted. After the intron step, `mitom_pipeline` rotates again so that the first exon of the orienting gene is at the start of the sequence; the rotation keeps `<`/`>` incomplete ends (they used to be lost, which gave table2asn ERRORs for trans-spliced pieces).

### 从头识别的方法（MFannot 移植）与早期验证（2026-10-09 / 10）/ De novo method and earlier validation

`mito_intron.py`（及 2026-11-08 从中拆出的边界规则 `mito_intron_rules.py`）= MFannot（Beck & Lang，https://github.com/BFL-lab/Mfannot ；数据 https://github.com/BFL-lab/MFannot_data ，GPL-3.0；Lang et al. 2023, doi:10.3389/fpls.2023.1222186）intron 流程的 Python 移植，本文件按 GPL-3.0 发布：

`mito_intron.py` (and `mito_intron_rules.py`, the boundary rules split out of it on 2026-11-08) is a Python port of the intron pipeline of MFannot (Beck & Lang, https://github.com/BFL-lab/Mfannot ; data https://github.com/BFL-lab/MFannot_data , GPL-3.0; Lang et al. 2023, doi:10.3389/fpls.2023.1222186); this file is released under GPL-3.0:

| 步骤 / step | MFannot 对应 / MFannot counterpart | 实现 / implementation |
|---|---|---|
| 蛋白库 / protein library | protein_collections/identified | 随附 `models/mfannot_proteins/`（313 个基因家族，基因名取自 `mt cox1;` 表头） / ships `models/mfannot_proteins/` (313 gene families; gene names taken from the `mt cox1;` header) |
| 预筛 / pre-screening | flip + blast + SelectProteinForExo | ORF ≥40 aa（首个起始密码子起）→ blastp（BLOSUM62，E ≤1e-10）→ 按基因名分组 → 同链相隔 ≥142 nt 的 ≥2 个 ORF 才进入下一步 / ORF >= 40 aa (from the first start codon) -> blastp (BLOSUM62, E <= 1e-10) -> grouped by gene name -> only >= 2 ORFs on the same strand >= 142 nt apart go on to the next step |
| 基因模型 / gene model | FindExonsInHypProtArray / SortAndGroupC4Report | Exonerate protein2genome（minintron 142、maxintron 5000、intronpenalty -35、`models/neutral.pssm`），分数 ≥200 的比对片段按 MFannot 规则串联，可跨越 >5 kb、内含其他基因的 intron / Exonerate protein2genome (minintron 142, maxintron 5000, intronpenalty -35, `models/neutral.pssm`); alignment pieces with score >= 200 are chained by MFannot's rules, and may span introns of > 5 kb that contain other genes |
| 类型 / type | RNAweasel | Infernal + Rfam RF00028 / RF00029（`models/rfam_group_I_II_introns.cm`） |
| 边界 / boundaries | Adjust_all_intronic_junctions | ±6 nt 同移，group I（外显子末端 U、intron 末端 ωG）≥19，group II（GUGYG…AY）≥40 / ±6 nt shifts together; group I (exon end U, intron end ωG) >= 19, group II (GUGYG…AY) >= 40 |
| 整合核心结果 / integration with the core result | — | 基因起止优先采用核心的同框 CDS 端点；无则按同源蛋白延伸到起始/终止（含不完整终止 T/TA） / gene start / end prefer the core's in-frame CDS ends; otherwise extended to a start / stop codon by the homologous protein (including an incomplete T/TA stop) |

与 MFannot 的差异（均写在脚本开头）：类型判定用 Rfam；边界移位以同源蛋白判断（MFannot 以起始模型判断，Exonerate 起始模型错时会拒绝正确边界）；跨环状原点的基因在旋转副本上建模（跨原点的内含子只写在 CDS join 中，table2asn 不接受两段式 intron 特征）；反式剪接基因写成部分 CDS 片段；未移植 MFannot 的小外显子 HMM 步骤。

Differences from MFannot (all written at the top of the script): type calls use Rfam; boundary shifts are judged by the homologous protein (MFannot judges by the start model, and rejects the correct boundary when Exonerate's start model is wrong); a gene across the circular origin is modelled on a rotated copy (an intron across the origin is written only in the CDS join, since table2asn does not accept a two-part intron feature); trans-spliced genes are written as partial CDS pieces; MFannot's small-exon HMM step is not ported.

**验证**（`benchmark/`，2026-10-09）：NCBI RefSeq 中全部带 intron 标注的动物线粒体基因组 85 条（刺胞 82、海绵 2、扁盘 1；真值 138 个 intron = 126 个 intron 特征 + 12 个由 CDS join 推出但缺 intron 特征的），阴性对照环节动物 20 条（含围沙蚕属 4 条）。输入模拟 MitoZ 输出（去掉 intron、拆分基因只留外显子 1）。

**Validation** (`benchmark/`, 2026-10-09): all 85 animal mitogenomes with intron annotation in NCBI RefSeq (cnidarians 82, sponges 2, placozoan 1; truth: 138 introns = 126 intron features + 12 inferred from CDS joins that lack an intron feature), and 20 annelids as negative controls (including 4 *Perinereis*). The input simulates MitoZ output (introns removed, split genes keeping only exon 1).

| | 完全一致 / exact | ≤6 nt | >6 nt | 漏检 / missed | 多报 / extra | 阴性假阳性 / false positives in negatives |
|---|---|---|---|---|---|---|
| 真 MFannot（Docker nbeck/mfannot） / real MFannot (Docker nbeck/mfannot) | 65 | 13 | 7 | 53 | 0 | 0 |
| mito_intron.py | **97** | 23 | 14 | **4** | 1 | **0** |
| mito_intron.py，排除本物种蛋白 / mito_intron.py, own-species proteins excluded | 97 | 23 | 14 | 4 | 1 | 0 |

MFannot 找到的 85 个 intron 中 83 个与本脚本坐标完全相同；MFannot 漏检主要是 >5 kb 的 nad5 大 intron（38 个中漏 38 个，本脚本找到 37 个）。与 NCBI 不一致的 41 处中 19 处 MFannot 与本脚本完全相同（NCBI 边界/终止密码子写法差异）。已知错误：扁盘动物 *Polyplacotoma* 反式剪接 COX1 产生一个 17.7 kb 假 intron（注释带 "UNCONFIRMED … check manually"）；海绵 *Plakinastrella* 两个 intron 间 9 nt 小外显子被合并。**intron 模型仍需人工核对。**

Of the introns MFannot found, nearly all have exactly the same coordinates as this script; the MFannot misses are mainly the large nad5 introns of > 5 kb. Of the places where this script disagrees with NCBI, a substantial part have identical coordinates in MFannot and this script (differences in how NCBI writes the boundary / stop codon). Known errors: in the placozoan *Polyplacotoma*, trans-spliced COX1 produces a 17.7 kb false intron (annotated with "UNCONFIRMED … check manually"); in the sponge *Plakinastrella*, the 9-nt small exon between two introns is merged. **The intron model still needs manual checking.**

**环节动物**（2026-10-09，NCBI 全部环节动物线粒体全基因组 900 条，去重后 733 条，遗传密码表 5）：环节动物线粒体 intron 主要是 cox1 的 **group II** intron（Rfam RF00029），与刺胞动物的 group I 不同。

**Annelids** (2026-10-09; all NCBI annelid complete mitogenomes, deduplicated, genetic code table 5): the introns of annelid mitochondria are mainly **group II** introns of cox1 (Rfam RF00029), unlike the group I introns of cnidarians.
- NCBI → 本脚本：NCBI 有 18 个 intron（3 个 intron 特征 + 15 个只体现在 CDS join 中），找到 15 个；漏 3 个（*Decemunciger* ND4 265 nt 断开处，MFannot 亦未检出）。

  NCBI -> this script: most of NCBI's introns (a few as intron features, most only in the CDS join) are found; a few are missed (*Decemunciger* ND4 at a 265-nt break, which MFannot does not detect either).
- 本脚本 → NCBI：21 个预测；3 个对应 NCBI intron 特征，12 个对应 NCBI CDS join，6 个在 NCBI 完全无注释的记录中（*Glycera fallax* ×2、臭海蛹 *Travisia forbesii*、*Chaetopterus dewysee*、*Notomastus* sp. B、*Podarkeopsis*）。

  This script -> NCBI: the predictions correspond to NCBI intron features, to NCBI CDS joins, or lie in records with no NCBI annotation at all (*Glycera fallax* x2, *Travisia forbesii*, *Chaetopterus dewysee*, *Notomastus* sp. B, *Podarkeopsis*).
- 与真 MFannot：13 个坐标完全一致（含上述无注释记录中的 6 个里的 5 个，以及臭海蛹 *T. sanrikuensis* / *T. forbesii*）；*Decemunciger* ND1 MFannot 判为 2 个 intron（中间 19 nt 小外显子，本脚本合为 1 个）；*Chaetopterus* 5 个约 300 nt、类型未定的预测 MFannot 未报，**未确证**。

  Against real MFannot: the coordinates are identical for most (including most of those in the records with no NCBI annotation, and *Travisia sanrikuensis* / *T. forbesii*); for *Decemunciger* ND1 MFannot calls 2 introns (with a 19-nt small exon between), this script merges them into 1; for *Chaetopterus*, 5 predictions of about 300 nt, type undetermined, are not reported by MFannot and are **unconfirmed**.
- 臭海蛹属 NCBI 线粒体全基因组 2 条：*T. sanrikuensis* LC677172（cox1 intron 9153..10034，NCBI 9152..10033）、*T. forbesii* OZ479805（NCBI 无注释；cox1 intron 7346..8597，与 MFannot 一致）。

  The two NCBI *Travisia* mitogenomes: *T. sanrikuensis* LC677172 (cox1 intron 9153..10034, NCBI 9152..10033) and *T. forbesii* OZ479805 (no NCBI annotation; cox1 intron 7346..8597, the same as MFannot).

### v5（2026-10-10）：NCBI 提交合规 + 证据门槛 / NCBI submission compliance + evidence thresholds

- **NCBI 校验**：`ncbi_validate.py` 把 GenBank 转为 NCBI 提交格式（.fsa + 5 列 .tbl），用 NCBI `table2asn -V vb` 校验（`--only_ours` 只统计本脚本写入的特征）。

  **NCBI validation**: `ncbi_validate.py` converts a GenBank file to the NCBI submission format (.fsa + 5-column .tbl) and validates it with NCBI `table2asn -V vb` (`--only_ours` counts only the features written by this script).
- **合规**：新 CDS 带 `/product`（沿用核心，否则标准名）；T/TA 不完整终止带 `/transl_except=(pos:…,aa:TERM)` 与标准 note；3′ 端延伸到阅读框内第一个终止密码子（允许与相邻基因重叠，`--stop_rule first`，默认）。

  **Compliance**: new CDSs carry `/product` (inherited from the core, otherwise the standard name); T/TA incomplete stops carry `/transl_except=(pos:…,aa:TERM)` and the standard note; the 3′ end is extended to the first stop codon in the reading frame (overlap with the adjacent gene is allowed, `--stop_rule first`, default).
- **证据门槛**（不依赖 NCBI 原注释）：拼接后 CDS 有内部终止、起始密码子非法、未定型 intron 剪接分数 < 14（2026-11-09 起改为接头连续 + 剪接必要，见上文“无类型内含子的判定”）、或与已有的另一标准线粒体蛋白基因 CDS 冲突（内含子编码 ORF 如 heg 除外）→ 不写入 GenBank，只在 `intron_report.tsv`（`confidence` / `written` / `evidence` 列）中列出。

  **Evidence thresholds** (independent of the original NCBI annotation): an intron is not written to the GenBank file, and is only listed in `intron_report.tsv` (columns `confidence` / `written` / `evidence`), when the spliced CDS has an internal stop, the start codon is not valid, an untyped intron has a splice score < 14 (since 2026-11-09 replaced by continuous junction + splicing necessity, see "Introns of no known type" above), or it conflicts with the CDS of another existing standard mitochondrial protein gene (except intron-encoded ORFs such as heg).
- **小外显子**：多个同源蛋白在外显子交界处共同显示缺失时，在 intron 内搜索插入后对所有同源蛋白比对都改善的短外显子（MFannot AnnotateMiniExons 思路）。

  **Small exons**: when several homologous proteins together show a deletion at an exon junction, the intron is searched for a short exon whose insertion improves the alignment to all the homologs (the idea of MFannot AnnotateMiniExons).

| 指标 / metric | 刺胞/海绵/扁盘 85 条（138 intron） / 85 Cnidaria / Porifera / Placozoa records (138 introns) | 环节动物 733 条（NCBI 18 intron） / 733 annelid records (18 NCBI introns) |
|---|---|---|
| 本脚本写入特征的 NCBI 校验 ERROR / NCBI validation ERROR of the features written by this script | **0**（v3：202） | **0** |
| intron 与 NCBI 完全一致 / ≤6 nt / >6 nt / 漏检 / introns vs NCBI: exact / <= 6 nt / > 6 nt / missed | 98 / 23 / 6 / 11 | 找到 11，漏检 7 / found, missed |
| CDS 与 NCBI 完全一致 / CDS identical to NCBI | 99/126（3′ 端 126/126） / 3′ end | — |
| 阴性记录假阳性 / false positives in negative records | 0 | 写入 16 个：11 个有 NCBI 支持，5 个在 NCBI 完全无注释的记录中（均与真 MFannot 坐标一致） / written: some with NCBI support, some in records with no NCBI annotation at all (all with coordinates identical to real MFannot) |

与 NCBI 不一致的 28 处 intron 边界，用独立证据（同源蛋白比对得分、group I ωG 规则）判定：13 处两种边界蛋白相同（写法差异），14 处本脚本更好（多为石珊瑚 COX1，NCBI 边界不符合 ωG 规则），1 处 NCBI 更好。NCBI 原记录自身在 table2asn 下也有错误（刺胞 13 个 `MultiIntervalIntron`，环节 1 个 `AbuttingIntervals`），本脚本不照抄 NCBI 坐标。环节动物漏检 7 个：*Decemunciger* ND4 ×3（MFannot 亦漏）、*Chaetopterus* ×4（约 300 nt、未定型、剪接位点弱、MFannot 未报，NCBI 仅以 CDS join 表示；是否为 intron **未确证**，按门槛不写入）。

Intron boundaries that differ from NCBI were judged by independent evidence (homologous-protein alignment score, the group I ωG rule): in some, the two boundaries give the same protein (a difference in writing); in most, this script is better (mostly stony-coral COX1, where NCBI's boundary does not fit the ωG rule); in a minority, NCBI is better. NCBI's own records also have errors under table2asn (cnidarian `MultiIntervalIntron`, annelid `AbuttingIntervals`), so this script does not copy NCBI coordinates. Annelid introns missed: *Decemunciger* ND4 (MFannot misses them too) and *Chaetopterus* (about 300 nt, untyped, weak splice sites, not reported by MFannot, and NCBI represents them only by a CDS join; whether they are introns is **unconfirmed**, and by the thresholds they are not written).

**近缘参考（借鉴 MitoFinder）**：`--ref_gb <近缘物种 GenBank>`（`mitom run --intron_ref_gb`）把参考 CDS 翻译加入同源蛋白库。MFannot 自带库以真菌/原生生物为主，动物极少、无环节动物，近缘参考可提高环节动物等类群的检出。MitoFinder 本身的 intron 注释只靠 BLASTX 相似性合并片段、不搜边界（其 README 明言），故不替换现有方法。

**Close-relative reference (learned from MitoFinder)**: `--ref_gb <GenBank of a close species>` (`mitom run --intron_ref_gb`) adds the translated reference CDSs to the homologous-protein library. MFannot's built-in library is mainly fungi / protists, with few animals and no annelids, so a close-relative reference improves detection in annelids and other clades. MitoFinder's own intron annotation only merges fragments by BLASTX similarity and does not search for boundaries (its README says so), so it does not replace the current method.

### v6（2026-10-10）：结构模型扩充 + 裸序列验证 / Structure-model expansion + bare-sequence validation

- **类型模型**：Rfam RF00028/RF00029 + 7 个 group II D1–D4 模型（RF01998, RF01999, RF02001, RF02003–05, RF02012；随附）+ Sellés Vidal et al. 2025 的 14 个 group I 亚型模型（IA1–IE3，基于 42,459 个 group I intron 训练，doi:10.1093/bioadv/vbaf020；无许可声明故不随附，首次运行下载到 `~/.cache/mito_intron/`，`--no_subtype_cms` 关闭）。注释写出亚型（如 IB1）。NCBI 真值 intron 判型：刺胞 138 个由 74 → 109 个，环节 27 个不变（18 个 group II）；273 条编码序列 0 误判。

  **Type models**: Rfam RF00028/RF00029 + 7 group II D1–D4 models (RF01998, RF01999, RF02001, RF02003–05, RF02012; shipped) + the 14 group I subtype models (IA1–IE3) of Sellés Vidal et al. 2025 (trained on 42,459 group I introns, doi:10.1093/bioadv/vbaf020; not shipped because no licence is stated, downloaded to `~/.cache/mito_intron/` on first run, `--no_subtype_cms` turns it off). The annotation writes the subtype (e.g. IB1). Typing of the NCBI truth introns: more cnidarian introns are typed, the annelid ones are unchanged (group II); no mistyping on the coding sequences tested.
- **起始/不完整 CDS**：首密码子非起始密码子时上游（再下游）搜索起始密码子；阅读框延伸出序列末端时写为 `<`/`>` 不完整 CDS（NCBI 允许），不再因"无起始密码子"丢弃。

  **Start / incomplete CDS**: if the first codon is not a start codon, a start codon is searched for upstream (then downstream); when the reading frame extends beyond the end of the sequence, an incomplete CDS with `<`/`>` is written (allowed by NCBI) instead of being discarded for "no start codon".
- **借鉴 MitoFinder**：`--ref_gb` 近缘参考；外显子总长 ≤ 3×蛋白长+33。

  **Learned from MitoFinder**: `--ref_gb` close-relative reference; total exon length <= 3 x protein length + 33.
- **其他工具实测/调研**：MITOS2（mamba，refseq89m）在 Metridium 上把 cox1/nad5 注成 -a/-b 片段、外显子边界偏差最多 174 nt、无 intron 特征；在臭海蛹 cox1 上只注出第二外显子。MitoFinder `--allow-intron` 仅按 BLASTX 相似性合并、不搜边界（其 README）。MitoZ、MITOS2、GeSeq、PGA、CPGAVAS2 等均无 intron 模型；目前只有 MFannot/RNAweasel 建模 intron（调研见 benchmark 记录）。

  **Other tools, tested / surveyed**: MITOS2 (mamba, refseq89m) annotated cox1/nad5 of *Metridium* as -a/-b fragments, with exon boundaries off by up to 174 nt and no intron features; on *Travisia* cox1 it annotated only the second exon. MitoFinder `--allow-intron` merges only by BLASTX similarity and does not search for boundaries (its README). MitoZ, MITOS2, GeSeq, PGA and CPGAVAS2 have no intron model; at present only MFannot/RNAweasel model introns (survey in the benchmark records).

| 验证 / test | 结果 / result |
|---|---|
| 刺胞/海绵/扁盘 85 条（138 intron，模拟 MitoZ 输入） / 85 Cnidaria / Porifera / Placozoa records (138 introns, simulated MitoZ input) | 完全一致 98 / ≤6 nt 23 / >6 nt 6 / 漏检 11（同 v5）；高可信 126、中 3、低 4（v5：105 / 22 / 6）；本脚本写入特征 table2asn ERROR 0 / exact / <= 6 nt / > 6 nt / missed (same as v5); confidence high / medium / low; table2asn ERROR of the features written by this script |
| 环节动物 741 条（**裸序列**，去掉全部注释，不用 MitoZ） / 741 annelid records (**bare sequences**, all annotation removed, no MitoZ) | NCBI→本脚本：27 个 intron 中完全一致 6 / ≤6 nt 6 / >6 nt 4 / 漏检 11（4 个 intron 伸到片段末端无第二外显子、3 个 *Decemunciger* ND4 与 4 个 *Chaetopterus* 约 300 nt 未定型，MFannot 均未检出）；本脚本→NCBI：写入 21 个，**0 个与 NCBI 矛盾**（8 intron 特征、8 CDS join、5 个在 NCBI 无注释记录中且与 MFannot 一致），全部高可信 / NCBI -> this script: of 27 introns, exact / <= 6 nt / > 6 nt / missed (the missed ones: introns extending to the end of a fragment with no second exon, *Decemunciger* ND4, and the untyped *Chaetopterus* introns of about 300 nt; MFannot detected none of these); this script -> NCBI: those written, **none contradicts NCBI** (intron features, CDS joins, and records with no NCBI annotation that agree with MFannot), all high confidence |

### v7（2026-10-10）：提速与稳健性 / Speed-up and robustness

- **瓶颈实测**（cProfile）：v6 的时间几乎全在 cmsearch，且只集中在命中的少数模型（IB3 12 s、RF00028 10 s）；与 intron 长度无关（只扫两端 1.5 kb 不更快，判型也不变）。

  **Bottleneck measured** (cProfile): in v6 nearly all the time was spent in cmsearch, and only on the few models that hit (IB3, RF00028); it is independent of intron length (scanning only 1.5 kb at each end is no faster, and the typing does not change).
- **两步判型**：先用 Infernal `--rfam` 快速过滤，仍未定型的再全灵敏度搜索；同一记录所有 intron 一次调用；只搜正链（`--toponly`）、`--threads` 多线程。165 个 NCBI intron 判型数与大类与全灵敏度完全一致（147/165）。

  **Two-step typing**: Infernal `--rfam` as a fast filter first, then a full-sensitivity search for those still untyped; one call for all the introns of a record; only the plus strand (`--toponly`), `--threads` multithreading. The typing counts and major classes of the NCBI introns tested are identical to full sensitivity.
- **小外显子两步搜索**：先单同源蛋白、交界窗口粗筛，再对前 5 个候选完整验证。

  **Two-step small-exon search**: coarse screening with a single homologous protein in the junction window, then full verification of the top 5 candidates.
- **稳健性**：外部程序超时（`--timeout`，默认 900 s）不再卡死；单条记录出错不影响整批（原样保留并在报告中注明）；`--jobs N` 多文件并行（mitom_pipeline 传 `--jobs`）。

  **Robustness**: external programs that time out (`--timeout`, default 900 s) no longer hang; an error in one record does not affect the batch (the record is kept as it was and noted in the report); `--jobs N` runs files in parallel (mitom_pipeline passes `--jobs`).
- **效果**：单记录 *Porites* 42 → 6 s、*Plakinastrella* 103 → 8.5 s、*Travisia forbesii* 4.2 → 1.2 s；全部 846 条验证记录 12 分钟；坐标与 v6 完全相同，判型 +1（NC_040136 nad5，IB1）；测试 464 → 107 s。

  **Effect**: large speed-ups on single records (*Porites*, *Plakinastrella*, *Travisia forbesii*) and on the whole validation set; coordinates identical to v6, one more intron typed (NC_040136 nad5, IB1); the test suite runs faster.

依赖：`conda create -n mitointron -c conda-forge -c bioconda exonerate blast infernal biopython`。

Dependencies: `conda create -n mitointron -c conda-forge -c bioconda exonerate blast infernal biopython`.

```bash
mitom intron -i final_gb_dir -o out_dir --genetic_code 4 --bin /root/miniconda3/envs/mitointron/bin
python3 -m pytest -q tests          # Metridium / Nematostella / Platynereis，排除本物种蛋白
```

## 核基因标记 / Nuclear markers: rDNA and histones (mito_nuc.py)

`mitom nuc annotate` 在 contig（基因组浅测序 / SPAdes 组装、GetOrganelle 组装的 rDNA 单元、PCR 片段、GenBank 记录均可）中注释核 rDNA
（18S、ITS1、5.8S、ITS2、28S）和组蛋白（H1、H2A、H2B、H3、H4，核遗传密码 1），思路同 MitoFinder：用 NCBI 完整注释作参考，
把参考的基因边界通过比对逐列投影到 contig 上。
环状序列（GetOrganelle 结果中名字以 `(circular)` 结尾的路径、GenBank 中 topology 为 circular 的记录）起点任意，可能落在基因中间：
先转到最长的无基因区段（rDNA 单元的 IGS）中点开始再注释，`nuc_annotation.gb` 的定义行注明新起点在原序列中的位置。

`mitom nuc annotate` annotates nuclear rDNA (18S, ITS1, 5.8S, ITS2, 28S) and histones (H1, H2A, H2B, H3, H4, nuclear genetic code 1) in contigs (genome-skimming / SPAdes assemblies, rDNA units assembled by GetOrganelle, PCR fragments and GenBank records all work). The idea is the same as MitoFinder: the complete NCBI annotation is used as the reference, and the reference gene boundaries are projected column by column onto the contig through the alignment. A circular sequence (a path whose name ends with `(circular)` in a GetOrganelle result; a record with circular topology in GenBank) has an arbitrary start that may fall inside a gene: it is first rotated to start at the midpoint of the longest gene-free stretch (the IGS of an rDNA unit) before annotation, and the definition line of `nuc_annotation.gb` states where the new start lies in the original sequence.

```bash
mitom nuc annotate -i scaffolds.fasta -o nuc_out --threads 32            # 默认 --genes rdna,histone
mitom nuc annotate -i my_seqs.gb -o nuc_out --genes rdna
mitom nuc annotate -i get_org_out/anonym.K127.complete.graph1.1.path_sequence.fasta -o nuc_out --genes rdna --organism "Genus species"
# 参考库（随附 models/nuc/，2026-10-23 NCBI）重建或加入自己的注释（如 Geneious 导出的 GenBank）：
mitom nuc fetch -o nuc_refs
mitom nuc build --gb nuc_refs/rdna.gb nuc_refs/histone.gb my_refs.gb -o my_lib
mitom nuc annotate -i scaffolds.fasta -o nuc_out --lib my_lib
```

| 输出 / output | 内容 / contents |
|---|---|
| `nuc_features.tsv` | 每个特征一行：contig、长度、SPAdes 覆盖度、基因、起止、链、5'/3' 是否不完整、最佳参考及物种、一致度、边界依据（aligned 比对得到 / projected 无比对外推）、`aligned_nt`（参考基因内比对到该特征的 nt，越少越要人工核对）；`copy_ratio`（contig 覆盖度 / 单拷贝覆盖度，后者取 SPAdes contig ≥ 1 kb 的长度加权中位数）；组蛋白另有 `class` / `class_evidence`：`replication-dependent`（串联基因簇拷贝，系统发育常用的 H3/H4：H3 第 31/87/89/90 位 = A S V M，或 copy_ratio ≥ 5）、`replacement`（替换变体：H3.3 = S A I G，或含内含子）、`single-copy`（copy_ratio < 2，如 H2A.X）、空 = 证据不足；序列证据优先于覆盖度 / one row per feature: contig, length, SPAdes coverage, gene, start, end, strand, whether the 5'/3' end is incomplete, best reference and species, identity, boundary basis (aligned = from an alignment; projected = extrapolated without alignment), `aligned_nt` (nt of the reference gene aligned to this feature; the fewer, the more it needs manual checking); `copy_ratio` (contig coverage / single-copy coverage, the latter being the length-weighted median of SPAdes contigs >= 1 kb); histones also have `class` / `class_evidence`: `replication-dependent` (copies in tandem gene clusters, the H3/H4 commonly used in phylogenetics: H3 positions 31/87/89/90 = A S V M, or copy_ratio >= 5), `replacement` (replacement variants: H3.3 = S A I G, or with an intron), `single-copy` (copy_ratio < 2, e.g. H2A.X), empty = insufficient evidence; sequence evidence takes precedence over coverage |
| `nuc_annotation.gb` | 有特征的 contig 的 GenBank（rRNA / misc_RNA / CDS，不完整端为 `<` `>`，CDS 带 translation、transl_table 1） / GenBank file of the contigs that have features (rRNA / misc_RNA / CDS, incomplete ends as `<` `>`, CDS with translation and transl_table 1) |
| `submission/nuc.fsa`、`nuc.tbl`、`nuc.sqn` | NCBI 提交文件：定义行按 NCBI 规范自动生成（如 "... 18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence"），来源信息取自输入 GenBank（voucher、isolate、geo_loc_name、collection_date 等）或 `--organism` / `--src`（NCBI source table TSV）；内部说明不写入；rRNA 在负链的序列翻成正向；≥10 N 在 .sqn 中为缺口（100 N = 长度未知）；内部终止的组蛋白写成 gene /pseudo。装有 table2asn 时自动 `-V vb -Z` 校验，结果在 `validation.txt`、`nuc.val`、`nuc.dr`。提交前需加提交者模板（`table2asn -t template.sbt`），NCBI 要求 collection_date 与 geo_loc_name / NCBI submission files: the definition line is generated automatically following NCBI conventions (e.g. "... 18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence"); source information is taken from the input GenBank (voucher, isolate, geo_loc_name, collection_date, etc.) or from `--organism` / `--src` (NCBI source table TSV); internal notes are not written; an rRNA on the minus strand is written in the forward direction; >= 10 N becomes a gap in .sqn (100 N = length unknown); a histone with an internal stop is written as gene /pseudo. When table2asn is installed, `-V vb -Z` validation is run automatically, with the results in `validation.txt`, `nuc.val`, `nuc.dr`. Before submission add the submitter template (`table2asn -t template.sbt`); NCBI requires collection_date and geo_loc_name |
| `genes/<基因>.fasta` | 每个拷贝的序列（组蛋白另有 `.faa`） / sequence of each copy (histones also have `.faa`); the file name is `genes/<gene>.fasta` |
| `consensus/<基因>.fasta`、`consensus.tsv` | 每个基因一条序列（GenBank 输入含多个物种时每个物种一条，`sample` 列为物种名）：多拷贝 rDNA 在组装中散成许多 contig 时，片段先在 scaffold 缺口（≥10 N）处切开，高覆盖片段按重叠拼接（拷贝变异分叉处接覆盖度高的一支）、按参考排序（ITS 没有参考，用样本自己的完整拷贝），再按覆盖度加权逐位投票；与共识一致度 < 99 %（`--consensus_identity`）的片段（其他物种、另一种 rDNA、污染）不用，`consensus.tsv` 列出用了哪些 contig、丢了哪些及原因；没有片段覆盖的区段为 N。组蛋白各拷贝是旁系同源，不合并，按 `class` 取一条：基因簇拷贝优先，其次未定、单拷贝，替换变体最后；同类中完整、覆盖度高者优先。SPAdes 常把基因簇拆成 < 150 nt 的高覆盖碎片：此时报告基因簇碎片（标为不完整），而不用完整的单拷贝变体顶替（MJ20260820041 环节动物 6 个组装：旧规则 3 个取到 H3.3、2 个取到基因簇 H3、3 个 H2A 取到 H2A.X）；`consensus.tsv` 的 `class` 列标明所取拷贝的类别 / one sequence per gene (when a GenBank input contains several species, one per species, with the species name in the `sample` column): when multi-copy rDNA is scattered over many contigs in the assembly, fragments are first cut at scaffold gaps (>= 10 N); high-coverage fragments are joined by overlap (at a fork of copy variation, the branch with higher coverage is taken), ordered by the reference (ITS has no reference, so the sample's own complete copy is used), then voted position by position weighted by coverage; fragments whose identity to the consensus is < 99 % (`--consensus_identity`) (other species, another kind of rDNA, contamination) are not used, and `consensus.tsv` lists which contigs were used and which were dropped and why; stretches covered by no fragment are N. Histone copies are paralogs, are not merged, and one is taken by `class`: gene-cluster copies first, then undetermined, then single-copy, replacement variants last; within a class, complete and high-coverage ones first. SPAdes often splits a gene cluster into < 150 nt high-coverage fragments: then the gene-cluster fragment is reported (marked incomplete) instead of being replaced by a complete single-copy variant (in the 6 MJ20260820041 annelid assemblies the old rule took H3.3 or H2A.X variants in several); the `class` column of `consensus.tsv` states the class of the copy taken |

大组装先用每属一个 rDNA 单元、每门每种组蛋白一个蛋白做快速预筛，只把命中的区域（±20 kb）送入精细注释：
SRR16188826 SPAdes scaffolds（148 万条、7.2 亿 nt）32 线程 2 分 18 秒。

A large assembly is first pre-screened quickly with one rDNA unit per genus and one protein of each histone per phylum, and only the hit regions (±20 kb) are sent to detailed annotation (timing for an SRR16188826 SPAdes assembly: Chinese text).

准确度（`benchmark/nuc_bench.py`，排除同属参考；边界完全一致 / ±5 nt / ±20 nt）：
建库用的 NCBI 378 个物种完整 rDNA 单元 83.2 / 92.7 / 95.2 %（barrnap 0.9 `--kingdom euk` 1.8 / 42.1 / 79.9 %）；
另取不在库里的 71 个物种完整单元作独立测试 42.3 / 82.6 / 90.6 %（barrnap 2.6 / 58.0 / 71.8 %，完全一致率低主要因为占 21 条的一个甲螨研究
有自己的端点约定）。完整组蛋白 CDS 整个模型（起止、外显子、基因名）完全一致：建库集 677 个 97.5 %，独立集 172 个 97.7 %。
另以 RefSeq 精编 rRNA（人、大鼠、猕猴、兔、牛、果蝇、秀丽隐杆线虫）为金标准（`benchmark/nuc_gold_bench.py`）：558 个金标准边界中
NCBI 提交者注释了 9.3 %（其中完全一致 96.2 %），mitom 排除同属参考时全部边界完全一致 85.8 %、±1 nt 94.3 %，在提交者注释过的 53 个边界上
79.2 %（低于提交者）；用完整参考库时 100 %，但此时同种的精编记录就在库里，不是独立检验。库中 28 条精编 rRNA（`build --curated`）比对到边界时
直接决定该边界。
建库时用其他研究的参考交叉检查每条参考的边界，被多数研究否定的边界（如某蜱类研究 28S 3' 端早约 530 nt）不参与投票
（`models/nuc/rdna_qc.tsv`）。边界跟随参考提交者的约定（各研究对 5.8S、28S 端点的约定不完全一致）；参考稀少的类群（海绵、刺胞、粘孢子虫、扁形动物）边界多靠外推，
建议加入本类群的已注释记录。方法与调优记录：`knowledge/topics/nuclear_markers.md`。

When the library is built, the boundaries of each reference are cross-checked against the references of other studies; a boundary rejected by the majority of studies (e.g. the 28S 3' end of a tick study, about 530 nt too early) does not vote (`models/nuc/rdna_qc.tsv`). Boundaries follow the convention of the reference submitter (studies differ in their conventions for the 5.8S and 28S ends); for clades with few references (sponges, cnidarians, myxozoans, flatworms) the boundaries are mostly extrapolated, so annotated records of your own clade should be added. Method and tuning record: `knowledge/topics/nuclear_markers.md`.

Accuracy is measured with `benchmark/nuc_bench.py`, excluding same-genus references (boundaries exact / ±5 nt / ±20 nt), against barrnap 0.9 `--kingdom euk`: on the complete rDNA units of the 378 NCBI species used to build the library, and on an independent test set of 71 complete units of species not in the library (the low exact-match rate there is mainly because one acarid study that accounts for 21 of them has its own endpoint convention); and for the complete histone CDS (the whole model: start / end, exons, gene name) on the library set (677) and the independent set (172). It was also tested against RefSeq-curated rRNAs (human, rat, rhesus, rabbit, cow, *Drosophila*, *C. elegans*) as the gold standard (`benchmark/nuc_gold_bench.py`): NCBI submitters annotated only a part of the 558 gold-standard boundaries; mitom is evaluated with same-genus references excluded, and also on the boundaries the submitters annotated (where it is lower than the submitters); with the full reference library it is 100 %, but then the curated record of the same species is in the library, so this is not an independent test. The 28 curated rRNAs in the library (`build --curated`) decide a boundary directly when they are aligned to it.

rRNA 端点按与实验测定的成熟末端（RefSeq 精编人类等 rRNA）同源的位置确定：先用 BLAST 投影找到基因，再用 Rfam 模型
（RF01960 / RF00002 / RF02543，随附 `models/nuc/rfam/`）精修每个端点——先 profile HMM（快，命中离模型端 ≤3 列且离投影 ≤30 nt 才采用），
其余用协方差模型（±150 nt，命中跑出窗口再 ±600 nt）；模型不能确认且贴着序列末端的端点标为不完整；同种精编参考（≥99.5 %）决定的端点保留。
`evidence` 列：hmm / cm（模型）、curated / aligned / projected（参考）。`--no_cm` 关闭。
独立金标准（RefSeq 精编 rRNA 映射到 8 个物种 336 条记录，558 个边界，排除同属参考，`benchmark/nuc_gold_bench.py`）：精确 96.2 %、±1 nt 99.6 %、
未注出 0.4 %（改动前 85.8 / 96.8 / 1.6 %）；精确不一致的 21 个是秀丽隐杆线虫（模型外侧核苷酸数与哺乳类差 1，排除同属后库中无线虫精编；
线虫 RefSeq 端点为推定值）。CM 金标准（`benchmark/nuc_cm_gold.py`）多鳞虫 + 沙蚕 146 个边界：100 %（与本软件同一模型，循环），
NCBI 提交者 0 %。按提交者注释算（`nuc_bench.py`）一致率下降（38.6 % 精确），原因是各研究的约定不同：如 *Junonia* 18S 5' 端提交者早 4 nt、5.8S 早 6 nt（长型），
mitom 与人类精编同源位置一致。方法与证据：`knowledge/topics/nuclear_markers.md` 第 10、11 节。

rRNA endpoints are set at the positions homologous to the experimentally determined mature ends (RefSeq-curated rRNAs of human etc.): the gene is first found by BLAST projection, then each endpoint is refined with the Rfam models (RF01960 / RF00002 / RF02543, shipped in `models/nuc/rfam/`) - first the profile HMM (fast; a hit is used only when it is <= 3 columns from the model end and <= 30 nt from the projection), otherwise the covariance model (±150 nt; if the hit runs out of the window, ±600 nt); an endpoint that the models cannot confirm and that lies against the sequence end is marked incomplete; an endpoint decided by a curated reference of the same species (>= 99.5 %) is kept. `evidence` column: hmm / cm (model), curated / aligned / projected (reference). `--no_cm` turns it off.

Independent gold standard (RefSeq-curated rRNAs mapped to 336 records of 8 species, 558 boundaries, same-genus references excluded, `benchmark/nuc_gold_bench.py`): the exact mismatches are *C. elegans* (the number of nucleotides outside the model differs by 1 from mammals; after excluding the genus the library has no curated nematode; nematode RefSeq endpoints are putative). CM gold standard (`benchmark/nuc_cm_gold.py`), scale worms + sandworms: the same model as this software, so the comparison is circular; the submitters' annotation does not match it. Agreement computed against the submitters' annotation (`nuc_bench.py`) falls because studies differ in their conventions: e.g. in *Junonia* the submitter's 18S 5' end is 4 nt early and the 5.8S 6 nt early (long form); mitom agrees with the homologous positions of the curated human rRNAs. Method and evidence: `knowledge/topics/nuclear_markers.md` sections 10 and 11.

低配运行（2 线程，`/usr/bin/time`，2026-10-26）：3 条 contig 4.6 s / 132 MB（1 线程 5.1 s / 101 MB）；61 条 GenBank 23 s / 269 MB；
200 条 rDNA 记录 3 分 6 秒 / 533 MB（优化前 24 分 20 秒 / 2.7 GB）；5.6 Mb 混合输入（走预筛）3 分 25 秒 / 608 MB。
30 个 SPAdes 组装（0.005–1.7 Gb，8 线程，2026-10-30）：6 秒至 2.2 分钟，≤ 1.6 GB（改动前 16 线程最长 20 分钟、16 GB；barrnap 0.9 8 线程 2 秒至 4.6 分钟）。

Resource use (time and memory) was measured on small runs with 2 threads (2026-10-26) and on 30 SPAdes assemblies of 0.005–1.7 Gb with 8 threads (2026-10-30), including a comparison with barrnap 0.9; the figures are in the Chinese text.

NCBI 提交校验（table2asn 1.28 `-V vb -Z`，2026-10-25）：多鳞虫 + 沙蚕 50 条、23007 3 条 0 条消息；独立 rDNA 898 条本软件写的特征 0 ERROR
（剩 40 个 rRNA 跨越长度未知缺口的 WARNING，与 NCBI 现有记录写法相同；1 个 ERROR 是 NCBI 原记录 source 文字括号不配对）；组蛋白 2012 条 0 ERROR
（修复前 214：负链 CDS 5'/3' partial 标记互换、partial 端未延到序列末端、H1 假内含子、CDS 按 /gene 名连错拷贝），剩 11 个起点不在序列端的 5' partial WARNING（NCBI 标为假基因或 PCR 片段）。

NCBI submission validation (table2asn 1.28 `-V vb -Z`, 2026-10-25): scale worms + sandworms 50 records and 23007 3 records gave no messages; on independent rDNA records, the features written by this software have no ERROR (what remains are WARNINGs for rRNAs spanning gaps of unknown length, written the same way as in existing NCBI records; the one ERROR is an unmatched parenthesis in the source text of the original NCBI record); histone records have no ERROR (before the fix: 5'/3' partial marks of negative-strand CDSs swapped, partial ends not extended to the sequence end, false H1 introns, CDSs joined across the wrong copy by `/gene` name), with a few 5' partial WARNINGs whose start is not at the sequence end (NCBI marks these as pseudogenes or PCR fragments).

精度（`benchmark/nuc_fp_bench.py`：每个报出的特征是否与记录自身的同名注释重叠，2026-10-25）：独立 rDNA 898 条报出 18S 927 / 5.8S 890 / 28S 895，
未被记录注释的 28S 1 个（鸟类 IGS 中 1.5 kb 类 28S 片段），重复 0（之前 5 个重复、6 个假 28S）；200 条完整线粒体基因组 0 个特征；
2032 条组蛋白记录中 807 个"未注释"的组蛋白经逐条核对 805 个是真基因（提交者只写了 gene 特征或无特征）、名称错 0。
多鳞虫科 / 沙蚕科（NCBI 2026-10-25，排除同属参考）：含 rDNA 的 50 条全部注出；18S 3' 端都在保守的 ...GATCATTA、28S 3' 端都在 ...GATTTGT，
与 NCBI 不一致处多为提交者把 28S 延到序列末端（含 IGS）或 18S 3' 端偏 60–140 nt；5.8S 3' 端与 KY7538xx 等差 7–8 nt（约定差异，未确证）。
沙蚕科 NCBI 没有分段注释的完整单元（只有整段 misc_RNA 的 PP239203 *Perinereis wilsoni*），用 Geneious 注释的 23007-1 对照：
各边界相差 0–3 nt，28S 3' 端 mitom 在保守端 ...GATTTGT，Geneious 早 271 nt。

Polynoidae / Nereididae (NCBI 2026-10-25, same-genus references excluded): all the records that contain rDNA were annotated; the 18S 3' ends all end at the conserved ...GATCATTA and the 28S 3' ends at ...GATTTGT; where they disagree with NCBI, it is mostly that the submitter extended the 28S to the end of the sequence (including the IGS) or put the 18S 3' end off; the 5.8S 3' end differs from KY7538xx etc. (a convention difference, unconfirmed). Nereididae have no NCBI complete unit with segmented annotation (only PP239203 *Perinereis wilsoni* with a single whole-length misc_RNA), so 23007-1 annotated in Geneious was used for comparison: the boundaries differ by only a few nt; at the 28S 3' end mitom is at the conserved end ...GATTTGT, and Geneious is far earlier.

Precision (`benchmark/nuc_fp_bench.py`: whether each reported feature overlaps the record's own same-name annotation, 2026-10-25): on independent rDNA records the features reported are almost all matched by the record's own annotation (one 28S-like fragment of 1.5 kb in the IGS of a bird was not annotated by its record; duplicate and false 28S calls were removed); no feature is reported on complete mitogenomes; of the histones "unannotated" in the records, nearly all were checked one by one and are real genes (the submitter wrote only a gene feature, or no feature), with no wrong names.


## GenBank / EMBL 转 NCBI 提交文件 / Submission files (mito_tbl.py, `mitom gb2tbl`)

把任何已注释的 GenBank / EMBL 转成 NCBI 的 5 列特征表 `.tbl` 和带来源修饰符的 `.fsa`（可全部合并、每个文件一对或每条记录一对），并可用 table2asn 生成 `.sqn`。与同类工具比较（2026-10-30）后吸收了它们的长处：GB2sequin（CHLOROBOX）的 EMBL 输入、来源信息表（`--source_table`，即 table2asn 的 `.src`）和产物名表（`--product_table`）；NCBI EDirect gbf2tbl / gbf2fsa 的管道输入、把 DEFINITION 作为 FASTA 标题、每行 70 个碱基。另外：缺 `/product` 的线粒体基因按基因名或 `/anticodon` 补上（CDS 缺产物名是 NCBI 的 ERROR MissingProteinName）；去掉 INSDC 之外的限定词（Geneious、SnapGene、ApE 留下的）和非 ASCII 文字；跨原点、不完整端点照 NCBI 写法；读不了的文件、无序列的记录报告并跳过。来源修饰符名逐个用 table2asn 1.28 核实过（`voucher`、`title`、`definition` 它不认，故分别映射或另作处理）。

Converts any annotated GenBank / EMBL file into the NCBI 5-column feature table `.tbl` and the `.fsa` with source modifiers (all merged, one pair per file, or one pair per record), and can generate `.sqn` with table2asn. After comparing with similar tools (2026-10-30), their strengths were taken over: from GB2sequin (CHLOROBOX), EMBL input, a source information table (`--source_table`, i.e. table2asn's `.src`) and a product name table (`--product_table`); from NCBI EDirect gbf2tbl / gbf2fsa, pipe input, the DEFINITION as the FASTA title, and 70 bases per line. In addition: a mitochondrial gene lacking `/product` gets one from the gene name or `/anticodon` (a CDS with no product name is the NCBI ERROR MissingProteinName); qualifiers outside INSDC (left by Geneious, SnapGene, ApE) and non-ASCII text are removed; features across the origin and incomplete ends are written as NCBI does; unreadable files and records without a sequence are reported and skipped. The source modifier names were each verified with table2asn 1.28 (it does not recognise `voucher`, `title` and `definition`, so they are mapped or handled separately).

**压力测试**（17 个刁钻样例，每个都经 table2asn 判定；`benchmark/results/gb2tbl_stress_2026-10-30/`）：

**Stress test** (17 tricky cases, each judged by table2asn; `benchmark/results/gb2tbl_stress_2026-10-30/`):

| 样例 / case | 考什么 / what it tests | MitoM gb2tbl ERROR | NCBI gbf2tbl + gbf2fsa ERROR |
|---|---|---|---|
| c01–c03、c09、c13、c14 | 普通记录、内含子、脊椎动物负链、`100^101` 位点、超出末端的特征、怪 ID / 小写 / IUPAC / ordinary records, introns, vertebrate minus strand, a `100^101` site, a feature beyond the end, odd ID / lowercase / IUPAC | 0 | 0 |
| c04 / c05 | 跨原点基因（正链 / 负链） / gene across the origin (plus / minus strand) | 0 | 2 / 2（SeqLocOrder） |
| c06 | 不完整端点 / incomplete ends | 0 | 2（NoStop、SuspiciousFrame） |
| c07 / c17 | 非 ASCII 文字 / Latin-1 编码的文件 / non-ASCII text / a Latin-1 encoded file | 0 / 0 | 6 / 4（NonAsciiAsn） |
| c08 | 编辑软件导出的私有限定词 / private qualifiers exported by editing software | 0 | 4（IllegalDbXref、GeneIdMismatch） |
| c10 | 没有 source、没有物种名 / no source, no species name | 0 | 0 |
| c15 | 真实的 ENA EMBL 文件（AF178678） / a real ENA EMBL file (AF178678) | 0 | 读不了，没有输出 / unreadable, no output |
| c16 | 编辑软件导出：只有基因名、没有产物名 / export from editing software: gene names only, no product names | 0 | 31（MissingProteinName、TrnaCodonWrong 等） |
| c11 / c12 | 没有序列的记录 / 截断的文件 / a record with no sequence / a truncated file | 报告并跳过 / reported and skipped | 输出坏记录 / 没有输出 / outputs a broken record / no output |

15 个有效样例 gb2tbl 全部 0 ERROR，gbf2tbl + gbf2fsa 共 51 个；后者的 FASTA 不带 `[mgcode=]`，每条记录另有 13–14 个遗传密码警告。GB2sequin 只有网页版，未实测；Wolflab gb2tbl.py、Loman genbank_to_tbl.py 是 Python 2 脚本，只处理单条记录，不支持不完整端点和多段位置（读代码判断）。

In all 15 valid cases gb2tbl gives 0 ERROR, against a total of ERRORs with gbf2tbl + gbf2fsa; the latter's FASTA has no `[mgcode=]`, and each record also has genetic-code warnings. GB2sequin is web-only and was not tested; the Wolflab gb2tbl.py and Loman genbank_to_tbl.py are Python 2 scripts that handle only a single record and support neither incomplete ends nor multi-part locations (judged by reading the code).

## 流程图 / Flowcharts

`mitom run` 的 12 步（与 README 第 3 节的步骤表、终端的 `[Step k/12]` 一致；`mitom annotate` 不做定向、第二轮、对比与输出整理，共 8 步）。
The 12 steps of `mitom run` (as the step table in README section 3; `mitom annotate` skips reorient, round 2, compare and output layout: 8 steps):
[中文 SVG](flowcharts/mitom_workflow.zh.svg) · [English SVG](flowcharts/mitom_workflow.en.svg) ·
[中文 Mermaid](flowcharts/mitom_workflow.zh.mmd) · [English Mermaid](flowcharts/mitom_workflow.en.mmd)。
重新生成 / regenerate：`python3 scripts/render_workflow_diagrams.py`
