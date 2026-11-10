# 用转录组注释 / 核对线粒体基因组：学到的知识 / RNA-seq for mitogenome annotation: what we learned

日期 2026-11-10/11，机器 cody。数据与数字的出处：`benchmark/results/rna_mappers_2026-11-10/README.md`（比对软件），
`benchmark/results/annelid_rna_2026-11-10/README.md`（环节动物 54 对 / 23 对干净数据），
`benchmark/results/rnaseq_2026-11-10/README.md`（Polydora ORFan）。代码：`mito_rna.py`（`mitom rna`）。

## 1. 生物学基础 / biology

- 动物线粒体基因组整体转录成长的多顺反子前体，在 tRNA 处切开（tRNA punctuation，Ojala et al. 1981, (doi:10.1038/290470a0)
  doi:10.1038/290470a0）；成熟 mRNA 加 poly(A)，补全不完整终止密码子 T / TA。
- 所以测序深度是连续的（前体填满间隔），**深度不能告诉你基因在哪里结束**；能告诉你末端的是 poly(A) 位点： (mito_rna.py:67)
  reads 比到基因组后，末尾一串不比对的 A（反向链为开头一串 T）。
- 表达量：用基因平均深度 / 13 个蛋白基因深度中位数（相对深度）。Selifanova et al. 2023 用这种比较判断 (benchmark/results/rnaseq_2026-11-10/README.md:11)
  *Polydora* ORFan 是否转录（我们在其数据上复现，`rnaseq_2026-11-10`）。

## 2. 数据选择最重要 / data selection matters most

| 规则 / rule | 为什么 / why（出处） |
|---|---|
| 转录组必须来自线粒体基因组所属的同一基因组项目，最好同一标本；配对只用 ncbi-genome-inventory 的关联（`relation_derived` × `linked_assemblies_derived` / `specimen_derived`），**绝不按物种名配** | 按物种名配的 2 对差 8 %（实为别的物种），快速进化的 ATP8、ND6、ND4L、ND3 几乎比不上去，看起来像“不表达” |
| 比对后算 reads 与基因组的差异（edit rate）；同项目 > 2 % 标为 FLAG（污染 / 标错 / 隐存种），结果不用 | 同一标本 30 对最大 1.08 %；标记的 2 对是两种蚯蚓（5.0 %、7.7 %） |
| 学规则先用干净数据：同一标本、差异 < 1 %、PolyA 建库、蛋白基因深度 ≥ 100×、poly(A) reads ≥ 1,000、13 基因齐全、无重复序列 | 54 对中 23 对满足；异常的先放一边，免得把数据问题当成规律 |
| 3′ 端证据要 PolyA（oligo-dT）建库；随机引物 / total RNA 几乎没有 poly(A) reads | 按建库方式的确认率见 annelid README |
| 清洗：fastp 只去接头即可；不要 `--trim_poly_x`（会切掉 poly(A)），不要 `-x` | poly(A) 是证据本身 |
| 只取一部分 reads 也够：前 800 万对就有上千倍线粒体深度（局域网 100 Mbit/s 时这样做） | 27 对用的就是前 800 万对 |

## 3. 比对 / mapping

- 用 **bwa mem**（`mitom rna` 默认）。转录组个体常与组装个体不同（差 0.1–3.5 %），bwa mem 最稳；minimap2 -ax sr (benchmark/results/rna_mappers_2026-11-10/README.md:71)
  在 3 % 差异、61 bp 时丢很多 reads。
- 不要用：端到端比对（bowtie2 默认、segemehl：不软剪切，看不到 poly(A)）；剪接比对（STAR、HISAT2、minimap2 (benchmark/results/rna_mappers_2026-11-10/README.md:71)
  splice：动物线粒体没有剪接内含子，报出的都是假内含子）。
- 只比对到线粒体基因组，MAPQ ≥ 20。数字见 `rna_mappers_2026-11-10/README.md`。 (benchmark/results/rna_mappers_2026-11-10/README.md:71)

## 4. 假信号与偏差 / artifacts（不处理会误判）

| 现象 / artifact | 识别 / how to recognise | 处理 / handling in mito_rna |
|---|---|---|
| oligo-dT 内部引发 | poly(A) 位点前 10 nt 里有 ≥ 7 个连续 A 或共 ≥ 9 个 A（1,047 个真末端中只有 1 个满足）；*Erpobdella* 最强的“末端”7,080 reads 就是基因组上 9 个 A | 标 `internal_priming`，不计 |
| 序列在记录中重复 | 基因 25-mer 在别处也出现（*Themiste* ATP8、ND6 各两份，95–100 %）：reads MAPQ 0 被过滤，深度 = 0 | `duplicated_share` > 0.5：不判断表达，注释说明 |
| 短基因深度低 | 深度随基因长度上升（39 种中 36 种）；ATP8 只有中位数的约 0.07 | 不是注释错误；判断短 ORF 要和 ATP8 / ND4L / ND3 比 |
| 读到 3′ 端的 reads 带接头 / 引物 | A 串后接 SMART oligo-dT 引物序列（未修剪的 reads） | 规则只看软剪切开头的 A 串，不受影响 |
| tRNA 3′ 端的 CCA | tRNA 末端的“poly(A)” | 不计 |
| 物种不符 | edit rate > 2 %（同项目） | FLAG |

## 5. 干净数据学到的规律 / patterns from the 23 clean pairs

- poly(A) reads 停在哪里：蛋白基因 3′ 端 41.1 %，rRNA 3′ 端 47.4 %，其他只有 11.4 %——注释与转录组高度一致。 (benchmark/results/annelid_rna_2026-11-10/README.md:63)
- MitoM 注释为不完整终止（T / TA）的 63 个基因：全部 poly(A) 就在注释末端（0–2 nt）。 (benchmark/results/annelid_rna_2026-11-10/README.md:63)
- 表达排序在物种间可重复（Spearman 中位数 0.82）：COX1–3 > ATP6、CYTB > ND5、ND2、ND4、ND1、ND6 > ND4L、ND3 > ATP8。 (benchmark/results/annelid_rna_2026-11-10/README.md:63)
  所以相对深度可以用来判断 ORF 像不像真基因。
- ND4L 通常没有自己的 poly(A) 末端（20/23）：与 ND4 重叠，同一转录本。ATP8 也常与下游重叠。没有末端证据不等于错。 (benchmark/results/annelid_rna_2026-11-10/README.md:63)
- 完整终止密码子附近 146 个强位点： (benchmark/results/annelid_rna_2026-11-10/README.md:63)

| 情况 / case | 位点 / sites | 转录本末端 / transcript end |
|---|---|---|
| 不与下游重叠 | 71 | 终止密码子后几个到约 30 nt（短 3′ UTR） |
| 不与下游重叠 | 14 | 正好在终止密码子 |
| 与下游 tRNA 重叠 | 27 | 在 tRNA 起点切开，留 T / TA |
| 与下游 tRNA 重叠 | 23 | 保留完整终止密码子 |
| 与下游蛋白基因重叠 | 8 | 在下游基因起点切开，留 T / TA |

- **转录组的真正价值**：基因组序列本身定不了的末端。终止密码子与下游 tRNA 重叠时，两种结局都常见，只有 poly(A) (benchmark/results/annelid_rna_2026-11-10/README.md:36)
  能逐个确定；MitoM 标为 "stop codon not determined" 的末端（*W. pigra* ND6、*Sabella* ND1、ND3）被 poly(A) 补上，
  长度都与同科中位数一致。

## 6. 现行规则（`mito_rna.py --refine`）/ current rules

1. 确认：poly(A) ≥ 3 reads 在注释 3′ 端 ±10 nt；或完整终止密码子后 11–30 nt、中间无其他基因（短 3′ UTR）。
2. 缩短：原末端未被确认、基因末 90 nt 内有 ≥ 10 reads 的位点，截断后框内留 T / TA。已确认的末端从不移动
   （*Travisia* COX3 的教训：内部位点会切掉保守 C 端 …YWWGS）。
3. 延伸：基因没有终止密码子时，延伸到其后 90 nt 内留下框内 T / TA 的位点，中间无终止密码子、无其他基因。
4. rRNA 3′ 端：只移到 30 nt 内、占附近 poly(A) reads ≥ 90 % 的位点，原末端无 reads。
5. 每个修正都用同科已注释蛋白的长度核对（目前 5 处全部一致）。

## 9. 内含子：用剪接连接 reads 检验 / introns checked with junction reads

方法：不用剪接比对软件（会报假内含子），只检验 MitoM 预测的内含子：在 reads 中直接数 30 nt 探针——剪接后的
外显子|外显子，与未剪接的外显子|内含子、内含子|外显子（两条链）。数字见 annelid README 的内含子表。

- 4 个 MitoM 内含子全部被证实、剪接位点精确：*Chaetopterus dewysee* COX1（320 nt，此前“未确证”）、*Glycera fallax*
  COX1 两个、*Travisia forbesii* COX1；剪接后连接 reads 4,895–102,944 条，远多于未剪接的 (benchmark/results/annelid_rna_2026-11-10/README.md:93)
- *Travisia* ATP8 的 "1 nt 内含子"（MitoM 用 join 跨过一个移码）：剪接连接 0 条、未剪接 533 / 677 条——RNA 保留了这个
  碱基，不是内含子；该 ATP8 需人工核查（移码或位置不对） (benchmark/results/annelid_rna_2026-11-10/README.md:93)

## 10. 外部资料核查 / checked claims from outside material (2026-11-11)

- “MitoZ-RNA / MToolBox 用 poly(A) 定转录边界”：未找到依据；MitoZ 无 RNA-seq 输入，MToolBox 是人类变异注释 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6582343/)
- “MFannot 的 RNA-seq 模块（HISAT2 + StringTie）”：MFannot 论文说现有细胞器注释软件都不用 RNA 数据 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10352661/)
- PMGA 存在（植物线粒体注释网站），所称 AUES、AEWDG 算法查不到：未确证
- HISAT2 + StringTie 用于动物线粒体：我们的比对评测中 HISAT2 最差，前体连续转录使 StringTie 拼不出基因断点 (benchmark/results/rna_mappers_2026-11-10/README.md:71)
- 值得用的：剪接连接验证内含子（上节）；Iso-Seq 全长转录本（环节动物仅 3 个组装有，未确证其线粒体覆盖）；RNA 编辑在动物线粒体少见，零星编辑位点不影响 2 % 的物种核对：未确证

## 7. 未解决 / open

- 终止密码子与下游 tRNA 重叠时由 poly(A) 决定末端：规律已学到（上表），规则尚未写入；计划先在 23 对干净数据上定参数， (mito_rna.py:301)
  再拿其余 29 对作独立检验。
- *Travisia* ATP8 的 1 nt join 不被 RNA 支持，需人工核查 (benchmark/results/annelid_rna_2026-11-10/README.md:93)
- *Protula*、*Arenicola* rrnL 内部 190–330 nt 处的强位点：原因不明。 (mito_rna.py:72)
- rrnS 内部弱位点（≤ 70 reads）：可能是降解中间体，未使用。 (mito_rna.py:70)
- 两种蚯蚓（5–8 %）与 *Themiste* 重复：需要分类学判断。 (mito_rna.py:74)

## 8. 工作上的教训 / lessons for the workflow

- 优先用用户自己的工具（ncbi-genome-inventory 配对、organise-sra 下载 / locate），功能缺失报 issue (benchmark/results/annelid_rna_2026-11-10/README.md:93)
  （ncbi-genome-inventory #1、#2；organise_sra #4），不自己绕开。
- 结论之前先看 reads 差异：低深度可能是物种不符或序列重复，而不是不表达。 (benchmark/results/annelid_rna_2026-11-10/README.md:27)
- 每个数字出自哪一批数据要写清（登录号、run、项目），见 `annelid_rna_2026-11-10/datasets.tsv`。 (benchmark/results/annelid_rna_2026-11-10/select_pairs.py:23)
