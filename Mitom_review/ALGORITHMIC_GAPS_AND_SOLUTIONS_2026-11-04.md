# MitoM 算法 / 精度不足与解决方案（只记不修）

**日期：** 2026-11-04  
**仓库：** `/media/cody/18T_R1/yzz/Mitom` @ `1912f7b`（MitoM 1.0.0）  
**约束：** 本文件只记录问题与建议改法，**不改代码**；交给 Claude Code 按现有风格实施。  
**勿重复工程缺陷清单：** Codex R2–R9、`fork`、`/dev/null`、`.tbl`/`pseudo\tNone`、假零 ERROR 等见 `VERIFICATION_REPORT.md`、`NEW_BUGS_2026-11-04.md`。本文件侧重**专业性、算法逻辑、精准度**；其中 P0 有少量与工程正确性重叠的项（直接影响精度主张可信度）。

**背景材料：** `docs/METHODS.md`、`knowledge/experiments/rrna_ceiling/`、`knowledge/topics/taxon_specific_mtdna.md`、`docs/reviews/2026-11-03_publication_recommendations.md`。

---

## 总览

| 优先级 | ID | 主题 | 建议动作 |
|---|---|---|---|
| P0 | S1 | 遗传密码细粒度 + import/run 统一 | 覆盖表 + 共用 `code_of_lineage` |
| P0 | S2 | 验证失败不可变「0 ERROR」/ 零产物 exit 0 | returncode + 产物检查；FATAL 计入；无产物 exit≠0 |
| P0 | S3 | `--topology linear` 禁止 circular 预处理 | 对齐 R5 |
| P1 | S4 | PCG 内含子感知前移 | 常有内含子基因首轮允许间隙；反式剪接单独路径 |
| P1 | S5 | 起始/终止证据分层 | 多候选 + 理由；P11 限分布内；T/TA 可配置 |
| P1 | S6 | rRNA 软边界，少追 exact-NCBI | 区间+不确定度；主指标改检出/可用性 |
| P1 | S7 | 难类群 tRNA | 默认流水线补缺投票；无臂 CM；多拷贝策略 |
| P2 | S8 | NUMT / 核污染 | 深度/共线/多拷贝降权，标 note 不默认删 |
| P2 | S9 | 多分体 / 线性 / 类群特异完备性 | 套件级汇总；类群特异基因清单 |
| P2 | S10 | 不完整组装 | 禁假环化；partial + 缺失列表 |
| P3 | S11 | 每基因置信度 TSV | evidence / alternatives / flag |
| P3 | S12 | 评测去泄漏 + 真值分层 | 数据身份清单；A/B/C 真值；冻结后再改宣传数字 |

**建议实施顺序：** S1–S3 → S4–S7 → S8–S10 → S11–S12（冻结版本后独立评测，再改 README「零 ERROR / 84%」等数字）。

---

## P0

### S1 — 遗传密码：细粒度表 + 路径统一

**问题**  
- 门默认粗：`Trematoda→Platyhelminthes→9`（常见应为 **21**）；`Tunicata→Chordata→2`（常见应为 **13**）；Hemichordata 默认 9（常见 24/33）。  
- `import --clade auto` 与 `run --clade auto` 可对同一被囊记录给出不同 code（D1）。  
- `mito_classify` blastx 写死 `-query_gencode 5`（N5），脊椎/被囊查询框可能偏。

**方案**  
1. 增加亚门/纲级覆盖：`Trematoda→21`、`Tunicata→13`；Hemichordata 按参考命中或密码子投票在 24/33 间选。  
2. `import` 与 `run` **共用**同一套 `code_of_lineage` / 分类决策，禁止两套默认。  
3. blastx：先粗分门，再按候选 code 重跑或改用不依赖错误 gencode 的比对策略；结果写入 `clade_decision.tsv`。  
4. 回归：被囊、吸虫、半索各至少 1 条 fixture，断言 `import`/`run` code 一致且符合覆盖表。

### S2 — 验证与退出码可信（支撑精度宣传）

**问题**  
table2asn 失败仍可报 0 ERROR；FATAL/REJECT 不进阻断；零最终产物仍 exit 0（R2/R3/R7 等）。宣传「2497 零 ERROR」不可信。

**方案**  
1. 读 `returncode`；无 `.sqn`/`.val` 不得记 0 error。  
2. FATAL/REJECT 计入阻断计数。  
3. 无 `final_annotation` 产物 → 非 0 exit。  
4. 修完后**独立重跑**验证集，再更新 README 数字。

### S3 — linear 拓扑不作 circular 预处理

**问题**  
`--topology linear` 仍可能先按环状裁末端重复（R5）。

**方案**  
显式 `linear` 时跳过所有 circular 去重/ORIGIN 扩展；缺 `topology=` 的短 contig 默认 linear，勿默认 circular（N9/D6）。

---

## P1

### S4 — PCG：内含子感知前移

**问题**  
核心 `protein2dna` 无内含子；大内含子/反式剪接靠后置，首轮易缺外显子或假 ORF。留一内含子 299/309，扁盘/反式剪接仍弱。

**方案**  
1. 对已知常有内含子的基因（如刺胞 nad5、相关 cox1）首轮用允许间隙的比对（如 `protein2genome`），或「无内含子命中过短 → 自动开 intron 模式」。  
2. 反式剪接：检测两段 ORF + 同一蛋白 HSP → 标 `trans-spliced`，禁止硬拼单一假 ORF。  
3. 与现有 `mito_intron` 位点库衔接：首轮已切 exon 时，后置只校验/补亚型，避免双重改写冲突。

### S5 — 起始/终止：证据分层，少硬套单一准则

**问题**  
P1/P3/P11 大量拟合 NCBI/腹足习惯；非环节「最上游起始」会错；重叠与 T/TA 门级启发式外推宽。

**方案**  
1. 每个 CDS 保留并报告：`homology_span`、`orf_span`、`punctuated_by_tRNA`、`chosen` + 简短理由。  
2. P11 仅在「有模型且类群在训练分布内」使用；否则保守规则 + `low_confidence`。  
3. T/TA：仅当邻居 tRNA 真伸入且读框允许；门级阈值改为可配置表，写入 `curation_report.tsv`。  
4. 论文/评测区分：裸核心 vs 统一 curate 后（见发表建议 §4）。

### S6 — rRNA：软边界，承认 NCBI 天花板

**问题**  
两端 exact≈48–52%（环节更低）；近缘 NCBI 自洽上限约 60–70%；3′ 常与邻居重叠或远早于邻居。

**方案**  
1. 主指标改为：检出成功、与邻居不冲突、overlap±N nt；两端 exact 作附录。  
2. 输出不确定区间（HMM 匹配起止 vs 邻居裁剪起止）；圈图可画软边界。  
3. 文档引用 `rrna_ceiling`，避免暗示「可逼近 95% exact」。

### S7 — 难类群 tRNA

**问题**  
无臂/极端 tRNA（线虫、棘头等）召回弱；多工具对照在 `mitom trna`，默认流水线不够。

**方案**  
1. 默认流程对**缺失 isotype**：MiTFi → ML →（可选）ARWEN **只补缺、不覆盖高分**。  
2. 无臂结构单独 CM；多拷贝取结构分最高，其余 `pseudo` 或丢弃并写报告。  
3. 工具失败勿当成「0 tRNA」（D7：`check=False`）。

---

## P2

### S8 — NUMT / 核污染

**问题**  
几乎无过滤；WGS / `mitom nuc` 邻域或组装污染可出假 CDS/tRNA。

**方案**  
1. 有 `--infq`/深度：深度≈核水平且密码子/GC 异常 → 可疑。  
2. 无 reads：与核 scaffold 共线，或同基因多拷贝且仅一拷贝在典型 mt 邻域 → 降权。  
3. 默认不删，写 `note=possible NUMT` / misc_feature；summary 计可疑数。

### S9 — 多分体 / 线性 / 类群特异完备性

**问题**  
单分子 + 标准 13 PCG + 22 tRNA 假设过强；水母线性/多分体、钙质海绵、Globodera、无 tRNA 类群等覆盖不足。

**方案**  
1. Medusozoa/Calcarea 等固化线性策略；多分体允许多 FASTA 为「基因组套件」，跨 contig 汇总完备性。  
2. 完备性检查改用**类群特异基因清单**，勿对无 tRNA 类群一律报缺 22 个 tRNA。

### S10 — 不完整组装

**问题**  
假环化、基因被切半、高 N → 漏基因或错误拓扑。

**方案**  
1. `<10 kb`、端点基因截断、高 N：禁止按 circular 裁剪。  
2. 输出 `partial` + 缺失基因列表；summary 明确 incomplete。

---

## P3

### S11 — 每基因置信度 TSV

**问题**  
缺「为何这样定界 + 备选 + 分数」，人工复核成本高。

**方案**  
增加（或扩展现有 reports）逐基因表：  
`accession, gene, start, end, strand, score, evidence, alternatives, flag`  
人工只盯 `flag≠ok`。证据标签建议区分：homology / structure / neighbor-punctuation / NCBI-habit / model-P11。

### S12 — 评测独立性与真值分层

**问题**  
P11 全量训后再跑同批、refseq89m/内含子库/tRNA CM 与评测接触未统一清单；「约定」与「生物学」混淆。

**方案**  
1. 按 `publication_recommendations.md`：每条 accession 记录是否进过蛋白库 / 分类参考 / tRNA CM / P11 / rRNA / 内含子库。  
2. 真值分层：  
   - **A** NCBI 习惯一致率（现有可比基准）  
   - **B** 人工/文献核实难类群子集  
   - **C**（可选）RNA-seq/长读证实 3′  
3. 冻结软件版本后再跑独立留科/时间留出；主文强调 B（+C）；A 与全库数字标注非独立。  
4. 未完成 S2 前，不更新「零 ERROR」宣传。

---

## 与现有审查文件的关系

| 文件 | 内容 |
|---|---|
| `VERIFICATION_REPORT.md` | 本机核实 R1–R9 与共识工程项 |
| `NEW_BUGS_2026-11-04.md` | N/D/E 工程缺陷（只审不修） |
| `docs/reviews/2026-11-03_publication_recommendations.md` | 发表与评测方法建议 |
| **本文件** | 算法/精度差距 + **可执行解决方案**与实施顺序 |

Claude Code 实施时：工程项跟 NEW_BUGS / VERIFICATION；算法项跟本文件 S1–S12；发表数字跟 publication_recommendations，且依赖 S2/S12 完成后再改文档。
