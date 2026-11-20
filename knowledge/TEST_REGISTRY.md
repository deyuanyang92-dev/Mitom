# 测试登记表 / Test registry

**做任何测试前先查这里；测完立即在对应分组末尾加一行。** 目的：两台电脑、多个会话不重复测试。
**Check here before running any test; add a row right after it.** Two machines, many sessions: no repeated tests.

- machine：`deyuan` = 另一台电脑（结果路径 `/media/deyuan/18T_B/yzz/mitom_work/`），`cody` = 本机（`/media/cody/18T_R1/`），`?` = 来源未写。
- 每行：日期、机器、问题、数据、比较的工具与版本、关键数字（照抄来源）、结论 / 动作（commit）、详细记录位置。
- 详细结果仍放 `benchmark/results/<名称_日期>/README.md`（含 run.sh）；这里只放一行索引。
- 另一台电脑的测试若只在本地、未提交到仓库，这里看不到：测完请 commit + push（至少 README 和汇总表）。
- “PCG exact” = 蛋白基因两端与 NCBI 完全一致。真值默认是 NCBI RefSeq 注释。

建立：2026-11-08，从 knowledge/benchmark_findings.md、benchmark/results/*、docs/benchmark_*.md、docs/METHODS.md、
knowledge/experiments/*、Mitom_review/*、CHANGELOG、git log 汇总（数字照抄来源，未重新运行）。

## 1. 准确度基准（对 NCBI / 其他工具）/ Accuracy benchmarks

| date | machine | question | data | tools & versions | key result | conclusion / action (commit) | where recorded |
|---|---|---|---|---|---|---|---|
| 2026-10-20 | ? | Metazoa 优化用样本集 | 173 RefSeq NC_，30 门 | MitoM（MitoZ 路线 mitom_v1） | 集合建立；错误类型驱动 P1/P3/P7/D1 规则 | 07ccd14 | benchmark/sets/mito_bench_2026-10.tsv; docs/benchmark_metazoa_2026-10.md |
| 2026-10-21 | deyuan | MitoM vs 3 工具，环节动物 | 40 环节动物 RefSeq（15 目 29 科，code 5） | MitoM（MitoZ 3.6 + MITOS2 核心）、MitoZ 3.6、MITOS2 2.1.10、MitoFinder 1.4.2 | PCG exact（517）80.1 / 63.2 / 53.6 / 48.4 %；rRNA exact 30.0/17.5/20.0/12.5 %；明显错误记录 16/36/39/34；中位 100/244/353/206 s | P3 终止修正、D1 去重、ND5/ND4L 上游起始；PCG exact 383→407→414/517 | docs/benchmark_annelida_2026-10.md; benchmark/results/annelida_2026-10/ |
| 2026-10-21 | ? | `--if_reannotation_2 auto` vs `yes` | 6 条开发记录（6 门） | MitoM b4c1b93 | 0 个特征不同；NC_026995.1 293.8→169.0 s | auto 设为默认（commit ?） | knowledge/experiments/reannotation/ab_dev6.tsv |
| 2026-10-22 | deyuan | 4 工具小集对比 | 24 RefSeq（8 门×3） | MitoM（MitoZ 路线）、MitoZ 3.6、MITOS2 2.1.10、MitoFinder 1.4.2 | 17 条共同：PCG exact 82.1/78.4/50.5/59.2 %；rRNA exact 37.5/25.0/15.6/9.4 %；有错误记录 7/16/17/17；中位 175/249/385/204 s | 7455075 | docs/benchmark_metazoa_2026-10.md; benchmark/results/metazoa_mini_2026-10/ |
| 2026-10-22 | deyuan | 原生核心 vs MitoZ 路线（24） | 同上 | mitom_native / mitom / MitoZ | PCG exact 247/240/174；中位 109/175/249 s；内存 0.32/0.73/0.49 GB | 原生核心为默认（5bac38a, b3b1336） | benchmark/results/metazoa_mini_2026-10/report_native.md |
| 2026-10-22 | deyuan | 原生核心 vs MitoZ 路线（173） | 173 RefSeq，30 门 | mitom_native vs mitom_v1 | PCG exact 79.0 / 66.7 %；漏 21/39；NCBI 错误记录 11/113；25/30 门领先 | 9e888d0 | benchmark/results/metazoa_173_2026-10/ |
| 2026-10-22 | deyuan | 起止规则留一法调参 | 255 RefSeq | MitoM native | PCG 两端 2186→2261；P1 移动 ≤15 nt 31 对/12 错 | MAX_START_MOVE 15、MAX_PCG_OVERLAP 30（695d570） | docs/benchmark_six_tools_2026-10-22.md |
| 2026-10-22 | deyuan | 5 工具（含 MITOS1）+ 棘皮 | 24 集 + 6 棘皮；MITOS1 8 条 | MitoM、MitoZ 3.6、MITOS2 2.1.10、MitoFinder 1.4.2、MITOS1 1.1.7 | 17 条：PCG exact 84.9/78.4/50.9/59.2 %（留一 82.6）；棘皮 94.9 vs 91.5 %；MITOS1 4/8 崩溃 | 7455075 | docs/benchmark_six_tools_2026-10-22.md; benchmark/results/six_tools_2026-10-22/ |
| 2026-10-27 | deyuan | 完整流程 v30 vs MITOS2（大集） | 1050（24 门）/669 环节/577 未训练科/188 脊椎/13 锚定 | MitoM v30 vs MITOS2 2.1.10 | PCG exact 84.0/59.4、81.2/51.4、80.9/53.4、88.9/72.7、87.3/47.9 %；rRNA 48.9/16.9 …；ERROR 记录全部 0 | 6575782, 6ab18fd | benchmark/results/full_v30_2026-10-27/; docs/METHODS.md |
| 2026-11-07 | ? | 审查修复 3–5 轮前后真实数据对比 | 49 NCBI 记录（8 门），FASTA+GenBank，8 次运行 | MitoM eaa8ef5 vs e930614 | 坐标 0 差异；ERROR 0 | 无回归（a70f8eb） | Mitom_review/FIX_STATUS_2026-11-05.md |
| ? | ? | 真实 MitoZ 结果经 curate | 31 条鳞沙蚕 MitoZ 注释 | MitoZ vs mito_curate | PCG 5′/3′ 同 NCBI 253/211→265/254；table2asn ERROR 128→4 | 软体动物数据未测 | docs/METHODS.md（curation） |

## 2. MITOS2 / MitoZ / import / curate

| date | machine | question | data | tools & versions | key result | conclusion / action (commit) | where recorded |
|---|---|---|---|---|---|---|---|
| 2026-10-20 | ? | MITOS2 蛋白基因边界准确度 | 6 基因组（76 位点，含 NC_000931.1） | MITOS2 2.1.x | 5′ 一致 66 %，3′ 78 %，两端 51 % | 用于边界方法设计 | knowledge/experiments/mitos2_pcg/; knowledge/topics/pcg_boundary_methods.md |
| 2026-10-28 | deyuan | **curate 对各工具结果的作用** | 24 RefSeq（三者都完成 18） | MitoM、MitoZ 3.6、MITOS2 2.1.10，各自 + `mitom import` | PCG exact：MITOS2 50.0→56.5 %，MitoZ 75.7→78.7 %，MitoM 58.7→83.0 %；rRNA ≤3 nt 24–29→62–65 %；ERROR 记录 18/18→0/18 | curate 不能补漏注基因；MITOS2 端点仍偏（48 个 >30 nt） | benchmark/results/curation_cmp_2026-10-28/README.md |
| 2026-11-08 | cody | MITOS2 调用 / 导入 / curate 实际演示（**与上一行重复，仅单条**） | NC_000931 去注释 | MITOS2 2.1.10（`--engine mitos2 --mitos_args`）、2.1.9 单独运行 + import | PCG exact 4/13→6/13（MitoM 自身 10/13）；rRNA 0/2→2/2；ERROR 2→0；295 s +8 s vs 35 s；curate 保留 MITOS2 越过下游 tRNA 5–20 nt 的完整终止（cox1/atp8/cox3/nad6） | `--mitos_args`、runmitos.py 支持、docs/MITOS2.md（766653c） | docs/MITOS2.md |
| 2026-11-08 | cody | 冻结评测流程冒烟测试（freeze.sh → run_frozen.sh → stats） | NC_000931，干净 worktree | mitom_native、mitom_blind（只给序列） | 两者都 10/13 蛋白基因两端一致、rRNA/tRNA 全对；脏工作区时 run_frozen 拒绝运行；manifest 记下 MitoZ 自动找到的是 3.5（FastMitoAssembler 环境） | 流程可用（benchmark/frozen/） | benchmark/frozen/PLAN.md |
| 2026-11-08 | cody | 现有 RefSeq 能否组成独立测试集 | `/media/cody/18T_R1/NCBI/Server1/annelida_ref_curate_250320.gb`（170 条环节动物 NC_） | frozen_set.py | 169/170 记录号已在仓库出现；所有属、科见过；物种未见过 0 条 | 测试集改为 RefSeq+GenBank、物种未见过为主；属未见过另报 | benchmark/frozen/PLAN.md §2 |

## 3. tRNA

| date | machine | question | data | tools & versions | key result | conclusion / action (commit) | where recorded |
|---|---|---|---|---|---|---|---|
| 2026-10-20 | ? | tRNA 各工具 vs NCBI | 85 刺胞/海绵/扁盘；40 环节 | MITOS2、ARWEN、ARAGORN、tRNAscan-SE 2 | MITOS2 96.7 %/FP 21、97.3 %/FP 3；ARWEN/ARAGORN FP 348/462、79/117；tRNAscan 61.2/46.1 % | MITOS2 为默认，其余只作参考 | docs/METHODS.md; benchmark/trna_sets.tsv |
| 2026-10-20 | ? | tRNAscan-SE `-O` | NC_000931.1 | tRNAscan-SE | 不加 `-g` 时 Trp(TCA) 判为 Sup | 记录 | knowledge/experiments/trnascan/ |
| 2026-10-27 | deyuan | MitoM-ML vs MITOS2（未训练类群） | 1,043 无脊椎（565 未调参）；188 脊椎；501 环节 | MitoM-ML+MiTFi vs MITOS2 | 召回/精确/两端：无脊椎 98.55/99.04/86.48 vs 96.59/98.79/84.63 %；脊椎 98.38/98.55/92.66 vs 98.31/98.24/92.61 | MitoM-ML 默认；补搜 CPU 11.9→1.0 h（4e9b99e, 6ab18fd） | docs/METHODS.md |
| 2026-10-27 | deyuan | tRNA 锚定集 | 13 锚定基因组 | native / v27 / MITOS2 | 召回 98.61/97.57/96.88 % | — | benchmark/results/full_v30_2026-10-27/ |
| 2026-10-31 | ? | Leu/Ser 编号 | 规则测试 | MITOS2 2.1.10 规则 | 所有类群按 MITOS2（与 NCBI 脊椎动物相反），有意 | 655a75c | docs/METHODS.md |

| 2026-11-19 | cody | MitoM-ML 环状起点（旋转不变）；MiTFi 重复扫描能否跳过 | tests/data 8 环状 + 3 issues 记录；4 记录 219 个旋转 | trna_ml / mito_trna（exp/trna-circ-mitfi-2026-11-19），Infernal 1.1，mitotools | 旋转不一致 8/219 → 0/219（MITOM_TRNA_CIRC=1）；未旋转 8/8 位置相同；线性不加尾产生 11 个片段调用（未采用）；MiTFi 在 2/10 记录跳过，每条省 ~2 s / 57 s，core.tsv 10/10 相同；golden 两种设置均相同（intron_fix.log 为基线已有）；pytest 856 通过 | CIRC_FIX 实验开关默认关（保留集准确率需在 deyuan 测）；MiTFi 按需默认开 | benchmark/results/trna_circ_2026-11-19/README.md |
## 4. rRNA

| date | machine | question | data | tools & versions | key result | conclusion / action (commit) | where recorded |
|---|---|---|---|---|---|---|---|
| 2026-10-21 | ? | 按门 HMM vs CM vs BLAST | 93 开发记录，186 rRNA | CM、HMM、BLAST | 留种找到 cm 165、hmm 153、blast 145、cm+hmm 175/186 | HMM + CM 兜底（4d0846b, 1fba3f2） | knowledge/experiments/rrna_ref/ |
| 2026-10-21 | ? | EMMA HMM vs MITOS2 CM | 22 门 184 特征 | EMMA vs CM | 83 vs 158/184 | 不采用 EMMA | knowledge/experiments/emma_hmm/ |
| 2026-10-21 | ? | CM glocal rRNA 端点 | 91 记录 | MITOS2 CM | mt12 两端 ≤3 nt 25/91 | 无单一规则 | knowledge/experiments/rrna_ends/ |
| 2026-10-26 | deyuan | rRNA 端点能否到 95 % | 1,050 + 669，1,719 留一 | MitoM vs MitoZ vs MITOS2 | NCBI 近缘自身两端一致仅 72 %；两端 exact 48.3/32.5 %；共识 ≤3 nt 94.0/95.7 % | 对 NCBI 不可达（NCBI 不一致）（5368d07） | knowledge/experiments/rrna_ceiling/README.md |
| 2026-10-26 | ? | 文献实验 rRNA 端点 | 13 基因组 32 端点 | MitoM、NCBI、MitoZ、MITOS2 | 16 个 paper 端点 exact 12/14/2/2 | 仅记录 | rrna_ceiling/README.md §6 |
| 2026-10-26 | ? | RNA-seq 定 rRNA 端点 | ENA 5 条 | RNA-seq | 普通 RNA-seq 无单碱基分辨率 | 放弃 | rrna_ceiling/README.md |
| 2026-11-04/05 | cody | 去掉 rRNA 的 NC_000931 用 `rrna-hmm search` | NC_000931 | mito_rrna_hmm | 补回 2 个（覆盖 0.99/0.97） | — | Mitom_review/VERIFICATION_REPORT.md |

## 5. 内含子 / Introns

| date | machine | question | data | tools & versions | key result | conclusion / action (commit) | where recorded |
|---|---|---|---|---|---|---|---|
| 2026-10-09 | ? | 从头内含子 vs MFannot | 85 记录 138 内含子；20 环节阴性 | mito_intron.py vs MFannot（Docker） | exact 97 vs 65；漏 4 vs 53；阴性 FP 0 | 采用 MFannot 移植 | docs/METHODS.md |
| 2026-10-29 | deyuan | 内含子金标准 + 端到端 | 193 记录 309 内含子；76 无内含子 | `mitom run` | 309/309；留一 299；0 假阳性；ERROR 0 | 位点库（f270b92） | benchmark/results/intron_gold_2026-10-29/ |
| 2026-10-30 | ? | 起点在基因内的旋转测试 | 270 记录 × 5 种旋转 | mitom_lo | near3 232→295；near5 190(+12 假)→295 | 738f84d, 66e0e09, 2dae5b3 | benchmark/results/intron_rotation_2026-10-30/ |
| 2026-10-30 | ? | MJ20260820041 cox1 内含子 | 22 MJ 线粒体基因组 | MitoM | 18 条加了 cox1 内含子；0 ERROR | 已核实 | docs/METHODS.md |
| 2026-11-09 | cody | NCBI 中金标准以外带内含子的动物线粒体基因组能否注出 | NCBI esearch intron[FKEY] ≥10 kb 动物 291 条中不在金标准/位点库的 134 条（刺胞 122、苔藓虫 5、海绵 4、鸟 3） | `mitom run`（默认与留一）f04d839，mitotools；eval_introns.py；blastp refseq89m 比较拼接产物 | 刺胞 189 个：exact 119 / ≤6 39 / >6 30 / 漏 1（非可信内含子）/ 多 1（NCBI 漏注）；位置不同时 MitoM 产物更好 48、相同 19、NCBI 更好 2；海绵 6/6；苔藓虫 50 个只注出 12（Cupuladriidae ~250 nt 短内含子漏，23 个有低置信候选未写）；留一与默认相同 | 刺胞/海绵可用；苔藓虫短内含子是缺口 | benchmark/results/intron_ncbi_2026-11-09/README.md |
| 2026-11-09 | cody | 无类型内含子改按“接头连续 + 剪接必要”判定，是否更通用、不增多报 | 金标准 270 条（NCBI 重新下载，309 个金标准内含子）、金标准以外 134 条、按门分层 362 条无内含子 RefSeq；留一，只重跑内含子步骤 | mito_intron.py --untyped_rule continuity vs splice（f04d839 之后），mitotools；eval_introns.py | 金标准新旧相同 299/4/0/6、多报 2；苔藓虫 12→31/50（多报 2 个为 NCBI 漏注的 cox3 内含子）；刺胞/海绵不变；无内含子 362 条多报 0→0；golden GenBank 不变 | 新规则设为默认 | benchmark/results/intron_algo_2026-11-09/README.md |
| 2026-10-30 | ? | NCR 注释 | 22 MJ；20 NCBI Hesionidae | `--ncr` | 26 NCR，0 ERROR | a50402d | docs/METHODS.md |

## 6. 拓扑 / 方向 / 类群 / Topology, orientation, clade

| date | machine | question | data | tools & versions | key result | conclusion / action (commit) | where recorded |
|---|---|---|---|---|---|---|---|
| 2026-10-30 | ? | 全部首尾相接 vs 线状注释 | 93 记录 | MitoM | PCG 差 ≤2 点；+22 % 时间 | 不全部相接，只查端点 | docs/METHODS.md |
| 2026-10-30 | ? | `--topology auto` vs 标签 | 93 + 20 基因内切开 + 20 缺口 | MitoM | 基因内切开 20/20 环状，PCG 78.9 vs 72.9 % | auto 默认（de5ec0c, c2fe1c3） | docs/METHODS.md |
| 2026-10-31 | ? | 无标签判断环/线 | 91 基因组 × 4 = 360 序列 | MitoM a50402d/2dae5b3/cc8e82f | 正确 64.4→85.8→84.3 % | 去掉 21 tRNA 规则（cc8e82f） | benchmark/results/topology_unlabelled_2026-10-31/ |
| 2026-10-31 | ? | `--clade auto` | 238 RefSeq（29 门）+ 片段 / 突变 / N | mito_classify.py | 密码 93.3 %；门 对/不定/错 63.0/32.4/4.6 % | 不定时用 Animal（79b4b8d） | docs/METHODS.md |
| ? | ? | 起始密码子模型 P11 | 1,683 留一记录 | 梯度提升树 | exact 起点 89.1→90.7 % | fc56b1e | CHANGELOG.md |

## 7. NCBI 提交 / table2asn / gb2tbl

| date | machine | question | data | tools & versions | key result | conclusion / action (commit) | where recorded |
|---|---|---|---|---|---|---|---|
| 2026-10-20 | ? | table2asn 对遗传密码 / transl_except / T-TA 的行为 | 手工 fsa/tbl 用例 | table2asn | 需 transl_except 才去掉 NoStop | 采用 | knowledge/experiments/table2asn/ |
| 2026-10-22 | deyuan | `/transl_except` 坐标迁移 | 173 记录 | MitoM 评测 | ERROR 记录 110→37 | 修评测 | knowledge/benchmark_findings.md |
| 2026-10-27 | deyuan | 全部基准输出 table2asn | 2,497 记录 | MitoM v30 | ERROR 记录 0 | — | benchmark/results/full_v30_2026-10-27/ncbi_errors.txt |
| 2026-10-30 | deyuan | gb2tbl vs NCBI gbf2tbl | 17 手工用例 | `mitom gb2tbl` vs gbf2tbl/gbf2fsa | 15 个有效用例 ERROR 0 vs 51 | d3247dc, 171154d | benchmark/results/gb2tbl_stress_2026-10-30/ |
| 2026-11-04/05 | cody | table2asn 失败被当作 0 ERROR | 本地复现 | table2asn 1.28 | 退出码不可用；Critical 未读取 | 需 `.val`；Critical 计入 | Mitom_review/FIX_STATUS_2026-11-05.md |

## 8. 稳健性 / 安装 / 审查复现 / Robustness, install, reviews

| date | machine | question | data | tools & versions | key result | conclusion / action (commit) | where recorded |
|---|---|---|---|---|---|---|---|
| 2026-10-20 | ? | 安装难度 | Docker Linux、macOS | install | 表已记录 | f1c4201 | docs/benchmark_metazoa_2026-10.md |
| 2026-10-30 | ? | 线程库 | 4 记录 | MitoM | 47 min→34 s（128 核） | 738f84d | CHANGELOG.md |
| 2026-10-31 | ? | 每命令模糊测试 | 19 命令 × 24 种坏输入 | MitoM | 无 traceback，exit 2 | 7b0522d, 79b4b8d | docs/METHODS.md |
| 2026-10-31 | ? | 坐标 / 旋转不变 | 5 记录 × 12 旋转 × 反向互补；55 记录 | MitoM | 特征序列不变 | 已修 | docs/METHODS.md |
| 2026-10-31 | ? | 续跑 / 锁 | 重跑、Ctrl-C、kill -9、并发 | MitoM | 结果相同；同目录第二次运行被拒 | `.mitom.lock` | docs/METHODS.md |
| 2026-10-31 | ? | 每记录线程数 | 148 记录 | MitoM | 57 min vs 约 4 min | 8462e88 | CHANGELOG.md |
| 2026-11-04/05 | cody | 模型审查问题能否复现 | 1912f7b 最小复现 + 真实数据 | MitoM 1.0.0 | R2–R9 复现（R1 不成立）；`curate -i tests/data` ERROR 9→0 | fix/review-2026-11-05 | Mitom_review/VERIFICATION_REPORT.md |
| 2026-11-05 | ? | Cursor 复审 | 1912f7b→eaa8ef5 | 审查 | F1–F5 | a04228b | Mitom_review/FIX_REVIEW_2026-11-05.md |
| 2026-11-06 | cody | Codex 第二意见 | 未提交修复 | Codex | 3 项（import 删输入、线状 tRNA 跨端、超时算成功） | 已修 fc7ed54 | Mitom_review/SECOND_OPINION_2026-11-06.md |
| 2026-11-06 | cody | 全新用户安装 | 空 HOME，PATH=/usr/bin:/bin | install.sh、doctor、`mitom run` NC_000933 | exit 0，1 min 40 s；发现旧链接、scikit-learn 不符 | 87ebf24, 17947a4 | Mitom_review/FIX_STATUS_2026-11-06_install.md |
| 2026-11-08 | cody | MITOS2 环境 | 本机 env `mitos2`（2.1.9）、新装 `mitos`（2.1.10） | runmitos | 2.1.9 只有 runmitos.py；不激活环境直接运行缺 cmsearch | 两种名字都认（766653c） | docs/MITOS2.md §6 |
| 2026-11-08 | cody | 审查复现：`--mitos_args` 粘连 / 缩写选项；2.1.9 发现；各子命令读 `.gb.gz` 目录 | 合成参数；本机 mitos2（2.1.9）环境；只含 NC_000931.gb.gz 的目录 × 11 个子命令 | MitoM 766653c→本次 | `-c2`、`--cod 2` 能绕过（已拒）；doctor 找不到 runmitos.py（已修）；只有 gb2tbl 不认 .gb.gz（已修），其余 10 个正常 | Mitom_review/FIX_STATUS_2026-11-08.md | tests/test_review_2026_11_08.py |
| 2026-11-08 | cody | 结构重构 1–6 步前后输出是否相同 | scripts/golden.py 6 例（NC_000931 GenBank、NC_006902 内含子、NC_001636+NC_008164 两门、NC_010217 纯 FASTA、KY753833 MitoZ import、OW387151 curate）；NC_000931 run.log；3 条 rDNA contig `nuc annotate` | MitoM 9fffd02→bf80dd9，mitotools 环境 | 每步后 6 例文本输出逐字相同；run.log / 终端仅耗时不同；nuc 输出相同；pytest 796→807 通过 | 重构不改结果；以后改结构先跑 `scripts/golden.py`（refactor/structure） | CHANGELOG.md Unreleased；tests/data/golden/ |
| 2026-11-08 | cody | 第二次拆分（mito_curate → 4 个规则模块，mito_intron_ref → junction / models；删除旧名转发；setup.py）前后输出是否相同 | scripts/golden.py 6 例（同上一行） | MitoM f606ccf→本提交，mitotools 环境 | 6 例文本输出逐字相同；pytest 813 通过；wheel 含全部模块；pip install -e . 可用 | 拆分不改结果 | CHANGELOG.md Unreleased |
| 2026-11-09 | cody | `mitom annotate` 改用 run 的阶段函数后结果是否相同 | NC_000931 `annotate`（新 golden 例 annotate_gb）；其余 6 例 | MitoM 4b79cca→本提交，mitotools 环境 | all.final.gbf、05.curated、06.ncbi_submission 逐字相同；新增 original_annotation/、clade_decision.tsv、metadata_transfer_log.tsv、body_source 列，日志改为 run.log；7 例 golden 相同；pytest 813 通过 | annotate 与 run 只剩一份后处理代码 | CHANGELOG.md Unreleased；tests/data/golden/annotate_gb |
| 2026-11-09 | cody | GitHub issue #6 #7 #8 #9 能否复现、修复后结果 | NCBI MZ052209.1、OZ287881.1、OZ346613.1（tests/data/issues，golden 例 run_issues）；NC_001636+NC_008164 两门（新布局） | MitoM 8e1a9f6→本提交，mitotools 环境 | #7 当前代码不复现（1 个 trnP 7832..7894，22 tRNA）；#8 两条 linear 组装记录改判环状、cox1 起始、source 1..L、OZ287881 NCR1 6,121 bp 连成一段，NCBI ERROR 0；#6 日志路径指向真实报告；#9 新布局顶层 summary + README + 5 个目录，genes.fna 45 条、proteins.faa 39 条；classic 布局 7 例 golden 除 run_params、intron_fix.log 路径外逐字相同；pytest 816 通过 | 四个 issue 已修复 / 核实 | CHANGELOG.md Unreleased；tests/data/golden/run_issues |
| 2026-11-09 | cody | 结构审阅后的小改（start_run 按阶段表计步、去掉 hasattr；annotate 帮助；文档数字只写一处）是否影响测试 | 全部单元测试（不跑 golden：两条命令执行的阶段与终端行未变） | MitoM 2571869→本提交，mitotools 环境 | pytest 817 通过（新增 1 个阶段表测试），2 跳过 | 无行为变化 | CHANGELOG.md Unreleased |
| 2026-11-09 | cody | 校正结果改用显式标志决定是否生成提交文件、日志标签改名后，输出是否不变 | 8 个 golden 例；全部单元测试 | MitoM cb9f111→本提交，mitotools 环境 | golden 8 例逐字相同；pytest 817 通过，2 跳过 | 无行为变化 | CHANGELOG.md Unreleased |
| 2026-11-09 | cody | 环节动物线粒体内含子注释成功率（含 NCBI 无注释记录） | NCBI 环节动物 8–100 kb 全部 1,242 条：有内含子注释的 18 条（15 条在金标准、3 条为其副本）；无任何注释的 335 条 | 金标准留一结果（601cc54，intron_algo_2026-11-09）；`mitom run` 默认（mito_bench_run mitom_native）601cc54+ORF 步骤，mitotools；tblastn refseq89m 72 个环节动物蛋白 | 有注释：16 条 24 个内含子 exact 24/24，多报 1（KT989324 NCBI 漏注）；无内含子 RefSeq 环节动物 37 条多报 0；无注释 335 条：12 条写出 20 个内含子（燐虫科 4 条 11 个、小头虫科 3 条 3 个），剪接后蛋白 20/20 完整；燐虫科/小头虫科 23 条无漏检迹象；6 处可能漏检（5 个整基因未注出、Helobdella nad2） | 环节动物内含子可用；分化基因识别是缺口 | benchmark/results/intron_annelid_2026-11-09/README.md |
| 2026-11-09 | cody | 未知功能 ORF（mito_orf）能否注出、是否乱报，与 MITOS2 比较 | 才女虫科 NCBI 46 条（≥10 CDS）；无内含子 RefSeq 362 条（26 门）；MitoM 注释的无注释环节动物 335 条；MITOS2 vs MitoM 3 条 Polydorini | mito_orf.py（本提交），ORF 库 16,073 条 RefSeq；MITOS2 refseq89m；mitotools | 才女虫科：Polydorini ~355 密码子 ORF 4/4 写出，另 8 个保守 ORF 读进 ND5/trnK 给 WARNING；362 条只写 1 个（Aurelia orf340；第一版 7 个跨门假同源为串联重复，修正后 0）；335 条写 16、保守未写 9；MITOS2 0 个；单条 NC_061377 无近缘 RefSeq 时 WARNING，--orf_ref 后写出；golden NC_010217 补上海绵 atp9 | --orf yes 设为默认 | benchmark/results/orf_2026-11-09/README.md |
| 2026-11-10 | cody | ORF 检出率（已知 ORF 真值）及按失败例改进后的效果 | 260 条 RefSeq 去掉 NCBI 注释的 342 个非标准 CDS（MutS 125、hypothetical 123、atp9 57、归巢内切酶 10），逐条留一；同上 362 条无内含子、才女虫科 46、无注释环节动物 335 | mito_orf.py（本提交：内含子开放、多同源替换下限、三位密码子模式、Fisher 合并），ORF 库重建 96,187 ORF；长度基线（ORFfinder 式）；mitotools | 写出 175→228/342，WARNING 47，漏 67（多为无同源的特有 hypothetical）；atp9 57/57、MutS 123/125、归巢内切酶 6/10；只看长度 ≥75 密码子召回 326 但 362 条中写出 357 个，MitoM 写出 5 个（均有同源与编码模式证据）；密码子偏好区分力太弱未采用 | 新判定设为默认；剩余缺口记 issue #10 | benchmark/results/orf_2026-11-09/README.md |
| 2026-11-10 | cody | 按文献/GitHub 调研改进 ORF：tRNA 边界结束、Pfam 家族、RNAcode 是否有用 | 同上 342 个已知 ORF、362 条无内含子、才女虫科 46、无注释环节动物 335；Selifanova 2023 的 *Polydora* ORFan；OQ078742 `mitom run` | mito_orf.py（本提交），Pfam 14 个 HMM（hmmsearch --cut_ga，HMMER 3.4），RNAcode 0.3.1（rnacode_test 环境，线粒体密码子换标准同义）；mitotools | 写出 228→233/342，WARNING 45；MutS 123 个、归巢内切酶 8、polB 3 由 Pfam 命名；论文三个 ORFan 在 4 条记录全部写出（单条运行也写出）；362 条写出 7（Aurelia polB 由 Pfam 证实）；RNAcode 与两项检验合用无增益，未采用；table2asn 对 3′ 不完整 CDS 只给 WARNING | 新规则默认；RNAcode 不引入 | benchmark/results/orf_2026-11-09/README.md；rnacode/README.md |
| 2026-11-10 | cody | 结构审阅（fb526c3）后：ORF/NCR 独立成第 8 步、ORF 改子进程、annotate 旋转 golden | 8 个原 golden 例；新例 annotate_rotated（NC_000931 起点挪 5,000 nt 的 FASTA）；全部单元测试 | MitoM fb526c3→本提交，mitotools | 原 8 例输出逐字相同；annotate_rotated 转回 cox1 起点（COX2 1599、ATP8 3634，与原记录一致）；pytest 全部通过 | 结构修正，无输出变化 | CHANGELOG.md Unreleased；tests/data/golden/annotate_rotated |
| 2026-11-10 | cody | 漏注蛋白基因（atp8 等）能否用近缘谱 HMM 补回、会不会乱报 | 189 条记录（每门 ≤15 + 86 条有 atp8 的缨鳃虫目）逐基因去掉再找（留一属）；200 条无 atp8 注释的扁形动物/线虫；335 条无注释环节动物（MitoM 注释） | mito_pcg_hmm.py（本提交），参考库 RefSeq 16,073 + NCBI 环节动物 1,242 条（140,824 蛋白），MAFFT 7.526、HMMER 3.4；mitotools | 留一属 2,412/2,448 位置正确（atp8 187/189，atp6 186/188）；阴性 5/200（1 个丝虫为误报，4 个自由生活扁虫很可能是漏注的真 atp8）；335 条中缺 atp8 的 38 条找到 29（龙介虫 19/26），同属基因排列一致（Hydroides 都在 ND4L 前）；golden 输出不变 | --pcg_hmm yes 默认 | benchmark/results/pcg_hmm_2026-11-10/README.md |
| 2026-11-10 | cody | ORF 串联重复屏蔽改为阵列 ≥30 nt、跨原点起点修复后的 ORF 结果 | 同 ORF 四套（342 已知 ORF、362 无内含子、才女虫科 46、335 环节动物） | mito_orf.py（本提交）；mitotools | 已知 ORF 写出 233→232；362 条写出 7→8；才女虫科 12 写出不变；335 条 OZ210882.1 原先使 ORF 步骤崩溃，修复后写出 34 | 修复 | benchmark/results/orf_2026-11-09/README.md；CHANGELOG |
| 2026-11-10 | cody | atp8/atp6 按真正的门测找回率；无注释环节动物中未找回的原因 | 25 个类群（每门 ≤15 条注释了该基因的记录，缨鳃虫目单列）415 次搜索，留一属；335 条环节动物中 atp8 未找到 9、atp6 6 条逐条核查 | mito_pcg_hmm.py de0bbc0；mitotools | atp8 187/196，atp6 218/219；漏的多为近缘普遍丢失 atp8 而自身保留的扁虫/线虫（先验的代价）；未找回的 atp8 中 3 条、atp6 中 5 条为不完整/组装问题记录（跳过），6 条龙介虫（Protula 等 4 属）只有弱命中、该属无参考；去掉不完整记录 atp8 29/35 | 无代码变化 | benchmark/results/pcg_hmm_2026-11-10/README.md |
| 2026-11-10 | cody | 基因组计划（DToL/ERGA/VGP，HiFi）完整线粒体基因组上 MitoM 能否注出全部 13 个蛋白基因 | NCBI 6,309 条（5,560 条无注释）中每门 ≤25 条成环 ≥12 kb：221 条、15 门 | mitom run 留一（mitom_lo，核心去本记录、pcg_hmm 去同属），MitoM de0bbc0+；mitotools | 13 个全 180/221；节肢/脊索/棘皮/苔藓/纽形等 87 条全齐；atp8 以外 12 个基因在 218 条中只缺 1 个（海绵 nad4l）；缺 atp8 主要为已知丢失（线虫 22、扁虫 6）与分化类群（龙介虫 2、线形 4、无腔 2、软体 2）；HMM 补回 atp8 13、atp6 1；项目有注释的 15 条完全一致；缺 atp8 以外基因的 4 条逐条核对：海绵 Lycopodina 为 3 条染色体的多分体（nad4l 在 291 bp 小环）、线虫 Allodiplogaster 同一组装两个环互补、Kudoa 与同属同样丢失 5 个基因——MitoM 均正确；星虫 Themiste 未定 | 无代码变化（mito_bench_run 增加按属排除） | benchmark/results/ebp_mito_2026-11-10/README.md |
| 2026-11-10 | cody | 多分体线粒体基因组识别（mito_multipart） | 基因组计划 6,309 条中 70 个组装有 ≥2 条线粒体记录；长度不同的 16 组 37 条用 mitom run 注释（mitom_native） | mito_multipart.py（本提交）；mitotools | 合起来才完整 6 组（蠼螋、千足虫、线虫两个环；凤蝶 3 种两段）；一条已全 6 组；拷贝 1 组；合起来也不完整 2 组（Lagis 片段、海绵 Mycale 缺环，tblastn 确认）；write_bare 丢 DBLINK 已修 | 第 9 步默认运行 | benchmark/results/ebp_mito_2026-11-10/README.md |
| 2026-11-10 | cody | 转录组能否为线粒体注释提供证据（ORF 表达、3′ 端） | Polydora cf. ciliata RNA-seq SRR26830859（2.78 亿对 61 bp，随机引物；organise-sra 下载）比对到 MitoM 注释的 OQ078742 | mito_rna.py（本提交），minimap2 2.28 -ax sr；mitotools | 690,537 reads（论文 672,141）；ORFan-544/427 相对深度 0.41/0.37（ND4L 0.23、ND2 0.20），ORFan-1060 0.16 且只有 C 端转录（五段深度 8,1,8,91,63），ND3 后区段 0.01——复现论文三点结论；文库判为反向链特异；polyA reads 仅 10 条，3′ 端证据不可用（如实提示） | mitom rna 新增 | benchmark/results/rnaseq_2026-11-10/README.md |
| 2026-11-10 | cody | 转录组比对到线粒体基因组用哪个软件（poly(A) 软剪切、无剪接、个体分化、跨原点） | 模拟：3 个环节动物线粒体基因组 × 分化 0/1/3 % × 2×61/2×150，各 20 万对，带 poly(A) 与前体片段；真实：*Erpobdella* SRR10997416、*Polydora* SRR26830859 前 400 万对 | minimap2 2.31（sr / splice:sr），bwa 0.7.19，bwa-mem2，bowtie2 2.5.5（端到端 / local），HISAT2 2.2.3（± 剪接），STAR 2.7.10b（± 剪接），BBMap 40.05，Subread，segemehl；env mapbench | bwa mem 位置正确 98.2 %（3 %、61 bp 时 96.7 %），poly(A) 找回 0.96；旧默认 minimap2 -ax sr 93.5 %（76.3 %）；*Erpobdella*（差 3.5 %）蛋白基因深度 bwa 289 对 minimap2 80；内部引发位点规则（前 10 nt 连续 A ≥ 7 或 A ≥ 9：1,047 个真 3′ 端中 1 个） | mito_rna 默认 bwa mem（--mapper），internal_priming 列（本提交） | benchmark/results/rna_mappers_2026-11-10/README.md |
| 2026-11-11 | cody | 环节动物全基因组项目线粒体基因组 + 同项目转录组：MitoM 注释能否被转录组确认 / 修正 | ncbi-genome-inventory Annelida 的关联配对 54 对（同一标本 30、同项目不同标本 22、未知 2）；reads：deyuan112 前 800 万对 27、NCBI 整个 run 27 | mitom run（Annelida-segmented-worms），mito_rna bwa mem + --refine（本提交）；organise-sra locate / download（direct） | 物种核对：2 对差 > 2 %（两种蚯蚓，标记）；52 对中蛋白基因 3′ 端被 poly(A) 确认 510/673；修正 2 个（W. pigra ND6 177→152、Sclerolinum ND2 339→331，均等于同科中位长度）；Travisia COX3 误截促成“已确认末端不移动” | --refine 规则（本提交）；refined_gb/ 52 条 | benchmark/results/annelid_rna_2026-11-10/README.md |
| 2026-11-11 | cody | 同上 52 对的异常值分析（基因深度 robust z、未确认末端、孤立 poly(A) 位点、rRNA 末端） | 同上 | mito_rna --refine（本提交：重复序列、3′ UTR 确认、无终止密码子时延伸、rRNA 3′ 端） | 蛋白基因 3′ 端：末端 464、3′ UTR 15、修正 4、重复 5、无证据 185（ND4L 50）/673；rRNA 确认 88、修正 1 /104；新修正 Sabella ND1 300→305、ND3 +TA、rrnL −11 nt；Themiste ATP8/ND6 序列重复。更正：上一行 510/673 含 rRNA 注释，实为 464 | 本提交 | benchmark/results/annelid_rna_2026-11-10/README.md |
| 2026-11-11 | cody | 独立模块 mito_rnaseq_annot.py（学到的转录组规则做成算法）在 54 对上：学习集 23 对 / 测试集 31 对 | 同上（ncbi-genome-inventory 配对） | mito_rnaseq_annot.py（本提交），bwa mem | 内含子确认 4、移码拼接标记 3、未注释连接噪声 0（规则：≥ 20 % 侧翼深度）；重叠终止修正 8 / 8，3′ 端修正 3 / 1；所有可核对修正与同科中位长度差 ≤ 2 aa（10/10、8/8）；同标本蛋白基因内 ≥ 80 % RNA 插入缺失 12 处待 DNA 核实 | 本提交 | benchmark/results/annelid_rna_2026-11-10/README.md |
| 2026-11-11 | cody | MITOS2 vs MitoM vs MitoM + 转录组（mito_rnaseq_annot）在同 52 条环节动物基因组项目线粒体上 | 同上（学习 23 / 测试 29） | MITOS2 2.1.9（mito_bench_run，refseq89m）；MitoM（mitom_out）；mito_rnaseq_annot | 蛋白长度与同科中位数 ≤ 2 aa：测试集 175/292、231/295、233/295；3′ 端 poly(A) 支持 183、252、254；已证实内含子精确 0/4、4/4、4/4 | 本提交 | benchmark/results/annelid_rna_2026-11-10/README.md |
| 2026-11-11 | cody | MITOS2 经 MitoM 修正（mitos2_cur）/ 全部步骤（mitos2_full，本提交新增工具）后是否等于 MitoM | 同上 52 条 | mito_bench_run.py mitos2_cur、mitos2_full | 测试集蛋白 ≤ 2 aa：175 → 193 → 193，MitoM 231；3′ 端 poly(A)：183 → 212 → 212，MitoM 252；内含子 0 → 0 → 1，MitoM 4/4；差距在核心的 ND2、COX1、ND4、ND5 边界 | 修正补不回核心差距 | benchmark/results/annelid_rna_2026-11-10/README.md |
| 2026-11-11 | cody | MitoM 的优势是算法还是数据：只用 MITOS2 的数据（refseq89m）再跑；剩余偏差的来源 | 同上 52 条；Hermodice OZ026631 转录组共识校正 | mitom run --pcg_hmm no --orf no --rrna_hmm no，MITOM_START_MODEL=no，MITOM_EXCLUDE_LEVEL=all（本提交） | 测试集：同数据 MitoM 230/294、252 与完整 231/295、252 几乎相同（MITOS2 175/292、183），内含子 4/4 从头预测；偏差 118 中 99 个为近缘长度不一；7 个严重截短落在 ≥ 80 % RNA 插入缺失处，Hermodice 校正 10 处后 COX3 135→259、CYTB 0→366、ND4 353→450、ND6 151→158 | 优势来自算法；组装碱基错误为主要剩余误差 | benchmark/results/annelid_rna_2026-11-10/README.md |
| 2026-11-12 | cody | 转录组与基因组不一致的插入缺失：组装错误还是 RNA 层面？同一标本 HiFi 核实 | 5 条同标本记录（Hermodice、Sabellastarte、Spirobrachia、Lepidonotus、Ridgeia），HiFi 5 run 33.8 GB（organise-sra 直连） | hifi_site_check.py（本提交），minimap2 map-hifi，bwa mem | ≥ 80 % RNA 位点 21/21 为组装错误（HiFi 转录组版 412–2,233 reads vs 基因组版 0–103）；Spirobrachia ND2 部分位点为 DNA 异质性；Ridgeia 的 HiFi run 只有核内拷贝（编辑率 10 %），不可判；臭海蛹 COX3 验证为 RNA 层面 | 规则：同标本 ≥ 80 % = 组装错误；部分比例需 DNA | benchmark/results/annelid_rna_2026-11-10/README.md |
| 2026-11-12 | cody | 用 MitoHiFi 复现 Hermodice OZ026631（登录序列的 10 处组装错误是否由 MitoHiFi 产生） | HiFi ERR12760787（同标本） | MitoHiFi 3.2.1（~/mybiosofts/MitoHiFi，env mitohifi_env），hifiasm 0.19.5，MitoFinder 1.4.2；compare_mitohifi.py（本提交） | 近缘参考 PP999651 招募 0 条 reads；自身作参考：环状 15,861 bp，10/10 位点为 reads 版本，与登录序列差 18 indel；MitoHiFi 在 MitoFinder 处中断（31.7 kb 双拷贝 contig、未注释参考无 tRNA） | 登录错误不能用当前 MitoHiFi 复现 | benchmark/results/annelid_rna_2026-11-10/README.md |
| 2026-11-12 | cody | P9 修复：核心读码框错时取另一框的完整 ORF，而非移码拼接 | 67 条环节动物基因组项目线粒体（mitom run 前后对比）；鸟、龟 RefSeq 6 条（ND3 真移码） | mito_curate_pcg.other_frame_orf（本提交）；全部单元测试 861 通过 | 873 个蛋白基因中只变 1 个：Travisia ATP8 49 aa 拼接 → 53 aa 单框；龟 ND3 仍在第 174 位拼接；Sabellastarte ND3、Spirobrachia ND6（组装错误）仍为拼接 | 本提交 | knowledge/benchmark_findings.md |

## 9. 核标记 / mitom nuc

| date | machine | question | data | tools & versions | key result | conclusion / action (commit) | where recorded |
|---|---|---|---|---|---|---|---|
| 2026-10-25 | ? | rDNA 边界 vs barrnap | 378 种；71 独立；组蛋白 677+172 | mito_nuc vs barrnap 0.9 | exact 83.2 vs 1.8 %；组蛋白 97.5/97.7 % | 8a0debe | docs/METHODS.md |
| 2026-10-25/26 | ? | RefSeq rRNA 金标准 | 336 记录 558 边界 | mito_nuc | exact 85.8→96.2 % | Rfam 端点细化 | docs/METHODS.md |
| 2026-10-25 | ? | 假阳性与 table2asn | 898 rDNA；2,012 组蛋白；200 线粒体 | table2asn 1.28 | 0 ERROR | 已修 | docs/METHODS.md |
| 2026-10-26 | ? | 大组装运行时间 | SRR16188826（148 万 contig） | mito_nuc | 200 记录 24 min→3 min | 预筛 | docs/METHODS.md |
| ? | ? | 30 个 SPAdes 组装（≤1.7 Gb） | 30 组装 | mito_nuc | 6 s–2.2 min（8 线程） | 低复杂度屏蔽 | CHANGELOG.md |

## 10. 其他 / Other

| date | machine | question | data | tools & versions | key result | conclusion / action (commit) | where recorded |
|---|---|---|---|---|---|---|---|
| 2026-10-26 | ? | 蛋白基因第二拷贝 | 鹭、鸬鹚等记录 | MitoM | 多得 9 个 NCBI 有的基因、1 个 NCBI 没有 | COPY_IDENTITY 95 % | CHANGELOG.md |
| 2026-11-19 | cody | 线粒体覆盖度怎么算最快最准：using_bwa_draw_coverage（bwa+samtools depth）/ MitoZ / GetOrganelle（bowtie2）/ 流式计数；环状原点是否丢深度 | 模拟：NC_001636.1 15.5 kb，wgsim 2×150，旋转半圈使跨原点 reads 有真值；20 万 / 100 万对，SNP 0 / 5 %，1 / 8 线程 | bwa 0.7.19，bwa-mem2，minimap2 2.31，bowtie2 2.5.5，samtools（env mapbench）；mito_coverage.py | 原点深度：samtools depth 无环处理 0.90–0.94 倍真值，延长参考后 1.00；流式计数比 BAM+sort+depth 快约 25–30 %；5 % 差异下 bwa r 0.9985 > minimap2 0.9954 > bowtie2 0.9938（minimap2 最快 2.5 s，bowtie2 最慢 3–5 倍）；pyasm 是字节码库，未比 | `mitom coverage`（默认 bwa mem + 流式计数 + 环状延长折回；--mapper minimap2 可选）；仅模拟数据，真实数据未测 | benchmark/results/coverage_2026-11-19/README.md |
| 2026-11-19 | cody | tRNA 判别器（caller / bounds）改为 JSON 树 + numpy 推理（无 scikit-learn、无 pickle）后，概率是否与 scikit-learn 相同；整个流程输出是否不变 | models/trna_ml/caller_train.tsv.gz（82,448 行）、bounds_train.tsv.gz（129,310 行）；golden 全部案例 | scikit-learn 1.7.2；mito_gbdt.py；tests/test_trna_json.py；scripts/golden.py --jobs 2；pytest 856 通过 | decision_function 最大差 0.0，predict_proba 最大差 2.2e-16（两张表所有行；含缺失值构造），预测类别 0 差异；golden 仅提示 reports/intron_fix.log 为参考中没有的新文件（与本改动无关），其余相同；屏蔽 sklearn 后推理可用；预发布 JSON 随仓库提供，首次运行不再训练 | exp/trna-json-2026-11-19 | tests/test_trna_json.py |
| 2026-11-20 | cody | 外部审核报告（mitom (1).md）的代码缺陷逐条核对 + 7 处修复是否改变真实结果 | annelida_ref_curate_250320.gb（170 条 RefSeq）、GCA_001317505.1 gbff；pytest 881 通过；golden 通过 | mitotools 环境 | complete_claim：170 条中 4 条旧逻辑误用第 1 条声明；mito_check 前后 170 条结果相同；.gbff 旧代码崩溃、新读 75 条；元数据 O(N²) 已在 9ba101d 修好（二次调用 0.0001 s）；code 14/33 无本地数据，NCBI 几乎不用，未改 | 7 处修复（placed[6]、complete_claim 逐记录、rRNA fill、genes.fna 去 pseudo、.gbff、缩写终止子、intron 工具检查） | benchmark/results/review_check_2026-11-20/README.md |
| 2026-11-20 | cody | issue #12：`mitom run --infq` 的 coverage.tsv / depth.tsv.gz / summary 深度列是否正确 | NC_000931 旋转半圈后 wgsim 60000 对 2×150（跨原点 reads）；`mitom run --infq` | minimap2，mitotools | 平均深度 1152.54（理论 60000×300/15619=1152.5）、mapped_reads 120000、breadth_1x 100 %、min 5；用时 56 s | 通过；issue #11 镜像 md5 bb6325e2… 与 REF_MD5 一致 | tests/test_mito_coverage.py、tests/test_organize_outputs.py |
| 2026-11-20 | cody | 强内含子模型能否胜过核心的弱冲突 CDS（issue #15，MW592986.1 atp6）；会不会改变其他记录 | intron_ncbi_2026-11-09 的 134 条 NCBI 记录，内含子步骤旧/新代码 | intron_only.py，mitotools | 134 条中 1 条报告不同（MW592986.1 atp6 low→high，位置同 NCBI 16489..17094），其余 133 条相同 | 通过；mito_intron.outweighed_conflicts | benchmark/results/weak_conflict_2026-11-20/README.md |

## 未测 / 未定 / Not yet tested or open

| item | source |
|---|---|
| MITOS2 导入后 curate：完整终止越过下游 tRNA ≤ 20 nt 被保留（按核心调定），对 MITOS2 结果是否截成 T/TA 未在 24 条集上测 | mito_curate.py OVERLAP_FALLBACK；本表 §2 2026-11-08 |
| 软体动物数据上的 curate 规则未测 | docs/METHODS.md（curation） |
| rRNA：13 条环节 RefSeq 端点自相矛盾；牡蛎分裂 16S、被囊 rrnL、刺胞 R1 填充 | knowledge/benchmark_findings.md; knowledge/experiments/rrna_ceiling/README.md |
| 起始密码子：NCBI 无蛋白证据选下游起始；双壳类 N 端与 12S；线虫 atp8；线虫/棘头 MiTFi 漏 tRNA（ARWEN 43 对/20 错，未采用） | knowledge/benchmark_findings.md |
| NC_027103.1（Rhodactis 20 kb nad5）、NC_061743.1 所有工具都错 | knowledge/benchmark_findings.md |
| P11 对已 curate 记录再 curate 会改变起点（2 例，待定） | docs/METHODS.md |
| 混合门嵌合体无警告；小门的门判断常错（--clade auto） | docs/METHODS.md（clade auto） |
| 拓扑：只缺非编码区/tRNA 的 contig 判为环状；myxozoa 判为线状 | benchmark/results/topology_unlabelled_2026-10-31/README.md |
| 内含子：Chaetopterus 约 300 nt 预测、Decemunciger ND4 未确证 | docs/METHODS.md; benchmark/results/intron_gold_2026-10-29/gold_README.md |
| 注释步骤超时、[organism=.] 占位待用户决定；Docker 全新安装未测 | Mitom_review/FIX_REVIEW_2026-11-05.md; Mitom_review/FIX_STATUS_2026-11-06_install.md |
| 2026-11-20 | cody | 同源投票+打乱零模型能否找回 Exonerate 漏掉的短外显子（issue #15，19 条 bryozoan 漏检内含子） | intron_ncbi_2026-11-09 的 134 条 NCBI 记录 | intron_only.py，eval_introns.py，mitotools | 全链搜索：exact 148→146、missed 22→25（新外显子边界偏 ~40 nt，并给已正确的 cox1/cox2 加假外显子）；仅无模型基因：与 head 完全相同（148/22/3） | 无增益；`--short_exons` 保持默认关闭 | benchmark/results/short_exons_2026-11-20/README.md |
| 2026-11-20 | cody | 短外显子：边界/框架拟合 + 增益守卫（--short_exons=all）（issue #15） | intron_ncbi_2026-11-09 的 134 条 NCBI 记录 | intron_only.py，eval_introns.py，mitotools | exact 148→151、missed 22→19、fp 3，无回退（3 条 cob 变 exact） | 有净增益；金标准/无内含子集未跑，默认仍关闭 | benchmark/results/short_exons_2026-11-20/README.md |
| 2026-11-20 | cody | Codex 深度审查内含子算法（gpt-6-astra）+ 修 finding 1/4（short6） | 同上 134 条 | second-opinion.sh，intron_only.py，eval_introns.py | 与 short5 完全相同（151/47/31/19，fp 3）；未修：locate 查询坐标映射、多新外显子联合拟合、nad4 供体 +42 nt、cox3 末端假阳性；3 条“漏检”实为 NCBI 的 nad3 移码缺口 | 两处防御性修复；其余待做 | ~/.cache/second-opinion/20261120-182924.M0ZA/review.md |
| 2026-11-20 | cody | 边界精度：nad5 3′端偏 −36 nt（BK069890–BK071243 系列）是算法错还是 NCBI 坐标？ | 同上 134 条；BK069890.1 nad5 | intron_only.py，eval_introns.py，homolog_scorer，mitotools | REFINE_3 24→40 无任何变化（151/47/31/19/3）；BK069890.1 nad5 拼接蛋白同源得分：预测 1186，真值位置（−36）1112，+15 起出现终止密码子 | 同源证据支持预测而非 NCBI；非算法缺陷，REFINE_3 已还原为 24，代码无改动 | 本行 |
| 2026-11-20 | cody | Codex 深度审查（核心注释 + RNA/核rDNA）的 8+4 条高危发现修复，逐条先写失败测试 | 单元测试（合成数据） | pytest 全套，second-opinion deep | 见 CHANGELOG Unreleased；全套 pytest 通过（见提交） | 本次提交 | CHANGELOG.md |
