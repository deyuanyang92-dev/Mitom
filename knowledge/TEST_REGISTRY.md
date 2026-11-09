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
