# MitoM 软件深度评估与优化建议报告

> **评估对象**：`deyuanyang92-dev/Mitom`（v0.0.1 测试版，含核心源码、打包与环境配置、基准评测体系）
> **评估维度**：代码架构、算法可靠性、性能与资源利用、后期维护与合规性、分优先级优化路线图

---

## 执行摘要

**MitoM 是一款实证根基极扎实、工程完成度远超同类科研脚本的动物线粒体基因组注释系统。** 它成功从早期的“MitoZ + MITOS2 封装脚本”演进为单 Conda 环境、自研原生核心（`mito_core.py`）的一站式工具链，在蛋白基因起止边界精度、tRNA 召回/精确率、内含子建模及 NCBI `table2asn` 零报错提交合规性上，全面领先 MitoZ、MITOS2 与 MitoFinder。

与此同时，随着功能版图迅速扩张（新增机器学习 tRNA 判别、GBDT 起始密码子模型、未知 ORF 选择检验、多基因座核标记 `mitom nuc`、RNA-seq 证据校验等），软件正处在**“研究型多脚本集合”向“长期可维护标准软件包”跨越的关键拐点**。当前最核心的风险集中在：
1. **架构层**：40+ 顶层 `.py` 平铺无包命名空间、跨步骤通过自由文本 `/note` 正则传参、多线程后调用 `fork()` 进程池的死锁隐患；
2. **算法与代码层**：部分边界条件与静默失效细节（如 `trna_ml.py` 跨原点坐标与聚类未取模、`mito_pcg_hmm.py` 死参数 `--k`、`mito_start.py` 哨兵值冲突、`mito_intron.py` 运行时联网无校验拉取外部模型）；
3. **性能层**：批量多记录 GenBank 场景存在两处隐蔽的 $O(N^2)$ 解析瓶颈，以及 tRNA 双引擎对 22 个 pan-CM 的重复扫描；
4. **维护与合规层**：多许可证声明不一致（CC0 与 GPL-3.0 混用冲突）、打包清单漏洞（Docker/Bioconda 漏排或误收文件）以及仓库体积膨胀。

| 评估维度 | 综合评级 | 核心亮点 | 主要短板 |
|---|---|---|---|
| **代码架构** | 中偏上（B+） | 单环境原生核心收敛彻底；`RUN_STAGES` + `RunContext` 流程清晰；`mito_debug.py` 可观测性极佳 | 根目录平铺 40+ 模块；步骤间靠 GenBank `/note` 正则耦合；多线程后 `fork` 风险；核基因引发范围膨胀 |
| **算法可靠性** | 高（A-） | 防泄露评测体系达到顶刊水准；多策略互补显著领先竞品；统计检验原理扎实 | 若干边界/环状取模缺陷；门级启发式规则与 ML 决策面耦合；内含子亚型模型无哈希校验联网拉取 |
| **性能表现** | 优秀（A-） | 单条中位 70 s / 峰值 182 MB 全面领先竞品；`mito_env.py` 线程风暴治理到位 | 大批量多记录 GenBank 终稿写出存在 $O(N^2)$ 重复解析；tRNA 双引擎重复扫描 pan-CM；`mito_pcg_hmm` 缺 HMM 缓存 |
| **后期维护与合规** | 中（B-） | 有 `golden.py` 回归守护、`test_structure.py` 架构门禁、`frozen/` 冻结评测与文档自动生成 | CC0 与 GPL-3.0 声明冲突且缺 GPL 文本；7 套构建/依赖配置难同步；`.dockerignore` 漏排审查目录；基准脚本残留本机硬编码路径 |

---

## 一、代码架构评估

### 1.1 架构亮点（值得长期保留的设计）

1. **彻底的“去外部重依赖”原生核心收敛**
   * 项目从早期依赖 MitoZ + MITOS2（3 个 Conda 环境、约 8 GB 磁盘）重构为单一 `mitotools` 环境的自研核心 `mito_core.py`（仅需 1 个环境 + MITOS2 `refseq89m` 参考数据），将 MitoZ/MITOS2 降级为可选 `--engine` 或 `mitom import` 导入源。这大幅降低了安装门槛，也让跨原点搜索、旁系同源屏蔽与边界精修实现了完全自主可控。
2. **清晰的阶段编排与上下文隔离**
   * `mitom_pipeline.py` 已收敛为薄入口，核心工作流拆分至 `mitom_run.py` 的 `RUN_STAGES` / `ANNOTATE_STAGES` 阶段表，通过 `RunState` 与 `RunContext`（`mitom_context.py`）在各步骤间传递状态，消除了早期大量流程级全局变量。
   * 后处理步骤（`mitom_steps.py`）统一使用 `StepResult` 接口与 `@optional_step` 装饰器实现故障降级；重型/高风险步骤（`mito_curate.py`、`mito_intron.py`、`mito_intron_ref.py`、`mito_orf.py`、`mito_pcg_hmm.py`）以独立子进程运行实现崩溃隔离。
3. **罕见优秀的可观测性与防呆工程**
   * **结构化日志与自动诊断**：`mito_common.run_tool()` 统一托管全部外部工具调用并输出固定语法（`[STEP]`、`[CMD]`、`[RECORD]`、`[ENV]`、`[DONE]`），配合 `mito_debug.py`（`mitom debug`）能直接给出失败步骤、失败记录的首条错误与一键复现命令。
   * **原子写与进程锁**：`mito_common.write_gb()` / `write_text_atomic()` 采用同目录临时文件 + `os.replace` 原子替换（并统一修复 Biopython 将无值 `/pseudo` 写成 `/pseudo=""` 的历史坑）；`mitom_prepare.lock_output()` 结合 PID 与 `/proc/<pid>/stat` 启动时间戳防止并发写坏同一输出目录，且能安全回收陈旧锁。
   * **大文件拆分有序**：`mito_curate.py` 按规则族拆为 Mixin（`mito_curate_{base,origin,pcg,rna}.py`），`mito_intron_ref.py` 拆为 `sites/junction/models`，`mito_nuc.py` 拆为 `lib/rdna/histone`，职责边界较为清晰。

### 1.2 核心架构问题与技术债

1. **顶层平铺 40+ 模块，缺少标准 Python 包命名空间**
   * 所有 `.py` 文件全部裸放在仓库根目录，混杂了领域核心（`mito_*.py`）、流程编排（`mitom_*.py`）、无前缀核心（`trna_ml.py`、`ncbi_validate.py`）与兼容垫片（`batch_mitoz.py`、`setup.py`）。
   * 负面连锁效应明显：`setup.py` 靠 `HERE.glob("*.py")` 自动发现模块；`packaging/bioconda/build.sh` 用 `cp mitom *.py "$dest/"` 拷贝文件，会把 `setup.py` 和开发遗留脚本一并打进正式发行包；`benchmark/` 内大量脚本仍需 `sys.path.insert` 甚至 `benchmark/mito_names.py` 这种 `importlib` 重定向垫片。
2. **跨阶段依赖自由文本 `/note` 正则传参（隐式强耦合）**
   * `mito_core.py` 将关键内部元数据（参考蛋白比对区间 `ref ... score ... residues a-b/n`、MiTFi 的 `E` 值与 `copy` 标记）编码进 GenBank 特征的 `/note` 文本中；随后 `mito_curate_rna.py`（`N1`）、`mito_core.py`（`end_evidence`）、`mito_intron.py`（`CORE_REF_RE`）与 `mito_intron_ref.py`（`_intron_note`）各自用正则表达式逆向解析，而 `mito_curate_rna.py` 在 `N1` 阶段又会用 `INTERNAL_NOTE` 将这些注释删掉。
   * `mitom_steps.py` 同样依赖解析子进程的 `stderr` 文本（`summary_line`、`[mito_orf] WARNING `）与 TSV 报告来回传结构化信息。这种“以人类可读文本充当内部 IPC 协议”的模式对步骤顺序与文案微调高度敏感，一旦上游日志/注释格式改动，下游规则（如 P11 起始模型、T1 弱 tRNA 清理、内含子核心参考回溯）会静默失效。
3. **多线程之后显式 `fork()` 进程池（POSIX 死锁与兼容性隐患）**
   * `mitom_steps.py`（`run_plots`）、`mito_intron.py`、`mito_intron_ref.py`、`mito_intron_sites.py` 均显式使用 `mp.get_context("fork").Pool(...)`（利用子进程继承模块级 `_JOB` 全局字典避免 pickle）。
   * 但在 `mitom run` 流程中，前序阶段（`mitom_stage.py`、`mitom_snakemake.py`）已使用 `ThreadPoolExecutor` 与后台守护线程（`watch_files`），并导入了 `numpy`/`matplotlib`。在多线程进程中调用 `fork()` 极易触发锁状态死锁或 macOS 崩溃（Python 3.12+ 已对此抛出 `DeprecationWarning`）。
4. **同类机器学习子系统的推理架构割裂**
   * `mito_start.py` 采用了极佳的工程实践：训练期将 `HistGradientBoostingClassifier` 导出为纯 JSON 树，推理期仅用 `numpy` 计算（零 `pickle`、零 `scikit-learn` 运行时依赖）。
   * 反观 `trna_ml.py` 仍采用“分发压缩训练表 $\rightarrow$ 首次运行现场训练 $\rightarrow$ `pickle` 缓存”模式，迫使整个软件在 `env/mitotools.yml` 和 Bioconda 中强锁 `scikit-learn==1.7.2`。
5. **功能边界扩张（Scope Creep）与外部包脆弱 Monkey-Patch**
   * `mito_nuc*.py`（面向 GB 级 WGS 组装的核 rDNA / 组蛋白注释）与 `mito_rnaseq_annot.py` 与线粒体注释核心耦合在同一顶层目录，导致依赖树、CLI 和文档持续膨胀。
   * `scripts/runmitos_fast.py` 通过在运行时 Monkey-Patch `mitfi.util.execute`、`mitfi.cmsearch` 和 `mitfi.remove_empty`（先 dry-pass 记录命令并删除新建空文件，再并发回放）强绑定 `mitos=2.1.10` 内部实现，升级极为脆弱。

---

## 二、算法可靠性评估

### 2.1 方法论核心优势

1. **科研级防数据泄露评测体系**：`trna_ml` 按科级哈希划分并排除 `refseq89m` 已有记录；`train_start_model.py` 基于留一核心输出做科级 5 折 `GroupKFold`；`frozen_stats.py` 采用按科聚类的 bootstrap 估计 95% 置信区间；`nuc_cm_gold.py` 使用独立的 Rfam 协方差模型后验概率构建金标准，规避了 NCBI 历史标注不一致的标签噪声。
2. **多源互补的算法设计**：
   * **蛋白基因**：`blastx`（含弱命中敏感重试）+ `exonerate protein2dna` + `coding_frame` 移码纠偏 + 旁系同源屏蔽 + `ADJACENT` 基因序兜底（`atp8`/`nad4l`）+ `mito_pcg_hmm.py` 近缘谱 HMM（含基于近缘基因存在率的贝叶斯先验几率修正与 Bonferroni $\times 3$ 校正）。
   * **tRNA**：MitoM-ML（pan CM + 220 个 clade CM + ARWEN/ARAGORN + GBDT 判别与边界模型）+ MiTFi 二次校验（`p < 0.9` 与次级拷贝需确认）+ MiTFi 缺失家族补救 + 脊椎动物末端切换，在 565 条未见科测试集上实现 98.75% 召回 / 99.06% 精确 / 88.70% 两端精确，全面超越 MITOS2。
   * **内含子**：MFannot 风格从头建模（`mito_intron.py`，支持跨原点旋转建模、微外显子插入、同源连续性与剪接增益检验）+ 近缘位点迁移（`mito_intron_ref.py`，结合 60 nt 接点探针与金标准规范滑移），在 193 条金标准记录上默认 309/309 精确、留一 299/309 精确且 76 条无内含子负对照 0 假阳性。
   * **未知 ORF**：`mito_orf.py` 从碱基组成零模型 $E_{\text{len}}$、多同源序列最小生成树（MST）替换下界 $p_{\text{open}}$、局部核苷酸比对密码子 1/2/3 位点变异模式检验、串联重复屏蔽到 Pfam `--cut_ga` 域注释，显著优于单纯长度阈值法。

### 2.2 已定位的具体算法隐患与代码缺陷

| 优先级 | 位置 | 缺陷 / 隐患描述 | 潜在影响 |
|---|---|---|---|
| **高** | `trna_ml.py`: `end_pairs()`、`_loci()`、`call()` | 1. `end_pairs()` 中 `sub = g2[s - 1:e] if e <= 2 * L and s >= 1 else ""`，当 `s = 1` 向左偏移 `-1` 变为 `s = 0` 时直接返回空串并被丢弃，丢失环状基因组原点左侧 `-1 nt` 合法边界。<br>2. `_loci()` 与 `call()` 中的 `_ov(a, b)` 未做环状取模，且 `drop_wrap_duplicates()` 仅对同模型（`x["model"] == h["model"]`）去重；若 `pan` 模型命中 `1..68` 而 `clade` 模型命中 `L..L+67`，两者既不去重也不聚类，可能产生跨原点重复调用。<br>3. `scan()` 对线性序列仍拼接 `seq[:WRAP]`，虽最终丢弃越界调用，但已影响 `n_strong_genome` 等全局上下文特征。 | 跨原点 tRNA 边界偏移或产生重复假阳性；线性短序列上下文特征轻微失真 |
| **高** | `mito_intron.py`: `subtype_cms()` | 首次运行且无本地缓存时，直接联网请求 GitHub API（`LaraSellesVidal/Group1IntronDatabase` `main` 分支）下载 `.cm` 模型，**无版本固定、无 SHA256 校验、无开源许可证确认**。离线节点会阻塞至超时（30 s）后靠 `.failed` 缓和 24 h。 | 供应链安全风险、跨时间不可复现风险、离线集群首次卡顿 |
| **中** | `mito_start.py`: `candidates()`、`StartModel.best()`；`train_start_model.py` | 1. `candidates()` 中 `c > b - 30` 在极短/截断 CDS（$b - a < 30$）时会连原始起点 `a` 一并过滤，返回空列表；独立调用 `StartModel.best([])` 会抛出 `ValueError`。<br>2. 用 `MISSING = -999.0` 作为缺失哨兵，且 `export()` 丢弃了 `HistGradientBoostingClassifier` 原生 `missing_go_to_left` 属性；当 `up_gap`/`core_off` 为较大负数时会与缺失值碰撞。<br>3. `train_start_model.py` 依赖 sklearn 私有属性 `model._predictors`，且 `start_model.json` 未校验上游 P1 规则版本，P1 改动会导致 `is_p1`/`cand` 特征分布静默漂移。 | API 调用边界崩溃；大重叠基因误判为缺失特征；未来升级 sklearn 或调整 P1 规则后模型静默退化 |
| **中** | `mito_curate.py`: `Curator.run()`；`mito_taxa.py` | 1. `Curator.run()` 在单轮 `for f in CDS` 循环内串联执行 `origin_cut -> frameshift -> frame_repair -> stop_trim -> pcg -> across_origin`，而 `pcg()` 会实时读取相邻 CDS 坐标；基因 $i$ 判断下游边界时，基因 $i+1$ 的移码/修剪尚未执行，存在隐蔽的遍历顺序依赖。<br>2. `mito_taxa.py` 中硬编码了多条门级特例（环节动物 `ND5,ND4L`、脊索动物 `abbreviated_overlap`、线虫/节肢 `no_upstream_atg` 等），全部依赖字符串前缀匹配，既拟合了 NCBI 历史标注惯例（非严格生物学真相），又与下游 P11 起始模型共享决策面，增加了排查难度。 | 相邻重叠基因边界受遍历顺序微弱影响；手工规则与 ML 模型耦合难调优 |
| **中** | `mito_pcg_hmm.py`: `search_record()` | CLI 定义了 `--k` 参数并传入 `search_record(..., k=a.k)`，但函数体内完全未使用 `k`（硬编码 `max(SIZES)` 即 100，且调用 `prior()` 未传 `k`）。 | 用户通过 `--k` 调参完全无效（死参数 Bug） |
| **中** | `mito_orf.py`: `third_position_p()` | 对密码子 1/2/3 位点变异计数 `d` 直接采用 2 自由度 $\chi^2$ 渐近公式 $p = \exp(-X/2)$，未设最小变异总数阈值（如 $n \ge 15$）。当局部核苷酸比对仅匹配极少数变异（例如 $n = 7, d = [1, 0, 6]$）时，$\chi^2$ 渐近近似会高估显著性。 | 极低变异或近缘短片段比对下产生过度乐观的 $p_{\text{pos3}}$ |
| **低** | `mito_rnaseq_annot.py`: `Tap._read()`、`single_frame_orf()` | 1. `single_frame_orf()` 直接跳过负链（`if f.location.strand == -1: return None`），仅对正链移码基因生效。<br>2. `Tap._read()` 解析 `SA:Z:` 补充比对计算内含子跨度时未处理跨原点环状坐标，`step_introns()` 第一循环访问 `depth[a]`/`depth[b]` 未做 `% len(depth)` 保护。 | 负链错误移码基因无法检出替代单框 ORF；跨原点剪接位点可能遗漏或越界 |
| **低** | `mito_nuc.py`: `FastaIndex` | `FastaIndex` 打开底层文件描述符 `self.fd = os.open(...)`，但未实现上下文管理器（`__enter__`/`__exit__`），若 `annotate()` 中途抛异常会导致 fd 泄漏。 | 异常路径下文件句柄泄漏 |

---

## 三、性能与资源利用评估

### 3.1 实测性能优势
* **端到端资源领先**：在 24 条跨 8 门 RefSeq 基准中，MitoM 单条中位耗时 **70 s**（MitoFinder 204 s、MitoZ 249 s、MITOS2 385 s），中位峰值内存 **182 MB**、最大峰值 **270 MB**（优于 MitoFinder 的 283 MB、MitoZ 的 719 MB 与 MITOS2 的 1,929 MB）。
* **底层性能意识成熟**：
  * `mito_env.py` 在导入最早期即把 `OMP_NUM_THREADS`、`OPENBLAS_NUM_THREADS` 等默认压至 `1`，彻底根治多进程并发时 NumPy/OpenBLAS/sklearn 抢核引发的负载爆炸（128 核实测从 47 min 降至 34 s）；
  * `mito_common.read_fasta()` 绕过 `Bio.SeqIO` 启动期对 `numpy` 的隐式导入（单核心进程常驻内存从 40 MB 降至 18 MB）；
  * Infernal 内存矩阵限 128 MB、内含子 `cmsearch` 用 `--cpu 0` 避免多线程复制 DP 矩阵；
  * `mito_nuc.py` 通过 `FastaIndex`（`os.pread`）、并行 `makeblastdb` + `blastdb_aliastool`、低复杂度过滤（`-dust yes`）及 `MAX_*_REGIONS` 上限，将 1.7 Gb 级全基因组组装的核基因提取压到 2.2 min 以内。

### 3.2 关键性能瓶颈与优化空间

1. **【严重】批量多记录 GenBank 终稿写出存在两处 $O(N^2)$ 解析瓶颈**
   * 在 `mitom_run.py` 的 `stage_write_final()` 中，每个样本都串行调用 `finalize_genbank_output()`。当输入为一个包含 $N$ 条记录的多记录 GenBank 文件（如 `mitom fetch` 产出的 600+ 条 `Annelida.gb`）时，所有样本的 `metadata_src_gb` 指向**同一个大文件**。
   * 但在 `mitom_output.py` 的 `finalize_genbank_output()` 内：
     * Step 3 每次都调用 `transfer_gb_metadata_file()`，重新读取整个大文件并调用 `build_metadata_index(src_text)`；
     * Step 3.5 每次都调用 `build_organism_index([metadata_src_gb])`，重新用 `SeqIO.parse` 完整解析整个大文件。
   * 对 $N = 636$ 条记录，这会触发 $636 \times 636 \approx 40.4$ 万次 GenBank 记录解析！只需按 `metadata_src_gb` 路径在阶段级预构建并缓存索引（每个唯一源文件只解析 1 次），即可将 $O(N^2)$ 直接降为 $O(N)$。
2. **【高】tRNA 阶段双引擎对 22 个 pan-CM 重复扫描**
   * `mito_trna.core_trnas()` 在 `use_ml=True` 时，先由 `annotate_trnas_ml()` 调用 `trna_ml.scan()` 跑完 22 个 pan-CM + clade-CM，紧接着无条件调用 `annotate_trnas()`（MiTFi）再把 22 个 pan-CM 完整 `cmsearch` 一遍（外加敏感补搜与 JVM 启动）。
   * 实际上，MiTFi 仅在三种情况下才真正被用到：① ML 存在 `p < 0.9` 或次级拷贝需要 `trna_second_opinion`；② 属于 `Vertebrata` 需要 `mitfi_ends`；③ 存在未检出的标准家族需要 `rescue_trnas`。按需触发 MiTFi，或直接共享同一份 pan-CM `cmsearch` 输出，预计可砍掉 tRNA 阶段 25%–45% 的算力开销。
3. **【中】`mito_pcg_hmm.py` 缺乏跨记录的近缘 HMM 缓存**
   * 在批量注释同一科/属（如 200 条环节动物均缺失 `atp8`）时，`mito_pcg_hmm.py` 为每条记录都分别对 10/30/100 条近缘参考重跑 MAFFT + `hmmbuild`（600 次冗余建库），而选出的参考序列集几乎完全相同。按 `(gene, tuple(ref_accs))` 在批次临时目录内缓存 `.hmm` 文件，可大幅提速。
4. **【中】后处理多处仍为单线程串行执行**
   * `run_submission()`（`mitom_steps.py`）对每条记录串行执行 `table2asn`；`mito_curate.py` 与 `mito_pcg_hmm.py` 串行处理文件列表。在数百条记录的批处理中应复用 `--jobs` 并发度。
5. **【低】热点循环微瓶颈与首跑并发冗余**
   * `mito_multipart.py` 的 `kmers()` 在 `for i in range(len(s)):` 内层循环里逐碱基调用 `str.maketrans("ACGT", "TGCA")`，应提取为模块级常量 `_COMP`。
   * `trna_ml.classifiers()` 虽已用 `os.replace` 保证原子落盘，但未加文件锁：若用户未跑 `mitom setup` 直接启动 `--jobs 8`，8 个进程会同时从头训练 GBDT 模型，造成首跑 CPU 与内存尖峰。
   * `mito_start.py` 的 `StartModel.score()` 每次推理均重新构造列名映射字典 `cols`，应在 `__init__` 中预计算。

---

## 四、后期维护、打包与合规性评估

### 4.1 许可证与合规风险（公开发表前务必理顺）
* **仓库内部许可证声明自相矛盾**：
  * `pyproject.toml` 与 `packaging/bioconda/meta.yaml` 明确声明复合许可证 `CC0-1.0 AND GPL-3.0-only AND MIT`（因为 `mito_intron.py`、`mito_intron_rules.py` 与 `models/mfannot_proteins` 源自 GPL-3.0 的 MFannot，`models/mitfi` 源自 MIT 的 MITOS2）。
  * 但 `README.md` 徽章与声明段落、`CITATION.cff` 及根目录 `LICENSE` 文件均仅声明 `CC0-1.0`，且仓库未附带 GPL-3.0 完整许可证文本（违反 GPL-3.0 再分发必须附带许可证全文的条款）。
* **外部模型版权边界不明**：`mito_intron.py` 运行时拉取的 `Group1IntronDatabase` 在 GitHub 上无明确开源许可声明，将其纳入默认自动下载链路存在潜在版权争议。

### 4.2 仓库体积、打包与依赖治理
1. **仓库严重膨胀，科研过程产物与生产代码混杂**：完整仓库解压达 90.1 MB、1,437 个文件，大量体积来自 `benchmark/results/` 下直接提交的 `.gb`、大 `.tsv`/`.tsv.gz` 和图片，以及双轨审查目录（根目录 `Mitom_review/` 16 个文件 + `docs/reviews/`），拉慢克隆速度并干扰外部贡献者判断当前状态。
2. **Docker 与 Bioconda 打包漏洞**：
   * `.dockerignore` 排除了 `benchmark`、`knowledge`、`packaging`、`tests`，但**漏排了 `Mitom_review`、`docs/reviews`、`docs/images`**，导致生产镜像夹带历史审查文档；
   * `packaging/bioconda/build.sh` 用 `cp mitom *.py "$dest/"` 会误把 `setup.py` 和废弃兼容脚本 `batch_mitoz.py` 打入生产包。
3. **多处依赖清单漂移风险**：项目同时维护 `pyproject.toml`、`setup.py`、`install.sh`、`Dockerfile`、`env/mitotools.yml`、`env/mitos.yml`、`packaging/bioconda/meta.yaml` 共 7 处构建/依赖文件。目前 `pyproject.toml` 未锁 `scikit-learn` 版本，而 `env/mitotools.yml` 与 `meta.yaml` 强锁 `1.7.2`。
4. **基准脚本残留本机硬编码路径**：`benchmark/trna_ml/scripts/` 下的 `common.py`、`run_baselines.py`、`run_mitfi.py` 硬编码了 `/opt/micromamba/envs/invertrnax/bin`、`/root/miniconda3/envs/mitos/bin` 以及本机中文路径 `/mnt/n/yzz-分类工作流/temp_scripts/Mitoz-annotate`，影响外部可复现性。

### 4.3 维护体系的突出优点
* `scripts/golden.py`（全流程金标准回归测试）+ `scripts/make_options_doc.py`（自动从 `-h` 生成 108 KB 的 `docs/OPTIONS.md` 并由 `test_docs.py` 校验同步）+ `tests/test_structure.py`（强制检查无 `global`、无 `sys.path` 篡改、统一用 `write_gb` 和 `run_tool`）+ `benchmark/frozen/`（论文预注册冻结评测机制）构成了极高水准的工程护城河。

---

## 五、分优先级优化建议路线图

### P0：公开发布 / 投稿前必须修复

1. **消除 `stage_write_final` 的 $O(N^2)$ 性能陷阱**
   * 在 `mitom_run.py` 的 `stage_write_final()` 入口，按唯一 `metadata_src_gb` 路径预先构建 `metadata_index`（`exact_dict, clean_dict, norm_dict`）和 `organism_index`，以参数传入 `finalize_genbank_output()`，将多记录 GenBank 批处理的元数据回写复杂度从 $O(N^2)$ 降至 $O(N)$。
2. **修复已确认的代码缺陷**
   * `trna_ml.py`：修复 `end_pairs()` 在 `s = 0` 时的环状坐标回绕；使 `_loci()` / `call()` 中的 `_ov()` 支持环状偏移 `(-L, 0, L)`，避免跨原点重复调用。
   * `mito_start.py` & `train_start_model.py`：保护 `candidates()` 永不丢弃原始起点 `a`、`StartModel.best()` 入口判空；将 `-999.0` 哨兵替换为 `np.nan` 并导出 `missing_go_to_left`；在 `start_model.json` 中记录 P1 规则版本号。
   * `mito_pcg_hmm.py`：将未使用的 `--k` 参数正确透传至 `prior()` 与近缘筛选，或移除冗余 CLI 选项。
   * `mito_orf.py`：为 `third_position_p()` 增加最小变异数门限（如 `sum(d) < 15` 时使用精确检验或保守回退），防止低计数下 $\chi^2$ 假阳性。
   * `mito_multipart.py`：将 `str.maketrans("ACGT", "TGCA")` 移出 `kmers()` 内层循环。
3. **理顺许可证与模型分发合规**
   * 在根目录补充 `LICENSE.GPL-3.0` 与 `THIRD_PARTY_LICENSES.md`，将 `README.md` 徽章/声明、`CITATION.cff` 统一更新为与 `pyproject.toml` 一致的复合许可说明（`CC0-1.0 AND GPL-3.0-only AND MIT`）。
   * 针对 `mito_intron.py` 的 `subtype_cms()`：默认改为仅使用本地自带 Rfam 模型（或由 `mitom setup --with-intron-subtypes` 显式下载、锁定 commit SHA 并校验 SHA256），禁止在注释运行期静默联网拉取未校验文件。
4. **修复打包清单与容器卫生**
   * `.dockerignore` 补齐 `Mitom_review`、`docs/reviews`、`docs/images`；
   * `packaging/bioconda/build.sh` 明确排除 `setup.py` 和 `batch_mitoz.py`。

### P1：显著提升稳定性、性能与可维护性

5. **统一 GBDT 推理范式，解除 `scikit-learn==1.7.2` 强锁**
   * 将 `mito_start.py` 的“纯 JSON 树 + `numpy` 推理”架构复用到 `trna_ml.py` 的 `caller` 与 `bounds` 两个分类器。收益一举三得：彻底移除运行时对 `scikit-learn==1.7.2` 的强版本锁；消除 `pickle` 反序列化风险；消灭首跑并发重复训练缓存的冷启动开销。
6. **消除多线程后的 `fork()` 风险**
   * 将 `mitom_steps.py`（`run_plots`）与 `mito_intron*.py` 中的 `mp.get_context("fork")` 替换为 `spawn`（配合显式传参取代 `_JOB` 全局共享）或独立子进程调度，保证 Python 3.12+ 与跨平台长期稳定。
7. **优化 tRNA 与后处理并发性能**
   * 在 `mito_trna.core_trnas()` 中实现 pan-CM `cmsearch` 输出复用，或对“非脊椎动物 + 全部 22 个家族检出且 $p \ge 0.9$”的基因组直接跳过冗余的完整 MiTFi 扫描；
   * 在 `mito_pcg_hmm.py` 中对 `(gene, ref_ids)` 建库结果做会话级缓存，并为 `run_submission`（`table2asn`）和 `mito_curate.py` 引入 `--jobs` 并行。
8. **解耦基于 `/note` 正则的内部隐式协议**
   * 在 `mito_core.py` 写出 `<prefix>.core.tsv`（目前已生成）或 sidecar JSON 时，保存结构化的核心比对区间与 tRNA 置信度；让 `mito_curate.py`、`mito_intron.py` 优先读 sidecar（缺失时再回退读 `/note`），解除步骤间脆弱的文本耦合。
   * 在 `Curator.run()` 中拆分为“Pass 1：全 CDS 移码/截断修复 (`origin_cut/frameshift/frame_repair/stop_trim`) $\rightarrow$ Pass 2：全 CDS 邻域边界精修 (`pcg/across_origin/pcg_check`)”，消除相邻基因遍历顺序依赖。

### P2：中长期架构演进

9. **重组为标准 `src/mitom/` Python 包并清理仓库**
   * 将 40+ 顶层 `.py` 归拢到 `src/mitom/` 子包（`core/`、`curate/`、`intron/`、`trna/`、`nuc/`、`rna/`、`cli/`），保留轻量的 `mitom` 入口；
   * 将 `benchmark/results/` 大体积评测产物迁至独立仓库或 Zenodo，清理归档 `Mitom_review/` 与 `docs/reviews/`，消除 `benchmark/trna_ml/scripts/` 中的本机硬编码路径，并尽量复用生产端函数防止双轨漂移。
10. **收敛门级启发式规则与产品边界**
    * 定期评估 `mito_taxa.py` 的门级特例规则与 `mito_start.py`（P11）的重叠度，在输出报告中清晰区分“生物序列证据”与“NCBI 历史标注惯例”；
    * 明确 `mito_nuc`（核 rDNA/组蛋白）与 `mito_rnaseq_annot.py` 的产品定位，若持续扩展建议拆为独立插件包，保持 MitoM 核心的简洁聚焦。