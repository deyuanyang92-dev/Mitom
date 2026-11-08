# MitoM 代码工艺评审（优雅性 / 模块化 / 可维护性）

**日期：** 2026-11-08  
**角色：** 高级软件架构视角（不评科学精度、不列功能 bug）  
**范围：** `feat/mitos2-args` 附近 HEAD；根目录约 36 个 `.py`，合计约 2.7 万行  

---

## 总评

**等级：可接受，但偏乱（acceptable → messy）**  

有拆分意图和共享工具纪律，**不是无法导航的泥球**；也**谈不上优雅**。真正的模块边界仍被兼容门面、巨型编排函数、模块全局状态、subprocess-CLI、非包布局打穿。  
**好维护吗？短期能改、长期贵。** 改一处旁路或 CLI 选项时，认知负担明显高于同类成熟工具（MITOS2 本体另论；相对 MitoFinder 这类更小的流水线，本仓编排层偏重）。

---

## 优雅吗？

**否。** 局部有干净抽象（`StepResult`/`optional_step`、`write_gb`、`mito_env` 工具表），整体观感是：

- 扁平脚本农场，不是可安装的包  
- 注释/docstring 极密（利于科学追溯），但无法补偿「找改点」成本  
- 命名 `mito_*` / `mitom_*` 有意图，实现上互相再导出、互相 `sys.path.insert`  
- 同一概念双轨：`run_tool` vs 裸 `subprocess`；`write_gb` vs 残留 `SeqIO.write`

优雅的标志是：**改一处意图清晰、测试不绑死私有符号、导入图无环**。当前都不稳。

---

## 模块化吗？

**半模块化：文件级有边界，系统级漏水。**

### 做得好的

| 层 | 文件 | 评价 |
|---|---|---|
| 共享 I/O / CLI | `mito_common` | 有正统入口（`write_gb`/`run_tool`/`cli_guard`） |
| 工具发现 | `mito_env` | 横切干净，宜保留 |
| 流水线拆片 | `mitom_inputs` / `stage` / `orient` / `steps` / `output` / `log` | 命名与职责大致对齐 |
| 步骤契约 | `mitom_steps.StepResult` | 小而清晰 |

### 假边界 / 漏水

1. **`mitom_pipeline.py`（~1678 LOC）**  
   上帝文件：再导出大量符号（含私有名），`workflow_run` / `build_parser` 各数百行。测试大量 `import mitom_pipeline as bm`，**把坏结构钉死了**。

2. **步骤用「再起一个 Python CLI」耦合**  
   `mitom_steps` spawn `mito_curate.py` 等。有隔离，但 API = 命令行字符串，错误通道扁平，难做进程内 mock/DI。

3. **`mito_common` 杂货袋（~630 LOC）**  
   I/O、几何、翻译、CLI、`CoreFeature`、还 lazy import `mitom_inputs` → 循环依赖风险。

4. **巨型叶子**  
   `mito_nuc` ~2305、`mito_intron_ref` ~2042、`mito_intron` ~1653、`mito_curate` ~1419：既是库又是完整 CLI。

5. **无正式包**  
   无 `__init__.py` / `pyproject.toml`；bioconda 复制 `*.py`；到处 `sys.path.insert(0, HERE)`。

6. **模块级全局当 IPC**  
   `LAST_PLAN`、`LAST_DEDUP_STATS`、`ORIENTATION_FALLBACK`、`CIRC_L`、`EXCLUDE_REFS` 等——并发与单测都脆。

---

## 好维护吗？

| 维度 | 判断 |
|---|---|
| 新人定位「注释主路径」 | 中等（有步骤日志与目录约定） |
| 改遗传密码 / 拓扑 / 验证 | 贵（多入口、多回退、测试绑 facade） |
| 加一个后处理步骤 | 中等（有 `optional_step` 模式） |
| 包成 pip/conda 正经库 | 差（需先包化） |
| 安全重构算法外结构 | 中等偏难（测试依赖私有再导出） |

**结论：** 维护成本已被审查轮次推高的「正确性补丁」摊薄了结构债；再堆功能而不收边界，会进入「每修一个审查项多一层特例」的模式。

---

## 严格建议（按优先级）

### P0 — 立刻停止恶化

1. **禁止扩大 `mitom_pipeline` 再导出面**；新代码直接 import 真实模块。  
2. 测试逐步脱离 `import mitom_pipeline as bm` 摸私有符号；facade 白名单压到少量公共 API。

### P1 — 三个月内收益最大

3. **`workflow_run` → 薄编排**（`stage_input` → annotate → reorient → round2 → post → finalize）；引入 `RunContext` dataclass，**新状态不进模块全局**。  
4. **`build_parser` → `mitom_cli.py`**。  
5. **正式包化**：`pyproject.toml` + console_scripts；过渡期根脚本只 `main()`；去掉散落 `sys.path.insert`。  
6. **统一契约**：新代码禁止裸 `SeqIO.write(..., genbank)`；`mito_nuc`/`mito_plot` 迁 `run_tool`；消灭 `mito_intron` 同名吞失败的 `run_tool`。

### P2 — API 化后处理

7. `mito_curate` / intron 抽出 `run(...) -> Result`；`mitom_steps` 直接调用；CLI 变薄；subprocess 作可选隔离后端而非唯一路径。

### P3 — 拆叶子（只切文件，不重写算法）

8. `mito_intron_ref`、`mito_nuc`、`mito_curate` 按 I/O | 搜索 | 合并 / 提交 切开，保留公开函数名。

---

## 不要重写的部分

- `mito_env` 发现顺序与 `TOOLS`/`REQUIRES`  
- `write_gb` 原子写 + `/pseudo` 约定、`cli_guard`、`run_tool` 日志格式  
- `StepResult` / `optional_step`、编号输出目录与 `OUTPUT_LAYOUT`  
- `mito_taxa` 表驱动与拼写兼容  
- **算法内核**（core / intron / tRNA / HMM / curate 规则）——结构上只切边界，不换科学实现  
- 多子命令产品表面与 `mitom debug`

---

## 一句话给负责人

把 MitoM 当「能发表的方法 + 能交付的软件」时：科学与补丁密度已经很高；**架构上差的是包边界和编排薄度，不是再多写一个审查修复。** 下一阶段最大杠杆是：**停 facade 扩张 → 包化 → `workflow_run` 变薄**，而不是重写生物学。
