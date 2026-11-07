# MitoM 系统复查：缺陷 / 稳健性 / 优化 / 期刊标准（只记不改）

**日期：** 2026-11-08  
**分支：** `feat/mitos2-args` @ `766653c`（相对 `release/v0.0.1` / `a70f8eb` 多 2 笔：`.gb.gz` import、`--mitos_args`）  
**工作区：** `docs/MITOS2.md` 有未提交修改（curation 对比表改写）  
**约束：** 只审查、不改代码。

---

## 0. 总判决

| 问题 | 判决 |
|---|---|
| 上一轮遗留的 import `.gb.gz` | **已修**（`6dfca2c`） |
| 新功能 `--mitos_args` | **有真 bug**：粘连短选项可绕过拦截并覆盖 `-c/-o/-r` |
| 工程稳健性（相对 11-05/06 审查） | 主路径明显好转；部分旁路仍脆 |
| 是否达到 Bioinformatics / NAR 软件文常见门槛 | **尚未**。更接近「可试用的测试版 / 预印本工具说明」 |

---

## 1. 仍存在的 Bug（已核对代码）

### B1 — 高：`--mitos_args="-c2"` 等粘连短选项绕过拦截

- **位置：** `mito_engines.mitos_extra_args`、`mitos2_command`
- **机制：** 只拒绝整 token ∈ `MITOS_OWN_OPTS`（`-c`、`--code`、`-o`…）。`-c2`、`-rrefseq63m`、`-o/tmp/x` 被放行。
- **命令顺序：** MitoM 先写 `-c 5 -o … -r …`，再把 `extra` 接在后面。argparse 惯例「后者覆盖前者」。
- **本机探针：** `mitos_extra_args("-c2")` → `['-c2']`；`mitos2_command(..., extra=['-c2'])` 末尾带 `-c2`。
- **后果：** 用户以为用了 `--genetic_code`，MITOS2 实际可能用另一套密码/输出目录/参考版本；结果难查。
- **建议：** 拒绝以 `-c/-o/-r/-R/-i` 开头的粘连 token，或对 short opts 做与 runmitos 同构解析后再判；加回归测试「后置覆盖」。

### B2 — 中：仅有 `runmitos.py`（bioconda mitos 2.1.9）时发现不一致

- **位置：** `mito_env.mitos_bin` 已认 `runmitos.py`；`TOOLS` / `find("runmitos")` / doctor 版本行 / `mito_trna` 仍偏 `runmitos`。
- **后果：** `--engine mitos2` 能跑，doctor / `--debug` / 部分 tRNA 路径仍可能报「runmitos not found」，日志缺可靠 MITOS2 版本。
- **建议：** 统一走 `mitos_bin()`；engine=mitos2 时强制记版本到 `[ENV]`。

### B3 — 低～中：旁路仍不认 `.gb.gz`

- `import` 主路径已修。
- `iter_gb_files`、`_is_gb_file`、`mito_check` 等未先 `normalize`/`prepare_input` 的路径仍看 `p.suffix`，目录里只有 `.gb.gz` 时可能当空输入。
- 不阻塞 `mitom run`（会先解压），但单独子命令仍脆。

---

## 2. 不够稳健（非必现崩溃，但生产上会踩）

| ID | 主题 | 说明 |
|---|---|---|
| R1 | 失败语义 | 有 NCBI ERROR 时步骤仍可 `ok=True`（注释有效、错误在报告里）——对「可提交」叙事仍易误读 |
| R2 | `--mitos_args` 语义护栏 | 合法参数如关蛋白/tRNA、改 `--circrot` 可空注释或打乱方向，现无 WARNING |
| R3 | 外部工具失败像「0 结果」 | 多数路径已修；自定义 `mitos_args` 失败时需确认日志够醒目 |
| R4 | 测试干扰 / 环境 | 历史曾有全量 pytest 顺序干扰；需保持隔离 |
| R5 | 依赖版本漂移 | doctor 已查 scikit-learn；用户 PATH 上旧 `mitom` 仍是常见踩坑（install 文档已写） |
| R6 | 算法层 S4–S12 | 内含子前移、rRNA 软边界、NUMT、多分体、置信度 TSV、评测去泄漏——未做，属精度债 |

---

## 3. 值得优化（高价值优先）

1. **堵住 B1**，并统一 MITOS2 发现/版本（B2）。  
2. **旁路 `.gz` 一视同仁**：`iter_gb_files` / check / curate 入口统一 `_plain_suffix` + `_read_text`。  
3. **每基因置信度 / 备选坐标 TSV**（原 S11）：人工复核成本立刻下降，也利于论文「可追溯」。  
4. **评测去泄漏 + 冻结版本重跑**（原 S12 / 发表建议）：否则 84%、零 ERROR 不宜作主结论。  
5. **rRNA 改报软边界**（S6）：避免被 exact-NCBI 天花板卡死叙事。  
6. **PCG 内含子感知前移**（S4）：刺胞大内含子等难例。  

低优先级：NUMT、多分体套件、`--mitos_args` 危险组合 WARNING。

---

## 4. 是否达到 Bioinformatics / NAR 软件文标准

对照该类期刊常见要求（可复现、独立评测、与对照工具公平比较、明确局限）：

| 要求 | 现状 |
|---|---|
| 可用软件 + 文档 + 安装 | **基本具备**（install、doctor、USAGE、MITOS2 指南、CHANGELOG） |
| 方法可复述 | **部分**：METHODS 有；工具版本/完整命令需再固化 |
| 独立测试集 / 训练泄漏审计 | **不足**（P11 全量、refseq89m、内含子库接触未统一清单） |
| 冻结版本上重跑主数字 | **不足**（v30 等旧报告 ≠ 当前 HEAD 自动成立） |
| 与 MITOS2/MitoZ 等公平对照 | **有材料**；需统一协议与版本号后写进主文 |
| 正确性可信（验证不假零 ERROR） | **主路径已明显改善**；宣传「零 ERROR」仍建议独立复核后再写 |
| 局限与真值分层 | **文档有免责**；论文需把 NCBI 习惯 vs 生物学边界写清 |

**结论：**  
- **工程测试版 / 组内使用 + 预印本：** 可以。  
- **Bioinformatics Application Note / NAR web server·software 投稿级：** 尚差「冻结版本 + 独立评测 + 泄漏审计 + 对照协议」四块；B1 应在对外版本前修掉。

建议对外叙事：先标 **0.0.x 测试版**，修好 B1/B2，冻结一版做独立重测，再谈正式发表数字。

---

## 5. 与既有审查文件关系

| 文件 | 关系 |
|---|---|
| `FIX_REVIEW_2026-11-05.md` F1–F5 | 后续提交已修；本轮未见回退 |
| `SECOND_OPINION_2026-11-06.md` | import 保护输入、线性 tRNA、timeout——已合入 |
| `ALGORITHMIC_GAPS_AND_SOLUTIONS_2026-11-04.md` S4–S12 | 仍是精度/发表方向，未实施 |
| `docs/reviews/2026-11-03_publication_recommendations.md` | 期刊门槛细则仍适用 |
| **本文件** | 2026-11-08 系统复查；含新功能 B1/B2 |

---

## 6. 建议给 Claude Code 的顺序（若继续修）

1. B1 `--mitos_args` 粘连短选项 + 测试  
2. B2 MITOS2 发现/版本统一  
3. B3 旁路 `.gb.gz`  
4. 冻结版本 → 独立评测（发表路径）  
