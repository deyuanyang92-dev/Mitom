# Mitom_review · 本机核实汇总报告

**核实日期：** 2026-11-04  
**仓库副本：** `/media/cody/18T_R1/yzz/Mitom`  
**核实 commit：** `1912f7bf5c9efe3d5887e4584036d1c5f2ae8a41`（`main`，含 `Mitom_review/` 与 `docs/reviews/`）  
**审查原文基线：** 多份报告针对 `655a75c`（2026-10-31）；本机 HEAD 在其后，逻辑问题多数仍在。

**审查材料：**
| 文件 | 来源 |
|---|---|
| `Mitom_review/MitoM_code_review_Grok4.7_2026-11-04.md` | Grok 4.7 |
| `Mitom_review/MitoM_code_review_GLM5.3_2026-11-04.md` | GLM 5.3 |
| `Mitom_review/MitoM_code_review_KimiK3_2026-11-04.md` | Kimi K3 |
| `Mitom_review/MitoM_代码审查报告_2026-11-03.md` + 第二轮 | 早期静态审查 |
| `Mitom_review/两套审查对比与我的盲区_2026-11-03.md` | 对比说明 |
| `docs/reviews/2026-11-03_codex_release_review.md` | Codex 发布前审查（R1–R9） |
| `docs/reviews/2026-11-03_publication_recommendations.md` | 发表建议（非缺陷清单） |
| `docs/reviews/2026-10-21/*` | 旧轮 rRNA/CLI/curate 修复记录 |

**方法：** 不凭经验；对照本机二进制 + 帮助输出 + 最小复现 / 真实 fixtures。复现脚本输出：`/media/cody/18T_R1/yzz/_verify_mitom_reviews/`。

---

## 0. 本机软件版本（核对基准）

| 软件 | 本机实测 | 出处 |
|---|---|---|
| MitoM | **1.0.0** | `./mitom -V`；commit `1912f7b` |
| Biopython | **1.84** | `python3 -c 'import Bio'`；`mitom doctor` |
| nhmmer | **HMMER 3.4** | `mitotools` / `BUSCO_v5.8.2` env（`mito_env.find` 当前指向 BUSCO env） |
| table2asn | **1.28.1094** | `/home/cody/miniforge3/envs/mitotools/bin/table2asn` |
| Python（运行 mitom） | 能 `import Bio` 的解释器（launcher 可 re-exec） | `mitom` |

克隆：GitHub 直连 TLS 失败；**HTTP 代理 `127.0.0.1:7890` 可用**，7892 不可用。

---

## 1. 总判决

- **Codex R1（CRLF）在本机 git clone 上不成立**（工作区 LF，`./mitom -V` 正常）。Codex 原文已写明边界：属当时工作副本，不能代表 GitHub 源码包。
- **Codex R2–R9 在本机均复现为真**（R2/R5/R7 动态；R3/R4/R6/R8/R9 代码+局部复现）。
- **三模型共识的工程项**（`fork`、`/dev/null`、`/pseudo\tNone`）本机为真。
- **docs 发表建议**合理：特别点出 R2–R4 会削弱「零 ERROR」宣传可信度——本机已证实验证失败可被计成 0 ERROR。
- **真实数据：** `mitom curate -i tests/data` → table2asn **ERROR 9→0**；剥掉 rRNA 的 `NC_000931` 上 `rrna-hmm search` 成功加入 2 条（coverage 0.99/0.97）。

**发布结论（与 Codex 一致）：当前不具备「验证可信 + 失败可感知」的放行依据；应先修 R2–R7 与 `/pseudo` 提交表，再谈零 ERROR 宣传。**

---

## 2. Codex R1–R9 核实表

| ID | 声称 | 本机判决 | 证据 |
|---|---|---|---|
| R1 | CRLF 使 `./mitom` / `install.sh` 不能跑 | **假（对本 clone）** | `file`/`\\r\\n` 检查均为 LF；`./mitom -V` → `MitoM 1.0.0`；`bash install.sh -h` OK。边界：Windows `autocrlf` 工作区仍可能中招。 |
| R2 | table2asn 失败仍报 0 ERROR | **真** | `validate_record(rec, …, "/bin/false")` → `[]`，无 `.sqn`；源码不读 `returncode`。 |
| R3 | FATAL/REJECT 不进阻断计数 | **真** | `mitom_steps.py:157` 只统计 `("ERROR","WARNING")`；`n["ERROR"]` 不含 FATAL。 |
| R4 | gb2tbl 超时记到空 SeqID | **真** | `mito_tbl.py:561` `per[""]["ERROR"] += 1`；报告读 `val.get(e.sid)`。 |
| R5 | `--topology linear` 前已裁末端重复 | **真** | 12100 nt（末 100=首 100）+ `topology="linear"` → 输出 **12000 nt** 且 `topology=linear`；日志先按 circular 删除。 |
| R6 | `prepare_import` 在 `lock_output` 前 | **真** | `mitom_pipeline.py` 源码顺序：`err = prepare_import` 先于 `err = lock_output`。 |
| R7 | 零最终产物仍 exit 0 | **真** | `mitom annotate` 对 2000×N FASTA：`[DONE] exit=0`，无 final GenBank。 |
| R8 | 规范化路径碰撞静默丢输入 | **真** | 同目录 `same.fa`(A) + `same.fa.gz`(B) → 仅保留 **B** 内容，无碰撞警告。 |
| R9 | 复用 `-o` 旧比较表/图混入 | **未完整动态复跑** | 代码路径与 Codex 描述一致（清理不全）；标为 **很可能真**，未再造三类旧文件全流程。 |

早期静态审查「漏掉 Codex 6 个 P1」的自评（`两套审查对比…`）**成立**：本机也确认那些项主要靠动态复现才能钉死。

---

## 3. 三模型（2026-11-04）共识 / 高优先级抽样

| 主题 | 来源 | 判决 | 证据 |
|---|---|---|---|
| `multiprocessing.get_context("fork")` | GLM 1.1 / Kimi 1.2 | **真** | `mitom_steps.py` 绘图 Pool |
| nhmmer `-o /dev/null` | GLM 1.2 / Kimi 1.3 | **真** | `mito_rrna_hmm.py:242` |
| `.tbl` 写出 `pseudo\tNone` | Grok 1.5 / GLM 3.1 / Kimi 1.6 | **真** | `feature_lines` → `\t\t\tpseudo\tNone` |
| Chordata 硬编码 code 2 | Grok 1.2 / GLM 3.2 / Kimi 3.1 | **真（代码）** | `PHYLUM_CODE["Chordata"]=2`，`CLADE_CODES["Chordata"]=2`；未跑被囊动物全注释 |
| 吸虫应用 code 9 而非 21 | Grok 1.1 | **真（代码）** | `Platyhelminthes: 9`；无 Trematoda→21 |
| 部分失败 exit 0 | 多份 / Codex R7 | **真** | 见 R7；与「仅全部 annotation_failed 才返回 1」一致 |
| Windows stale-lock `os.kill` | GLM/Kimi | **未在 Linux 复现** | 代码存在；本机 Linux 不触发「杀错进程」语义，标 **平台相关、未测** |

未逐条复现：Grok 内含子/原点旋转等高密度算法项、Kimi `mito_intron` 目录规范化等——需单独用例；报告中不标「已否证」。

---

## 4. docs 意见核对

| docs 主张 | 核对 |
|---|---|
| Codex：不具备发布放行依据 | **同意**（R2–R7 本机为真） |
| 发表建议：先修正确性，再冻结评测；「2497 零 ERROR」须独立复核 | **同意**——R2–R4 直接削弱该数字可信度 |
| 发表建议：旧 v30 / 全量 P11 不能当独立泛化 | **方法正确**；本轮未重跑千级评测 |
| 2026-10-21 A/B/C_fixes（rRNA/CLI/curate） | 相关符号仍在（`p95_len`、`startswith("added")`、`MAX_ORIGIN_OVERLAP`、`hint` 列等）；本机 curate/rrna-hmm 实跑通过 |

---

## 5. 实际数据测试

| 测试 | 结果 |
|---|---|
| `mitom curate -i tests/data` | ERROR **9→0**，WARNING 8→18；`validation_compare.tsv` 已存 |
| `rrna-hmm search` 于剥 rRNA 的 `NC_000931` | 2 rRNA added；coverage 0.99 / 0.97 |
| pytest 抽样（ncbi_validate / topology / input_hardening / curate / rrna_hmm） | **173 passed**（12 s） |
| `mitom doctor` | 必需工具 OK；`runmitos`/`arwen` 可选缺失 |

日志目录：`/media/cody/18T_R1/yzz/_verify_mitom_reviews/`（含 `repro_codex_r.txt`、`curate_head/`、`pytest_smoke.log`）。

---

## 6. 建议修复优先级（核实后）

1. **P0：** R2（table2asn 退出码/产物）、R3（FATAL/REJECT 计数）、R4（超时 SeqID）、`/pseudo` 无值限定符  
2. **P0：** R5（`--topology linear` 禁止先按 circular 裁剪）、R7（零产物勿 exit 0）  
3. **P1：** R6 锁序、R8 路径碰撞、R9 输出清理  
4. **P1：** `/dev/null` 与 `fork`（跨平台）；Chordata/Trematoda 遗传密码表  
5. **发表：** 按 `publication_recommendations.md` 冻结版本后再报「零 ERROR」与千级一致率  

---

## 7. Claude Code 终核

本机 `claude` 2.1.312 经代理 7890 调用后返回：

`You've hit your weekly limit · resets Nov 5, 12pm (Asia/Shanghai)`

（日志：`/media/cody/18T_R1/yzz/_verify_mitom_reviews/claude_final.log`）。额度恢复后用本报告做终核即可；本轮核实结论不依赖该步。
