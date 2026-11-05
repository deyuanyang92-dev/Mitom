# Mitom_review 终核

**日期：** 2026-11-04  
**对照报告：** `Mitom_review/VERIFICATION_REPORT.md`（commit `1912f7b`，MitoM 1.0.0）  
**说明：** 本机 `claude` 2.1.312 经代理 7890 调用后先有 403，后为周限额：`You've hit your weekly limit · resets Nov 5, 12pm (Asia/Shanghai)`（`_verify_mitom_reviews/claude_final.log`）。未能用 Claude Code 官方会话终核；本文件为独立抽查终核。额度恢复后再跑一次更稳妥。

---

## 判决

**同意** Cursor / `Mitom_review/VERIFICATION_REPORT.md` 总判决：

- Codex **R1 对本机 clone 为假**；**R2–R8 为真**（R9 代码一致、未全流程复跑）。
- 三模型共识项（`fork`、`/dev/null`、`pseudo\tNone`、Chordata=2、Platyhelminthes=9）本机为真。
- **当前不具备「验证可信 + 失败可感知」的发布放行依据**；应先修 R2–R7 与 `/pseudo`。

---

## 抽查证据（≥8 条）

| # | 报告项 | 本机 |
|---|---|---|
| 1 | R2 table2asn 失败仍空 msgs | `validate_record(..., "/bin/false")` → `[]` |
| 2 | R3 只计 ERROR/WARNING | `mitom_steps.py` 可见；FATAL 不进阻断计数 |
| 3 | R5 linear 仍先按 circular 裁 | 核实日志：12100→12000，`topology=linear` |
| 4 | R6 prepare_import 先于 lock | `mitom_pipeline.py` 源码顺序成立 |
| 5 | R7 零产物 exit 0 | 核实日志：`[DONE] exit=0`，无 final GenBank |
| 6 | `pseudo\tNone` 进 .tbl | `feature_lines` 探针复现 |
| 7 | `get_context("fork")` | `mitom_steps.py` |
| 8 | nhmmer `-o /dev/null` | `mito_rrna_hmm.py` |
| 9 | Chordata=2 / Platy=9 | `mito_taxa.PHYLUM_CODE` |
| 附 | 真数据 curate / rrna-hmm | `_verify_mitom_reviews`：ERROR 降为 0；剥 rRNA 后 2 条 added |

---

## 报告需注意处（非否定总判决）

1. 早期误写的「10-21 修复已全部落地 → 可放行」**不适用于**本轮 Mitom_review；10-21 项仍在，但 **Codex R2–R7 是发布阻断项**。
2. R9 标「很可能真」可保留；不必升为已动态钉死。
3. Claude Code 403 未清：用户需 `claude` 登录/配额恢复后再终核一次更稳妥。

---

## 开放项

| 项 | 建议 |
|---|---|
| 修 R2–R7 + pseudo | **仍开放（P0）** |
| R8/R9、fork、/dev/null、遗传密码表 | **仍开放（P1）** |
| 千级「零 ERROR」宣传 | **冻结版本并独立复核前勿用** |
| R1 CRLF | **对本 clone 可关闭**；Windows autocrlf 边界保留一句即可 |
| Claude Code 鉴权 | **仍开放**（403） |

---

## 建议下一步

1. 按报告 §6 修 P0（R2–R7、pseudo）。  
2. 鉴权恢复后：`claude -p` 只读 `Mitom_review/VERIFICATION_REPORT.md` + 本文件做第二终核。  
3. 修完后用 `_verify_mitom_reviews/repro_codex_r.txt` 同类探针回归。
