# 说明书整理交接（给后续 Claude 核对）

日期：2026-11-10。机器：cody。未提交。

## 目标

使用说明对标 GetOrganelle / PhyloSuite：分栏、不挤；每命令有输入 / 输出 / 示例；特色；benchmark 只链表；内部数据集可重建。

## 本地文件（请先读这些）

| 文件 | 角色 |
|---|---|
| `docs/README.md` | 文档总目录（新） |
| `docs/USAGE_zh.md` | 使用手册（重排，约 632 行） |
| `docs/EXAMPLES.md` | 跟 `tests/data/` 的例子（新） |
| `docs/DATASETS.md` | `models/` 制备（新） |
| `README.md` | 落地页导航 + 精度表（数字只在此与 METHODS） |
| `docs/METHODS.md` | 方法与评测；步数 / ORF 进程说明已改 |

相关代码注释：`mitom_run.py`（start_run 12/8）、`mitom_log.py`（样例 `/12`）、`tests/test_mitom_debug.py`（fixture `/12`；旧版日志 `1/11` 保留）。

未纳入本次：`benchmark/results/rna_mappers_2026-11-10/`（无关未跟踪目录）。

## 已与代码对齐（核对入口）

```bash
python3 -c "import mitom_run as m; print(len(m.RUN_STAGES), len(m.ANNOTATE_STAGES))"
# 期望：12 8
```

- `RUN_STAGES` 含独立 `stage_orf_ncr`（第 8 步）。
- `ANNOTATE_STAGES` 8 步。
- ORF：`mitom_steps.orf_cmd` → 子进程 `mito_orf.py`。
- 有效 clade 名：`Annelida-segmented-worms`（不是 `Annelida`）。
- `mitom compare -o` 是 **TSV 文件路径**。
- `mitom gb2tbl` 目录输入默认输出名 `submission`（不是目录名）。
- `benchmark/train_start_model.py` 需要 `--set` 与 `--out`。
- 核库：`mitom nuc build --gb … -o models/nuc`。

## 请重点复查

1. USAGE §4 与 README 步骤表是否与 `RUN_STAGES` 一一对应（勿再改成 11）。
2. EXAMPLES 例 5 / 6 的 clade、金标准 classic 布局说明是否够清楚。
3. DATASETS 各脚本命令能否按脚本 `-h` / 头注释跑通。
4. 数字只写一次：USAGE 里是否还有应挪进表的实测句（minimap2 用时、gb2tbl ERROR 条数等）。
5. 是否重跑 `python3 scripts/render_workflow_diagrams.py` 同步 SVG 字幕（脚本字幕已是 12；图文件可能仍旧）。
6. 新页英标是否够用（CLAUDE.md 要求中英；现为中文为主 + 短英标）。

## 自测已跑

`pytest -q tests/test_mitom_debug.py` → 6 passed（改 fixture 后）。
