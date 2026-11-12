# pstack 系统建议摘要（61b06af）

日期：2026-11-12。只读建议 + 已修 USAGE 步数漂移。未提交。

## 结论

编排够用（run 12 / annotate 8）。不要再拆 `mito_intron.py` / `mito_nuc.py`。对外短板是说明书与冻结评测。

## P0

1. 重做并提交文档分栏：`docs/README` + 精简 `USAGE` + `EXAMPLES` + `DATASETS`（此前本地整理未进 git，已丢）。
2. 新 benchmark 先写 `knowledge/TEST_REGISTRY.md`；数字只进一张表。

## P1

3. 推进 `benchmark/frozen/PLAN.md`。
4. USAGE 写清 `mitom rna` vs `mito_rnaseq_annot.py`（后者不是 run 的一步）。

## P2 / 不做

5. 可选：`RunState` 上显式 Mode；重渲流程图 SVG。
6. 不做：再拆大流程文件；取消子进程隔离。

## 本轮已改（未提交）

- `docs/USAGE_zh.md` §4：补第 8 步 ORF/NCR，改为 12 步；`run` / `compare` 文案对齐。
- `mitom_run.py`：`start_run` 注释改为 run 12 / annotate 8。

证据：`python3 -c "import mitom_run as m; print(len(m.RUN_STAGES), len(m.ANNOTATE_STAGES))"` → `12 8`。

## 后续 / Follow-up (2026-11-12, Claude)

- P0 第 1 条：文档分栏没有丢失，完整保存在 `stash@{0}`（2026-11-10），已恢复并提交（f228eff），随后中英双语（c7c7902）；
  stash 已删。The docs hub was not lost: it was in `stash@{0}`, restored (f228eff), made bilingual (c7c7902).
- 本轮所改的 `USAGE_zh.md` §4 与 `mitom_run.py` 注释已随 f228eff 提交。The edits listed above went in with f228eff.
- P1 第 4 条已做：USAGE § 5.7 有 `mitom rna` 与 `mito_rnaseq_annot.py` 的对照表。Done: USAGE § 5.7 table.
- 仍未做：冻结评测（`benchmark/frozen/PLAN.md`）、P2。Still open: the frozen evaluation, P2.
