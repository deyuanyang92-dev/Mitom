# RELEASE_READINESS_2026-11-08 处理结果 / Fix status (Claude Code, 2026-11-08, cody)

依据 / source: `RELEASE_READINESS_2026-11-08.md`（审 `feat/mitos2-args` @ 766653c）。回归测试：`tests/test_review_2026_11_08.py`。

| # | 审查项 | 核实 | 状态 | 改动 |
|---|---|---|---|---|
| 1 | `--mitos_args="-c2"` 绕过拦截 | 属实。另外 runmitos 用 argparse，缩写（`--cod 2`、`--outd x`、`--refd`、`--lin`）同样能绕过 | 已修 | `mito_engines._own_opt`：精确、`=`、粘连短选项（`-c2 -o/tmp -rX -RX -iX`）、长选项缩写都拒绝；实测 `mitom run --mitos_args="--best -c2"` 运行前报错，不建输出目录 |
| 2 | 只有 `runmitos.py`（2.1.9）时 doctor / debug 找不到 | 属实（doctor 用 `find("runmitos")`） | 已修 | 回退放进 `mito_env.find()`：所有调用者（doctor、debug、mito_trna、mitos_bin）一致；本机 2.1.9 环境实测找到 `runmitos.py` |
| 3 | 旁路子命令不认 `.gb.gz` 目录 | 部分属实：只含 `NC_000931.gb.gz` 的目录逐个实测 check、compare、curate、gb2fasta、validate、plot、dedup、reorient、intron 都正常（`prepare_input` 解压）；**只有 gb2tbl** 报 “no GenBank / EMBL file”（它自己列文件、不走 prepare_input） | 已修 | `mito_tbl`：`.gz` 文件按去掉 `.gz` 的后缀识别并解压（最多 3 层），输出名不带 `.gb`；非 GenBank 的 `.gz` 忽略 |
| — | 失败被当成步骤成功、`--mitos_args` 危险组合无警告、S4–S12、冻结评测 / 独立测试集 | 未处理 | 待定 | 属设计 / 评测工作，需用户决定优先级 |

运行中的错误与原因见 `docs/MITOS2.md` 第 6 节。
