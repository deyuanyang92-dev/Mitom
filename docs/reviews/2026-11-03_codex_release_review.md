# MitoM 发布前代码审查（Codex，2026-11-03）

审查基线：`655a75c`，版本 `1.0.0`，当前 `/mnt/p/Mitom` 工作区。仅记录问题，未修改业务代码、安装脚本、模型或现有测试。工作区已有的大量 tracked 差异在忽略行末空白后为空（`git diff --ignore-space-at-eol --stat`），本次新增内容仅为此审查记录。

结论：当前不具备发布放行依据。以下包含 8 项已复现的源码问题，以及 1 项当前工作副本的换行符阻断。P1 表示应在正式发布前处理，P2 表示特定输入或重复运行下的正确性问题。

审查重点：输入转换与拓扑、注释与最终输出流程、异常与退出状态、NCBI 验证、分组运行与输出锁、安装打包。另检查了核心 CDS 定位、坐标变换、tRNA 模型加载、内含子处理和边界修整的相关代码及现有测试；这不等于已独立证明所有类群的生物学注释准确性。

## R1 · P1 · 当前副本的 CRLF 使 Linux 命令和安装入口不能运行

位置：`mitom:1`、`install.sh:20`、`packaging/bioconda/build.sh:3`；Dockerfile 使用当前目录作为构建输入。

当前这些文件均为 `i/lf w/crlf`，即 Git 索引内是 LF，工作文件为 CRLF。实际执行：

```text
./mitom --version
/usr/bin/env: ‘python3\r’: No such file or directory

bash install.sh --help
install.sh: line 20: set: pipefail\r: invalid option name
```

影响：直接运行命令、安装以及从当前目录构建的 Docker / Bioconda 包均受影响。`bash -n` 可以通过，但不能检出 `set -euo pipefail\r` 的运行时问题。

边界：这是当前工作副本的问题，不能据此声称 Git 提交或从 GitHub 生成的源码包也含 CRLF。正式分发物需单独验证。

## R2 · P1 · table2asn 非零退出、没有验证产物时，仍报告零错误

位置：`ncbi_validate.py:119`、`ncbi_validate.py:126`、`mito_tbl.py:558`；结果进入 `mitom_steps.py:149` 的提交汇总。`mito_nuc.py:1864` 也有未检查退出码的相同模式。

`run_tool()` 默认 `check=False`，这里既不检查返回码，也不确认本次验证产生了有效产物；缺少 `.val` 时循环读取空列表，返回“没有消息”。

实际复现：把验证程序设为 `/bin/false`，`ncbi_validate.validate_record()` 返回 `[]`，没有生成 `.sqn`。随后通过 `run_submission()` 运行同一失败程序，得到：

```text
StepResult(ok=True)
1 records checked: 1 without NCBI ERROR, 0 with ERROR
record  ERROR  WARNING
sample  0      0
```

影响：程序崩溃、参数错误等基础设施失败被解释成验证成功，`summary.tsv` 的 `ncbi_errors=0` 和“可提交”说明不可信。`mito_tbl.validate()` 同样忽略返回码，且没有清理旧 `.val` / `.sqn`，重用目录时还可能读取上次产物。核标记路径的相同模式由代码检查确认，未运行真实核标记外部工具链。

## R3 · P1 · FATAL / REJECT 验证消息没有进入阻断错误计数

位置：`mitom_steps.py:157`、`mitom_steps.py:161`；`mito_tbl.py:663`、`mito_tbl.py:664`。

解析器支持 `ERROR|REJECT|FATAL`，但提交汇总只计数与字符串 `ERROR` 完全相等的消息。`gb2tbl` 也只读取 `Counter()['ERROR']`。

故障注入复现：让模拟验证程序以 0 退出，并写入一条符合现有解析器格式的 FATAL 消息：

```text
validation_messages.tsv:
sample  FATAL  SEQ_INST.BadDeltaSeq  ...

submission_summary.tsv:
sample  0  0
```

同时 `StepResult.ok=True`，摘要显示该记录没有 NCBI ERROR。此复现验证的是严重级别处理逻辑，消息本身为人工构造，未冒充真实 table2asn 输出。

影响：致命或拒绝级别的结果也被输出为 `ncbi_errors=0`。即使修复 R2 的程序退出码检查，这个计数问题仍需独立处理。

## R4 · P1 · gb2tbl 超时被记到空 SeqID，所有真实记录仍显示验证通过

位置：`mito_tbl.py:560`、`mito_tbl.py:561`、`mito_tbl.py:663`。

超时分支执行 `per['']['ERROR'] += 1`，而报告仅查询 `val[e.sid]`。空 ID 上的错误既未传播到该批次内的记录，也未作为独立失败行输出。

实际复现：模拟 `run_tool()` 抛出 `subprocess.TimeoutExpired`，执行 `gb2tbl --validate`。虽先打印了超时提示，随后却出现：

```text
[gb2tbl] table2asn: 0 ERROR(s) in 0 record(s); .sqn files to upload
gb2tbl_report.tsv: sample 的 NCBI_ERROR=0、NCBI_WARNING=0
退出码：0
```

影响：机器读取报告无法区分“验证通过”和“从未完成验证”，并被告知存在可上传的 `.sqn`。没有 SeqID 的其他验证错误也会落入同一个空键而丢失。

## R5 · P1 · --topology linear 生效前已经裁剪末端重复序列

位置：`mitom_inputs.py:759`、`mitom_inputs.py:763`、`mitom_inputs.py:774`。

FASTA 输入先按描述推断拓扑并裁剪末端重叠，之后才应用用户指定的 `--topology linear`。没有拓扑标签的长序列在前一步默认按 circular 处理。

实际复现：固定随机种子生成 12,000 nt，再追加其前 100 nt；FASTA 不标拓扑，调用：

```python
split_or_sanitize_fastas([fa], out, True, 15, topology='linear')
```

结果：输入 12,100 nt，输出 12,000 nt，最终拓扑却为 `linear`。日志也先称按 circular 删除了 100 nt。

影响：用户明确要求线性注释，序列仍被不可恢复地缩短；对于具有真实末端重复的线性分子，丢失的可能是生物学序列，而非组装重叠。此问题不依赖任何外部注释程序。

## R6 · P1 · import 在取得输出锁之前覆盖正在使用的输入

位置：`mitom_pipeline.py:1506`、`mitom_pipeline.py:1523`；实际写入位于 `mitom_pipeline.py:346`。

`main()` 先调用 `prepare_import()`，将本次输入写进 `<out>/00.imported/imported.gb`，然后才调用 `lock_output()`。

实际复现：在输出目录放入活进程 PID 对应的 `.mitom.lock` 和代表当前运行输入的 `00.imported/imported.gb`，再启动第二个 `mitom import`。

结果：第二个命令最终返回 2，提示另一个进程正在写入，但 `imported.gb` 已经被第二个输入覆盖。

影响：锁没有保护导入准备阶段。第一个任务后续读取注释、序列或元数据时，可能读取第二个任务的数据，产生错误结果或失败；被拒绝的命令仍污染了活动任务。

## R7 · P1 · 没有生成最终记录也可能返回成功

位置：`mitom_pipeline.py:854`、`mitom_pipeline.py:443`；最终写入异常仅被加入 `summary.final_write_failed`。

`workflow_run()` 的退出判断只关注“全部核心注释失败”，没有判断全部最终写入失败、最终产物为空等情况。`workflow_annotate()` 也有类似判断，并且 `samples=[]` 时直接绕过失败条件。

复现一：使用有效的、已有 CDS 注释的 GenBank 输入，走现有 input 引擎；注入 `finalize_genbank_output()` 写入失败。实际生成 `reports/final_write_failed.txt` 并列出 sample，`final_annotation/` 没有 GenBank 文件，但返回 0，并打印 `MitoM finished`。这个用例模拟文件写入失败，未修改源函数。

复现二：真实命令 `mitom annotate` 读取一条 2,000 个 N 的 FASTA，指定 `--clade Mollusca --genetic_code 5`，关闭后处理和绘图。输入被清洗阶段全部跳过，零最终记录，仍退出 0，日志为 `[DONE] exit=0`。该用例未替换注释或退出逻辑。

影响：批处理、调度器和发布验收可能把没有产出的运行判为成功。此处记录的是“零产出仍成功”，并未把所有可选步骤的失败都要求成主流程失败。

## R8 · P2 · 输入规范化后的路径碰撞会静默丢失整份输入

位置：`mitom_inputs.py:142`、`mitom_inputs.py:144`、`mitom_inputs.py:146`。

规范化将 `same.fa.gz` 写成 `00.input/same.fa`，同时也把原目录中的 `same.fa` 写到同一目标，没有唯一化或碰撞检查。EMBL 改成 `.gb` 后也可能产生同类碰撞。

实际复现：目录中同时包含 `same.fa` 的 record_A 和 `same.fa.gz` 的 record_B，两者序列不同。`normalize_input()` 返回的目录只有一个 `same.fa`，内容为 record_B；record_A 丢失且没有碰撞提示。

影响：后续去重和注释根本看不到被覆盖的记录，无法在去重报告里解释缺失原因。属于规范化写入覆盖，不是有意的序列去重。

## R9 · P2 · 复用输出目录时，旧比较表、原始注释和图片会混入新结果

位置：`mitom_output.py:413`、`mitom_output.py:447`；`mitom_steps.py:355`、`mitom_steps.py:372`。

工作目录清理和最终整理仅删除部分固定产物，没有清理 `original_vs_mitom.tsv`、`original_annotation/`、`plots/`。FASTA 运行跳过比较步骤，绘图关闭或记录集合缩小时，也不会删除对应旧产物。

实际复现：构造旧的这三类文件，顺序调用当前流程实际使用的 `clear_previous_work()`、`new_reports()`、FASTA 比较跳过分支和 `organize_outputs()`，新最终 GenBank 已写入，但上述旧文件全部保留。

影响：例如先运行一批 GenBank，再对另一批 FASTA 复用 `-o`，会留下上一批比较结果、原始注释与图。输出 README 还可能继续列出这些存在的目录或文件，使它们看起来属于本次运行。

## 验证记录与限制

- 最小复现脚本位于 `/tmp/mitom_codex_review_repro.py`，日志 `/tmp/mitom-codex-repro.log`。脚本只写临时目录，通过直接调用或故障注入验证上述行为，没有更改被审查函数。
- 全量现有测试命令：`PYTHONDONTWRITEBYTECODE=1 python3 -m pytest -q -p no:cacheprovider --basetemp=/tmp/mitom-codex-review-tests`。结果为 **707 passed / 2 failed / 1 skipped / 4 errors**，耗时 533.53 秒，日志 `/tmp/mitom-codex-pytest.log`。4 项 packaging 初始化错误均来自 R1 的 `build.sh` CRLF；不是 4 个独立缺陷。
- 两项 failed 为 `test_mitom_orbicella_only_ncbi_trnas`、`test_mitom_replaces_mitoz_trnas_keeps_pcgs`。输出报告均为 `mitom failed; input kept`；独立复现获得实际错误 `[Errno 30] Read-only file system: '/root/.local/share/mitom'`，属于本会话默认模型缓存目录不可写。使用项目支持的 `MITOM_DATA=/tmp/mitom-codex-ml-data` 单独重跑两项，**2 passed in 28.61s**，日志 `/tmp/mitom-codex-ml-pytest.log`；本报告不把这两项原始失败算作新注释缺陷。此次重跑调用了实际 cmsearch / cmalign 和模型，并未替换预测结果。
- 对顶层 Python 模块及 `mitom` 执行 `ruff check --select F --no-cache`，发现 18 项 F401 未使用导入；本轮未发现该组规则下的未定义名字问题。F401 为维护规范问题，未计入上面的 8 项逻辑缺陷。
- 当前 PATH 没有直接提供 blastx、cmsearch、table2asn、snakemake，但 `mito_env` 能自动发现其他 Conda 环境内的 blastx、cmsearch、cmalign 和 table2asn，以及 `/opt/mitos_refs`。真实完整流程的跨类群注释精度、完整 Docker 构建和 NCBI 提交文件有效性没有在本环境完成验证；不能用故障注入或单元测试代替这些验证。
- 本报告不把“未观察到其他问题”表述为“没有其他 bug”，也不作正式发布认证。
