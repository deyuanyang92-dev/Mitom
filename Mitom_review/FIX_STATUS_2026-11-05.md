# MitoM 审查项处理结果（2026-11-05，Claude Code）

**分支：** `fix/review-2026-11-05`（基于 `1912f7b`）
**依据：** `VERIFICATION_REPORT.md`、`NEW_BUGS_2026-11-04.md`、`ALGORITHMIC_GAPS_AND_SOLUTIONS_2026-11-04.md`
**回归测试：** `tests/test_review_2026_11_05.py`（每个已修项至少 1 个用例）；全量 `pytest tests`

状态：**已修** = 代码已改并有测试；**不成立** = 核实后报告结论不对；**不改** = 属实但属于设计选择或收益低，原因见表。

## 一、改之前在本机核实时的新发现

| 发现 | 说明 |
|---|---|
| table2asn 1.28 的严重级别 `Critical` 原来完全没被读取 | `MSG_RE` 只认 INFO/WARNING/ERROR/REJECT/FATAL。实测越界 feature 报 `Critical: valid [SEQ_FEAT.Range]`，被静默丢掉。现在按 REJECT 计入 ERROR 列 |
| table2asn 的**退出码不能用来判断是否崩溃** | 实测：有校验错误时，单记录 exit 0，多记录 exit 2，`.val`/`.sqn` 都照常写出。所以 R2 不能按报告建议的"读 returncode"来修，否则会把正常结果全记成失败。现在只按"没有 `.val`"判定为未校验 |
| `pseudo\tNone` 的真实触发点 | `mito_curate.py` 用内存里的记录做"after"验证，P10 把 `/pseudo` 设成了 `[None]`；`mito_intron_ref` 的 `trans_splicing` 也是这样 |
| 全 N FASTA 用 `annotate` 仍 exit 0 | 0 个样本时没走到最终产物检查；现与 `run` 一致，exit 2 |

## 二、逐项状态

### P0（验证可信度、退出码）

| ID | 状态 | 改动 |
|---|---|---|
| R2 / D11 / N6 / N15 | 已修 | `ncbi_validate.validate_record`、`mito_tbl.validate`、`mito_nuc.write_submission`：没有 `.val`（崩溃、程序缺失）记为 `MITOM.Table2asnFailed` ERROR；提交步骤 `ok=False` |
| R3 | 已修 | `mitom_steps.run_submission`：REJECT/Critical/FATAL 计入 ERROR 列（`BLOCKING`） |
| R4 | 已修 | gb2tbl 超时或失败时，该对文件的每个 SeqID 各记 1 个 ERROR（原来记到空 SeqID 下，没有任何行会读到） |
| R5 | 已修 | `--topology linear` 时不再按环状裁掉末端重复（端到端实测：12100 nt 进，12100 nt 出） |
| R7 | 已修 | `run`/`annotate` 没有写出任何最终 GenBank → exit 1；`annotate` 没有可注释序列 → exit 2 |
| D12 | 部分 | 记录有 NCBI ERROR 时步骤仍算成功（注释结果本身有效，错误已写进报告）；只有未能校验时才 `ok=False` |
| N1 / N2 / D8 | 已修 | `feature_lines`：值为 None 或空时只写 key |

### 遗传密码

| ID | 状态 | 说明 |
|---|---|---|
| D1 | 已修 | `auto_clade`：CDS 没有 `transl_table` 时，用 `code_of_lineage` 取参考库里同类群的密码（被囊动物 13） |
| D2 | 已修 | `mito_intron.run_exonerate` 用 `mito_core.exonerate_code()`（code 24/33 传 64 字母表） |
| N12（吸虫→21） | 不改 | 参考库 34 条吸虫记录：32 条 `transl_table=9`，2 条为 14。现在按数据决定（9），与 GenBank 提交习惯一致；9 和 21 只差 ATA 一个密码子 |
| N13（被囊→13） | 已修 | 即 D1；`run` 路径原本就给 13 |
| D15（半索动物→24/33） | 不成立 | 肠鳃纲用 code 9（参考库 3/3）；24/33 只用于羽鳃类 |
| N5（blastx `-query_gencode 5`） | 不改 | 只用于找门；密码另由比对上的密码子判定（`_codons`） |
| D5（分类失败回落到 Annelida/5） | 不改 | 这是 `--clade auto` 出现之前的默认行为，会给出 WARNING 并写入 `clade_decision.tsv` |

### 安装

| ID | 状态 | 说明 |
|---|---|---|
| D3 | 已修 | Miniforge 下载名用 `MacOSX` |
| D13 / E23 | 已修 | 不管 `--envs` 怎么给，都一定建 `mitotools`（wrapper 依赖它） |
| D14 | 已修 | 不带 `-y` 的重试改从 `</dev/null` 读，不会卡在等键盘输入 |

### 输出目录、输入

| ID | 状态 | 说明 |
|---|---|---|
| R6 | 已修 | `import` 先加锁，再写 `00.imported/` |
| R8 | 已修 | 解压或转换后重名（`same.fa` + `same.fa.gz`）→ InputError |
| R9 / N14 | 已修 | 复用 `-o` 时删除上一轮的结果（OUTPUT_LAYOUT、annotation-summary.tsv、all.final.gbf、README.txt）；`-i` 在 `-o` 之内时绝不删除输入 |
| N8 | 已修（部分） | 锁的持有进程属于其他用户时（`PermissionError`）视为仍在运行，原来会被当成陈旧锁抢走。PID 复用不处理：报错信息已给出锁文件路径 |
| D19 | 已修 | atexit 释放锁时容错 |
| N9 / D6 | 已修（GenBank） | LOCUS 行没有 topology 且 < 10 kb → linear（与 FASTA 规则一致）。≥ 10 kb 且无标签的 FASTA 先按环状注释、再判定，这是基于基准测试的设计（CHANGELOG），不改 |
| E19 | 已修 | 双层 gzip 也能解开 |
| E1 / E2 / E9 | 不成立（CLI 层面） | 所有命令都先经过 `prepare_input` / `normalize_input` 解压 `.gz`；只有直接调用内部函数时才会出现 |
| D4（CRLF 分组） | 不成立 | `read_text()` 会把 `\r\n` 转成 `\n`，实测能正确切成 2 块 |
| N11 | 已修（主要路径） | core、curate、trna、import、dedup、nuc 写 GenBank 改用 `write_gb`（原子写入） |
| N4 | 已修 | Step 4 重新计数失败时记 WARNING，不再静默 |
| N10 | 已修 | transfer_metadata / replace_organism 有文件失败 → exit 1 |
| N7 | 已修 | intron CM 下载的文件名只取 basename，且必须以 `.cm` 结尾 |

### 分组、注释步骤

| ID | 状态 | 说明 |
|---|---|---|
| N3 | 已修 | 失败组的注释结果不并入（它的 reports 和日志保留） |
| D17 | 已修 | `_concat_tsv` 表头不同时按列名合并 |
| E15 | 已修 | `mito_rescue.present_genes`：能识别 cox1/COI/rrnS 等命名 |
| E16 | 已修 | `reorient_check` 同样把基因名规范化 |
| E22 | 已修 | 无法平移到新原点的 feature 改为丢弃并记 WARNING（原来保留旧坐标，位置是错的） |
| E26 | 已修 | rescue 的 MITOS2 cmsearch 线程数 = CPU 数 ÷ 外层并行数（`mitos2_command(jobs=)`） |
| D7 | 已修 | `mito_trna` 预测器非零退出时打印 WARNING，不再像"0 个 tRNA" |
| E3 | 已修 | 提交用的 `.tbl` 里，NCR `/note` 去掉 "(MitoM)" |
| E11 | 不改 | 只在非标准写法 `complement(join(1..50,951..1000))` 时出现；NCBI 标准写法实测结果正确 |
| E4 | 不成立 | curate 选邻居时只看 CDS/tRNA/rRNA，不看 misc_feature（NCR）；跨原点的情况由 `spans()` 处理 |
| E20 | 不改 | ≥ L/2 的间隙按设计不当作内含子（防止跨原点误判） |
| D9 | 不改 | 没有 CDS 的记录不需要 `[mgcode=]` |
| D10 | 不改 | 越界 feature 保留用户原始数据并给出警告；table2asn 的 Critical 现在会计入 ERROR |
| D16 | 不改 | plot 的 ORIGIN_EXT 只用于把 reads 映射到原点附近，与 core 用途不同 |
| D18 | 不改（加报错） | BLAST 本身不支持带空格的路径（makeblastdb 实测失败），改为给出明确的报错 |
| E21 | 不改 | 下载时只追加写入，被截断的只可能是最后一条，这种情况已处理 |

### 导入、元数据

| ID | 状态 | 说明 |
|---|---|---|
| E5 | 已修 | MitoZ 输出目录之外的 GenBank 也会导入（MitoZ 目录内仍只取 `*.result/*.gbf`） |
| E24 | 不改 | 同一基因组同时给了 MitoZ 和 MITOS2 结果属于用户的输入选择；导入后 dedup 会按序列去重 |
| E6 | 不改 | 文档已写明 MitoZ 结果默认按 circular，线性 contig 需加 `--topology linear` |
| E7 | 不改 | 与 core 一致：`/organism` 先用序列名，再由 transfer_metadata / replace_organism 替换 |
| E12 | 已修 | 重名改名时跳过已存在的 `id_N` |
| E13 | 已修 | accession 先取完整 id，再取 LOCUS name |
| E14 | 已修 | smart 匹配也会更新 LOCUS 长度（`_with_length`） |
| E17 | 不成立 | 同一轮运行里样本名唯一 |
| E8 / E18 / E25 | 不改 | 名字全是符号这类输入不现实；样本名在 `split_or_sanitize_fastas` 里另有去重 |
| `fork` / Windows 锁 | 不改 | 不支持原生 Windows（`install.sh` 直接退出）；macOS 上 fork 可用 |
| `/dev/null` | 已修 | 改用 `os.devnull` |

### 算法建议 S1–S12

S1–S3 已由上面的项覆盖。S4–S12（内含子前移、起止证据分层、rRNA 软边界、NUMT、多分体、置信度 TSV、评测去泄漏）是研究和发表层面的方向，不是缺陷，本轮没有实施。其中 S12 应在冻结版本后、更新 README 的"零 ERROR"数字之前完成。

## 三、验证（2026-11-05，mitotools 环境 Python 3.11 / Biopython 1.88，table2asn 1.28.1094）

| 检查 | 结果 |
|---|---|
| 全量 `pytest tests`，修复前（1912f7b 的 worktree） | 709 passed，2 failed，3 skipped |
| 全量 `pytest tests`，修复后 | **736 passed**（新增 27 个），2 failed，3 skipped |
| 那 2 个失败：`test_mito_trna.py::test_mitom_orbicella_only_ncbi_trnas`、`::test_mitom_replaces_mitoz_trnas_keeps_pcgs` | 修复前就有；单独运行 `tests/test_mito_trna.py` 时两边都通过，只在全量运行时失败（测试之间的顺序干扰），与本轮改动无关，尚未处理 |
| `mitom curate -i tests/data`，修复前后对比 | 两边相同：ERROR 9→0，WARNING 8→18（`Critical` 与 table2asn 判定的改动不改变这批数据的计数） |
| `mitom run --clade auto`，NC_006902 + NC_010217（海绵 + 珊瑚，两个分组） | exit 0，两组合并，两条记录均 0 ERROR |
| 全 N FASTA：`mitom annotate` | exit 2（修复前为 0） |
| 末尾 100 nt 重复首端的 12100 nt 序列，`--topology linear` | 输出 12100 nt linear（修复前为 12000） |
| 同一个 `-o` 先跑 GenBank、再跑 FASTA | `original_vs_mitom.tsv`、`original_annotation/` 被删除，用户自己的 `my_notes.txt` 保留 |
| 实测 table2asn | 有校验错误时，单记录 exit 0，多记录 exit 2，`.val` 都会写出：退出码不能用来判断崩溃 |
| 实测 makeblastdb | 路径带空格时直接失败（D18 改为给出明确报错） |
