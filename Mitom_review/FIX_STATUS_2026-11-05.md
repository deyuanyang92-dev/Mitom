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

---

# 第二轮：三份模型报告（Grok 4.7 / GLM 5.3 / Kimi K3）和 11-03 静态审查的逐条核实

方法：每份报告交给一个子代理，逐条读代码、做最小复现（不改仓库），给出判定和补丁；补丁由 Claude Code 合并、写回归测试、跑全量测试。判定：**已修** / **不成立** / **设计**（有注释、文档或测试锁定该行为）/ **不改**（原因见表）。

## Grok 4.7（算法质疑最多的一份）

| ID | 判定 | 要点 |
|---|---|---|
| 1.1 吸虫→21 | 不成立 | 参考库 32/34 条吸虫记录用 9 |
| 1.2 Chordata→2 | 已修 | `import --clade auto`（第一轮）；显式 `--clade Chordata` 也改用记录自身的密码（GLM 3.2）；`Urochordata` 同义名→13 |
| 1.3 / 1.4 单独运行 curate 时固定用 code 5 | 已修 | 无 `--genetic_code` 且 CDS 无 `transl_table` 时，按分类取密码（`mito_tbl.record_code`），并写入 `/transl_table`；验证用的 `.fsa` 带 `[mgcode=]`。实测 NC_000933：14/14 条翻译与 NCBI 一致（原来 13 条不一致） |
| 1.6 code 14 | 不成立 | 已有的 14 会保留 |
| 1.7 半索动物 | 已修（范围窄） | 羽鳃类 Rhabdopleuridae→24、Cephalodiscidae→33；肠鳃纲仍为 9 |
| 1.8 trnL/trnS 编号 | 设计 | 按 MITOS2 规则（655a75c 有文档和测试） |
| 1.9 裸 `atp` | 已修，比报告说的更严重 | N1 把 `atp` 改名为 ATP6 后，D1 当作重复的 ATP6 **删掉了**（实测 atp8 被删）。现从同义词表中去掉裸 `atp` |
| 1.10 单碱基 transl_except | 不成立 | 单碱基合法的只有 TERM，会被重写 |
| 1.11 `/codon_start` | 已修，范围比报告更大 | `codon_start≠1` 的 CDS 原来会被移起点、当作移码拼接、丢失 translation；现在跳过边界规则并记录 |
| 1.12 plastid | 已修 | `plastid:apicoplast` 不再写成 chloroplast |
| 1.13 下载检索式 | 已修 | 加上 `"complete mitochondrial genome"[Title]`（实测漏掉 406 条后生动物记录） |
| 1.14 Radopholus | 不改 | 参考库无此属记录；gb2tbl 已提示密码为推定 |
| 1.15 "subunit 1" | 已修 | 阿拉伯数字的 COX 产物名可识别（名称快照 190 行变化，全部是这两条规则） |
| 2.1 两个终止→换读码框 | 不改（不确定） | 机制能复现，但这是 NC_086642 的启发式，读码框里有 2 个终止的基因本身已损坏；没有基准测试证据前不改 |
| 2.2 内含子端延伸用外接矩形 | 已修 | 邻居是跨原点 join 或本基因位于邻居的内含子里时，原来会跳过该邻居；现在按每一段判断。报告给的原因（`fs>fe`）不对，负链子项不成立 |
| 2.3 参考内含子 3′ 端 | 设计 | 与代码注释和测试一致 |
| 2.4 负链核 rDNA | 不成立 | blastn 负链命中是 qstart<qend、sstart>send，已有测试锁定 |
| 2.5 tRNA 重叠长度少 1 | 已修 | 只在重叠恰好等于一半时有影响 |
| 2.6 atp8/nad4L 回退 ORF 不能跨原点 | 已修，另有一处 bug | 另外发现负链跨原点时会返回附加拷贝里的坐标 |
| 2.7 `--stop_rule homolog` | 设计 | 非默认选项，docstring 写明 |
| 2.8 ARWEN 反密码子为空 | 已修 | 只补反密码子；模型特征保持不变 |
| 3.1 跨原点切成两段的基因 | 已修（非标准输入） | 两个独立 feature（不是 join）时取 5′ 段 |
| 3.2 内含子 note 被删 | 已修 | 只删当前 CDS 自己的跨原点 note |
| 3.3 RefSeq COMMENT | 设计 | dedup 规则 (2)，有测试 |
| 3.5 第 2 轮覆盖第 1 轮 | 不确定 | 没有证据表明会丢基因 |
| 3.7 FASTA 与线性 mtDNA 类群规则 | 设计 | 文档写明 |
| 3.10 Windows | 不支持原生 Windows | |
| 3.11 去重计数 | 已修 | 截断的记录不再计入重复 |
| 3.12 合并文件非原子 | 已修 | |
| 3.13 跨原点的 `/transl_except` 没有平移 | 不改 | 比较时只检查有无，不影响结果 |

## GLM 5.3

| ID | 判定 | 要点 |
|---|---|---|
| 2.1 O2 负链反了 | **不成立** | 负链按转录方向 `parts[0]` 就是原点前那段；现有测试锁定。**报告建议的改法会保留 3′ 越界并丢掉基因，没有采用** |
| 2.2 tRNA `group()` 只和第一个成员比 | 已修 | 夹在两个预测之间的预测会并入同一组（改变 tRNA 结果，见 A/B） |
| 2.3 输入 tRNA 跨原点取外接矩形 | 已修，比报告说的更严重 | 原来变成 1..16000，会在冲突中挤掉一个真实 tRNA |
| 3.2 显式 `--clade Chordata` | 已修 | GenBank 输入用记录自身的密码（被囊 13、头索 5） |
| 3.3 defline | 部分 | 非 ASCII 字符 table2asn 会报 FATAL，现已计入；`[organism=.]` 占位未处理（需你决定验证是否与 gb2tbl 保持一致） |
| 1.4 locale | 不成立 | LC_ALL=C / POSIX 下 Python 3.11 会自动启用 UTF-8 模式，实测不受影响 |
| 1.5 tarfile | 已修（11-03 A） | |
| 1.6 注释阶段无超时 | 不改（需你决定） | 给 MitoZ/MITOS2 外壳加超时会留下孤儿子进程，需要按进程组 kill |
| 1.7 / 1.8 / 1.9 / 4.1 / 4.2 / 4.4 | 已修 | 异常信息写入日志；样本 logger 关闭；所有 GenBank 写出都改用 `write_gb`；第 2 轮失败的日志保留；exonerate 失败给出警告；snakemake 写 sample.log |
| 4.5 / 4.6 | 设计 | |

## Kimi K3

| ID | 判定 | 要点 |
|---|---|---|
| 1.5 `/pseudo=""` | **已修，影响每次 `mitom run`** | 最终 ID 改写和重定向用的是原始 `SeqIO.write`，会把 `/pseudo` 写成 `/pseudo=""`；全部改为 `write_gb` |
| 1.4 | 已修 | 目录里有 `.gb.gz` 时，`mito_intron.py` 单独运行会崩溃 |
| 1.8 / 1.13 / 1.14 | 已修 | 程序报错信息保留；失败程序的部分输出不再解析；CM 下载先写 `.part`，有超时和完成标记 |
| 1.11 | 已修 | `reverse_complement_record` 也改为丢弃无法变换的 feature 并警告 |
| 2.1 | 设计 | `tests/test_mito_intron.py:495` 锁定 |
| 2.2 | 已修 | 与 Grok 1.3 合并 |
| 2.3 | 部分 | 未标注时按环状是设计；FASTA 头里写了 `topology=linear` 现在会被遵守 |
| 2.4 / 2.5 | 已修 / 只改注释 | |
| 2.6 MITOS2 片段顺序 | 不成立 | MITOS2 2.1.9 先按 qstart 排序再标字母 |
| 2.7 mito_ncr 单独运行 | 已修 | 输入检查、`cli_guard` |
| 3.2 atp9 | 已修 | `--engine mitos2` 不再丢 atp9（import 本来就不丢） |
| 3.3 / 3.4 | 已修 | 内含子步骤翻译以 M 开头；`n^1` 位点不再写成倒置区间 |
| 3.5 | 不成立 | |
| 4.1–4.5 | 已修 | 02 阶段 / 第 2 轮失败的日志保留；`mitom debug` 能识别 02；curate / rescue / rrna_hmm 的失败写成 `[RECORD]` 行；同名 failed_logs 不再互相覆盖；blastx 改走 `run_tool` |
| 1.1/1.2/1.10/1.12/4.2(a)(c) | 原生 Windows 不支持 | |
| 1.9 / 1.15 / 1.17 / 4.6 | 不改 / 不成立 | |

## 11-03 静态审查

| ID | 判定 | 要点 |
|---|---|---|
| 问题 3 | 已修 | `--overwrite_new_body yes` 会原地重写 body 文件，改为原子写入 |
| 问题 4 | 已修 | 文件在运行中变短时，分片循环会卡死（已复现），现在报错 |
| 问题 7 | 已修 | 生产代码里实际的未用 import 是 17 个（不是 198 个），已删除 |
| 新发现 A | 已修，前提有误 | 3.8.17 / 3.9.17 / 3.10.12 / 3.11.4 起已有 filter；更早的小版本上 `trna_ml` 会崩溃，现在统一走安全解压 |
| 新发现 C | 已修 | 3 处死代码删除 |
| 问题 5 / 6 / 8 / 9 | 不成立 | |
| R1 CRLF | 已修 | 新增 `.gitattributes`（代码和脚本固定 LF，tsv 原样保留） |

## 测试中顺序相关的 2 个失败

根因：`test_mito_engines.py::test_setup_downloads_checks_and_writes_the_config` 里，`mito_env.setup()` 会把 `os.environ["MITOS_REFDIR"]` 设成测试用的假参考目录；而对一个原本不存在的变量调用 `monkeypatch.delenv`，pytest 什么也不记录，测试结束后无法恢复。之后的 tRNA 测试在子进程里读到了这个假目录。已修，全量测试 0 失败。

## 第二轮验证

| 检查 | 结果 |
|---|---|
| 全量 `pytest tests` | **743 passed，0 failed**，3 skipped；新增回归测试共 36 个（`tests/test_review_2026_11_05.py`） |
| 名称快照 `tests/data/name_expected.tsv.gz` | 已重新生成；190 行变化，逐行核对，全部来自 Grok 1.9（裸 `atp`）和 1.15（"subunit 1/2/3"）两条规则 |
| 真实数据新旧对比：从 NCBI 下载 49 条记录（内含子金标准中的 25 条：刺胞 16 + 环节/扁盘/海绵；mini 基准 24 条，覆盖 8 个门），裸序列 `mitom run --clade auto`，旧代码 `7cb43d0` 对比新代码 | 49 条全部写出，exit 0；所有 CDS/tRNA/rRNA/intron 的坐标**完全相同（0 处差异）**；文本只多了 2 行 tRNA `/note="anticodon: …"`（Grok 2.8）；NCBI ERROR 新旧都是 0 |
| 说明 | 改变注释结果的几项（Grok 2.2 内含子端、2.6 跨原点 ORF、GLM 2.2 tRNA 分组、GLM 2.3 输入 tRNA）在这 49 条记录上没有触发。完整的内含子金标准（309 个内含子，需要另一台机器上的 `mitom_work/intron_gold_2026-10-29/gold/gb`）尚未重跑 |

## 第三轮：Cursor 复查（`FIX_REVIEW_2026-11-05.md`，2026-11-06）

| ID | 判定 | 处理 |
|---|---|---|
| F1 显式 `--clade Chordata` 混合密码 | 属实 | 记录间密码不同（如脊椎 2 + 被囊 13）时现在报 ERROR，提示用 `--clade auto`（按密码分组注释）或 `--genetic_code`，不再静默全部用 2。`auto_clade` 现在也能读 `.gb.gz`（`import` 在解压前调用它） |
| F2 双层 gzip 文件名 | 属实 | `a.fa.gz.gz` → `a.fa`（结尾的 `.gz` 全部去掉，单文件和目录输入都一样） |
| F3 失败组的 tsv 进总报告 | 属实 | `reports/*.tsv` 只合并成功的组；失败组的日志、.json、failed_logs 仍保留（N3 的原意） |
| F4 tRNA 分组桥接 | 属实，报告中的例子不成立 | 150..230 与 1..70 不重叠，不会桥接；真实情况是相邻的同氨基酸 tRNA（1..70、71..140）中间夹一条 36..105。现在新预测必须与组内每一条都有重叠才加入；GLM 2.2 的夹心情况仍是一组 |
| F5 环上 `n^1` | 属实 | 环状记录写成 `join(n..n,1..1)`；线性记录仍取末端 1 nt |

验证：新增 6 个回归测试；全量 `pytest` **753 passed，0 failed**，3 skipped。
