# Claude Code 修复复查（只记不改）

**日期：** 2026-11-05  
**范围：** `1912f7b` → `eaa8ef5`（`fix/review-2026-11-05` 已合入）  
**对照：** `FIX_STATUS_2026-11-05.md`、当前源码（未改代码）  
**结论：** 多数审查项的修法与 table2asn「有 .val 才算校验过」一致。下面是**新引入或修完仍会踩中的逻辑问题**。

## 仍会错的

| ID | 级别 | 问题 | 会怎样 |
|---|---|---|---|
| F1 | 高 | 显式 `--clade Chordata` 只给**整个输入一个密码** | 脊椎+被囊（或文昌鱼）混在一起时，`auto_clade` 因多种密码报错，这段报错被丢掉，全部退回密码 **2**。被囊仍按脊椎密码注释 |
| F2 | 高 | 双层 gzip 只剥掉一层文件名 | `a.gb.gz.gz` 解开后仍叫 `a.gb.gz`（内容已是普通 GenBank）。目录输入认后缀，这条记录会被当成没有 |
| F3 | 中 | 失败组的 tsv 仍并进总报告 | 注释 gb 不合并，但 `reports/*.tsv` 用了含失败组的列表。总表里会出现没有最终 gb 的行 |
| F4 | 中 | tRNA 分组变成「连坐」 | 两个不重叠的同氨基酸 tRNA，中间夹一条重叠预测，会合成一组，只留一条 |
| F5 | 中 | 环上 `n^1` 位点被夹成末端 1 nt | `clean_record` 把贴着序列末端的零长度位点收成 `n-1..n`，原点另一侧那 1 nt 没有 |

## F1 — 显式 Chordata

- `mitom_pipeline.resolve_tools`：没写 `--genetic_code` 时先令 Chordata=2，再调 `auto_clade`；只有「一个门且一个密码且无报错」才改掉。
- `mitom_inputs.auto_clade`：记录之间密码不同就返回错误字符串（脊椎 2、被囊 13、文昌鱼 5）。
- 调用处只判断 `not err`，**不打印这句错误**。
- `--clade auto` 会按条分组；显式 `--clade Chordata` 不会。
- 测试 `test_explicit_chordata_uses_the_tunicate_code` 只有一条被囊，测不到混合文件。
- `import` 在 `resolve_tools` 之前**不**做 `normalize_input`。`import --clade Chordata -i 被囊.gb.gz` 时 `auto_clade` 不认 `.gz`，同样留在密码 2。`run` 会先解压，单层 `.gb.gz` 没这个问题。

## F2 — 双层 gzip 文件名

- `normalize_input`：内容最多解 3 层 gzip，文件名只在结尾是 `.gz` 时切掉 **3 个字符**。
- `a.fa.gz.gz` → `00.input/a.fa.gz`，里面已是普通 FASTA/GenBank。
- 目录输入靠 `_is_gb_file` / `_is_fasta_file` 看后缀，`.gz` 对不上，整份输入变空。
- 单文件有时还能靠内容猜格式；**目录不行**。这正是 E19 要修的情况。

## F3 — 失败组报告

- `merge_outputs`：失败组不复制 `final_annotation` 等（对）。
- 合并 `reports` 时用的是 `outs_all`（含失败组）。`.tsv` 会 `_concat_tsv` 进总 `reports/`。
- 失败组写到一半的校正表/汇总，会和成功组排在一张表里。

## F4 — tRNA `group()`

- `mito_trna.group`：和组里**任意一条**重叠 ≥ 一半就加入。
- 例：trnF 在 1..70 和 200..270，中间一条 150..230 与两边都够重叠 → 三合一，冲突消解只留一条。
- 测试 `test_trna_call_between_two_others_joins_their_group` 只锁了「本来就是同一条」的夹心，没有两条真基因被桥接的情况。

## F5 — `n^1`

- `mito_tbl.clean_record`：`start==end` 且在 `0..n` 时改成 `[start-1, start+1]`，再夹到序列内。
- 位点在长度 n 上（环的 `n^1`）变成 `n-1..n`，不是跨原点的 2 nt。
- 旧写法会得到越过序列末端的区间；新写法合法，但环上这个位点少了一侧。

## 看过、先不单列成缺陷

- table2asn：不以退出码判崩溃、没有 `.val` 才记失败。与 1.28 的实测一致（有校验错误时单条 exit 0、多条 exit 2，`.val` 仍在）。
- `orf_before` 正链跨原点：落在附加段里的坐标会减回 `n`；跨过接点的由 `location()` 写成 join。现有测试不传 `n`，负链没有把邻居挪进附加段。49 条对比没触发，先不当已证实的错。
- `codon_start` 不是 1 的 CDS 跳过边界规则，但 `translation()` 仍按 `codon_start` 翻译。和注释一致。
- S4–S12（内含子前移、rRNA 软边界、NUMT 等）本轮未做，与 `FIX_STATUS` 相符。
