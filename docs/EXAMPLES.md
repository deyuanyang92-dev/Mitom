# 跟做示例

每个例子固定四块：**输入 → 命令 → 先看的输出 → 对照**。数据在 `tests/data/`。请先完成 [USAGE 安装](USAGE_zh.md#2-安装)。

金标准目录 `tests/data/golden/*` 多用经典布局（`final_annotation/`、`reports/`）；下面默认写出分组布局（`annotation/`、`qc/`）。对照金标准时加 `--output_layout classic`，或只比注释内容。

数字只在 [BENCHMARK.md](BENCHMARK.md) 与 [METHODS.md](METHODS.md)。

---

## 例 1 · GenBank → 完整注释（环节动物）

| | |
|---|---|
| **输入** | `tests/data/NC_000931.gb`（RefSeq，带分类） |
| **期望** | 蛋白基因 + rRNA + tRNA 齐全；`summary.tsv` 里 NCBI ERROR = 0；圈图与提交文件齐全 |

```bash
mitom run -i tests/data/NC_000931.gb -o out_annelid
```

| 先看 | 内容 |
|---|---|
| `out_annelid/summary.tsv` | 基因数、是否 cox1 起始、NCBI ERROR |
| `out_annelid/annotation/*.gb` | 最终注释 |
| `out_annelid/submission/` | `.fsa` / `.tbl` / `.sqn` / `.val` |
| `out_annelid/annotation/plots/` | 圈图 PNG / PDF |
| `out_annelid/qc/` | 各步报告 |
| `out_annelid/logs/run.log` | 完整日志 |

对照：`tests/data/golden/run_annelid_gb/`。

---

## 例 2 · 无分类 FASTA（自动判门）

| | |
|---|---|
| **输入** | 由 GenBank 转出的裸 FASTA |
| **期望** | `qc/clade_decision.tsv` 写出门与遗传密码；注释可完成 |

```bash
mitom gb2fasta -i tests/data/NC_010217.gb -o fasta_in --split
mitom run -i fasta_in -o out_fasta
```

| 先看 | 内容 |
|---|---|
| `out_fasta/qc/clade_decision.tsv` | 门 / 密码怎么判的 |
| `out_fasta/annotation/` | 最终 GenBank |

对照：`tests/data/golden/run_fasta_blind/`。小门或片段很短时请手写 `--clade`（见 [BENCHMARK § 类群](BENCHMARK.md#2-类群与遗传密码---clade-auto)）。

---

## 例 3 · 多门混合一次跑

| | |
|---|---|
| **输入** | 两个不同门的 GenBank 放同一目录 |
| **期望** | 终端按组分别注释；结果合并进同一个 `-o` |

```bash
mkdir -p two_phyla
cp tests/data/NC_001636.gb tests/data/NC_008164.gb two_phyla/
mitom run -i two_phyla -o out_two
```

先看：`out_two/summary.tsv`（两条都有结果）、`out_two/qc/clade_decision.tsv`。

对照：`tests/data/golden/run_two_phyla/`。

`annotate` / `import` 一次只能一个门；多门请用 `run`。

---

## 例 4 · 刺胞内含子

| | |
|---|---|
| **输入** | `tests/data/NC_006902.gb` |
| **期望** | `intron_report.tsv` / `intron_ref_report.tsv` 有候选；注释里有 intron / join |

```bash
mitom run -i tests/data/NC_006902.gb -o out_coral
```

| 先看 | 内容 |
|---|---|
| `out_coral/qc/intron_report.tsv` | 从头识别 |
| `out_coral/qc/intron_ref_report.tsv` | 近缘位点迁移 |
| `out_coral/annotation/*.gb` | 是否写出 intron / join |

对照：`tests/data/golden/run_coral_intron/`。评测数字：[BENCHMARK § 内含子](BENCHMARK.md#4-内含子)。

---

## 例 5 · 只整理已有注释（`import`）

| | |
|---|---|
| **输入** | `tests/data/KY753833_mitoz.gbf`（MitoZ 风格，环节动物） |
| **期望** | **保留**原基因；只做定向、校正、提交、作图 |

```bash
mitom import -i tests/data/KY753833_mitoz.gbf -o out_import \
  --clade Annelida-segmented-worms --genetic_code 5
```

先看：`out_import/annotation/`、`out_import/qc/curation_report.tsv`、`out_import/submission/`。

对照：`tests/data/golden/import_mitoz/`（经典布局）。

| | `run` | `import` |
|---|---|---|
| 原基因 | 丢掉，全部重注 | **保留** |
| 适用 | 要从序列重做 | 已有 MitoZ / MITOS2 / 手注结果 |

---

## 例 6 · 只要一轮注释（`annotate`）

| | |
|---|---|
| **输入** | 同例 1 |
| **期望** | 工作目录布局（`04.final_gb/` 等），无 round-2 整理 |

```bash
mitom annotate -i tests/data/NC_000931.gb -o out_ann
```

对照：`tests/data/golden/annotate_gb/`。日常完整流程请用 `run`。

---

## 例 7 · 从 NCBI 拉取再注释

| | |
|---|---|
| **输入** | 类群名 + 邮箱（或 `NCBI_EMAIL`） |
| **输出** | `<名称>.gb`、`fetch_report.tsv`、`fetch.log` |

```bash
mitom fetch --taxon Annelida -o fetch_out --email you@example.com --refseq_only
mitom run -i fetch_out/Annelida.gb -o out_fetch
```

`--dry_run` 只看检索式与条数，不下载。同一目录可重复运行，只补新记录。

---

## 例 8 · 核 rDNA / 组蛋白

| | |
|---|---|
| **输入** | 全基因组组装 FASTA（或已注释 GenBank） |
| **期望** | `nuc_features.tsv`、`consensus/`、可选 `submission/` |

```bash
mitom nuc annotate -i your_scaffolds.fasta -o nuc_out
mitom nuc annotate -i your_scaffolds.fasta -o nuc_out --genes rdna
```

准确度与用时：[METHODS 核基因节](METHODS.md#核基因标记--nuclear-markers-rdna-and-histones-mito_nucpy)。

---

## 例 9 · 圈图两种样式

```bash
mitom plot -i out_annelid/annotation -o maps_circos --style circos
mitom plot -i out_annelid/annotation -o maps_mitom --style mitom
```

样例图：`docs/images/map_circos.png`、`map_mitom.png`、`map_introns.png`。有测序数据时加 `--fq1` / `--fq2`。

---

## 例 10 · 排错

```bash
mitom doctor
mitom debug out_annelid
mitom run -i tests/data/NC_000931.gb -o out_debug --debug
```

更多：[USAGE § 排错](USAGE_zh.md#8-排错)。
