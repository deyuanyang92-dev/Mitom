# 跟做示例 / Worked examples

每个例子固定四块：**输入 → 命令 → 先看的输出 → 对照**。数据在 `tests/data/`。请先完成 [USAGE 安装](USAGE_zh.md#2-安装)。

Each example has four fixed parts: **input → command → what to look at first → reference**. Data are in `tests/data/`. Finish the [USAGE install](USAGE_zh.md#2-安装) first.

金标准目录 `tests/data/golden/*` 多用经典布局（`final_annotation/`、`reports/`）；下面默认写出分组布局（`annotation/`、`qc/`）。对照金标准时加 `--output_layout classic`，或只比注释内容。

The golden folders `tests/data/golden/*` mostly use the classic layout (`final_annotation/`, `reports/`); the examples below show the grouped layout (`annotation/`, `qc/`). To compare with a golden folder, add `--output_layout classic`, or compare only the annotation content.

数字只在 [BENCHMARK.md](BENCHMARK.md) 与 [METHODS.md](METHODS.md)。
Numbers are only in [BENCHMARK.md](BENCHMARK.md) and [METHODS.md](METHODS.md).

---

## 例 1 · GenBank → 完整注释（环节动物） / Example 1: GenBank to full annotation (annelid)

| | |
|---|---|
| **输入 / Input** | `tests/data/NC_000931.gb`（RefSeq，带分类 / with taxonomy） |
| **期望 / Expect** | 蛋白基因 + rRNA + tRNA 齐全；`summary.tsv` 里 NCBI ERROR = 0；圈图与提交文件齐全 / all protein genes, rRNAs, tRNAs; NCBI ERROR = 0 in `summary.tsv`; maps and submission files present |

```bash
mitom run -i tests/data/NC_000931.gb -o out_annelid
```

| 先看 / Look at | 内容 / Contents |
|---|---|
| `out_annelid/summary.tsv` | 基因数、是否 cox1 起始、NCBI ERROR / gene count, cox1 start, NCBI ERROR |
| `out_annelid/annotation/*.gb` | 最终注释 / final annotation |
| `out_annelid/submission/` | `.fsa` / `.tbl` / `.sqn` / `.val` |
| `out_annelid/annotation/plots/` | 圈图 PNG / PDF / circular maps |
| `out_annelid/qc/` | 各步报告 / per-step reports |
| `out_annelid/logs/run.log` | 完整日志 / full log |

对照：`tests/data/golden/run_annelid_gb/`。
Reference: `tests/data/golden/run_annelid_gb/`.

---

## 例 2 · 无分类 FASTA（自动判门） / Example 2: FASTA without taxonomy (phylum auto-detected)

| | |
|---|---|
| **输入 / Input** | 由 GenBank 转出的裸 FASTA / plain FASTA converted from GenBank |
| **期望 / Expect** | `qc/clade_decision.tsv` 写出门与遗传密码；注释可完成 / `qc/clade_decision.tsv` gives phylum and genetic code; annotation completes |

```bash
mitom gb2fasta -i tests/data/NC_010217.gb -o fasta_in --split
mitom run -i fasta_in -o out_fasta
```

| 先看 / Look at | 内容 / Contents |
|---|---|
| `out_fasta/qc/clade_decision.tsv` | 门 / 密码怎么判的 / how phylum / code were decided |
| `out_fasta/annotation/` | 最终 GenBank / final GenBank |

对照：`tests/data/golden/run_fasta_blind/`。小门或片段很短时请手写 `--clade`（见 [BENCHMARK § 类群](BENCHMARK.md#2-类群与遗传密码--clade-auto)）。
Reference: `tests/data/golden/run_fasta_blind/`. For a small phylum or a very short fragment, set `--clade` by hand (see [BENCHMARK § clade](BENCHMARK.md#2-类群与遗传密码--clade-auto)).

---

## 例 3 · 多门混合一次跑 / Example 3: several phyla in one run

| | |
|---|---|
| **输入 / Input** | 两个不同门的 GenBank 放同一目录 / two GenBank files of different phyla in one folder |
| **期望 / Expect** | 终端按组分别注释；结果合并进同一个 `-o` / groups are annotated separately; results merged into one `-o` |

```bash
mkdir -p two_phyla
cp tests/data/NC_001636.gb tests/data/NC_008164.gb two_phyla/
mitom run -i two_phyla -o out_two
```

先看：`out_two/summary.tsv`（两条都有结果）、`out_two/qc/clade_decision.tsv`。
Look at: `out_two/summary.tsv` (both records have results) and `out_two/qc/clade_decision.tsv`.

对照：`tests/data/golden/run_two_phyla/`。
Reference: `tests/data/golden/run_two_phyla/`.

`annotate` / `import` 一次只能一个门；多门请用 `run`。
`annotate` / `import` take one phylum at a time; use `run` for several.

---

## 例 4 · 刺胞内含子 / Example 4: cnidarian introns

| | |
|---|---|
| **输入 / Input** | `tests/data/NC_006902.gb` |
| **期望 / Expect** | `intron_report.tsv` / `intron_ref_report.tsv` 有候选；注释里有 intron / join / both reports list candidates; the annotation has intron / join |

```bash
mitom run -i tests/data/NC_006902.gb -o out_coral
```

| 先看 / Look at | 内容 / Contents |
|---|---|
| `out_coral/qc/intron_report.tsv` | 从头识别 / de novo detection |
| `out_coral/qc/intron_ref_report.tsv` | 近缘位点迁移 / site transfer from relatives |
| `out_coral/annotation/*.gb` | 是否写出 intron / join / whether intron / join is written |

对照：`tests/data/golden/run_coral_intron/`。评测数字：[BENCHMARK § 内含子](BENCHMARK.md#4-内含子)。
Reference: `tests/data/golden/run_coral_intron/`. Benchmark numbers: [BENCHMARK § introns](BENCHMARK.md#4-内含子).

---

## 例 5 · 只整理已有注释（`import`）
*Example 5: tidy an existing annotation only (`import`)*

| | |
|---|---|
| **输入 / Input** | `tests/data/KY753833_mitoz.gbf`（MitoZ 风格，环节动物 / MitoZ style, annelid） |
| **期望 / Expect** | **保留**原基因；只做定向、校正、提交、作图 / the original genes are **kept**; only orientation, curation, submission files and maps |

```bash
mitom import -i tests/data/KY753833_mitoz.gbf -o out_import \
  --clade Annelida-segmented-worms --genetic_code 5
```

先看：`out_import/annotation/`、`out_import/qc/curation_report.tsv`、`out_import/submission/`。
Look at: `out_import/annotation/`, `out_import/qc/curation_report.tsv`, `out_import/submission/`.

对照：`tests/data/golden/import_mitoz/`（经典布局）。
Reference: `tests/data/golden/import_mitoz/` (classic layout).

| | `run` | `import` |
|---|---|---|
| 原基因 / Original genes | 丢掉，全部重注 / dropped, all re-annotated | **保留** / **kept** |
| 适用 / Use when | 要从序列重做 / redo from the sequence | 已有 MitoZ / MITOS2 / 手注结果 / you have MitoZ / MITOS2 / manual results |

---

## 例 6 · 只要一轮注释（`annotate`） / Example 6: one annotation round only (`annotate`)

| | |
|---|---|
| **输入 / Input** | 同例 1 / same as Example 1 |
| **期望 / Expect** | 工作目录布局（`04.final_gb/` 等），无 round-2 整理 / working-folder layout (`04.final_gb/` etc.), no round-2 tidy-up |

```bash
mitom annotate -i tests/data/NC_000931.gb -o out_ann
```

对照：`tests/data/golden/annotate_gb/`。日常完整流程请用 `run`。
Reference: `tests/data/golden/annotate_gb/`. Use `run` for the normal full workflow.

---

## 例 7 · 从 NCBI 拉取再注释 / Example 7: fetch from NCBI, then annotate

| | |
|---|---|
| **输入 / Input** | 类群名 + 邮箱（或 `NCBI_EMAIL`） / taxon name + e-mail (or `NCBI_EMAIL`) |
| **输出 / Output** | `<名称>.gb`、`fetch_report.tsv`、`fetch.log` / `<name>.gb`, `fetch_report.tsv`, `fetch.log` |

```bash
mitom fetch --taxon Annelida -o fetch_out --email you@example.com --refseq_only
mitom run -i fetch_out/Annelida.gb -o out_fetch
```

`--dry_run` 只看检索式与条数，不下载。同一目录可重复运行，只补新记录。
`--dry_run` only shows the search query and record count, without downloading. The same folder can be re-run; only new records are added.

---

## 例 8 · 核 rDNA / 组蛋白 / Example 8: nuclear rDNA / histones

| | |
|---|---|
| **输入 / Input** | 全基因组组装 FASTA（或已注释 GenBank） / whole-genome assembly FASTA (or annotated GenBank) |
| **期望 / Expect** | `nuc_features.tsv`、`consensus/`、可选 `submission/` / `nuc_features.tsv`, `consensus/`, optional `submission/` |

```bash
mitom nuc annotate -i your_scaffolds.fasta -o nuc_out
mitom nuc annotate -i your_scaffolds.fasta -o nuc_out --genes rdna
```

准确度与用时：[METHODS 核基因节](METHODS.md#核基因标记--nuclear-markers-rdna-and-histones-mito_nucpy)。
Accuracy and run time: [METHODS, nuclear markers section](METHODS.md#核基因标记--nuclear-markers-rdna-and-histones-mito_nucpy).

---

## 例 9 · 圈图两种样式 / Example 9: two map styles

```bash
mitom plot -i out_annelid/annotation -o maps_circos --style circos
mitom plot -i out_annelid/annotation -o maps_mitom --style mitom
```

样例图：`docs/images/map_circos.png`、`map_mitom.png`、`map_introns.png`。有测序数据时加 `--fq1` / `--fq2`。
Example maps: `docs/images/map_circos.png`, `map_mitom.png`, `map_introns.png`. Add `--fq1` / `--fq2` when you have reads.

---

## 例 10 · 排错 / Example 10: troubleshooting

```bash
mitom doctor
mitom debug out_annelid
mitom run -i tests/data/NC_000931.gb -o out_debug --debug
```

更多：[USAGE § 排错](USAGE_zh.md#8-排错)。
More: [USAGE § troubleshooting](USAGE_zh.md#8-排错).
