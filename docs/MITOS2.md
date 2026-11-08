# MITOS2 结果的批量整理 / Batch post-processing of MITOS2 results

MITOS2 是常用的线粒体注释工具，但它不输出 GenBank、不统一方向和起点、不修整基因端点，结果通常要在 Geneious 里逐条手工调整。
MitoM 把这些步骤自动化，适合大批量数据：

MITOS2 is widely used, but it writes no GenBank file, does not orient the genome, and does not curate gene ends, so
results are usually fixed by hand in Geneious. MitoM automates these steps for many genomes at once:

| 步骤 / Step | MITOS2 | MitoM |
|---|---|---|
| GenBank 文件 / GenBank file | 无 / no | `final_annotation/*.gb`（Geneious 可直接打开 / opens in Geneious） |
| 统一方向，cox1 开头 / plus strand, cox1 first | 无 / no | 自动 / automatic |
| 端点与基因名按 NCBI 规范修整 / ends and names curated to NCBI rules | 无 / no | 自动，每处修改有记录 / automatic, every change logged |
| NCBI 提交文件与检查 / NCBI files and table2asn check | 无 / no | `.tbl .fsa .sqn` |
| 环状图 / genome map | 有 / yes | 有 / yes |

示例均为实际运行（2026-11-08）：*Platynereis dumerilii* NC_000931，15,619 bp，遗传密码 5，去掉注释后作为 `Pdum.fasta`。
All examples are real runs (2026-11-08) on NC_000931 without its annotation (`Pdum.fasta`).

---

## 1. 已有 MITOS2 结果 → GenBank / Existing MITOS2 result → GenBank

```bash
mitom import -i mitos2_result/ --fasta Pdum.fasta -o out --genetic_code 5 --clade Annelida-segmented-worms
```

- `-i`：MITOS2 结果目录（含 `result.bed`），或包含多个结果目录的上一级目录 / a MITOS2 result folder, or a folder of them.
- `--fasta`：MITOS2 用的序列，序列名与 `result.bed` 第一列一致 / the genome MITOS2 ran on, same sequence names.
- `--genetic_code`：必填，与 MITOS2 的 `-c` 相同 / required, same as MITOS2 `-c`.
- 只要转换、不修整 / conversion only: 加 / add `--curate no --plot no`.

本例的基因数、NCBI ERROR 和用时见第 5 节的表；每处修整见 `curated/curation_report.tsv`。
Gene counts, NCBI ERROR and time for this example: the table in section 5; every edit: `curated/curation_report.tsv`.

## 2. 由 MitoM 调用 MITOS2 / Run MITOS2 from MitoM

```bash
mitom run -i Pdum.fasta -o out --engine mitos2 --genetic_code 5 --mitos_args="--best"
```

MITOS2 只负责找基因，后续步骤同第 1 节。安装：`mamba env create -f env/mitos.yml`（MITOS2 2.1.10）。这次调用的用时只写在下表。
MITOS2 finds the genes. Everything after that is section 1. Install with `mamba env create -f env/mitos.yml` (MITOS2 2.1.10). The time for this run is only in the table.

| 步骤 / step | 用时 / time |
|---|---|
| `mitom run --engine mitos2`（本例 / this example） | 261 s |

| MITOS2 参数 / option | MitoM |
|---|---|
| `-c` | `--genetic_code` |
| `-r`, `--refdir` | `--refdir`（默认自动找 refseq89m / found automatically） |
| `--linear` | 按拓扑自动设置 / set from the topology |
| 其他 / any other | `--mitos_args="--best --finovl 30"`（必须带 `=` / write it with `=`） |

常用 / useful: `--best`（每个基因只留最好的拷贝 / best copy only）、`--finovl N`（基因间最大重叠，默认 50 / max overlap, default 50）、
`--evalue X`（蛋白 BLAST 阈值，默认 2 / protein BLAST cutoff, default 2）、`--ncbicode`（NCBI 起止密码子 / NCBI start/stop codons）。
完整列表 / full list: `runmitos -h`.

## 3. 输出 / Output

| 文件 / File | 用途 / Use |
|---|---|
| `final_annotation/<id>.gb` | 最终结果：已定向、修整、NCBI 基因名 / final: oriented, curated, NCBI names |
| `final_annotation_all.gb` | 同上合并为一个文件 / the same in one file |
| `original_annotation/<id>.gb` | MITOS2 结果原样转成的 GenBank（坐标不变）/ MITOS2 result as GenBank, coordinates unchanged |
| `original_vs_mitom.tsv` | 每个基因修整前后对比 / each gene before vs after |
| `summary.tsv` | 每条记录：基因数、是否 cox1 开头、NCBI ERROR 数 / per record: gene counts, cox1 first, NCBI errors |
| `reports/curation_report.tsv` | 每处修改：规则、旧/新位置、原因 / every change: rule, old/new location, reason |
| `ncbi_submission/` | `.fsa .tbl`（BankIt）、`.sqn`、`submission_summary.tsv`（ERROR 为 0 可提交 / 0 ERROR = ready） |
| `plots/<id>.png/.pdf` | 环状图 / genome map |
| `reports/`, `run.log` | 各步骤报告、全部命令 / step reports, all commands |

`curation_report.tsv`（实测 / real）:

```text
gene rule old           new           reason
rrnS R1   10826..11624  10828..11617  trimmed to flanking genes and filled to them
rrnL R1   11660..12876  11681..12852  trimmed to flanking genes and filled to them
ND5  P1   7270..8987    7279..8959    start codon nearest the preceding gene without overlap
ND5  P3   7270..8987    7279..8959    abbreviated stop T/TA immediately before downstream tRNA
ND1  P3   13105..14035  13105..14029  abbreviated stop T/TA immediately before downstream tRNA
ND3  P1   14165..14503  14159..14503  start codon nearest the preceding gene without overlap
```

不转入 GenBank 的 MITOS2 特征：OH / OL、内含子（MitoM 另行识别内含子）。
Not carried over: MITOS2 OH / OL and introns (MitoM finds introns itself).

## 4. Geneious 修改后再修整 / Curate after editing in Geneious

```bash
mitom curate -i edited_gb/ -o curated --genetic_code 5            # 修整端点和名称 / curate ends and names
mitom gb2tbl -i curated/ -o ncbi --genetic_code 5 --validate      # NCBI 提交文件 / submission files
```

| 文件 / File | 用途 / Use |
|---|---|
| `curated/<name>.gb` | 修整后的 GenBank / curated GenBank |
| `curated/curation_report.tsv` | 每处修改，不同意的回 Geneious 改回 / every change; revert in Geneious if you disagree |
| `curated/validation_compare.tsv` | table2asn 修整前后 ERROR/WARNING 对比 / errors before vs after |

curate 只改端点和名称，不增删基因。伸入下游 tRNA 不超过 20 nt 的完整终止密码子保留不改；MITOS2 常这样结尾（本例 cox1、atp8、cox3、nad6），请在 Geneious 中核对。
curate changes ends and names only; it never adds or removes genes. A complete stop codon running up to 20 nt into the
next tRNA is kept; MITOS2 often ends genes this way (here cox1, atp8, cox3, nad6); check these in Geneious.

## 5. 本例与 NCBI 比较 / This genome vs NCBI

本例与 NCBI 注释比较 / this genome vs the NCBI annotation:

| | 蛋白基因两端一致 / PCG both ends | rRNA 两端 / rRNA ends | tRNA | NCBI ERROR | 用时 / time |
|---|---|---|---|---|---|
| MITOS2 | 4/13 | 0/2 | 22/22 | 2 | 295 s |
| MITOS2 + MitoM 整理 / curation | 6/13 | 2/2 | 22/22 | 0 | +8 s |
| MitoM 自身注释 / own annotation | 10/13 | 2/2 | 22/22 | 0 | 35 s |

24 条 RefSeq（8 门，MITOS2 2.1.10）中三者都完成的 18 条，整理前后 / 18 of 24 RefSeq records (8 phyla), before and after curation
（[benchmark/results/curation_cmp_2026-10-28](../benchmark/results/curation_cmp_2026-10-28/README.md)）:

| | 蛋白基因两端一致 / PCG exact | rRNA ≤ 3 nt | tRNA 找到 / found | 有 NCBI ERROR 的记录 / records with ERROR | 端点偏 > 30 nt / ends > 30 nt off |
|---|---|---|---|---|---|
| MITOS2 | 50.0 % | 23.5 % | 95.9 % | 18 / 18 | 56 |
| MITOS2 + MitoM 整理 / curation | 56.5 % | 64.7 % | 95.9 % | 0 / 18 | 48 |
| MitoM 自身注释 / own annotation | 83.0 % | 61.8 % | 99.5 % | 0 / 18 | 10 |

整理让结果全部通过 NCBI 检查，并改 rRNA 端点。它补不回 MITOS2 漏掉的基因，蛋白基因端点仍会偏。多类群开发期基准（不是冻结评测）的数字在 [METHODS.md](METHODS.md) 的「完整流程与 MITOS2 的对比」表里，这里不抄。
Curation makes every record pass the NCBI check and moves rRNA ends. It does not add genes MITOS2 missed, and many protein ends stay off. The metazoan development-benchmark numbers (not the frozen evaluation) are in the "full run vs MITOS2" table in [METHODS.md](METHODS.md), not copied here.

“一致”指与 NCBI 写法一致，不等于生物学上正确；MitoM 的规则来自 NCBI 记录。找基因的能力两者接近，差别主要在端点。请人工核对结果。
“Exact” means as NCBI writes it, not proven biologically right; MitoM's rules were learned from NCBI records. Gene
finding is similar; the gain is mainly in gene ends. Check the results by hand.

## 6. 常见问题 / Troubleshooting

| 现象 / Symptom | 解决 / Fix |
|---|---|
| `runmitos.py`: `No such file or directory: 'cmsearch'` | 先激活 MITOS2 环境 / activate the MITOS2 env (`conda activate mitos`) |
| `mitom doctor`: `runmitos not found`（MITOS2 2.1.9） | 已支持 `runmitos.py` / `runmitos.py` now accepted; 或装 / or install 2.1.10 |
| `--mitos_args "--best"`: `expected one argument` | 写成 / write `--mitos_args="--best"` |
| `--mitos_args: -c2 is set by MitoM` | `-c -o -r -R -i --linear --json`（含 `-c2`、`--cod` 等写法）由 MitoM 设置，用 `--genetic_code` 等 / set by MitoM, use `--genetic_code` etc. |
| `MITOS2 results: give --genetic_code` | 加 / add `--genetic_code` |
| `genome sequence 'X' not found` | FASTA 序列名改成与 `result.bed` 一致 / rename to match `result.bed` |
