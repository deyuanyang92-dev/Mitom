# 内部数据集与模型制备 / Internal datasets

安装时 `mitom setup` 下载 refseq89m 并准备 tRNA 模型。`models/` 另有 MitoM 自建库。日常用户不必重建。维护者按下面命令更新。

`mitom setup` downloads refseq89m and prepares the tRNA models at install time. `models/` also holds MitoM's own libraries. Ordinary users need not rebuild them; maintainers update them with the commands below.

Numbers and eval settings: [METHODS.md](METHODS.md). This page only says how to rebuild files.

## 总览 / Overview

| 文件 / File | 用途 / Purpose | 制备脚本 / Built by |
|---|---|---|
| refseq89m（`mitom setup`） | 核心注释的参考蛋白、tRNA / rRNA 模型 / reference proteins, tRNA / rRNA models for the core | MITOS2 官方包 / official MITOS2 package |
| `models/ref_taxa.tsv.gz` | `--clade auto`：参考基因组的门与遗传密码 / `--clade auto`: phylum and genetic code of reference genomes | `scripts/build_ref_taxa.py` |
| `models/pcg_refs.faa.gz` | 漏基因时按近缘建 HMM / HMMs from relatives for missed genes | `scripts/build_pcg_refs.py` |
| `models/orf_library.fna.gz` | 未知 ORF 的同源 / homologues for unknown ORFs | `scripts/build_orf_library.py` |
| `models/orf_families.hmm.gz` | ORF 家族命名（Pfam 等）/ ORF family names (Pfam etc.) | 随附；无单独重建脚本（`mito_orf.py` 直接读取）/ shipped; no rebuild script (read directly by `mito_orf.py`) |
| `models/intron_sites.*` | 内含子位点迁移 / intron site transfer | 金标准核实后打包（见 METHODS 内含子节）/ packed after gold-standard verification (METHODS intron section) |
| `models/rrna_hmm/` | 按门的 12S / 16S HMM / per-phylum 12S / 16S HMMs | `scripts/fetch_rrna_refs.py` + `mito_rrna_hmm.py build` |
| `models/start_model/` | 起始密码子选择 / start-codon choice | `benchmark/train_start_model.py` |
| `models/nuc/` | 核 rDNA / 组蛋白参考 / nuclear rDNA / histone references | `mito_nuc` 的建库子命令（METHODS 核基因节）/ build subcommand of `mito_nuc` (METHODS nuclear section) |
| `models/mfannot_proteins/` | 内含子从头识别的蛋白库 / protein library for de novo intron detection | MFannot 数据（GPL）/ MFannot data (GPL) |
| `models/rfam_group_I_II_introns.cm` | 内含子类型 / intron types | Rfam RF00028 / RF00029 |

## 1. 参考分类表 `ref_taxa.tsv.gz` / Reference taxonomy table

**输入 / Input**

- 本机已装好的 refseq89m（`--refdir`）
  refseq89m installed on this machine (`--refdir`)
- 一份含分类信息的 RefSeq 线粒体 GenBank，例如：
  a RefSeq mitochondrial GenBank file with taxonomy, for example:

```bash
mitom fetch --taxon Metazoa --refseq_only -o refpool --email you@example.com
```

**命令 / Command**

```bash
python3 scripts/build_ref_taxa.py \
  --refdir ~/mitos_refs/refseq89m \
  --gb refpool/Metazoa.gb \
  -o models/ref_taxa.tsv.gz
```

**输出 / Output**

- `models/ref_taxa.tsv.gz`：每条参考登录号一行，含谱系与 `transl_table`
  one row per reference accession, with lineage and `transl_table`
- 屏幕列出 refseq89m 有、但 GenBank 文件里没有的登录号（这些命中不计入自动分门）
  the screen lists accessions in refseq89m but not in the GenBank file (hits on them are not counted in automatic phylum assignment)

## 2. 蛋白基因参考 `pcg_refs.faa.gz` / Protein-gene references

**输入 / Input**：已注释的动物线粒体 GenBank（可用 ORF 建库时下的 RefSeq 包，再加需要的类群）。
Annotated animal mitochondrial GenBank files (the RefSeq set downloaded for the ORF library, plus any clades you need).

```bash
python3 scripts/build_pcg_refs.py refpool/Metazoa.gb annelida.gb \
  --out models/pcg_refs.faa.gz --per_genus 2
```

**输出 / Output**：`models/pcg_refs.faa.gz`。表头为 `>登录号|基因|遗传密码|谱系`。每个属每个基因最多 `--per_genus` 条。
`models/pcg_refs.faa.gz`. Header: `>accession|gene|genetic code|lineage`. At most `--per_genus` sequences per gene per genus.

## 3. ORF 库 `orf_library.fna.gz` / ORF library

从 NCBI 拉 RefSeq 动物线粒体全基因组，抽出非标准 CDS，以及未注释区段里 ≥ 40 个密码子的开放读框。
Fetches RefSeq animal mitochondrial genomes from NCBI and extracts non-standard CDSs and open reading frames of at least 40 codons in unannotated stretches.

```bash
python3 scripts/build_orf_library.py work_orf --out models/orf_library.fna.gz
```

**输出 / Output**

- `models/orf_library.fna.gz`
- `work_orf/chunk_*.gb`（可断点续跑 / resumable）

需要网络（NCBI E-utilities）。库是核苷酸（`.fna.gz`），不是蛋白 `.faa`。
Needs network access (NCBI E-utilities). The library is nucleotide (`.fna.gz`), not protein `.faa`.

## 4. rRNA HMM 按门参考 / Per-phylum rRNA HMM references

**下载参考记录 / Download reference records**

```bash
python3 scripts/fetch_rrna_refs.py -o refs_rrna --all-metazoa --per-phylum 150 --seed 7
# 或只要几个门：
python3 scripts/fetch_rrna_refs.py -o refs_rrna --phyla Annelida Mollusca
```

可选 `--set benchmark/sets/mito_bench_2026-10.tsv`，把评测集登录号排除出建库集。
Optional `--set benchmark/sets/mito_bench_2026-10.tsv` keeps the benchmark accessions out of the build set.

**建库 / Build**

```bash
python3 mito_rrna_hmm.py build --gb refs_rrna/*.gb -o models/rrna_hmm
```

**输出 / Output**

- 下载目录（如 `refs_rrna/`）：`<门>.gb`、`refs_<门>.txt`、`fetch_log.tsv`
  download folder (e.g. `refs_rrna/`): `<phylum>.gb`, `refs_<phylum>.txt`, `fetch_log.tsv`
- 建库目录 `models/rrna_hmm/`：`<门>.rrnS.hmm`、`<门>.rrnL.hmm`、`manifest.tsv`
  build folder `models/rrna_hmm/`: `<phylum>.rrnS.hmm`, `<phylum>.rrnL.hmm`, `manifest.tsv`

## 5. 起始密码子模型 / Start-codon model

`--set` 与 `--out` 必填。`--set` 为 `基准表.tsv:该次 leave-out 运行根目录`。
`--set` and `--out` are required. `--set` is `benchmark table.tsv:root folder of that leave-out run`.

```bash
python3 benchmark/train_start_model.py \
  --set benchmark/sets/mito_bench_2026-10.tsv:/path/to/mitom_lo_root \
  --out models/start_model
# 写出 models/start_model/start_model.json
```

训练与交叉验证口径见 METHODS 校正表 P11 行。这里不复述百分数。
Training and cross-validation setup: METHODS curation table, row P11. Percentages are not repeated here.

## 6. 内含子位点库 / Intron site library

位点库来自逐条核实的 NCBI 带内含子记录（金标准流程见 METHODS 内含子节与 `benchmark/results/`）。维护步骤：
The library comes from NCBI intron-bearing records verified one by one (gold-standard procedure: METHODS intron section and `benchmark/results/`). Maintenance steps:

1. 更新金标准坐标与类型判定。
   Update the gold-standard coordinates and type calls.
2. 按仓库既有打包方式写出 `models/intron_sites.tsv.gz` 与 `models/intron_sites.fa.gz`。
   Write `models/intron_sites.tsv.gz` and `models/intron_sites.fa.gz` the way the repository already packs them.
3. 跑内含子金标准（见 `knowledge/TEST_REGISTRY.md`），再提交。
   Run the intron gold standard (see `knowledge/TEST_REGISTRY.md`), then commit.

不要把被测记录同时放进位点库再报“默认一致”的数字。
Do not leave a tested record in the site library and then report a "default" match number.

## 7. 核基因参考 `models/nuc/` / Nuclear-gene references

```bash
mitom nuc fetch -o nuc_refs          # 下载 NCBI 参考（可选）
mitom nuc build --gb nuc_refs/rdna.gb nuc_refs/histone.gb -o models/nuc
```

可追加自己的注释 GenBank。准确度表在 METHODS 核基因节。边界 QC 在 `models/nuc/rdna_qc.tsv`。
You can add your own annotated GenBank files. The accuracy table is in the METHODS nuclear section. Boundary QC is in `models/nuc/rdna_qc.tsv`.

## 8. 用户侧数据（不是内部库） / User-side data (not internal libraries)

| 目的 / Purpose | 命令 / Command |
|---|---|
| 拉某类群线粒体全基因组 / Fetch mitogenomes of a clade | `mitom fetch --taxon Annelida -o fetch_out` |
| 去重后再注释 / Deduplicate, then annotate | `mitom dedup -i fetch_out/Annelida.gb -o dedup_out` |
| 给 ORF 步骤额外近缘参考 / Extra close references for the ORF step | `mitom run ... --orf_ref near.gb` |
| 给内含子步骤额外近缘参考 / Extra close references for the intron step | `mitom run ... --intron_ref_gb near.gb` |

## 检查 / Checks

改库之后至少跑：
After changing a library, run at least:

```bash
mitom doctor
pytest -q tests/test_structure.py tests/test_step_argv.py
# 相关领域再跑 knowledge/TEST_REGISTRY.md 里对应行，不要重复未改动的大评测
```
