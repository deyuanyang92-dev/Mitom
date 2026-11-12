# 内部数据集与模型制备 / Internal datasets

安装时 `mitom setup` 下载 refseq89m 并准备 tRNA 模型。`models/` 另有 MitoM 自建库。日常用户不必重建。维护者按下面命令更新。

Numbers and eval settings: [METHODS.md](METHODS.md). This page only says how to rebuild files.

## 总览

| 文件 | 用途 | 制备脚本 |
|---|---|---|
| refseq89m（`mitom setup`） | 核心注释的参考蛋白、tRNA / rRNA 模型 | MITOS2 官方包 |
| `models/ref_taxa.tsv.gz` | `--clade auto`：参考基因组的门与遗传密码 | `scripts/build_ref_taxa.py` |
| `models/pcg_refs.faa.gz` | 漏基因时按近缘建 HMM | `scripts/build_pcg_refs.py` |
| `models/orf_library.fna.gz` | 未知 ORF 的同源 | `scripts/build_orf_library.py` |
| `models/orf_families.hmm.gz` | ORF 家族命名（Pfam 等） | 随附；无单独重建脚本（`mito_orf.py` 直接读取） |
| `models/intron_sites.*` | 内含子位点迁移 | 金标准核实后打包（见 METHODS 内含子节） |
| `models/rrna_hmm/` | 按门的 12S / 16S HMM | `scripts/fetch_rrna_refs.py` + `mito_rrna_hmm.py build` |
| `models/start_model/` | 起始密码子选择 | `benchmark/train_start_model.py` |
| `models/nuc/` | 核 rDNA / 组蛋白参考 | `mito_nuc` 的建库子命令（METHODS 核基因节） |
| `models/mfannot_proteins/` | 内含子从头识别的蛋白库 | MFannot 数据（GPL） |
| `models/rfam_group_I_II_introns.cm` | 内含子类型 | Rfam RF00028 / RF00029 |

## 1. 参考分类表 `ref_taxa.tsv.gz`

**输入**

- 本机已装好的 refseq89m（`--refdir`）
- 一份含分类信息的 RefSeq 线粒体 GenBank，例如：

```bash
mitom fetch --taxon Metazoa --refseq_only -o refpool --email you@example.com
```

**命令**

```bash
python3 scripts/build_ref_taxa.py \
  --refdir ~/mitos_refs/refseq89m \
  --gb refpool/Metazoa.gb \
  -o models/ref_taxa.tsv.gz
```

**输出**

- `models/ref_taxa.tsv.gz`：每条参考登录号一行，含谱系与 `transl_table`
- 屏幕列出 refseq89m 有、但 GenBank 文件里没有的登录号（这些命中不计入自动分门）

## 2. 蛋白基因参考 `pcg_refs.faa.gz`

**输入**：已注释的动物线粒体 GenBank（可用 ORF 建库时下的 RefSeq 包，再加需要的类群）。

```bash
python3 scripts/build_pcg_refs.py refpool/Metazoa.gb annelida.gb \
  --out models/pcg_refs.faa.gz --per_genus 2
```

**输出**：`models/pcg_refs.faa.gz`。表头为 `>登录号|基因|遗传密码|谱系`。每个属每个基因最多 `--per_genus` 条。

## 3. ORF 库 `orf_library.fna.gz`

从 NCBI 拉 RefSeq 动物线粒体全基因组，抽出非标准 CDS，以及未注释区段里 ≥ 40 个密码子的开放读框。

```bash
python3 scripts/build_orf_library.py work_orf --out models/orf_library.fna.gz
```

**输出**

- `models/orf_library.fna.gz`
- `work_orf/chunk_*.gb`（可断点续跑）

需要网络（NCBI E-utilities）。库是核苷酸（`.fna.gz`），不是蛋白 `.faa`。

## 4. rRNA HMM 按门参考

**下载参考记录**

```bash
python3 scripts/fetch_rrna_refs.py -o refs_rrna --all-metazoa --per-phylum 150 --seed 7
# 或只要几个门：
python3 scripts/fetch_rrna_refs.py -o refs_rrna --phyla Annelida Mollusca
```

可选 `--set benchmark/sets/mito_bench_2026-10.tsv`，把评测集登录号排除出建库集。

**建库**

```bash
python3 mito_rrna_hmm.py build --gb refs_rrna/*.gb -o models/rrna_hmm
```

**输出**

- 下载目录（如 `refs_rrna/`）：`<门>.gb`、`refs_<门>.txt`、`fetch_log.tsv`
- 建库目录 `models/rrna_hmm/`：`<门>.rrnS.hmm`、`<门>.rrnL.hmm`、`manifest.tsv`

## 5. 起始密码子模型

`--set` 与 `--out` 必填。`--set` 为 `基准表.tsv:该次 leave-out 运行根目录`。

```bash
python3 benchmark/train_start_model.py \
  --set benchmark/sets/mito_bench_2026-10.tsv:/path/to/mitom_lo_root \
  --out models/start_model
# 写出 models/start_model/start_model.json
```

训练与交叉验证口径见 METHODS 校正表 P11 行。这里不复述百分数。

## 6. 内含子位点库

位点库来自逐条核实的 NCBI 带内含子记录（金标准流程见 METHODS 内含子节与 `benchmark/results/`）。维护步骤：

1. 更新金标准坐标与类型判定。
2. 按仓库既有打包方式写出 `models/intron_sites.tsv.gz` 与 `models/intron_sites.fa.gz`。
3. 跑内含子金标准（见 `knowledge/TEST_REGISTRY.md`），再提交。

不要把被测记录同时放进位点库再报“默认一致”的数字。

## 7. 核基因参考 `models/nuc/`

```bash
mitom nuc fetch -o nuc_refs          # 下载 NCBI 参考（可选）
mitom nuc build --gb nuc_refs/rdna.gb nuc_refs/histone.gb -o models/nuc
```

可追加自己的注释 GenBank。准确度表在 METHODS 核基因节。边界 QC 在 `models/nuc/rdna_qc.tsv`。

## 8. 用户侧数据（不是内部库）

| 目的 | 命令 |
|---|---|
| 拉某类群线粒体全基因组 | `mitom fetch --taxon Annelida -o fetch_out` |
| 去重后再注释 | `mitom dedup -i fetch_out/Annelida.gb -o dedup_out` |
| 给 ORF 步骤额外近缘参考 | `mitom run ... --orf_ref near.gb` |
| 给内含子步骤额外近缘参考 | `mitom run ... --intron_ref_gb near.gb` |

## 检查

改库之后至少跑：

```bash
mitom doctor
pytest -q tests/test_structure.py tests/test_step_argv.py
# 相关领域再跑 knowledge/TEST_REGISTRY.md 里对应行，不要重复未改动的大评测
```
