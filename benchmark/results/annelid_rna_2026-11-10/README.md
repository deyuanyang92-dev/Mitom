# Annelid genome-project mitogenomes annotated with their own RNA-seq (2026-11-10/11, cody)

问题：NCBI 上环节动物全基因组项目（多为 DToL / ERGA 长读长组装）的线粒体基因组，几乎都没有注释。用 MitoM 注释后，
再用同一基因组项目的转录组核对、修正，能得到一批多准确的注释？
Question: annelid genome-project mitogenomes (long-read assemblies, mostly unannotated): annotate with MitoM, then check
and correct with the same project's RNA-seq.

## Data / 数据（每一对的全部登录号在 `datasets.tsv`）

1. 清单：`ncbi-genome-inventory Annelida`（用户工具）。配对只用它推导的关联：线粒体基因组 `relation_derived`
   （in_assembly / same_biosample）对应 RNA run 的 `linked_assemblies_derived` 与 `specimen_derived`
   （`select_pairs.py`）。不按物种名配对。候选 run 全部列在 `candidates.tsv`。
   Pairs only through the inventory's links (`select_pairs.py`); never by species name.
2. 每个线粒体基因组选一个 run：同一标本 > 同项目不同标本 > 同项目、标本未知；再选 PolyA 建库、双端、数据量大的。
3. 获取：`organise-sra locate` 先查本机与局域网服务器（cody + deyuan112），已有的从 deyuan112 取前 800 万对
   （局域网 100 Mbit/s），其余 `organise-sra download` 经 direct-download 直连下载（整个 run）。来源写在
   `results.tsv` 的 `reads_source`。
4. 注释：`mitom run --clade Annelida-segmented-worms`（MitoM c1d29e5 之后的版本）；转录组证据：`mito_rna.py`
   （bwa mem，`../rna_mappers_2026-11-10/`）；写回注释：`mito_rna.py --refine --run <run>`。

| 配对级别 / link level | 对数 / pairs | reads 与线粒体基因组差异中位数 / median edit rate | 最大 / max | 标记 / flagged |
|---|---|---|---|---|
| 同一标本 / same specimen | 30 | 0.35 % | 1.08 % | 0 |
| 同项目不同标本 / same project, other specimen | 22 | 0.45 % | 7.69 % | 2 |
| 同项目、标本未知 / same project, specimen unknown | 2 | 0.81 % | 0.81 % | 0 |

## Species check / 物种核对

同一基因组项目的转录组与线粒体基因组应是同一物种：reads 差异 > 2 % 标为 FLAG（污染或样品标错，也可能是隐存种），
其表达和 3′ 端不使用。标记的 2 对都是蚯蚓（正蚓科）、同项目不同标本：*Aporrectodea caliginosa* OZ374031 /
ERR13621440、*Lumbricus terrestris* OX457054 / ERR12642427（差异见 `results.tsv`）。这两个种都是已知的隐存种复合体，
需要分类学判断。
Same project should mean same species: > 2 % divergence is flagged (contamination, mislabelled sample, or a cryptic
species); the two flagged pairs are earthworms of known species complexes.

## Results / 结果（52 对通过物种核对，52 种 / 52 pairs pass）

| 建库 / library | 对数 / pairs | 蛋白基因 3′ 端被 poly(A) 确认 / PCG 3′ ends confirmed |
|---|---|---|
| PolyA | 47 | 467/608 |
| RANDOM | 4 | 38/52 |
| cDNA_oligo_dT | 1 | 5/13 |
| 合计 / total | 52 | 510/673 |

3′ 端修正（`--refine`；原注释末端无 poly(A)、内部有强位点、截断后留下框内不完整终止 T / TA）/ corrected 3′ ends:

| 线粒体基因组 / mitogenome | 物种 / species | run | 基因 / gene | 密码子 / codons | poly(A) reads | 近缘同科中位长度 / family median (aa) |
|---|---|---|---|---|---|---|
| CM084515.1 | *Whitmania pigra*（同一标本 / same specimen） | SRR26513850 | ND6 | 177 → 152（原注释 "stop codon not determined"） | 18 | 152 |
| OZ404008.1 | *Sclerolinum annulatum* | ERR14224601 | ND2 | 339 → 331 | 107 | 331 |

- 规则如何得来：初版也把 *Travisia forbesii* COX3 截短 6 个密码子（内部位点 11 reads），会切掉保守的 C 端
  …YWWGS；该基因原注释末端已有 3 条 poly(A) reads。于是加了“已被确认的末端不移动”，之后只剩上表 2 个，长度都等于
  同科中位数。The first version also cut Travisia COX3 (losing the conserved ...YWWGS); a confirmed end is now never
  moved.
- 内部引发位点（基因组 A 串上的 oligo-dT）标出不计：52 对共 204 个。Internal-priming sites flagged: 204.
- 表达：短基因测得少是建库长度偏差（ATP8 最短，相对深度最低），不说明注释错；判断短 ORF 要和 ATP8 / ND4L / ND3 比。
  *Glycera fallax* OZ210882 的两个大 ORF 都有转录（`results.tsv` 的 `orfs` 列）。
  Short genes get few reads (library length bias); ORFs of Glycera fallax are transcribed.

## Files / 文件

- `select_pairs.py`, `run_rna.sh`, `finish.sh`, `summarise.py`
- `datasets.tsv`：每对的全部登录号（线粒体基因组、BioSample、组装、RNA run、Experiment、BioProject、BioSample、配对级别）
- `candidates.tsv`：所有关联的 RNA run
- `results.tsv`：每对一行，含物种核对、3′ 端、修正、ORF、reads 来源
- `refined_gb/`：52 条通过核对的注释（MitoM + 转录组证据，note 中写明 run）。Refined annotations of the 52 pairs.

Work data on cody: `/media/cody/18T_R1/yzz/_realtest/annelid_rna_2026-11-10/`（`rna_out_v1_species_match/`：
早先按物种名配对的结果，已作废 / an earlier species-name pairing, superseded）。

## Run errors / 运行中的问题

| 问题 / problem | 原因 / cause | 处理 / fix |
|---|---|---|
| 两对 reads 差 8 % | 我自写脚本按物种名配对，绕开了 inventory 的关联 | `select_pairs.py` 只用 inventory 关联；ncbi-genome-inventory issue #2 建议直接输出配对表 |
| GoaT 0 个在测物种 | 经 Clash 代理连接被关闭（直连正常） | ncbi-genome-inventory issue #1 |
| locate 找到 0 个 | `find -xdev` 跳过了各数据盘（独立挂载） | organise_sra issue #4（分支 feat/locate） |
| 局域网拉取只拿到 1 个文件 | `while read` 循环里的 ssh 读走了 stdin | `ssh -n` |
| 等待脚本不结束 | `pgrep -f` 匹配到启动它的命令本身 | 去掉等待，下载完成后再运行 |
