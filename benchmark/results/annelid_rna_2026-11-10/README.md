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

| 蛋白基因 3′ 端 / PCG 3′ ends (673) | 数目 / n |
|---|---|
| poly(A) 正好在注释末端（±10 nt）/ at the annotated end | 464 |
| 完整终止密码子后 11–30 nt 结束（短 3′ UTR，中间无其他基因）/ after a complete stop (short 3′ UTR) | 15 |
| 末端被修正 / end moved | 4 |
| 序列在记录中重复，深度不可测 / duplicated in the record | 5 |
| 无 poly(A) 证据 / no evidence | 185 |

无证据最多的是 ND4L（52 种中 50 种）：ND4L 与 ND4 共用一个转录本，自己没有 poly(A)，这是正常的；其次 ATP8、ND6
（表达低、短）。rRNA 3′ 端：104 个中 88 个被确认，1 个修正。
Most unconfirmed ends are ND4L (shares its transcript with ND4) and short, weakly expressed genes.

修正 / corrections（`--refine`；每个都与同科已注释蛋白的中位长度一致）:

| 线粒体基因组 / mitogenome | 物种 / species | 基因 / gene | 改动 / change | poly(A) reads | 同科中位 / family median (aa) |
|---|---|---|---|---|---|
| CM084515.1 | *Whitmania pigra*（同一标本） | ND6 | 177 → 152 codons，止于 TA | 18 | 152 |
| OZ404008.1 | *Sclerolinum annulatum* | ND2 | 339 → 331，止于 TA | 107 | 331 |
| OZ204559.1 | *Sabella spallanzanii* | ND1 | 300 → 305，延伸到 T | 17 | 305 |
| OZ204559.1 | *Sabella spallanzanii* | ND3 | 115 + TA（原无终止密码子） | 29 | 115 |
| OZ204559.1 | *Sabella spallanzanii* | rrnL | 3′ 端移 11 nt | 1,978 | — |

*Sabella* 的 ND1、ND3 和 *W. pigra* 的 ND6，MitoM 原本都标为 "stop codon not determined"：转录组正好补上了 MitoM
定不了的末端。MitoM had left three of these ends open ("stop codon not determined"); the RNA-seq closes them.

## Outliers / 异常值分析（如何得到上面的规则）

| 异常 / outlier | 例子 / examples | 原因 / cause | 处理 / rule |
|---|---|---|---|
| 某基因几乎没有 reads（邻近基因很高） | *Themiste lageniformis* OZ077499 的 ATP8、ND6（各出现两次） | 序列在记录中重复（95–100 % 的 25-mer 出现两次），reads MAPQ 0 被过滤 | `duplicated_share` > 0.5：不判断表达，注释加说明 |
| 其他低表达基因（ND6、ND3、ATP8、ND4L，7–226×） | *Pterobdella*、*Branchellion*、*Siboglinum* … | 无重复，真实低表达 + 短片段建库偏差 | 不是注释错误；短 ORF 要与短基因比 |
| poly(A) 在注释末端后 11–30 nt | 27 个，全部有完整终止密码子 TAA / TAG | 短 3′ UTR | 中间无其他基因、≥ 3 reads：算确认 |
| poly(A) 在注释末端之后、基因没有终止密码子 | *Sabella* ND1（GTT 结尾）、ND3 | MitoM 未定末端 | 延伸到留下框内 T / TA 的位点（中间无终止密码子、≥ 10 reads） |
| poly(A) 在注释末端之内 | *Travisia* COX3（原末端已有 3 reads） | 降解 / 加工中间体 | 已确认的末端不移动 |
| rRNA 内部强位点 | *Sabella* rrnL −11 nt（95 % reads）；*Protula*、*Arenicola* rrnL −190 至 −330 nt | 前者为真 3′ 端；后者不明 | 30 nt 内、占 ≥ 90 % 才移动；远处的只报告 |
| 同科一致的高表达（如管栖须腕动物 ND2 17–34 倍） | *Lamellibrachia*、*Oasisia*、*Ridgeia* | 生物学（科内一致） | 不处理 |
| 同项目 reads 差异 > 2 % | 两种蚯蚓 | 污染 / 标错 / 隐存种 | FLAG，不用 |

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
