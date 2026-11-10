# Missing protein genes by profile HMM of the nearest relatives (mito_pcg_hmm.py, 2026-11-10)

问题：atp8、atp6 等高度分化的基因，注释核心（blastx + exonerate，参考 refseq89m）常常漏掉，例如龙介虫科
（Serpulidae）的 atp8：NCBI 上 9 条有注释的龙介虫只有 2 条注释了 atp8。用近缘种的谱 HMM（HMMER）补找，能找回多少，会不会乱报？
Question: highly divergent genes (atp8, atp6) are missed by the core (e.g. the atp8 of tube worms, annotated in 2 of 9
NCBI serpulids). How many does a profile HMM of the nearest relatives recover, and does it invent genes?

## Method / 方法（`mito_pcg_hmm.py`；`mitom run --pcg_hmm yes`，默认，第 6 步补漏中、MITOS2 补漏之后）

只对注释完仍缺的蛋白基因运行（补充，不替代主流程）。Only genes still missing after the core run through it.

1. 参考：`models/pcg_refs.faa.gz`（`scripts/build_pcg_refs.py`：16,073 条 RefSeq + 1,242 条 NCBI 环节动物，140,824 条
   蛋白，每属每基因最多 2 条，带分类谱系）。按谱系共享层级选最近的；层级相同时按记录自身 cox1（无则 cob ……）的 3-mer
   相似度排（同为缨鳃虫目时缨鳃虫科排在须腕动物之前）；无分类信息（FASTA）时只按相似度。
2. 最近 10、30、100 条各建一个谱 HMM（MAFFT + hmmbuild），在没有外显子覆盖的 ORF（≥ 20 个密码子）中 `hmmsearch --max`；
   E 值乘 3（三次尝试，Bonferroni）。
3. 接受：E ≤ 0.01 × 近缘先验比值 min(1, p/(1−p))，p = 最近 30 条已注释记录中有该基因的比例（平滑）；模型覆盖 ≥ 40 %。
   近缘普遍没有的基因（如寄生扁虫、色矛纲线虫的 atp8）需要强得多的证据。
4. 位置：最靠近比对起点的起始密码子（先上游）到终止密码子（与下一基因重叠 ≤ 30 nt 允许）；之后由校正步骤调整。

## Results / 结果

留一属（同属参考全部排除）/ leave-genus-out，189 条记录（每门最多 15 条 + 86 条有 atp8 的缨鳃虫目），逐个基因去掉再找：

| 基因 / gene | 次数 / searched | 位置正确 / correct | 位置错 / wrong place | 未找到 / missed |
|---|---|---|---|---|
| atp8 | 189 | **187** | 0 | 2 |
| atp6 | 188 | 186 | 0 | 2 |
| nad4l / nad6 / nad1 / cox2 | 188–189 each | all | 0 | 0 |
| cox3 / nad2 / nad4 / cob / nad3 | 188–189 each | 186–187 | 0–2 | 0–3 |
| cox1 | 188 | 172 | 11 | 5 |
| nad5 | 187 | 179 | 8 | 0 |
| 全部 / all 13 | 2,448 | **2,412（98.5 %）** | 23 | 13 |

正确的 2,412 个中起点完全一致 1,817、终点完全一致 1,876。cox1 / nad5 位置错的多为带内含子的基因只找到一个外显子（内含子
步骤之后补）。atp8 未找到的 2 个：一个无腔动物；*Spirobranchus* sp. PX970954（留一属后同科只剩 1 条参考，E = 0.03，差一点）。

误报 / false positives：200 条没有注释 atp8 的扁形动物、线虫记录中报出 5 个 atp8。1 个是丝虫 *Mansonella*（色矛纲，无
atp8：误报）；4 个是自由生活扁虫（三肠目 *Parakontikia*、*Amaga*，多肠目 *Hoploplana*、*Planocera*）。参考库里已有 3 种
自由生活扁虫注释了 atp8（*Macrostomum*、*Stenostomum*、*Syndesmis*），三肠目的命中又在同一位置，很可能是 NCBI 漏注的真 atp8。

无注释环节动物 335 条（`mitom run` 注释，无真值；`annelid_unannotated335_pcg_hmm.tsv`）：缺 atp8 的 38 条中找到 29 条——
龙介虫 19/26、小头虫 5/5、才女虫 3/3；atp6 找到 3/9。基因排列一致：*Hydroides* 9 条都在 ND4L 之前（与 NCBI 唯一注释了
atp8 的 *Serpula uschakovi* 相同），*Spirobranchus corniculatus* 5 条都在 ATP6 之前，小头虫在 ND4L / ND4 前，才女虫在 ND6 后。

- Leave-genus-out, 2,448 gene searches: 2,412 at the right place (atp8 187/189). False positives in 200 flatworm /
  nematode records: 1 sure (a filarial nematode); 4 in free-living flatworms, which have atp8 in three annotated
  relatives, most likely real unannotated atp8.
- 335 unannotated annelids: atp8 found in 29 of 38 records missing it (tube worms 19/26), at the same gene neighbours
  within each genus as the one NCBI-annotated serpulid atp8.

### atp8 / atp6 按门 / by phylum（`eval_atp.py`，`atp_by_phylum.tsv`；每门最多 15 条注释了该基因的记录，留一属）

上表的抽样按 Metazoa 下一级分组（Spiralia 等），环节、纽形、腕足、苔藓动物几乎没有覆盖；这里按真正的门抽样，只测 atp8、atp6。
The table above samples by the rank below Metazoa; here atp8 / atp6 are sampled by real phylum.

| 类群 / group | atp8 正确 / correct | atp6 正确 / correct |
|---|---|---|
| Annelida | 15/15 | 15/15 |
| Annelida: Sabellida | 14/15（not found 1） | 14/15（not found 1） |
| Arthropoda | 15/15 | 15/15 |
| Brachiopoda | 4/5（not found 1） | 6/6 |
| Bryozoa | 8/8 | 9/9 |
| Chordata | 15/15 | 15/15 |
| Cnidaria | 15/15 | 15/15 |
| Echinodermata | 14/15（not found 1） | 15/15 |
| Entoprocta | 2/2 | 2/2 |
| Gastrotricha | 1/1 | 1/1 |
| Hemichordata | 6/6 | 6/6 |
| Kinorhyncha | 2/2 | 3/3 |
| Mollusca | 15/15 | 15/15 |
| Nematoda | 13/15（not found 2） | 15/15 |
| Nematomorpha | 2/2 | 2/2 |
| Nemertea | 15/15 | 15/15 |
| Onychophora | 4/4 | 4/4 |
| Phoronida | 1/1 | 1/1 |
| Placozoa | — | 3/3 |
| Platyhelminthes | 0/3（not found 3） | 15/15 |
| Porifera | 15/15 | 15/15 |
| Priapulida | 2/2 | 2/2 |
| Rotifera | — | 2/2 |
| Tardigrada | 3/3 | 3/3 |
| Xenacoelomorpha | 6/7（not found 1） | 10/10 |
| **合计 / total** | **187/196** | **218/219** |

未找到的 atp8 / atp8 missed：扁形动物 3 个（*Macrostomum*、*Stenostomum*、*Syndesmis*：少数有 atp8 的扁虫，近缘大多没有，
先验使门槛变严）、线虫 2 个（旋毛虫 *Trichinella*、鞭虫 *Trichuris*，同理）；*Spirobranchus* sp.、*Lingula*、*Stichopus*、
*Paratomella* 各 1 个。atp6 只漏 1 个（须腕动物 *Siboglinum ekmani*）。The misses are mostly the few flatworms / nematodes
that keep atp8 while their relatives lost it: the prior is the price paid for the low false-positive rate.

### 无注释环节动物中没找回的原因 / Why some annelid genes were not found

335 条中 atp8 未找到 9 条、atp6 未找到 6 条，逐条核查（长度、线状 / 环状、N、其余基因）：

| 原因 / reason | atp8 | atp6 | 记录 / records |
|---|---|---|---|
| 记录不完整或组装问题（线状短片段、UNVERIFIED、缺多个基因）：跳过 / incomplete record: skipped | 3 | 5 | *Lagis*、*Maldane*、*Trochonerilla*、*Trilobodrilus*、*Armandia*（13.7 kb 标 complete 但缺 cox1 / cox2 / cob）；*Hydroides dirampha*（线状、UNVERIFIED、有 N，ND3 与 ND4L 间只有 12 nt） |
| 基因组完整，只缺这个基因，有弱命中（E 0.01–0.3）/ complete, only this gene missing, weak hits | 6 | 0 | 龙介虫 *Protula* 3 种、*Ficopomatus*、*Pileolaria*、*Marifugia*：基因排列与找回 atp8 的 *Hydroides* / *Spirobranchus* 不同，参考库中这些属都没有 atp8 注释（NCBI 的 *Ficopomatus enigmaticus* 也没有）——丢失还是太分化，现有数据定不了 |
| 高度重排的基因组 / rearranged genome | 0 | 1 | 水蛭 *Helobdella europaea*（同时缺 nad5） |

去掉不完整记录后：atp8 找回 29/35（龙介虫 19/25），atp6 3/4。Excluding incomplete records: atp8 29/35 (tube worms
19/25), atp6 3/4. The 6 remaining tube worms need annotated atp8 of their genera (e.g. a colleague's verified genomes).

## How the rule was reached / 规则的由来（每处都由失败例找出）

- v1 → v2：跨原点的命中坐标未折回（cox1 位置错 43 → 11）；加入近缘先验（阴性 3/80 → 1/200）。
- v3：`mito_orf` 的串联重复屏蔽把 3 份 CTT（Leu-Leu-Leu）也屏蔽了，龙介虫 atp8 被打上 X；改为阵列 ≥ 30 nt 才屏蔽
  （ORF 结果几乎不变，见 `../orf_2026-11-09/`）。
- v4：谱系相同时按 accession 排序，须腕动物挤掉了缨鳃虫——改为按 cox1 相似度；50 个密码子的基因由哪套参考找到不确定——
  10 / 30 / 100 三个模型都试、E × 3；HMMER 的快速过滤器在几百个短 ORF 上无必要且会漏掉弱命中——`--max`。atp8 181 → 187。

## Files / 文件

`eval_pcg_hmm.py`, `eval_atp.py`, `run.sh`; `leave_genus_out.tsv` (one row per search), `atp_by_phylum.tsv`,
`annelid_unannotated335_pcg_hmm.tsv`.
Work data on cody: `/media/cody/18T_R1/yzz/_realtest/pcg_hmm_2026-11-10/`.
