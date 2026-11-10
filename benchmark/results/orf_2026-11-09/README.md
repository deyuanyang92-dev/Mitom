# ORFs of unknown function in animal mitogenomes (mito_orf.py, 2026-11-09)

问题：13 个蛋白基因以外的编码区（MutS、atp9、F-ORF、内含子 ORF、才女虫科的未知 ORF ……），MitoM 能找回多少、会不会乱写？
与只看长度的 ORF 查找（NCBI ORFfinder / Geneious 式）和 MITOS2 比较。方法见 `docs/METHODS.md`“未知功能的 ORF”与
`mito_orf.py` 文件头。
Question: how many protein genes beyond the 13 does MitoM recover, without writing chance ORFs? Compared with ORF
finding by length (ORFfinder / Geneious style) and MITOS2. Method: docs/METHODS.md, mito_orf.py header.

## Data / 数据

| 数据 / set | 记录 / records | 真值 / truth |
|---|---|---|
| 已知 ORF / known ORFs (`make_truth.py`) | 260 条 RefSeq，去掉 NCBI 注释的非标准 CDS 后运行 / RefSeq records run without their non-standard CDS | 342 个 NCBI 注释的非标准 CDS（MutS 125、hypothetical 123、atp9 57、归巢内切酶 10 …）/ annotated CDS |
| 无内含子 / intron-free (`../intron_algo_2026-11-09/negatives_set.tsv`) | 362 RefSeq，26 门 / phyla | 未知；写出的应很少且有依据 / unknown; writes should be few and explained |
| 才女虫科 / Spionidae (NCBI, >= 10 CDS) | 46 | Polydorini trnQ–trnA ORF |
| 无注释环节动物 / unannotated annelids (`mitom run`, ../intron_annelid_2026-11-09/) | 335 | 未知 / unknown |

每条已知 ORF 记录单独运行，参考库中去掉它自己（留一）。Each known-ORF record run alone, itself left out of the library.

## Results / 结果

| 方法 / method | 已知 ORF 写出 / known ORFs written (342) | 只给 WARNING / WARNING only | 无内含子 362 条写出 / written in 362 intron-free |
|---|---|---|---|
| 只看长度 ≥ 75 密码子（ORFfinder 式）/ length >= 75 codons | 326 | — | 357（118 条记录 / records） |
| 只看长度 ≥ 150 密码子 / length >= 150 codons | 214 | — | 71 |
| 第一版 / first version (f2e73e1) | 175 | 69 | 1 |
| 第二版 / second version (d13f41c) | 228 | 47 | 5 |
| tRNA 边界、Pfam、两项各 0.05 / tRNA ends, Pfam, both tests at 0.05 (80adc30) | 233 | 45 | 7 |
| **当前 / now**（串联重复阵列 ≥ 30 nt 才屏蔽、跨原点起点修复 / repeat arrays >= 30 nt, origin fix） | **232** | **46** | **8** |
| MITOS2 | 0（不注释额外 ORF / no extra ORFs） | — | 0 |

已知 ORF 按类别 / known ORFs by kind（写出 / WARNING / 漏 — written / WARNING / missed）：

| atp9 | MutS | 归巢内切酶 / homing endonuclease | hypothetical / ORF（含 F-ORF）| 其他 / other |
|---|---|---|---|---|
| 57 / 0 / 0 | 123 / 0 / 2 | 7 / 0 / 3 | 28 / 41 / 54 | 18 / 4 / 5 |

Pfam 命名（真值中写出的）：MutS 123 个全部为 MutS_V（“MutS-like protein”），归巢内切酶 8 个为 LAGLIDADG_1，3 个 polB（NCBI
有的称 “replication helicase subunit”）为 DNA_pol_B_2。Pfam names among the known ORFs written: all 123 MutS, 8 homing
endonucleases (LAGLIDADG_1), 3 polB (DNA_pol_B_2).

其他数据 / other sets：才女虫科 46 条写出 12——Selifanova 等 2023（*Genome Biol Evol*，PMC10721130）报道的三个 *Polydora*
ORFan 在 4 条记录中全部写出（下表）；无注释环节动物 335 条写出 34（修复跨原点起点后，原先使该步骤崩溃的 OZ210882.1 也能运行）。

| Selifanova et al. 2023（*P.* cf. *ciliata*） | MitoM（OQ078742，单条运行 / run alone） | 结尾 / end |
|---|---|---|
| ORFan-544 | orf188，1783..>2346 | 读进 trnK，3′ 不完整 / runs into trnK, open 3' end |
| ORFan-1060 | orf353，4680..5739 | 终止密码子 / stop codon |
| ORFan-427 | orf146，11977..>12414 | 读进 trnR，3′ 不完整 / runs into trnR, open 3' end |

（单条运行的坐标是 MitoM 自己的注释，与 NCBI 记录上的相差 9 nt；table2asn 对两个 3′ 不完整的 CDS 只给 WARNING。）

- **从第一版改进的地方**（每处都由真值中的失败例找出）：内含子不再排除（归巢内切酶、内含子中的 MutS：0 → 6、103 → 123）；
  多个同源合并计替换数（逐位置碱基种类数，或最小生成树长度的一半，取大者）——近缘种单独不够，多个可以累加；编码读框改用
  三个密码子位置的完整模式（第三位最多、第二位最少），选重叠读框时按此判断（F-ORF 原先选成了反义链）；两项检验用 Fisher
  法合并（每项 ≤ 0.05，合并 ≤ 1e-4）。第一版之前还修过两处：控制区串联重复造成的跨门假同源（屏蔽重复）；蛋白比对造成的
  第三位偏差（改用核苷酸比对计位置）。
- **只看长度不可行**：≥ 75 个密码子时候选里已有 326/342 个真 ORF，但同时在 362 条无内含子记录中写出 357 个 ORF；MitoM 的
  候选步骤就是这一步，第二步按证据判断。
- **362 条无内含子记录中写出的 7 个**：*Aurelia limbata* orf340（Pfam 定为 B 族 DNA 聚合酶，即钵水母线粒体 polB，Kayal 等
  2012）及其旁边的 orf104；*Pseudanodonta complanata* orf79（与 20 种淡水蚌已注释的 F-ORF 同源，NCBI 此条未注释）；玻璃海绵
  *Oopsacas minuta* orf110；*Iheyomytilidicola lauensis* orf148；双壳类 *Saxidomus purpurata* orf90、*Teredo bartschi* orf57。
  都有同源和编码模式证据；双壳类常有额外 ORF，后两个是否真实未知。
- **规则怎么定的**（`rnacode/README.md`）：RNAcode（比较基因组编码潜力）单独用会把很多控制区判为编码；与本方法两项检验
  合用几乎没有增益，未采用。两项检验（保持开放、编码模式）都是必要条件，各取 0.05：非编码序列同时通过约 0.0025，362 条记录
  683 个有同源的候选预计约 2 个，实际 7 个中至少 5 个很可能真实。
- **335 条无注释环节动物写出的 30 个**（多出的是 *Polydora* 的 ORFan-544 / -427 同源）：Polydorini ORF（3 条 *Polydora* 加 *P. caeca*）；*Glycera fallax* orf728 与
  *Podarkeopsis* orf778 互为同源，都在 cox1 的 ~2.4 kb II 类内含子中，应是内含子编码的逆转录酶 / 成熟酶；*Hydroides*
  orf412；*Spirobranchus corniculatus* 5 条同种记录中各 3 个（orf154 / orf170 相邻且互为重复，同源包括基因组内的重复拷贝，
  是否真实不确定）。
- **漏检的 64 个**：多为没有任何同源的类群特有 hypothetical protein（黏孢子虫、扁盘动物、钵水母、栉水母），其真实性本身未
  证实；较长的给 WARNING。密码子使用偏好（与本记录 13 个基因比，GeneMark 式）也测过：在 362 条无内含子记录中只放进 14 个
  随机 ORF 时只能找回 29 % 的真 ORF，区分力太弱，未采用。
- **单个基因组**：NC_061377 单独运行时近缘属不在 RefSeq，orf357 只给 WARNING；`--orf_ref` 给出 OQ078742 后写出
  （`NC_061377_with_orf_ref.tsv`）。MITOS2 在 NC_061377、OQ078742、MW528029 中都没有注释它。

- What changed since the first version, each found from failures on the truth: introns searched (homing endonucleases
  0 -> 6, MutS 103 -> 123); substitutions summed over many homologs (per-position base kinds or half the MST); the
  frame test uses the full codon pattern (position 3 most, 2 least), which also picks the right strand of overlapping
  ORFs; the two tests Fisher-combined.
- Length alone: at >= 75 codons the candidates hold 326/342 true ORFs, but 357 ORFs would be written in 362 intron-free
  records. The 5 MitoM writes there all have homologs and the coding pattern.
- Misses: mostly lineage-specific hypothetical proteins without any homolog; long ones are WARNINGs. Codon usage was
  tested and is too weak to use.

## Files / 文件

`run.sh`; `make_truth.py`, `score_recall.py`, `length_baseline.py`; `recall_truth.tsv`, `recall_scored.tsv`,
`length_baseline.tsv`; `*_orf_report.tsv` (one per set); `NC_061377_with_orf_ref.tsv`;
`OQ078742_single_run_orf_report.tsv` (`mitom run`); `rnacode/` (the RNAcode comparison). Work data on cody:
`/media/cody/18T_R1/yzz/_realtest/orf_2026-11-09/`.
