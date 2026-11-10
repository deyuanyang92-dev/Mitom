# Complete genome-project mitogenomes: does MitoM find every gene? (2026-11-10)

问题：地球生物基因组计划类项目（DToL、ERGA、VGP ……，NCBI 标题 “genome assembly, organelle: mitochondrion”，多为 PacBio
HiFi + MitoHiFi）组装的线粒体基因组完整、成环。在完整的基因组上，MitoM 能否注释出全部 13 个蛋白基因？缺的是漏注、真丢失还是组装问题？
Question: genome-project mitogenomes (DToL, ERGA, VGP ..., mostly PacBio HiFi) are complete and circular. Does MitoM
find all 13 protein genes in them, and is any gene missing a miss, a real loss or an assembly problem?

## Data / 数据

`fetch.py`：NCBI 上全部 6,309 条（2026-11-10；`records.tsv` 逐条列出门、长度、拓扑、项目自己注释的基因数）。其中 5,560 条
**没有任何注释**（只有序列），所以项目注释只能在少数记录上作对照。`make_set.py`：每门最多 25 条成环、≥ 12 kb 的记录，
属尽量不重复，共 221 条、15 个门（`set.tsv`）。`mitom run` 留一（`mito_bench_run.py --tools mitom_lo`：核心参考去掉该记录，
蛋白基因 HMM 的参考去掉同属，`$MITOM_EXCLUDE_GENUS`），MitoM de0bbc0 之后版本。

All 6,309 NCBI records; most have no annotation at all. Test set: up to 25 circular records >= 12 kb per phylum (221,
15 phyla), `mitom run` leave-one-out (the record out of the core's references, its genus out of the HMM references).

## Results / 结果（`summarise.py` → `per_record.tsv`）

| 门 / phylum | 记录 / records | 13 个全 / all 13 | 缺的基因（次数）/ missing genes (n) | HMM 补回 / added by HMM |
|---|---|---|---|---|
| Annelida | 25/25 | 22 | atp8 2, cob 1, cox3 1 | atp6 1, atp8 1 |
| Arthropoda | 25/25 | 25 | — | — |
| Chordata | 25/25 | 25 | — | — |
| Cnidaria | 25/25 | 24 | atp6 1, atp8 1, cox3 1, nad2 1, nad4 1 | — |
| Mollusca | 25/25 | 23 | atp8 2 | atp8 3 |
| Nematoda | 25/25 | 3 | atp8 22, cob 1, cox2 1, nad1 1, nad2 1, nad3 1 | — |
| Porifera | 25/25 | 24 | nad4l 1 | — |
| Echinodermata | 16/16 | 16 | — | atp8 1 |
| Bryozoa | 8/8 | 8 | — | atp8 4 |
| Nematomorpha | 6/6 | 2 | atp8 4 | atp8 2 |
| Platyhelminthes | 6/6 | 0 | atp8 6 | — |
| Nemertea | 4/4 | 4 | — | — |
| Xenacoelomorpha | 3/3 | 1 | atp8 2 | atp8 1 |
| Ctenophora | 2/2 | 2 | — | — |
| Hemichordata | 1/1 | 1 | — | — |

- **全部 13 个基因**：180/221 条。节肢、脊索、棘皮、苔藓、纽形、栉水母、半索动物 87 条全部齐全。
- **atp8 以外的 12 个基因**：除下表 4 条外，217 条 × 12 = 2,604 个基因全部注出；下表中 3 条是生物学原因（多分体、基因丢失），MitoM 正确。
- **缺 atp8**：线虫 22、扁形动物 6——色矛纲线虫与寄生扁虫已知丢失 atp8；其余为高度分化的类群——龙介虫 2（*Serpula
  vermicularis*、*Marifugia*）、线形动物 4、无腔动物 2、软体动物 2（*Limaria*、*Helminthope*）。
- **蛋白基因 HMM 补回**：atp8 13 个（苔藓 4、软体 3、线形 2、棘皮 1、环节 1、无腔 1 …），atp6 1 个。
- **项目自己有注释的 15 条**：与 MitoM 完全一致（没有项目有而 MitoM 缺的，也没有 MitoM 多出的）。

长读长完整组装不该“缺”基因，所以每个意外缺失都与同属 / 同种的其他记录核对（`check/`）：

| 记录 / record | 缺 / missing | 核查 / what the relatives show | 结论 / verdict |
|---|---|---|---|
| OZ017794.1 海绵 *Lycopodina hypogea*（31,099 nt） | nad4l | 另一研究把它的线粒体基因组发表为 3 条染色体（PP657140–PP657142）：1 号 31,099 nt（与 OZ017794 等长）含 12 个基因 + atp9，**nad4l 单独在 291 bp 的 2 号环上**，3 号 140 bp 无基因 | 多分体线粒体基因组，该记录本就不含 nad4l：MitoM 正确 / multipartite mtDNA: correct |
| OZ057410.1 / OZ057411.1 线虫 *Allodiplogaster pararmata* | 410 缺 cob、cox2、nad1、nad2、nad3；411 缺 cox1、cox3 | 两条来自同一个体、同一个 HiFi 组装（nxAllPara2.1）。410 有完整 cox1、cox3、nad4、nad4l、nad6，411 有完整 atp6、cob、cox2、nad1、nad2、nad3、nad5，各带另一部分的残片；合起来 12 个基因全（色矛纲线虫无 atp8；同属 RefSeq 两种也是 12 个） | 两个环的多分体线粒体基因组：MitoM 正确 / two-circle mtDNA: correct |
| OZ199841.1 黏孢子虫 *Kudoa neothunni*（NCBI 归刺胞动物） | atp6、atp8、cox3、nad2、nad4 | RefSeq 同属 *K. hexapunctata*、*K. septempunctata* 缺的正是同样 5 个 | 真实的基因丢失：MitoM 正确 / real gene loss: correct |
| OZ077499.1 星虫 *Themiste lageniformis* | cob、cox3 | ATP8–trnY–trnQ–ND6 在记录中重复两次，中间 1.5 kb 中 cob 只剩 N 端 34 aa（E = 1e-16，模型 9 %）；该属在 NCBI 没有其他完整线粒体基因组 | 未定：组装问题，或另一个环未提交 / unresolved |

Of the four records missing protein genes other than atp8, three are biology (two multipartite mitochondrial
genomes, one gene loss shared by the genus) and MitoM is right; one (*Themiste*) is unresolved.

## Parts of one genome / 多分体线粒体基因组（`mito_multipart.py`，`mitom run` 第 9 步）

上面两个例子（*Lycopodina*、*Allodiplogaster*）说明：逐条看会“缺”基因，其实是一个基因组的几个分子。`mito_multipart.py`
把同一组装（COMMENT “The assembly …”）或同一 BioSample（DBLINK）的记录放在一起：各自缺基因、合起来 ≥ 12 个蛋白基因的，报告为
“parts of one genome”（多分体，或一个基因组的几段）；序列相同的（INSDC 与 RefSeq 拷贝）只算一次。
检验：6,309 条中 70 个组装提交了 ≥ 2 条线粒体记录，多数是 INSDC + RefSeq 同一序列；长度不同的 16 组（37 条，`multipart_groups.tsv`）
用 `mitom run` 注释后分析（`multipart_report.tsv`）：

| 判定 / verdict | 组 / groups | 物种 / species |
|---|---|---|
| 合起来才完整 / complete only together | 6 | 蠼螋 *Labia minor*（7 + 6 个基因）、千足虫 *Choneiulus palmatus*（9 + 4）、线虫 *Allodiplogaster pararmata*（5 + 10，缺 atp8 为色矛纲常态）——各两个环；凤蝶 *Troides* 3 种——各两条线状片段（一个基因组分两段） |
| 一条已全，其余为额外拷贝或片段 / one record complete | 6 | 蚜虫 *Sitobion fragariae*（13 + 5）、蜈蚣 *Scolopendra cretica*（13 + 8 + 5）、小蜂 *Ixodiphagus* 等 |
| 同一序列 / copies | 1 | 甲虫 *Dendarus foraminosus* |
| 合起来也不完整 / incomplete together | 2 | *Lagis koreni*（两条线状片段，8 个基因）；海绵 *Mycale flammula*（4 个小环 6 个基因：atp6、cox2、cox3、nad3 在序列中 tblastn 找不到同源，cob / nad4 只有与 ND2 / ND5 的交叉相似——组装缺环） |

A group whose records each lack genes but hold >= 12 of the 13 together is reported as parts of one genome: 6 of the
16 multi-record genome-project assemblies (two-circle mtDNA of an earwig, a millipede and a nematode; three
butterflies assembled in two pieces). The benchmark helper `write_bare` lost DBLINK (BioSample) when it copied a
record; fixed (MitoM itself keeps it: golden run_issues).

## Files / 文件

`fetch.py`, `make_set.py`, `summarise.py`; `records.tsv` (all 6,309), `set.tsv` (the 221), `per_record.tsv`;
`multipart_groups.tsv`, `multipart_report.tsv`; `check/`: the second *Allodiplogaster* circle (OZ057411.1, `mitom run`) and the three *Lycopodina* chromosomes (PP657140-PP657142).
Work data on cody: `/media/cody/18T_R1/yzz/_realtest/ebp_mito_2026-11-10/`.
