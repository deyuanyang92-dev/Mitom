# Transcript evidence from RNA-seq (mito_rna.py, `mitom rna`, 2026-11-10)

问题：用转录组能否给线粒体注释提供证据（ORF 是否表达、基因 3′ 端）？用公开数据复现论文结论来检验。
Question: can RNA-seq give evidence for a mitochondrial annotation (are ORFs expressed, where do genes end)? Checked
by reproducing a published study.

## Data / 数据

*Polydora* cf. *ciliata* 转录组 SRR26830859（Selifanova et al. 2023, *Genome Biol Evol*, PMC10721130：Illumina
NovaSeq，双端 61 bp，随机引物 cDNA，2.78 亿对 reads；用 organise-sra / kingfisher 下载，md5 校验）。参考：MitoM 对 OQ078742
（同一物种的线粒体基因组）的注释（`../orf_2026-11-09/`），其中三个 ORFan 已由 `mito_orf` 写出（orf188、orf353、orf146）。
只比对到线粒体基因组（minimap2 2.28 `-ax sr`），MAPQ ≥ 20：690,537 条 reads 采用（论文 Bowtie2：672,141 条）。

## Results / 结果（`OQ078742_rna_features.tsv`）

| 区段 / feature | 相对深度 / relative depth | 覆盖比例 / covered share | 论文 / the paper |
|---|---|---|---|
| ND4L / ND2（表达最低的蛋白基因 / least-expressed genes） | 0.23 / 0.20 | 0.70 / 0.70 | — |
| orf188 = ORFan-544 | **0.41** | 1.00 | 有转录，水平与 NAD4L、NAD2 相当 / transcribed like NAD4L, NAD2 |
| orf146 = ORFan-427 | **0.37** | 0.78 | 同上 / same |
| orf353 = ORFan-1060 | 0.16 | **0.32** | 只有部分转录：C 端 / only partly: the C-terminal part |
| NCR2（ND3 后 631 nt）/ after ND3 | **0.01** | 0.00 | 第 4 个区段（ND3 下游）无稳定产物 / no stable product |
| 12S / 16S | 125 / 57 | — | — |

orf353 分五段的深度 / depth along orf353 in fifths：8、1、8、**91、63**（第 212–353 位氨基酸有转录；论文：C 端 213–240、
261–282 位高度保守且有 RNA-seq 覆盖）。文库自动判为**反向链特异**（read1 反义，正向比例 0.00）。

PolyA：全部 reads 中只有 10 条带 polyA 尾（随机引物总 RNA、61 bp reads），这套数据给不出 3′ 端证据；工具如实提示。
polyA 位点的识别由单元测试（人工 SAM）检验，需要 polyA 富集的数据再做实测。

- All three conclusions of the paper are reproduced: ORFan-544 and ORFan-427 are transcribed at the level of the
  least-expressed genes, ORFan-1060 only in its C-terminal part, the region after ND3 not at all. The library is
  recognised as reverse-stranded. Poly(A) reads are too few in this random-primed library for 3' ends (reported).

## Notes / 说明

- 覆盖是连续的（多顺反子前体），所以只凭覆盖定不了边界；相对深度、覆盖比例、polyA 位点才是证据。
- 蛋白基因中位深度 < 10 时不下表达结论（50 万对 reads 的子集：比对上 1,242 条，提示 reads 太少）。
- 尚未接入 `mitom run`；目前作为 `mitom rna` 对注释结果单独运行。

## Files / 文件

`run.sh`; `OQ078742_rna_features.tsv`, `OQ078742_rna_polya_sites.tsv`. Work data on cody:
`/media/cody/18T_R1/yzz/_realtest/rnaseq_2026-11-10/` (the 644 MB curl subset of the first attempt was deleted after
the organise-sra download: `DELETED.txt`).
