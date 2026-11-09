# Introns: NCBI animal mitogenomes not in the gold standard (2026-11-09)

Question: does MitoM annotate the introns of every animal mitogenome in NCBI, beyond the 193 records of the intron
gold standard (`../intron_gold_2026-10-29/`, which is also MitoM's intron-site library)?

问题：NCBI 中所有带内含子的动物线粒体基因组，金标准（同时是 MitoM 的内含子位点库）以外的记录，MitoM 能否注释出内含子？

## Data / 数据

NCBI esearch 2026-11-09 (`make_set.py`): `txid33208[Organism:exp] AND mitochondrion[filter] AND intron[FKEY] AND
10000:100000[SLEN]`, every record in `ncbi_records.tsv` with its status:

| 状态 / status | 记录 / records |
|---|---|
| 已在金标准 / in the gold standard | 89 |
| 金标准 NC_ 记录的 INSDC 副本 / INSDC copy of a gold NC_ record | 66 |
| 本次测试的 NC_ 记录的 INSDC 副本 / INSDC copy of a tested NC_ record | 2 |
| 本次测试 / tested (`set.tsv`) | 134：刺胞 Cnidaria 122、苔藓虫 Bryozoa 5、海绵 Porifera 4、脊索 Chordata 3 |

None of the 134 is in the site library. Records with an intron only as a CDS join (no intron feature) are not found by
this query (the gold standard's annelid records are mostly of that kind). Shorter records with intron features
(fragments, < 10 kb) were not tested.

134 条都不在位点库中。只用 CDS join 表示内含子、没有 intron 特征的记录（金标准中的环节动物多为这种）检索不到；< 10 kb 的片段未测。

## Run / 运行

`run.sh`: `mitom run` on the bare record (source feature only; the record's clade and genetic code), 24 records at a
time x 2 threads, MitoM f04d839, mitotools environment, machine cody; scored with `benchmark/eval_introns.py` against
NCBI's intron features (not checked placements). Leave-one-out (`mitom_lo`: the record left out of the core's
references and the site library) gives the same counts as the default (`mitom_native`), as expected for records the
library does not hold. `check_products.py` splices each CDS as NCBI and as MitoM do, translates it and scores it with
blastp against the refseq89m proteins of the gene (`products_native.tsv`; `better` = no internal stop and the higher
bit score).

## Results / 结果

Rule of these runs: the one before 2026-11-09 (`--intron_rule splice`). With the rule that replaced it (continuity
and necessity of the splice, `../intron_algo_2026-11-09/`) the bryozoan introns found rise from 12 to 31 of 50; the
other groups are unchanged.

`mitom_native_summary.tsv` (truth = NCBI intron features):

| 类群 / group | NCBI 内含子 / NCBI introns | 完全一致 / exact | ≤ 6 nt | > 6 nt | 漏检 / missed | 多报 / extra |
|---|---|---|---|---|---|---|
| 刺胞 / Cnidaria | 189 | 119 | 39 | 30 | 1 | 1 |
| 海绵 / Porifera | 6 | 5 | 1 | 0 | 0 | 0 |
| 苔藓虫 / Bryozoa | 50 | 11 | 0 | 1 | 38 | 0 |
| 鸟类 *Goura* nad3（非内含子 / not introns） | 3 | — | — | — | 3 | 0 |

Where MitoM's and NCBI's placements differ in Cnidaria (`products_native.tsv`):

| 位置差 / placement difference | MitoM 产物更好 / MitoM's product better | 相同产物 / same product | NCBI 更好 / NCBI better |
|---|---|---|---|
| ≤ 6 nt（39） | 20 | 17 | 2 |
| > 6 nt（30） | 28 | 2 | 0 |

- **Cnidaria**: all introns found but one. Of the > 6 nt cases, 28 are NCBI errors and 2 give the same protein
  (NCBI's intron feature starts 200 / 721 nt inside the first exon of its own CDS join; MitoM matches the join: OQ740730.1, PX980026.1). The errors: 18 sea-anemone records of 2025-2026
  (BK069890-903, BK071241-243, PV232310.1) end the nad5 intron 36 nt late, which drops 12 conserved codons from nad5; BK068676-82 give
  only a 247-aa piece of nad5. The missed one (OR906318.1, a "trnY-AUA intron" of 594 nt splitting the tRNA across the
  origin, inside the nad5 region) is not a credible intron. The extra one is a cox1 intron NCBI did not annotate
  (BK068680.1: NCBI's cox1 reads through it, blastp covers 49 % of its product; MitoM's spliced cox1 is complete).
  The two "NCBI better" cases are cox1 introns 1-4 nt off (HG423144.1, OL372279.1; 6-9 bits).
- **Porifera**: all found, including the two group II introns of *Acanthella acuta* (OM729606.2, cox1 and cox2).
- **Bryozoa** (Jenkins et al. 2022, Sci Rep, "Unprecedented frequency of mitochondrial introns in colonial
  bilaterians"): the long introns (500-2,200 nt; *Exechonella*, *Parantropora*, PP990758.1) are found at NCBI's
  position. The short ones (240-290 nt, up to 18 per genome in Cupuladriidae: *Cupuladria* MW592987.1,
  *Discoporella* MW592989.1) are missed: the de novo step finds 23 of the 38 missed at or near NCBI's position, but
  they have no group I / II signature and weak splice sites, so they stay low-confidence and are not written
  (`intron_report.tsv`); NCBI's spliced products are the better ones there. No false introns in these records.
- ***Goura*** (pigeons): NCBI's 1-nt nad3 "intron" is the bird nad3 frameshift; MitoM writes none (correct).

- 刺胞动物：除 1 个外全部注出；位置差 > 6 nt 的 30 个中 28 个是 NCBI 写错、2 个蛋白相同（2025–2026 年 18 条海葵记录 nad5 内含子晚结束 36 nt，丢 12 个保守密码子）；漏的 1 个不是可信的内含子；多报的 1 个是 NCBI 漏注的 cox1 内含子。
- 海绵：全部注出（含 *Acanthella acuta* 的 cox1、cox2 II 类内含子）。
- 苔藓虫：长内含子注出；Cupuladriidae 的短内含子（240–290 nt）大多漏掉——从头识别找到了其中大多数，但没有 I/II 类特征，置信度低而未写入。这是目前唯一明显的缺口。

## Files / 文件

`ncbi_records.tsv`, `set.tsv`, `mitom_native_introns.tsv` / `_summary.tsv`, `mitom_lo_introns.tsv` / `_summary.tsv`,
`products_native.tsv`, `make_set.py`, `check_products.py`, `run.sh`. The GenBank records and MitoM outputs (5 MB +
runs) are in `/media/cody/18T_R1/yzz/_realtest/intron_ncbi_2026-11-09/` on cody; `make_set.py --fetch` rebuilds them.
