# Introns of no known type: decided by biology, not by thresholds fitted to data (2026-11-09)

问题：Rfam 认不出 I / II 类的内含子，什么证据下才写入？旧规则（I 类剪接分数 ≥ 14、模型比上同源蛋白 ≥ 80 %）是从
基准数据挑出来的阈值，在新数据（苔藓虫，`../intron_ncbi_2026-11-09/`）上漏掉大半。
Question: when is an intron Rfam cannot type written? The old rule (group I splice score >= 14, >= 80 % homolog
coverage) was fitted to benchmark data and missed most introns of new data (bryozoans).

## Biology / 生物学依据

- 动物线粒体里认不出类型的内含子几乎都是退化的 II 类（环节动物 cox1；苔藓虫）。Jenkins et al. 2022（*Sci Rep*
  12:10889，PMC9240083）：Cupuladriidae 的 35 个内含子 221–519 nt，没有 II 类 V/VI 结构域，起点没有 GUGYG、16 个
  终点没有 AY，作者只能按相邻外显子的同源比对定边界；22 个经 PCR 证实存在。I 类剪接规则对它们不适用。
- An intron removes DNA only: the gene's homologs run on across the junction of the spliced protein (**continuity**).
- An intron interrupts the coding sequence: without splicing the gene does not cover its homolog in place
  (**necessity**).
- Both hold for any intron type and any lineage; the old thresholds were stand-ins for them.

## Rule / 规则（mito_intron.py `--untyped_rule continuity`, `mitom run --intron_rule continuity`, default）

An untyped intron is written when
1. `junction_continuity`: more than half of the homologs (the model's, the core's reference, up to nine of the
   library) align on both sides of the junction, and the median of the homolog residues missing / added there is
   <= 1 (the codon the junction splits). A discontinuous junction is first searched +-30 nt at each end for the
   placement that makes it continuous;
2. `splicing_gain` >= 20 BLOSUM62 points (~9 bits; random sequence scores below 0): the spliced protein beats both
   unspliced alternatives (the exon before dropped and the exon after read on upstream in frame to the first stop, and
   the mirror case).

Typed introns are not subject to these (their structure is the evidence). The coverage gate is dropped under this rule
(a divergent homolog aligns in part; a false junction shows as discontinuity). `--intron_rule splice` gives the old
rule. Also new: a gene whose annotated CDS has an internal stop is an intron candidate (imported annotations, other
engines); a CDS that still has one after the intron step is listed in `intron_report.tsv` for checking. The report
gains the columns `junction_span`, `junction_missing_aa`, `junction_extra_aa`, `splice_score`, `splicing_gain`.

## Data / 数据

| 数据集 / set | 记录 / records | 真值 / truth |
|---|---|---|
| 金标准 / gold standard (`build_gold.py`, `gold_set.tsv`) | 270：193 条带内含子 + 76 条无内含子 + KT989324.1 | 309 个逐个核实的内含子 / verified introns |
| 金标准以外的 NCBI 记录 / NCBI records outside it (`../intron_ncbi_2026-11-09/`) | 134 | NCBI 内含子特征 / NCBI intron features |
| 无内含子 / intron-free (`build_negatives.py`, `negatives_set.tsv`) | 362 RefSeq，按门分层，26 个门 / 362 RefSeq, 26 phyla | 无内含子 / none |

Every record: `mitom run` leave-one-out (`MITOM_EXCLUDE_REF`) with `--intron_fix no --keep_intermediate yes`
(`mito_bench_run.py --tools mitom_prein`), then only the intron step (`intron_only.py`: mito_intron.py +
mito_intron_ref.py, as `mitom run`), scored with `benchmark/eval_introns.py`. The intron-only rerun with the old rule
gives exactly the full leave-one-out runs (gold 299 / 4 / 0 / 6, 2 extra; NCBI set as `../intron_ncbi_2026-11-09/`).

## Candidates / 候选（`candidates_*.tsv`: every model written, MITOM_INTRON_WRITE_ALL=1, de novo step only）

Typed candidates: 469 true, 2 "false" that are introns NCBI does not annotate (KT989324.1 cox1, group II, see the
gold-standard report; BK068680.1 cox1, group I, NCBI's cox1 reads through it). Untyped candidates:

| 证据 / evidence | 真 / true | 假 / false |
|---|---|---|
| 接头不连续 / junction not continuous | 2（MW592987.1 nad2、nad5，位置偏 > 6 nt） | 10（海绵 nad2/nad5、丝盘虫 nad6 ×2、刺胞 nad1/nad6、软体 cox2 ×3、轮虫 nad1） |
| 连续但不必要 / continuous, not needed (gain < 20) | 0 | 1（*Kuphus* NC_063762.1 nad3：2.5 kb "内含子"包住 nad4，原位 nad3 完整） |
| 连续且必要 / continuous and needed | 25（gain ≥ 96） | 2 → 实为 NCBI 漏注的 cox3 内含子（gain 34、36，下） |

The two "false" untyped candidates that pass both are a third cox3 intron in both cupuladriids (MW592987.1
18194..18458, MW592989.1 12024..12282): NCBI's cox3 stops inside the intron (`...WHFVDVVWLFL SII*`); spliced, it ends
`...WHFVDVVWLFL FVSIYWWGK*`, the conserved COX3 C-terminus (human `...WHFVDVVWLFLYVSIYWWGS`). Their gain is small
because the last exon is 9 codons long.

## Results / 结果（`*_base_*` = old rule, `*_cont3_*` = new rule）

| 数据集 / set | 规则 / rule | 完全一致 / exact | ≤ 6 nt | > 6 nt | 漏检 / missed | 多报 / extra |
|---|---|---|---|---|---|---|
| 金标准 309 / gold 309 | 旧 / old | 299 | 4 | 0 | 6 | 2 |
| | 新 / new | 299 | 4 | 0 | 6 | 2 |
| 苔藓虫 50（NCBI） / bryozoans 50 | 旧 / old | 11 | 0 | 1 | 38 | 0 |
| | 新 / new | 23 | 7 | 1 | 19 | 2 |
| 刺胞 189 + 海绵 6（NCBI） / cnidarians + sponges | 旧 = 新 / old = new | 124 | 40 | 30 | 1 | 1 |
| 无内含子 362 / intron-free 362 | 旧 / old | — | — | — | — | 0 |
| | 新 / new | — | — | — | — | 0 |

The gold extras are the two probably real introns of `../intron_gold_2026-10-29/` (KT989324.1 cox1, Trichoplax rrnL);
the bryozoan extras are the cox3 introns above. Golden regression: all GenBank outputs identical; intron_report.tsv
gains the five columns, run_params.json the `intron_rule` key.

## Still missed / 仍未解决

- 16 bryozoan introns next to short, divergent exons (28-81 aa; nad4, nad4L, nad6, cox2 / cytb end exons): exonerate
  does not align these exons at all (each report must reach a raw score of 200 to be chained), so no model has the
  intron. Needs a homology-constrained exon search in the known gap with a significance test, not a lower threshold.
- MW592986.1 atp6: a strong, continuous model is refused because the core annotated a weak ATP8 (score 84) inside
  atp6 exon 2 (the real atp8 is at 3757..3855). Needs a general rule for a model against a weak conflicting CDS.
- MW592987.1 nad2 / nad5: junctions 6-7 homolog residues short with no continuous placement within +-30 nt (perhaps a
  real deletion next to the junction); left unwritten.

## Files / 文件

`build_gold.py`, `build_negatives.py`, `intron_only.py`, `label_candidates.py`; `*_set.tsv`, `candidates_*.tsv`,
`*_summary.tsv`, `*_introns.tsv`. Work data on cody: `/media/cody/18T_R1/yzz/_realtest/intron_algo_2026-11-09/` and
`../intron_ncbi_2026-11-09/` root.
