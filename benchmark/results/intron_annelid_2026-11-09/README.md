# Introns in annelid mitogenomes (2026-11-09)

问题：环节动物线粒体基因组的内含子，MitoM 能注释出多少？
Question: how many introns of annelid mitogenomes does MitoM annotate?

## Data / 数据

`scan_annelid.py`：NCBI 上全部 1,242 条 8–100 kb 的环节动物线粒体记录（2026-11-09），其中 18 条有内含子（intron
特征，或 CDS join 中 ≥ 142 nt 的间隔）。15 条已在金标准中（`../intron_gold_2026-10-29/`），另外 3 条是金标准中 NC_
记录的 INSDC 副本（OR123446 = NC_087808，OR123442 = NC_087806，KT726958 = NC_028710），所以金标准已包含 NCBI 上所有带
内含子的环节动物线粒体基因组。< 8 kb 的片段（多为 cox1 片段）没有测。逐条结果见 `annelid_records.tsv`。

All 1,242 NCBI annelid mitochondrial records of 8-100 kb: 18 have introns; 15 are in the gold standard and the other
3 are INSDC copies of its NC_ records. The gold standard therefore holds every annelid mitogenome with introns at
NCBI. Fragments < 8 kb were not tested.

## Results / 结果

Leave-one-out (`MITOM_EXCLUDE_REF`: the record out of the core's references and of the intron-site library), MitoM
601cc54, `--intron_rule continuity` (default); the runs of `../intron_algo_2026-11-09/` (`gold_cont3_introns.tsv`,
group Annelida).

| | 记录 / records | 内含子 / introns | 完全一致 / exact | ≤ 6 nt | 漏检 / missed | 多报 / extra |
|---|---|---|---|---|---|---|
| 有内含子 / with introns | 16 | 24 | 24 | 0 | 0 | 1 |
| 无内含子 / intron-free (`../intron_algo_2026-11-09/negatives_set.tsv`, RefSeq Annelida) | 37 | 0 | — | — | — | 0 |

- 24 个内含子：cox1 18 个、nad1 3 个、nad4 3 个，长 265–2,468 nt。Rfam 定为 II 类的 7 个，II 类但模型不完整的 7 个
  （`group II (?)`），未定型的 1 个（*Chaetopterus qiani* cox1，按接头连续 + 剪接必要写入），其余 9 个由位点库转移
  （`mito_intron_ref.py`）。
- 多报的 1 个是 *Glycera unicornis* KT989324.1 的 cox1 II 类内含子：NCBI 没有注释，但金标准报告中判断它很可能是真实
  的（`../intron_gold_2026-10-29/`）。
- 24 introns: cox1 18, nad1 3, nad4 3, 265-2,468 nt; 7 Rfam group II, 7 group II with an incomplete model, 1 untyped
  (*Chaetopterus qiani* cox1, written by junction continuity + splicing necessity), 9 placed by site transfer. The one
  extra is the unannotated, probably real group II cox1 intron of *Glycera unicornis* KT989324.1.

## NCBI 上没有任何注释的记录 / Records NCBI holds without any annotation

上面只统计了 NCBI 注释了内含子的记录。1,242 条中有 335 条**没有任何注释**（只有序列；多鳞虫科 50、龙介虫科 26、小头虫
科 17、才女虫科 14 ……），其中的内含子 NCBI 看不到，所以“只有 18 条有内含子”低估了实际情况（燐虫科、小头虫科尤其多）。
这 335 条用 `mitom run` 默认参数注释（遗传密码 5，`unannotated_introns.py` → `unannotated_introns.tsv`），没有真值，
所以按两条证据核对：写出的内含子剪接后蛋白是否完整；没有写出内含子的基因是否偏短或缺失，以及缺的部分是否在记录的别处
（`missed_intron_check.py`：用 refseq89m 的 72 个环节动物的同名蛋白做 tblastn → `missed_intron_check.tsv`）。

335 NCBI annelid records have no annotation at all, so their introns are invisible at NCBI. `mitom run` on them (no
truth): the introns written were checked for a complete spliced protein, genes without an intron for being short or
missing, and the missing part for lying elsewhere on the record (tblastn of the 72 annelid refseq89m proteins).

| | 记录 / records | 内含子 / introns |
|---|---|---|
| 写出内含子 / introns written | 12 | 20（cox1 13、nad2 2、nad4 2、nad6 2、cob 1） |
| 其中燐虫科 / Chaetopteridae | 4（*Mesochaetopterus* 2 种各 4 个：cox1、nad2、nad4、nad6；*Phyllochaetopterus* sp. cox1、cob；*Chaetopterus dewysee* cox1） | 11 |
| 其中小头虫科 / Capitellidae | 3（*Notomastus* sp.、*Dasybranchus caducus*、*Notomastus* sp. B：cox1） | 3 |
| 剪接后蛋白完整、无内部终止 / spliced protein complete, no internal stop | 12 | 20 / 20 |
| 基因偏短或缺失、缺的部分在别处 ≥ 142 nt（可能漏检）/ short or missing gene, rest >= 142 nt away (possibly missed) | 6 | — |

- 燐虫科、小头虫科的 23 条无注释记录：除 3 条 *Capitella capitata*（PV742346 只有 9.7 kb；PV742351 用 *C. teleta* 的
  COX1/COX2 做 tblastn 找不到任何同源，记录本身缺这两个基因）外，13 个基因都完整、无内部终止，没有漏检内含子的迹象。
- 6 处可能的漏检中，5 处是整段基因没有注释出来（龙介虫科 *Hydroides* / *Protula* atp6、*Trilobodrilus* /
  *Trochonerilla* nad5、*Thysanocardia* nad2），属于高度分化基因的识别问题；1 处是 *Helobdella europaea* OZ317477 的
  nad2（只注释出 159 aa，其余部分在 ≥ 142 nt 之外）：内含子步骤找到过一个 8.9 kb 的候选，剪接后有内部终止，未写出，
  在 `intron_report.tsv` 中列为待查。这条水蛭基因组还缺 atp6、nad5，cox1 / cob / nad1 都偏短，需要人工核查。
- 2 条燐虫的 cox2 在内含子步骤后被标为“仍有内部终止”，后面的边界校正（P7）已改正：提示文字已注明以最终 NCBI 检查为准。
- In the 23 unannotated Chaetopteridae / Capitellidae records all 13 genes are complete without internal stops, except
  in 3 *Capitella capitata* records that lack the genes (no COX1 / COX2 homolog at all): no sign of a missed intron.
  Of the 6 possible misses elsewhere, 5 are whole divergent genes the annotation did not find; 1 is *Helobdella
  europaea* nad2 (159 aa; an 8.9 kb candidate intron had an internal stop and is listed for checking).

局限 / Limits：注释过的记录只有 16 条、24 个内含子，以 cox1 II 类为主；无注释记录没有真值，只能按蛋白完整性核对。
Annotated: only 16 records / 24 introns; unannotated: no truth, checked by protein completeness only.

## Files / 文件

`scan_annelid.py`, `annelid_records.tsv`, `unannotated_introns.py` / `.tsv`, `missed_intron_check.py` / `.tsv`.
Work data on cody: `/media/cody/18T_R1/yzz/_realtest/orf_2026-11-09/annelid/` (1,242 records) and `annelid_unannot/`
(the 335 runs, `mito_bench_run.py --tools mitom_native`, MitoM 601cc54 + this commit's ORF step).
