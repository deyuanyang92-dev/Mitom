# ORFs of unknown function in animal mitogenomes (mito_orf.py, 2026-11-09)

问题：线粒体基因组里 13 个蛋白基因以外的编码区（如才女虫科报道的未知 ORF），MitoM 能否注释、会不会乱报？与 MITOS2 比较。
Question: can MitoM annotate protein-coding regions beyond the 13 genes (e.g. the unknown ORF of spionid polychaetes)
without writing chance ORFs, and how does it compare with MITOS2?

## Method / 方法（`mito_orf.py`; `mitom run --orf yes`, default）

ORF 是否是基因，取决于它是否受选择、一直保持开放；只看长度判断不了。富含 AT 的序列终止密码子少，80–110 个密码子的 ORF
随机就常见，而短的真蛋白又不比随机 ORF 长。所以按以下原理判定：
An ORF is a gene when selection has kept it open; length alone cannot tell (AT-rich sequence holds 80-110-codon
chance ORFs; short real proteins are no longer). Evidence, from first principles:

1. Length beyond chance: `E_len = 2 F p_start (1 - p_stop)^L` from the base composition of the record's unannotated
   sequence (F nt); `E_len <= 0.01`.
2. Kept open: a homolog (blastp E <= 1e-5) that is itself an ORF - in the other records of the run, the ORF library
   (`models/orf_library.fna.gz`: 120,248 ORFs in unannotated sequence and 371 non-standard CDS of 16,073 RefSeq animal
   mitogenomes, `scripts/build_orf_library.py`) or `--orf_ref` relatives. The d nucleotide differences between the
   codons of the two copies (blastp HSP) are a lower bound of the substitutions; each makes a stop with probability f
   (from the ORF's codons) in non-coding sequence: `p_open = (1 - f)^d <= 0.01`.
3. Coding frame: in a local nucleotide alignment of the two (which does not know the frame, unlike the protein
   alignment that favours differences at position 3), the differences fall mostly at codon position 3: binomial
   `p_pos3 <= 0.01` against 1/3.
4. Tandem repeats (copies >= 85 % identical, period 3-150 nt) are masked first: an array of a stop-free unit is an
   ORF of any length and aligns with any other array (control regions of birds, fishes, insects - the false hits of
   the first version).
5. ORFs overlapping a stronger one (other frame or strand) are its shadow and are dropped.

`conserved` (2 + 3) and `known_family` (homolog in the MFannot protein families; atp9 gets its own name) are written
as CDS; `long_unconfirmed` (1 only), `standard_fragment` (homolog of one of the 13 genes) and ORFs whose frame runs on
into a gene without a stop are WARNINGs in the log, rows in `orf_report.tsv` and text in the NCR /note - not silent.

## Results / 结果

| 数据 / set | 记录 / records | 写出 / written | 保守但未写出 / conserved, not written | 较长未确认 / long, unconfirmed | 13 基因片段 / fragment of the 13 |
|---|---|---|---|---|---|
| 才女虫科 / Spionidae (NCBI, >= 10 CDS) | 46 | 4 | 8 | 1 | 0 |
| 无内含子 RefSeq，26 门 / intron-free RefSeq, 26 phyla | 362 | 1 | 0 | 18 | 4 |
| 无注释环节动物（MitoM 注释）/ unannotated annelids (annotated by MitoM) | 335 | 16 | 9 | 29 | 21 |

- **才女虫科**：Polydorini（*Polydora*、*Boccardiella*）trnQ 与 trnA 之间的 ~355 个密码子 ORF，在 4 条有注释的记录中
  全部写出（NCBI 都没有注释）；在 335 条无注释记录中的 3 条 *Polydora* 中也写出（另 1 条 *P. caeca* 写出 orf137）。
  同一组属里另有两处很强的编码信号（ND5 前 ~140 个密码子、trnE 与 trnK 之间 ~190 个密码子），读框没有终止就读进 ND5
  或 trnK/16S，所以没有写出，只给 WARNING（ND5 那处与 ND5 同读框，可能是 ND5 的起点更靠上游）。
- **362 条无内含子记录**：只写出 1 个，*Aurelia limbata* orf340，与同属 *A. coerulea* 的 ORF 同源（钵水母线粒体基因组
  有额外的 ORF）。第一版写出的 7 个跨门“同源”全部是控制区串联重复，加入重复屏蔽和核苷酸比对检验后消失。
- **335 条无注释环节动物**：写出的 16 个是 3 条 *Polydora* 的 ~355 密码子 ORF、*P. caeca* orf137、*Spirobranchus
  corniculatus* 两个相邻、互为重复的 ORF（5 条同种记录 × 2）、*Hydroides minax* orf81 和 *S. giganteus* orf152。保守但
  未写出的 9 个：6 个是龙介虫科 *Hydroides* 与 COX1 / ATP6 同读框的 ORF（这些基因的起点可能更靠上游），3 个是
  *Polydora* ND5 前那一处，都给 WARNING。
- **MITOS2** 在三条 Polydorini 记录（NC_061377、OQ078742、MW528029）中都没有注释这个 ORF（只有 13 个基因、tRNA、
  rRNA 和 OH/OL）。MitoM 单条运行时，OQ078742 和 MW528029 写出 orf353；NC_061377 的近缘属不在 RefSeq，只给 WARNING，
  加 `--orf_ref`（OQ078742）后写出（`NC_061377_with_orf_ref.tsv`）。
- **Golden 回归**：*Plakinastrella* NC_010217（只给 FASTA）中注释核心漏掉的海绵 atp9 由本步骤补上（known_family）。

- Spionidae: the ~355-codon Polydorini ORF between trnQ and trnA is written in all 4 annotated records (NCBI
  annotates none) and in 3 unannotated *Polydora*; two more strong coding signals read on into ND5 (in frame) or
  trnK/16S without a stop and are WARNINGs.
- 362 intron-free records: 1 written (*Aurelia limbata* orf340, homolog in *A. coerulea*). The 7 cross-phylum "homologs"
  of the first version were tandem repeats of control regions; repeat masking and the nucleotide-alignment test
  removed them.
- MITOS2 annotates none of these ORFs. A single genome whose relatives are not in RefSeq gets a WARNING; with
  `--orf_ref` relatives it is written.

## Limits / 局限

- 证据来自同源比较：一个类群如果在 RefSeq、同一次运行和 `--orf_ref` 中都没有分化足够的近缘序列，真实 ORF 只能给出
  WARNING。短 ORF（< 60 个密码子）要求同源有足够多的替换，才能证明受选择。
- 读框没有终止就读进下一个基因的 ORF 不写出（终点无法确定）。
- Evidence needs a homolog diverged enough to show selection; short ORFs need many substitutions; ORFs whose frame
  runs into the next gene are not written (no end).

## Files / 文件

`run.sh`; `*_orf_report.tsv` (one per set); `NC_061377_with_orf_ref.tsv`. Work data on cody:
`/media/cody/18T_R1/yzz/_realtest/orf_2026-11-09/`.
