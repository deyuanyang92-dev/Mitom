# 强内含子模型 vs 弱冲突 CDS / A strong intron model against a weak conflicting CDS (2026-11-20, cody; issue #15)

问题 / Problem: MW592986.1 atp6 — a continuous model (exonerate score 528, 98 % of the homolog) was refused ("model overlaps
existing ATP8 CDS") because the core had called a weak ATP8 (score 84, 53/55 aa) inside atp6 exon 2; the real atp8 is
at 3757..3855.

规则 / Rule (`mito_intron.outweighed_conflicts`): the conflict is waived, and `apply_model` replaces the weak CDS, only when
the model is already "high" (typed introns, no internal stop) and every conflicting CDS (a) has a core score in its
note, (b) has score x 3 <= the model's exonerate score, (c) lies >= 80 % inside the model's exons. A CDS without a core
score (the submitter's annotation) is never outweighed.

测试 / Test: the intron step (`intron_only.py`, leave-one-out) of the 134 NCBI records of `intron_ncbi_2026-11-09`,
old vs new code, same input.

| | 记录 / records | 报告不同 / reports differing |
|---|---|---|
| old vs new | 134 | 1 (MW592986.1) |

MW592986.1: atp6 16489..17094, low / not written -> high / written; spliced CDS = NCBI's intron, ATP8 call replaced.
The other 133 reports are identical.

Not done: the 16 bryozoan introns next to short divergent exons (needs a new search, see issue #15).
