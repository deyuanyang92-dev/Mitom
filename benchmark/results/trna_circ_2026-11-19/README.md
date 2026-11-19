# tRNA: circular origin in MitoM-ML, MiTFi only when needed (2026-11-19, cody)

结果数字只写在下表；中文和英文段落都指向下表。 / The numbers are in the table only.

| check / 检查 | data / 数据 | before / 改前 | after / 改后 |
|---|---|---|---|
| rotations whose MitoM-ML calls differ from the unrotated record / 旋转后 tRNA 调用与原记录不同的旋转数 (`rotation_check.py`) | NC_000931.1, NC_001636.1, MZ052209.1, NC_010217.1; each tRNA put at position 1, across the origin (+30) and its end at position 20 | 8 / 219 (all MZ052209.1) | 0 / 219 (`MITOM_TRNA_CIRC=1`) |
| unrotated records, calls (positions) old vs CIRC_FIX / 未旋转记录 | 8 circular records of tests/data | - | same positions in 8/8; p changed in 3 calls (NC_001636 trnD 0.944 -> 0.977, NC_006902 trnM 0.971 -> 0.964, NC_010217 trnK 0.708 -> 0.726) |
| scan with the last 300 nt prepended, hits / 前接 300 nt 扫描 | NC_000931.1 | 316 hits | 316 hits, same E, bits, anticodon (-Z kept) |
| linear record without the appended copy / 线性不加尾 | NC_000931.1 opened at 44 positions, linear | 0 fragment calls | 11 fragment calls of 28-62 nt (rejected, not kept) |
| MiTFi skipped / 跳过 MiTFi | 10 records of tests/data (`core_timing.tsv`, mito_core annotate, 2 threads) | MiTFi 1.8-18.5 s per record | skipped in 2/10 (NC_000931, NC_001636: -1.8 / -2.1 s of ~57 s); core.tsv identical 10/10 |
| golden (9 cases) | scripts/golden.py --jobs 2 | - | default and `MITOM_TRNA_CIRC=1`: all compared files identical; only "new output intron_fix.log" in 8 cases, which the base commit 9ba101d gives too |
| pytest | tests (no frozen / mito_sample) | - | 856 passed, 2 skipped |

中文：(1) 环状起点问题可复现：串联重复 trnP 的重叠特征和起点处 trnS1 的端点随旋转改变。修复（`MITOM_TRNA_CIRC=1`，默认关闭）：重叠按环计算、端点候选可跨起点、扫描时前接最后 300 nt。线性序列不加尾会产生片段调用，未采用。未训练类群的保留集（benchmark/trna_ml）数据只在 deyuan 机器上，cody 无法测准确率；只有起点附近的 tRNA 特征改变。
(2) MiTFi 只有在 p < 0.9 的保留调用、脊椎动物、或缺少的 tRNA 家族时才影响结果；其余情况跳过，结果在逻辑上相同（`MITOM_TRNA_ALWAYS_MITFI=1` 恢复旧行为）。节省很少，见上表。

English: (1) The origin problem reproduces: in MZ052209.1 the overlap feature of the tandem trnP copies and the ends of trnS1 at the origin change with the rotation. Fix (`MITOM_TRNA_CIRC=1`, off by default): circular overlaps, end pairs across the origin, the last 300 nt searched before the start. No appended copy for linear sequences gives fragment calls and was not kept. The held-out unseen-family set (benchmark/trna_ml) is only on deyuan, so accuracy was not measured on cody; only features of tRNAs near the origin change.
(2) MiTFi changes the result only through a kept call with p < 0.9, a vertebrate or a missing tRNA family; otherwise it is skipped and the result is the same by construction (`MITOM_TRNA_ALWAYS_MITFI=1` restores the old behaviour). The saving is small, see the table.

Scripts: `benchmark/trna_ml/scripts/rotation_check.py`, `linear_wrap_check.py`, `circ_fix_same.py`; `timing2.sh` + `mitfi_t.py` here.
