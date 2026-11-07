# 冻结评测方案 / Frozen evaluation plan

论文中的准确度数字只取自这一次评测：版本冻结、测试集与 MitoM 的训练和调参数据无重叠、指标事先确定、只跑一次。
The paper's accuracy numbers come from this one run only: frozen version, a test set disjoint from everything MitoM
was trained or tuned on, metrics fixed in advance, run once.

## 1. 规则 / Rules

1. 运行前先提交本方案、测试集和冻结清单（manifest）；提交时间就是预注册时间。/ Commit this plan, the set and the manifest before running.
2. 看到任何结果后，不再修改代码、参数、测试集或指标。发现的问题留到下一个版本，下一轮评测用新测试集。/ No change after the first result; fixes go to the next version and a new set.
3. 失败、超时都计入（算作漏注）；所有结果都报告，包括 MitoM 落后的指标。/ Failures count; every result is reported, losses included.
4. 结果登记到 `knowledge/TEST_REGISTRY.md`，存放于 `benchmark/results/frozen_<日期>/`。/ Register the run; store it there.

## 2. 测试集 / Test sets

| 集 / set | 内容 / content | 何时 / when |
|---|---|---|
| A 独立属 / unseen genus（主 / primary） | RefSeq 动物线粒体全基因组；记录号从未出现在仓库任何训练、调参或评测文件中，**所在属也未出现过**；每属 1 条、每科 ≤2 条、每门 ≤40 条 / accession and genus never seen; ≤1 per genus, ≤2 per family, ≤40 per phylum | 现在 / now |
| B 时间外推 / temporal | 冻结日期之后发表的记录（`mitom fetch --query '... AND ("<冻结日期>"[PDAT] : "3000"[PDAT])'`），同样做泄漏审计 / records published after the freeze date | 冻结后 ≥ 3 个月 / ≥ 3 months later |

“见过”由 `frozen_set.py` 从仓库自动收集：`benchmark/`（所有集合、tRNA 训练拆分、历次结果）、`models/`（refseq89m 记录表、rRNA HMM 参考、内含子位点）、
`knowledge/experiments/`、`tests/data/` 中出现过的任何记录号；以及这些记录的属和科。tRNA 分支模型（clade CM）覆盖的类群单独标出（`clade_cm_rank`，1 = 同属）。
不在 refseq89m 中，也就同时保证 MITOS2 和 MitoM 的蛋白参考库里没有被测记录自身。

"Seen" is collected from the repository (any accession in benchmark/, models/, knowledge/experiments/, tests/data/, plus
their genera and families); tRNA clade-model coverage is reported separately. Not in refseq89m means neither MITOS2's
nor MitoM's protein references hold the record itself.

## 3. 冻结 / Freeze

`freeze.sh` 记录：git commit 与 tag（要求无未提交修改）、测试集 md5、`models/` 与 refseq89m 每个文件的 md5、
`mitom doctor`、各 conda 环境的完整导出（mitotools、mitos、mitoz、mitofinder）、各程序版本。
`run_frozen.sh` 在 HEAD 不是该 commit、有未提交修改或测试集被改动时拒绝运行。

## 4. 工具与命令 / Tools and commands

全部经 `benchmark/mito_bench_run.py` 运行（完整命令写入 `commands.log` 和每条记录的 `log.txt`）。输入为去掉注释的序列。
Run through mito_bench_run.py; every command is logged.

| 工具 / tool | 运行方式 / how | 给的信息 / given |
|---|---|---|
| `mitom_native` | `mitom run`（默认参数） | NCBI 记录的类群和遗传密码 / clade and code of the record |
| `mitom_blind` | `mitom run -i in.fa`，不给 `--clade`、`--genetic_code` | 只有序列 / sequence only |
| `mitos2` | `runmitos`（2.1.10，refseq89m，默认参数，full 模式） | 遗传密码 / code |
| `mitos2_cur` | `mitos2` 的结果经 `mitom import`（只做整理） | 同上 / same |
| `mitoz`、`mitoz_cur` | MitoZ 3.6 `annotate`；及其经 `mitom import` | 类群、密码；MitoZ 不支持的门记为 unsupported |
| `mitofinder`（可选） | MitoFinder 1.4.2 `-a`，参考 = 池中最近缘记录（排除同种） | 参考基因组 / a reference |

主比较 / primary comparison：`mitom_native` vs `mitos2`（两者信息相同）。

## 5. 指标（事先确定）/ Metrics (fixed in advance)

真值 = NCBI RefSeq 注释；匹配规则见 `mito_bench_eval.py`（重叠 ≥ 50 %，同名同链）。
| 指标 / metric | 定义 / definition |
|---|---|
| **pcg_exact（主 / primary）** | 蛋白基因两端与 NCBI 完全一致的比例 / both ends identical |
| pcg_le3、pcg_found | 两端差 ≤ 3 nt；找到 / within 3 nt; found |
| rrna_le3 | rRNA 两端差 ≤ 3 nt |
| trna_recall、trna_precision | tRNA 召回 / 精确率 |
| err_free | 无明显错误的记录比例（缺基因、内部终止、链错、重复、名称错、端点差 > 30 nt、table2asn ERROR、失败） |

区间：按科重抽样的 bootstrap（2,000 次，95 %），配对差值（同一组重抽样）；`frozen_stats.py`。
另报告分层：门、科是否见过、clade_cm_rank、NCBI 注释来源（`annotation_method`：真值由 MITOS2 产生的记录对 MITOS2 有利）、拓扑。
Intervals: family-clustered bootstrap, paired differences. Strata as listed.

## 6. 步骤 / Steps

```bash
conda activate mitotools
# 1. 候选池 / pool (all RefSeq metazoan mitogenomes)
mitom fetch --taxon Metazoa --refseq_only -o pool                       # -> pool/Metazoa.gb
# 2. 泄漏审计与选集 / leak audit and selection
python3 benchmark/frozen/frozen_set.py audit  --pool pool/Metazoa.gb -o frozen/pool_audit.tsv
python3 benchmark/frozen/frozen_set.py select --pool pool/Metazoa.gb -o frozen/test_accs.txt      # + test_accs.audit.tsv
python3 benchmark/make_large_set.py --pool pool/Metazoa.gb --out data --tsv frozen/set.tsv --accs frozen/test_accs.txt
# 3. 冻结 / freeze  (commit + tag first: git tag v0.1.0-frozen)
bash benchmark/frozen/freeze.sh frozen/set.tsv benchmark/results/frozen_<date>
git add benchmark/results/frozen_<date> && git commit -m "Frozen evaluation: set and manifest"   # = 预注册 / pre-registration
# 4. 一次运行 / the one run (resumable: rerun the same command after an interruption)
bash benchmark/frozen/run_frozen.sh benchmark/results/frozen_<date> data 8 2
# 5. 提交结果，登记 / commit the results, add the registry row
```

预计 / estimate：400 条，MITOS2 约 6 min/条，8 并行约 5 h；MitoM 约 1 min/条。建议在另一台电脑（deyuan）上运行。

## 7. 局限（论文中写明）/ Limits (state them in the paper)

- 真值是 NCBI 注释，本身来自不同软件、规则不一；“一致”≠“正确”。MitoM 的端点规则学自 NCBI 记录，所以更像 NCBI 写法。
  / NCBI annotation is the truth; agreement is not correctness; MitoM's rules were learned from NCBI style.
- 需要一个人工核实的子集（如每门 2–3 条，用转录组 / 文献核对端点）才能谈生物学准确度；现有内含子金标准已参与训练，不能用。
  / A hand-verified subset is needed for biological accuracy; the intron gold set was used in training.
- 默认参数比较；不对其他工具调参。/ Default parameters for every tool.
