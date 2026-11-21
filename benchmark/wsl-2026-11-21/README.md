# WSL2 性能测试记录 (2026-11-21)

对应 issue: #16(批量慢)、#17(内存)、#18(安装大)、#19(下载慢)。本目录给出原始日志与复现方法,便于据此做优化。

## 环境
- MitoM 0.0.1, git dfafa9b;Python 3.11.18(conda env mitotools,mamba 安装)
- WSL2 (Linux 6.18, glibc 2.35), 8 CPU, 可用内存约 3GB(.wslconfig memory=3GB)
- 输入在 ext4;对照一次放在 /mnt/d(9p)

## 数据
- `data/ann10.gb`:10 条环节动物(Annelida)线粒体基因组 GenBank 记录,368KB
- `data/accessions.txt`:对应的 10 个 accession(可从 NCBI 重新取得)

## 复现命令(`$WORK` 为任意工作目录)
```bash
/usr/bin/time -v -o t.txt mitom run -i data/ann10.gb -o $WORK/out                 # 基线,默认 jobs
/usr/bin/time -v -o t.txt mitom run -i data/ann10.gb -o $WORK/out --jobs 4 --threads 2
/usr/bin/time -v -o t.txt mitom run -i data/ann10.gb -o $WORK/out --jobs 8 --threads 1
```

## 结果(10 条)
| 配置 | 总耗时 | 单进程峰值 RSS | 说明 |
|---|---|---|---|
| 基线(ext4) | 11m44s | 317MB | time -v 显示 CPU 利用 384%,默认已多核 |
| 基线(/mnt/d 9p) | 12m15s | 317MB | 9p 不是主要瓶颈 |
| jobs4 threads2 | 10m10s | 317MB | |
| jobs8 threads1 | 8m11s | 317MB | 仅 1.43x |

注意:`time -v` 只报告主进程峰值,不含 snakemake 子进程的总内存,`--jobs 8` 的总内存尚未采样(待测)。

## 基线步骤耗时(总约 683s,见 logs/baseline_jobs1/run.log)
| 步骤 | 耗时 | 占比 |
|---|---|---|
| STEP 2 annotate(核心 snakemake) | 604.7s | 88%,约 60s/条 |
| STEP 6 add_missing | 约 59s | 其中 intron 54.5s(mito_intron 33s、intron_ref 21s);pcg_hmm 4.7s |
| orf_ncr | 11.8s | |
| curate | 2.8s | |
| ncbi | 2.6s | |

## 目录
- `logs/<配置>/run.log` 主流程日志(含各步耗时);`intron_fix.log`;`run_params.json` 实际参数;`stdout.log`;`time_v.txt` 完整 `time -v`
- `ext4_vs_9p_result.md`、`jobs_result.md` 汇总

## 建议的优化方向与待补数据
1. 优先看 annotate 核心步骤(88%):还缺单条 `--debug` 重跑,拆分 60s/条里哪个工具最慢。
2. add_missing 中 intron 约 54s/10 条,可考虑缓存参考库或并行。
3. `--jobs` 加速仅 1.43x,可能受 3GB 内存或核心步骤自身多线程影响,待补总内存采样。
4. 本机无法测试 16GB 内存场景。

## 单条 --debug 拆分(logs/single_record_debug/,data/ann1.gb = ann10.gb 第一条)
命令:`mitom run -i data/ann1.gb -o $WORK --debug --keep_intermediate yes`。总 1m41s,CPU 282%,峰值 RSS 317MB。`run.log` 中每条外部命令都有 `[CMD rc=0 N.Ns]`。

| 步骤 | 耗时 |
|---|---|
| STEP 2 annotate | 57.8s |
| STEP 6 add_missing | 26.2s(introns 21.3s,pcg_hmm 4.8s) |
| 其余步骤合计 | 约 3s |
| 启动与环境检查 | 约 10s(未计入步骤,从日志时间戳估算) |

annotate 内部各工具累计(同一条记录):
| 工具 | 调用 | 累计耗时 |
|---|---|---|
| blastx | 13 次(逐基因,串行) | 26.3s |
| cmsearch | 3 次 | 14.4s |
| cmalign | 22 次 | 10.2s(每次 0.1 到 1.0s) |
| nhmmer | 2 次 | 2.0s |
| blastp / tblastn / aragorn | 各 1 次 | 1.9s / 1.5s / 1.9s |
| exonerate | 47 次 | 0.6s |

结论:这三个工具合计约 51s,占 annotate 的 88%。
- blastx 对 13 个蛋白基因逐个调用,每次 0.7 到 4.3s,耗时在串行调用;可考虑合并为一次查询或并行。
- cmsearch 3 次共 14.4s,是单次最重的工具;可考虑限制搜索范围或复用结果。
- cmalign 22 次小调用,合并后可省启动开销。
- 注:`[CMD]` 里 python 包装调用的 85s 是嵌套外层,不要与上表相加。

## --jobs 8 --threads 1 总内存(logs/jobs8_total_memory/)
每 2s 对整棵进程树(同一 session)求和 RSS,`samples.tsv` 列:t_s、total_rss_MB、n_proc、mem_avail_MB。10 条 ann10.gb,总耗时 9m36s(含采样)。

| 指标 | 数值 |
|---|---|
| 总 RSS 峰值 | 1434MB(t=510s,19 个进程) |
| 系统可用内存最低 | 1008MB(WSL 约 2.9GB) |
| `time -v` 报告的峰值 | 317MB(只含主进程,低估) |

说明:总 RSS 把共享页重复计入,真实占用略低于 1434MB。结论:`--jobs 8` 在 3GB 下没有 OOM,但总内存约为 `time -v` 数字的 4.5 倍;评估 16GB/4GB 机器时应以总 RSS 为准,不能用 `time -v`。
