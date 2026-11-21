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
