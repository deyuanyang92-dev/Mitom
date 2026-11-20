# 外部审核报告核对 / Check of an external review report (2026-11-20, cody)

Report: `~/mitom (1).md` (not in the repo). Each claim was read against the code of branch exp/review-fixes-2026-11-19.

| 条目 item | 结果 result |
|---|---|
| P0-1 `placed[5]` -> `[6]` | 真；已修 / real, fixed |
| P0-2 `complete_claim` 只读第 1 条 | 真；170 条 RefSeq 中 4 条会误用第 1 条的 complete genome 声明；已改为逐记录（缓存）/ real, 4 of 170 records wrong, fixed |
| P0-6 线性 rRNA fill 到序列末端 | 真；只填到相邻基因 / real, fixed |
| P0-8 缩写终止子末密码子 | 真；170 条真实记录结果不变，合成记录可触发 / real; real data unchanged |
| P1-3 `.gbff` | 真；旧代码崩溃，新读 75 条 / real, fixed |
| P1-7 `genes.fna` 含 pseudo | 真；已修 / real, fixed |
| P2-1 `Path("")` | 真；已修 / real, fixed |
| 元数据 O(N²) / metadata O(N²) | 报告过时：9ba101d 已缓存 / stale, already fixed |
| P0-5 线性 tRNA | 部分过时：CIRC_FIX 已有 `circular` / partly stale |
| mitom_groups 线程风暴 | 说重：n_jobs 已按样本数分配 / overstated |
| P0-4 code 14/33 | 理论上对；本地无数据，NCBI 几乎不用 14/33，未改 / not changed, no data |
| P0-9 span 大内含子 | 2026-10-31 已按最大空隙修过，极端情形未复现 / not reproduced |

第二批 / second batch (same day): 
| P0-3 /pseudo、/transl_except | 真；pcg_check / origin_join_check 先用 apply_excepts，导入的 /pseudo CDS 不再 curate、不写 /translation，O1 保留内部例外，带例外的 CDS 跳过 P7/P8 / real, fixed |
| P1-1 KeepDeNovo | 部分：合并与单独转移模型都有内部终止子时保留 de novo；只有合并模型有终止子时仍用转移模型（de novo 多出的内含子多半是错的）/ partly, fixed |
| P1-2 pcg_hmm | `--k` 已传给 prior()（报告过时）；改了错误的帮助文字，补 /translation；`taken` 环形取模未验证 / partly stale |
| P2-4 .mitom.lock | 真；临时文件 + os.link 原子建锁；48 进程同时抢锁，旧代码 3/3 失败、新代码 3/3 通过 / real, fixed |

未处理 / not done: P2-5 fork -> spawn (risk, no evidence), P1-2 `taken` modulo n, P0-4 code 14/33 (no data).
