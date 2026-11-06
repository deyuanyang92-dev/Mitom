# 安装与自检（2026-11-06）

## 发现的问题

| 问题 | 原因 | 修复 |
|---|---|---|
| 本机 `mitom --version` 显示 0.9.0，`mitom debug` 不存在 | `~/.local/bin/mitom` 是 10-23 旧安装留下的符号链接，指向 `~/mybiosofts/Mitom`（PR #2 时的副本）；之后没有重跑 `install.sh` | `mitom doctor` 检查 PATH 上的 `mitom` 运行哪个文件夹，不是本文件夹或是旧式符号链接时标 `!!`、给出修复命令并返回非零 |
| 旧命令用 `mitgard` 环境的 python3 运行，没有 scikit-learn，默认的 MitoM-ML tRNA 失败；旧 doctor 仍说 "All required programs found" | 启动器只检查 Biopython | 启动器要求 Biopython + numpy + scikit-learn，优先 `mitotools` 环境；doctor 检查 scikit-learn/numpy（缺失为 `!!`），并显示与验证版本 1.7.2 的差异 |
| `mitotools` 环境是 scikit-learn 1.6.1，`env/mitotools.yml` 要求 1.7.2 | 环境创建后 yml 改过，未再 `install.sh` | 重跑 `install.sh`（env update） |
| README 没有更新步骤 | — | 加 `git pull && bash install.sh`；求助时附 `mitom doctor` 和 `mitom debug out` 的输出 |
| Docker 全新安装测试无法进行 | 本机 Docker 守护进程未运行（需要 sudo） | 改用空 HOME、只含系统 PATH 的全新用户模拟安装 |

## 全新用户安装测试（空 HOME，PATH 只有 /usr/bin:/bin，无 conda，2026-11-06）

| 检查 | 结果 |
|---|---|
| `bash install.sh` | exit 0：自动装 Miniforge、建 mitotools（scikit-learn 1.7.2）、下载 refseq89m、写配置 |
| `mitom doctor` | exit 0，必需项全部 OK |
| `mitom run -i NC_000933.gb -o out` | exit 0，1 分 40 秒；`mitom debug out` 正确读出各步骤 |
| 发现的手动步骤 | 新系统上 `~/.local/bin` 不在 PATH，需要手改 `~/.bashrc` → `install.sh` 现在自动写入（只写一次；zsh/macOS 也写 `~/.zshrc`） |
| 测试中发现的新错误 | 新的 doctor 检查在安装过程中（PATH 尚未生效）报 "no mitom on PATH"，使 `install.sh` exit 1 → `~/.local/bin/mitom` 指向本文件夹时不算问题；修复后重装 exit 0 |
