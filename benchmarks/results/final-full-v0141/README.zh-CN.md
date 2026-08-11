# v0.14.1 完整数据发布门禁

[English](README.md)

状态：**36/36 个首轮候选结果全部通过完整文件 `cmp`**。

- 输入：11,230,392个位点、412个样本；每次固定保留5,425,725个位点。
- 负载：七个Original兼容过滤参数加`--recode --recode-INFO-all --stdout`。
- 场景：BGZF+TBI、BGZF+automatic CSI、Plain VCF、BCF自适应stream。
- 线程：`1 2 4 8 12 16 24 28 32`；全部严格串行运行。
- Original：复用v0.12.1锁定的VCFtools 0.1.17科学golden及单次时间；
  v0.14.1没有重跑Original。
- 自动CSI：bcftools 1.24 / HTSlib 1.24。
- 候选：本地可选PGO二进制SHA-256
  `0d585fe4e146feac1b63dd93a0cd238b2cd5b1baaa1d51a21a367c8f4b0b714e`。

首轮wall time（秒）：

| 输入 | Original | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF + TBI | 2267.88 | 203.72 | 146.86 | 75.48 | 44.36 | 33.47 | 34.06 | 32.43 | 30.93 | 32.84 |
| BGZF + automatic CSI | 2267.88 | 206.13 | 248.91 | 132.26 | 98.61 | 89.97 | 89.53 | 89.16 | 87.79 | 86.22 |
| Plain VCF | 2092.91 | 204.35 | 97.04 | 58.69 | 44.39 | 42.13 | 43.62 | 41.89 | 42.52 | 41.69 |
| BCF自适应stream | 1943.47 | 549.11 | 466.67 | 466.99 | 120.01 | 62.82 | 50.46 | 44.36 | 43.67 | 45.00 |

这些是首轮单次值，不是均值；automatic CSI包含建索引时间。严格预算下的
BCF在1–4线程受串行记录解码限制；本机28/32线程边际不能解释为运行时上限
或更高线程服务器的扩展承诺。

永久提交的精简证据包括：

- [`asset-validation.tsv`](asset-validation.tsv)：输入、索引和golden的
  大小/SHA-256门禁；
- [`manifest.tsv`](manifest.tsv)：候选、负载、环境和策略；
- [`all-runs.tsv`](all-runs.tsv)：每个Original/候选任务一行；
- [`summary.tsv`](summary.tsv)：wall、加速比、CPU、峰值RSS、backend和
  exact状态；
- [`run-v0141-full-release-matrix.sh`](../../run-v0141-full-release-matrix.sh)：
  可断点续跑的发布驱动。

大输入、锁定Original golden、每任务日志保留在基准主机本地。每个57–59 GB
候选输出只有在逐字节通过后才删除；最终scratch目录为空。
