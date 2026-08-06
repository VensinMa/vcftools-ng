# vcftools-ng v0.13.0 W03–W10：230k 固定性能矩阵

日期：2026-08-06  
状态：v0.13.1 发布证据；本次发布采用固定 23k 正确性门禁和 230k 性能矩阵。  
范围：W03、W04、W05、W07、W08、W09、W10；W11/W12 不在本轮优化范围内。

## 门禁与复现边界

- 固定输入：`osmanthus412.flags.23chr_10k.vcf`
- 记录数：230,000；样本数：412；输入大小：3,201,584,833 bytes
- 输入 SHA-256：`e62f8f06617371229825bb4747777469ac173c56aa4bc03cec469b29bbee1be2`
- Original：VCFtools 0.1.17；只使用已经锁定的一次基线和永久 oracle，没有重新运行
- vcftools-ng：每个场景、每个线程运行 3 次，表格使用 wall time 中位数
- 线程：1、4、8、16、32
- 230k 门禁：225/225 PASS；每次结果均与 Original oracle 逐字节一致
- 23k 全矩阵门禁：285/285 PASS；包含 GATK 多等位 W10
- 实际参与性能测试的候选二进制 SHA-256：`fd23cd1b67ea85c5617e6adcd4f1b1b3b54b72829c6ef705da6af34702e64e17`
- 修改前冻结候选 SHA-256：`732b4d9c2866d4302dd4db7a9f4d2ee10ea34b08a6b835280d7ed429984da59f`

原始记录：

- `original-runs.tsv`
- `baseline-current/runs.tsv`
- `final-mmap/runs.tsv`
- `oracles/SHA256SUMS`
- 运行器：`../../run-v0130-w03-w10-230k.sh`

## 最终 wall time 中位数（秒）

| 场景 | Original | 1 | 4 | 8 | 16 | 32 |
|---|---:|---:|---:|---:|---:|---:|
| W03 positions 1% shuffled | 0.98 | 0.26 | 0.14 | 0.13 | 0.13 | 0.12 |
| W03 positions 1% sorted | 1.00 | 0.26 | 0.15 | 0.14 | 0.13 | 0.12 |
| W03 positions 50% shuffled | 2.53 | 1.37 | 0.40 | 0.25 | 0.22 | 0.16 |
| W03 positions 50% sorted | 2.53 | 1.37 | 0.41 | 0.26 | 0.22 | 0.16 |
| W04 exclude 1% duplicates/absent | 4.00 | 2.46 | 0.73 | 0.39 | 0.31 | 0.20 |
| W04 exclude 50% duplicates/absent | 2.55 | 1.38 | 0.40 | 0.28 | 0.21 | 0.16 |
| W05 keep 100% counts | 4.00 | 2.47 | 0.67 | 0.41 | 0.31 | 0.19 |
| W05 keep 25% counts | 3.91 | 1.18 | 0.35 | 0.23 | 0.18 | 0.14 |
| W05 keep 50% counts | 3.94 | 1.61 | 0.45 | 0.27 | 0.22 | 0.15 |
| W07 window pi non-overlap | 4.06 | 2.70 | 0.68 | 0.41 | 0.31 | 0.20 |
| W07 window pi overlap | 4.16 | 2.68 | 0.66 | 0.38 | 0.32 | 0.20 |
| W08 Tajima's D 100 kb | 4.03 | 2.69 | 0.67 | 0.39 | 0.32 | 0.20 |
| W09 site FST large pair | 4.55 | 2.58 | 0.70 | 0.39 | 0.33 | 0.22 |
| W09 site FST small pair | 4.09 | 0.87 | 0.28 | 0.19 | 0.15 | 0.14 |
| W10 window FST biallelic | 4.55 | 2.60 | 0.71 | 0.41 | 0.33 | 0.21 |

## 相对 Original VCFtools 0.1.17 的加速倍率

| 场景 | 1 | 4 | 8 | 16 | 32 |
|---|---:|---:|---:|---:|---:|
| W03 positions 1% shuffled | 3.77x | 7.00x | 7.54x | 7.54x | 8.17x |
| W03 positions 1% sorted | 3.85x | 6.67x | 7.14x | 7.69x | 8.33x |
| W03 positions 50% shuffled | 1.85x | 6.32x | 10.12x | 11.50x | 15.81x |
| W03 positions 50% sorted | 1.85x | 6.17x | 9.73x | 11.50x | 15.81x |
| W04 exclude 1% duplicates/absent | 1.63x | 5.48x | 10.26x | 12.90x | 20.00x |
| W04 exclude 50% duplicates/absent | 1.85x | 6.37x | 9.11x | 12.14x | 15.94x |
| W05 keep 100% counts | 1.62x | 5.97x | 9.76x | 12.90x | 21.05x |
| W05 keep 25% counts | 3.31x | 11.17x | 17.00x | 21.72x | 27.93x |
| W05 keep 50% counts | 2.45x | 8.76x | 14.59x | 17.91x | 26.27x |
| W07 window pi non-overlap | 1.50x | 5.97x | 9.90x | 13.10x | 20.30x |
| W07 window pi overlap | 1.55x | 6.30x | 10.95x | 13.00x | 20.80x |
| W08 Tajima's D 100 kb | 1.50x | 6.01x | 10.33x | 12.59x | 20.15x |
| W09 site FST large pair | 1.76x | 6.50x | 11.67x | 13.79x | 20.68x |
| W09 site FST small pair | 4.70x | 14.61x | 21.53x | 27.27x | 29.21x |
| W10 window FST biallelic | 1.75x | 6.41x | 11.10x | 13.79x | 21.67x |

## 相对本轮修改前冻结候选的 wall time 变化

负值代表更快。

| 场景 | 1 | 4 | 8 | 16 | 32 |
|---|---:|---:|---:|---:|---:|
| W03 positions 1% shuffled | -46.9% | -56.2% | -56.7% | -58.1% | -60.0% |
| W03 positions 1% sorted | -48.0% | -53.1% | -50.0% | -58.1% | -60.0% |
| W03 positions 50% shuffled | -14.9% | -25.9% | -30.6% | -40.5% | -50.0% |
| W03 positions 50% sorted | -14.4% | -24.1% | -25.7% | -40.5% | -50.0% |
| W04 exclude 1% duplicates/absent | -8.9% | -6.4% | -17.0% | -31.1% | -37.5% |
| W04 exclude 50% duplicates/absent | -14.3% | -25.9% | -22.2% | -43.2% | -48.4% |
| W05 keep 100% counts | -9.2% | -14.1% | -8.9% | -29.5% | -40.6% |
| W05 keep 25% counts | -16.9% | -28.6% | -30.3% | -48.6% | -53.3% |
| W05 keep 50% counts | -13.4% | -22.4% | -25.0% | -42.1% | -50.0% |
| W07 window pi non-overlap | -0.7% | -16.0% | -21.2% | -31.1% | -37.5% |
| W07 window pi overlap | -1.1% | -17.5% | -19.1% | -30.4% | -39.4% |
| W08 Tajima's D 100 kb | -2.2% | -16.2% | -15.2% | -25.6% | -37.5% |
| W09 site FST large pair | -8.8% | -14.6% | -17.0% | -25.0% | -31.2% |
| W09 site FST small pair | -20.2% | -33.3% | -44.1% | -55.9% | -54.8% |
| W10 window FST biallelic | -8.8% | -13.4% | -22.6% | -28.3% | -36.4% |

## 实现结论

本轮保留的主要改动：

1. positions/exclude-positions 直接文本路径及按染色体 cursor。
2. 样本投影只解析选中样本 GT，100% keep 退化为 identity projection。
3. π、Tajima's D、site/window FST 使用紧凑逐位点贡献和固定顺序归约。
4. 常见二倍体 GT（`0/0`、`0/1`、`1/1`、`./.`）以及 `GT` 位于 FORMAT 首列时使用专用解析路径。
5. 两群体二等位 FST 使用栈上定长累积器。
6. Plain VCF 自适应零拷贝：多线程使用只读 `mmap`；单线程仅在 positions 或样本投影等 I/O 主导任务中启用，纯 π/Tajima 路径继续使用 `pread`。

只读映射会让 Linux 的峰值 RSS 统计接近被访问的输入文件大小；这些主要是可回收的文件映射页，不是等量匿名堆内存。230k 输入下映射场景约 3,056–3,064 MiB RSS，单线程纯 π/Tajima `pread` 路径约 133 MiB。后续 230 万及完整数据测试必须同时观察 wall time、major faults、存储类型和系统内存压力，不能只按 RSS 数字否决映射路径。

未保留的实验：两群体多等位 FST 的 thread-local 扁平 scratch。23k W10 多等位在 1/4 线程回退 1.65%/4.80%，高线程变化与未改动的二等位场景噪声同量级，因此已撤销。

## 后续更大规模门禁

只有在需要确认阶段性收益时再运行固定 230 万真实子集。该门禁应同时验证：

- 23 条染色体每条 100,000 位点；
- 1/4/8/16/32 线程；
- Original oracle 逐字节一致性；
- SSD Plain VCF 的 `mmap` 收益能否随数据规模保持；
- 可回收映射页、major page faults 与系统内存压力；
- 若用于机械硬盘，必须单独评估并发 page-fault 是否造成寻道回退。
