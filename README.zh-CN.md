# vcftools-ng

[English](README.md) | [简体中文](README.zh-CN.md)

vcftools-ng 是 VCFtools 0.1.17 的实验性高性能、输出兼容后继实现。

**最新正式版：**
[v0.12.1 — 融合位点统计与可扩展精确重编码](https://github.com/VensinMa/vcftools-ng/releases/tag/v0.12.1)

完整 11,230,392 位点的七场景第一轮发布门禁已经通过：
1/2/4/8/16/32 线程的 42/42 个候选输出均与 VCFtools 0.1.17 逐字节
一致且更快。第 2–5 轮暂缓，下面的全量结果均明确为单次结果。

vcftools-ng 使用自适应、有序的输入分片和有界流水线：

- Plain VCF 使用按完整记录边界对齐的并行字节范围；
- BGZF VCF + TBI/CSI 和 BCF + CSI 使用独立有序索引区域；
- 缺少索引的本地 BGZF VCF/BCF 默认使用有效线程预算调用 bcftools
  自动构建 CSI；
- 输入、计算和有序输出可以重叠，结果顺序保持确定；
- `--freq`、`--freq2`、`--counts`、`--missing-site`、
  `--site-depth`、`--site-mean-depth` 和 `--site-quality` 可以进入
  直接融合文本路径，避免构建中间 `bcf1_t`；
- 乱序 shard 的待提交输出带有背压，内存不会随剩余文件无限增长。

## 推荐安装方式

Linux x86_64 便携包解压后即可运行，包内包含 vcftools-ng、用于自动
构建 CSI 的 bcftools、HTSlib 以及所需的非 glibc 运行库：

```bash
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.12.1/vcftools-ng-v0.12.1-linux-x86_64.tar.gz
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.12.1/vcftools-ng-v0.12.1-linux-x86_64.tar.gz.sha256
sha256sum -c vcftools-ng-v0.12.1-linux-x86_64.tar.gz.sha256
tar -xzf vcftools-ng-v0.12.1-linux-x86_64.tar.gz
./vcftools-ng-v0.12.1-linux-x86_64/bin/vcftools-ng --version
```

便携包面向 glibc 2.17 或更高版本，在 CentOS 7 兼容的
manylinux2014 基线上构建。无需安装 CMake、编译器、Conda、系统
HTSlib 或系统 bcftools。解压后请保持 `bin`、`lib` 和 `libexec`
目录的相对位置不变。

## 兼容性原则

VCFtools 0.1.17 是兼容性 oracle。只有通过完整文件 `cmp` 的工作负载
才称为“精确兼容”；数值接近或行数相同不能替代逐字节验证。

“精确兼容”和“性能更快”是两个独立结论：

- 精确兼容表示记录的真实工作负载与 Original golden 完全一致；
- 性能提升表示该参数所在的真实工作负载在基准主机上实测更快；
- 这不等于承诺所有硬件、文件系统、异常输入和任意参数组合都更快。

详细参数状态、证据等级、输入覆盖和兼容边界见
[参数兼容性与优化状态](docs/parameter-compatibility.md)。

## 已兼容并有真实性能证据的 Original 参数

主要参数组包括：

- 位点统计：`--freq`、`--freq2`、`--counts`、`--missing-site`、
  `--site-depth`、`--site-mean-depth`、`--site-quality`；
- 个体/HWE：`--depth`、`--missing-indv`、`--het`、`--hardy`；
- 群体统计：`--site-pi`、`--window-pi`、`--TajimaD`、
  `--weir-fst-pop` 及其窗口参数；
- LD/PCA：`--geno-r2`、LD 窗口参数、`--pca`、`--pca-no-norm`；
- 转换：`--recode`、`--recode-bcf`、`--recode-INFO-all`、
  `--stdout`；
- 常用位点、深度、缺失、MAF/MAC/HWE、染色体、位置、BED、
  thinning、非参考等位基因、FILTER/INFO Flag、基因型和样本过滤。

四种双文件 diff 输出都通过了完整文件门禁，其中
`--diff-site-discordance` 有独立真实性能基准；另外三种 diff 输出尚
未单独声明性能加速比。

## vcftools-ng 新增参数

下列参数在 Original VCFtools 0.1.17 中不存在：

| 参数 | 作用 |
|---|---|
| `--threads N`、`-t N` | 设置总 CPU 预算；未指定时检测调度器/CPU affinity |
| `--batch-size N` | 调整通用有界流水线批大小 |
| `--input FILE` | 自动识别 VCF/BGZF VCF/BCF 的输入别名 |
| `--compat exact` | 显式选择当前唯一的精确兼容模式 |
| `--input-backend auto\|stream\|plain\|indexed` | 自动选择或强制输入后端 |
| `--bcftools FILE` | 指定自适应策略判定值得构建 CSI 时使用的 bcftools |
| `--recode-vcf-gz` | 并行生成确定性的 BGZF 压缩 VCF |

特别说明：Original 0.1.17 **没有** `--recode-vcf-gz`。该扩展通过将
解压后的内容与 Original `--recode` 逐字节比较来验证；不存在可供比较的
Original 压缩输出。

一次扫描同时生成多个统计输出也是 vcftools-ng 扩展；Original 通常需要
多次独立运行。

## Original 已知缺陷与有意边界

README 不会把 Original 的错误隐藏在“兼容”表述中：

- Original 的 BCF→VCF 路径会错误重写带引号逗号/等号的结构化 header，
  并具有特殊的缺失 GT（`-1/-1`）和 Character/String FORMAT
  尾字节/NUL 行为；`--compat exact` 为了逐字节一致会重现这些输出；
- Original 在 BCF 输入、基因型 masking 和 `--recode-bcf` 同时使用时
  会损坏基因型，vcftools-ng 会拒绝该组合；
- Original PCA 在存在缺失基因型时存在未定义/错位行为，精确 PCA 要求
  完整基因型（使用 `--max-missing 1`）；
- Original 的 `--non-ref-af-any` 单独使用时不生效，以及 diff 对部分
  缺失二倍体的行为，在精确模式中作为可观察兼容语义保留；
- 当前快速路径不支持多倍体 GT。

未来如果增加标准修正版行为，必须使用明确的新模式，不能静默改变
`--compat exact`。

## 使用示例

```bash
./bin/vcftools-ng \
  --bcf input.bcf \
  --threads 16 \
  --min-alleles 2 --max-alleles 2 \
  --minGQ 10 --minQ 30 --min-meanDP 7 \
  --max-missing 0.9 --maf 0.1 \
  --recode --recode-INFO-all \
  --out results/filtered
```

Plain VCF 不能建立 CSI/TBI，因为索引使用 BGZF 虚拟偏移。已有有效
`.csi`/`.tbi` 会被独立验证并保留；vcftools-ng 不覆盖用户已有索引。
v0.12.2 只把索引作为自适应加速手段：BGZF完整重编码在1线程时流式读取，
2线程及以上复用或构建索引；BCF完整重编码即使存在CSI也选择流式路径；
BGZF/BCF区域查询才优先复用或构建索引；紧凑型全文件统计从4线程起复用
已有索引，但不为单次扫描临时构建索引。`--no-auto-index` 已移除。
`--input-backend stream|indexed` 仍可用于高级诊断和显式覆盖。

## v0.12.1 完整数据第一轮发布门禁

全量发布工作负载使用七个真实项目过滤参数，并输出保留全部 INFO 的完整
VCF。32 CPU 主机上的 wall time（秒）如下：

| 场景 | Original | 1 线程 | 2 线程 | 4 线程 | 8 线程 | 16 线程 | 32 线程 |
|---|---:|---:|---:|---:|---:|---:|---:|
| BGZF + TBI | 2267.88 | 387.00 | 204.35 | 118.08 | 77.83 | 53.01 | 47.01 |
| BGZF + 自动 CSI | 2267.88 | 552.17 | 290.43 | 170.28 | 129.52 | 106.91 | 100.88 |
| BGZF + 禁止自动索引 | 2267.88 | 304.39 | 304.35 | 255.95 | 255.01 | 257.05 | 262.18 |
| Plain VCF | 2092.91 | 287.87 | 294.76 | 101.53 | 71.58 | 49.78 | 52.88 |
| BCF + CSI | 1943.47 | 321.28 | 323.02 | 162.08 | 109.14 | 58.58 | 42.08 |
| BCF + 自动 CSI | 1943.47 | 459.21 | 387.58 | 198.52 | 126.74 | 68.56 | 52.86 |
| BCF + 禁止自动索引 | 1943.47 | 319.97 | 161.73 | 109.67 | 56.70 | 41.24 | 40.41 |

相对 VCFtools 0.1.17 的加速倍率：

| 场景 | 1 线程 | 2 线程 | 4 线程 | 8 线程 | 16 线程 | 32 线程 |
|---|---:|---:|---:|---:|---:|---:|
| BGZF + TBI | 5.86× | 11.10× | 19.21× | 29.14× | 42.78× | 48.24× |
| BGZF + 自动 CSI | 4.11× | 7.81× | 13.32× | 17.51× | 21.21× | 22.48× |
| BGZF + 禁止自动索引 | 7.45× | 7.45× | 8.86× | 8.89× | 8.82× | 8.65× |
| Plain VCF | 7.27× | 7.10× | 20.61× | 29.24× | 42.04× | 39.58× |
| BCF + CSI | 6.05× | 6.02× | 11.99× | 17.81× | 33.18× | 46.19× |
| BCF + 自动 CSI | 4.23× | 5.01× | 9.79× | 15.33× | 28.35× | 36.77× |
| BCF + 禁止自动索引 | 6.07× | 12.02× | 17.72× | 34.28× | 47.13× | 48.09× |

42/42 个候选输出全部通过完整文件 `cmp`，单次实测加速范围为
4.11×–48.24×。自动 CSI 场景每次都从无 sidecar 的独立路径开始，因此
时间包含一次完整索引构建；正常使用时，成功生成的索引会保留并在后续
命令中复用。低线程只运行一次时构建索引可能不划算，高线程 BGZF 则可
明显受益。本次负载中，BGZF 从 2 线程开始自动 CSI 快于无索引，
BCF 则在所有测试线程数下均为无索引更快。这些 v0.12.1 数据促成了
v0.12.2 自适应策略：BCF完整文件重编码会自动选择流式路径，
`--no-auto-index` 不再接受。这些是单轮门禁数据，第 2–5 轮完成前不声明
最终均值或严格的相邻线程单调扩展。

## 从源码构建

源码构建需要 CMake、C++20 编译器、HTSlib、LAPACK、zlib 和 POSIX
threads：

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DHTSLIB_ROOT=/path/to/htslib
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## 测试与复现

- 固定开发测试流程：
  [docs/benchmark-workflow.md](docs/benchmark-workflow.md)
- 固定发布流程：
  [docs/release-workflow.md](docs/release-workflow.md)
- 三场景开发门禁：
  [benchmarks/run-development-gate.sh](benchmarks/run-development-gate.sh)
  （BGZF VCF + TBI、Plain VCF、BCF 自适应流式全扫描路径）
- v0.12.1 完整七场景发布驱动：
  [benchmarks/run-v0121-full-release-matrix.sh](benchmarks/run-v0121-full-release-matrix.sh)
- 参数兼容矩阵：
  [docs/parameter-compatibility.md](docs/parameter-compatibility.md)
- 版本历史：
  [docs/VERSION_HISTORY.md](docs/VERSION_HISTORY.md)

大体积真实输入和 Original golden 保留在基准主机本地；其大小、SHA-256、
命令、环境、时间、CPU、RSS 和精简结果表会随版本记录提交到 master。
