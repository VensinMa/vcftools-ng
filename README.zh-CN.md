# vcftools-ng

[English](README.md) | [简体中文](README.zh-CN.md)

vcftools-ng 是 VCFtools 0.1.17 的实验性高性能、输出兼容后继实现。

**最新正式版：**
[v0.11.3](https://github.com/VensinMa/vcftools-ng/releases/tag/v0.11.3)

**发布候选版：** v0.12.1 — 融合位点统计与可扩展精确重编码。固定三场景开发
门禁已经通过，完整七场景发布测试正在进行。

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
| `--no-auto-index` | 禁止自动构建 CSI |
| `--bcftools FILE` | 指定自动索引使用的 bcftools |
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
使用 `--no-auto-index` 可保留无索引流式路径。

## v0.12.1 开发门禁结果

在 230 万真实位点、412 个样本和七个真实过滤参数下，固定 Original
基线未重新运行。三种输入的 1/2/4/8/16/32 线程共 18 个输出全部逐字节
通过：

| 场景 | Original | 1线程 | 2线程 | 4线程 | 8线程 | 16线程 | 32线程 |
|---|---:|---:|---:|---:|---:|---:|---:|
| BGZF + TBI | 392.73s | 74.04s | 40.77s | 23.13s | 15.48s | 10.56s | 8.44s |
| Plain VCF | 353.46s | 54.27s | 54.10s | 15.80s | 11.32s | 8.77s | 8.45s |
| BCF + CSI | 327.22s | 65.70s | 65.58s | 33.15s | 22.38s | 11.80s | 8.26s |

正式发布前还必须完成 1123 万位点、七场景、五次重复的完整矩阵。

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
- v0.12.1 完整七场景发布驱动：
  [benchmarks/run-v0121-full-release-matrix.sh](benchmarks/run-v0121-full-release-matrix.sh)
- 参数兼容矩阵：
  [docs/parameter-compatibility.md](docs/parameter-compatibility.md)
- 版本历史：
  [docs/VERSION_HISTORY.md](docs/VERSION_HISTORY.md)

大体积真实输入和 Original golden 保留在基准主机本地；其大小、SHA-256、
命令、环境、时间、CPU、RSS 和精简结果表会随版本记录提交到 master。
