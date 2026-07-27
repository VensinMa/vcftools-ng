# vcftools-ng

[English](README.md) | [简体中文](README.zh-CN.md)

vcftools-ng 是 VCFtools 0.1.17 的实验性高性能、输出兼容后继实现。

**最新正式版：**
[v0.12.4 — 标准可复现运行日志](https://github.com/VensinMa/vcftools-ng/releases/tag/v0.12.4)

v0.12.4 默认生成 `PREFIX.log`，记录完整命令、输入输出、过滤参数、线程
分配、资源使用、CSI/TBI 验证、自适应索引决定及原因、耗时、警告和最终
状态。`--log-file FILE` 可自定义路径，`--no-log-file` 只关闭日志文件；
科学数据 stdout 始终保持纯净。

v0.12.3 引入的完整彩色终端手册继续保留，可通过 `vcftools-ng --help`
查看全部支持参数、示例、输出后缀、组合限制和自适应后端策略。

完整 11,230,392 位点的四场景发布矩阵已经通过：
1/2/4/8/16/32 线程、每个配置五次重复，共 120/120 个 vcftools-ng
输出均与保留的 VCFtools 0.1.17 golden 逐字节一致。v0.12.2 复用了
v0.12.1 已锁定并重新校验哈希的 Original 时间与 golden，本轮没有重跑
Original。v0.12.4 不改变科学过滤、统计、输入调度或输出字节，因此继承
这组已锁定的正式发布矩阵，不将旧数据重新标成新跑基准。启用标准日志后
另外完成了 230 万位点三场景门禁，18 个场景/线程组合全部逐字节通过。

vcftools-ng 使用自适应、有序的输入分片和有界流水线：

- Plain VCF 在 1–2 线程时流式读取，3 线程及以上使用按完整记录边界
  对齐的并行字节范围；
- BGZF 完整重编码在 1 线程时流式读取，2 线程及以上复用或构建索引；
- BCF 完整重编码即使存在 CSI 也使用流式路径，选择性区域查询才复用或
  构建索引；
- 输入、计算和有序输出可以重叠，结果顺序保持确定；
- `--freq`、`--freq2`、`--counts`、`--missing-site`、
  `--site-depth`、`--site-mean-depth` 和 `--site-quality` 可以进入
  直接融合文本路径，避免构建中间 `bcf1_t`；
- 乱序 shard 的待提交输出带有背压，内存不会随剩余文件无限增长。

## 推荐安装方式

Linux x86_64 便携包解压后即可运行，包内包含 vcftools-ng、用于自动
构建 CSI 的 bcftools、HTSlib 以及所需的非 glibc 运行库：

```bash
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.12.4/vcftools-ng-v0.12.4-linux-x86_64.tar.gz
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.12.4/vcftools-ng-v0.12.4-linux-x86_64.tar.gz.sha256
sha256sum -c vcftools-ng-v0.12.4-linux-x86_64.tar.gz.sha256
tar -xzf vcftools-ng-v0.12.4-linux-x86_64.tar.gz
./vcftools-ng-v0.12.4-linux-x86_64/bin/vcftools-ng --version
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
| `--log-file FILE` | 自定义运行日志路径并覆盖写入 |
| `--no-log-file` | 关闭日志文件，但保留终端 stderr 信息 |

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

## 命令行帮助

运行 `vcftools-ng -h` 或 `vcftools-ng --help` 可以查看完整终端手册，
其中包含基础示例、所有已支持参数的用途、输出文件后缀、参数组合限制和
自适应输入/索引策略。

交互式终端会为分组标题、参数、示例和重要提示自动着色；通过管道或重定向
保存时仍输出无 ANSI 转义码的纯文本：

```bash
vcftools-ng --help
vcftools-ng --help > vcftools-ng-help.txt
NO_COLOR=1 vcftools-ng --help
CLICOLOR_FORCE=1 vcftools-ng --help
```

## 标准运行日志

普通运行默认生成 `PREFIX.log`，延续 Original VCFtools 的输出前缀习惯，
并增加完整的复现元数据：

| 调用方式 | 日志行为 |
|---|---|
| `--out subset` | 覆盖写入 `subset.log` |
| `--out results/sample` | 覆盖写入 `results/sample.log` |
| 未指定 `--out` | 覆盖写入 `out.log` |
| `--log-file FILE` | 覆盖写入用户指定路径 |
| `--no-log-file` | 不生成日志文件，但终端信息仍写入 stderr |

终端与日志文件由同一个中心 logger 输出完全相同的诊断信息。使用
`--recode --stdout` 时，VCF 字节只进入 stdout；日志始终只进入 stderr
和日志文件，不会污染 VCF。

日志记录完整命令、工作目录、开始/结束时间、输入格式/大小/存储类型、
已有 TBI/CSI 及验证结果、自适应索引决定和原因、建索引线程数与耗时、
实际输入 backend、各阶段线程分配、全部输出和过滤参数、样本/位点数、
输出大小、wall/CPU 时间、峰值 RSS、警告、错误和最终退出状态。已有有效
索引即使本次策略决定不使用，也会被明确记录并保留，绝不删除或覆盖。

```bash
vcftools-ng --gzvcf input.vcf.gz --threads 24 \
  --recode --out subset
# 数据：subset.recode.vcf
# 日志：subset.log

vcftools-ng --gzvcf input.vcf.gz --counts \
  --log-file logs/counts-run.log --out counts
```

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

### 直接生成 BGZF VCF

`--recode-vcf-gz` 是 vcftools-ng 新增的输出参数，可以将过滤后的记录
直接写为 BGZF VCF，不会先生成未压缩中间文件，也不要求同时指定
`--recode`。

| 项目 | 行为 |
|---|---|
| 支持的输入 | `--vcf`、`--gzvcf`、`--bcf` 或自动识别的 `--input` |
| 只输出压缩 VCF | 使用 `--recode-vcf-gz`，不加 `--recode` |
| 输出文件名 | `PREFIX.recode.vcf.gz`，`PREFIX` 来自 `--out` |
| INFO 字段 | 添加 `--recode-INFO-all` 以保留输入中的全部 INFO |
| 压缩方式 | 使用有效 `--threads` 预算的确定性 BGZF 压缩 |
| 输出索引 | v0.12.4 不会自动创建 |
| 兼容性 | 解压内容与 Original `--recode` 做逐字节比较 |

```bash
vcftools-ng \
  --gzvcf input.vcf.gz \
  --threads 24 \
  --min-alleles 2 --max-alleles 2 \
  --minQ 40 --minGQ 20 --minDP 5 --maxDP 30 \
  --min-meanDP 10 --max-missing 0.9 --maf 0.1 \
  --recode-vcf-gz --recode-INFO-all \
  --out subset
```

该命令生成：

```text
subset.recode.vcf.gz
```

`--recode` 与 `--recode-vcf-gz` 不互斥。如果确实同时需要未压缩和
BGZF 两种文件，可以一起指定：

```bash
vcftools-ng \
  --gzvcf input.vcf.gz \
  --threads 24 \
  --recode --recode-vcf-gz --recode-INFO-all \
  --out subset
```

第二条命令只执行一次过滤、解码和 VCF 格式化，同时写出：

```text
subset.recode.vcf
subset.recode.vcf.gz
```

同时写两种格式会增加输出 I/O 和压缩工作。`--stdout` 仅支持单独的
`--recode`，不能与 `--recode-vcf-gz` 组合。

在已经测试的线程数下，压缩输出字节保持确定性；解压后的 VCF 已经与
Original VCFtools 0.1.17 的 `--recode` 输出通过完整文件逐字节比较。
Original 本身没有 `--recode-vcf-gz` 参数。目前尚未记录严格同条件的
`--recode` 与 `--recode-vcf-gz` 性能对照。

如果下游需要索引访问，可在过滤完成后创建 TBI：

```bash
bcftools index --tbi --threads 24 subset.recode.vcf.gz
```

Plain VCF 不能建立 CSI/TBI，因为索引使用 BGZF 虚拟偏移。已有有效
`.csi`/`.tbi` 会被独立验证并保留；vcftools-ng 不覆盖用户已有索引。
v0.12.2 只把索引作为自适应加速手段：BGZF完整重编码在1线程时流式读取，
2线程及以上复用或构建索引；BCF完整重编码即使存在CSI也选择流式路径；
BGZF/BCF区域查询才优先复用或构建索引；紧凑型全文件统计从4线程起复用
已有索引，但不为单次扫描临时构建索引。`--no-auto-index` 已移除。
`--input-backend stream|indexed` 仍可用于高级诊断和显式覆盖。

## 继承的 v0.12.2 完整数据五轮发布矩阵

全量发布工作负载使用七个真实项目过滤参数，并输出保留全部 INFO 的完整
VCF。下表为 vcftools-ng 在 32 CPU 主机上严格串行运行五次的 wall time
均值（秒）。Original 为 v0.12.1 已锁定的单次基线，本轮没有重跑。
v0.12.4 的标准日志功能不改变科学执行路径或输出字节，因此继承这组矩阵。

| 场景 | Original | 1 线程 | 2 线程 | 4 线程 | 8 线程 | 16 线程 | 32 线程 |
|---|---:|---:|---:|---:|---:|---:|---:|
| BGZF + TBI | 2267.88 | 309.28 | 204.91 | 112.35 | 73.93 | 52.05 | 42.32 |
| BGZF + 自动 CSI | 2267.88 | 308.21 | 281.03 | 167.50 | 128.60 | 107.16 | 99.59 |
| Plain VCF | 2092.91 | 285.45 | 285.22 | 89.31 | 63.71 | 47.84 | 50.61 |
| BCF 自适应流式 | 1943.47 | 325.66 | 164.14 | 110.57 | 57.26 | 41.30 | 40.20 |

相对 VCFtools 0.1.17 的加速倍率：

| 场景 | 1 线程 | 2 线程 | 4 线程 | 8 线程 | 16 线程 | 32 线程 |
|---|---:|---:|---:|---:|---:|---:|
| BGZF + TBI | 7.33× | 11.07× | 20.19× | 30.68× | 43.57× | 53.59× |
| BGZF + 自动 CSI | 7.36× | 8.07× | 13.54× | 17.63× | 21.16× | 22.77× |
| Plain VCF | 7.33× | 7.34× | 23.44× | 32.85× | 43.75× | 41.35× |
| BCF 自适应流式 | 5.97× | 11.84× | 17.58× | 33.94× | 47.05× | 48.35× |

120/120 个输出全部通过完整文件 `cmp` 和 SHA 校验，五次均值加速范围为
5.97×–53.59×。自动 CSI 场景每轮都包含一次新建索引的成本；BCF
完整重编码在所有线程档位均自动选择流式路径，区域查询仍会使用 CSI。

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
- v0.12.2 完整四场景发布驱动：
  [benchmarks/run-v0122-full-release-matrix.sh](benchmarks/run-v0122-full-release-matrix.sh)
- v0.12.4 技术记录：
  [docs/versions/v0.12.4.md](docs/versions/v0.12.4.md)
- v0.12.4 启用日志后的三场景门禁：
  [benchmarks/results/development-v0124-logging-final/README.md](benchmarks/results/development-v0124-logging-final/README.md)
- 参数兼容矩阵：
  [docs/parameter-compatibility.md](docs/parameter-compatibility.md)
- 版本历史：
  [docs/VERSION_HISTORY.md](docs/VERSION_HISTORY.md)

大体积真实输入和 Original golden 保留在基准主机本地；其大小、SHA-256、
命令、环境、时间、CPU、RSS 和精简结果表会随版本记录提交到 master。
