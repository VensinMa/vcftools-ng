# vcftools-ng

[English](README.md) | [简体中文](README.zh-CN.md)

vcftools-ng 是 VCFtools 0.1.17 的实验性高性能、输出兼容后继实现。

**最新正式版：**
[v0.13.0 — 事务式 BGZF 输出与热路径加速](https://github.com/VensinMa/vcftools-ng/releases/tag/v0.13.0)

v0.13.0 已加入科学输出事务式发布、日志与资源限制
加固，并默认输出更节省磁盘的 BGZF VCF。新的 ordered pipeline 会在运行
开始时一次编译不变的执行决策，将 DP/过滤样本扫描融合为一次，以 batch
连续块提交 VCF 文本，复用每个 BGZF worker 的压缩状态，并让 Plain VCF
worker 跨 shard 复用文件描述符。文本输入在严格总线程预算内采用偏向输入
阶段的自适应分配。便携构建已锁定 bcftools 1.24 与 HTSlib 1.24；完整数据
第一轮发布门禁已经 108/108 通过。v0.13.0 正式基准只采用这轮已验证的
单次时间；未完成的后续重复不会纳入任何公开计时结果。

全部科学结果先写入目标目录中的私有暂存文件，完成 flush、close 和写入
状态检查后才统一发布；运行失败会清理暂存文件并保留原有正式结果。日志
镜像写入失败只会关闭文件镜像，不会破坏 stderr 或科学输出。所有浮点参数
拒绝 NaN/Inf。未指定 `--threads` 时会对调度器、CPU affinity、cgroup 和
硬件限制求交，并将自动值限制为最多 128；阶段规划还会服从文件描述符
上限。显式线程值不受自动 128 上限限制，但仍服从实际资源分配上限。

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
- v0.13.0 将同一个直接文本内核扩展到 Plain/BGZF VCF 重编码，
  七个常用过滤参数与统计、重编码可以共享一次扫描；不支持的过滤、选择、
  输入格式或高级分析会在发布任何输出前自动回退到通用兼容流水线；
- 乱序 shard 的待提交输出带有背压，内存不会随剩余文件无限增长。

### 性能结论的适用范围

加速比只对实际测试的工作负载成立。v0.13.0 直接内核的性能数据适用于文档
指定的七参数位点统计/VCF 重编码组合，不能笼统外推到样本或位点选择、个体
归约、窗口 pi、Tajima's D、FST、LD、PCA、diff 或所有存储系统。这些操作
在已记录范围内仍保持精确兼容，也可能受益于公共流水线，但需要各自的扩展性
实测。固定的开发、候选版本和完整数据证据规则见
[代表性工作负载基准矩阵](docs/benchmark-workload-matrix.zh-CN.md)。

## 推荐安装方式

Linux x86_64 便携包解压后即可运行。`bin` 中只放 `vcftools-ng`；用于
自动构建 CSI 的私有 bcftools 放在 `libexec`。v0.13.0 便携构建同时
打包 bcftools 1.24、HTSlib 1.24 和所需的非 glibc 运行库：

```bash
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.13.0/vcftools-ng-v0.13.0-linux-x86_64.tar.gz
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.13.0/vcftools-ng-v0.13.0-linux-x86_64.tar.gz.sha256
sha256sum -c vcftools-ng-v0.13.0-linux-x86_64.tar.gz.sha256
tar -xzf vcftools-ng-v0.13.0-linux-x86_64.tar.gz
./vcftools-ng-v0.13.0-linux-x86_64/bin/vcftools-ng --version
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
| `--threads N`、`-t N` | 设置输入/计算/I/O 共享 worker 预算；自动检测各资源上限并最多取 128 |
| `--batch-size N` | 调整通用有界流水线批大小 |
| `--input FILE` | 自动识别 VCF/BGZF VCF/BCF 的输入别名 |
| `--compat exact` | 显式选择当前唯一的精确兼容模式 |
| `--input-backend auto\|stream\|plain\|indexed` | 自动选择或强制输入后端 |
| `--bcftools FILE` | 指定自适应策略判定值得构建 CSI 时使用的 bcftools |
| `--recode-vcf-gz` | 并行生成确定性的 BGZF 压缩 VCF |
| `--recode-vcf` | 在默认 `--recode` 改为 BGZF 后，显式生成未压缩 VCF |
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

独立的 [Original VCFtools 0.1.17 已知问题档案](docs/original-vcftools-0.1.17-known-issues.md)
持续记录每项问题的触发条件、Original 实际输出、vcftools-ng 处理策略和
回归证据。

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
实际输入 backend、各阶段线程分配、执行内核/组件、高层阶段耗时、全部输出
和过滤参数、样本/位点数、输出大小、wall/CPU 时间、峰值 RSS、警告、错误
和最终退出状态。已有有效
索引即使本次策略决定不使用，也会被明确记录并保留，绝不删除或覆盖。

```bash
vcftools-ng --gzvcf input.vcf.gz --threads 24 \
  --recode --out subset
# 数据：subset.recode.vcf.gz
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

### VCF 重编码输出

从 v0.13.0 源码候选版开始，`--recode` 默认直接写出 BGZF VCF；
`--recode-vcf-gz` 是同一输出的显式别名，两者都不会生成未压缩中间
文件。只有明确需要普通 VCF 文件时才使用新增的 `--recode-vcf`。

| 项目 | 行为 |
|---|---|
| 支持的输入 | `--vcf`、`--gzvcf`、`--bcf` 或自动识别的 `--input` |
| 默认压缩输出 | `--recode` 或 `--recode-vcf-gz` |
| 显式普通 VCF | `--recode-vcf` |
| 输出文件名 | `PREFIX.recode.vcf.gz`，`PREFIX` 来自 `--out` |
| INFO 字段 | 添加 `--recode-INFO-all` 以保留输入中的全部 INFO |
| 压缩方式 | 使用有效 `--threads` 预算的确定性 BGZF 压缩 |
| 输出索引 | 不会自动创建 |
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

如果确实同时需要未压缩和 BGZF 两种文件，可同时指定
`--recode-vcf` 与 `--recode-vcf-gz`：

```bash
vcftools-ng \
  --gzvcf input.vcf.gz \
  --threads 24 \
  --recode-vcf --recode-vcf-gz --recode-INFO-all \
  --out subset
```

第二条命令只执行一次过滤、解码和 VCF 格式化，同时写出：

```text
subset.recode.vcf
subset.recode.vcf.gz
```

同时写两种格式会增加输出 I/O 和压缩工作。为保持兼容，
`--recode --stdout` 仍输出普通 VCF stdout，不能再组合文件型 BGZF 输出。

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

## v0.13.0 完整数据第一轮性能

全量发布工作负载使用七个真实项目过滤参数并保留全部 INFO。下表是32 CPU
主机上通过门禁的第一轮 application wall time（秒），属于单次实测而非
多轮均值。Original 继续使用 v0.12.1 已锁定并校验哈希的基线，本轮没有
重跑。4种输入、3种输出/存储方式和9个线程档位共108/108项完整内容验证
全部通过。

### SSD 输出未压缩 VCF（`--recode-vcf`）

| 输入 | Original | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 2267.88 | 243.46 | 166.44 | 85.78 | 44.09 | 34.27 | 28.50 | 27.69 | 29.28 | 28.50 |
| BGZF VCF + 自动 CSI | 2267.88 | 241.42 | 264.60 | 140.10 | 100.80 | 91.84 | 85.66 | 82.55 | 87.44 | 85.75 |
| Plain VCF | 2092.91 | 241.25 | 126.93 | 69.34 | 42.25 | 41.53 | 41.99 | 40.01 | 40.39 | 40.70 |
| BCF 自适应 | 1943.47 | 317.31 | 159.73 | 160.44 | 82.43 | 48.10 | 43.13 | 35.91 | 35.65 | 37.02 |

生成等价格式未压缩 VCF 时，相对 VCFtools 0.1.17 的加速倍率：

| 输入 | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 9.32× | 13.63× | 26.44× | 51.44× | 66.18× | 79.57× | 81.90× | 77.45× | 79.57× |
| BGZF VCF + 自动 CSI | 9.39× | 8.57× | 16.19× | 22.50× | 24.69× | 26.48× | 27.47× | 25.94× | 26.45× |
| Plain VCF | 8.68× | 16.49× | 30.18× | 49.54× | 50.39× | 49.84× | 52.30× | 51.81× | 51.42× |
| BCF 自适应 | 6.12× | 12.17× | 12.11× | 23.58× | 40.40× | 45.06× | 54.12× | 54.52× | 52.50× |

### SSD 输出 BGZF VCF（`--recode` 或 `--recode-vcf-gz`）

| 输入 | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 1229.88 | 528.12 | 270.56 | 150.06 | 114.85 | 91.27 | 72.80 | 72.04 | 71.30 |
| BGZF VCF + 自动 CSI | 1221.37 | 622.72 | 324.19 | 205.69 | 171.43 | 149.57 | 128.60 | 127.71 | 128.74 |
| Plain VCF | 1018.63 | 509.09 | 257.79 | 134.33 | 104.59 | 87.78 | 66.65 | 67.24 | 67.14 |
| BCF 自适应 | 1019.48 | 511.56 | 265.43 | 160.21 | 123.11 | 98.45 | 83.54 | 80.16 | 81.39 |

相对 Original VCFtools 0.1.17 普通 VCF 工作流的端到端加速倍率：

| 输入 | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 1.84× | 4.29× | 8.38× | 15.11× | 19.75× | 24.85× | 31.15× | 31.48× | 31.81× |
| BGZF VCF + 自动 CSI | 1.86× | 3.64× | 7.00× | 11.03× | 13.23× | 15.16× | 17.64× | 17.76× | 17.62× |
| Plain VCF | 2.05× | 4.11× | 8.12× | 15.58× | 20.01× | 23.84× | 31.40× | 31.13× | 31.17× |
| BCF 自适应 | 1.91× | 3.80× | 7.32× | 12.13× | 15.79× | 19.74× | 23.26× | 24.24× | 23.88× |

### HDD 同一物理磁盘输入并输出 BGZF VCF

| 输入 | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 1322.12 | 524.01 | 270.81 | 152.87 | 115.65 | 95.05 | 73.03 | 72.35 | 71.47 |
| BGZF VCF + 自动 CSI | 1225.84 | 626.17 | 326.51 | 209.94 | 174.60 | 150.75 | 129.77 | 130.50 | 128.87 |
| Plain VCF | 1076.97 | 823.98 | 826.74 | 869.95 | 883.69 | 909.75 | 997.62 | 1029.14 | 1059.79 |
| BCF 自适应 | 1065.09 | 514.18 | 265.73 | 160.99 | 123.46 | 98.92 | 84.01 | 80.37 | 81.93 |

SSD BGZF 倍率比较的是同一过滤任务的完成时间，但输出编码不同：Original
写普通 VCF，vcftools-ng 写 BGZF；解压后的科学内容等价。59.43 GB 普通
VCF 压缩后为10.20 GB，减少82.8%。自动 CSI 行包含新建索引耗时。由于
保留的 Original 基线在 SSD 上测得，HDD 场景不报告会混入存储设备差异的
倍率。Plain VCF 在机械硬盘同盘读写时受到 I/O 上限约束，高线程反而变慢；
该结果被如实保留。

## 从源码构建

源码构建需要 CMake、C++20 编译器、HTSlib、LAPACK、zlib 和 POSIX
threads：

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DHTSLIB_ROOT=/path/to/htslib
cmake --build build -j
cmake --install build --prefix /path/to/install-prefix
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
- v0.13.0 输入/输出/存储完整数据发布门禁：
  [benchmarks/results/v0130-input-output-storage/README.md](benchmarks/results/v0130-input-output-storage/README.md)
- 参数兼容矩阵：
  [docs/parameter-compatibility.md](docs/parameter-compatibility.md)
- 版本历史：
  [docs/VERSION_HISTORY.md](docs/VERSION_HISTORY.md)

大体积真实输入和 Original golden 保留在基准主机本地；其大小、SHA-256、
命令、环境、时间、CPU、RSS 和精简结果表会随版本记录提交到 master。
