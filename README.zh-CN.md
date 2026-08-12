# vcftools-ng

[English](README.md) · [英文详细技术文档](TECHNICAL_REFERENCE.md) · [中文详细技术文档](TECHNICAL_REFERENCE.zh-CN.md)

vcftools-ng 是面向常用 VCF 过滤、统计、群体遗传分析和格式重编码任务的
高性能 VCFtools 0.1.17 兼容实现。

**推荐版本：** [v0.14.2 — Portable I/O Recovery](https://github.com/VensinMa/vcftools-ng/releases/tag/v0.14.2)

它保留 VCFtools 风格的命令，并在明确声明兼容的范围内保证科学结果一致；
内部使用 HTSlib、按工作负载选择的专用解析器、有界并行流水线和自适应输入
后端。实际速度取决于参数、输入格式、存储设备、样本数和线程数，因此下表是
真实工作负载的实测证据，不是对所有命令的统一加速承诺。

## 哪些任务更快？

以下参数组均存在于 Original VCFtools 0.1.17，并且已经通过真实数据逐字节门禁
和性能测试。

| 工作负载 | 主要兼容参数 | 代表性实测结果 |
|---|---|---:|
| 一次扫描生成位点统计 | `--freq`、`--freq2`、`--counts`、`--missing-site`、`--site-depth`、`--site-mean-depth`、`--site-quality` | 六项同时输出：8/16 线程 **38.68x / 46.11x** |
| 个体与 HWE 统计 | `--depth`、`--missing-indv`、`--het`、`--hardy`、`--site-quality` | 单项在 8/16 线程为 **13.02x–21.17x**；五项一次扫描在 16 线程为 **39.33x** |
| 位点和样本选择 | `--positions`、`--exclude-positions`、`--keep`、`--remove`、`--indv`、`--remove-indv` | 固定 23 万位点矩阵在 32 线程为 **8.17x–27.93x** |
| π、Tajima 与 FST | `--site-pi`、`--window-pi`、`--window-pi-step`、`--TajimaD`、`--weir-fst-pop`、FST 窗口参数 | 代表性 23 万位点任务在 32 线程为 **20.15x–29.21x** |
| LD 与 PCA | `--geno-r2` 及其窗口参数、`--pca`、`--pca-no-norm` | 已有结果一致的专用内核；加速依赖数据，PCA 提升可能较小 |
| 过滤并重编码 | 常用位点、基因型、FILTER/INFO、区间和样本过滤与 `--recode` | v0.14.2 完整数据、等价普通 VCF 输出：索引 BGZF 输入在 1–32 线程为 **10.58x–80.28x** |
| 格式转换与差异比较 | `--recode-bcf` 和支持的 `--diff*` 输出 | BCF 转换在 16 线程为 **18.89x**；位点 discordance 为 **15.14x** |

群体/选择结果是三次运行中位数；v0.14.2 完整数据重编码是 11,230,392 位点的
单轮发布证据，v0.13.0/v0.14.1/v0.14.2 三版本共 81 行便携包 A/B 全部通过
完整文件一致性门禁。不同工作负载之间不要直接外推，详细条件见
[基准适用范围](docs/benchmark-workload-matrix.zh-CN.md)和
[v0.14.2 完整数据 A/B 证据](benchmarks/results/full-unified-v0142-ab/README.zh-CN.md)。

### 完整数据过滤并输出普通 VCF

固定 11,230,392 位点工作负载相对 Original VCFtools 0.1.17 的加速倍率：

| 输入 | 1 线程 | 8 线程 | 16 线程 | 28 线程 | 32 线程 |
|---|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 10.58x | 55.69x | 71.59x | 78.64x | 80.28x |
| Plain VCF | 9.77x | 64.38x | 65.30x | 63.96x | 65.98x |
| BCF 自适应流式路径 | 3.52x | 19.79x | 37.16x | 43.90x | 38.54x |

此表是 v0.14.2 的 `--recode-vcf` 结果，因此输出容器与 Original 的未压缩
VCF 等价。Original 没有重跑，倍率分母是相同输入和过滤命令下保留并经哈希
验证的 v0.12.1 基线。默认 BGZF 输出占用空间小得多，但压缩成本不同；相关结果见
[SSD/HDD 输出表](benchmarks/results/v0130-input-output-storage/README.md)。
性能差异在 5% 以内视为基本持平，但科学输出始终必须完全一致。

## 安装

Linux x86_64 便携包解压即用，要求 glibc 2.17 或更新版本，并已在干净的
CentOS 7 和 Ubuntu 20.04 容器中验证。`bin` 中只有 `vcftools-ng`，私有运行
工具位于 `libexec`，不会污染用户的 PATH。

```bash
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.14.2/vcftools-ng-v0.14.2-linux-x86_64.tar.gz
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.14.2/vcftools-ng-v0.14.2-linux-x86_64.tar.gz.sha256
sha256sum -c vcftools-ng-v0.14.2-linux-x86_64.tar.gz.sha256
tar -xzf vcftools-ng-v0.14.2-linux-x86_64.tar.gz
./vcftools-ng-v0.14.2-linux-x86_64/bin/vcftools-ng --help
```

不需要编译器、CMake、Conda 环境、系统 HTSlib 或系统 bcftools。解压后的
`bin`、`lib`、`libexec` 目录需要放在一起。

## 快速使用

### 过滤并输出压缩 VCF（推荐）

```bash
vcftools-ng --gzvcf input.vcf.gz --threads 24 \
  --min-alleles 2 --max-alleles 2 --remove-indels \
  --minQ 40 --minGQ 20 --minDP 5 --maxDP 30 \
  --min-meanDP 10 --max-missing 0.9 --maf 0.1 \
  --recode --recode-INFO-all --out filtered
```

生成 `filtered.recode.vcf.gz` 和 `filtered.log`。压缩结果是 BGZF VCF；当前
版本不会自动给新输出建立索引，如有需要再执行：

```bash
bcftools index --tbi --threads 24 filtered.recode.vcf.gz
```

### 一次扫描生成多种统计结果

```bash
vcftools-ng --gzvcf input.vcf.gz --threads 16 \
  --freq --counts --missing-site --site-depth --site-mean-depth \
  --site-quality --out cohort
```

推荐把兼容的输出放在一次命令中：Original 通常每个统计命令扫描一次输入，
vcftools-ng 可以共享解析和基因型解码。

### 个体统计

```bash
vcftools-ng --gzvcf input.vcf.gz --threads 16 \
  --depth --missing-indv --het --hardy --out cohort
```

### 窗口 π、Tajima's D 和 FST

```bash
vcftools-ng --gzvcf input.vcf.gz --threads 16 \
  --window-pi 100000 --window-pi-step 10000 --TajimaD 100000 \
  --weir-fst-pop population1.txt --weir-fst-pop population2.txt \
  --fst-window-size 100000 --fst-window-step 10000 --out diversity
```

窗口分析要求同一染色体内位置非递减，并且染色体区段连续。每个群体文件每行
放一个样本 ID。

## 与 Original VCFtools 的重要差异

| 行为 | vcftools-ng 默认值 | 如何修改 |
|---|---|---|
| 文件形式的 `--recode` | 确定性 BGZF VCF：`PREFIX.recode.vcf.gz` | 用 `--recode-vcf` 输出未压缩的 `PREFIX.recode.vcf` |
| 输入后端/索引 | 根据格式、工作负载、存储和有效线程自适应 | 高级诊断：`--input-backend stream\|plain\|indexed` |
| 线程数 | 综合调度器、CPU affinity、cgroup 和硬件限制，自动模式最高 128 | 用 `--threads N` / `-t N` 指定有效资源范围内的预算 |
| 运行日志 | 默认写入 `PREFIX.log`，终端信息继续保留 | `--log-file FILE` 或 `--no-log-file` |
| 多种输出 | 能兼容的分析共享同一次扫描 | 在一条命令中列出全部所需输出参数 |
| 运行失败 | 科学结果先暂存，全部成功后事务式发布 | 无需设置；失败时保留原有目标文件 |

`--recode-vcf-gz` 是显式生成压缩 VCF 的参数，Original 0.1.17 中不存在；
`--recode-vcf` 也是 vcftools-ng 新参数。使用 `--recode --stdout` 时仍输出普通
VCF 到 stdout，所有日志只进入 stderr 和日志文件，不会污染数据流。

自适应索引策略绝不覆盖已有的有效 TBI/CSI。全文件 BGZF 重编码通常在两个及
以上有效线程时选择索引并行；缺少输入索引时，只在预期有收益时构建。全文件
BCF 即使已有 CSI，通常仍选择更快的流式路径；限制染色体/坐标的 BCF 查询才
优先利用索引。日常使用建议保持 `--input-backend auto`。

## 已优化的 Original 兼容参数

下面只按功能组简要列出已有逐字节与性能证据的 Original 参数，不表示所有参数
笛卡尔组合都已经测试。

- 统计：`--freq`、`--freq2`、`--counts`、`--missing-site`、深度与质量统计、
  个体统计、HWE、π、Tajima、FST、LD 和 PCA。
- 位点过滤：等位数、SNP/indel、QUAL、平均深度、缺失率、MAF/MAC/HWE、
  non-reference、染色体/位置/BED、FILTER/INFO Flag 和 `--thin`。
- 基因型与样本过滤：`--minGQ`、`--minDP`、`--maxDP`、基因型 FT、
  `--keep`、`--remove`、`--indv`、`--remove-indv`。
- 输出与比较：`--recode`、`--recode-bcf`、`--recode-INFO-all`、`--stdout`
  和已支持的 `--diff*` 参数组。

精确参数名、验证证据、不支持的组合、多倍体边界和刻意继承的 Original 历史
行为见[参数兼容矩阵](docs/parameter-compatibility.md)。

## 常用 vcftools-ng 特有参数

| 参数 | 默认状态 | 用途 |
|---|---|---|
| `--threads N`、`-t N` | 自动，最高 128 | 设置全流程共享 CPU 预算 |
| `--input FILE` | — | 自动识别 VCF、BGZF VCF 或 BCF |
| `--recode-vcf-gz` | 关闭 | 显式生成确定性 BGZF VCF |
| `--recode-vcf` | 关闭 | 显式生成未压缩 VCF |
| `--log-file FILE` | `PREFIX.log` | 自定义运行日志路径 |
| `--no-log-file` | 关闭 | 只关闭日志文件，不关闭 stderr |
| `--input-backend ...` | `auto` | 专家级后端覆盖/诊断 |
| `--corrected-depth-arithmetic` | 关闭 | 用检查溢出的 64 位深度求和替代 Original 的 32 位回绕行为 |

运行 `vcftools-ng --help` 可查看完整彩色终端手册、参数组合、输出后缀、默认值
和示例。

## 推荐设置

- 日常使用保持 `--input-backend auto`。
- 调度器/cgroup 配置正确时可以不写 `--threads`；否则应填写作业实际获得的
  CPU 数，而不是服务器的总核心数。
- 优先使用 `--recode` 输出 BGZF，避免产生超大的普通 VCF，机械硬盘上尤其
  重要。只有下游明确要求未压缩 VCF 时才使用 `--recode-vcf`。
- 需要保留全部 INFO 注释时添加 `--recode-INFO-all`。
- 默认保留日志以便复现，并把兼容的统计输出合并到同一次运行。

## 文档与证据入口

科学兼容指文档覆盖的工作负载输出与保留的 Original VCFtools 0.1.17 golden
逐字节一致。由于 Original 没有压缩 VCF 输出参数，BGZF 结果按解压后的 VCF
字节比较。不满足专用内核条件的输入会回退到通用精确路径或明确拒绝，不会
静默进入不兼容的快速内核。

| 想了解的内容 | 独立文档 |
|---|---|
| 原首页长版的全部内容 | [中文详细技术文档](TECHNICAL_REFERENCE.zh-CN.md) |
| 精确参数名、优化状态和兼容边界 | [参数兼容矩阵](docs/parameter-compatibility.md) |
| 发布历史与各版本变化 | [发布历史](docs/VERSION_HISTORY.md)和[逐版本记录](docs/versions/README.md) |
| 逐字节/哈希门禁和固定测试制度 | [基准与一致性流程](docs/benchmark-workflow.md) |
| 每项性能结论由哪个工作负载支持 | [代表性工作负载矩阵](docs/benchmark-workload-matrix.zh-CN.md) |
| 多版本完整结果和原始基准记录 | [基准档案](benchmarks/README.md) |
| v0.14.2 完整数据 A/B 与结果门禁 | [v0.14.2 证据](benchmarks/results/full-unified-v0142-ab/README.zh-CN.md) |
| 从源码构建和验证命令 | [构建与验证](TECHNICAL_REFERENCE.zh-CN.md#从源码构建) |
| 输入/索引调度和 capability planner | [自适应输入后端](docs/architecture/adaptive-input-backends.md)和 [QueryPlan](docs/architecture/query-plan.md) |
| Release 如何测试、打包与发布 | [发布流程](docs/release-workflow.md) |
| Original 行为中哪些被保留、修正或拒绝 | [Original VCFtools 已知问题记录](docs/original-vcftools-0.1.17-known-issues.md) |

[文档总索引](docs/README.md)把这些入口集中在一个页面，并继续链接到原始 TSV
数据和复现脚本。

许可证：LGPL-3.0-or-later，见 [LICENSE](LICENSE)。
