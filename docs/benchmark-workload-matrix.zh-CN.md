# 代表性工作负载基准矩阵

[English](benchmark-workload-matrix.md)

性能结论只对实际测试的工作负载成立。某个输出、过滤组合、输入格式、
存储设备或样本选择密度得到的加速比，不能表述为所有 vcftools-ng 命令的
统一保证。

## v0.13.0 的实际范围

v0.13.0 的直接文本内核并非只支持统计。对于满足条件的 Plain/BGZF VCF，
它可以把七参数过滤、位点局部统计和 VCF 重编码融合到一次扫描。七个过滤
参数是 `--min-alleles`、`--max-alleles`、`--minQ`、`--minGQ`、
`--min-meanDP`、`--max-missing` 和 `--maf`；可融合的位点输出是
`--freq`、`--freq2`、`--counts`、`--missing-site`、`--site-depth`、
`--site-mean-depth` 和 `--site-quality`。

样本/位点选择、FILTER/INFO/FT 选择、尚未融合的过滤、BCF、个体归约、
窗口统计、FST、LD、PCA 和 diff 使用通用兼容流水线。它们可以受益于公共
输入、调度和输出优化，但不能直接套用七参数融合内核的加速比。

每次运行日志记录 `Execution kernel`、`Execution components`、输入后端、
线程分配和高层阶段耗时。基准结果必须同时记录这些字段，让退出 fast path
的情况可以直接识别。

## 开发基准矩阵

矩阵按计算结构覆盖代表性工作负载，不做全部参数的笛卡尔积：

| 编号 | 工作负载 | 必测变化 |
|---|---|---|
| W01 | 无过滤位点计数 | `--counts` |
| W02 | 生产环境七参数过滤计数 | 固定真实过滤组合 |
| W03 | 位点保留 | 约 1%/50%，已排序/打乱列表 |
| W04 | 位点剔除 | 约 1%/50%，列表包含重复和不存在位点 |
| W05 | 样本保留 | 约 25%/50%/100% 样本并输出 counts |
| W06 | 样本选择后重编码 | 保留 50% 样本并使用真实过滤组合 |
| W07 | 窗口 pi | 生产环境重叠窗口和非重叠窗口 |
| W08 | Tajima's D | 生产环境窗口大小 |
| W09 | site FST | 真实群体中样本最少和最多的群体对 |
| W10 | window FST | 真实群体对、重叠窗口，另加允许多等位场景 |
| W11 | 输出密集重编码 | 能支持时分别输出 Plain VCF、BGZF VCF 和 BCF |
| W12 | 通用多等位路径 | 不限制二等位的位点统计 |

当前桂花项目生产参数为：pi 和 window FST 使用 100 kb 窗口、10 kb 步长，
Tajima's D 使用 100 kb 窗口。该值来自当前项目配置；以后修改时应记录为
新的 profile，不能静默改变固定基准。

开发阶段使用 23,000 个真实位点，线程 `1 4 8 16 32`，最多重复三次。
候选版本使用标准 2,300,000 位点真实子集，线程 `1 2 4 8 16 32`，至少
重复三次。11,230,392 位点最终门禁只选四类：W02、W06、W07/W08 之一和
W10。

## 一致性与固定基线

- Original VCFtools 0.1.17 的 golden 和时间只生成一次，保存哈希并长期
  复用；日常优化不重新运行 Original。
- 第一轮执行完整逐字节门禁，通过后后续重复才只用于计时。
- 文本结果使用 `cmp`；BGZF 解压后与 Original VCF oracle 比较；BCF 使用
  现有的规范化兼容流程。
- exact 模式下窗口和群体统计仍要求逐字节一致，数值容差不能代替门禁。
- 结果记录重复次数、wall、CPU、RSS、输出大小、kernel、backend、样本数、
  输入/保留位点数以及输出行数或窗口数。

## 阶段耗时的解释

内置计时保持低开销，只记录高层阶段：

- `input/index planning`：格式/索引检查和输入 source 建立；
- `pipeline setup`：位点选择器、样本映射和输出构建；
- `ordered input/compute/commit`：通用流水线中重叠的读取、解码、分析和有序提交；
- `fused scan/filter/output`：直接融合文本内核；
- `output finalization`：最终 flush、关闭和 writer 校验。

解析、reducer、压缩和 ordered-commit 等待的细分归因应放在专门 profiling
构建中。正式 Release 基准不能启用逐记录计时，否则计时开销会扭曲性能。
