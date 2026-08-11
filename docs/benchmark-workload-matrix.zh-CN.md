# 代表性工作负载基准矩阵

[English](benchmark-workload-matrix.md)

性能结论只对实际测试的工作负载成立。某个输出、过滤组合、输入格式、
存储设备或样本选择密度得到的加速比，不能表述为所有 vcftools-ng 命令的
统一保证。

## v0.14.1 的实际范围

v0.13.0 的直接文本内核并非只支持统计。对于满足条件的 Plain/BGZF VCF，
它可以把七参数过滤、位点局部统计和 VCF 重编码融合到一次扫描。七个过滤
参数是 `--min-alleles`、`--max-alleles`、`--minQ`、`--minGQ`、
`--min-meanDP`、`--max-missing` 和 `--maf`；可融合的位点输出是
`--freq`、`--freq2`、`--counts`、`--missing-site`、`--site-depth`、
`--site-mean-depth` 和 `--site-quality`。

v0.13.1 进一步让满足条件的 Plain VCF 位点保留/剔除、样本投影、窗口 pi、
Tajima's D、site FST 和 window FST 进入直接文本执行族。这些路径共享选中
样本 GT 解码、紧凑逐位点贡献、确定性有序归约和自适应只读映射范围。
v0.13.2将十参数生产过滤、FILTER/INFO/FT和共享site pi加入该执行族；无索引
BGZF加入有界解压/计算重叠；满足条件的LD、精确PCA和有索引BCF discordance
使用专用内核。不满足条件的形态仍回退通用兼容流水线。

v0.14.1将每次运行编译为一个不可变能力计划，去掉GT-only分析中未使用的
FORMAT字段工作，优化确定性LD/PCA后处理存储，并加固异常输入、故障和倍性
边界；上述未正式发布的v0.13.2工作全部纳入v0.14.1。

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
v0.13.1锁定230,000位点SSD/NVMe W03-W10矩阵；v0.13.2为新增内核补充稳定
230k A/B，避免把23k中亚秒级
启动噪声当作吞吐结论。更大候选版本使用标准2,300,000位点真实子集；
v0.14.1稳定本地扩展线程集合为`1 2 4 8 12 16 24 28 32`。其
11,230,392位点发布门禁在四种代表输入场景（BGZF+TBI、BGZF+automatic
CSI、Plain VCF、BCF自适应stream）下对同一七参数精确重编码负载测试上述
九个线程数。

可复用驱动是
[`benchmarks/run-workload-matrix.sh`](../benchmarks/run-workload-matrix.sh)。
把
[`benchmarks/workload-matrix-profile.example.sh`](../benchmarks/workload-matrix-profile.example.sh)
复制到仓库外，填写锁定输入、选择列表、群体文件、oracle 和结果目录。
驱动在运行前校验 `ORACLE_ROOT/SHA256SUMS`，不会自动生成或覆盖 Original
oracle。

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

v0.13.1 已提交矩阵、精简计时、oracle/输入哈希和固定runner见
[`benchmarks/results/workload-matrix-230k-v0130/RESULTS.zh-CN.md`](../benchmarks/results/workload-matrix-230k-v0130/RESULTS.zh-CN.md)。
v0.13.2九类精确门禁、oracle哈希和A/B摘要见
[`benchmarks/results/v0132-development-gate/README.md`](../benchmarks/results/v0132-development-gate/README.md)。
v0.14.1完整数据发布驱动和精简结果分别位于
[`benchmarks/run-v0141-full-release-matrix.sh`](../benchmarks/run-v0141-full-release-matrix.sh)
与`benchmarks/results/final-full-v0141/`。
