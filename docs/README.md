# vcftools-ng documentation index

[Project quick start](../README.md) · [中文快速开始](../README.zh-CN.md) ·
[Detailed English reference](../TECHNICAL_REFERENCE.md) ·
[中文详细技术文档](../TECHNICAL_REFERENCE.zh-CN.md)

This page separates user guidance from engineering and release evidence. The
root README is intentionally short; no historical benchmark, gate, or design
record was discarded when it was shortened.

## User documentation / 用户文档

- [Quick start and recommended settings / 快速开始与推荐设置](../README.md)
  ([中文](../README.zh-CN.md))
- [Parameter compatibility and optimization status / 参数兼容与优化状态](parameter-compatibility.md)
- Full colored command manual / 完整彩色终端手册：`vcftools-ng --help`
- [Original VCFtools 0.1.17 known issues / Original 已知问题](original-vcftools-0.1.17-known-issues.md)

## Releases / 发布记录

- [Release history / 发布历史](VERSION_HISTORY.md)
- [Per-version engineering records / 逐版本工程记录](versions/README.md)
- [Release notes / Release 说明](releases/v0.14.2.md)
- [Release, packaging, and publication workflow / 测试、打包与发布流程](release-workflow.md)

## Exactness gates and benchmarks / 一致性门禁与基准

- [Fixed benchmark and exactness workflow / 固定基准与一致性流程](benchmark-workflow.md)
- [Representative workload and claim matrix](benchmark-workload-matrix.md)
  ([中文](benchmark-workload-matrix.zh-CN.md))
- [All retained benchmark records / 全部基准档案](../benchmarks/README.md)
- [v0.14.2 complete-data portable A/B](../benchmarks/results/full-unified-v0142-ab/README.md)
  ([中文](../benchmarks/results/full-unified-v0142-ab/README.zh-CN.md))
- [v0.14.1 complete-data release gate](../benchmarks/results/final-full-v0141/README.md)
  ([中文](../benchmarks/results/final-full-v0141/README.zh-CN.md))
- [v0.13.1 locked 230k W03-W10 matrix](../benchmarks/results/workload-matrix-230k-v0130/RESULTS.md)
  ([中文](../benchmarks/results/workload-matrix-230k-v0130/RESULTS.zh-CN.md))
- [v0.13.0 input/output/storage matrix](../benchmarks/results/v0130-input-output-storage/README.md)

Every retained result directory includes, as applicable, input/oracle hashes,
raw TSV timing and resource data, logs, exactness status, and the resumable
driver. Original baselines are generated once and hash-locked rather than
rerun during routine optimization.

## Architecture / 架构

- [Adaptive input and index backends / 自适应输入与索引后端](architecture/adaptive-input-backends.md)
- [Immutable capability QueryPlan](architecture/query-plan.md)
- [Detailed implementation, build, and verification reference](../TECHNICAL_REFERENCE.md)
  ([中文](../TECHNICAL_REFERENCE.zh-CN.md))

## Build and verification / 构建与验证

- [Build from source](../TECHNICAL_REFERENCE.md#build-from-source)
  ([中文](../TECHNICAL_REFERENCE.zh-CN.md#从源码构建))
- [CTest and exact differential verification](../TECHNICAL_REFERENCE.md#verify)
  ([中文](../TECHNICAL_REFERENCE.zh-CN.md#测试与复现))
- [Portable Linux x86_64 packaging evidence](../benchmarks/results/packaging-v0141/README.md)

Performance numbers are workload-specific. A speedup from one input format,
output format, storage device, sample density, or thread count must not be
presented as a guarantee for another workload. Scientific exactness always
requires byte/hash equality; the project's 5% tie band applies only to timing.
