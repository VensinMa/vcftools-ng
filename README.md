# vcftools-ng

[简体中文](README.zh-CN.md) · [Detailed technical reference](TECHNICAL_REFERENCE.md) · [中文技术文档](TECHNICAL_REFERENCE.zh-CN.md)

A high-performance, output-compatible successor to VCFtools 0.1.17 for common
VCF filtering, statistics, population analyses, and recoding workloads.

**Recommended release:** [v0.14.3 — Strict CPU Budget Recovery](https://github.com/VensinMa/vcftools-ng/releases/tag/v0.14.3)

vcftools-ng keeps VCFtools-style commands and exact scientific output where
compatibility is claimed, while using HTSlib, workload-specific parsers,
bounded parallel pipelines, and adaptive input backends. Performance depends
on the command, input format, storage, sample count, and thread count; the
tables below report representative measured workloads rather than a universal
speedup claim.

## What gets faster?

The following parameter families exist in Original VCFtools 0.1.17 and have
real-data byte-comparison gates plus measured performance evidence.

| Workload | Main compatible parameters | Representative measured result |
|---|---|---:|
| Site statistics in one scan | `--freq`, `--freq2`, `--counts`, `--missing-site`, `--site-depth`, `--site-mean-depth`, `--site-quality` | Six outputs: **38.68x / 46.11x** at 8/16 threads |
| Individual and HWE statistics | `--depth`, `--missing-indv`, `--het`, `--hardy`, `--site-quality` | Individual outputs: **13.02x–21.17x** at 8/16 threads; five outputs together: **39.33x** at 16 threads |
| Position and sample selection | `--positions`, `--exclude-positions`, `--keep`, `--remove`, `--indv`, `--remove-indv` | **8.17x–27.93x** at 32 threads on the locked 230k matrix |
| Diversity and FST | `--site-pi`, `--window-pi`, `--window-pi-step`, `--TajimaD`, `--weir-fst-pop`, `--fst-window-size`, `--fst-window-step` | **20.15x–29.21x** at 32 threads for representative 230k window/site workloads |
| LD and PCA | `--geno-r2` and LD controls; `--pca`, `--pca-no-norm` | Exact optimized kernels; gains are workload-dependent and PCA gains can be modest |
| Filtering and recode | Common site, genotype, FILTER/INFO, interval, and sample filters with `--recode` | v0.14.2 full-data exact gate: **10.58x–80.28x** at 1–32 threads for indexed BGZF input and equivalent plain-VCF output |
| Conversion and comparison | `--recode-bcf`; `--diff*` outputs | BCF conversion: **18.89x** at 16 threads; site discordance: **15.14x** |

The population/selection figures are three-run medians. The v0.14.2
complete-data recode figures are single-repeat release evidence over
11,230,392 records; the 81-row v0.13.0/v0.14.1/v0.14.2 portable A/B passed
complete-file exactness gates throughout.
See [benchmark scope and methodology](docs/benchmark-workload-matrix.md) and
the [v0.14.2 complete-data A/B evidence](benchmarks/results/full-unified-v0142-ab/README.md)
before comparing results across different workloads.

### Complete-data filtering and plain-VCF recode

Speedup over retained VCFtools 0.1.17 goldens for the common filtering/recode
workload (11,230,392 records):

| Input | 1 thread | 8 threads | 16 threads | 28 threads | 32 threads |
|---|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 10.58x | 55.69x | 71.59x | 78.64x | 80.28x |
| Plain VCF | 9.77x | 64.38x | 65.30x | 63.96x | 65.98x |
| BCF adaptive stream | 3.52x | 19.79x | 37.16x | 43.90x | 38.54x |

These v0.14.2 rows use `--recode-vcf` so the output container is equivalent
to Original's uncompressed VCF. Original was not rerun: the denominator is
the retained hash-validated v0.12.1 baseline for the identical input and
filter command. Default BGZF output is much smaller but has a different
compression cost; use the dedicated
[SSD/HDD output tables](benchmarks/results/v0130-input-output-storage/README.md)
for that comparison. Differences within 5% are treated as performance ties,
never as permission for scientific-output differences.

## Install

The portable Linux x86_64 archive runs after extraction. It targets glibc
2.17+ and has been checked on clean CentOS 7 and Ubuntu 20.04 containers.
`bin` contains only `vcftools-ng`; its private runtime tools remain under
`libexec`.

```bash
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.14.3/vcftools-ng-v0.14.3-linux-x86_64.tar.gz
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.14.3/vcftools-ng-v0.14.3-linux-x86_64.tar.gz.sha256
sha256sum -c vcftools-ng-v0.14.3-linux-x86_64.tar.gz.sha256
tar -xzf vcftools-ng-v0.14.3-linux-x86_64.tar.gz
./vcftools-ng-v0.14.3-linux-x86_64/bin/vcftools-ng --help
```

No compiler, CMake, Conda environment, system HTSlib, or system bcftools is
required. Keep the extracted `bin`, `lib`, and `libexec` directories together.

## Quick start

### Filter and write compressed VCF (recommended)

```bash
vcftools-ng --gzvcf input.vcf.gz --threads 24 \
  --min-alleles 2 --max-alleles 2 --remove-indels \
  --minQ 40 --minGQ 20 --minDP 5 --maxDP 30 \
  --min-meanDP 10 --max-missing 0.9 --maf 0.1 \
  --recode --recode-INFO-all --out filtered
```

This writes `filtered.recode.vcf.gz` and `filtered.log`. The compressed output
is BGZF VCF; vcftools-ng does not automatically index a newly written output:

```bash
bcftools index --tbi --threads 24 filtered.recode.vcf.gz
```

### Produce several statistics in one scan

```bash
vcftools-ng --gzvcf input.vcf.gz --threads 16 \
  --freq --counts --missing-site --site-depth --site-mean-depth \
  --site-quality --out cohort
```

Combining compatible outputs is strongly recommended: Original normally scans
the input once per command, while vcftools-ng can share parsing and genotype
decoding.

### Individual statistics

```bash
vcftools-ng --gzvcf input.vcf.gz --threads 16 \
  --depth --missing-indv --het --hardy --out cohort
```

### Window pi, Tajima's D, and FST

```bash
vcftools-ng --gzvcf input.vcf.gz --threads 16 \
  --window-pi 100000 --window-pi-step 10000 --TajimaD 100000 \
  --weir-fst-pop population1.txt --weir-fst-pop population2.txt \
  --fst-window-size 100000 --fst-window-step 10000 --out diversity
```

Window analyses require nondecreasing positions and contiguous chromosome
segments. Use one sample ID per line in each population file.

## Important differences from Original VCFtools

| Behavior | vcftools-ng default | How to override |
|---|---|---|
| `--recode` file output | Deterministic BGZF VCF: `PREFIX.recode.vcf.gz` | Use `--recode-vcf` for uncompressed `PREFIX.recode.vcf` |
| Input backend/index | Adaptive by format, workload, storage, and effective threads | Advanced diagnostics: `--input-backend stream\|plain\|indexed` |
| Threads | Detect scheduler, affinity, cgroup, and hardware limits; automatic maximum 128 | `--threads N` / `-t N` restricts the complete process tree to at most N runnable CPUs; I/O-waiting workers may overlap without increasing CPU capacity |
| Run log | `PREFIX.log` is written and terminal diagnostics remain enabled | `--log-file FILE` or `--no-log-file` |
| Multiple outputs | Compatible analyses share one scan | List all required output parameters in one command |
| Failed output | Scientific files are staged and published transactionally | No override; pre-existing destinations are preserved on failure |

`--recode-vcf-gz` is an explicit alias for compressed VCF output and does not
exist in Original 0.1.17. `--recode-vcf` is also a vcftools-ng extension. With
`--recode --stdout`, VCF remains uncompressed on stdout and diagnostics never
enter the data stream.

Adaptive indexing never overwrites an existing valid TBI/CSI sidecar. For a
full-file BGZF recode, indexed work is usually selected at two or more effective
threads; a profitable missing input index may be built. Full-file BCF normally
uses the faster streaming path even if CSI exists, while coordinate-restricted
BCF queries can use the index. Leave `--input-backend auto` unchanged for normal
use.

## Optimized compatible parameters

The compact list below groups Original parameter names with exact-output and
performance evidence; it is not a claim that every Cartesian combination has
been tested.

- Statistics: `--freq`, `--freq2`, `--counts`, `--missing-site`,
  `--site-depth`, `--site-mean-depth`, `--site-quality`, `--depth`,
  `--missing-indv`, `--het`, `--hardy`, pi, Tajima, FST, LD, and PCA options.
- Site filters: allele count, SNP/indel, QUAL, mean depth, missingness,
  MAF/MAC/HWE, non-reference frequency/count, chromosome/position/BED,
  FILTER/INFO Flag, and thinning options.
- Genotype/sample filters: `--minGQ`, `--minDP`, `--maxDP`, genotype FT,
  `--keep`, `--remove`, `--indv`, and `--remove-indv`.
- Output/comparison: `--recode`, `--recode-bcf`, `--recode-INFO-all`,
  `--stdout`, and the supported `--diff*` family.

For the exact names, evidence, unsupported combinations, polyploid boundary,
and inherited Original quirks, use the
[parameter compatibility matrix](docs/parameter-compatibility.md).

## vcftools-ng-only options worth knowing

| Option | Default | Purpose |
|---|---|---|
| `--threads N`, `-t N` | Automatic, capped at 128 | Strict shared CPU budget for the complete pipeline, including BGZF output compression |
| `--input FILE` | — | Auto-detect VCF, BGZF VCF, or BCF |
| `--recode-vcf-gz` | Off | Explicit deterministic BGZF VCF output |
| `--recode-vcf` | Off | Explicit uncompressed VCF output |
| `--log-file FILE` | `PREFIX.log` | Select run-log path |
| `--no-log-file` | Off | Disable only the log file, not stderr |
| `--input-backend ...` | `auto` | Expert backend override/diagnostic |
| `--corrected-depth-arithmetic` | Off | Checked 64-bit site-depth arithmetic instead of Original's 32-bit-wrap behavior |

Run `vcftools-ng --help` for the complete colored terminal manual, combinations,
output suffixes, defaults, and examples.

## Recommended settings

- Keep `--input-backend auto` unless diagnosing a backend.
- Omit `--threads` when scheduler/cgroup limits are correct; otherwise set it
  to the CPU allocation actually granted to the job, not the server total.
- `--threads N` is a process-tree CPU ceiling, not a pthread-count promise.
  On Linux the process and automatic bcftools children are restricted to at
  most N allowed CPUs. Input/output workers may overlap while blocked on I/O
  or ordered queues, but cannot execute on more than N cores. Explicit
  values may exceed 128 on large servers but are intersected with the actual
  scheduler/cgroup/affinity allocation; only automatic detection is capped at
  128. Resource-planning invariants are tested through 65,536 logical threads.
- Prefer `--recode` (BGZF) to avoid very large plain-VCF writes, especially on
  HDD storage. Use `--recode-vcf` only when another tool requires plain VCF.
- Add `--recode-INFO-all` when all input INFO annotations must be retained.
- Keep the default log for reproducibility and request compatible statistics
  together in one invocation.

## Documentation and evidence

Scientific compatibility means complete output bytes match retained Original
VCFtools 0.1.17 goldens for the documented workload. Compressed VCF is checked
after decompression because Original has no equivalent compressed-output
parameter. Unsupported shapes fall back to the general exact pipeline or are
rejected explicitly; they do not silently enter an incompatible fast kernel.

| Interested in | Document |
|---|---|
| Everything previously on the long homepage | [Detailed technical reference](TECHNICAL_REFERENCE.md) |
| Exact parameter names, optimized status, and boundaries | [Parameter compatibility matrix](docs/parameter-compatibility.md) |
| Release history and what changed in each version | [Release history](docs/VERSION_HISTORY.md) and [per-version records](docs/versions/README.md) |
| Byte/hash gates and the fixed test policy | [Benchmark and exactness workflow](docs/benchmark-workflow.md) |
| Which benchmark supports which performance claim | [Representative workload matrix](docs/benchmark-workload-matrix.md) |
| Complete multi-version and raw benchmark results | [Benchmark archive](benchmarks/README.md) |
| v0.14.2 full-data A/B and exact-output gates | [v0.14.2 evidence](benchmarks/results/full-unified-v0142-ab/README.md) |
| v0.14.3 strict CPU-budget design and fair 230k A/B | [v0.14.3 record](docs/versions/v0.14.3.md) |
| Build from source and verification commands | [Build and verify](TECHNICAL_REFERENCE.md#build-from-source) |
| Input/index scheduling and capability planning | [Adaptive input backends](docs/architecture/adaptive-input-backends.md) and [QueryPlan](docs/architecture/query-plan.md) |
| How a release is tested, packaged, and published | [Release workflow](docs/release-workflow.md) |
| Original behaviors that are retained, corrected, or rejected | [Original VCFtools known-issue ledger](docs/original-vcftools-0.1.17-known-issues.md) |

The [documentation index](docs/README.md) provides the same routes in one
place, including links to raw TSV data and reproduction scripts.

License: LGPL-3.0-or-later. See [LICENSE](LICENSE).
