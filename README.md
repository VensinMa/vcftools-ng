# vcftools-ng

[English](README.md) | [简体中文](README.zh-CN.md)

Experimental high-performance, output-compatible successor to VCFtools 0.1.17.

**Latest release:** [v0.12.2 — Workload-Adaptive Indexing and Exact Recode
Scaling](https://github.com/VensinMa/vcftools-ng/releases/tag/v0.12.2)

The full 11,230,392-record release matrix passed in four representative input
scenarios at 1/2/4/8/16/32 threads. All 120 vcftools-ng outputs (five repeats
per configuration) were byte-identical to the retained VCFtools 0.1.17
goldens. The hash-locked Original timings and goldens from v0.12.1 were reused;
Original was not rerun for v0.12.2.

The current implementation uses adaptive ordered input shards and a bounded,
order-preserving pipeline. Plain VCF streams at 1–2 threads and uses aligned
byte ranges at 4+ threads. Full-file BGZF recode streams at one thread, then
reuses or builds an index at 2+ threads. Full-file BCF recode streams even
when CSI exists; selective BGZF/BCF queries still reuse or build an index.
Persistent
compute workers and the ordered committer overlap work on up to three batches.
For eligible unfiltered site statistics, an adaptive fused path avoids
`bcf1_t` construction. It currently covers `--freq`, `--freq2`, `--counts`,
`--missing-site`, `--site-depth`, `--site-mean-depth`, and `--site-quality`.
Plain VCF workers parse aligned byte ranges directly, while indexed BGZF
workers query ordered tabix regions. Worker output is bounded and committed
in exact input order.

## Install (recommended)

The portable Linux x86_64 archive is ready to run after extraction. It
includes vcftools-ng, bcftools for automatic CSI construction, HTSlib, and
the required non-glibc runtime libraries.

```bash
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.12.2/vcftools-ng-v0.12.2-linux-x86_64.tar.gz
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.12.2/vcftools-ng-v0.12.2-linux-x86_64.tar.gz.sha256
sha256sum -c vcftools-ng-v0.12.2-linux-x86_64.tar.gz.sha256
tar -xzf vcftools-ng-v0.12.2-linux-x86_64.tar.gz
./vcftools-ng-v0.12.2-linux-x86_64/bin/vcftools-ng --version
```

The archive requires Linux x86_64 with glibc 2.17 or newer. It is built on a
CentOS 7-compatible manylinux2014 baseline and is tested in clean CentOS 7
and Ubuntu 20.04 containers. It therefore also covers newer CentOS Stream,
Rocky Linux, AlmaLinux, Ubuntu, and Debian releases on x86_64. No CMake,
compiler, Conda environment, system HTSlib, or system bcftools is required.
Keep the extracted `bin`, `lib`, and `libexec` directories together.

## Version records

Every development version has a local record of supported parameters,
compatibility fixtures, retained-record counts, byte-comparison results,
original/8-thread/16-thread benchmarks, speedups, CPU, RSS, and reproduction
commands where those measurements were recorded:

- [Version history](docs/VERSION_HISTORY.md)
- [Per-version archive and update policy](docs/versions/README.md)
- [Machine-readable benchmark table](docs/versions/benchmarks.tsv)
- [New-version record template](docs/versions/TEMPLATE.md)

## Compatibility contract

VCFtools 0.1.17 is the compatibility oracle. Development correctness and
performance are evaluated on a real 2,300,000-record subset: 100,000 records
from each of 23 chromosomes. Every supported output must pass `cmp` against
the original program. The v0.11.2 final input-backend matrix additionally ran
the complete 11,230,392-record dataset in seven input/index scenarios, with
original plus 1/2/4/8/16/32-thread vcftools-ng runs repeated five times. All
245 outputs were byte-identical. v0.11.3 then ran its unfiltered `--counts`
fast path on the standard 2,300,000-record subset in the same seven scenarios
at 1/2/4/8/16/32 threads. All 42 candidate outputs were byte-identical and
every candidate was faster than its scenario's original run. The final
11,230,392-record gate then reused the hash-validated five-run Original
oracles from v0.11.2 and ran v0.11.3 five times at every thread count:
210/210 new outputs were byte-identical and faster than Original.

v0.12.1 generated the retained Original oracles for BGZF VCF, Plain VCF, and
BCF. v0.12.2 hash-validated and reused those oracles and locked timings, then
ran the same seven-filter exact-recode workload in four representative
scenarios. All 120 outputs across five repeats and 1/2/4/8/16/32 threads
passed complete-file `cmp` and SHA validation.

### Supported parameter status

The command-line surface is divided into three evidence classes:

1. Original 0.1.17 parameters with byte-identical real-data gates and measured
   optimization;
2. Original parameters with exact output gates but limited standalone
   performance evidence (three of the four two-file diff outputs);
3. vcftools-ng-only execution/output extensions.

`--recode-vcf-gz` belongs to the third group: Original VCFtools 0.1.17 does
**not** provide this option. Conversely, `--freq2` and `--stdout` are Original
parameters and are included in the compatibility surface.

The complete parameter-by-parameter matrix, performance evidence, input
coverage, compatibility boundaries, and intentionally inherited Original bugs
are documented in
[Parameter compatibility and optimization status](docs/parameter-compatibility.md).

This is not yet a complete replacement for every VCFtools option. A parameter
or interaction must not be treated as compatible until it has a differential
test. Existing gates validate documented real workloads, not the full
Cartesian product of every parameter and every input encoding.

Known Original defects are not hidden by the compatibility claim:

- exact BCF-to-VCF mode reproduces Original's malformed structured-header,
  missing-GT, and Character/String FORMAT bytes;
- Original corrupts BCF output when BCF input is combined with genotype
  masking, so vcftools-ng rejects that `--recode-bcf` combination;
- Original PCA is undefined/misaligned with missing genotypes, so exact PCA
  requires complete genotypes (`--max-missing 1`).

These behaviors and the intentionally preserved `--non-ref-af-any` and
partially missing diff quirks are detailed in the parameter matrix. A future
standards-correct mode must be explicit and must not silently alter
`--compat exact`.

## Build from source

Source builds are intended for developers or platforms not covered by the
portable archive. They require CMake, a C++20 compiler, HTSlib, LAPACK, zlib,
and POSIX threads.

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DHTSLIB_ROOT=/path/to/htslib
cmake --build build -j
```

## Run

```bash
./build/vcftools-ng \
  --bcf tests/fixtures/osmanthus412.23chr_100k.bcf \
  --threads 16 \
  --min-alleles 2 --max-alleles 2 --remove-indels \
  --minQ 40 --minGQ 20 --minDP 5 --maxDP 30 \
  --min-meanDP 10 --max-missing 0.9 --maf 0.1 \
  --recode --recode-INFO-all \
  --out results/subset
```

BCF backend performance is workload-dependent: indexed regions are valuable
for selective queries, while streaming is faster for a full-file
filter/recode. v0.12.2 therefore treats indexing only as an adaptive
acceleration technique. Plain VCF is never indexed; one-thread BGZF recode
streams; multi-thread BGZF recode reuses or builds an index; full-file BCF
recode streams even when CSI exists; and selective BGZF/BCF queries reuse or
build an index. Compact full-scan statistics reuse an existing index from
four threads but do not build a one-use index. Existing sidecars are never
overwritten. `--input-backend stream|indexed` remains an advanced override,
and `--bcftools FILE` selects the executable used when adaptive CSI
construction is profitable.

## v0.12.2 full-data five-repeat performance

The release workload applies seven real-project filters and writes a complete
VCF with all INFO fields. Values are mean wall-clock seconds from five
strictly serial vcftools-ng repeats on the 32-CPU host. Original values are
the hash-locked single-run v0.12.1 baselines and were not rerun.

| Input path | Original | 1 thread | 2 threads | 4 threads | 8 threads | 16 threads | 32 threads |
|---|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 2267.88 | 309.28 | 204.91 | 112.35 | 73.93 | 52.05 | 42.32 |
| BGZF VCF + automatic CSI | 2267.88 | 308.21 | 281.03 | 167.50 | 128.60 | 107.16 | 99.59 |
| Plain VCF | 2092.91 | 285.45 | 285.22 | 89.31 | 63.71 | 47.84 | 50.61 |
| BCF adaptive stream | 1943.47 | 325.66 | 164.14 | 110.57 | 57.26 | 41.30 | 40.20 |

Speedup over VCFtools 0.1.17:

| Input path | 1 thread | 2 threads | 4 threads | 8 threads | 16 threads | 32 threads |
|---|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 7.33× | 11.07× | 20.19× | 30.68× | 43.57× | 53.59× |
| BGZF VCF + automatic CSI | 7.36× | 8.07× | 13.54× | 17.63× | 21.16× | 22.77× |
| Plain VCF | 7.33× | 7.34× | 23.44× | 32.85× | 43.75× | 41.35× |
| BCF adaptive stream | 5.97× | 11.84× | 17.58× | 33.94× | 47.05× | 48.35× |

All 120 outputs were byte-identical. Automatic-CSI rows include fresh CSI
construction in every repeat. The observed mean speedup range was
5.97×–53.59×. The BCF adaptive path avoids indexed full-file traversal at
every tested thread count; selective region queries still use CSI.

Plain VCF cannot use CSI/TBI because those formats store BGZF virtual
offsets. It remains on the parallel aligned-byte-range adapter and never
invokes automatic indexing. Existing `.csi` and `.tbi` files are loaded and
validated independently. A valid sidecar is always preserved; if only
unusable or stale sidecars exist, automatic mode warns and falls back instead
of overwriting user files.

The input architecture, remaining ordered-BGZF work, and 1–32-thread
development gate are documented here:

- [Adaptive input backend design](docs/architecture/adaptive-input-backends.md)
- [Input backend benchmark driver](benchmarks/run-input-backend-matrix.sh)
- [Final full-dataset benchmark driver](benchmarks/run-final-full-matrix.sh)
- [Final full-dataset summary](benchmarks/results/final-full-v0112/summary.tsv)
- [v0.11.3 subset summary](benchmarks/results/adaptive-v0113-subset-final/summary.tsv)
- [v0.11.3 full-data summary](benchmarks/results/final-full-v0113/summary.tsv)
- [v0.12.2 full-data matrix](benchmarks/results/final-full-v0122/README.md)
- [v0.12.2 technical record](docs/versions/v0.12.2.md)

## Verify

```bash
ctest --test-dir build --output-on-failure
```

The differential test performs unfiltered and combined-filter comparisons on
all 2.3 million subset records, then verifies single-thread and multithread
determinism.

Run the standard 8/16-thread statistics benchmark with:

```bash
./benchmarks/run-subset.sh
```

Run the filtered VCF recode benchmark with:

```bash
./benchmarks/run-recode-subset.sh
```

Run the sample-selection plus six-output benchmark with:

```bash
./benchmarks/run-sample-subset.sh
```

Run the MAC/HWE filtering benchmark with:

```bash
./benchmarks/run-genetics-filter-subset.sh
```

Run the BED interval-filter benchmark with:

```bash
./benchmarks/run-bed-subset.sh
```

Run the ordered thinning benchmark with:

```bash
./benchmarks/run-thin-subset.sh
```

Run the v0.10 individual statistics benchmark with:

```bash
./benchmarks/run-v010-statistics.sh
```

Run the v0.10 population/LD/PCA/conversion/diff benchmark with:

```bash
./benchmarks/run-v010-advanced.sh
```

Run the non-reference AF/AC benchmark with:

```bash
./benchmarks/run-non-ref-subset.sh
```

## Benchmark dataset

See [data/README.md](data/README.md) for provenance, formats, validation, and
conversion timings. The 2.3-million-record compatibility fixture is documented
in [tests/fixtures/README.md](tests/fixtures/README.md).
