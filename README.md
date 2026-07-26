# vcftools-ng

[English](README.md) | [简体中文](README.zh-CN.md)

Experimental high-performance, output-compatible successor to VCFtools 0.1.17.

**Latest release:** [v0.11.3](https://github.com/VensinMa/vcftools-ng/releases/tag/v0.11.3)

**Release candidate:** v0.12.1 — Fused Site Statistics and Scalable Exact
Recode. Its permanent
[three-scenario development gate](docs/benchmark-workflow.md) has passed; the
full seven-scenario release matrix is in progress.

· all 210 full-data adaptive-counts candidate runs byte-identical and faster
than VCFtools 0.1.17

The current implementation uses adaptive ordered input shards and a bounded,
order-preserving pipeline. Plain VCF uses aligned byte ranges; BGZF VCF plus
TBI/CSI and BCF plus CSI use independent indexed regions; inputs without a
sidecar are automatically indexed with bcftools when possible, and inputs
without a usable direct backend fall back to an HTSlib stream. Persistent
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
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.11.3/vcftools-ng-v0.11.3-linux-x86_64.tar.gz
curl -LO https://github.com/VensinMa/vcftools-ng/releases/download/v0.11.3/vcftools-ng-v0.11.3-linux-x86_64.tar.gz.sha256
sha256sum -c vcftools-ng-v0.11.3-linux-x86_64.tar.gz.sha256
tar -xzf vcftools-ng-v0.11.3-linux-x86_64.tar.gz
./vcftools-ng-v0.11.3-linux-x86_64/bin/vcftools-ng --version
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

BCF plus CSI is currently the fastest input. BGZF VCF plus TBI/CSI and plain
VCF also use parallel input adapters and continue scaling from 8 to 16
threads. A local BGZF VCF or BCF without an index now runs
`bcftools index --csi --threads N` automatically and then selects the indexed
adapter. `N` is the explicit vcftools-ng thread count or the automatically
detected scheduler/affinity count. Use `--no-auto-index` to retain the
streaming behavior, `--input-backend stream` to force it, or
`--bcftools FILE` to select the executable. Optional `--chr`, `--from-bp`,
and `--to-bp` selections are pushed into indexed shards and rechecked by the
compatibility filter.

For unfiltered `--counts`, v0.11.3 uses a lower-overhead adaptive fused
backend. The 2.3-million-record validation on the 32-CPU host measured:

| Input path | Original | 1 thread | 2 threads | 4 threads | 8 threads | 16 threads | 32 threads |
|---|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 59.72 s | 18.59 s | 10.51 s | 8.41 s | 4.45 s | 2.89 s | 1.79 s |
| BGZF VCF + automatic CSI | 60.02 s | 18.43 s | 10.42 s | 10.52 s | 15.85 s | 14.57 s | 13.45 s |
| BGZF VCF, no automatic index | 60.18 s | 18.50 s | 10.46 s | 10.44 s | 10.40 s | 10.40 s | 10.42 s |
| Plain VCF | 31.21 s | 12.69 s | 6.82 s | 4.10 s | 3.22 s | 2.56 s | 2.21 s |
| BCF + CSI | 41.56 s | 12.64 s | 12.61 s | 7.16 s | 5.04 s | 3.22 s | 2.49 s |
| BCF + automatic CSI | 41.52 s | 12.75 s | 12.74 s | 14.48 s | 8.80 s | 5.24 s | 4.04 s |
| BCF, no automatic index | 41.46 s | 12.68 s | 12.58 s | 11.25 s | 5.78 s | 3.37 s | 3.36 s |

Automatic-CSI timings include first-run index construction. With
`--no-auto-index`, BGZF is deliberately limited by one ordered compressed
stream; extra requested threads do not add overhead or change output.

The final five-repeat 11.23-million-record means were:

| Input path | Original | 1 thread | 2 threads | 4 threads | 8 threads | 16 threads | 32 threads |
|---|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 296.69 s | 97.20 s | 51.76 s | 41.67 s | 21.46 s | 12.69 s | 8.37 s |
| BGZF VCF + automatic CSI | 297.06 s | 96.92 s | 51.44 s | 49.10 s | 77.88 s | 69.39 s | 64.92 s |
| BGZF VCF, no automatic index | 293.28 s | 96.34 s | 51.09 s | 49.03 s | 48.99 s | 48.96 s | 49.03 s |
| Plain VCF | 163.95 s | 90.53 s | 47.27 s | 27.24 s | 20.98 s | 20.56 s | 20.97 s |
| BCF + CSI | 199.39 s | 60.36 s | 60.17 s | 36.46 s | 25.84 s | 14.96 s | 10.90 s |
| BCF + automatic CSI | 198.45 s | 60.54 s | 60.05 s | 71.56 s | 43.76 s | 25.07 s | 18.66 s |
| BCF, no automatic index | 198.44 s | 61.32 s | 62.20 s | 55.47 s | 28.14 s | 16.50 s | 16.02 s |

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
