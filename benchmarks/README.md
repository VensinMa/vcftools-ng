# Subset benchmark

The permanent optimization workflow is documented in
[`docs/benchmark-workflow.md`](../docs/benchmark-workflow.md). Daily
development uses only BGZF+TBI, Plain VCF, and the BCF adaptive streaming
full-scan path through `run-development-gate.sh`, reusing hash-locked Original
outputs and timings. The BCF row uses the default `auto` policy and verifies
that it selects the streaming backend.
The complete four-scenario matrix is reserved for an explicitly authorized
release candidate: BGZF+TBI, BGZF+automatic-CSI, Plain VCF, and adaptive BCF.

## Policy

Until the implementation reaches a mature compatibility milestone, benchmarks
use only the real 2,300,000-record subset. It contains exactly 100,000 records
from each of chromosomes `chr1` through `chr23` and all 412 samples. The full
11,230,392-record dataset is not part of the current benchmark loop.

All reported `vcftools-ng` outputs were checked byte for byte with `cmp`
against files produced by VCFtools 0.1.17.

## v0.11 adaptive input backends

The counts workload isolates input, parsing, shared GT decoding, and ordered
output on the same 2,300,000 records. Every row was compared byte for byte
with an original VCFtools 0.1.17 output.

| Input | Original | 8 threads | 16 threads | Speedup 8/16 | CPU 8/16 | RSS KB 8/16 | Exact |
|---|---:|---:|---:|---:|---:|---:|---|
| BGZF VCF + TBI | 68.23 | 16.06 | 11.25 | 4.25× / 6.06× | 535% / 894% | 409,476 / 541,000 | PASS |
| BGZF VCF, no index | 68.22 | 53.01 | 53.12 | 1.29× / 1.28× | 140% / 145% | 24,768 / 25,020 | PASS |
| Plain VCF | 36.15 | 11.11 | 8.23 | 3.25× / 4.39× | 685% / 1156% | 483,716 / 795,148 | PASS |
| BCF + CSI | 43.89 | 5.15 | 3.32 | 8.52× / 13.22× | 458% / 832% | 268,756 / 514,536 | PASS |
| BCF, no index | 43.78 | 5.80 | 3.52 | 7.55× / 12.44× | 430% / 846% | 37,392 / 47,356 | PASS |

The indexed BGZF and Plain adapters retain substantial 8-to-16-thread
scaling. Unindexed BCF still benefits from HTSlib's parallel BGZF
decompression. Unindexed BGZF VCF remains flat because one thread performs
text record parsing after decompression.

Raw evidence, backend logs, SHA256 comparisons, environment metadata, and
summary values are stored in
[`benchmarks/results/input-vnext-counts`](results/input-vnext-counts).
The complete regression log is
[`benchmarks/results/v011-ctest.txt`](results/v011-ctest.txt).

## v0.11.1 automatic CSI construction

When a local BGZF VCF/BCF lacks TBI/CSI, the default behavior now builds CSI
with the effective vcftools-ng thread count and then uses indexed regions.
The real-subset first-run measurement includes both indexing and counts:

| Engine | Original | 8 threads | 16 threads | Speedup 8/16 | Exact |
|---|---:|---:|---:|---:|---|
| BGZF VCF, starting without index | 68.22 | 27.84 | 22.69 | 2.45× / 3.01× | PASS |

Raw timings, logs, and hashes are stored in
[`benchmarks/results/v0111-auto-index`](results/v0111-auto-index). CSI-only
construction measurements are in
[`benchmarks/results/csi-index-build`](results/csi-index-build), and the
cumulative regression log is
[`benchmarks/results/v0111-ctest.txt`](results/v0111-ctest.txt).

## v0.11.2 validated and protected indices

v0.11.2 explicitly validates CSI and TBI rather than treating file existence
as validity. It also confirms that Plain VCF never enters the index-building
path.

| Case | Original | 8 threads | 16 threads | Speedup 8/16 | Exact |
|---|---:|---:|---:|---:|---|
| BGZF, build and validate CSI, then counts | 68.22 | 27.82 | 23.17 | 2.45× / 2.94× | PASS |
| Plain VCF, no index attempt | 36.15 | 11.18 | 8.29 | 3.23× / 4.36× | PASS |

Evidence is stored in
[`benchmarks/results/v0112-index-validation`](results/v0112-index-validation);
the cumulative test log is
[`benchmarks/results/v0112-ctest.txt`](results/v0112-ctest.txt).

## Input-backend matrix

The post-v0.10 input-engine work treats file encodings as separate execution
paths. Its design and correctness constraints are recorded in
[`docs/architecture/adaptive-input-backends.md`](../docs/architecture/adaptive-input-backends.md).

The benchmark driver covers:

- BGZF VCF with TBI;
- the same BGZF bytes without an index sidecar;
- plain VCF;
- BCF with CSI;
- the same BCF bytes without an index sidecar.

It always runs original VCFtools 0.1.17 and byte-compares each vcftools-ng
result. Its deliberately unindexed historical cases pass `--no-auto-index`;
all ordinary vcftools-ng runs retain automatic CSI construction. Development
runs default to 8 and 16 threads:

```bash
NG=./build/vcftools-ng \
ORIGINAL=/path/to/vcftools-0.1.17/bin/vcftools \
SUBSET_BGZF=tests/fixtures/osmanthus412.23chr_100k.vcf.gz \
SUBSET_PLAIN=/path/to/osmanthus412.23chr_100k.vcf \
SUBSET_BCF=tests/fixtures/osmanthus412.23chr_100k.bcf \
./benchmarks/run-input-backend-matrix.sh
```

The final scalability gate can extend the same command without changing its
correctness oracle:

```bash
THREAD_LIST="1 2 4 8 16 32 64 128 256" REPEATS=3 \
CACHE_STATE=warm \
./benchmarks/run-input-backend-matrix.sh benchmarks/results/input-final
```

Thread points above the scheduler allocation are omitted by the caller. Final
reports must also state storage type and whether each repetition used a cold
or warm page cache; the script records the available host and block-device
metadata plus the caller-supplied `CACHE_STATE`, but does not flush the
operating-system cache. `EXPECTED_RECORDS` defaults to 2,300,000 for the real
subset and is stored in the manifest.

## v0.10 individual and HWE statistics

The five new outputs are `--depth`, `--missing-indv`, `--het`, `--hardy`,
and `--site-quality`. Original VCFtools requires five independent scans;
vcftools-ng can produce all five from one shared GT/DP decode.

| Workload | Original | 8 threads | 16 threads | Speedup 8/16 | Exact |
|---|---:|---:|---:|---:|---|
| Five separate outputs | 335.08 s | 23.01 s | 19.72 s | 14.56× / 16.99× | PASS |
| Five outputs, one scan | 335.08 s | 8.63 s | 8.52 s | 38.83× / 39.33× | PASS |

The combined 8-thread run used 428% CPU and 68,344 KB peak RSS; the
16-thread run used 457% CPU and 63,576 KB peak RSS. Reproduce the original,
8-thread, and 16-thread measurements and complete-file comparisons with:

```bash
./benchmarks/run-v010-statistics.sh
```

## First vertical-slice results

One invocation produced all five outputs: allele frequencies, allele counts,
site missingness, site depth, and mean site depth.

| Input | Threads | Wall time | CPU | Peak RSS | Exact result |
|---|---:|---:|---:|---:|---|
| BCF | 8 | 7.55 s | 398% | 64.3 MiB | PASS |
| BCF | 16 | 6.56 s | 417% | 64.6 MiB | PASS |
| BGZF VCF | 8 | 56.29 s | 140% | 64.8 MiB | PASS |
| BGZF VCF | 16 | 56.23 s | 138% | 65.1 MiB | PASS |

The original VCFtools 0.1.17 needs one scan per output. Its five sequential
subset runs took 67.96, 67.84, 67.47, 74.25, and 75.09 seconds, or 352.61
seconds in total. Against that five-output workload, the 16-thread BCF run is
53.8 times faster.

The BCF scaling curve flattens after eight threads because input decoding and
ordered output formatting become dominant. The nearly flat BGZF VCF results
show that text record parsing is the principal bottleneck for that input
format; more worker threads alone cannot remove it.

## Combined-filter differential case

The real subset was also checked with:

```text
--min-alleles 2 --max-alleles 2 --remove-indels --minQ 40
--minGQ 20 --minDP 5 --maxDP 30 --min-meanDP 10
--max-missing 0.9 --maf 0.1
```

Both implementations retained 129,855 of 2,300,000 records. All five output
files were byte-identical. The 16-thread BCF implementation generated all five
files in one 4.68-second scan.

## Filtered VCF recode

The second vertical slice adds `--recode` and `--recode-INFO-all`. Genotypes
removed by `--minGQ`, `--minDP`, or `--maxDP` have only their GT changed to
missing, matching VCFtools 0.1.17; the remaining FORMAT values are retained.

The INFO-all golden was generated by VCFtools 0.1.17 from the BGZF VCF subset
in 156.14 seconds. The output contains 129,855 records and is 1.4 GB.

| Input | Threads | Wall time | CPU | Peak RSS | Exact result |
|---|---:|---:|---:|---:|---|
| BCF | 8 | 5.42 s | 526% | 70.4 MiB | PASS |
| BCF | 16 | 4.58 s | 609% | 71.4 MiB | PASS |

The 8-thread run is 28.8 times faster than the original; the 16-thread run is
34.1 times faster. The default `--recode` behavior, which replaces INFO with
`.`, also passed the full 2.3-million-record byte comparison.

## Reproduce

```bash
./benchmarks/run-subset.sh
```

The statistics script runs 8 and 16 threads, compares every generated file with the
VCFtools 0.1.17 golden output, records timing under `benchmarks/results/`, and
removes the duplicate generated files after validation.

The recode script follows the same 8/16-thread policy:

```bash
./benchmarks/run-recode-subset.sh
```

## Site-selection filters

The third vertical slice adds chromosome, closed-range, and position-list
selection. These checks run before GT/DP/GQ decoding.

| Selection | NG threads | Original | NG | Exact result |
|---|---:|---:|---:|---|
| positions + exclude-positions + not-chr | 8 | 39.80 s | 3.39 s | PASS |
| indexed chr + from-bp + to-bp | 8 | 48.88 s | 0.24 s | PASS |
| indexed chr + from-bp + to-bp | 16 | 48.88 s | 0.19 s | PASS |

The first case retained 16 records after scanning the complete input. The
indexed cases read and retained the 29,816 records in `chr7:1-2000000`
directly from the CSI. They are about 204 and 257 times faster than the
original 48.88-second full scan. Every comparison covered the entire recoded
VCF byte for byte.

### Rejected scheduler experiment

A persistent worker-pool prototype was measured against the existing
batch-local scheduler with all other variables held constant. It increased
the 8/16-thread recode times from 5.46/4.58 seconds to 6.12/5.56 seconds
because the per-batch condition-variable barrier roughly doubled system time.
The prototype was removed. A future pipeline implementation must overlap
input, compute, and output instead of only changing worker lifetime.

## Sample selection

The fourth vertical slice adds file-based and inline sample inclusion and
exclusion. The compatibility workload retains 17 of 412 samples, applies the
combined site/genotype filters, and produces all five statistics plus an
INFO-all recoded VCF in one scan.

| Threads | Wall time | CPU | Peak RSS | Six exact outputs |
|---:|---:|---:|---:|---|
| 8 | 4.08 s | 490% | 63.6 MiB | PASS |
| 16 | 3.73 s | 524% | 63.4 MiB | PASS |

VCFtools 0.1.17 requires six scans for these outputs. The measured individual
runs total 623.49 seconds sequentially, versus one 4.08/3.73-second scan:
approximately 153 and 167 times faster. The recoded VCF header, retained
sample columns, allele counts, missingness, depth statistics, GT masking, MAF,
and call-rate decisions are all covered by byte comparisons.

Reproduce with:

```bash
./benchmarks/run-sample-subset.sh
```

## Minor-allele count and exact HWE filtering

The fifth vertical slice adds `--mac`, `--max-mac`, and `--hwe`. MAC is the
minimum count across all alleles, including REF, exactly as in VCFtools
0.1.17. The HWE test ports the original Wigginton exact-test recurrence and
uses only fully called diploid genotypes remaining after genotype and sample
filters.

The compatibility workload retains 17 samples and 57,595 of 2,300,000 sites
with `--mac 4 --max-mac 12 --hwe 0.001`. The original BGZF VCF run took
104.85 seconds. Reproduce the 8/16-thread BCF measurements with:

```bash
./benchmarks/run-genetics-filter-subset.sh
```

| Threads | Wall time | CPU | Peak RSS | Exact result |
|---:|---:|---:|---:|---|
| 8 | 4.01 s | 498% | 63.8 MiB | PASS |
| 16 | 3.64 s | 538% | 63.3 MiB | PASS |

This is approximately 26.1 and 28.8 times faster than the original run.

## BED interval filtering

The sixth vertical slice adds `--bed` and `--exclude-bed`. BED records are
indexed by chromosome; each lookup uses a binary search and occurs before
genotype/depth decoding. The exact VCFtools overlap rule is retained,
including the longest REF/ALT allele span.

The include workload retains 79,437 sites. The original BGZF VCF scan took
40.95 seconds, while the first 8-thread BCF calibration took 3.50 seconds.
The complementary `--exclude-bed` path retained 2,220,563 sites and also
matched all 2,220,564 output lines byte for byte.

Reproduce the standard 8/16-thread include measurements with:

```bash
./benchmarks/run-bed-subset.sh
```

| Threads | Wall time | CPU | Peak RSS | Exact result |
|---:|---:|---:|---:|---|
| 8 | 3.50 s | 460% | 62.5 MiB | PASS |
| 16 | 3.39 s | 489% | 63.6 MiB | PASS |

These runs are approximately 11.7 and 12.1 times faster than the original
include scan.

## Ordered site thinning

The seventh vertical slice adds `--thin`. VCFtools applies thinning after all
other site and genotype filters, so the result depends on the previous site
that survived those filters. `vcftools-ng` computes independent filters in
parallel and applies this final state transition in input order when each
batch is committed.

The real-subset workload combines the BED selector with `--thin 10000`,
retains 433 sites, and matches all 434 count-file lines. The original took
39.94 seconds; the first 8-thread BCF calibration took 3.51 seconds.

Reproduce the standard 8/16-thread measurements with:

```bash
./benchmarks/run-thin-subset.sh
```

| Threads | Wall time | CPU | Peak RSS | Exact result |
|---:|---:|---:|---:|---|
| 8 | 3.49 s | 461% | 63.3 MiB | PASS |
| 16 | 3.36 s | 484% | 62.5 MiB | PASS |

These runs are approximately 11.4 and 11.9 times faster than the original.

## Non-reference AF and AC filters

The eighth vertical slice adds the all-ALT and any-ALT lower/upper bounds for
non-reference allele frequency and count. The real-subset workload combines
all eight options with sample and BED selection, retains 22,168 sites, and
matches all 22,169 count-file lines. The original took 41.02 seconds; the
first 8-thread BCF calibration took 3.47 seconds.

VCFtools 0.1.17 contains a compatibility quirk: `--non-ref-af-any` used
without a non-`any` AF limit does not filter. The new engine intentionally
retains that behavior in exact mode and tests it against the unfiltered
2.3-million-site golden.

Reproduce the standard 8/16-thread measurements with:

```bash
./benchmarks/run-non-ref-subset.sh
```

| Threads | Wall time | CPU | Peak RSS | Exact result |
|---:|---:|---:|---:|---|
| 8 | 3.47 s | 462% | 62.3 MiB | PASS |
| 16 | 3.36 s | 484% | 63.0 MiB | PASS |

These runs are approximately 11.8 and 12.2 times faster than the original.

## FILTER/INFO/FT compatibility and overlapped pipeline

The ninth vertical slice adds site FILTER, INFO Flag, and genotype FT
filtering. A 23,000-record fixture derived from the real 23-chromosome subset
covers missing, PASS, single, multiple, and order-sensitive flags. Site
filtering, genotype counts, missingness, GT masking, FT retention, and recode
files are byte-identical to VCFtools 0.1.17.

The scheduler is now a bounded reader → compute → ordered-commit pipeline:

- at most three batches are in flight;
- batches are split into 64–256-record slices;
- persistent workers reuse decode scratch space;
- thinning and output remain single-owner and input ordered;
- reader, worker, or writer failure cancels every lane without deadlock.

The default batch size changed from 8,192 to 2,048 after real-subset tuning.
For the 1.4 GB filtered recode workload, this reduced peak RSS from 191 MiB
to about 55 MiB and improved wall time from 3.73/3.79 seconds to 3.36/3.41
seconds at 8/16 threads.

Against the removed batch-barrier scheduler:

| Workload | Threads | Barrier | Pipeline | Improvement |
|---|---:|---:|---:|---:|
| sample selection + six outputs | 8 | 4.09 s | 3.11 s | 24.0% |
| sample selection + six outputs | 16 | 3.76 s | 3.13 s | 16.8% |
| 1.4 GB filtered recode | 8 | 5.72 s | 3.37 s | 41.1% |
| 1.4 GB filtered recode | 16 | 4.67 s | 3.37 s | 27.8% |

## v0.10 individual, population, LD/PCA, conversion, and diff

v0.10 keeps the 2,300,000-record real-subset gate and adds individual
statistics, ordered population/window statistics, genotype LD, PCA,
deterministic parallel BGZF conversion, and two-file comparison.

The five new individual outputs require 335.08 seconds as five original
scans. A single shared vcftools-ng scan takes 8.63/8.52 seconds at 8/16
threads, a 38.83×/39.33× speedup. Site π, FST, LD, and PCA each remain
byte-identical and take approximately 3.0–3.5 seconds.

The full 4.3 GB BCF output is compared as a compressed byte stream:

| Implementation | Wall | CPU | Peak RSS | Speedup | Exact |
|---|---:|---:|---:|---:|---|
| Original 0.1.17 | 425.51 s | 99% | 5,980 KB | 1.00× | oracle |
| vcftools-ng, 8 threads | 38.81 s | 868% | 57,024 KB | 10.96× | PASS |
| vcftools-ng, 16 threads | 22.53 s | 1677% | 61,392 KB | 18.89× | PASS |

Two-file site discordance compares the complete 2,300,001-line output:

| Implementation | Wall | CPU | Peak RSS | Speedup | Exact |
|---|---:|---:|---:|---:|---|
| Original 0.1.17 | 136.14 s | 99% | 6,728 KB | 1.00× | oracle |
| vcftools-ng, 8 threads | 8.97 s | 393% | 8,660 KB | 15.18× | PASS |
| vcftools-ng, 16 threads | 8.99 s | 392% | 8,824 KB | 15.14× | PASS |

The diff workload is limited by reading two BCF streams, hence the flat
8-to-16-thread scaling. Raw evidence is under `benchmarks/results/v010-*`.
Reproduce the release measurements with:

```bash
./benchmarks/run-v010-statistics.sh
./benchmarks/run-v010-advanced.sh
```
