# Parameter compatibility and optimization status

This document describes the v0.13.1 command-line surface. VCFtools 0.1.17
source and output files are the compatibility oracle.

The terms below are deliberately separate:

- **Exact** means the documented workload has passed complete-file `cmp`
  against an Original 0.1.17 golden. It does not mean that every possible
  parameter combination, malformed file, ploidy, or VCF/BCF encoding has been
  tested.
- **Performance verified** means the parameter was exercised in a real-data
  workload that was faster than Original on the benchmark host. It is measured
  evidence, not a promise that every filesystem, dataset, and machine will be
  faster.
- **Extension** means the option does not exist in Original VCFtools 0.1.17.

## Original-compatible parameters with performance evidence

All options in this section exist in Original 0.1.17, have byte-comparison
coverage, and execute through a lower-overhead or parallel vcftools-ng path.
Performance is measured by workload group rather than by claiming an
independent speedup for every possible combination.

| Workload group | Original-compatible parameters | Evidence |
|---|---|---|
| Site frequency, missingness, depth, and quality | `--freq`, `--freq2`, `--counts`, `--missing-site`, `--site-depth`, `--site-mean-depth`, `--site-quality` | Real 2.3-million-record complete-file goldens. The v0.11.4 fused BGZF path passed at 8/16 threads. Six outputs in one scan measured 38.68×/46.11× against six Original scans; `--freq2` measured 7.20×/8.55×. |
| Individual and HWE statistics | `--depth`, `--missing-indv`, `--het`, `--hardy` | Real 2.3-million-record exact gates and individual/combined benchmarks in v0.10.0. |
| Diversity and population statistics | `--site-pi`, `--window-pi`, `--window-pi-step`, `--TajimaD`, `--weir-fst-pop`, `--fst-window-size`, `--fst-window-step` | v0.13.1 direct Plain VCF paths passed the 23k and 230k exact matrices at 1/4/8/16/32 threads. At 32 threads, window pi, Tajima's D, site FST, and window FST measured 20.15×–29.21× over Original in the locked 230k workloads. |
| LD and PCA | `--geno-r2`, `--ld-window`, `--ld-window-min`, `--ld-window-bp`, `--ld-window-bp-min`, `--min-r2`, `--pca`, `--pca-no-norm` | Real-data exact gates and LD/PCA benchmarks in v0.10.0. |
| VCF/BCF recode | `--recode`, `--recode-bcf`, `--recode-INFO-all`, `--stdout` | Complete VCF/BCF byte gates and conversion benchmarks. File-based `--recode` changes only the container to BGZF; decompressed VCF bytes retain the Original-compatible path, while `--recode-vcf` selects a directly comparable plain file. v0.13.1 also gates selected-sample plus seven-filter direct BGZF recode. `--recode --stdout` remains plain VCF. |
| Common numeric/site filters | `--min-alleles`, `--max-alleles`, `--remove-indels`, `--keep-only-indels`, `--minQ`, `--min-meanDP`, `--max-meanDP`, `--max-missing`, `--max-missing-count`, `--maf`, `--max-maf`, `--mac`, `--max-mac`, `--hwe` | Exact filtered statistics/recode gates and v0.2/v0.5/v0.11.4 real-data benchmarks. |
| Chromosome, position, and interval filters | `--chr`, `--not-chr`, `--from-bp`, `--to-bp`, `--positions`, `--exclude-positions`, `--bed`, `--exclude-bed`, `--thin` | Complete-file exact gates and v0.3/v0.6/v0.7 benchmarks. v0.13.1 additionally gives `--positions`/`--exclude-positions` a direct Plain VCF path; 1%/50%, sorted/shuffled, duplicate, and absent-position 230k fixtures passed byte gates and measured 8.17×–20.00× at 32 threads. |
| Non-reference filters | `--non-ref-af`, `--max-non-ref-af`, `--non-ref-af-any`, `--max-non-ref-af-any`, `--non-ref-ac`, `--max-non-ref-ac`, `--non-ref-ac-any`, `--max-non-ref-ac-any` | Complete-file exact gates and v0.8 benchmark. |
| FILTER/INFO Flag filters | `--keep-filtered`, `--remove-filtered`, `--remove-filtered-all`, `--keep-INFO`, `--remove-INFO` | Annotated real-data exact gate and v0.9 filtered-recode benchmark. `--keep-INFO`/`--remove-INFO` are Flag filters, matching this Original path. |
| Genotype filters | `--minGQ`, `--minDP`, `--maxDP`, `--remove-filtered-geno`, `--remove-filtered-geno-all` | Exact GT-masking and FT gates; filtered recode/statistics benchmarks in v0.2, v0.9, and v0.11.4. |
| Sample filters | `--keep`, `--remove`, `--indv`, `--remove-indv` | Six complete outputs passed the original real-data gate. v0.13.1 direct sample projection passed 25%/50%/100% 230k counts gates; 32-thread speedup measured 21.05×–27.93×. Selected-sample plus seven-filter BGZF recode is covered by the 23k byte gate. |

The Original input/output selectors `--vcf`, `--gzvcf`, `--bcf`, and `--out`
are also supported. Correctly matched VCF, BGZF VCF, and BCF inputs feed the
optimized adapters. This is not a claim that every parameter above has been
tested in a full Cartesian product with every input format.

## Original-compatible parameters with exactness but limited performance evidence

Two-file comparison uses Original's `--diff`, `--gzdiff`, or `--diff-bcf`
input selector with:

- `--diff-site`
- `--diff-indv`
- `--diff-site-discordance`
- `--diff-indv-discordance`

All four outputs have complete-file byte gates. A real 2.3-million-record
performance benchmark exists for `--diff-site-discordance` (15.18×/15.14× at
8/16 threads). The other three share the new merge/decode engine but have not
been benchmarked separately, so no independent speedup claim is made for
them. Diff filtering is currently limited to chromosome/position/BED and
sample selection.

## vcftools-ng-only extensions

These options are additions, not Original-compatible parameter names:

| Extension | Purpose and validation |
|---|---|
| `--threads N`, `-t N` | Sets the shared input/compute/I/O CPU budget. When omitted, vcftools-ng intersects scheduler, affinity, cgroup, and hardware limits and caps automatic selection at 128. |
| `--batch-size N` | Tunes the bounded generic pipeline batch size. |
| `--input FILE` | Auto-detected input alias for `--vcf`/`--gzvcf`/`--bcf`. |
| `--compat exact` | Explicitly selects the only currently implemented compatibility mode. |
| `--input-backend auto\|stream\|plain\|indexed` | Selects or forces an input adapter. `auto` is the default adaptive speed policy. |
| `--bcftools FILE` | Selects the bcftools executable used when the adaptive policy decides CSI construction is profitable. |
| `--recode-vcf-gz` | Writes deterministic, parallel BGZF-compressed VCF. Original 0.1.17 has no compressed-VCF recode option. Validation compares the decompressed bytes with Original `--recode`; compressed bytes are deterministic across tested thread counts, but cannot be compared with a nonexistent Original artifact. |
| `--recode-vcf` | Writes an uncompressed VCF file after v0.13.0 changed the default `--recode` file output to BGZF. `--recode --stdout` remains plain VCF for compatibility. |
| `--log-file FILE` | Overrides the default `PREFIX.log` run-log path. The log is overwritten and receives the same diagnostics as stderr. |
| `--no-log-file` | Disables the run-log file without disabling stderr diagnostics. |

Original VCFtools creates `PREFIX.log` by default. vcftools-ng now preserves
that output-prefix behavior while adding structured execution, index-policy,
resource, and exit metadata. The log is operational evidence, not a scientific
result artifact, and is not compared byte-for-byte with Original.

`--no-auto-index` was removed in v0.12.2. Index construction and index use are
now independent adaptive decisions. `--input-backend stream` and
`--input-backend indexed` remain advanced diagnostic/override controls, but
normal use should leave the default `auto` mode enabled.

Producing several supported statistics in one scan is also a vcftools-ng
extension. Original normally requires separate invocations.

`--help`, `-h`, and `--version` exist in both programs but their text/version
output is naturally program-specific.

## Deliberate boundaries and Original defects

The maintained issue ledger is
[`original-vcftools-0.1.17-known-issues.md`](original-vcftools-0.1.17-known-issues.md).
It records minimal triggers, observed Original behavior, and whether
vcftools-ng retains or rejects each case.

Exact mode reproduces observable Original output, including some Original
quirks. Where Original produces corrupt or undefined output, vcftools-ng
rejects the unsafe combination instead of presenting it as compatible.

- Original's BCF-to-VCF path rewrites structured header records incorrectly
  when quoted values contain commas or equals signs. It also has unusual
  missing-GT (`-1/-1`) and Character/String FORMAT last-byte/NUL behavior.
  vcftools-ng intentionally reproduces those emitted bytes in exact
  BCF-to-VCF mode. This is compatibility with the oracle, not a claim that
  the resulting header/text is standards-correct.
- Original corrupts genotypes for BCF input combined with genotype masking
  and `--recode-bcf`. vcftools-ng rejects this combination. Exact raw
  `--recode-bcf` therefore currently requires BCF input, no sample subsetting,
  and no genotype masking.
- Original's PCA path is undefined/misaligned with missing genotypes.
  vcftools-ng exact mode requires complete genotypes and reports that
  `--max-missing 1` is needed.
- Original's `--non-ref-af-any`-alone no-op behavior and partially missing
  diploid diff behavior are intentionally retained because they are observable
  compatibility semantics.
- Polyploid GT is outside the current exact fast-path contract and is rejected.

If a future corrected mode is added, it must use an explicit mode/name and
must not silently change `--compat exact` output.

## Current development caveats

- The v0.12.1 direct text parser has real 2.3-million-record gates for its
  normal diploid GT/DP/QUAL workloads plus synthetic permutations for field
  order and missing values. Malformed records, extreme numeric overflow,
  zero-sample files, and unusual polyploid encodings are not claimed as
  compatible.
- The retained differential suite proves documented workloads, not every
  option × input-format × option-interaction combination.
- The v0.12.1 full-data first-repeat gate passed 42/42 configurations. Its
  deferred repeats 2–5 are not yet part of the performance claims.
