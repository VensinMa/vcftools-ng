# Parameter compatibility and optimization status

This document describes the v0.14.2 command-line surface. VCFtools 0.1.17
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

v0.14.1 introduced the immutable `QueryPlan`, retained by v0.14.2. The plan
decides required FORMAT fields, fused-kernel eligibility, general-pipeline
decode, and the logged fallback reason before scientific outputs are opened.
This changes execution cost, not the parameter compatibility categories below.

| Workload group | Original-compatible parameters | Evidence |
|---|---|---|
| Site frequency, missingness, depth, and quality | `--freq`, `--freq2`, `--counts`, `--missing-site`, `--site-depth`, `--site-mean-depth`, `--site-quality` | Real 2.3-million-record complete-file goldens. The v0.11.4 fused BGZF path passed at 8/16 threads. Six outputs in one scan measured 38.68×/46.11× against six Original scans; `--freq2` measured 7.20×/8.55×. |
| Individual and HWE statistics | `--depth`, `--missing-indv`, `--het`, `--hardy` | Real 2.3-million-record exact gates and individual/combined benchmarks in v0.10.0. |
| Diversity and population statistics | `--site-pi`, `--window-pi`, `--window-pi-step`, `--TajimaD`, `--weir-fst-pop`, `--fst-window-size`, `--fst-window-step` | v0.13.1 direct Plain VCF paths passed the 23k and 230k exact matrices at 1/4/8/16/32 threads. v0.13.2 removes chromosome-string/map work, uses integer window indexing, applies an exact integer difference array to overlapping pi, and lets `--site-pi` share the fused scan. FST/Tajima floating accumulation order is unchanged. Window analyses require nondecreasing positions and contiguous chromosome segments; unsorted inputs are rejected transactionally rather than reproducing Original's duplicated/misleading window rows. |
| LD and PCA | `--geno-r2`, `--ld-window`, `--ld-window-min`, `--ld-window-bp`, `--ld-window-bp-min`, `--min-r2`, `--pca`, `--pca-no-norm` | v0.13.2 uses bit-plane/popcount LD only for eligible complete biallelic diploid data and SoA exact PCA with fixed per-cell site order. Phased GT is decoded by allele value, while multiallelic, haploid, and incomplete sites cannot silently enter a mathematically incompatible specialized kernel. The fixed 23k 1/4/8/16/32 gate is byte-identical. |
| VCF/BCF recode | `--recode`, `--recode-bcf`, `--recode-INFO-all`, `--stdout` | Complete VCF/BCF byte gates and conversion benchmarks. File-based `--recode` changes only the container to BGZF; decompressed VCF bytes retain the Original-compatible path, while `--recode-vcf` selects a directly comparable plain file. v0.13.1 also gates selected-sample plus seven-filter direct BGZF recode. `--recode --stdout` remains plain VCF. |
| Common numeric/site filters | `--min-alleles`, `--max-alleles`, `--remove-indels`, `--keep-only-indels`, `--minQ`, `--min-meanDP`, `--max-meanDP`, `--max-missing`, `--max-missing-count`, `--maf`, `--max-maf`, `--mac`, `--max-mac`, `--hwe` | Exact filtered statistics/recode gates and v0.2/v0.5/v0.11.4 real-data benchmarks. v0.13.2 fuses all listed filters into the direct text plan, including Original's exact HWE calculation. |
| Chromosome, position, and interval filters | `--chr`, `--not-chr`, `--from-bp`, `--to-bp`, `--positions`, `--exclude-positions`, `--bed`, `--exclude-bed`, `--thin` | Complete-file exact gates and v0.3/v0.6/v0.7 benchmarks. v0.13.1 additionally gives `--positions`/`--exclude-positions` a direct Plain VCF path; 1%/50%, sorted/shuffled, duplicate, and absent-position 230k fixtures passed byte gates and measured 8.17×–20.00× at 32 threads. |
| Non-reference filters | `--non-ref-af`, `--max-non-ref-af`, `--non-ref-af-any`, `--max-non-ref-af-any`, `--non-ref-ac`, `--max-non-ref-ac`, `--non-ref-ac-any`, `--max-non-ref-ac-any` | Complete-file exact gates and v0.8 benchmark. |
| FILTER/INFO Flag filters | `--keep-filtered`, `--remove-filtered`, `--remove-filtered-all`, `--keep-INFO`, `--remove-INFO` | v0.13.2 direct FILTER/INFO Flag paths pass annotated 23k gates at 1/4/8/16/32 threads. `--keep-INFO`/`--remove-INFO` remain Flag-only and reject non-Flag header declarations. |
| Genotype filters | `--minGQ`, `--minDP`, `--maxDP`, `--remove-filtered-geno`, `--remove-filtered-geno-all` | v0.13.2 fuses GQ/DP/FT masking into the same sample pass and passes counts, missing-site, and recode byte gates at 1/4/8/16/32 threads. |
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

All four outputs have complete-file byte gates. v0.13.2 additionally
parallelizes indexed BCF site/individual discordance by contig and commits in
the first header's order; unsupported shapes retain the serial exact path. A
real 2.3-million-record
performance benchmark exists for `--diff-site-discordance` (15.18×/15.14× at
8/16 threads). The other three share the new merge/decode engine but have not
been benchmarked separately, so no independent speedup claim is made for
them. Diff filtering is currently limited to chromosome/position/BED and
sample selection.

## vcftools-ng-only extensions

These options are additions, not Original-compatible parameter names:

| Extension | Purpose and validation |
|---|---|
| `--threads N`, `-t N` | Sets the process-tree CPU budget. On Linux vcftools-ng and automatic bcftools children inherit an affinity set of at most N CPUs; waiting I/O workers may overlap. When omitted, scheduler, affinity, cgroup, and hardware limits are intersected and automatic selection is capped at 128. |
| `--batch-size N` | Tunes the bounded generic pipeline batch size. |
| `--input FILE` | Auto-detected input alias for `--vcf`/`--gzvcf`/`--bcf`. |
| `--compat exact` | Explicitly selects the only currently implemented compatibility mode. |
| `--input-backend auto\|stream\|plain\|indexed` | Selects or forces an input adapter. `auto` is the default adaptive speed policy. |
| `--bcftools FILE` | Selects the bcftools executable used when the adaptive policy decides CSI construction is profitable. |
| `--recode-vcf-gz` | Writes deterministic, parallel BGZF-compressed VCF. Original 0.1.17 has no compressed-VCF recode option. Validation compares the decompressed bytes with Original `--recode`; compressed bytes are deterministic across tested thread counts, but cannot be compared with a nonexistent Original artifact. |
| `--recode-vcf` | Writes an uncompressed VCF file after v0.13.0 changed the default `--recode` file output to BGZF. `--recode --stdout` remains plain VCF for compatibility. |
| `--log-file FILE` | Overrides the default `PREFIX.log` run-log path. The log is overwritten and receives the same diagnostics as stderr. |
| `--no-log-file` | Disables the run-log file without disabling stderr diagnostics. |
| `--corrected-depth-arithmetic` | With `--site-depth` and/or `--site-mean-depth`, uses checked 64-bit sums instead of Original unsigned-32-bit wrapping. Disabled by default; exact outputs retain OVI-013. Corrected artifacts are separately hash-locked extension evidence. |

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
- Polyploid GT is outside the current exact statistical contract. Commands
  that require GT semantics reject it; site-only filtering plus raw VCF
  recode preserves the sample field as opaque text (OVI-012).

The direct VCF parser is not limited to `0/0`, `0/1`, `1/1`, `./.`, `0|0`,
`0|1`, and `1|1`. The real DeepVariant fixture also contains `./0` and `./1`;
the real GATK fixture contains phased calls, haploid `.`, and multiallelic
diploids using allele indices 2-4. A locked synthetic Original gate further
covers reverse heterozygotes, both orientations of partial missingness,
phased partial missingness, called haploids, and multi-digit allele indices.
See the [GT compatibility audit](../benchmarks/results/v0132-development-gate/GT_COMPATIBILITY.md).
The `--missing-site` gate also preserves Original's phased trailing-missing
haploid interpretation (OVI-011) without changing allele counts or filters.

### GT grammar x analysis semantics

"Exact" below means a dedicated Original 0.1.17 byte gate exists. "Boundary"
means the command does not silently enter a diploid-only mathematical kernel;
it is skipped, rejected, or otherwise limited as stated.

| GT grammar | Counts | Missingness | VCF recode | HWE | FST | LD | PCA |
|---|---|---|---|---|---|---|---|
| Haploid | Exact | Exact | Preserved | Not claimed | Fully-diploid boundary | Complete-diploid boundary | Complete-diploid boundary |
| Diploid, unphased | Exact | Exact | Exact | Exact for biallelic complete calls | Exact | Exact for eligible biallelic sites | Exact for complete eligible sites |
| Diploid, phased | Exact | Exact, including OVI-011 | Exact | Allele/dosage compatible | Exact | Phase-independent dosage | Phase-independent dosage |
| Partially missing diploid | Exact | Exact with legacy quirks | Exact | Command-specific boundary | Command-specific boundary | Incomplete-call boundary | Rejected by exact PCA |
| Multiallelic diploid | Exact | Exact | Exact | Biallelic boundary | Generic exact FST path | Biallelic boundary | Eligibility boundary |
| Ploidy >= 3 | Semantic commands reject | Semantic commands reject | Site-only opaque passthrough; masking rejects | Reject | Reject | Reject | Reject |

The test-only environment variable `VCFTOOLS_NG_TEST_PARSER` accepts `auto`,
`generic`, or `specialized`. Release gates compare generic and specialized
scalar results against each other and against Original-locked DeepVariant,
GATK, and synthetic GT artifacts. It is development infrastructure, not a
supported user-facing tuning option.

If a future corrected mode is added, it must use an explicit mode/name and
must not silently change `--compat exact` output.

## Current development caveats

- v0.14.1 adds directed malformed/truncated input, numeric overflow, NaN/Inf,
  zero-selected-sample, no-contig-header, oversized-record, GT grammar, and
  polyploid boundary tests. These fixtures define documented behavior; they
  are not a claim that genotype-semantic analyses support ploidy above two.
- The retained differential suite proves documented workloads, not every
  option × input-format × option-interaction combination.
- The v0.14.1 full-data first-repeat release gate passed 36/36 configurations
  over four input scenarios and nine thread counts. It is explicitly
  single-run evidence; adaptive follow-up repeats are not yet presented as
  means.
