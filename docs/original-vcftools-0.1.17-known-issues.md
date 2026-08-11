# Original VCFtools 0.1.17 known issues / 已知问题档案

This file records behavior in Original VCFtools 0.1.17 that is corrupt,
undefined, standards-inconsistent, or operationally surprising. It prevents
“exact compatibility” from silently being interpreted as “the Original
behavior is correct.”

本文单独记录 Original VCFtools 0.1.17 中已确认的错误、未定义行为、偏离
格式规范或明显不合理的历史语义。“精确兼容”只表示输出与 oracle 一致，
不表示 Original 的行为一定正确。

## Policy / 处理原则

- **Exact-compatible quirk / 兼容保留：** if output is valid and the behavior
  is observable, `--compat exact` may reproduce it byte for byte.
- **Unsafe Original bug / 拒绝继承：** if Original emits corrupt or undefined
  output, vcftools-ng rejects the combination with an explicit error.
- **Corrected behavior / 修正版：** a future standards-correct alternative
  must use an explicit option or compatibility mode; it must never silently
  change an established exact golden.
- Every new entry should include a minimal reproducer, Original version and
  binary hash where practical, observed output, vcftools-ng decision, and a
  regression/golden reference.

## Confirmed entries / 已确认条目

### OVI-001 — BCF structured-header corruption

- **Class:** confirmed Original defect; exact text quirk retained where safe.
- **Trigger:** BCF input converted to VCF when quoted structured-header values
  contain commas or equals signs.
- **Original behavior:** fields are split/reprinted incorrectly, producing a
  malformed structured header.
- **vcftools-ng:** exact BCF-to-VCF mode reproduces the observed oracle bytes
  for covered goldens. This is not described as standards-correct output.
- **Evidence:** BCF conversion cases in `tests/differential.sh` and the
  compatibility matrix.

### OVI-002 — BCF missing-GT and Character/String FORMAT byte quirks

- **Class:** confirmed Original defect/legacy encoding behavior; retained only
  for exact BCF-to-VCF compatibility.
- **Trigger:** BCF input containing missing GT or Character/String FORMAT
  values.
- **Original behavior:** missing GT may be printed as `-1/-1`; Character/String
  values can expose last-byte or embedded-NUL artifacts.
- **vcftools-ng:** covered exact mode reproduces the bytes. VCF-input goldens
  are preferred when a valid textual oracle is required.
- **Evidence:** `tests/fixtures/README.md` and BCF conversion differential
  cases.

### OVI-003 — Corrupt BCF after genotype masking

- **Class:** confirmed unsafe Original defect; deliberately rejected.
- **Trigger:** BCF input plus a genotype filter such as `--minGQ`, `--minDP`,
  `--maxDP`, or genotype FT removal, together with `--recode-bcf`.
- **Original behavior:** emitted BCF genotypes are corrupt/untrustworthy.
- **vcftools-ng:** rejects this combination. Exact raw `--recode-bcf` requires
  BCF input without genotype masking or sample subsetting.
- **Evidence:** option validation in `src/main.cpp` and the BCF recode gate.

### OVI-004 — PCA with missing genotypes

- **Class:** confirmed undefined/misaligned Original behavior; deliberately
  rejected in exact mode.
- **Trigger:** PCA input contains missing genotypes.
- **Original behavior:** matrix cells and sample/site contributions can become
  undefined or misaligned.
- **vcftools-ng:** exact PCA requires complete genotypes and directs users to
  filter with `--max-missing 1`.
- **Evidence:** PCA compatibility gate and `docs/versions/v0.10.0.md`.

### OVI-005 — `--non-ref-af-any` alone is a no-op

- **Class:** confirmed surprising logic; retained as observable exact
  semantics.
- **Trigger:** `--non-ref-af-any` is used without the corresponding maximum
  form in the affected Original branch.
- **Original behavior:** the requested lower-bound filter does not take
  effect.
- **vcftools-ng:** `--compat exact` retains the no-op behavior for the covered
  parameter interaction. A corrected implementation would need an explicit
  mode.
- **Evidence:** non-reference filter cases in `tests/differential.sh`.

### OVI-006 — Partially missing diploid diff semantics

- **Class:** confirmed surprising comparison behavior; retained for exact
  compatibility.
- **Trigger:** two-file diff encounters diploid genotypes such as `0/.`,
  `./1`, or comparable partial missingness.
- **Original behavior:** missing/called and discordance decisions do not follow
  a simple allele-wise missingness rule.
- **vcftools-ng:** preserves the observed branch decisions in exact diff
  outputs.
- **Evidence:** flag/diff semantic cases in `tests/differential.sh`.

### OVI-007 — Filtered GT is always rewritten as `./.`

- **Class:** valid VCF but lossy and surprising legacy behavior; retained for
  exact compatibility.
- **Trigger:** a genotype is rejected by `--minGQ` (and equivalently by other
  genotype-masking filters) during VCF recode.
- **Original behavior:** the GT becomes unphased diploid missing `./.` even if
  the input was phased (`1|1`), haploid (`1`), or a single missing allele
  (`.`). Other FORMAT values remain present.
- **vcftools-ng:** the direct recode kernel deliberately writes `./.`. It does
  not preserve original phase or ploidy for a filtered GT.
- **Minimal oracle:** `tests/fixtures/fast-recode.vcf` with
  `--minGQ 10 --recode --recode-INFO-all`.
- **Golden:** `tests/golden/fast-recode.info.vcf`, generated by Original
  VCFtools 0.1.17 and checked at 1/4/8/16/32 threads.

### OVI-008 — Spurious warnings for valid comma-rich INFO descriptions

- **Class:** confirmed parser/warning defect; not reproduced.
- **Trigger:** valid VCF structured INFO header lines whose quoted
  `Description` contains commas, including GATK `MLEAC` and `MLEAF` lines.
- **Original behavior:** emits `Expected at least 2 parts in INFO entry`
  warnings even though the header is valid. In the covered VCF-input case,
  the header and scientific record output are still preserved.
- **vcftools-ng:** accepts these valid quoted descriptions without the
  misleading warning. Exact scientific-output compatibility does not require
  duplicating diagnostic text.
- **Evidence:** `osmanthus205.gatk.23chr_1k` Original oracle stderr retained as
  `tests/golden/gatk205-seven-filter.{counts,recode}.stderr`; the source and
  oracle identities are recorded in the fixture provenance file.

### OVI-009 — Corrupt raw BCF recode with non-canonical FORMAT encodings

- **Class:** confirmed unsafe Original defect; exact-compatibility gap found by
  the v0.13.0 23k workload matrix.
- **Trigger:** BCF input containing the legacy Character/String FORMAT
  encodings covered by `osmanthus412.flags.23chr_1k.bcf`, followed by raw
  `--recode-bcf --recode-INFO-all`; genotype masking is not required.
- **Original behavior:** exits with status zero but emits an unreadable BCF.
  HTSlib reports `Invalid FORMAT type 15` at `chr1:5330` and aborts with a BCF
  read error. The reproducing Original 0.1.17 binary SHA-256 is
  `8950bcdc1900e6c86df93c39502d46752ec5bdaa01426b86f56cbe94c14fae15`.
- **vcftools-ng v0.13.0:** emits valid standards-readable BCF rather than the
  corrupt Original bytes. This trigger must not be advertised as exact BCF
  compatibility; a follow-up should either reject it in exact mode or expose
  corrected output through an explicit non-exact mode.
- **Benchmark decision:** W11 BCF uses the canonical normalized 23k BCF for
  which Original and vcftools-ng outputs match byte for byte. The corrupt
  trigger and smoke artifacts are retained locally under
  `benchmarks/results/workload-matrix-23k-v0130/smoke-bcf-copy.*`.

### OVI-010 — `--geno-r2` segfault after writing output with RNC FORMAT

- **Class:** confirmed Original process-lifetime defect; not reproduced.
- **Trigger:** the 23k DeepVariant-derived VCF containing the
  `GT:DP:AD:GQ:PL:RNC:FT` layout and its `RNC` Character FORMAT declaration,
  with `--geno-r2 --ld-window 10 --ld-window-bp 100000`.
- **Original behavior:** writes the `.geno.ld` artifact, then terminates with
  SIGSEGV (exit status 139) instead of completing successfully. Adding the
  explicit biallelic filters does not prevent the crash. The reproducing
  Original 0.1.17 binary SHA-256 is
  `8950bcdc1900e6c86df93c39502d46752ec5bdaa01426b86f56cbe94c14fae15`.
- **vcftools-ng:** completes normally. The v0.13.2 development gate locks the
  already Original-validated v0.13.1 LD artifact as its successful reference,
  rather than treating a crashing command as a valid oracle generator.
- **Evidence:** `benchmarks/run-v0132-development-gate.sh` and the retained
  local 23k LD A/B records under `benchmarks/results/v0132-ld-bitset-ab/`.

### OVI-011 — Phased trailing-missing GT is counted as haploid in `.lmiss`

- **Class:** standards-inconsistent but valid-output legacy behavior; retained
  only for exact site-missingness output.
- **Trigger:** textual diploid GT whose separator is `|` and second allele is
  missing, such as `0|.` or `.|.`.
- **Original behavior:** `--missing-site` subtracts the second chromosome and
  also subtracts its missing count, interpreting the trailing missing allele
  as the BCF haploid sentinel. Thus `0|.` contributes `N_DATA=1,N_MISS=0` and
  `.|.` contributes `N_DATA=1,N_MISS=1`. Allele counts and `--max-missing`
  still treat the same three-character GT as diploid, so Original is
  internally inconsistent across outputs.
- **vcftools-ng:** exact `.lmiss` reproduces this observable rule without
  changing allele counts, ploidy-dependent analyses, or missingness filters.
- **Minimal oracle:** `tests/fixtures/genotype-forms.vcf`.
- **Goldens:** `tests/golden/genotype-forms.{frq.count,lmiss}`, generated by
  Original VCFtools 0.1.17 and checked by `tests/adaptive_counts.sh`.

### OVI-012 — Polyploid GT blocks site-only raw VCF recode

- **Class:** unnecessary Original architectural limitation; corrected safe
  passthrough extension.
- **Trigger:** a VCF contains ploidy above two and requests only site-level
  selection/filtering such as `--positions` or `--minQ` plus raw VCF recode.
- **Original behavior:** genotype parsing runs even though no GT semantics are
  needed, then exits with `Polyploidy found, and not supported by vcftools`
  before writing variant records.
- **vcftools-ng:** the fused site-only path treats sample fields as opaque text
  and preserves them byte-for-byte. Counts, missingness, genotype masking,
  HWE, FST, LD, PCA, and other GT-semantic paths still reject ploidy above two
  before committing scientific outputs.
- **Minimal fixture:** `tests/fixtures/fast-counts-polyploid.vcf`, covering
  `0/0/1`, `0|1|1`, `./1/2`, `0/./2`, and `10/10/11`.
- **Regression:** `tests/adaptive_counts.sh` checks site-only passthrough at
  1/32 threads and semantic rejection for every record. Because Original
  produces no successful artifact, passthrough goldens are explicitly marked
as vcftools-ng extension evidence rather than exact Original goldens.

### OVI-013 — Site-depth arithmetic wraps at 32 bits

- **Class:** confirmed numeric limitation retained by exact compatibility.
- **Trigger:** `--site-depth` or `--site-mean-depth` with DP values whose sum
  or sum of squares exceeds an unsigned 32-bit integer. Degenerate sites with
  zero or one reported DP also expose Original's NaN variance formatting.
- **Original behavior:** accumulates both `SUM_DEPTH` and `SUMSQ_DEPTH` in
  `unsigned int`. For two samples with `DP=100000`, `SUMSQ_DEPTH` is the
  wrapped value `2820130816`; the corresponding variance is
  `-1.71799e+10`. A site with no DP emits `-nan` for mean and variance, and a
  site with one DP emits `-nan` variance.
- **vcftools-ng policy:** exact mode deliberately preserves these output
  bytes. `--corrected-depth-arithmetic` is an explicit non-exact extension
  that uses checked 64-bit sums and sum-of-squares; it fails rather than
  silently overflowing that representation.
- **Minimal oracle:** `tests/fixtures/numeric-edge.vcf`.
- **Goldens:** `tests/golden/numeric-edge.{ldepth,ldepth.mean}`, generated by
  Original VCFtools 0.1.17 and checked in both fused and generic
  implementations by `tests/reliability_regression.sh`. Corrected extension
  artifacts are separately hash-locked and are not labeled Original goldens.

## Candidate-entry template / 新条目模板

```text
OVI-NNN — Short name
Class: confirmed defect | undefined | surprising legacy behavior | unverified
Trigger:
Original command and binary SHA-256:
Observed result:
Why it is unreasonable/incorrect:
vcftools-ng policy: retain | reject | explicit corrected mode
Minimal fixture and golden:
Regression test:
```

Unverified suspicions must not be promoted to this confirmed section until a
minimal Original 0.1.17 reproduction is preserved.
