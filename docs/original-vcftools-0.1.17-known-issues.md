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
