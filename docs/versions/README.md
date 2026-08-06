# Version archive policy

This directory is the permanent local archive for each `vcftools-ng`
development version.

## Compatibility oracle and fixtures

- Oracle: VCFtools 0.1.17.
- Main fixture: `tests/fixtures/osmanthus412.23chr_100k`, containing
  2,300,000 real records, 100,000 from each of `chr1`–`chr23`, and 412
  samples.
- Flag fixture: `tests/fixtures/osmanthus412.flags.23chr_1k`, containing
  23,000 records derived from the same real data and deterministic
  FILTER/INFO/FT annotations.
- GATK fixture: `tests/fixtures/osmanthus205.gatk.23chr_1k`, containing 23,000
  records and all 205 samples from an independent GATK callset. Its source,
  artifact, and Original-oracle hashes are locked in the fixture provenance
  and `tests/golden/SHA256SUMS`.
- Exactness: complete output files are compared with `cmp`; golden hashes are
  stored in `tests/golden/SHA256SUMS`.
- v0.11.2 completed a five-repeat, seven-scenario benchmark on all
  11,230,392 records. Later routine development returns to the 2.3-million
  subset until the next final-stage full-data gate.
- v0.12.1 completed a new full-data first-repeat gate for its seven-filter
  exact-recode workload: 42/42 candidate configurations passed. It was
  published under the explicitly staged workflow and its recorded values
  remain labeled `full_11.23m_r1`, not means.
- v0.12.2 reused the hash-locked v0.12.1 Original baselines and completed a
  four-scenario, five-repeat full-data matrix: 120/120 vcftools-ng outputs
  passed. Original was not rerun.
- v0.12.3 changes terminal help, documentation, and portable-package
  validation only. It inherits the v0.12.2 exactness and performance matrix;
  no full-data timing is relabeled as a v0.12.3 measurement.
- v0.12.4 adds standard `PREFIX.log` run records and explicit index-decision
  evidence. Scientific output bytes remain under the existing Original/golden
  gates. Its final development gate passed all 18 three-scenario/thread
  combinations against the locked 2.3-million-record Original goldens, and
  its portable package passed CentOS 7 and Ubuntu 20.04 verification. With
  release-owner approval, the full-data performance matrix is inherited from
  v0.12.2 rather than rerun.
- v0.13.0 completed a 108/108 first-repeat full-data release gate across four
  input scenarios, three same-device output/storage scenarios, and
  1/2/4/8/12/16/24/28/32 threads. The published matrix uses the validated
  first repeat only; incomplete follow-up records are excluded and the values
  are never labeled as means.

## Required record for every future version

Each version file must contain:

1. version and implementation scope;
2. newly supported parameters;
3. compatibility semantics and known original quirks;
4. fixture, sample count, input format, retained record count, and outputs;
5. full byte-comparison result;
6. original VCFtools 0.1.17 wall time;
7. `vcftools-ng` 8-thread and 16-thread wall time, CPU, and peak RSS;
8. speedup versus original;
9. exact reproduction command or benchmark script;
10. SHA/golden evidence and regression-test result.

Use `未记录` rather than estimating a historical value that was not measured.
Starting with the next version after v0.9.0, original, 8-thread, and 16-thread
measurements are mandatory in the same benchmark report.

## Archive files

- Human-readable index: `docs/VERSION_HISTORY.md`
- Machine-readable benchmark table: `docs/versions/benchmarks.tsv`
- Per-version records: `docs/versions/v*.md`
- New-version template: `docs/versions/TEMPLATE.md`
- Raw timing files: `benchmarks/results/*.time.txt`
- Detailed benchmark narrative: `benchmarks/README.md`
