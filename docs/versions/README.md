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
- v0.13.1 adds a locked 230,000-record matrix for position/sample selection,
  window pi, Tajima's D, and FST. Its 225/225 performance runs and 285/285
  23,000-record development runs are byte-identical to retained Original
  oracles. The v0.13.0 full-data recode matrix remains historical evidence and
  is not relabeled as a v0.13.1 measurement.
- v0.13.2 adds a fixed 23,000-record nine-family gate and locked 230,000-record
  A/B records for fused production filters, FILTER/INFO/FT, no-index BGZF, LD,
  exact PCA, and indexed BCF diff. It does not infer a new full-data result.
- v0.14.1 incorporates that development candidate, adds an immutable
  capability plan and the final correctness/memory/post-scan work, then passes
  a 36/36 first-repeat complete-data release gate across four input scenarios
  and nine thread counts. The v0.12.1 Original timings and scientific goldens
  are reused only after complete size/SHA-256 validation; Original is not
  rerun. Follow-up repeats remain separate from the publication gate.
- v0.14.2 runs a unified real-portable A/B against v0.13.0 and v0.14.1,
  restores private libdeflate, adds BCF-aware strict stream planning, and
  selects aligned pread for Plain VCF above 8 GiB. The initial matrix passes
  81/81 output gates; final post-fix Plain rows retain the same oracle hash.
  Performance differences within 5% are treated as effectively tied.
- v0.14.3 fixes the v0.14.2 whole-workflow CPU oversubscription defect. Linux
  process-tree affinity now enforces `--threads N`; waiting I/O/output pools
  may overlap within those N CPUs. The locked 2.3-million exact suite passes,
  and fair 230k v0.14.2 A/B deltas at 4/8/16/24/32 CPUs all remain within 5%.

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
- v0.14.1 implementation and decision audit:
  `docs/versions/v0.14.1-resume-plan.md`
- v0.14.2 portable A/B and release correction:
  `docs/versions/v0.14.2.md`
- v0.14.3 strict CPU-budget correction:
  `docs/versions/v0.14.3.md`
- New-version template: `docs/versions/TEMPLATE.md`
- Raw timing files: `benchmarks/results/*.time.txt`
- Detailed benchmark narrative: `benchmarks/README.md`
