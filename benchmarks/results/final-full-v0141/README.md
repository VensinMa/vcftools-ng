# v0.14.1 complete-data release gate

[简体中文](README.zh-CN.md)

Status: **36/36 first-repeat candidate outputs PASS complete-file `cmp`**.

- Input: 11,230,392 records, 412 samples.
- Retained: 5,425,725 records in every run.
- Workload: seven Original-compatible filters plus `--recode
  --recode-INFO-all --stdout`.
- Scenarios: BGZF+TBI, BGZF+automatic CSI, Plain VCF, BCF adaptive stream.
- Threads: `1 2 4 8 12 16 24 28 32`; strictly serial execution.
- Original: locked v0.12.1 VCFtools 0.1.17 scientific goldens and one-run
  timings; zero Original reruns for v0.14.1.
- Automatic CSI: bcftools 1.24 / HTSlib 1.24.
- Candidate: local opt-in PGO binary SHA-256
  `0d585fe4e146feac1b63dd93a0cd238b2cd5b1baaa1d51a21a367c8f4b0b714e`.

Single-run wall seconds:

| Input | Original | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF + TBI | 2267.88 | 203.72 | 146.86 | 75.48 | 44.36 | 33.47 | 34.06 | 32.43 | 30.93 | 32.84 |
| BGZF + automatic CSI | 2267.88 | 206.13 | 248.91 | 132.26 | 98.61 | 89.97 | 89.53 | 89.16 | 87.79 | 86.22 |
| Plain VCF | 2092.91 | 204.35 | 97.04 | 58.69 | 44.39 | 42.13 | 43.62 | 41.89 | 42.52 | 41.69 |
| BCF adaptive stream | 1943.47 | 549.11 | 466.67 | 466.99 | 120.01 | 62.82 | 50.46 | 44.36 | 43.67 | 45.00 |

These are first-repeat values, not means. Automatic CSI includes index
construction. Strict-budget BCF is record-decode limited at 1-4 threads; the
observed local 28/32-thread edges are not runtime caps or scaling claims.

Compact permanent evidence:

- [`asset-validation.tsv`](asset-validation.tsv): input/index/golden
  size and SHA-256 gate;
- [`manifest.tsv`](manifest.tsv): candidate, workload, environment, and policy;
- [`all-runs.tsv`](all-runs.tsv): one row per Original/candidate run;
- [`summary.tsv`](summary.tsv): wall, speedup, CPU, peak RSS, backend, exactness;
- [`run-v0141-full-release-matrix.sh`](../../run-v0141-full-release-matrix.sh):
  resumable release driver.

Large inputs, retained Original goldens, per-run logs, and transient 57-59 GB
candidate outputs stay on the benchmark host. Successful scratch outputs were
removed only after byte comparison; the final scratch directory is empty.
