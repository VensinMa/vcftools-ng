# vcftools-ng version history

All exactness claims use VCFtools 0.1.17 as the oracle. Unless stated
otherwise, validation uses the real 2,300,000-record, 23-chromosome,
412-sample subset. Times are wall-clock seconds.

| Version | Main addition | Original | 8 threads | 16 threads | Exact |
|---|---|---:|---:|---:|---|
| [v0.1.0](versions/v0.1.0.md) | Five statistics in one scan | 352.61 | 7.55 | 6.56 | PASS |
| [v0.2.0](versions/v0.2.0.md) | Common filters and VCF recode | 156.14 | 5.42 | 4.58 | PASS |
| [v0.3.0](versions/v0.3.0.md) | Position/chromosome selection and CSI query | 48.88 | 0.24 | 0.19 | PASS |
| [v0.4.0](versions/v0.4.0.md) | Sample inclusion/exclusion | 623.49 | 4.08 | 3.73 | PASS |
| [v0.5.0](versions/v0.5.0.md) | MAC and exact HWE | 104.85 | 4.01 | 3.64 | PASS |
| [v0.6.0](versions/v0.6.0.md) | BED include/exclude | 40.95 | 3.50 | 3.39 | PASS |
| [v0.7.0](versions/v0.7.0.md) | Ordered `--thin` | 39.94 | 3.49 | 3.36 | PASS |
| [v0.8.0](versions/v0.8.0.md) | Non-reference AF/AC filters | 41.02 | 3.47 | 3.36 | PASS |
| [v0.9.0](versions/v0.9.0.md) | FILTER/INFO/FT and overlapped pipeline | 156.14 | 3.37 | 3.37 | PASS |

The workload in each row is the version's representative compatibility
benchmark; workloads differ between rows. Consult the per-version page before
comparing versions directly. v0.9.0 also reran the six-output sample workload
at 3.11/3.13 seconds for 8/16 threads.

## Cumulative supported surface at v0.9.0

- Inputs: VCF, BGZF VCF, BCF.
- Outputs: `--freq`, `--counts`, `--missing-site`, `--site-depth`,
  `--site-mean-depth`, `--recode`, `--recode-INFO-all`.
- Position/site selection: chromosome, closed range, position lists, BED,
  ordered thinning.
- Samples: file-based and inline inclusion/exclusion.
- Site/genotype filters: allele number/type, quality, depth, missingness,
  MAF/MAC/HWE, non-reference AF/AC, FILTER, INFO Flag, and FORMAT/FT.
- Execution: bounded reader → compute → ordered-commit pipeline with 8/16
  thread compatibility gates.

Machine-readable values are in [benchmarks.tsv](versions/benchmarks.tsv).

