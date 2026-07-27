# v0.12.4 standard-logging final development gate

Status: **PASS**.

## Workload and oracle

- real 2,300,000-record, 23-chromosome, 412-sample subset;
- seven real-project filters plus full INFO-preserving VCF recode;
- BGZF VCF + TBI, Plain VCF, and BCF adaptive stream;
- 1/2/4/8/16/32 threads;
- retained VCFtools 0.1.17 goldens and timings;
- Original reruns: zero;
- complete-file `cmp` plus SHA-256 for every candidate.

All 18/18 candidate outputs were byte-identical and retained 1,121,342 of
2,300,000 sites.

## Results

Wall-clock seconds:

| Scenario | Original | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 392.73 | 61.55 | 40.99 | 23.47 | 15.72 | 10.63 | 8.62 |
| Plain VCF | 353.46 | 59.34 | 58.63 | 17.72 | 12.85 | 9.55 | 9.93 |
| BCF adaptive stream | 327.22 | 66.83 | 33.57 | 22.46 | 11.62 | 8.02 | 7.81 |

Speedup over the locked Original baseline:

| Scenario | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 6.38× | 9.58× | 16.73× | 24.98× | 36.95× | 45.56× |
| Plain VCF | 5.96× | 6.03× | 19.95× | 27.51× | 37.01× | 35.60× |
| BCF adaptive stream | 4.90× | 9.75× | 14.57× | 28.16× | 40.80× | 41.90× |

## Logging assertions

The committed logs demonstrate:

- BGZF VCF validates and uses the existing TBI;
- Plain VCF records that CSI/TBI is not applicable and selects stream at
  1–2 threads or aligned ranges at 3+;
- BCF validates and preserves the existing CSI but selects stream for
  full-file recode and records `Index used: no`;
- terminal stderr and the selected log file use the same central logger;
- every successful run records resource use and `Exit status: success`.

The exact binary exercised by this gate has SHA-256:
`46e27de325bba7cdc02ea8ed6031a647e9bc6fc59b7c112366a1a77da62f491e`.

Machine-readable evidence:

- `summary.tsv`
- `all-runs.tsv`
- `manifest.tsv`
- `runs/*.tsv`
- `runs/*.time.txt`
- `logs/*.stderr.txt`
