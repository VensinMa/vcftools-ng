# v0.13.0 full-data input/output/storage release gate

Status: **single-repeat release gate PASS (108/108); published timings use
the validated first repeat only**.

This gate used all 11,230,392 records and retained 5,425,725 records with the
seven-filter workload recorded in `manifest.tsv`. Original VCFtools 0.1.17 was
not rerun: the driver hash-validated and reused the retained v0.12.1 VCF
goldens and timings. Every v0.13.0 output passed complete-content SHA-256,
record-count, and retained-record-count gates. Plain VCF bytes matched the
corresponding Original output directly; BGZF output was validated after
decompression and was deterministic across thread counts.

The host has 32 logical Intel Core i9-14900KF CPUs and 125 GiB RAM. SSD input
and output use `/dev/nvme3n1p2`; HDD input and output use the same rotational
`/dev/sda1` filesystem. Values below are first-repeat application wall seconds,
not multi-run means. `summary.tsv` also records durable wall time, flush time,
CPU utilization, peak RSS, output bytes, backend, and exactness.

## SSD uncompressed VCF output

| Input | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 243.46 | 166.44 | 85.78 | 44.09 | 34.27 | 28.50 | 27.69 | 29.28 | 28.50 |
| BGZF VCF + automatic CSI | 241.42 | 264.60 | 140.10 | 100.80 | 91.84 | 85.66 | 82.55 | 87.44 | 85.75 |
| Plain VCF | 241.25 | 126.93 | 69.34 | 42.25 | 41.53 | 41.99 | 40.01 | 40.39 | 40.70 |
| BCF adaptive | 317.31 | 159.73 | 160.44 | 82.43 | 48.10 | 43.13 | 35.91 | 35.65 | 37.02 |

## SSD BGZF VCF output

| Input | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 1229.88 | 528.12 | 270.56 | 150.06 | 114.85 | 91.27 | 72.80 | 72.04 | 71.30 |
| BGZF VCF + automatic CSI | 1221.37 | 622.72 | 324.19 | 205.69 | 171.43 | 149.57 | 128.60 | 127.71 | 128.74 |
| Plain VCF | 1018.63 | 509.09 | 257.79 | 134.33 | 104.59 | 87.78 | 66.65 | 67.24 | 67.14 |
| BCF adaptive | 1019.48 | 511.56 | 265.43 | 160.21 | 123.11 | 98.45 | 83.54 | 80.16 | 81.39 |

## Speedup over Original for SSD BGZF VCF output

| Input | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 1.84× | 4.29× | 8.38× | 15.11× | 19.75× | 24.85× | 31.15× | 31.48× | 31.81× |
| BGZF VCF + automatic CSI | 1.86× | 3.64× | 7.00× | 11.03× | 13.23× | 15.16× | 17.64× | 17.76× | 17.62× |
| Plain VCF | 2.05× | 4.11× | 8.12× | 15.58× | 20.01× | 23.84× | 31.40× | 31.13× | 31.17× |
| BCF adaptive | 1.91× | 3.80× | 7.32× | 12.13× | 15.79× | 19.74× | 23.26× | 24.24× | 23.88× |

These ratios use the locked Original plain-VCF wall time as the end-to-end
workflow baseline. The filtering task and decompressed scientific content are
equivalent, but the output encoding and artifact size differ.

## Same-HDD input and BGZF VCF output

| Input | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 1322.12 | 524.01 | 270.81 | 152.87 | 115.65 | 95.05 | 73.03 | 72.35 | 71.47 |
| BGZF VCF + automatic CSI | 1225.84 | 626.17 | 326.51 | 209.94 | 174.60 | 150.75 | 129.77 | 130.50 | 128.87 |
| Plain VCF | 1076.97 | 823.98 | 826.74 | 869.95 | 883.69 | 909.75 | 997.62 | 1029.14 | 1059.79 |
| BCF adaptive | 1065.09 | 514.18 | 265.73 | 160.99 | 123.46 | 98.92 | 84.01 | 80.37 | 81.93 |

## Reference speedup over Original for same-HDD BGZF output

| Input | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 1.72× | 4.33× | 8.37× | 14.84× | 19.61× | 23.86× | 31.05× | 31.35× | 31.73× |
| BGZF VCF + automatic CSI | 1.85× | 3.62× | 6.95× | 10.80× | 12.99× | 15.04× | 17.48× | 17.38× | 17.60× |
| Plain VCF | 1.94× | 2.54× | 2.53× | 2.41× | 2.37× | 2.30× | 2.10× | 2.03× | 1.97× |
| BCF adaptive | 1.82× | 3.78× | 7.31× | 12.07× | 15.74× | 19.65× | 23.13× | 24.18× | 23.72× |

These requested reference ratios divide the locked Original SSD plain-VCF
time by the same-HDD vcftools-ng BGZF time. They therefore include both
output-encoding and storage-device differences and are not a controlled
same-device comparison.

The 59.43 GB Plain VCF result compressed to a deterministic 10.20 GB BGZF
VCF, an 82.8% size reduction. The Plain-VCF-on-HDD rows show the expected
same-device I/O ceiling: reading the 122.91 GB uncompressed input and writing
compressed output becomes slower at high concurrency. The result is retained
rather than hidden because it defines the storage-dependent scaling boundary.

## Storage-aware usage recommendation

- For high-performance SSD/NVMe storage with ample capacity, `--recode-vcf`
  is the preferred minimum-wall-time mode because it avoids output
  compression. The tradeoff is the 59.43 GB Plain VCF artifact.
- For conventional or slower HDD storage, keep the default `--recode` or use
  the equivalent `--recode-vcf-gz`, and prefer BGZF VCF input with a valid
  TBI/CSI index. The 10.20 GB BGZF result reduces output traffic by 82.8%,
  while compressed indexed input avoids the severe same-disk Plain VCF I/O
  ceiling visible above.

This recommendation is specific to the measured workload and storage classes;
filesystem caching, shared-storage contention, compression ratio, and later
access patterns can move the crossover point.

## Speedup over Original for SSD Plain VCF output

| Input | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF VCF + TBI | 9.32× | 13.63× | 26.44× | 51.44× | 66.18× | 79.57× | 81.90× | 77.45× | 79.57× |
| BGZF VCF + automatic CSI | 9.39× | 8.57× | 16.19× | 22.50× | 24.69× | 26.48× | 27.47× | 25.94× | 26.45× |
| Plain VCF | 8.68× | 16.49× | 30.18× | 49.54× | 50.39× | 49.84× | 52.30× | 51.81× | 51.42× |
| BCF adaptive | 6.12× | 12.17× | 12.11× | 23.58× | 40.40× | 45.06× | 54.12× | 54.52× | 52.50× |

The uncompressed ratios compare identical output encoding. The BGZF ratios
compare end-to-end completion of the same filtering task against Original's
plain-VCF output; VCFtools 0.1.17 has no `--recode-vcf-gz` option. HDD ratios
are reference values against the SSD Original baseline, not same-device
measurements.

## Repeat policy and published run count

The release driver runs at most three repeats. A row whose first application
wall time exceeds 1,800 seconds runs once. Otherwise it runs twice; if either
run exceeds 600 seconds and the symmetric difference is below 10%, repeat
three is skipped. All faster or more variable rows run three times. Each skip
is retained as a machine-readable record, and the summary reports the actual
run count. For the final v0.13.0 publication, optional follow-up repeats were
stopped and excluded; every table in this report therefore uses the validated
first repeat (`runs=1`) consistently.

Reproduce with:

```bash
GATE_ONLY=1 benchmarks/run-v0130-output-storage-matrix.sh
GATE_ONLY=0 benchmarks/run-v0130-output-storage-matrix.sh
```
