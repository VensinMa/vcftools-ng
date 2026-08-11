# v0.14.2 unified complete-data A/B

[简体中文](README.zh-CN.md)

Status: **81/81 initial portable A/B outputs PASS** complete-file SHA-256,
size, and 5,425,725/11,230,392 retained-site gates. Original VCFtools was not
rerun. The comparison used the real v0.13.0, v0.14.1, and v0.14.2 portable
archives, one strictly serial run per row, one SSD, the same seven filters,
and uncompressed VCF output. Differences within 5% are treated as effectively
tied.

The initial A/B exposed two release-build regressions rather than a slower
parser: v0.14.1's portable HTSlib lacked libdeflate, and mapping the complete
122.9 GB Plain VCF produced about 100 GB RSS while defeating explicit-read
readahead. v0.14.2 bundles libdeflate 1.25 and maps only cacheable Plain inputs
up to 8 GiB; larger files use bounded aligned `pread` workers. The final Plain
rows are therefore in `final-plain.tsv` and supersede the initial candidate
Plain rows in `all-runs.tsv`.

Application wall seconds:

| Input | Version | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 28 | 32 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF + TBI | v0.13.0 | 360.87 | 222.47 | 114.02 | 59.59 | 44.76 | 35.91 | 31.43 | 28.83 | 28.52 |
| BGZF + TBI | v0.14.1 | 358.89 | 220.90 | 113.88 | 59.32 | 44.26 | 35.79 | 30.87 | 29.84 | 29.09 |
| BGZF + TBI | v0.14.2 | 214.43 | 151.02 | 77.44 | 40.72 | 34.21 | 31.68 | 26.07 | 28.84 | 28.25 |
| Plain VCF | v0.13.0 | 208.71 | 109.78 | 59.28 | 37.98 | 28.08 | 26.16 | 36.62 | 34.41 | 33.58 |
| Plain VCF | v0.14.1 | 208.14 | 100.72 | 61.84 | 42.61 | 37.08 | 37.85 | 38.67 | 36.58 | 36.13 |
| Plain VCF | v0.14.2 final | 214.12 | 97.63 | 53.83 | 32.51 | 27.42 | 32.05 | 29.65 | 32.72 | 31.72 |
| BCF stream | v0.13.0 | 473.81 | 238.56 | 239.38 | 122.12 | 71.12 | 60.20 | 49.21 | 47.17 | 47.29 |
| BCF stream | v0.14.1 | 665.13 | 472.36 | 472.72 | 124.28 | 68.71 | 62.06 | 63.75 | 65.19 | 77.96 |
| BCF stream | v0.14.2 | 551.39 | 470.46 | 237.79 | 98.22 | 67.59 | 52.30 | 45.45 | 44.27 | 50.43 |

v0.13.0 consumed about 142% CPU when asked for one BCF thread and 283% at two
threads. Those low-thread rows are retained as historical observations, not
strict-budget comparisons. At 32 threads, the host is at its local scheduling
and storage edge; the 24/28-thread rows are more representative of its BCF
knee. Application time and post-`sync -f` durable time remain separate in the
raw tables so a 57-59 GB writeback fluctuation is not described as a parser
regression.

Evidence:

- `all-runs.tsv`: initial 81-row portable A/B;
- `final-plain.tsv`: post-fix final-source Plain candidate rows;
- `asset-validation.tsv`: input, oracle, and release-archive identities;
- `runs/` and `logs/`: per-row timing metadata and captured diagnostics;
- `../../run-v0142-unified-full-ab.sh`: resumable driver.
