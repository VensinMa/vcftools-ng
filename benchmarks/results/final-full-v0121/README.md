# v0.12.1 full-data first-repeat release gate

Status: **PASS for all first-repeat gates; repeats 2–5 are pending**.

## Method

- Candidate source: `f33751479f875a962039d92bd4ba17753d5bef15`
- Oracle: VCFtools 0.1.17
- Input: 11,230,392 records, 23 chromosomes, 412 samples
- Retained output records: 5,425,725
- Workload: `--min-alleles 2 --max-alleles 2 --minGQ 10 --minQ 30
  --min-meanDP 7 --max-missing 0.9 --maf 0.1 --recode
  --recode-INFO-all --stdout`
- Scenarios: BGZF+TBI, BGZF automatic CSI, BGZF no automatic index,
  Plain VCF, BCF+CSI, BCF automatic CSI, and BCF no automatic index
- Threads: 1, 2, 4, 8, 16, and 32
- Execution: strictly serial; operating-system cache was not flushed
- Exactness: complete-file `cmp`, plus output size and SHA-256 recording

Original was run once for each actual input format. The corresponding format
oracle was shared by scenarios that differ only in vcftools-ng index policy.
All 42 candidate outputs were byte-identical. Automatic-CSI rows used a new
path without a sidecar and include index construction.

## First-repeat wall time

Seconds from the 32-CPU release host:

| Scenario | Original | ng1 | ng2 | ng4 | ng8 | ng16 | ng32 | ng32 speedup |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BGZF + TBI | 2267.88 | 387.00 | 204.35 | 118.08 | 77.83 | 53.01 | 47.01 | 48.24× |
| BGZF + automatic CSI | 2267.88 | 552.17 | 290.43 | 170.28 | 129.52 | 106.91 | 100.88 | 22.48× |
| BGZF, no automatic index | 2267.88 | 304.39 | 304.35 | 255.95 | 255.01 | 257.05 | 262.18 | 8.65× |
| Plain VCF | 2092.91 | 287.87 | 294.76 | 101.53 | 71.58 | 49.78 | 52.88 | 39.58× |
| BCF + CSI | 1943.47 | 321.28 | 323.02 | 162.08 | 109.14 | 58.58 | 42.08 | 46.19× |
| BCF + automatic CSI | 1943.47 | 459.21 | 387.58 | 198.52 | 126.74 | 68.56 | 52.86 | 36.77× |
| BCF, no automatic index | 1943.47 | 319.97 | 161.73 | 109.67 | 56.70 | 41.24 | 40.41 | 48.09× |

Complete speedup matrix:

| Scenario | ng1 | ng2 | ng4 | ng8 | ng16 | ng32 |
|---|---:|---:|---:|---:|---:|---:|
| BGZF + TBI | 5.86× | 11.10× | 19.21× | 29.14× | 42.78× | 48.24× |
| BGZF + automatic CSI | 4.11× | 7.81× | 13.32× | 17.51× | 21.21× | 22.48× |
| BGZF, no automatic index | 7.45× | 7.45× | 8.86× | 8.89× | 8.82× | 8.65× |
| Plain VCF | 7.27× | 7.10× | 20.61× | 29.24× | 42.04× | 39.58× |
| BCF + CSI | 6.05× | 6.02× | 11.99× | 17.81× | 33.18× | 46.19× |
| BCF + automatic CSI | 4.23× | 5.01× | 9.79× | 15.33× | 28.35× | 36.77× |
| BCF, no automatic index | 6.07× | 12.02× | 17.72× | 34.28× | 47.13× | 48.09× |

Every candidate configuration was faster than Original in this first repeat.
The observed speedup range was 4.11×–48.24×. Plain VCF and indexed BCF show
small single-repeat reversals between adjacent thread counts; no monotonic
scaling claim is made until repeats 2–5 are complete.

Automatic CSI was faster than no-index BGZF from two threads upward. No-index
BCF was faster than automatic CSI at every tested thread count. This report
does not change the default automatic-index policy; users can explicitly
select `--no-auto-index`.

## Retained oracle identities

| Original input format | Bytes | SHA-256 |
|---|---:|---|
| BGZF VCF | 59,434,159,204 | `7548416e01d4a318b81c5d1feb9429f60c7995205d66169242c3792af4c4fc14` |
| Plain VCF | 59,434,159,621 | `d4f2a15e8c5ad0cc12abf4a3ab308bb48f22adf8ec13776d72ceeaa1f8d402b8` |
| BCF | 57,211,771,106 | `dde9edd98d5d05aa885e0e2f78a9696cfe0802c472d4dfd81ccffb49339107f3` |

The actual goldens remain on the benchmark host. `manifest.tsv` records the
input/index identities and host environment; `all-runs.tsv` and `summary.tsv`
contain the compact timing, CPU, RSS, exactness, size, and hash evidence.

## Resume repeats 2–5

The release driver stopped after the first-repeat gate by request:

```bash
GATE_ONLY=1 ./benchmarks/run-v0121-full-release-matrix.sh
```

It can resume without repeating any completed run:

```bash
GATE_ONLY=0 ./benchmarks/run-v0121-full-release-matrix.sh
```

The driver reruns Original for repeats 2–5, runs every candidate configuration
four more times, compares every output, and then replaces the single-repeat
summary with five-run means.
