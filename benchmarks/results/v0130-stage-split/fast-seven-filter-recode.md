# v0.13.0 seven-filter direct-recode development check

This is a development measurement, not a release benchmark. It uses the
23,000-site, 412-sample flagged real-data fixture
`tests/fixtures/osmanthus412.flags.23chr_1k.vcf.gz` and the high-frequency
filter/recode workload:

```text
--min-alleles 2 --max-alleles 2
--minGQ 10 --minQ 30 --min-meanDP 7
--max-missing 0.9 --maf 0.1
--recode-vcf-gz --recode-INFO-all
```

All v0.12.4 and v0.13.0 files were byte-identical with SHA-256
`d8455a9642c220045d95705a74beb794faa9b6eda58ac22c792939f419236e14`.
The kept-site count was 11,761 of 23,000.

| Threads | v0.12.4 wall (s) | v0.13.0 wall (s) | Speedup | v0.13.0 CPU | v0.13.0 RSS (KiB) |
|---:|---:|---:|---:|---:|---:|
| 4 | 0.91 | 0.70 | 1.30x | 499% | 137376 |
| 8 | 0.56 | 0.41 | 1.37x | 977% | 140216 |
| 16 | 0.50 | 0.29 | 1.72x | 1540% | 248252 |
| 32 | 0.61 | 0.28 | 2.18x | 1773% | 264912 |

Each value is one warm-cache run on the 32-CPU development host, so it is
only a direction-finding comparison. No 2.3-million-site differential or
11.23-million-site release benchmark was run for this change.

The direct path parses GT/GQ/DP once, performs genotype masking and the seven
site filters in the same worker, reconstructs accepted VCF lines in shard
order, and feeds the existing deterministic BGZF compressor. The compatibility
gate also covers Plain VCF, stdout, INFO removal, simultaneous Plain/BGZF
outputs, and `--counts` combined with recode in one scan.
