# vcftools-ng v0.13.1 W03-W10 locked 230k performance matrix

[简体中文](RESULTS.zh-CN.md)

Date: 2026-08-06  
Status: v0.13.1 release evidence.  
Scope: W03, W04, W05, W07, W08, W09, and W10. W11/W12 were outside this optimization cycle.

## Exactness and reproducibility boundary

- Input: `osmanthus412.flags.23chr_10k.vcf`
- Records: 230,000; samples: 412; size: 3,201,584,833 bytes
- Input SHA-256: `e62f8f06617371229825bb4747777469ac173c56aa4bc03cec469b29bbee1be2`
- Original: locked one-run VCFtools 0.1.17 timings and permanent oracles; Original was not rerun
- vcftools-ng: three runs per case/thread; tables report median wall time
- Threads: 1, 4, 8, 16, and 32
- 230k gate: 225/225 PASS; every output byte-identical to its Original oracle
- 23k full matrix: 285/285 PASS, including the GATK multiallelic W10 case
- Benchmarked performance-candidate SHA-256: `fd23cd1b67ea85c5617e6adcd4f1b1b3b54b72829c6ef705da6af34702e64e17`
- Frozen pre-optimization candidate SHA-256: `732b4d9c2866d4302dd4db7a9f4d2ee10ea34b08a6b835280d7ed429984da59f`

Compact evidence is retained in `original-runs.tsv`, `baseline-current/runs.tsv`,
`final-mmap/runs.tsv`, `oracles/SHA256SUMS`, and
`../../run-v0130-w03-w10-230k.sh`.

## Median wall time (seconds)

| Workload | Original | 1 | 4 | 8 | 16 | 32 |
|---|---:|---:|---:|---:|---:|---:|
| W03 positions 1% shuffled | 0.98 | 0.26 | 0.14 | 0.13 | 0.13 | 0.12 |
| W03 positions 1% sorted | 1.00 | 0.26 | 0.15 | 0.14 | 0.13 | 0.12 |
| W03 positions 50% shuffled | 2.53 | 1.37 | 0.40 | 0.25 | 0.22 | 0.16 |
| W03 positions 50% sorted | 2.53 | 1.37 | 0.41 | 0.26 | 0.22 | 0.16 |
| W04 exclude 1% duplicates/absent | 4.00 | 2.46 | 0.73 | 0.39 | 0.31 | 0.20 |
| W04 exclude 50% duplicates/absent | 2.55 | 1.38 | 0.40 | 0.28 | 0.21 | 0.16 |
| W05 keep 100% counts | 4.00 | 2.47 | 0.67 | 0.41 | 0.31 | 0.19 |
| W05 keep 25% counts | 3.91 | 1.18 | 0.35 | 0.23 | 0.18 | 0.14 |
| W05 keep 50% counts | 3.94 | 1.61 | 0.45 | 0.27 | 0.22 | 0.15 |
| W07 window pi non-overlap | 4.06 | 2.70 | 0.68 | 0.41 | 0.31 | 0.20 |
| W07 window pi overlap | 4.16 | 2.68 | 0.66 | 0.38 | 0.32 | 0.20 |
| W08 Tajima's D 100 kb | 4.03 | 2.69 | 0.67 | 0.39 | 0.32 | 0.20 |
| W09 site FST large pair | 4.55 | 2.58 | 0.70 | 0.39 | 0.33 | 0.22 |
| W09 site FST small pair | 4.09 | 0.87 | 0.28 | 0.19 | 0.15 | 0.14 |
| W10 window FST biallelic | 4.55 | 2.60 | 0.71 | 0.41 | 0.33 | 0.21 |

## Speedup over Original VCFtools 0.1.17

| Workload | 1 | 4 | 8 | 16 | 32 |
|---|---:|---:|---:|---:|---:|
| W03 positions 1% shuffled | 3.77x | 7.00x | 7.54x | 7.54x | 8.17x |
| W03 positions 1% sorted | 3.85x | 6.67x | 7.14x | 7.69x | 8.33x |
| W03 positions 50% shuffled | 1.85x | 6.32x | 10.12x | 11.50x | 15.81x |
| W03 positions 50% sorted | 1.85x | 6.17x | 9.73x | 11.50x | 15.81x |
| W04 exclude 1% duplicates/absent | 1.63x | 5.48x | 10.26x | 12.90x | 20.00x |
| W04 exclude 50% duplicates/absent | 1.85x | 6.37x | 9.11x | 12.14x | 15.94x |
| W05 keep 100% counts | 1.62x | 5.97x | 9.76x | 12.90x | 21.05x |
| W05 keep 25% counts | 3.31x | 11.17x | 17.00x | 21.72x | 27.93x |
| W05 keep 50% counts | 2.45x | 8.76x | 14.59x | 17.91x | 26.27x |
| W07 window pi non-overlap | 1.50x | 5.97x | 9.90x | 13.10x | 20.30x |
| W07 window pi overlap | 1.55x | 6.30x | 10.95x | 13.00x | 20.80x |
| W08 Tajima's D 100 kb | 1.50x | 6.01x | 10.33x | 12.59x | 20.15x |
| W09 site FST large pair | 1.76x | 6.50x | 11.67x | 13.79x | 20.68x |
| W09 site FST small pair | 4.70x | 14.61x | 21.53x | 27.27x | 29.21x |
| W10 window FST biallelic | 1.75x | 6.41x | 11.10x | 13.79x | 21.67x |

## Implementation result

The release keeps direct positions/exclude-positions lookup, selected-sample GT
projection, compact ordered pi/Tajima/FST reductions, specialized common
diploid-GT parsing, a fixed two-population biallelic FST accumulator, and
adaptive zero-copy Plain VCF scanning.

For Plain VCF, multi-thread runs use a read-only mapping. One-thread runs use
the mapping only for I/O-dominated position or sample-projection workloads;
pure one-thread pi/Tajima paths retain `pread`. Relative to the frozen candidate,
the final 32-thread workloads improved by 31.2%-60.0%.

Linux peak RSS for mapped workloads approaches the accessed input size. These
pages are reclaimable file-backed mappings, not equivalent anonymous heap.
The 230k fixture reports approximately 3,056-3,064 MiB RSS for mapped paths and
approximately 133 MiB for one-thread pi/Tajima `pread` paths. Larger-data and
rotational-storage qualification must inspect major faults and memory pressure
alongside wall time.

The experimental thread-local flat multiallelic FST scratch was removed because
the 23k multiallelic case regressed by 1.65% and 4.80% at one and four threads;
its apparent high-thread improvement was indistinguishable from short-run noise.
