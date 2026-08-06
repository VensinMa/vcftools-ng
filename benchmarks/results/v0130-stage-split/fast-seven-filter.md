# v0.13.0 seven-filter fused-path development check

This is a development measurement, not a release benchmark. It uses the
230,000-site real-data fixture in
`/tmp/vcftools-ng-stage-split-230k` (10,000 sites from each of 23
chromosomes) and the high-frequency filter workload:

```text
--min-alleles 2 --max-alleles 2
--minGQ 10 --minQ 30 --min-meanDP 7
--max-missing 0.9 --maf 0.1 --counts
```

Every result matched SHA-256
`c6aa120a94e37bec0cd93215f12349d82bcb5bb2a2b294e65299c1def00fc0ba`.
The kept-site count was 118,305 of 230,000.

| Scenario | Threads | v0.12.4 wall (s) | v0.13.0 wall (s) | Speedup | v0.13.0 RSS (KiB) | Backend |
|---|---:|---:|---:|---:|---:|---|
| Plain VCF | 4 | 1.69 | 0.91 | 1.86x | 329684 | fast-counts-plain |
| Plain VCF | 8 | 0.96 | 0.55 | 1.75x | 383272 | fast-counts-plain |
| Plain VCF | 12 | 0.93 | 0.45 | 2.07x | 315436 | fast-counts-plain |
| Plain VCF | 16 | 0.78 | 0.46 | 1.70x | 325804 | fast-counts-plain |
| Plain VCF | 24 | 0.58 | 0.33 | 1.76x | 320968 | fast-counts-plain |
| Plain VCF | 28 | 0.57 | 0.32 | 1.78x | 307036 | fast-counts-plain |
| Plain VCF | 32 | 0.55 | 0.30 | 1.83x | 321324 | fast-counts-plain |
| BGZF VCF + TBI | 4 | 7.68 | 1.77 | 4.34x | 9048 | fast-counts-indexed-bgzf |
| BGZF VCF + TBI | 8 | 7.61 | 1.16 | 6.56x | 10068 | fast-counts-indexed-bgzf |
| BGZF VCF + TBI | 12 | 7.67 | 1.15 | 6.67x | 11628 | fast-counts-indexed-bgzf |
| BGZF VCF + TBI | 16 | 7.69 | 1.12 | 6.87x | 11372 | fast-counts-indexed-bgzf |
| BGZF VCF + TBI | 24 | 7.77 | 1.13 | 6.88x | 11608 | fast-counts-indexed-bgzf |
| BGZF VCF + TBI | 28 | 7.66 | 0.95 | 8.06x | 12452 | fast-counts-indexed-bgzf |
| BGZF VCF + TBI | 32 | 7.43 | 0.98 | 7.58x | 13744 | fast-counts-indexed-bgzf |
| BCF stream | 4 | 1.12 | 0.59 | 1.90x | 62260 | stream |
| BCF stream | 8 | 0.58 | 0.38 | 1.53x | 64216 | stream |
| BCF stream | 12 | 0.41 | 0.35 | 1.17x | 50644 | stream |
| BCF stream | 16 | 0.34 | 0.32 | 1.06x | 66712 | stream |
| BCF stream | 24 | 0.31 | 0.37 | 0.84x | 53188 | stream |
| BCF stream | 28 | 0.33 | 0.34 | 0.97x | 57740 | stream |
| BCF stream | 32 | 0.35 | 0.33 | 1.06x | 55208 | stream |

The v0.12.4 values are the locked direct-comparison baseline in
`comparison-v0124.tsv`. Each v0.13.0 value above is one warm-cache run, so
the table is suitable for direction finding only.
The 230,000-site BCF runs finish in roughly 0.3 seconds; their 24/28-thread
ratios are dominated by startup and measurement noise and are not used to
change the BCF production policy.

Rejected experiments:

- letting the pipeline reader replace one compute worker did not improve
  utilisation because it can block inside ordered `next_batch()`;
- asking input workers to steal one compute slice after every record chunk
  interrupted the text-parser hot path and regressed wall time;
- increasing in-flight batches from 3 to 8 had no stable benefit;
- increasing the shard look-ahead window from 1x to 4x input workers slightly
  helped Plain VCF but collapsed BGZF CPU utilisation to about 200%, so it was
  reverted.

The retained change extends the direct text-statistics kernel to the seven
filters above. It uses worker-local shard and output buffers deliberately:
the higher Plain-VCF RSS is accepted because it more than halves low-thread
wall time while preserving exact bytes.
