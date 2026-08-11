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
| [v0.10.0](versions/v0.10.0.md) | Individual/population statistics, LD/PCA, formats, diff | 335.08 | 8.63 | 8.52 | PASS |
| [v0.11.0](versions/v0.11.0.md) | Adaptive Plain/TBI/CSI ordered input shards | 68.23 | 16.06 | 11.25 | PASS |
| [v0.11.1](versions/v0.11.1.md) | Automatic CSI construction, including first-run cost | 68.22 | 27.84 | 22.69 | PASS |
| [v0.11.2](versions/v0.11.2.md) | Protected CSI/TBI validation, including first-run cost | 68.22 | 27.82 | 23.17 | PASS |
| [v0.11.3](versions/v0.11.3.md) | Adaptive direct text-to-count fusion | 296.69 | 21.46 | 12.69 | PASS |
| [v0.12.1](versions/v0.12.1.md) | Fused site statistics and scalable exact recode | 2267.88 | 77.83 | 53.01 | PASS |
| [v0.12.2](versions/v0.12.2.md) | Workload-adaptive indexing and exact recode scaling | 2267.88 | 73.93 | 52.05 | PASS |
| [v0.12.3](versions/v0.12.3.md) | Comprehensive colored terminal help | inherited | inherited | inherited | PASS |
| [v0.12.4](versions/v0.12.4.md) | Standard reproducible run logging | 392.73 | 15.72 | 10.63 | PASS |
| [v0.13.0](versions/v0.13.0.md) | Transactional BGZF output and hot-path acceleration | 2267.88 | 44.09 | 28.50 | PASS (first-repeat release gate) |
| [v0.13.1](versions/v0.13.1.md) | Adaptive zero-copy selection and population analytics | 4.16 | 0.38 | 0.32 | PASS (230k locked matrix) |
| [v0.13.2](versions/v0.13.2.md) | Exact fused filters, no-index BGZF, LD/PCA/diff kernels | inherited | 0.43 | 0.33 | PASS (23k exact, 230k A/B) |
| [v0.14.1](versions/v0.14.1.md) | Capability-planned exact analytics and hardened boundaries | 2267.88 | 44.36 | 34.06 | PASS (11.23m first-repeat release gate) |
| [v0.14.2](versions/v0.14.2.md) | Portable libdeflate, BCF-aware planning, and oversized-Plain pread | 2267.88 | 40.72 | 31.68 | PASS (81-row portable A/B plus final Plain gate) |

The workload in each row is the version's representative compatibility
benchmark; workloads differ between rows. Consult the per-version page before
comparing versions directly. v0.9.0 also reran the six-output sample workload
at 3.11/3.13 seconds for 8/16 threads.

## Cumulative supported surface in v0.14.2

- Inputs: VCF, BGZF VCF, BCF.
- Outputs: `--freq`, `--freq2`, `--counts`, `--missing-site`, `--site-depth`,
  `--site-mean-depth`, `--depth`, `--missing-indv`, `--het`, `--hardy`,
  `--site-quality`, `--site-pi`, `--window-pi`, `--TajimaD`,
  `--weir-fst-pop`, `--geno-r2`, `--pca`, `--pca-no-norm`,
  `--recode`, `--recode-bcf`, `--recode-vcf-gz`,
  `--recode-INFO-all`.
- Two-file comparison: site/individual membership and site/individual
  discordance.
- Position/site selection: chromosome, closed range, position lists, BED,
  ordered thinning.
- Samples: file-based and inline inclusion/exclusion.
- Site/genotype filters: allele number/type, quality, depth, missingness,
  MAF/MAC/HWE, non-reference AF/AC, FILTER, INFO Flag, and FORMAT/FT.
- Execution: adaptive Plain-range, indexed TBI/CSI, and streaming input
  adapters feeding the bounded compute → ordered-commit pipeline; 8/16
  thread compatibility gates and deterministic parallel BGZF compression.
  Index use is workload-adaptive: one-thread BGZF and full-file BCF recode
  stream, multi-thread BGZF recode can reuse/build an index, and selective
  BGZF/BCF queries retain indexed acceleration.
- Logging: terminal diagnostics are mirrored to `PREFIX.log` by default;
  `--log-file` overrides the path and `--no-log-file` disables only the file.
  Logs include explicit index decisions, resource use, outputs, filters, and
  success/failure status without entering scientific output streams.
- Fast path: eligible `--freq`, `--freq2`, `--counts`, `--missing-site`,
  `--site-depth`, `--site-mean-depth`, and `--site-quality` workloads on
  Plain/BGZF VCF directly parse text and can share one ordered scan. Indexed
  BGZF uses ordered tabix windows; requested concurrency is capped by CPU
  affinity and the file-descriptor budget. Worker output has bounded
  backpressure and recode workers use independent HTSlib output headers.
  v0.13.1 additionally covers direct position include/exclude, sample
  projection, window pi, Tajima's D, and site/window FST on eligible Plain VCF
  workloads, with adaptive read-only mapped ranges at profitable concurrency.
  v0.13.2 adds the production ten-filter set, FILTER/INFO/FT, shared site pi,
  and bounded no-index BGZF overlap; LD, exact PCA, and indexed BCF diff gain
  dedicated exact kernels where their eligibility rules are satisfied.
  v0.14.1 compiles field requirements, eligibility, fallback reason, and
  generic decode into one immutable plan; LD/PCA post-scan storage is further
  specialized without changing scientific accumulation or output order.
  v0.14.2 retains that plan, restores libdeflate in the portable archive,
  gives BCF stream decoding a format-aware strict budget, and automatically
  uses aligned pread instead of input-sized mmap for Plain VCF above 8 GiB.

v0.14.2 completed a same-host, same-output portable A/B across v0.13.0,
v0.14.1, and v0.14.2: three inputs, nine thread counts, and 81/81 exact
initial outputs. Final post-fix Plain outputs also matched the locked oracle.
Performance differences within 5% are reported as tied; application and
durable time remain separate for the 57-59 GB output workload.

v0.14.1 passed a complete 11,230,392-record first-repeat release matrix:

- four input scenarios at `1/2/4/8/12/16/24/28/32` threads;
- 36/36 candidate outputs byte-identical;
- observed single-run speedup range 3.54×–73.32×;
- all inputs, indexes, and retained Original goldens size/SHA-256 validated;
- zero Original reruns; follow-up adaptive repeats remain pending.

v0.12.2 passed a full-data, seven-filter exact-recode matrix inherited by
v0.12.3, whose runtime execution paths are unchanged:

- 11,230,392 input records; 5,425,725 retained records;
- four input scenarios at 1/2/4/8/16/32 threads, five repeats each;
- 120/120 candidate outputs byte-identical;
- observed mean speedup range 5.97×–53.59×;
- v0.12.1 Original timings and goldens reused after hash validation; zero
  Original reruns in v0.12.2.

The adaptive BCF stream averaged 40.20 seconds at 32 threads and used at most
168,952 KiB RSS. `--no-auto-index` has been removed because automatic mode now
selects indexing only when the workload is expected to benefit.

v0.11.0 also measured exact counts at:

- Plain VCF: 36.15/11.11/8.23 seconds for original/8/16;
- BCF plus CSI: 43.89/5.15/3.32 seconds;
- BGZF VCF without an index: 68.22/53.01/53.12 seconds.

The indexed and Plain adapters continue scaling from 8 to 16 threads. The
unindexed BGZF stream remains the next input bottleneck.

The v0.10 representative row above is the five-output shared-scan benchmark.
The same release also measured exact full BCF conversion at
425.51/38.81/22.53 seconds and exact two-file site discordance at
136.14/8.97/8.99 seconds for original/8/16 threads.

Machine-readable values are in [benchmarks.tsv](versions/benchmarks.tsv).
