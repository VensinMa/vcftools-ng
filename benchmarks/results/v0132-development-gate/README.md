# v0.13.2 development gate

This is the compact, reproducible record for the v0.13.2 development cycle.
It deliberately uses only the locked 23,000-record real-data fixtures for
correctness and the locked 230,000-record fixture for A/B performance work.
No 2.3-million or 11.23-million-record run is relabeled as v0.13.2 evidence.

## Exact gate

`benchmarks/run-v0132-development-gate.sh` checks 1/4/8/16/32 threads for:

- the real ten-filter counts workload;
- the shared exact HWE filter;
- site FILTER plus INFO Flag filtering and recode;
- genotype FT filtering with counts, missing-site, and recode in one scan;
- multi-output `--site-pi --counts` on the GATK-derived fixture;
- no-index BGZF decompression/compute overlap;
- bitset/popcount genotype LD;
- exact PCA SoA execution; and
- indexed BCF site-discordance diff.

All 45 executions passed byte comparison. `oracles/SHA256SUMS` locks eleven
scientific reference artifacts. Original VCFtools or v0.13.1 reference output
is generated only when `GENERATE_ORACLES=1`; normal development runs validate
the hashes and never rerun the slow reference.

The 23k wall values in `runs.tsv` are primarily correctness feedback and are
not used for headline performance claims because many executions are below
one second.

The separate [real-data GT audit](GT_COMPATIBILITY.md) enumerates every GT
token in both 23-chromosome fixtures. It confirms phased, partially missing,
haploid-missing, and multiallelic forms beyond the usual six examples and
documents the exact fast-path/fallback boundary.

The lightweight release test additionally forces generic and specialized
scalar parsers over the synthetic grammar and real DeepVariant/GATK layouts.
Five triploid spellings are transactionally rejected by GT-semantic commands,
while site-only `--minQ`/`--positions` VCF recode preserves them as opaque text
at 1 and 32 threads. This passthrough is a documented vcftools-ng extension;
Original 0.1.17 rejects even though the command does not need GT semantics.

After the final code and version update, the broader W01-W10 matrix also
passed **95/95** configurations (19 workloads x 1/4/8/16/32 threads, one final
repeat) against its locked oracles. Its compact record is retained at
`benchmarks/results/workload-matrix-23k-v0130/v0132-final-candidate-gate/runs.tsv`.

## 230k A/B findings

Stable A/B runs against the v0.13.1 code path established the following:

| Workload | Threads | v0.13.1 | v0.13.2 | Relative gain |
|---|---:|---:|---:|---:|
| Real ten-filter counts | 1 | 6.73 s | 3.21 s | 2.10x |
| Real ten-filter counts | 8 | 1.22 s | 0.43 s | 2.84x |
| Real ten-filter counts | 16 | 0.92 s | 0.33 s | 2.79x |
| Real ten-filter counts | 32 | 0.57 s | 0.20 s | 2.85x |
| Site FILTER counts | 1 | 6.70 s | 1.13 s | 5.93x |
| Site FILTER counts | 32 | 0.54 s | 0.11 s | 4.91x |
| Genotype FT missing | 1 | 6.73 s | 3.76 s | 1.79x |
| Genotype FT missing | 32 | 0.61 s | 0.23 s | 2.65x |
| No-index BGZF ten-filter counts | 1 | 7.60 s | 3.65 s | 2.08x |
| No-index BGZF ten-filter counts | 8 | 6.63 s | 0.72 s | 9.21x |
| No-index BGZF ten-filter counts | 32 | 6.66 s | 0.79 s | 8.43x |

Advanced reducer changes improved the locked 230k window pi, Tajima's D, and
window FST cases by approximately 1.12x-1.22x without changing floating-point
accumulation order. LTO added about 4%-7% at 16/32 threads in the measured
site-local workloads.

PGO was tested separately and is opt-in rather than silently required by the
normal build. Against the ordinary LTO build, its five-run means were:

| Workload | Threads | LTO | PGO | Relative gain |
|---|---:|---:|---:|---:|
| Counts | 1 | 2.358 s | 1.762 s | 1.34x |
| Counts | 8 | 0.318 s | 0.254 s | 1.25x |
| Counts | 16 | 0.240 s | 0.200 s | 1.20x |
| Counts | 32 | 0.156 s | 0.132 s | 1.18x |
| Ten-filter counts | 1 | 2.982 s | 2.388 s | 1.25x |
| Ten-filter counts | 8 | 0.406 s | 0.332 s | 1.22x |
| Ten-filter counts | 16 | 0.300 s | 0.258 s | 1.16x |
| Ten-filter counts | 32 | 0.188 s | 0.162 s | 1.16x |

The indexed BCF diff and pair-dense LD A/B summaries are retained in
`benchmarks/results/v0132-indexed-diff-ab/runs.tsv` and
`benchmarks/results/v0132-ld-bitset-ab/runs.tsv`. Detailed generated outputs
and stderr logs remain local and are intentionally ignored by Git.

## Safety gates

- Release/CTest excluding the pre-existing 2.3m differential test: 7/7 pass.
- Final W01-W10 23k candidate matrix: 95/95 pass.
- ASan+UBSan for the same seven tests: 7/7 pass.
- Directed ASan+UBSan coverage for the new filter, INFO/FT, BGZF, LD, PCA,
  diff, and multi-output paths: pass.
- Exact HWE filter: Original byte gate at 1/4/8/16/32 plus ASan+UBSan pass.
- INFO Flag filters: Original byte gate at 1/4/8/16/32; non-Flag INFO fields
  are rejected.
- Compressed-stage budget: 3 threads = 1 reader + 1 HTSlib I/O + 1 compute;
  4 threads = 1 reader + 2 HTSlib I/O + 1 compute.
