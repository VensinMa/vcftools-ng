# Capability query plan and shared decode

Status: implemented in v0.14.1.

## Decision

Every invocation is compiled once into an immutable `QueryPlan` before an
input backend or compute kernel is opened. The plan is the single source for:

- required GT, DP, GQ, FT, chromosome, allele-string, quality, and recode
  payloads;
- active frequency, genotype, and mean-depth filters;
- fused text-kernel eligibility;
- the first precise fallback reason written to stderr and `<out>.log`;
- the generic ordered pipeline's decode requirements;
- a compact capability mask copied into each fused text worker.

This replaces the former independent `can_use_fused_site_stats()` and generic
`compile_execution_plan()` rule sets. Backend-specific settings such as
population membership and output sinks remain in `FastSiteStatPlan`; they do
not independently decide which FORMAT fields to decode.

## Execution shape

```text
Options + detected input format
              |
              v
        immutable QueryPlan
          /             \
 eligible text VCF      general compatibility
          |                    |
 compact capability      shared HTSlib decode
 mask per worker          + BatchAnalysisPayload
          |                    |
 one record/sample scan   one record/sample scan
          \                    /
           deterministic ordered outputs
```

Both branches already feed every requested compatible analysis from the same
record and sample decode. Adding another queue of analysis lanes would copy
payloads and add synchronization without removing a scan, so v0.14.1 retains
the existing deterministic consumers. LD and PCA perform their required
post-scan pair/matrix computations from the shared dosage payload; that work
cannot be replaced by a site-local reducer.

## Exactness rules

- A kernel is selected before transactional scientific outputs are opened.
- Unsupported combinations fall back as a whole; a single invocation never
  publishes a mixture of fused and generic scientific files.
- The first rejected capability becomes the stable logged fallback reason.
- The compact fused capability mask is passed separately from
  `FastSiteStatPlan`, preserving the established hot-field layout.
- Worker results are committed by shard/batch ordinal, so output order and
  floating-point accumulation order remain deterministic.

## Performance evidence

Moving unused FORMAT-field selection into the plan removed DP/GQ/FT scans
from GT-only counts and FST workloads. On the locked 230k Plain VCF, a
ten-repeat interleaved A/B measured:

| Workload | Threads | Change |
|---|---:|---:|
| W01 counts | 8 / 16 / 32 | -21.93% / -26.98% / -17.54% |
| W10 window FST | 8 / 16 / 32 | -20.70% / -25.41% / -17.37% |
| W02 ten-filter counts | 1 | +2.03% |

W02 retains the complete specialized GT+DP+GQ parser; its measured change is
inside the 3% development threshold. A subsequent same-source PGO build
improved all tested W01/W02/recode/site-FST/window-FST combinations by
2.22-15.11% at 1, 8, 16, and 32 threads without changing output hashes.

Raw drivers:

- `benchmarks/run-v0141-candidate-230k.sh`
- `benchmarks/run-v0141-field-requirement-ab.sh`
- `benchmarks/build-pgo.sh`

PGO changes branch and code layout only. It does not enable `-march=native`,
fast-math, or a newer ISA baseline. Source builds keep reproducible opt-in PGO.
The public manylinux2014 portable archive uses Release+IPO/LTO because its
large real-data training fixture is intentionally not embedded in the source
archive; its portability and exactness are tested independently.

The same PGO binary was also checked on the locked 230k BGZF inputs. It
improved indexed W01/W02 by 6.43-12.67% / 7.26-11.17%. On the no-index
compressed stream it improved all cases except 32-thread W01, whose +1.44%
change is inside the 3% threshold; the remaining changes were -3.07% to
-12.57%.

## Deterministic post-scan consumers

The dosage payload remains shared by LD and PCA, but their post-scan storage
is specialized for the access pattern of each analysis:

- LD stores valid/heterozygous/homozygous-alt bit planes in one contiguous
  allocation per site. The pair loop keeps the established popcount and
  floating-point operation order.
- LD left sites are formatted in fixed 64-site blocks. Blocks are claimed in
  parallel and published by block ordinal, so pair and output order are
  unchanged while the 230k probe reduces about 230,000 stream/string slots to
  about 3,600.
- PCA precomputes normalized dosage in sample-major order. Every covariance
  cell still sums sites in input order before the unchanged legacy LAPACK
  ordering step.

The locked 230k W13 LD output passed 40/40 comparisons in each of the
allocation-layout and formatter A/B gates. The confirmed formatter result was
neutral at one thread and improved 8/16/32-thread wall time by
1.29%/2.13%/1.20%, with peak RSS lower by 3.8%-19.4%. The locked 64-sample PCA
output passed 40/40 comparisons; whole-run wall time improved by about 2% at
one and eight threads and was neutral at 16 and 32 threads. PCA finalization
itself fell from 0.251 to 0.115 seconds at one thread. These measurements are
local 230k development evidence, not full-data claims.

After these post-scan changes were frozen, PGO was retrained from the locked
23k mix. The final same-source Plain A/B passed 168/168 output gates. PGO
improved W01/W02 and large/site-window FST by about 6%-15%; recode, LD, and
PCA were either faster or inside the 3% noise threshold. The corresponding
BGZF+TBI and no-index BGZF W01/W02 A/B passed 96/96 gates and improved wall
time by about 3%-15%. Raw results are in `v0141-final-pgo-plain-ab/` and
`v0141-final-pgo-bgzf-ab/` under the locked 230k result root.
