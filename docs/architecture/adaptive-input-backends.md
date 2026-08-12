# Adaptive input backends and scalable ordered shards

Status: Plain-range and indexed-region adapters implemented in v0.11.0;
automatic CSI construction implemented in v0.11.1; protected explicit index
validation implemented in v0.11.2; direct text-to-count fusion implemented
in v0.11.3 for Plain VCF and indexed BGZF VCF; the same direct text kernel
extended to the seven-filter VCF-recode workload in v0.13.0. v0.14.1 shares
the strict stream-resource planner, failure propagation,
cancellation, reusable worker scratch, and parsed-record layout across the
eligible text paths. Ordered BGZF blocks remain a profile-gated future input
milestone rather than a release requirement.

Reference reviewed: `VensinMa/vcf2phylip` commit
`3a1d3ab24af2334bfa6bb61d5c3faaf45d5c4625`.

## Decision

`vcftools-ng` will not use one reader implementation for every file encoding.
It will inspect the real container and available index, then select an input
backend that exposes the greatest safe parallelism:

| Input | Preferred backend | Parallel unit | Expected scaling |
|---|---|---|---|
| Plain VCF | aligned byte ranges | complete-line byte span | storage and parse limited |
| BGZF VCF + TBI/CSI | indexed regions | ordered coordinate span | decompression and parse parallel |
| BCF full scan | ordered stream by default | bounded record batch | sequential BCF decode is usually fastest for full recode |
| Local BGZF VCF without index | adaptive stream or atomic CSI build | bounded batch or ordered coordinate span | depends on thread count and whether index cost is amortized |
| Selective BGZF VCF/BCF query | reuse or build a validated index | ordered coordinate span | avoids decoding unrelated regions |
| Ordinary gzip VCF | streaming fallback | bounded record batch | decompressor limited |

The backend is an input concern. General filtering, shared GT/DP decoding,
analysis payloads, ordered thinning, statistics, and writers consume one
common ordered-record contract. This is deliberately a semantic contract,
not one compulsory in-memory representation. Plain mmap/pread spans, indexed
region chunks, and streamed BGZF batches retain their storage-native ownership
so that abstraction does not introduce a copy into a queue of owned strings.
Eligible site-local Plain/BGZF workloads may fuse text parsing, filtering,
statistics, and VCF reconstruction inside the ordered shard worker;
unsupported options fall back before any output is published.

## What the vcf2phylip implementation demonstrates

The strong 8-to-16-thread scaling in two scenarios comes from removing the
single input feeder:

- Plain VCF is divided into proportional raw byte ranges. Every range is
  aligned to the next newline and read independently with `pread` when
  available. Parsing therefore scales without routing all text through one
  reader process.
- Indexed BGZF VCF is divided into contigs or coordinate windows. Each worker
  performs its own indexed query, decompression, and parsing. Results are
  returned in region order, and records are accepted only when their `POS`
  belongs to the worker's assigned window, preventing overlap duplicates.
- BGZF without an index uses multithreaded `bgzip` decompression but still
  sends one stdout stream through the parent chunker. That serial feeder
  explains why its 8- and 16-thread measurements are nearly identical.

Other transferable details are:

- detection from file bytes rather than filename suffix;
- scheduler and CPU-affinity-aware default worker discovery;
- a shared input/compute/I/O worker budget split between decompression and
  parsing;
- several tasks per worker for load balancing;
- bounded in-flight work and ordered result publication;
- an automatic fallback when an index is stale or a direct backend fails.

The Python process pool, one `tabix` subprocess per region, and per-region
temporary matrix files are implementation choices, not designs to copy into
the C++ engine.

## Historical bottleneck motivating the adapters

v0.10 opened one HTSlib stream, gave HTSlib at most four I/O threads, and had
one thread repeatedly calling `bcf_read`. Compute slices could overlap that
reader, but VCF text parsing and record duplication still passed through the
single reader. Full-file BCF also used this path; only a single selected BCF
chromosome could use an index iterator.

The later direct text adapters removed this ceiling for eligible Plain and
indexed BGZF workloads. The same pattern remains relevant when profiling a
fallback that still uses the generic stream. Historically, the symptoms were:

- BGZF VCF remains dominated by serial text parsing;
- BCF and inexpensive statistics flatten when the reader or ordered committer
  becomes dominant;
- adding compute workers cannot repair either bottleneck;
- `std::thread::hardware_concurrency()` alone is not a sufficient default on
  a scheduler-managed server.

## Common ordered-shard contract

Every direct backend satisfies the following logical `RecordShard` contract:

```text
RecordShard
  ordinal              stable total order, starting at zero
  source span          byte range, BGZF block range, or genomic region
  owned interval       the exact record-start interval owned by this shard
  records              records in original order
  first/last key       contig rank, POS, and within-position ordinal
  boundary metadata    data needed by thin/window/LD reducers
```

The contract does not require every backend to materialize that illustrative
structure. Plain input can expose non-owning views into mmap/pread storage;
indexed input can incrementally publish bounded region chunks; streamed input
can transfer ownership of bounded batches. The coordinator may execute shards
out of order, but exposes them to ordered consumers by `ordinal`. Memory is
bounded by a byte budget, not only a fixed number of batches. A slow early
shard applies backpressure instead of allowing hundreds of later shards to
accumulate in RAM.

The shared seam is responsible for:

- stable shard and chunk ordinals;
- one total CPU budget for reader, HTSlib I/O, and compute workers;
- bounded admission/backpressure;
- cancellation that wakes every reader, worker, and ordered waiter;
- first-failure propagation after all owned threads are joined;
- common stage-concurrency and shard metrics.

Backend implementations remain responsible for storage-specific byte/record
ownership and lifetime. A future refactor must demonstrate that it preserves
the zero-copy Plain path and bounded indexed publication before replacing the
current implementations.

Field requirements, kernel eligibility, and fallback logging are compiled by
the shared capability plan described in
[`query-plan.md`](query-plan.md).

Exact mode must preserve:

1. header contig order and physical record order;
2. duplicate records and the order of records sharing one position;
3. exactly one owner for each record;
4. VCFtools 0.1.17 filtering and formatting;
5. original site order for floating-point accumulations whose printed bytes
   depend on addition order.

If an index is absent, v0.11.1 first attempts atomic automatic CSI
construction for a local BGZF VCF/BCF. If indexing is disabled or fails, or
an existing index is stale or inconsistent with the header, automatic mode
falls back without changing output semantics.

## Backend designs

### Plain VCF: aligned byte ranges

1. Read and parse the header once and record the first data byte.
2. Divide the remaining file into substantially more ranges than active
   reader workers.
3. Move every interior boundary to the first byte after the next newline.
4. Use `pread` or independent file descriptors so workers do not share a file
   offset.
5. Parse complete lines with a worker-local header clone and scratch state.
6. publish shards in range order.

The implementation must test CRLF, a final line without newline, a record
larger than the target shard, empty contigs, duplicate positions, and malformed
rows at a boundary.

### Indexed BGZF VCF and indexed BCF

Use HTSlib APIs directly rather than starting a `tabix` process for every
task. A bounded pool of reader contexts owns independent file handles,
indices, headers, iterators, and decode scratch.

Regions follow header contig order. Coordinate windows are density-aware:
contig length is an initial estimate, and index linear bins or a bounded
profile refine the split so a dense chromosome does not become a straggler.
Queries may return a spanning record in adjacent windows; a shard owns only
records whose start position is inside its half-open interval. Record-start
ownership removes duplicates without losing long alleles or `INFO/END`
records.

For full-file exact mode, indexed parallelism is enabled only after a
validation pass establishes that region concatenation reproduces physical
file order. Otherwise the ordered BGZF backend is used.

### Adaptive CSI/TBI policy

Index construction and index use are speed decisions rather than unconditional
input preparation. The default `auto` policy classifies format, selection,
workload, and effective thread count before inspecting or building a sidecar:

| Input/workload | Default decision |
|---|---|
| Plain VCF | aligned ranges; never build CSI/TBI |
| BGZF VCF full-file recode, 1 thread | ordered stream |
| BGZF VCF full-file recode, 2+ threads | reuse TBI/CSI or build CSI |
| BCF full-file recode | ordered stream, even if CSI exists |
| BGZF VCF/BCF selective `--chr`/coordinate query | reuse or build an index |
| compact full-scan statistics, 1–3 threads | stream |
| compact full-scan statistics, 4+ threads | reuse a valid index, but do not build a one-use index |

When construction is selected, bcftools receives the effective vcftools-ng
thread count, whether explicitly requested or discovered from the scheduler,
CPU affinity, or hardware concurrency. The index is written to a
process-unique temporary sidecar and atomically published without replacing
an existing `INPUT.csi`.

An advisory input-file lock serializes vcftools-ng builders. Adaptive mode
falls back to the stream if bcftools or the target directory is unavailable;
forced indexed mode fails. `--input-backend stream` and
`--input-backend indexed` remain explicit diagnostic overrides.

CSI/TBI cannot index Plain VCF because their chunks use BGZF virtual offsets.
Plain VCF therefore remains on aligned byte ranges and does not enter the
automatic-index stage.

When the adaptive decision could use an index, every conventional sidecar is
explicitly loaded and checked for format, age, header compatibility, index
statistics, and a random-access probe on each non-empty contig. CSI is
preferred, but every reader receives the exact selected path, so an invalid
CSI cannot hide a valid TBI. Existing invalid or stale sidecars are protected:
adaptive mode falls back and never overwrites them. If the selected fast path
is streaming, unrelated sidecars are left untouched and are not part of the
critical path.

### BGZF without an index

Merely increasing `hts_set_threads` is insufficient because `bcf_read`
continues to parse records serially. The scalable path has three ordered
stages:

1. scan BGZF headers and issue independent compressed block reads;
2. decompress blocks in parallel and publish them by block ordinal;
3. stitch record boundaries, then send complete VCF lines or framed BCF
   records to parse workers.

VCF lines and BCF records can cross BGZF block boundaries, so block
decompression and record ownership must be separate steps. This backend is
the no-index fallback for genuine BGZF. Ordinary gzip cannot use it because
its Deflate stream has no independent BGZF block boundaries.

The current ordered HTSlib fallback still has one record feeder. Once that
feeder is saturated, waking more parse workers reduces throughput. v0.14.1
therefore caps fused compressed-stream compute workers at
four for GT-only/recode work and six when DP, GQ, FT, frequency, or mean-depth
logic is active. The cap is intersected with the strict total-thread plan, so
requests of 1-8 threads are unchanged. On the locked 230k no-index BGZF, this
reduced 16/32-thread W01 wall time by 7.83%/7.16% and W02 by 6.98%/5.48%.
Indexed-region and Plain backends do not use this ceiling.

### Ordinary gzip

Use the fastest available streaming decompressor and overlap it with parsing,
but report that decompression is a serial ceiling. Recommend BGZF plus TBI/CSI
or BCF for repeat high-throughput work. Exact automatic detection must not
mistake a `.vcf.gz` suffix for BGZF.

## Ordered analyses across shards

| Consumer | Safe parallel work | Ordered/boundary rule |
|---|---|---|
| Frequency, missingness, depth, HWE, quality | site-local decode and formatting | concatenate by shard ordinal |
| Individual depth/missingness | integer shard totals | merge in ordinal order |
| Heterozygosity | site contributions | retain original floating addition order |
| `--thin` | independent filters | boundary state carries the last accepted position per chromosome |
| π, Tajima's D, FST windows | site contributions and window membership | ordered span reducer handles overlapping windows |
| LD | shard-local pairs plus right halo | left site owns output; concatenate left-site order |
| PCA | cell/site contributions | fixed site order per cell in exact mode |
| Recode | worker formatting/compression | publish record or compressed-block ordinal |
| Diff | two independently sharded ordered streams | deterministic sorted merge join |

v0.11.0 first enabled the adapters for site-local statistics and recode, then
enabled them for `--thin`, windows, LD, and PCA after the complete
2.3-million-record differential gate preserved every output byte. Diff keeps
its independent two-stream implementation.

## Scaling to hundreds of CPUs

`--threads` supplies the process-tree CPU budget. On Linux the process affinity
is reduced to at most that many CPUs before worker pools or automatic index
children are created. Automatic
selection intersects every applicable limit rather than trusting only the
first one found:

1. scheduler allocation such as `SLURM_CPUS_PER_TASK`, `PBS_NP`, `NSLOTS`,
   or `LSB_DJOB_NUMPROC`;
2. process CPU affinity;
3. cgroup v1/v2 CPU quota;
4. online logical CPUs.

When `--threads` is omitted, the smallest applicable limit is capped at 128.
An explicit value is not subject to that automatic 128 ceiling, but it is
still reduced to the detected allocation limit.

The planner shares that CPU capacity between input/decompression and
parse/compute lanes, then applies file-descriptor and storage constraints.
I/O and ordered-output pools may contain more waitable pthread objects than N
so a blocked stage does not leave CPUs idle, but the complete process tree can
execute on at most N allowed CPUs. `--threads` is therefore not a promise about
the exact instantaneous pthread count. Concurrency is also bounded by storage
throughput, shard count, memory budget, and output bandwidth. Reporting both
requested and effective stage concurrency makes such limits visible.

Guidelines:

- create at least 4-8 ready shards per active worker, subject to a minimum
  useful byte/record span;
- allow more than 4,096 shards for very large files, but derive the cap from
  file size and memory rather than worker count alone;
- reuse persistent HTSlib reader contexts instead of opening a file per task;
- keep one persistent read-only descriptor per Plain VCF range worker and use
  `pread` for every assigned aligned range;
- compile immutable FORMAT/output requirements once per run and keep them out
  of the per-record hot path;
- batch accepted VCF lines at the ordered-commit seam and reuse worker-local
  BGZF compression contexts without changing compressed bytes;
- use worker-local scratch and mostly lock-free or sharded queues;
- keep completed payloads compact and spill large ordered text/BCF segments
  when the memory budget is reached;
- make NUMA placement and reader count tunable on multi-socket machines;
- distinguish compute scaling from page-cache or storage scaling in reports.

Hundreds of requested threads are not evidence of hundreds of useful reader
threads. The benchmark must record CPU utilisation and storage behaviour, and
the planner should stop increasing reader concurrency when profiling shows no
throughput gain.

## Implementation sequence and gates

1. **Complete in v0.11.0:** format/index inspection,
   scheduler/affinity CPU discovery, backend reporting, and a force-backend
   test option.
2. **Complete in v0.11.0:** plain VCF byte-range shards.
3. **Complete in v0.11.0:** TBI/CSI region shards for BGZF VCF and BCF.
4. **Complete in v0.11.1:** atomic, concurrency-safe automatic CSI
   construction using the effective thread budget.
5. **Complete in v0.11.2:** explicit CSI/TBI validation, protected invalid
   sidecars, and Plain VCF no-index policy.
6. **Complete in v0.11.3 for unfiltered `--counts`:** direct aligned-range
   VCF parsing and direct ordered-tabix GT counting without intermediate
   `bcf1_t` records; CPU-affinity and file-descriptor constrained workers.
7. **Complete in the v0.13.0 candidate:** compiled execution plans, fused
   DP/filter sample scans, batch VCF ordered commit, persistent Plain VCF
   descriptors, reusable deterministic BGZF compression contexts, an
   input-heavy text-parser allocation within the strict total CPU budget,
   and the seven high-frequency filters in the direct Plain/BGZF
   site-statistics and VCF-recode kernel.
8. **Complete in the v0.13.0 candidate:** replace complete in-memory indexed
   `RecordShard` payloads with bounded incremental
   `(shard_ordinal, chunk_ordinal)` publication.
9. **Complete in v0.14.1:** one strict stream CPU planner that
   accounts for the reader and HTSlib workers; a true serial one-thread path;
   common cancellation/failure gates; one parsed fixed-column/sample layout;
   worker-local reusable scratch; and precompiled small-population roles.
10. **Design decision locked in v0.14.1:** unify the semantic
    ordered-shard seam without forcing storage-native zero-copy adapters into
    one owned-string queue.
11. **Complete in v0.14.3:** enforce `--threads N` as a dynamic Linux CPU
    affinity budget inherited by bcftools children, explicitly account for
    hidden HTSlib/bcftools queue coordinators, force vendor BLAS runtimes to
    one thread, and allow I/O-waiting input/output pools to overlap without
    restoring the v0.14.2 CPU oversubscription bug.
12. Build one immutable capability/query plan for field requirements, fused
    eligibility, fallback reason, and logging.
13. Profile before implementing ordered BGZF block framing, adaptive batches,
    analysis lanes, LD cache blocking, SIMD, or further ordered-commit work.

For every scientific-output step, VCFtools 0.1.17 remains the oracle. Routine
development starts with generated fixtures and the two real 23k fixtures.
The 230k input is used for stabilized performance A/B; the 2.3-million-record
subset is an 8/16-thread late exact gate. Local scaling stops at the machine's
32 available CPUs. The complete 11.23-million-record benchmark is reserved
for a later final-stage gate after explicit approval.

## Benchmark matrix

Primary performance input cases:

- BGZF VCF with TBI;
- the same BGZF bytes through a path with no index;
- uncompressed VCF;

BCF remains in compatibility and release gates, but is not the primary target
for text-parser optimization.

Each report contains original VCFtools 0.1.17, vcftools-ng at every requested
thread count, wall/user/system time, CPU percentage, maximum RSS, speedup,
record count, byte-comparison status, selected backend, effective stage
concurrency, storage type, cold/warm cache state, and repetitions.

Early performance A/B uses `1,8,16,32`; stabilized planner scaling uses
`1,2,4,8,12,16,24,28,32`, bounded by process CPU affinity. Original and its
golden artifacts are reused during development rather than rerun for every
candidate.
