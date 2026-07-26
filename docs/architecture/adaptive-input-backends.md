# Adaptive input backends and scalable ordered shards

Status: Plain-range and indexed-region adapters implemented in v0.11.0;
automatic CSI construction implemented in v0.11.1; protected explicit index
validation implemented in v0.11.2; direct text-to-count fusion implemented
in v0.11.3 for Plain VCF and indexed BGZF VCF. Ordered BGZF blocks remain the
next input-engine milestone.

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
| BCF + CSI | indexed regions | ordered coordinate span | decompression and BCF decode parallel |
| Local BGZF VCF/BCF without index | build CSI, then indexed regions | index construction followed by ordered coordinate spans | first run includes indexing; subsequent runs use parallel regions |
| Ordinary gzip VCF | streaming fallback | bounded record batch | decompressor limited |

The backend is an input concern. Filtering, shared GT/DP decoding, analysis
payloads, ordered thinning, statistics, and writers continue to consume one
common ordered record stream.

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
- a total CPU budget split between decompression and parsing;
- several tasks per worker for load balancing;
- bounded in-flight work and ordered result publication;
- an automatic fallback when an index is stale or a direct backend fails.

The Python process pool, one `tabix` subprocess per region, and per-region
temporary matrix files are implementation choices, not designs to copy into
the C++ engine.

## Current vcftools-ng bottleneck

v0.10 opens one HTSlib stream, gives HTSlib at most four I/O threads, and has
one thread repeatedly call `bcf_read`. Compute slices can overlap that reader,
but VCF text parsing and record duplication still pass through the single
reader. Full-file BCF also uses this path; only a single selected BCF
chromosome can currently use an index iterator.

Consequently:

- BGZF VCF remains dominated by serial text parsing;
- BCF and inexpensive statistics flatten when the reader or ordered committer
  becomes dominant;
- adding compute workers cannot repair either bottleneck;
- `std::thread::hardware_concurrency()` alone is not a sufficient default on
  a scheduler-managed server.

## Common ordered-shard contract

Every direct backend produces `RecordShard` objects:

```text
RecordShard
  ordinal              stable total order, starting at zero
  source span          byte range, BGZF block range, or genomic region
  owned interval       the exact record-start interval owned by this shard
  records              records in original order
  first/last key       contig rank, POS, and within-position ordinal
  boundary metadata    data needed by thin/window/LD reducers
```

The coordinator may execute shards out of order, but exposes them to ordered
consumers by `ordinal`. Memory is bounded by a byte budget, not only a fixed
number of batches. A slow early shard applies backpressure instead of allowing
hundreds of later shards to accumulate in RAM.

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

`--threads` is a total CPU budget. Automatic selection checks, in order:

1. an explicit user value;
2. scheduler allocation such as `SLURM_CPUS_PER_TASK`, `PBS_NP`, `NSLOTS`,
   or `LSB_DJOB_NUMPROC`;
3. process CPU affinity;
4. online logical CPUs.

The planner assigns that budget to reader/decompress, parse/compute, and
compression lanes. It must not hard-cap total useful CPUs at 16 or 32.
Concurrency is nevertheless bounded by measured storage throughput, shard
count, memory budget, file-descriptor limit, and output bandwidth. Reporting
both requested and effective stage concurrency makes such limits visible.

Guidelines:

- create at least 4-8 ready shards per active worker, subject to a minimum
  useful byte/record span;
- allow more than 4,096 shards for very large files, but derive the cap from
  file size and memory rather than worker count alone;
- reuse persistent HTSlib reader contexts instead of opening a file per task;
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
7. Implement ordered BGZF block decompression and VCF/BCF record framing for
   no-index input.
8. Add cross-shard protocols for thinning and window statistics.
9. Add LD/PCA/diff shard consumers.
10. Profile and replace the single ordered committer with ordered segment
   publication where exactness permits.

For every step, VCFtools 0.1.17 remains the oracle. During development the
real 2.3-million-record subset is compared byte for byte. Routine local
scaling tests stop at the machine's 32 available CPUs. The complete
11.23-million-record benchmark is reserved for a later final-stage gate.

## Benchmark matrix

Required input cases:

- BGZF VCF with TBI;
- the same BGZF bytes through a path with no index;
- uncompressed VCF;
- BCF with CSI;
- the same BCF bytes through a path with no index.

Each report contains original VCFtools 0.1.17, vcftools-ng at every requested
thread count, wall/user/system time, CPU percentage, maximum RSS, speedup,
record count, byte-comparison status, selected backend, effective stage
concurrency, storage type, cold/warm cache state, and repetitions.

The development default remains 8 and 16 threads. Release-candidate local
scaling uses `1,2,4,8,16,32`, bounded by the process CPU affinity.
