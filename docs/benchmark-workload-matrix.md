# Representative workload benchmark matrix

[简体中文](benchmark-workload-matrix.zh-CN.md)

Performance claims are workload-specific. A speedup measured for one output,
filter set, input format, storage device, or selected-sample density must not
be presented as a guarantee for every vcftools-ng command.

## v0.14.2 scope

The v0.13.0 direct text kernel is not limited to statistics. For eligible
Plain or BGZF VCF input it fuses the common seven-filter workload with
site-local statistics and VCF recode. The eligible filters are
`--min-alleles`, `--max-alleles`, `--minQ`, `--minGQ`, `--min-meanDP`,
`--max-missing`, and `--maf`. Eligible site outputs are `--freq`, `--freq2`,
`--counts`, `--missing-site`, `--site-depth`, `--site-mean-depth`, and
`--site-quality`.

v0.13.1 additionally routes eligible Plain VCF position include/exclude,
sample projection, window pi, Tajima's D, site FST, and window FST through the
direct text family. Those paths share selected-sample GT decoding, compact
per-site contributions, deterministic ordered reduction, and adaptive
read-only mapped ranges. v0.13.2 extends that family to the ten-filter
production combination, FILTER/INFO/FT, and shared site pi. No-index BGZF gets
bounded decompression and compute overlap; eligible LD, exact PCA, and indexed
BCF discordance get dedicated kernels. Unsupported shapes still use the
general compatibility pipeline.

v0.14.1 compiles every invocation into one immutable capability plan, removes
unused FORMAT-field work from GT-only analyses, specializes deterministic LD
and PCA post-scan storage, and hardens malformed-input/failure/ploidy
boundaries. It incorporates the unreleased v0.13.2 work above.

v0.14.2 preserves those scientific semantics and corrects the portable and
storage layer: private libdeflate 1.25, format-aware strict BCF stream
planning, and aligned pread for Plain VCF above 8 GiB. Smaller Plain inputs
retain zero-copy mmap.

Every run log records `Execution kernel`, `Execution components`, input
backend, thread allocation, and high-level stage times. These fields must be
reported with benchmark results so a fast-path fallback is visible.

## Development matrix

The matrix covers computation shapes rather than every parameter Cartesian
product:

| ID | Workload | Required variants |
|---|---|---|
| W01 | Unfiltered site counts | `--counts` |
| W02 | Production seven-filter counts | exact production filter set |
| W03 | Site inclusion | approximately 1% and 50%; sorted and shuffled lists |
| W04 | Site exclusion | approximately 1% and 50%; duplicates and absent sites retained in the list fixture |
| W05 | Sample inclusion | approximately 25%, 50%, and 100% of samples with counts |
| W06 | Sample selection plus recode | 50% inclusion and production filters |
| W07 | Window pi | production overlapping and non-overlapping windows |
| W08 | Tajima's D | production window size |
| W09 | Site FST | smallest and largest real population pair |
| W10 | Window FST | a real pair, overlapping windows, plus one multiallelic case |
| W11 | Output-intensive recode | Plain VCF, BGZF VCF, and BCF output where supported |
| W12 | General multiallelic path | site statistics without the biallelic restriction |

The current Osmanthus production profile is 100 kb windows with a 10 kb step
for pi and window FST, and 100 kb for Tajima's D. This profile comes from the
current project configuration; changing it requires recording a new profile,
not silently changing the benchmark.

Development runs use 23,000 real records, threads `1 4 8 16 32`, and at most
three repeats. v0.13.1 locks a 230,000-record SSD/NVMe matrix for W03-W10;
v0.13.2 adds stable 230k A/B cases for its new kernels
so sub-second 23k startup noise is not presented as throughput. Larger
release-candidate runs use the standard 2,300,000-record real subset. The
stabilized local scaling set is `1 2 4 8 12 16 24 28 32`. v0.14.2 uses a
unified same-SSD portable A/B for BGZF+TBI, Plain VCF, and adaptive-stream
BCF. Application and post-`sync -f` durable time are recorded separately.

The reusable driver is
[`benchmarks/run-workload-matrix.sh`](../benchmarks/run-workload-matrix.sh).
Copy
[`benchmarks/workload-matrix-profile.example.sh`](../benchmarks/workload-matrix-profile.example.sh)
outside the repository and fill in the locked input, selection fixtures,
population files, oracle directory, and result directory. The driver validates
`ORACLE_ROOT/SHA256SUMS` before running and refuses to generate or replace an
Original oracle.

## Exactness and baselines

- Original VCFtools 0.1.17 goldens and timings are generated once, hashed,
  and retained. They are not regenerated during ordinary optimization.
- Each routine row runs once and is a complete byte/hash gate. Performance
  differences within 5% are treated as tied; this band never applies to
  scientific output bytes.
- Text outputs use `cmp`; BGZF output is decompressed and compared with the
  Original VCF text oracle; BCF uses the existing canonical compatibility
  procedure.
- Floating-point window and population outputs remain byte-exact in exact
  mode. A numerical tolerance is not a replacement for the gate.
- Results report run count, wall time, CPU, peak RSS, output bytes, kernel,
  backend, selected/total samples, input/kept sites, and produced rows or
  windows.

## Interpreting stage times

The built-in timings are deliberately high-level and low overhead:

- `input/index planning` includes format/index inspection and source setup;
- `pipeline setup` includes selectors, sample maps, and output construction;
- `ordered input/compute/commit` includes the general pipeline's overlapping
  read, decode, analysis, and ordered publication;
- `fused scan/filter/output` covers the direct text kernel;
- `output finalization` covers final flush, close, and writer validation.

Fine-grained parse, reducer, compression, and commit-wait attribution belongs
in dedicated profiling builds. Per-record timers must not be enabled in
release benchmarks because their overhead would distort the workload.

The committed v0.13.1 matrix, compact timings, oracle/input hashes, and exact
runner are in
[`benchmarks/results/workload-matrix-230k-v0130/RESULTS.md`](../benchmarks/results/workload-matrix-230k-v0130/RESULTS.md).
The v0.13.2 nine-family exact gate, oracle hashes, and A/B summary are in
[`benchmarks/results/v0132-development-gate/README.md`](../benchmarks/results/v0132-development-gate/README.md).
The v0.14.2 unified portable driver and result are
[`benchmarks/run-v0142-unified-full-ab.sh`](../benchmarks/run-v0142-unified-full-ab.sh)
and
[`benchmarks/results/full-unified-v0142-ab/`](../benchmarks/results/full-unified-v0142-ab/README.md).
