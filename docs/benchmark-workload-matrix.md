# Representative workload benchmark matrix

[简体中文](benchmark-workload-matrix.zh-CN.md)

Performance claims are workload-specific. A speedup measured for one output,
filter set, input format, storage device, or selected-sample density must not
be presented as a guarantee for every vcftools-ng command.

## v0.13.0 scope

The v0.13.0 direct text kernel is not limited to statistics. For eligible
Plain or BGZF VCF input it fuses the common seven-filter workload with
site-local statistics and VCF recode. The eligible filters are
`--min-alleles`, `--max-alleles`, `--minQ`, `--minGQ`, `--min-meanDP`,
`--max-missing`, and `--maf`. Eligible site outputs are `--freq`, `--freq2`,
`--counts`, `--missing-site`, `--site-depth`, `--site-mean-depth`, and
`--site-quality`.

Sample/site selection, FILTER/INFO/FT selection, unsupported filters, BCF,
individual reductions, window statistics, FST, LD, PCA, and diff use the
general compatibility pipeline. They can benefit from shared input,
scheduling, and output improvements, but the fused-kernel speedup does not
automatically apply to them.

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
three repeats. Release-candidate runs use the standard 2,300,000-record real
subset, threads `1 2 4 8 16 32`, and at least three repeats. The
11,230,392-record final gate uses four representative workloads: W02, W06,
either W07 or W08, and W10.

## Exactness and baselines

- Original VCFtools 0.1.17 goldens and timings are generated once, hashed,
  and retained. They are not regenerated during ordinary optimization.
- Repeat one is the byte-for-byte gate. Later repeats are timing-only after
  that gate succeeds.
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
