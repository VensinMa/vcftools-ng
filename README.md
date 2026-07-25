# vcftools-ng

Experimental high-performance, output-compatible successor to VCFtools 0.1.17.

The current implementation uses adaptive ordered input shards and a bounded,
order-preserving pipeline. Plain VCF uses aligned byte ranges; BGZF VCF plus
TBI/CSI and BCF plus CSI use independent indexed regions; inputs without a
sidecar are automatically indexed with bcftools when possible, and inputs
without a usable direct backend fall back to an HTSlib stream. Persistent compute
workers and the ordered committer overlap work on up to three batches. Each
worker reuses its decode scratch space; recoded records are committed in
exact input order.

## Version records

Every development version has a local record of supported parameters,
compatibility fixtures, retained-record counts, byte-comparison results,
original/8-thread/16-thread benchmarks, speedups, CPU, RSS, and reproduction
commands where those measurements were recorded:

- [Version history](docs/VERSION_HISTORY.md)
- [Per-version archive and update policy](docs/versions/README.md)
- [Machine-readable benchmark table](docs/versions/benchmarks.tsv)
- [New-version record template](docs/versions/TEMPLATE.md)

## Compatibility contract

VCFtools 0.1.17 is the compatibility oracle. Development correctness and
performance are evaluated on a real 2,300,000-record subset: 100,000 records
from each of 23 chromosomes. Every supported output must pass `cmp` against
the original program. The v0.11.2 final input-backend matrix additionally ran
the complete 11,230,392-record dataset in seven input/index scenarios, with
original plus 1/2/4/8/16/32-thread vcftools-ng runs repeated five times. All
245 outputs were byte-identical.

The current compatibility gate covers:

- Outputs: `--freq`, `--counts`, `--missing-site`, `--site-depth`,
  `--site-mean-depth`, `--depth`, `--missing-indv`, `--het`, `--hardy`,
  `--site-quality`, `--site-pi`, `--window-pi`, `--TajimaD`,
  `--weir-fst-pop`, `--geno-r2`, `--pca`, `--pca-no-norm`,
  `--recode`, `--recode-bcf`, `--recode-vcf-gz`,
  `--recode-INFO-all`
- Two-file outputs: `--diff-site`, `--diff-indv`,
  `--diff-site-discordance`, `--diff-indv-discordance`
- Site filters: `--min-alleles`, `--max-alleles`, `--remove-indels`,
  `--keep-only-indels`, `--minQ`, `--min-meanDP`, `--max-meanDP`,
  `--max-missing`, `--max-missing-count`, `--maf`, `--max-maf`,
  `--mac`, `--max-mac`, `--hwe`,
  `--non-ref-af`, `--max-non-ref-af`, `--non-ref-af-any`,
  `--max-non-ref-af-any`, `--non-ref-ac`, `--max-non-ref-ac`,
  `--non-ref-ac-any`, `--max-non-ref-ac-any`,
  `--chr`, `--not-chr`, `--from-bp`, `--to-bp`, `--positions`,
  `--exclude-positions`, `--bed`, `--exclude-bed`, `--thin`
- Site flag filters: `--keep-filtered`, `--remove-filtered`,
  `--remove-filtered-all`, `--keep-INFO`, `--remove-INFO`
- Genotype filters: `--minGQ`, `--minDP`, `--maxDP`,
  `--remove-filtered-geno`, `--remove-filtered-geno-all`
- Sample filters: `--keep`, `--remove`, `--indv`, `--remove-indv`
- Inputs: VCF, BGZF VCF, and BCF

This is not yet a complete replacement for every VCFtools option. Unsupported
options must not be treated as compatible until they have their own
differential tests.

Exact raw `--recode-bcf` currently requires BCF input and does not accept
sample subsetting or genotype masking; those cases lack a valid VCFtools
0.1.17 BCF oracle. `--recode-vcf-gz` is validated by byte-comparing its
decompressed stream with original `--recode`. Current diff filtering is
limited to chromosome/position/BED and sample selection.

## Build

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DHTSLIB_ROOT=/home/vensin/anaconda3/pkgs/htslib-1.23.1-h633afcb_0
cmake --build build -j
```

## Run

```bash
./build/vcftools-ng \
  --bcf tests/fixtures/osmanthus412.23chr_100k.bcf \
  --threads 16 \
  --min-alleles 2 --max-alleles 2 --remove-indels \
  --minQ 40 --minGQ 20 --minDP 5 --maxDP 30 \
  --min-meanDP 10 --max-missing 0.9 --maf 0.1 \
  --recode --recode-INFO-all \
  --out results/subset
```

BCF plus CSI is currently the fastest input. BGZF VCF plus TBI/CSI and plain
VCF also use parallel input adapters and continue scaling from 8 to 16
threads. A local BGZF VCF or BCF without an index now runs
`bcftools index --csi --threads N` automatically and then selects the indexed
adapter. `N` is the explicit vcftools-ng thread count or the automatically
detected scheduler/affinity count. Use `--no-auto-index` to retain the
streaming behavior, `--input-backend stream` to force it, or
`--bcftools FILE` to select the executable. Optional `--chr`, `--from-bp`,
and `--to-bp` selections are pushed into indexed shards and rechecked by the
compatibility filter.

Plain VCF cannot use CSI/TBI because those formats store BGZF virtual
offsets. It remains on the parallel aligned-byte-range adapter and never
invokes automatic indexing. Existing `.csi` and `.tbi` files are loaded and
validated independently. A valid sidecar is always preserved; if only
unusable or stale sidecars exist, automatic mode warns and falls back instead
of overwriting user files.

The input architecture, remaining ordered-BGZF work, 8/16-thread development
gate, and final 32-256-thread scaling matrix are documented here:

- [Adaptive input backend design](docs/architecture/adaptive-input-backends.md)
- [Input backend benchmark driver](benchmarks/run-input-backend-matrix.sh)
- [Final full-dataset benchmark driver](benchmarks/run-final-full-matrix.sh)
- [Final full-dataset summary](benchmarks/results/final-full-v0112/summary.tsv)

## Verify

```bash
ctest --test-dir build --output-on-failure
```

The differential test performs unfiltered and combined-filter comparisons on
all 2.3 million subset records, then verifies single-thread and multithread
determinism.

Run the standard 8/16-thread statistics benchmark with:

```bash
./benchmarks/run-subset.sh
```

Run the filtered VCF recode benchmark with:

```bash
./benchmarks/run-recode-subset.sh
```

Run the sample-selection plus six-output benchmark with:

```bash
./benchmarks/run-sample-subset.sh
```

Run the MAC/HWE filtering benchmark with:

```bash
./benchmarks/run-genetics-filter-subset.sh
```

Run the BED interval-filter benchmark with:

```bash
./benchmarks/run-bed-subset.sh
```

Run the ordered thinning benchmark with:

```bash
./benchmarks/run-thin-subset.sh
```

Run the v0.10 individual statistics benchmark with:

```bash
./benchmarks/run-v010-statistics.sh
```

Run the v0.10 population/LD/PCA/conversion/diff benchmark with:

```bash
./benchmarks/run-v010-advanced.sh
```

Run the non-reference AF/AC benchmark with:

```bash
./benchmarks/run-non-ref-subset.sh
```

## Benchmark dataset

See [data/README.md](data/README.md) for provenance, formats, validation, and
conversion timings. The 2.3-million-record compatibility fixture is documented
in [tests/fixtures/README.md](tests/fixtures/README.md).
