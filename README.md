# vcftools-ng

Experimental high-performance, output-compatible successor to VCFtools 0.1.17.

The current implementation uses HTSlib for VCF/BCF decoding and a bounded,
order-preserving pipeline. A reader, persistent compute workers, and the
ordered committer overlap work on up to three batches. Each worker reuses its
decode scratch space; recoded VCF records are formatted in workers and
committed in exact input order.

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

VCFtools 0.1.17 is the compatibility oracle. Before the implementation is
considered mature, correctness and performance are evaluated only on a real
2,300,000-record subset: 100,000 records from each of 23 chromosomes. Every
supported output must pass `cmp` against the original program. The complete
11,230,392-record dataset is intentionally reserved for a later milestone.

The current compatibility gate covers:

- Outputs: `--freq`, `--counts`, `--missing-site`, `--site-depth`,
  `--site-mean-depth`, `--recode`, `--recode-INFO-all`
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

BCF is currently the preferred high-performance input. HTSlib's textual VCF
record parsing remains largely serial, so adding compute threads does not
materially accelerate BGZF VCF input in this first implementation.
For a single `--chr` selection, BCF plus CSI uses indexed iteration; optional
`--from-bp` and `--to-bp` bounds are pushed into the index query and then
rechecked by the compatibility filter.

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

Run the non-reference AF/AC benchmark with:

```bash
./benchmarks/run-non-ref-subset.sh
```

## Benchmark dataset

See [data/README.md](data/README.md) for provenance, formats, validation, and
conversion timings. The 2.3-million-record compatibility fixture is documented
in [tests/fixtures/README.md](tests/fixtures/README.md).
