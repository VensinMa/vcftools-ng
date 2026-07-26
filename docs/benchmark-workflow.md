# Permanent optimization and benchmark workflow

This policy separates fast daily engineering from expensive release
qualification. It applies to v0.11.4 and every later optimization.

## 1. Lock the oracle once

The development workload is fixed to:

```text
--min-alleles 2 --max-alleles 2
--minGQ 10 --minQ 30
--min-meanDP 7 --max-missing 0.9 --maf 0.1
--recode --recode-INFO-all --stdout
```

The fixture is the real 2,300,000-record, 23-chromosome, 412-sample subset.
VCFtools 0.1.17 is run once for each distinct source format: BGZF VCF, Plain
VCF, and BCF. Its actual VCF outputs, not only checksums, are retained as
goldens. The lock records input bytes and SHA-256, TBI/CSI bytes and SHA-256,
golden bytes and SHA-256, Original wall time, workload, record count, and cache
policy.

An oracle may be regenerated only when the fixture, workload, Original
version, or compatibility contract intentionally changes. Regeneration is a
separate reviewed action; a development script must never silently rerun
Original or overwrite a golden.

The initial lock is
`benchmarks/results/v0114-real-filter-subset/baseline.lock.tsv`. Its three
actual golden VCFs occupy 34 GB and remain in the adjacent local `golden/`
directory. Git retains their identities and reproduction evidence, while the
files themselves remain excluded because of their size.

## 2. Daily correctness gate

Every code change first runs CTest and the small real-data SHA gates. A changed
recode/filter/input path then runs:

```bash
./benchmarks/run-development-gate.sh
```

This gate has exactly three scenarios:

1. BGZF VCF + valid TBI;
2. Plain VCF;
3. BCF + valid CSI.

It uses 1/2/4/8/16/32 threads by default. Before any candidate starts, it
checks every input, index, and retained golden against the lock. Original is
not executed. Every candidate output must pass `cmp`; SHA-256, bytes, wall
time, user/system CPU, CPU utilization, and RSS are then recorded. Candidate
VCFs are deleted after the exact comparison, because retaining every
12 GB duplicate would make iterative optimization impractical. The three
actual Original VCF goldens are never deleted.

For a focused experiment, `THREAD_LIST` may temporarily reduce the candidate
thread set. A change cannot be declared development-gate complete until all
six standard thread counts pass.

## 3. Development performance comparison

The locked Original time is the baseline for all later versions. If an oracle
was originally measured repeatedly and the runs show no material drift, the
lock stores the median (preferred) or mean and states which statistic was
used. The current first lock contains one measured Original run per format.

One candidate repeat is sufficient while tuning. A meaningful optimization
should improve the affected rows without regressing exactness, one/two-core
behavior, or peak memory beyond the documented tradeoff. CPU utilization is
reported rather than inferred from thread count.

## 4. Release qualification

The seven scenarios are not part of daily development:

- BGZF VCF + TBI;
- BGZF VCF + automatic CSI;
- BGZF VCF + `--no-auto-index`;
- Plain VCF;
- BCF + CSI;
- BCF + automatic CSI;
- BCF + `--no-auto-index`.

They run only after the user explicitly decides that a performance change is
large enough to prepare a release. Release qualification may rerun Original,
uses independent no-sidecar paths for every automatic-index run, validates the
first candidate repeat byte for byte, and then performs the agreed repeats.
It updates README, version history, environment report, portable Linux x86_64
archive, master, and GitHub Release together.

Automatic-index timing includes inspection, CSI construction, filtering, and
output. `--no-auto-index` rows measure the genuine sequential compressed
fallback. Existing valid or invalid user sidecars are never overwritten.

## 5. Storage and evidence rules

- Benchmark commands run serially; concurrent runs would contaminate CPU,
  cache, and disk measurements.
- Cache state is reported. The default is a warm/unspecified operating-system
  page cache; results must not claim cold-cache behavior unless caches were
  deliberately controlled.
- Linux rotational-storage detection and eligible page-cache prefetch are
  tested separately on an HDD. They do not replace the three standard SSD
  development scenarios.
- Compact TSV, hashes, logs, environment metadata, and the baseline lock are
  versioned. Very large inputs and actual goldens stay local but must remain
  present and hash-valid on the benchmark host.
- A failed `cmp`, missing golden, mismatched input/index hash, or changed
  workload stops the gate immediately.
