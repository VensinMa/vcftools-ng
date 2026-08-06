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
goldens. The lock records input bytes and SHA-256, relevant sidecar bytes and
SHA-256, golden bytes and SHA-256, Original wall time, workload, record count,
and cache policy.

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
3. BCF adaptive streaming full-scan path.

The BCF row is a speed-oriented full-file scenario, not an indexed-region
scenario. The driver uses the default `auto` policy and asserts
`Input backend: stream` in the diagnostic log. A neighbouring CSI is validated,
reported, and preserved, but deliberately not used, so its presence cannot
silently change the selected backend. BCF indexed access stays in targeted
region-query tests and the release-only four-scenario matrix uses one adaptive
BCF row.

It uses 1/2/4/8/16/32 threads by default. Before any candidate starts, it
checks every input, the required BGZF index, and every retained golden against
the lock. The neighbouring BCF CSI is not required by the baseline lock, but
when present its validation cost is part of the measured candidate run.
Original is not executed. Every candidate output must pass
`cmp`; SHA-256, bytes, wall time, user/system CPU, CPU utilization, and RSS are
then recorded. Candidate VCFs are deleted after the exact comparison, because
retaining every 12 GB duplicate would make iterative optimization impractical.
The three actual Original VCF goldens are never deleted.

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

The four scenarios are not part of daily development:

- BGZF VCF + TBI;
- BGZF VCF + automatic CSI;
- Plain VCF;
- BCF with the default adaptive policy.

They run only after the user explicitly decides that a performance change is
large enough to prepare a release. Release qualification reuses a retained
Original oracle without rerunning it when the Original version, input hashes,
workload, and compatibility contract are unchanged. A changed identity
requires one new Original run per input format. The driver uses independent
no-sidecar paths for every automatic-index run, validates the first candidate
repeat byte for byte, and then performs the agreed repeats. It updates README,
version history, environment report, portable Linux x86_64 archive, master,
and GitHub Release together.

The v0.13.0 storage matrix uses at most three candidate repeats. A row whose
first application-wall time exceeds 1,800 seconds runs only once. Otherwise,
every row runs twice after its first-repeat exactness gate. If either of those
two runs exceeds 600 seconds and their application-wall times differ by less
than 10% (symmetric percentage difference), repeat three is skipped and the
decision is retained as a machine-readable gate record. All other rows run
three times. Summary tables expose the actual repeat count; a one- or two-run
row is never reported as a three-run mean.

Automatic-index timing includes inspection, any CSI construction selected by
the adaptive policy, filtering, and output. Existing valid or invalid user
sidecars are never overwritten. The BCF row records the automatically selected
backend rather than splitting the same format into index-policy scenarios.

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
