# vcftools-ng vX.Y.Z

## Scope

Describe the implementation slice and architecture changes.

## Parameters added

- `--parameter`

## Exact semantics

List ordering rules, original quirks, missing-value behavior, and interactions
with existing filters.

## Validation

- Oracle: VCFtools 0.1.17
- Fixture:
- Input format:
- Samples:
- Records before/after:
- Outputs compared:
- `cmp` result:
- Golden files/SHA:
- CTest result:

## Benchmark

| Workload | Original | 8 threads | 16 threads | Speedup 8/16 | CPU 8/16 | RSS 8/16 | Exact |
|---|---:|---:|---:|---:|---:|---:|---|
| | | | | | | | |

## Reproduce

```bash
# Original VCFtools 0.1.17

# vcftools-ng 8 threads

# vcftools-ng 16 threads

# byte comparison
```

## Known limitations

- None recorded.

