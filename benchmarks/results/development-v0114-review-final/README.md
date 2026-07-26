# Final three-scenario pre-release gate

This run validates the reviewed runtime immediately before freezing the
v0.12.1 release candidate.

- Workload: seven real-project filters plus
  `--recode --recode-INFO-all --stdout`
- Fixture: 2,300,000 records, 23 chromosomes, 412 samples
- Scenarios: BGZF VCF + TBI, Plain VCF, and BCF + CSI
- Threads: 1, 2, 4, 8, 16, and 32
- Oracle: locked VCFtools 0.1.17 baseline; Original was not rerun
- Result: all 18 outputs passed complete-file `cmp` and SHA-256 validation

This gate includes the bounded fused-output collector and per-worker HTSlib
output-header fix found during final review. The full normal CTest suite
passed 3/3 in 123.87 seconds; the full ASan/UBSan suite passed 3/3 in 303.03
seconds.

`summary.tsv` contains wall time, CPU, RSS, speedup, and exactness.
`all-runs.tsv` and `manifest.tsv` retain the compact run and oracle identities.
