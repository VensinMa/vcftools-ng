# v0.11.4 development gate at `770aeee`

This is the first formal result produced by the permanent three-scenario
development workflow. It used the locked VCFtools 0.1.17 oracle and did not
rerun Original.

- Input scenarios: BGZF VCF + TBI, Plain VCF, and BCF + CSI
- Threads: 1, 2, 4, 8, 16, and 32
- Workload: the seven real-project filters plus
  `--recode --recode-INFO-all --stdout`
- Correctness: all 18 candidate VCFs passed byte-for-byte `cmp` and SHA-256
- Storage: candidate VCFs were removed only after comparison; the three actual
  Original golden VCFs remain retained under the locked baseline directory
- Cache policy: operating-system page cache was not flushed

`summary.tsv` is the comparison table. `all-runs.tsv` contains wall time,
user/system CPU, CPU utilization, RSS, bytes, SHA-256, and exactness for every
run. `manifest.tsv` identifies the source revision, binary version, workload,
lock, and benchmark policy.
