# v0.12.2 full-data four-scenario release matrix

Result: **PASS**.

- Candidate commit: `8980655d83b9b161fe558ce7d2d2a784289acc2a`
- Dataset: 11,230,392 records, 23 chromosomes, 412 samples
- Retained records: 5,425,725
- Scenarios: BGZF+TBI, BGZF+automatic CSI, Plain VCF, adaptive BCF
- Threads: 1/2/4/8/16/32
- Repeats: 5 per vcftools-ng configuration
- Exact outputs: 120/120 PASS
- Original runs in v0.12.2: 0

The Original VCFtools 0.1.17 timings and golden outputs were retained from
v0.12.1. The driver validated input, index, and oracle SHA-256 values before
reuse. Each vcftools-ng repeat ran strictly serially and performed a complete
`cmp`, SHA validation, and expected-backend assertion.

See:

- `summary.tsv` for five-repeat means, speedups, CPU, RSS, and backend;
- `all-runs.tsv` for every retained timing row;
- `manifest.tsv` for hashes, host, candidate commit, and timestamps;
- `../../run-v0122-full-release-matrix.sh` for the resumable driver.
