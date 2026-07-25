# v0.11.1 automatic CSI benchmark

The input is the 3,956,690,872-byte BGZF VCF fixture containing 2,300,000
real records. Each vcftools-ng run started from a different path without a
CSI or TBI sidecar. Wall time includes both `bcftools index --csi` and
vcftools-ng `--counts` processing.

The original baseline is the saved VCFtools 0.1.17 no-index measurement from
`benchmarks/results/input-vnext-counts`: 68.22 seconds. Both generated count
files matched its saved SHA-256,
`c2d9a316a7630f092eeec8ab89043f0a5dcc2b0bd579da95b2776bae566731cd`.

Raw timing and stderr files in this directory preserve the exact command
behavior, selected indexed backend, requested index thread count, and
2,300,000-record completion.
