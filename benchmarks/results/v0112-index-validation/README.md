# v0.11.2 validated-index benchmark

All cases use the real 2,300,000-record, 23-chromosome fixture and compare
the count-file SHA-256 with the saved VCFtools 0.1.17 oracle.

The BGZF runs started without a sidecar. Their wall time includes bcftools
CSI construction, explicit CSI validation, and counts processing. The Plain
VCF runs used the 25,173,128,179-byte uncompressed fixture. Logs confirm that
Plain input selected `plain-ranges`, never invoked automatic indexing, and
left no `.csi` or `.tbi`.

The original baselines are the saved v0.11 input-matrix measurements:
68.22 seconds for no-index BGZF VCF and 36.15 seconds for Plain VCF.
