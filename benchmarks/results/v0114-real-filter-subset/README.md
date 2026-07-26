# Locked real-filter baseline

This directory is the first locked baseline for the permanent development
gate. The workload is the seven real project filters followed by
`--recode --recode-INFO-all --stdout` on 2,300,000 records and 412 samples.

Daily optimization uses only BGZF VCF + TBI, Plain VCF, and the BCF adaptive
streaming full-scan path. The BCF row uses the default automatic policy and
requires it to select stream while ignoring a neighbouring CSI.
`baseline.lock.tsv` fixes the input, relevant index, actual Original VCF
output, and Original timing identities. The three actual golden VCFs remain
stored locally under `golden/`; they are intentionally excluded from Git
because they occupy 34 GB. `summary.tsv`, `all-runs.tsv`, the lock,
environment report, and logs are compact reproducibility evidence.

The seven-scenario run in this directory is the initial v0.11.4 investigation,
not a recurring development gate. It must not be repeated unless a release
candidate is explicitly authorized. Automatic-CSI and no-auto-index rows are
therefore retained as reference evidence only.

Candidate outputs are compared with `cmp` against the retained actual golden,
then removed to prevent every development iteration from consuming another
hundreds of gigabytes. Their size, SHA-256, timing, CPU, RSS, and exact status
remain in the run TSV.
