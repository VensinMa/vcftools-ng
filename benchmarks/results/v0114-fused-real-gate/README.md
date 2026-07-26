# v0.11.4 fused real-data gate

This gate closes the real-data evidence gap for the new direct text site-stat
path.

- Fixture: 2,300,000 records, 23 chromosomes, 412 samples
- Direct paths: BGZF VCF + TBI and Plain VCF
- Threads: 8 and 16
- Exactness: all eight candidate runs passed complete-file `cmp`
- Combined workload: `--freq --counts --missing-site --site-depth
  --site-mean-depth --site-quality`
- Separate workload: `--freq2`

The six existing Original goldens were reused. Because `--freq2` had only a
small synthetic golden, Original VCFtools 0.1.17 was run once to establish and
retain `tests/golden/subset-freq2.frq`:

- bytes: 86,146,636
- SHA-256:
  `34fb6cc9522a46504bd600c8f5ed7657745117df2f0f20a66139d50b875ad179`
- Original wall time: 68.16 seconds

The Original aggregate for the six BGZF scans was 393.74 seconds. Plain VCF
rows intentionally report `NA` for speedup because the retained Original
timings used BGZF input; their exact output and candidate runtime remain
recorded without making a cross-input performance claim.

Candidate artifacts were removed only after byte comparison. Actual goldens
remain stored locally and their identities remain in
`tests/golden/SHA256SUMS`.
