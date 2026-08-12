# v0.14.3 strict CPU-budget development gate

This directory records the locked 230,000-record production-filter A/B used
after correcting v0.14.2 CPU oversubscription. Both binaries were restricted
to the same `0..N-1` CPU affinity, run once at 4/8/16/24/32 CPUs, and produced
the same BGZF VCF SHA-256 at every row. Differences within 5% are predefined
as performance ties.

Files:

- `fair-230k.tsv`: wall time, CPU, RSS, and output SHA-256;
- `SHA256SUMS`: labeled input and v0.14.3 candidate binary identities (a
  manifest, not a directly executable `sha256sum --check` file).

The scientific release gate remains the complete 2.3-million-record CTest
differential. This 230k table is a performance regression gate, not an
Original baseline or a full-data benchmark.
