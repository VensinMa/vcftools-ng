# v0.13.1 Linux x86_64 portable-package validation

Date: 2026-08-06

Artifact:

- `vcftools-ng-v0.13.1-linux-x86_64.tar.gz`
- size: 5,312,940 bytes
- SHA-256: `489d34f5562ca02e815624a28a7f9395d64a44cb6fadcd482b802f7129a2abb4`
- build baseline: `quay.io/pypa/manylinux2014_x86_64`

Validation command:

```bash
packaging/linux-x86_64/build-portable.sh
packaging/linux-x86_64/test-portable.sh \
  dist/vcftools-ng-v0.13.1-linux-x86_64.tar.gz
```

Results:

| Container | Result |
|---|---|
| CentOS 7 | PASS |
| Ubuntu 20.04 | PASS |

Both clean containers verified:

- `vcftools-ng 0.13.1`;
- only `bin/vcftools-ng` is publicly exposed;
- private `libexec/bcftools.bin` reports bcftools 1.24 and HTSlib 1.24;
- no missing dynamic dependency;
- colored help and `NO_COLOR` behavior;
- automatic four-thread CSI construction;
- standard `PREFIX.log` completion metadata;
- deterministic BGZF output decompressed byte-for-byte against the retained
  real 23,000-record golden.
