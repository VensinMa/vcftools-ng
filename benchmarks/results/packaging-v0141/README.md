# v0.14.1 Linux x86_64 portable-package validation

Date: 2026-08-11

Artifact:

- `vcftools-ng-v0.14.1-linux-x86_64.tar.gz`;
- size: 5,332,319 bytes;
- SHA-256:
  `ce225f800bf0cede6151ad178b26c4d9c0a08ecac0f4225db09bd649dc21ecfa`;
- build baseline: `quay.io/pypa/manylinux2014_x86_64`.

Validation commands:

```bash
packaging/linux-x86_64/build-portable.sh
packaging/linux-x86_64/test-portable.sh \
  dist/vcftools-ng-v0.14.1-linux-x86_64.tar.gz
```

| Container | Result |
|---|---|
| CentOS 7 | PASS |
| Ubuntu 20.04 | PASS |

Both clean containers verified:

- `vcftools-ng 0.14.1`;
- only `bin/vcftools-ng` is publicly exposed;
- private `libexec/bcftools.bin` reports bcftools 1.24 and HTSlib 1.24;
- no missing dynamic dependency;
- colored help and `NO_COLOR` behavior;
- automatic four-thread CSI construction with the indexed direct kernel;
- standard `PREFIX.log` completion metadata;
- deterministic default BGZF output decompressed byte-for-byte against the
  retained real 23,000-record golden.

The public archive is Release+IPO/LTO on the CentOS 7-compatible baseline.
The local full-data performance gate used the separately identified opt-in
PGO build; the two artifacts share the same runtime source tree and exactness
contract, but their timing identities are not conflated.
