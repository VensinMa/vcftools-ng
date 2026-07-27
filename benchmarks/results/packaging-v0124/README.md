# v0.12.4 Linux x86_64 portable-package validation

Status: **PASS**.

## Artifact

- file: `vcftools-ng-v0.12.4-linux-x86_64.tar.gz`
- bytes: 5,226,503
- SHA-256:
  `873e95578baa63a7bb8cb21c389cc29b7021a71286c841091091bf92c9ac93a6`
- runtime baseline: Linux x86_64, glibc 2.17 or newer
- bundled tools: vcftools-ng 0.12.4, bcftools 1.23, HTSlib 1.23

## Build and clean-container verification

The archive was built with the CentOS 7-compatible
`quay.io/pypa/manylinux2014_x86_64:latest` toolchain. Both `centos:7` and
`ubuntu:20.04` passed:

- archive SHA-256;
- extraction and `vcftools-ng --version`;
- bundled bcftools/HTSlib version;
- no unresolved `ldd` dependency;
- forced ANSI-colored help and plain `NO_COLOR` help;
- automatic CSI construction with a four-thread budget;
- default standard log creation and index-decision fields;
- complete-file comparison of a real 23,000-record filtered recode against
  the retained Original golden.

Reproduce with:

```bash
./packaging/linux-x86_64/build-portable.sh
./packaging/linux-x86_64/test-portable.sh
```
