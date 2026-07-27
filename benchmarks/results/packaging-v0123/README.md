# v0.12.3 Linux x86_64 portable-package validation

Status: **PASS**.

## Artifact

- file: `vcftools-ng-v0.12.3-linux-x86_64.tar.gz`
- bytes: 5,191,796
- SHA-256:
  `7131735169972bd081854526b043f2da1c8e36aad137af3d42515a964a11fc72`
- runtime baseline: Linux x86_64, glibc 2.17 or newer
- bundled tools: vcftools-ng 0.12.3, bcftools 1.23, HTSlib 1.23

## Build and clean-container verification

The archive was built with the CentOS 7-compatible
`quay.io/pypa/manylinux2014_x86_64:latest` toolchain. Both `centos:7` and
`ubuntu:20.04` passed:

- archive SHA-256;
- extraction and `vcftools-ng --version`;
- bundled bcftools/HTSlib version;
- no unresolved `ldd` dependency;
- forced ANSI-colored help;
- plain help under `NO_COLOR=1`, even when `CLICOLOR_FORCE=1` is also set;
- automatic CSI construction with a four-thread budget;
- complete-file comparison of a real 23,000-record filtered recode against
  the retained Original golden.

Reproduce with:

```bash
./packaging/linux-x86_64/build-portable.sh
./packaging/linux-x86_64/test-portable.sh
```
