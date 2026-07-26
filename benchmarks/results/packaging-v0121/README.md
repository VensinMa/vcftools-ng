# v0.12.1 Linux x86_64 portable-package validation

Status: **PASS**.

## Artifact

- file: `vcftools-ng-v0.12.1-linux-x86_64.tar.gz`
- bytes: 5,187,794
- SHA-256:
  `b00e93877c8d97ecc09215e958f5e302b62b633ea4f8c1a9b871ab5e68174547`
- runtime baseline: Linux x86_64, glibc 2.17 or newer
- bundled tools: vcftools-ng 0.12.1, bcftools 1.23, HTSlib 1.23

## Build

The archive was built in
`quay.io/pypa/manylinux2014_x86_64:latest`. The image's default GCC 10
libstdc++ lacks floating-point `std::to_chars/from_chars`, so the packaging
script installs and activates CentOS SCL devtoolset-11. The corresponding
non-glibc C++ runtime is bundled. Runtime source code was not changed.

The package builder also rejects unresolved dependencies and no longer treats
the `rpm -qf` text for locally built libraries as an RPM package name.
Archive inspection confirmed that no erroneous temporary/license paths were
included.

## Clean-container verification

Both `centos:7` and `ubuntu:20.04` passed:

- archive SHA-256;
- extraction and `vcftools-ng --version`;
- bundled bcftools/HTSlib version;
- no unresolved `ldd` dependency;
- automatic CSI construction with the requested four-thread budget;
- complete-file comparison of a real 23,000-record filtered recode against
  the retained Original golden.

Reproduce with:

```bash
./packaging/linux-x86_64/test-portable.sh
```
