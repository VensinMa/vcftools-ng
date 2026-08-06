# v0.13.0 Linux x86_64 portable-package validation

Status: **PASS**.

## Artifact

- file: `vcftools-ng-v0.13.0-linux-x86_64.tar.gz`
- bytes: 5,280,148
- SHA-256:
  `cb743aa0663afaa78cfc91b4cc82e7d0ef818f09824ee872a403b8713f22a44e`
- runtime baseline: Linux x86_64, glibc 2.17 or newer
- bundled tools: vcftools-ng 0.13.0, bcftools 1.24, HTSlib 1.24

## Build and clean-container verification

The archive was built with the CentOS 7-compatible
`quay.io/pypa/manylinux2014_x86_64:latest` toolchain. Both `centos:7` and
`ubuntu:20.04` passed:

- archive SHA-256;
- extraction and `vcftools-ng --version`;
- exactly one public executable under `bin`;
- private bcftools/HTSlib 1.24 versions;
- relocatable `$ORIGIN/../lib` lookup for both private executables;
- no unresolved `ldd` dependency;
- forced ANSI-colored help and plain `NO_COLOR` help;
- automatic CSI construction with a four-thread budget;
- default standard log creation and index-decision fields;
- default BGZF VCF output; and
- complete decompressed comparison of a real 23,000-record filtered recode
  against the retained Original golden.

The first package attempt exposed a missing private-bcftools RPATH during the
direct dependency diagnostic. The build now links the helper with a relative
runtime path; the original failing clean-container test and the full portable
test both pass.

Reproduce with:

```bash
./packaging/linux-x86_64/build-portable.sh
./packaging/linux-x86_64/test-portable.sh
```
