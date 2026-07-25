# Portable Linux x86_64 release

`build-portable.sh` builds the release archive in the
`quay.io/pypa/manylinux2014_x86_64` container. This fixes the runtime baseline
at glibc 2.17 for CentOS 7 compatibility and prevents local Conda or system
libraries from leaking into the archive.

Requirements:

- Docker
- an x86_64 Linux build host

Build:

```bash
./packaging/linux-x86_64/build-portable.sh
```

Artifacts are written to `dist/`:

```text
vcftools-ng-v0.11.2-linux-x86_64.tar.gz
vcftools-ng-v0.11.2-linux-x86_64.tar.gz.sha256
```

The build pins and verifies the HTSlib 1.23 and bcftools 1.23 source
archives, builds vcftools-ng in the same container, recursively bundles
non-glibc shared-library dependencies, strips release binaries, and emits a
stable archive file order and metadata.
