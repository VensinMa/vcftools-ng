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

Validate the archive in clean CentOS 7 and Ubuntu 20.04 containers:

```bash
./packaging/linux-x86_64/test-portable.sh
```

Artifacts are written to `dist/`:

```text
vcftools-ng-vX.Y.Z-linux-x86_64.tar.gz
vcftools-ng-vX.Y.Z-linux-x86_64.tar.gz.sha256
```

The build installs the CentOS SCL GCC 11 toolchain, pins and verifies the
HTSlib 1.24 and bcftools 1.24 source archives, builds vcftools-ng in the same
container, recursively bundles non-glibc shared-library dependencies, strips
release binaries, and emits a stable archive file order and metadata.
Only `bin/vcftools-ng` is exposed as a user command; the bundled bcftools
binary remains private under `libexec`.
