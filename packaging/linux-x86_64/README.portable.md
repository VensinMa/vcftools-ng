# vcftools-ng @VERSION@ portable Linux x86_64

This archive is ready to run after extraction. `bin` contains only the public
`vcftools-ng` launcher. Private bcftools 1.24 support for automatic CSI
construction is stored under `libexec`; HTSlib 1.24 and the required non-glibc
runtime libraries, including private libdeflate 1.25, are also bundled.

Minimum platform:

- Linux x86_64
- glibc 2.17 or newer
- tested on CentOS 7 and Ubuntu

Run:

```bash
./bin/vcftools-ng --version
./bin/vcftools-ng --vcf input.vcf --threads 8 --counts --out result
```

Keep `bin`, `lib`, and `libexec` together. The launcher finds bundled files
relative to its own location, so the extracted directory can be moved.

For source, documentation, and licenses, see:
https://github.com/VensinMa/vcftools-ng
