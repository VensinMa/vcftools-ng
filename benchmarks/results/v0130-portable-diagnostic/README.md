# v0.13.0 portable-build diagnostic

This retained, non-release benchmark records one full-data run made with the
provisional glibc-2.17 portable package before the official matrix began.

- Input: BGZF VCF + TBI on the same NVMe device as the output.
- Output: uncompressed VCF.
- Threads: 1.
- Application wall time: 358.266018 seconds.
- Durable wall time: 358.941389 seconds.
- Result: byte-identical to the locked Original VCFtools 0.1.17 oracle.

The same source revision built natively completed the earlier development run
in 246.984843 seconds. The portable package deliberately targets an older
glibc/toolchain baseline and its HTSlib dependency set differs from the native
benchmark build. Mixing these values into one scaling table would therefore
be invalid. The release input/output/storage matrix uses the native Release
binary, matching the historical benchmark convention, while automatic CSI
construction is isolated through the bundled bcftools 1.24 wrapper.
