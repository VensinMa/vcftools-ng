# CSI index construction benchmark

Input: `tests/fixtures/osmanthus412.23chr_100k.vcf.gz`

- Records: 2,300,000
- BGZF VCF size: 3,956,690,872 bytes
- Tool: bcftools 1.23-1-g638bee2d, HTSlib 1.23-2-g6dd0d7d0
- Command: `bcftools index --csi --threads N --output OUTPUT INPUT`
- Runs were consecutive on the same local filesystem; operating-system caches were not flushed.
- The second one-thread run was placed after the 8- and 16-thread runs to check cache-order bias.

The warm one-thread repeat remained at 30.99 seconds, versus 11.78 seconds for
both 8 and 16 threads. The speedup therefore comes from parallel BGZF decoding,
not merely from the run order. Scaling saturates by eight requested threads on
this input and machine.

All four generated CSI files are 99,846 bytes and have SHA-256
`781d89e5199542bf2e1e91c104727236107e63eabecf0055303861a3c500ae6d`.
They are byte-identical.
