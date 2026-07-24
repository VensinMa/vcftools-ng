# Real-data compatibility fixture

`osmanthus412.23chr_100k` contains exactly 100,000 records from each of
`chr1` through `chr23`, for a total of 2,300,000 records and all 412 samples.
It is stored as both BCF plus CSI and BGZF VCF plus TBI.

The BCF and BGZF VCF headers were normalized to remove automatically-added
`bcftools_viewCommand` provenance lines. This makes the two encodings share an
identical VCF header while leaving every sample and variant record unchanged.
Both indexes report exactly 2,300,000 records.

`positions.keep.txt` and `positions.exclude.txt` contain real chr1–chr3
coordinates from this subset. They exercise the combined position inclusion,
position exclusion, and chromosome exclusion path.

`samples.keep.txt` and `samples.remove.txt` contain real sample IDs. Together
with inline `--indv` and `--remove-indv`, the compatibility case retains 17 of
the 412 samples.

For recode compatibility, VCFtools 0.1.17 reads the BGZF VCF fixture and
produces the golden VCF. `vcftools-ng` reads the equivalent BCF fixture for its
high-performance path. This is intentional: VCFtools 0.1.17's legacy BCF
reader writes embedded NUL padding from `Character` FORMAT fields into VCF
text, whereas its VCF reader and HTSlib produce identical valid record text.

The complete 11,230,392-record dataset is not used by the current regression
or benchmark loop.

`osmanthus412.flags.23chr_1k` is a small semantic fixture derived from the
same real records: the first 1,000 records of each chromosome, all 412
samples, and 23,000 records total. Deterministic annotations cover site
FILTER combinations, INFO Flags, and per-genotype FT combinations including
the order-sensitive `PASS;LowDP` and `LowDP;PASS` cases. No genotype or
variant allele was synthesized. Rebuild it with:

```bash
./tests/create-flagged-fixture.sh
```

VCFtools 0.1.17 reads the BGZF VCF encoding to create the flag goldens;
`vcftools-ng` reads the byte-equivalent BCF encoding.
