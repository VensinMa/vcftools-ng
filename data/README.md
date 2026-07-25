# Benchmark data

Simplified dataset prefix: `osmanthus412.snps`

## Dataset

- Samples: 412
- Variant records: 11,230,392
- Variants: filtered, biallelic SNPs
- First record: `chr1:5330 T>A`
- Last record: `chr23:22351121 C>T`

## Files

| File | Kind | Size | Storage |
|---|---|---:|---|
| `osmanthus412.snps.vcf.gz` | BGZF-compressed VCF | 18 GB | Symlink to the original file |
| `osmanthus412.snps.vcf.gz.tbi` | TBI index | 527 KB | Symlink to the original index |
| `osmanthus412.snps.bcf` | BCF | 21 GB | Generated locally |
| `osmanthus412.snps.bcf.csi` | CSI index | 471 KB | Generated locally |
| `osmanthus412.snps.vcf` | Uncompressed VCF | 115 GB | Generated locally |
| `osmanthus412.subset.vcf` | Uncompressed 2.3M-record VCF | 23.4 GiB | Direct BGZF decompression |

## Conversion

BCF:

```bash
bcftools view --threads 28 -Ob \
  -o osmanthus412.snps.bcf \
  osmanthus412.snps.vcf.gz
bcftools index --threads 28 osmanthus412.snps.bcf
```

Uncompressed VCF:

```bash
bcftools view --threads 28 -Ov \
  -o osmanthus412.snps.vcf \
  osmanthus412.snps.bcf
```

Measured conversion times on the development machine:

- BGZF VCF to BCF: 4 minutes 42.96 seconds
- BCF to uncompressed VCF: 3 minutes 27.99 seconds
- Raw line scan of uncompressed VCF: 41.53 seconds
- Real-subset BGZF decompression to plain VCF: 14.08 seconds

## gVCF

No gVCF was generated. This filtered cohort VCF contains variant sites, but not
the reference-confidence blocks required by gVCF (`<NON_REF>` alleles and
`INFO/END` intervals). Those blocks cannot be reconstructed losslessly from
this file. Creating a file with a `.g.vcf` suffix would not make it a valid
gVCF.

## Original source

```text
/home/vensin/workspace/Sweet_Osmanthus/05.variant_filter/02.vcftools_filter_snp_indel/412samples.SNP.biallelic.minGQ10.minQ30.meanDP6.maxmiss0.8.maf0.05.vcf.gz
```
