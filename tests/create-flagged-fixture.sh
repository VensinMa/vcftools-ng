#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
bcftools_bin=${BCFTOOLS:-/home/vensin/software/bcftools/bcftools}
bgzip_bin=${BGZIP:-/home/vensin/anaconda3/envs/hifiasm/bin/bgzip}
tabix_bin=${TABIX:-/home/vensin/anaconda3/envs/hifiasm/bin/tabix}
source_bcf=${1:-"$project_root/tests/fixtures/osmanthus412.23chr_100k.bcf"}
output_prefix=${2:-"$project_root/tests/fixtures/osmanthus412.flags.23chr_1k"}

mkdir -p "$project_root/tests/output"
work_dir=$(mktemp -d "$project_root/tests/output/flagged-fixture.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT

positions="$work_dir/positions.txt"
"$bcftools_bin" query -f '%CHROM\t%POS\n' "$source_bcf" |
    awk '
        $1 != chromosome {
            chromosome = $1
            count = 0
        }
        count < 1000 {
            print
            count++
        }
    ' >"$positions"

base_bcf="$work_dir/base.bcf"
"$bcftools_bin" view \
    --no-version \
    -R "$positions" \
    -Ob \
    -o "$base_bcf" \
    "$source_bcf"
"$bcftools_bin" index -f "$base_bcf"

samples="$work_dir/samples.txt"
"$bcftools_bin" query -l "$base_bcf" >"$samples"
sample_count=$(wc -l <"$samples")
annotation_vcf="$work_dir/annotation.vcf"

{
    printf '##fileformat=VCFv4.2\n'
    "$bcftools_bin" view -h "$base_bcf" | grep '^##contig='
    printf '##FILTER=<ID=q10,Description="Compatibility low-quality flag">\n'
    printf '##FILTER=<ID=Cluster,Description="Compatibility cluster flag">\n'
    printf '##INFO=<ID=Hotspot,Number=0,Type=Flag,Description="Compatibility hotspot flag">\n'
    printf '##INFO=<ID=Artifact,Number=0,Type=Flag,Description="Compatibility artifact flag">\n'
    printf '##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype-level filter">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
    while IFS= read -r sample; do
        printf '\t%s' "$sample"
    done <"$samples"
    printf '\n'
    "$bcftools_bin" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$base_bcf" |
        awk -v samples="$sample_count" '
            {
                mode = NR % 8
                if (mode == 0) {
                    filter = "PASS"
                    info = "."
                } else if (mode == 1 || mode == 5) {
                    filter = "q10"
                    info = "Hotspot"
                } else if (mode == 2 || mode == 6) {
                    filter = "Cluster"
                    info = "Artifact"
                } else if (mode == 3 || mode == 7) {
                    filter = "q10;Cluster"
                    info = "Hotspot;Artifact"
                } else {
                    filter = "."
                    info = "."
                }
                printf "%s\t%s\t.\t%s\t%s\t.\t%s\t%s\tFT",
                    $1, $2, $3, $4, filter, info
                for (sample = 1; sample <= samples; sample++) {
                    genotype_mode = (sample + NR) % 13
                    if (genotype_mode == 0)
                        ft = "LowDP"
                    else if (genotype_mode == 1)
                        ft = "MendelFail"
                    else if (genotype_mode == 2)
                        ft = "LowDP;MendelFail"
                    else if (genotype_mode == 3)
                        ft = "PASS;LowDP"
                    else if (genotype_mode == 4)
                        ft = "LowDP;PASS"
                    else if (genotype_mode == 5)
                        ft = "."
                    else
                        ft = "PASS"
                    printf "\t%s", ft
                }
                printf "\n"
            }
        '
} >"$annotation_vcf"

"$bgzip_bin" -f "$annotation_vcf"
"$tabix_bin" -f -p vcf "$annotation_vcf.gz"

"$bcftools_bin" annotate \
    --no-version \
    -a "$annotation_vcf.gz" \
    -c FILTER,INFO/Hotspot,INFO/Artifact,FORMAT/FT \
    -Ob \
    -o "$output_prefix.bcf" \
    "$base_bcf"
"$bcftools_bin" index -f "$output_prefix.bcf"
"$bcftools_bin" view \
    --no-version \
    -Oz \
    -o "$output_prefix.vcf.gz" \
    "$output_prefix.bcf"
"$tabix_bin" -f -p vcf "$output_prefix.vcf.gz"

"$bcftools_bin" index -n "$output_prefix.bcf"
