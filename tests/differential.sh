#!/usr/bin/env bash
set -euo pipefail

ng_binary=${1:?vcftools-ng binary required}
project_root=${2:?project root required}
fixture="$project_root/tests/fixtures/osmanthus412.23chr_100k.bcf"
golden="$project_root/tests/golden"
output_root="$project_root/tests/output"

mkdir -p "$output_root"
output=$(mktemp -d "$output_root/differential.XXXXXX")
trap 'rm -rf -- "$output"' EXIT

"$ng_binary" \
    --bcf "$fixture" \
    --threads 8 \
    --freq \
    --counts \
    --missing-site \
    --site-depth \
    --site-mean-depth \
    --out "$output/all"

cmp "$golden/subset-freq.frq" "$output/all.frq"
cmp "$golden/subset-counts.frq.count" "$output/all.frq.count"
cmp "$golden/subset-missing-site.lmiss" "$output/all.lmiss"
cmp "$golden/subset-site-depth.ldepth" "$output/all.ldepth"
cmp "$golden/subset-site-mean-depth.ldepth.mean" "$output/all.ldepth.mean"

filters=(
    --min-alleles 2
    --max-alleles 2
    --remove-indels
    --minQ 40
    --minGQ 20
    --minDP 5
    --maxDP 30
    --min-meanDP 10
    --max-missing 0.9
    --maf 0.1
)

"$ng_binary" \
    --bcf "$fixture" \
    --threads 8 \
    "${filters[@]}" \
    --freq \
    --counts \
    --missing-site \
    --site-depth \
    --site-mean-depth \
    --recode \
    --recode-INFO-all \
    --out "$output/filtered"

cmp "$golden/subset-filtered-freq.frq" "$output/filtered.frq"
cmp "$golden/subset-filtered-counts.frq.count" "$output/filtered.frq.count"
cmp "$golden/subset-filtered-missing-site.lmiss" "$output/filtered.lmiss"
cmp "$golden/subset-filtered-site-depth.ldepth" "$output/filtered.ldepth"
cmp \
    "$golden/subset-filtered-site-mean-depth.ldepth.mean" \
    "$output/filtered.ldepth.mean"
cmp \
    "$golden/subset-filtered-recode-info-all.recode.vcf" \
    "$output/filtered.recode.vcf"
rm -f -- "$output/filtered.recode.vcf"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 16 \
    "${filters[@]}" \
    --recode \
    --out "$output/recode-no-info"

cmp \
    "$golden/subset-filtered-recode.recode.vcf" \
    "$output/recode-no-info.recode.vcf"
rm -f -- "$output/recode-no-info.recode.vcf"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 1 \
    "${filters[@]}" \
    --freq \
    --counts \
    --missing-site \
    --site-depth \
    --site-mean-depth \
    --out "$output/single-thread"

cmp "$output/filtered.frq" "$output/single-thread.frq"
cmp "$output/filtered.frq.count" "$output/single-thread.frq.count"
cmp "$output/filtered.lmiss" "$output/single-thread.lmiss"
cmp "$output/filtered.ldepth" "$output/single-thread.ldepth"
cmp "$output/filtered.ldepth.mean" "$output/single-thread.ldepth.mean"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 8 \
    --positions "$project_root/tests/fixtures/positions.keep.txt" \
    --exclude-positions "$project_root/tests/fixtures/positions.exclude.txt" \
    --not-chr chr2 \
    --recode \
    --recode-INFO-all \
    --out "$output/position-selection"

cmp \
    "$golden/subset-position-selection-recode.recode.vcf" \
    "$output/position-selection.recode.vcf"
rm -f -- "$output/position-selection.recode.vcf"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 16 \
    --chr chr7 \
    --from-bp 1 \
    --to-bp 2000000 \
    --recode \
    --recode-INFO-all \
    --out "$output/region"

cmp \
    "$golden/subset-region-recode.recode.vcf" \
    "$output/region.recode.vcf"
rm -f -- "$output/region.recode.vcf"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 8 \
    --keep "$project_root/tests/fixtures/samples.keep.txt" \
    --indv W-DA-8 \
    --remove "$project_root/tests/fixtures/samples.remove.txt" \
    --remove-indv W-CP-5 \
    "${filters[@]}" \
    --freq \
    --counts \
    --missing-site \
    --site-depth \
    --site-mean-depth \
    --recode \
    --recode-INFO-all \
    --out "$output/samples"

cmp "$golden/subset-samples-freq.frq" "$output/samples.frq"
cmp "$golden/subset-samples-counts.frq.count" "$output/samples.frq.count"
cmp "$golden/subset-samples-missing-site.lmiss" "$output/samples.lmiss"
cmp "$golden/subset-samples-site-depth.ldepth" "$output/samples.ldepth"
cmp \
    "$golden/subset-samples-site-mean-depth.ldepth.mean" \
    "$output/samples.ldepth.mean"
cmp \
    "$golden/subset-samples-recode.recode.vcf" \
    "$output/samples.recode.vcf"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 16 \
    --keep "$project_root/tests/fixtures/samples.keep.txt" \
    --indv W-DA-8 \
    --remove "$project_root/tests/fixtures/samples.remove.txt" \
    --remove-indv W-CP-5 \
    "${filters[@]}" \
    --mac 4 \
    --max-mac 12 \
    --hwe 0.001 \
    --recode \
    --recode-INFO-all \
    --out "$output/mac-hwe"

cmp \
    "$golden/subset-mac-hwe-recode.recode.vcf" \
    "$output/mac-hwe.recode.vcf"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 8 \
    --bed "$project_root/tests/fixtures/regions.compatibility.bed" \
    --counts \
    --out "$output/bed-include"

cmp \
    "$golden/subset-bed-include-counts.frq.count" \
    "$output/bed-include.frq.count"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 16 \
    --exclude-bed "$project_root/tests/fixtures/regions.compatibility.bed" \
    --counts \
    --out "$output/bed-exclude"

cmp \
    "$golden/subset-bed-exclude-counts.frq.count" \
    "$output/bed-exclude.frq.count"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 16 \
    --bed "$project_root/tests/fixtures/regions.compatibility.bed" \
    --thin 10000 \
    --counts \
    --out "$output/bed-thin"

cmp \
    "$golden/subset-bed-thin-counts.frq.count" \
    "$output/bed-thin.frq.count"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 8 \
    --bed "$project_root/tests/fixtures/regions.compatibility.bed" \
    --keep "$project_root/tests/fixtures/samples.keep.txt" \
    --indv W-DA-8 \
    --remove "$project_root/tests/fixtures/samples.remove.txt" \
    --remove-indv W-CP-5 \
    --non-ref-af 0.1 \
    --max-non-ref-af 0.8 \
    --non-ref-af-any 0.2 \
    --max-non-ref-af-any 0.7 \
    --non-ref-ac 2 \
    --max-non-ref-ac 30 \
    --non-ref-ac-any 5 \
    --max-non-ref-ac-any 20 \
    --counts \
    --out "$output/non-ref"

cmp \
    "$golden/subset-non-ref-counts.frq.count" \
    "$output/non-ref.frq.count"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 16 \
    --non-ref-af-any 0.99 \
    --counts \
    --out "$output/non-ref-af-any-compat"

cmp \
    "$golden/subset-counts.frq.count" \
    "$output/non-ref-af-any-compat.frq.count"

flag_fixture="$project_root/tests/fixtures/osmanthus412.flags.23chr_1k.bcf"

"$ng_binary" \
    --bcf "$flag_fixture" \
    --threads 16 \
    --keep-filtered q10 \
    --keep-filtered PASS \
    --remove-filtered Cluster \
    --keep-INFO Hotspot \
    --remove-INFO Artifact \
    --recode \
    --recode-INFO-all \
    --out "$output/flags-site-info"

cmp \
    "$golden/flags-site-info.recode.vcf" \
    "$output/flags-site-info.recode.vcf"

"$ng_binary" \
    --bcf "$flag_fixture" \
    --threads 8 \
    --remove-filtered-all \
    --counts \
    --out "$output/flags-site-remove-all"

cmp \
    "$golden/flags-site-remove-all.frq.count" \
    "$output/flags-site-remove-all.frq.count"

flag_sample_filters=(
    --keep "$project_root/tests/fixtures/samples.keep.txt"
    --indv W-DA-8
    --remove "$project_root/tests/fixtures/samples.remove.txt"
    --remove-indv W-CP-5
)

"$ng_binary" \
    --bcf "$flag_fixture" \
    --threads 8 \
    "${flag_sample_filters[@]}" \
    --remove-filtered-geno-all \
    --counts \
    --missing-site \
    --recode \
    --recode-INFO-all \
    --out "$output/flags-ft-all"

cmp "$golden/flags-ft-all.frq.count" "$output/flags-ft-all.frq.count"
cmp "$golden/flags-ft-all.lmiss" "$output/flags-ft-all.lmiss"
cmp "$golden/flags-ft-all.recode.vcf" "$output/flags-ft-all.recode.vcf"

"$ng_binary" \
    --bcf "$flag_fixture" \
    --threads 16 \
    "${flag_sample_filters[@]}" \
    --remove-filtered-geno LowDP \
    --remove-filtered-geno MendelFail \
    --recode \
    --recode-INFO-all \
    --out "$output/flags-ft-specific"

cmp \
    "$golden/flags-ft-specific.recode.vcf" \
    "$output/flags-ft-specific.recode.vcf"

if "$ng_binary" \
    --bcf "$flag_fixture" \
    --threads 8 \
    --keep-INFO AC \
    --counts \
    --out "$output/invalid-info-type" >/dev/null 2>&1; then
    echo "--keep-INFO accepted a non-Flag INFO field" >&2
    exit 1
fi

echo "All 2,300,000-record outputs and 23,000-record flag semantics are byte-identical to VCFtools 0.1.17."
