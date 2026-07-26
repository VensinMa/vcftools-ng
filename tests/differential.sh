#!/usr/bin/env bash
set -euo pipefail

ng_binary=${1:?vcftools-ng binary required}
project_root=${2:?project root required}
fixture="$project_root/tests/fixtures/osmanthus412.23chr_100k.bcf"
vcf_fixture="$project_root/tests/fixtures/osmanthus412.23chr_100k.vcf.gz"
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

"$ng_binary" \
    --bcf "$fixture" \
    --threads 8 \
    --depth \
    --missing-indv \
    --het \
    --hardy \
    --site-quality \
    --out "$output/v010-statistics"

cmp "$golden/v010-depth.idepth" "$output/v010-statistics.idepth"
cmp "$golden/v010-missing-indv.imiss" "$output/v010-statistics.imiss"
cmp "$golden/v010-het.het" "$output/v010-statistics.het"
cmp "$golden/v010-hardy.hwe" "$output/v010-statistics.hwe"
cmp \
    "$golden/v010-site-quality.lqual" \
    "$output/v010-statistics.lqual"

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
    "$golden/subset-filtered-bcf-recode-info-all.recode.vcf" \
    "$output/filtered.recode.vcf"
rm -f -- "$output/filtered.recode.vcf"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 16 \
    "${filters[@]}" \
    --recode \
    --out "$output/recode-no-info"

cmp \
    "$golden/subset-filtered-bcf-recode.recode.vcf" \
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
    --gzvcf "$vcf_fixture" \
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
    --gzvcf "$vcf_fixture" \
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
    --gzvcf "$vcf_fixture" \
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
    --gzvcf "$vcf_fixture" \
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
flag_vcf_fixture="$project_root/tests/fixtures/osmanthus412.flags.23chr_1k.vcf.gz"

"$ng_binary" \
    --gzvcf "$flag_vcf_fixture" \
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
    --gzvcf "$flag_vcf_fixture" \
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
    --gzvcf "$flag_vcf_fixture" \
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

"$ng_binary" \
    --bcf "$flag_fixture" \
    --threads 16 \
    "${flag_sample_filters[@]}" \
    --minGQ 10 \
    --minDP 3 \
    --maxDP 60 \
    --remove-filtered-geno-all \
    --thin 100 \
    --depth \
    --missing-indv \
    --het \
    --hardy \
    --site-quality \
    --out "$output/flags-v010-filtered"

cmp \
    "$golden/flags-v010-filtered.idepth" \
    "$output/flags-v010-filtered.idepth"
cmp \
    "$golden/flags-v010-filtered.imiss" \
    "$output/flags-v010-filtered.imiss"
cmp \
    "$golden/flags-v010-filtered.het" \
    "$output/flags-v010-filtered.het"
cmp \
    "$golden/flags-v010-filtered.hwe" \
    "$output/flags-v010-filtered.hwe"
cmp \
    "$golden/flags-v010-filtered.lqual" \
    "$output/flags-v010-filtered.lqual"

"$ng_binary" \
    --bcf "$flag_fixture" \
    --threads 8 \
    --site-pi \
    --window-pi 100000 \
    --window-pi-step 50000 \
    --TajimaD 100000 \
    --out "$output/flags-v010-pi"

cmp \
    "$golden/flags-v010-site-pi.sites.pi" \
    "$output/flags-v010-pi.sites.pi"
cmp \
    "$golden/flags-v010-window-pi.windowed.pi" \
    "$output/flags-v010-pi.windowed.pi"
cmp \
    "$golden/flags-v010-tajima.Tajima.D" \
    "$output/flags-v010-pi.Tajima.D"

fst_populations=(
    --weir-fst-pop "$project_root/tests/fixtures/population-a.txt"
    --weir-fst-pop "$project_root/tests/fixtures/population-b.txt"
)

"$ng_binary" \
    --bcf "$flag_fixture" \
    --threads 16 \
    "${fst_populations[@]}" \
    --out "$output/flags-v010-site-fst"

cmp \
    "$golden/flags-v010-site-fst.weir.fst" \
    "$output/flags-v010-site-fst.weir.fst"

"$ng_binary" \
    --bcf "$flag_fixture" \
    --threads 8 \
    "${fst_populations[@]}" \
    --fst-window-size 100000 \
    --fst-window-step 50000 \
    --out "$output/flags-v010-window-fst"

cmp \
    "$golden/flags-v010-window-fst.windowed.weir.fst" \
    "$output/flags-v010-window-fst.windowed.weir.fst"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 8 \
    --site-pi \
    --window-pi 100000 \
    --window-pi-step 50000 \
    --TajimaD 100000 \
    --out "$output/v010-pi"

cmp "$golden/v010-site-pi.sites.pi" "$output/v010-pi.sites.pi"
cmp \
    "$golden/v010-window-pi.windowed.pi" \
    "$output/v010-pi.windowed.pi"
cmp "$golden/v010-tajima.Tajima.D" "$output/v010-pi.Tajima.D"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 8 \
    "${fst_populations[@]}" \
    --out "$output/v010-site-fst"

cmp \
    "$golden/v010-site-fst.weir.fst" \
    "$output/v010-site-fst.weir.fst"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 16 \
    "${fst_populations[@]}" \
    --fst-window-size 100000 \
    --fst-window-step 50000 \
    --out "$output/v010-window-fst"

cmp \
    "$golden/v010-window-fst.windowed.weir.fst" \
    "$output/v010-window-fst.windowed.weir.fst"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 8 \
    --thin 50000 \
    --geno-r2 \
    --ld-window 20 \
    --ld-window-bp 500000 \
    --min-r2 0.01 \
    --out "$output/v010-geno-r2"

cmp \
    "$golden/v010-geno-r2.geno.ld" \
    "$output/v010-geno-r2.geno.ld"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 16 \
    --keep "$project_root/tests/fixtures/population-a.txt" \
    --max-missing 1 \
    --thin 100000 \
    --pca \
    --out "$output/v010-pca"

cmp "$golden/v010-pca.pca" "$output/v010-pca.pca"

"$ng_binary" \
    --bcf "$fixture" \
    --threads 8 \
    --keep "$project_root/tests/fixtures/population-a.txt" \
    --max-missing 1 \
    --thin 100000 \
    --pca-no-norm \
    --out "$output/v010-pca-no-norm"

cmp \
    "$golden/v010-pca-no-norm.pca" \
    "$output/v010-pca-no-norm.pca"

for threads in 8 16; do
    "$ng_binary" \
        --bcf "$fixture" \
        --threads "$threads" \
        --chr chr1 \
        --to-bp 1000000 \
        --recode-bcf \
        --recode-INFO-all \
        --out "$output/v010-bcf-t$threads"

    cmp \
        "$golden/v010-chr1-1m.recode.bcf" \
        "$output/v010-bcf-t$threads.recode.bcf"
    rm -f -- "$output/v010-bcf-t$threads.recode.bcf"
done

for threads in 8 16; do
    "$ng_binary" \
        --gzvcf "$flag_vcf_fixture" \
        --threads "$threads" \
        --keep-filtered q10 \
        --keep-filtered PASS \
        --remove-filtered Cluster \
        --keep-INFO Hotspot \
        --remove-INFO Artifact \
        --recode-vcf-gz \
        --recode-INFO-all \
        --out "$output/v010-vcf-gz-t$threads"

    gzip -dc \
        "$output/v010-vcf-gz-t$threads.recode.vcf.gz" \
        > "$output/v010-vcf-gz-t$threads.recode.vcf"
    cmp \
        "$golden/flags-site-info.recode.vcf" \
        "$output/v010-vcf-gz-t$threads.recode.vcf"
    rm -f -- \
        "$output/v010-vcf-gz-t$threads.recode.vcf.gz" \
        "$output/v010-vcf-gz-t$threads.recode.vcf"
done

for mode in \
    diff-site \
    diff-indv \
    diff-site-discordance \
    diff-indv-discordance
do
    "$ng_binary" \
        --bcf "$flag_fixture" \
        --diff-bcf "$flag_fixture" \
        --"$mode" \
        --threads 8 \
        --out "$output/flags-v010-$mode"
done

cmp \
    "$golden/flags-v010-diff.sites_in_files" \
    "$output/flags-v010-diff-site.diff.sites_in_files"
cmp \
    "$golden/flags-v010-diff.indv_in_files" \
    "$output/flags-v010-diff-indv.diff.indv_in_files"
cmp \
    "$golden/flags-v010-diff.sites" \
    "$output/flags-v010-diff-site-discordance.diff.sites"
cmp \
    "$golden/flags-v010-diff.indv" \
    "$output/flags-v010-diff-indv-discordance.diff.indv"

"$ng_binary" \
    --bcf "$fixture" \
    --diff-bcf "$fixture" \
    --diff-site-discordance \
    --threads 16 \
    --out "$output/v010-diff-site-discordance"

cmp \
    "$golden/v010-diff-site-discordance.diff.sites" \
    "$output/v010-diff-site-discordance.diff.sites"

if "$ng_binary" \
    --bcf "$flag_fixture" \
    --threads 8 \
    --keep-INFO AC \
    --counts \
    --out "$output/invalid-info-type" >/dev/null 2>&1; then
    echo "--keep-INFO accepted a non-Flag INFO field" >&2
    exit 1
fi

echo "All v0.11 cumulative 2,300,000-record outputs, exact BCF conversions, and 23,000-record flag/diff semantics are byte-identical to VCFtools 0.1.17."
