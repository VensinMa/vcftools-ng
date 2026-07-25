#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ng_binary=${NG_BINARY:-"$project_root/build/vcftools-ng"}
original=${VCFTOOLS_ORIGINAL:-/home/vensin/anaconda3/envs/vcftools/bin/vcftools}
pca_original=${VCFTOOLS_PCA_ORACLE:-}
bcf="$project_root/tests/fixtures/osmanthus412.23chr_100k.bcf"
vcf="$project_root/tests/fixtures/osmanthus412.23chr_100k.vcf.gz"
golden="$project_root/tests/golden"
results="$project_root/benchmarks/results"
pop_a="$project_root/tests/fixtures/population-a.txt"
pop_b="$project_root/tests/fixtures/population-b.txt"

for path in "$ng_binary" "$original" "$bcf" "$vcf" "$pop_a" "$pop_b"; do
    if [[ ! -e "$path" ]]; then
        echo "Missing required path: $path" >&2
        exit 1
    fi
done

mkdir -p "$results" "$project_root/tests/output"
work=$(mktemp -d "$project_root/tests/output/v010-advanced.XXXXXX")
trap 'rm -rf -- "$work"' EXIT

timed() {
    local name=$1
    shift
    /usr/bin/time -v \
        -o "$results/$name.time.txt" \
        "$@" \
        >"$results/$name.stdout.txt" \
        2>"$results/$name.log"
}

echo "Validating combined pi/window-pi/Tajima output"
for threads in 8 16; do
    prefix="$work/pi-t$threads"
    timed "v010-ng${threads}-pi-combined" \
        "$ng_binary" \
        --bcf "$bcf" \
        --threads "$threads" \
        --site-pi \
        --window-pi 100000 \
        --window-pi-step 50000 \
        --TajimaD 100000 \
        --out "$prefix"
    cmp "$golden/v010-site-pi.sites.pi" "$prefix.sites.pi"
    cmp "$golden/v010-window-pi.windowed.pi" "$prefix.windowed.pi"
    cmp "$golden/v010-tajima.Tajima.D" "$prefix.Tajima.D"
done

echo "Validating site and window FST"
for threads in 8 16; do
    site_prefix="$work/site-fst-t$threads"
    window_prefix="$work/window-fst-t$threads"
    timed "v010-ng${threads}-site-fst" \
        "$ng_binary" \
        --bcf "$bcf" \
        --threads "$threads" \
        --weir-fst-pop "$pop_a" \
        --weir-fst-pop "$pop_b" \
        --out "$site_prefix"
    cmp "$golden/v010-site-fst.weir.fst" "$site_prefix.weir.fst"
    timed "v010-ng${threads}-window-fst" \
        "$ng_binary" \
        --bcf "$bcf" \
        --threads "$threads" \
        --weir-fst-pop "$pop_a" \
        --weir-fst-pop "$pop_b" \
        --fst-window-size 100000 \
        --fst-window-step 50000 \
        --out "$window_prefix"
    cmp \
        "$golden/v010-window-fst.windowed.weir.fst" \
        "$window_prefix.windowed.weir.fst"
done

echo "Validating LD and PCA"
for threads in 8 16; do
    ld_prefix="$work/ld-t$threads"
    pca_prefix="$work/pca-t$threads"
    timed "v010-ng${threads}-geno-r2" \
        "$ng_binary" \
        --bcf "$bcf" \
        --threads "$threads" \
        --thin 50000 \
        --geno-r2 \
        --ld-window 20 \
        --ld-window-bp 500000 \
        --min-r2 0.01 \
        --out "$ld_prefix"
    cmp "$golden/v010-geno-r2.geno.ld" "$ld_prefix.geno.ld"
    timed "v010-ng${threads}-pca" \
        "$ng_binary" \
        --bcf "$bcf" \
        --threads "$threads" \
        --keep "$pop_a" \
        --max-missing 1 \
        --thin 100000 \
        --pca \
        --out "$pca_prefix"
    cmp "$golden/v010-pca.pca" "$pca_prefix.pca"
done

if [[ -n "$pca_original" ]]; then
    timed "v010-original-pca" \
        "$pca_original" \
        --gzvcf "$vcf" \
        --keep "$pop_a" \
        --max-missing 1 \
        --thin 100000 \
        --pca \
        --out "$work/original-pca"
    cmp "$golden/v010-pca.pca" "$work/original-pca.pca"
else
    echo "VCFTOOLS_PCA_ORACLE is unset; using the stored 0.1.17+LAPACK golden"
fi

echo "Benchmarking exact full BCF conversion"
timed "v010-original-recode-bcf" \
    "$original" \
    --bcf "$bcf" \
    --recode-bcf \
    --recode-INFO-all \
    --out "$work/original-bcf"
for threads in 8 16; do
    prefix="$work/ng-bcf-t$threads"
    timed "v010-ng${threads}-recode-bcf" \
        "$ng_binary" \
        --bcf "$bcf" \
        --threads "$threads" \
        --recode-bcf \
        --recode-INFO-all \
        --out "$prefix"
    cmp "$work/original-bcf.recode.bcf" "$prefix.recode.bcf"
    rm -f -- "$prefix.recode.bcf"
done
sha256sum \
    "$work/original-bcf.recode.bcf" \
    >"$results/v010-original-recode-bcf.sha256"
rm -f -- "$work/original-bcf.recode.bcf"

echo "Benchmarking two-file site discordance"
timed "v010-original-diff-site-discordance" \
    "$original" \
    --gzvcf "$vcf" \
    --gzdiff "$vcf" \
    --diff-site-discordance \
    --out "$work/original-diff"
cmp \
    "$golden/v010-diff-site-discordance.diff.sites" \
    "$work/original-diff.diff.sites"
for threads in 8 16; do
    prefix="$work/ng-diff-t$threads"
    timed "v010-ng${threads}-diff-site-discordance" \
        "$ng_binary" \
        --bcf "$bcf" \
        --diff-bcf "$bcf" \
        --threads "$threads" \
        --diff-site-discordance \
        --out "$prefix"
    cmp "$work/original-diff.diff.sites" "$prefix.diff.sites"
done

echo "v0.10 advanced compatibility and 8/16-thread benchmarks: PASS"
