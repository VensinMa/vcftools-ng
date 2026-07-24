#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ng_binary=${1:-"$project_root/build/vcftools-ng"}
fixture=${2:-"$project_root/tests/fixtures/osmanthus412.23chr_100k.bcf"}
golden="$project_root/tests/golden"
result_dir="$project_root/benchmarks/results"

if [[ ! -x "$ng_binary" || ! -r "$fixture" ]]; then
    echo "Missing executable or 2,300,000-record BCF fixture" >&2
    exit 2
fi

sample_filters=(
    --keep "$project_root/tests/fixtures/samples.keep.txt"
    --indv W-DA-8
    --remove "$project_root/tests/fixtures/samples.remove.txt"
    --remove-indv W-CP-5
)
site_filters=(
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

mkdir -p "$result_dir" "$project_root/tests/output"
work_dir=$(mktemp -d "$project_root/tests/output/sample-benchmark.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT

compare_outputs() {
    local prefix=$1
    cmp "$golden/subset-samples-freq.frq" "$prefix.frq"
    cmp "$golden/subset-samples-counts.frq.count" "$prefix.frq.count"
    cmp "$golden/subset-samples-missing-site.lmiss" "$prefix.lmiss"
    cmp "$golden/subset-samples-site-depth.ldepth" "$prefix.ldepth"
    cmp \
        "$golden/subset-samples-site-mean-depth.ldepth.mean" \
        "$prefix.ldepth.mean"
    cmp \
        "$golden/subset-samples-recode.recode.vcf" \
        "$prefix.recode.vcf"
}

for threads in 8 16; do
    prefix="$work_dir/t$threads"
    timing="$result_dir/subset-samples-all-t$threads.time.txt"
    /usr/bin/time \
        -f $'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\ncpu=%P\nmax_rss_kb=%M' \
        -o "$timing" \
        "$ng_binary" \
        --bcf "$fixture" \
        --threads "$threads" \
        "${sample_filters[@]}" \
        "${site_filters[@]}" \
        --freq \
        --counts \
        --missing-site \
        --site-depth \
        --site-mean-depth \
        --recode \
        --recode-INFO-all \
        --out "$prefix"
    compare_outputs "$prefix"
    printf 'threads=%d: PASS (six byte-identical outputs)\n' "$threads"
    cat "$timing"
    rm -f -- "$prefix.frq" "$prefix.frq.count" "$prefix.lmiss" \
        "$prefix.ldepth" "$prefix.ldepth.mean" "$prefix.recode.vcf"
done
