#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ng_binary=${1:-"$project_root/build/vcftools-ng"}
fixture=${2:-"$project_root/tests/fixtures/osmanthus412.23chr_100k.bcf"}
golden="$project_root/tests/golden/subset-filtered-bcf-recode-info-all.recode.vcf"
result_dir="$project_root/benchmarks/results"

if [[ ! -x "$ng_binary" ]]; then
    echo "vcftools-ng binary is not executable: $ng_binary" >&2
    exit 2
fi

if [[ ! -r "$fixture" ]]; then
    echo "2,300,000-record BCF fixture is not readable: $fixture" >&2
    exit 2
fi

if [[ ! -r "$golden" ]]; then
    echo "VCFtools 0.1.17 recode golden is not readable: $golden" >&2
    exit 2
fi

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

mkdir -p "$result_dir" "$project_root/tests/output"
work_dir=$(mktemp -d "$project_root/tests/output/recode-benchmark.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT

for threads in 8 16; do
    prefix="$work_dir/t$threads"
    timing="$result_dir/subset-filtered-recode-info-all-t$threads.time.txt"

    /usr/bin/time \
        -f $'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\ncpu=%P\nmax_rss_kb=%M' \
        -o "$timing" \
        "$ng_binary" \
        --bcf "$fixture" \
        --threads "$threads" \
        "${filters[@]}" \
        --recode-vcf \
        --recode-INFO-all \
        --out "$prefix"

    cmp "$golden" "$prefix.recode.vcf"
    printf 'threads=%d: PASS (byte-identical)\n' "$threads"
    cat "$timing"
    rm -f -- "$prefix.recode.vcf"
done
