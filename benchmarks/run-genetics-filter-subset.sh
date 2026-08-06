#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ng_binary=${1:-"$project_root/build/vcftools-ng"}
fixture=${2:-"$project_root/tests/fixtures/osmanthus412.23chr_100k.bcf"}
golden="$project_root/tests/golden/subset-mac-hwe-recode.recode.vcf"
result_dir="$project_root/benchmarks/results"

if [[ ! -x "$ng_binary" || ! -r "$fixture" || ! -r "$golden" ]]; then
    echo "Missing executable, 2,300,000-record BCF fixture, or golden" >&2
    exit 2
fi

filters=(
    --keep "$project_root/tests/fixtures/samples.keep.txt"
    --indv W-DA-8
    --remove "$project_root/tests/fixtures/samples.remove.txt"
    --remove-indv W-CP-5
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
    --mac 4
    --max-mac 12
    --hwe 0.001
)

mkdir -p "$result_dir" "$project_root/tests/output"
work_dir=$(mktemp -d "$project_root/tests/output/genetics-filter.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT

for threads in 8 16; do
    prefix="$work_dir/t$threads"
    timing="$result_dir/subset-mac-hwe-recode-t$threads.time.txt"
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
    printf 'threads=%d: PASS (byte-identical recode)\n' "$threads"
    cat "$timing"
    rm -f -- "$prefix.recode.vcf"
done
