#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ng_binary=${1:-"$project_root/build/vcftools-ng"}
fixture=${2:-"$project_root/tests/fixtures/osmanthus412.23chr_100k.bcf"}
bed="$project_root/tests/fixtures/regions.compatibility.bed"
golden="$project_root/tests/golden/subset-bed-include-counts.frq.count"
result_dir="$project_root/benchmarks/results"

if [[ ! -x "$ng_binary" || ! -r "$fixture" || ! -r "$bed" ||
      ! -r "$golden" ]]; then
    echo "Missing executable, fixture, BED file, or golden" >&2
    exit 2
fi

mkdir -p "$result_dir" "$project_root/tests/output"
work_dir=$(mktemp -d "$project_root/tests/output/bed-benchmark.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT

for threads in 8 16; do
    prefix="$work_dir/t$threads"
    timing="$result_dir/subset-bed-counts-t$threads.time.txt"
    /usr/bin/time \
        -f $'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\ncpu=%P\nmax_rss_kb=%M' \
        -o "$timing" \
        "$ng_binary" \
        --bcf "$fixture" \
        --threads "$threads" \
        --bed "$bed" \
        --counts \
        --out "$prefix"
    cmp "$golden" "$prefix.frq.count"
    printf 'threads=%d: PASS (byte-identical BED-filtered counts)\n' \
        "$threads"
    cat "$timing"
    rm -f -- "$prefix.frq.count"
done
