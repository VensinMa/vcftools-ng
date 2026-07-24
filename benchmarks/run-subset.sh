#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ng_binary=${1:-"$project_root/build/vcftools-ng"}
fixture=${2:-"$project_root/tests/fixtures/osmanthus412.23chr_100k.bcf"}
golden="$project_root/tests/golden"
result_dir="$project_root/benchmarks/results"

if [[ ! -x "$ng_binary" ]]; then
    echo "vcftools-ng binary is not executable: $ng_binary" >&2
    exit 2
fi

if [[ ! -r "$fixture" ]]; then
    echo "2,300,000-record BCF fixture is not readable: $fixture" >&2
    exit 2
fi

mkdir -p "$result_dir" "$project_root/tests/output"
work_dir=$(mktemp -d "$project_root/tests/output/benchmark.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT

compare_outputs() {
    local prefix=$1
    cmp "$golden/subset-freq.frq" "$prefix.frq"
    cmp "$golden/subset-counts.frq.count" "$prefix.frq.count"
    cmp "$golden/subset-missing-site.lmiss" "$prefix.lmiss"
    cmp "$golden/subset-site-depth.ldepth" "$prefix.ldepth"
    cmp "$golden/subset-site-mean-depth.ldepth.mean" "$prefix.ldepth.mean"
}

for threads in 8 16; do
    prefix="$work_dir/t$threads"
    timing="$result_dir/subset-bcf-t$threads.time.txt"

    /usr/bin/time \
        -f $'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\ncpu=%P\nmax_rss_kb=%M' \
        -o "$timing" \
        "$ng_binary" \
        --bcf "$fixture" \
        --threads "$threads" \
        --freq \
        --counts \
        --missing-site \
        --site-depth \
        --site-mean-depth \
        --out "$prefix"

    compare_outputs "$prefix"
    printf 'threads=%d: PASS (byte-identical)\n' "$threads"
    cat "$timing"
    rm -f -- "$prefix.frq" "$prefix.frq.count" "$prefix.lmiss" \
        "$prefix.ldepth" "$prefix.ldepth.mean"
done
