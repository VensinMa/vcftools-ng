#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ng_binary=${1:-"$project_root/build/vcftools-ng"}
fixture=${2:-"$project_root/tests/fixtures/osmanthus412.23chr_100k.bcf"}
bed="$project_root/tests/fixtures/regions.compatibility.bed"
golden="$project_root/tests/golden/subset-non-ref-counts.frq.count"
result_dir="$project_root/benchmarks/results"

if [[ ! -x "$ng_binary" || ! -r "$fixture" || ! -r "$bed" ||
      ! -r "$golden" ]]; then
    echo "Missing executable, fixture, BED file, or golden" >&2
    exit 2
fi

filters=(
    --bed "$bed"
    --keep "$project_root/tests/fixtures/samples.keep.txt"
    --indv W-DA-8
    --remove "$project_root/tests/fixtures/samples.remove.txt"
    --remove-indv W-CP-5
    --non-ref-af 0.1
    --max-non-ref-af 0.8
    --non-ref-af-any 0.2
    --max-non-ref-af-any 0.7
    --non-ref-ac 2
    --max-non-ref-ac 30
    --non-ref-ac-any 5
    --max-non-ref-ac-any 20
)

mkdir -p "$result_dir" "$project_root/tests/output"
work_dir=$(mktemp -d "$project_root/tests/output/non-ref-benchmark.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT

for threads in 8 16; do
    prefix="$work_dir/t$threads"
    timing="$result_dir/subset-non-ref-counts-t$threads.time.txt"
    /usr/bin/time \
        -f $'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\ncpu=%P\nmax_rss_kb=%M' \
        -o "$timing" \
        "$ng_binary" \
        --bcf "$fixture" \
        --threads "$threads" \
        "${filters[@]}" \
        --counts \
        --out "$prefix"
    cmp "$golden" "$prefix.frq.count"
    printf 'threads=%d: PASS (byte-identical non-reference counts)\n' \
        "$threads"
    cat "$timing"
    rm -f -- "$prefix.frq.count"
done
