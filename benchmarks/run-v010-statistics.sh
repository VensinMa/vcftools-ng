#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ng_binary=${NG_BINARY:-"$project_root/build/vcftools-ng"}
original_binary=${VCFTOOLS_ORIGINAL:-/home/vensin/anaconda3/envs/vcftools/bin/vcftools}
vcf_fixture="$project_root/tests/fixtures/osmanthus412.23chr_100k.vcf.gz"
bcf_fixture="$project_root/tests/fixtures/osmanthus412.23chr_100k.bcf"
result_root="$project_root/benchmarks/results"
golden_root="$project_root/tests/golden"
output_root="$project_root/tests/output"

mkdir -p "$result_root" "$output_root"
run_output=$(mktemp -d "$output_root/v010-benchmark.XXXXXX")
trap 'rm -rf -- "$run_output"' EXIT

suffix_for() {
    case "$1" in
        depth) echo idepth ;;
        missing-indv) echo imiss ;;
        het) echo het ;;
        hardy) echo hwe ;;
        site-quality) echo lqual ;;
    esac
}

outputs=(depth missing-indv het hardy site-quality)

for output in "${outputs[@]}"; do
    suffix=$(suffix_for "$output")
    /usr/bin/time -v \
        -o "$result_root/v010-original-$output.time.txt" \
        "$original_binary" \
        --gzvcf "$vcf_fixture" \
        --"$output" \
        --out "$run_output/original-$output" \
        >"$result_root/v010-original-$output.stdout.txt" \
        2>"$result_root/v010-original-$output.log"
    cmp \
        "$golden_root/v010-$output.$suffix" \
        "$run_output/original-$output.$suffix"
    rm -f -- "$run_output/original-$output.$suffix"
done

for threads in 8 16; do
    for output in "${outputs[@]}"; do
        suffix=$(suffix_for "$output")
        /usr/bin/time -v \
            -o "$result_root/v010-ng${threads}-$output.time.txt" \
            "$ng_binary" \
            --bcf "$bcf_fixture" \
            --threads "$threads" \
            --"$output" \
            --out "$run_output/ng${threads}-$output" \
            >"$result_root/v010-ng${threads}-$output.stdout.txt" \
            2>"$result_root/v010-ng${threads}-$output.log"
        cmp \
            "$golden_root/v010-$output.$suffix" \
            "$run_output/ng${threads}-$output.$suffix"
        rm -f -- "$run_output/ng${threads}-$output.$suffix"
    done

    combined="$run_output/ng${threads}-combined"
    /usr/bin/time -v \
        -o "$result_root/v010-ng${threads}-combined.time.txt" \
        "$ng_binary" \
        --bcf "$bcf_fixture" \
        --threads "$threads" \
        --depth \
        --missing-indv \
        --het \
        --hardy \
        --site-quality \
        --out "$combined" \
        >"$result_root/v010-ng${threads}-combined.stdout.txt" \
        2>"$result_root/v010-ng${threads}-combined.log"
    cmp "$golden_root/v010-depth.idepth" "$combined.idepth"
    cmp "$golden_root/v010-missing-indv.imiss" "$combined.imiss"
    cmp "$golden_root/v010-het.het" "$combined.het"
    cmp "$golden_root/v010-hardy.hwe" "$combined.hwe"
    cmp "$golden_root/v010-site-quality.lqual" "$combined.lqual"
    rm -f -- "$combined".*
done

echo "v0.10 statistics: original, 8-thread, and 16-thread outputs are byte-identical."
