#!/usr/bin/env bash
set -euo pipefail

source_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

NG=${NG:-"$source_root/build/vcftools-ng"}
INPUT_BGZF=${INPUT_BGZF:-"$source_root/tests/fixtures/osmanthus412.23chr_100k.vcf.gz"}
INPUT_PLAIN=${INPUT_PLAIN:-"$source_root/data/osmanthus412.subset.vcf"}
THREAD_LIST=${THREAD_LIST:-"8 16"}
RESULT_ROOT=${1:-"$source_root/benchmarks/results/v0114-fused-real-gate"}

golden="$source_root/tests/golden"
work=$(mktemp -d "${TMPDIR:-/tmp}/vcftools-ng-fused-real.XXXXXX")
cleanup() {
    rm -rf -- "$work"
}
trap cleanup EXIT INT TERM

[[ -x "$NG" ]] || {
    printf 'vcftools-ng executable is missing: %s\n' "$NG" >&2
    exit 2
}
for path in \
    "$INPUT_BGZF" "$INPUT_BGZF.tbi" "$INPUT_PLAIN" \
    "$golden/subset-freq.frq" \
    "$golden/subset-freq2.frq" \
    "$golden/subset-counts.frq.count" \
    "$golden/subset-missing-site.lmiss" \
    "$golden/subset-site-depth.ldepth" \
    "$golden/subset-site-mean-depth.ldepth.mean" \
    "$golden/v010-site-quality.lqual"
do
    [[ -s "$path" ]] || {
        printf 'Required input or golden is missing: %s\n' "$path" >&2
        exit 2
    }
done

mkdir -p "$RESULT_ROOT"
summary="$RESULT_ROOT/summary.tsv"
printf 'workload\tinput\tthreads\toriginal_wall_s\tng_wall_s\tspeedup\tcpu_pct\tmax_rss_kb\texact\n' >"$summary"

run_case() {
    local workload=$1
    local input_kind=$2
    local input_path=$3
    local threads=$4
    local original_wall=$5
    local prefix="$work/${workload}-${input_kind}-t${threads}"
    local timing="$work/${workload}-${input_kind}-t${threads}.time"
    local input_option=--vcf
    if [[ "$input_kind" == bgzf_tbi ]]; then
        input_option=--gzvcf
    fi

    local output_args=()
    if [[ "$workload" == combined_six ]]; then
        output_args=(
            --freq --counts --missing-site --site-depth
            --site-mean-depth --site-quality
        )
    else
        output_args=(--freq2)
    fi

    /usr/bin/time -f '%e\t%U\t%S\t%P\t%M' -o "$timing" \
        "$NG" "$input_option" "$input_path" \
        --threads "$threads" \
        "${output_args[@]}" \
        --out "$prefix" \
        >"$work/${workload}-${input_kind}-t${threads}.stdout" \
        2>"$work/${workload}-${input_kind}-t${threads}.stderr"

    if [[ "$workload" == combined_six ]]; then
        cmp "$golden/subset-freq.frq" "$prefix.frq"
        cmp "$golden/subset-counts.frq.count" "$prefix.frq.count"
        cmp "$golden/subset-missing-site.lmiss" "$prefix.lmiss"
        cmp "$golden/subset-site-depth.ldepth" "$prefix.ldepth"
        cmp \
            "$golden/subset-site-mean-depth.ldepth.mean" \
            "$prefix.ldepth.mean"
        cmp "$golden/v010-site-quality.lqual" "$prefix.lqual"
        rm -f -- \
            "$prefix.frq" "$prefix.frq.count" "$prefix.lmiss" \
            "$prefix.ldepth" "$prefix.ldepth.mean" "$prefix.lqual"
    else
        cmp "$golden/subset-freq2.frq" "$prefix.frq"
        rm -f -- "$prefix.frq"
    fi

    local wall user system cpu rss
    IFS=$'\t' read -r wall user system cpu rss <"$timing"
    cpu=${cpu%\%}
    local speedup=NA
    if [[ "$original_wall" != NA ]]; then
        speedup=$(awk -v original="$original_wall" -v candidate="$wall" \
            'BEGIN {printf "%.4f", original / candidate}')
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
        "$workload" "$input_kind" "$threads" "$original_wall" \
        "$wall" "$speedup" "$cpu" "$rss" >>"$summary"
    printf 'PASS %s %s t%s wall=%ss speedup=%s%s\n' \
        "$workload" "$input_kind" "$threads" "$wall" "$speedup" \
        "$([[ "$speedup" == NA ]] && printf '' || printf 'x')"
}

for input_kind in bgzf_tbi plain_vcf; do
    input_path=$INPUT_BGZF
    if [[ "$input_kind" == plain_vcf ]]; then
        input_path=$INPUT_PLAIN
    fi
    for threads in $THREAD_LIST; do
        if [[ "$input_kind" == bgzf_tbi ]]; then
            run_case \
                combined_six "$input_kind" "$input_path" "$threads" 393.74
            run_case freq2 "$input_kind" "$input_path" "$threads" 68.16
        else
            run_case combined_six "$input_kind" "$input_path" "$threads" NA
            run_case freq2 "$input_kind" "$input_path" "$threads" NA
        fi
    done
done

printf 'FUSED REAL GATE PASS %s\n' "$RESULT_ROOT"
