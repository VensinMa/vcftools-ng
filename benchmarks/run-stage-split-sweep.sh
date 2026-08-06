#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
ng=${1:-"$project_root/build-v0130/vcftools-ng"}
fixture_dir=${2:-/tmp/vcftools-ng-stage-split-230k}
result_dir=${3:-"$project_root/benchmarks/results/v0130-stage-split"}
repeats=${REPEATS:-1}

totals=(4 8 12 16 24 28 32)
filters=(
    --min-alleles 2 --max-alleles 2
    --minGQ 10 --minQ 30 --min-meanDP 7
    --max-missing 0.9 --maf 0.1
    --counts --no-log-file
)

[[ -x "$ng" ]] || {
    echo "error: executable not found: $ng" >&2
    exit 2
}
for path in \
    "$fixture_dir/input.vcf" \
    "$fixture_dir/input.vcf.gz" \
    "$fixture_dir/input.vcf.gz.tbi" \
    "$fixture_dir/input.bcf"; do
    [[ -r "$path" ]] || {
        echo "error: stage-split fixture is missing: $path" >&2
        exit 2
    }
done

mkdir -p "$result_dir"
scratch=$(mktemp -d "${TMPDIR:-/tmp}/vcftools-ng-stage-split.XXXXXX")
cleanup() {
    rm -rf -- "$scratch"
}
trap cleanup EXIT INT TERM

summary="$result_dir/summary.tsv"
printf '%s\n' \
    $'scenario\ttotal_threads\tinput_threads\thts_io_threads\tcompute_threads\trepeat\twall_s\tuser_s\tsystem_s\tcpu_pct\tmax_rss_kb\tplanned_shards\toutput_sha256\texact' \
    >"$summary"

declare -A reference_hash=()

unique_candidates() {
    local total=$1
    shift
    local -A seen=()
    local value
    for value in "$@"; do
        ((value < 1)) && value=1
        ((value >= total)) && value=$((total - 1))
        [[ -n ${seen[$value]:-} ]] && continue
        seen[$value]=1
        printf '%s\n' "$value"
    done | sort -n
}

run_case() {
    local scenario=$1
    local total=$2
    local input_threads=$3
    local hts_io_threads=$4
    local compute_threads=$5
    local repeat=$6
    local prefix="$scratch/${scenario}-t${total}-i${input_threads}-h${hts_io_threads}-c${compute_threads}-r${repeat}"
    local -a input_args=()
    local -a environment=()

    case "$scenario" in
        plain)
            input_args=(--vcf "$fixture_dir/input.vcf" --input-backend plain)
            environment=(
                "VCFTOOLS_NG_TEST_INPUT_THREADS=$input_threads"
                "VCFTOOLS_NG_TEST_COMPUTE_THREADS=$compute_threads"
            )
            ;;
        bgzf_tbi)
            input_args=(--gzvcf "$fixture_dir/input.vcf.gz")
            environment=(
                "VCFTOOLS_NG_TEST_INPUT_THREADS=$input_threads"
                "VCFTOOLS_NG_TEST_COMPUTE_THREADS=$compute_threads"
            )
            ;;
        bcf_stream)
            input_args=(--bcf "$fixture_dir/input.bcf" --input-backend stream)
            environment=(
                "VCFTOOLS_NG_TEST_HTS_IO_THREADS=$hts_io_threads"
                "VCFTOOLS_NG_TEST_COMPUTE_THREADS=$compute_threads"
            )
            ;;
        *)
            echo "error: unknown scenario: $scenario" >&2
            exit 2
            ;;
    esac

    env "${environment[@]}" \
        /usr/bin/time -f $'%e\t%U\t%S\t%P\t%M' \
        -o "$prefix.time" \
        "$ng" "${input_args[@]}" --threads "$total" \
        "${filters[@]}" --out "$prefix" \
        >"$prefix.stdout" 2>"$prefix.stderr"

    local output="$prefix.frq.count"
    local hash
    hash=$(sha256sum "$output" | awk '{print $1}')
    local exact=PASS
    if [[ -z ${reference_hash[$scenario]:-} ]]; then
        reference_hash[$scenario]=$hash
    elif [[ ${reference_hash[$scenario]} != "$hash" ]]; then
        exact=FAIL
    fi
    local planned_shards
    planned_shards=$(awk '/^Planned input shards:/{print $4}' "$prefix.stderr")
    local wall user system cpu rss
    IFS=$'\t' read -r wall user system cpu rss <"$prefix.time"
    cpu=${cpu%%%}
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$scenario" "$total" "$input_threads" "$hts_io_threads" \
        "$compute_threads" "$repeat" "$wall" "$user" "$system" \
        "$cpu" "$rss" "$planned_shards" "$hash" "$exact" \
        >>"$summary"
    rm -f -- "$output" "$prefix.stdout" "$prefix.stderr" "$prefix.time"
}

for total in "${totals[@]}"; do
    mapfile -t text_inputs < <(
        unique_candidates "$total" \
            "$((total / 4))" "$((3 * total / 8))" \
            "$((total / 2))" "$((5 * total / 8))" \
            "$((3 * total / 4))" "$((4 * total / 5))" \
            "$((7 * total / 8))"
    )
    for input_threads in "${text_inputs[@]}"; do
        compute_threads=$((total - input_threads))
        for repeat in $(seq 1 "$repeats"); do
            run_case \
                plain "$total" "$input_threads" 0 \
                "$compute_threads" "$repeat"
            run_case \
                bgzf_tbi "$total" "$input_threads" 0 \
                "$compute_threads" "$repeat"
        done
    done

    mapfile -t hts_io_candidates < <(
        unique_candidates "$total" \
            1 "$((total / 16))" "$((total / 8))" \
            "$((total / 4))" "$((3 * total / 8))" \
            "$((total / 2))"
    )
    for hts_io_threads in "${hts_io_candidates[@]}"; do
        compute_threads=$((total - hts_io_threads))
        for repeat in $(seq 1 "$repeats"); do
            run_case \
                bcf_stream "$total" 0 "$hts_io_threads" \
                "$compute_threads" "$repeat"
        done
    done
done

best_body="$result_dir/.best.body.tsv"
awk -F '\t' '
    NR == 1 { next }
    {
        key = $1 SUBSEP $2
    }
    !(key in best) || $7 + 0 < best[key] {
        best[key] = $7 + 0
        line[key] = $0
        scenarios[$1] = 1
        totals[$2] = 1
    }
    END {
        for (scenario in scenarios) {
            for (total in totals) {
                key = scenario SUBSEP total
                if (!(key in line)) continue
                split(line[key], fields, "\t")
                print fields[1] "\t" fields[2] "\t" fields[3] "\t" fields[4] "\t" fields[5] "\t" fields[7]
            }
        }
    }
' "$summary" | sort -t $'\t' -k1,1 -k2,2n >"$best_body"
{
    printf '%s\n' \
        $'scenario\ttotal_threads\tinput_threads\thts_io_threads\tcompute_threads\twall_s'
    cat "$best_body"
} >"$result_dir/best.tsv"
rm -f -- "$best_body"

if awk -F '\t' 'NR > 1 && $14 != "PASS" { exit 1 }' "$summary"; then
    printf 'STAGE_SPLIT_SWEEP_PASS\n'
else
    echo 'error: stage split changed scientific output bytes' >&2
    exit 1
fi
