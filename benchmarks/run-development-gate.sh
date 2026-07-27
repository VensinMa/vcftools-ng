#!/usr/bin/env bash
set -euo pipefail

source_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

NG=${NG:-"$source_root/build/vcftools-ng"}
BASELINE_ROOT=${BASELINE_ROOT:-"$source_root/benchmarks/results/v0114-real-filter-subset"}
INPUT_BGZF=${INPUT_BGZF:-"$source_root/tests/fixtures/osmanthus412.23chr_100k.vcf.gz"}
INPUT_PLAIN=${INPUT_PLAIN:-"$source_root/data/osmanthus412.subset.vcf"}
INPUT_BCF=${INPUT_BCF:-"$source_root/tests/fixtures/osmanthus412.23chr_100k.bcf"}
THREAD_LIST=${THREAD_LIST:-"1 2 4 8 16 32"}
revision=${REVISION_LABEL:-"$(
    git -C "$source_root" rev-parse --short HEAD
)"}
RESULT_ROOT=${1:-"$source_root/benchmarks/results/development-$revision-$(date +%Y%m%d-%H%M%S)"}

readonly filter_args=(
    --min-alleles 2 --max-alleles 2
    --minGQ 10 --minQ 30
    --min-meanDP 7 --max-missing 0.9 --maf 0.1
)
readonly lock="$BASELINE_ROOT/baseline.lock.tsv"

[[ -x "$NG" ]] || {
    printf 'vcftools-ng executable is missing: %s\n' "$NG" >&2
    exit 2
}
[[ -s "$lock" ]] || {
    printf 'Locked Original baseline is missing: %s\n' "$lock" >&2
    exit 2
}

lock_value() {
    awk -F '\t' -v key="$1" \
        '$1 == key {print $2; exit}' "$lock"
}

validate_artifact() {
    local label=$1
    local path=$2
    local expected_bytes=$3
    local expected_sha256=$4
    [[ -f "$path" ]] || {
        printf 'Baseline artifact is missing: %s (%s)\n' \
            "$label" "$path" >&2
        exit 2
    }
    local actual_bytes
    actual_bytes=$(stat -c '%s' "$path")
    [[ "$actual_bytes" == "$expected_bytes" ]] || {
        printf 'Baseline size mismatch: %s expected=%s actual=%s\n' \
            "$label" "$expected_bytes" "$actual_bytes" >&2
        exit 2
    }
    local actual_sha256
    actual_sha256=$(sha256sum "$path" | cut -d' ' -f1)
    [[ "$actual_sha256" == "$expected_sha256" ]] || {
        printf 'Baseline SHA-256 mismatch: %s\n' "$label" >&2
        exit 2
    }
    printf 'LOCK PASS %-20s %s\n' "$label" "$actual_sha256"
}

INPUT_BGZF=$(realpath "$INPUT_BGZF")
INPUT_PLAIN=$(realpath "$INPUT_PLAIN")
INPUT_BCF=$(realpath "$INPUT_BCF")
BASELINE_ROOT=$(realpath "$BASELINE_ROOT")
NG=$(realpath "$NG")

bgzf_index=
if [[ -f "$INPUT_BGZF.tbi" ]]; then
    bgzf_index="$INPUT_BGZF.tbi"
elif [[ -f "$INPUT_BGZF.csi" ]]; then
    bgzf_index="$INPUT_BGZF.csi"
else
    printf 'BGZF development input has no TBI/CSI\n' >&2
    exit 2
fi
validate_artifact \
    bgzf-input "$INPUT_BGZF" \
    "$(lock_value bgzf_bytes)" "$(lock_value bgzf_sha256)"
validate_artifact \
    bgzf-index "$bgzf_index" \
    "$(lock_value bgzf_index_bytes)" \
    "$(lock_value bgzf_index_sha256)"
validate_artifact \
    plain-input "$INPUT_PLAIN" \
    "$(lock_value plain_bytes)" "$(lock_value plain_sha256)"
validate_artifact \
    bcf-input "$INPUT_BCF" \
    "$(lock_value bcf_bytes)" "$(lock_value bcf_sha256)"
validate_artifact \
    bgzf-golden "$BASELINE_ROOT/golden/gzvcf.vcf" \
    "$(lock_value golden_bgzf_bytes)" \
    "$(lock_value golden_bgzf_sha256)"
validate_artifact \
    plain-golden "$BASELINE_ROOT/golden/vcf.vcf" \
    "$(lock_value golden_plain_bytes)" \
    "$(lock_value golden_plain_sha256)"
validate_artifact \
    bcf-golden "$BASELINE_ROOT/golden/bcf.vcf" \
    "$(lock_value golden_bcf_bytes)" \
    "$(lock_value golden_bcf_sha256)"

mkdir -p "$RESULT_ROOT"/{logs,runs,scratch}
RESULT_ROOT=$(realpath "$RESULT_ROOT")

case_parameters() {
    backend_args=()
    case "$1" in
        bgzf_tbi)
            kind=gzvcf
            input=$INPUT_BGZF
            golden="$BASELINE_ROOT/golden/gzvcf.vcf"
            original_wall=$(lock_value original_bgzf_wall_s)
            ;;
        plain_vcf)
            kind=vcf
            input=$INPUT_PLAIN
            golden="$BASELINE_ROOT/golden/vcf.vcf"
            original_wall=$(lock_value original_plain_wall_s)
            ;;
        bcf_adaptive_stream)
            kind=bcf
            input=$INPUT_BCF
            golden="$BASELINE_ROOT/golden/bcf.vcf"
            original_wall=$(lock_value original_bcf_wall_s)
            ;;
        *)
            printf 'Unknown development scenario: %s\n' "$1" >&2
            exit 2
            ;;
    esac
}

for case_name in bgzf_tbi plain_vcf bcf_adaptive_stream; do
    case_parameters "$case_name"
    for threads in $THREAD_LIST; do
        [[ "$threads" =~ ^[1-9][0-9]*$ ]] || {
            printf 'Invalid thread count: %s\n' "$threads" >&2
            exit 2
        }
        run_id="$case_name-t$threads"
        metadata="$RESULT_ROOT/runs/$run_id.tsv"
        if [[ -s "$metadata" ]] &&
            awk -F '\t' \
                'NR == 2 && $11 == "PASS" {ok=1} END {exit !ok}' \
                "$metadata"; then
            printf 'SKIP completed %s\n' "$run_id"
            continue
        fi
        scratch="$RESULT_ROOT/scratch/$run_id"
        mkdir -p "$scratch"
        output="$scratch/output.vcf"
        timing="$RESULT_ROOT/runs/$run_id.time.txt"
        stderr="$RESULT_ROOT/logs/$run_id.stderr.txt"
        printf 'START %s %s\n' "$run_id" "$(date --iso-8601=seconds)"
        /usr/bin/time -f '%e\t%U\t%S\t%P\t%M' -o "$timing" \
            "$NG" "--$kind" "$input" --threads "$threads" \
            "${backend_args[@]}" \
            "${filter_args[@]}" \
            --log-file "$stderr" \
            --recode --recode-INFO-all --stdout \
            >"$output" 2>"$stderr"
        if [[ "$case_name" == bcf_adaptive_stream ]]; then
            grep -q '^Input backend: stream ' "$stderr" || {
                printf 'BCF adaptive-stream gate selected the wrong backend: %s\n' \
                    "$run_id" >&2
                exit 2
            }
        fi
        cmp "$golden" "$output"
        output_hash=$(sha256sum "$output" | cut -d' ' -f1)
        IFS=$'\t' read -r wall user system cpu rss <"$timing"
        cpu=${cpu%\%}
        printf 'case\tthreads\toriginal_wall_s\twall_s\tuser_s\tsystem_s\tcpu_pct\tmax_rss_kb\toutput_bytes\toutput_sha256\texact\n' \
            >"$metadata"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
            "$case_name" "$threads" "$original_wall" \
            "$wall" "$user" "$system" "$cpu" "$rss" \
            "$(stat -c '%s' "$output")" "$output_hash" \
            >>"$metadata"
        rm -rf -- "$scratch"
        printf 'DONE %s wall=%s exact=PASS\n' "$run_id" "$wall"
    done
done

all_runs="$RESULT_ROOT/all-runs.tsv"
first_metadata=$(
    find "$RESULT_ROOT/runs" -maxdepth 1 -name '*.tsv' -type f -print |
        sort | sed -n '1p'
)
head -1 "$first_metadata" >"$all_runs"
find "$RESULT_ROOT/runs" -maxdepth 1 -name '*.tsv' -type f -print0 |
    sort -z | xargs -0 -r -n1 tail -1 >>"$all_runs"

awk -F '\t' '
    BEGIN {
        OFS = "\t"
        print "case", "threads", "original_wall_s", "ng_wall_s", \
              "speedup", "cpu_pct", "max_rss_kb", "exact"
    }
    NR == 1 {next}
    {
        printf "%s\t%s\t%.2f\t%.2f\t%.4f\t%.1f\t%d\t%s\n", \
               $1, $2, $3, $4, $3 / $4, $7, $8, $11
    }
' "$all_runs" | {
    IFS= read -r header
    printf '%s\n' "$header"
    sort -t $'\t' -k1,1 -k2,2n
} >"$RESULT_ROOT/summary.tsv"

{
    printf 'key\tvalue\n'
    printf 'date\t%s\n' "$(date --iso-8601=seconds)"
    printf 'revision_label\t%s\n' "$revision"
    printf 'git_commit\t%s\n' "$(git -C "$source_root" rev-parse HEAD)"
    printf 'ng_version\t%s\n' "$("$NG" --version)"
    printf 'ng_binary_sha256\t%s\n' \
        "$(sha256sum "$NG" | cut -d' ' -f1)"
    printf 'git_tree_state\t%s\n' "$(
        if [[ -z "$(git -C "$source_root" status --short)" ]]; then
            printf clean
        else
            printf dirty
        fi
    )"
    printf 'baseline_root\t%s\n' "$BASELINE_ROOT"
    printf 'baseline_lock_sha256\t%s\n' \
        "$(sha256sum "$lock" | cut -d' ' -f1)"
    printf 'scenarios\tBGZF VCF + TBI;Plain VCF;BCF adaptive streaming full-scan path\n'
    printf 'bcf_development_policy\tdefault auto policy validates and preserves an existing CSI but must select stream for full-file recode\n'
    printf 'threads\t%s\n' "$THREAD_LIST"
    printf 'workload\t%s\n' "$(lock_value workload)"
    printf 'oracle_policy\tlocked; Original was not rerun\n'
    printf 'comparison\tcmp plus SHA-256\n'
    printf 'cache\toperating-system cache not flushed\n'
} >"$RESULT_ROOT/manifest.tsv"

printf 'DEVELOPMENT GATE PASS %s\n' "$RESULT_ROOT"
