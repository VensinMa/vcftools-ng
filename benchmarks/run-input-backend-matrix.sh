#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  NG=/path/to/vcftools-ng \
  ORIGINAL=/path/to/vcftools \
  SUBSET_BGZF=/path/to/subset.vcf.gz \
  SUBSET_PLAIN=/path/to/subset.vcf \
  SUBSET_BCF=/path/to/subset.bcf \
  ./benchmarks/run-input-backend-matrix.sh [result-directory]

Development defaults:
  THREAD_LIST="8 16"
  REPEATS=1
  WORKLOAD=counts
  EXPECTED_RECORDS=2300000
  CACHE_STATE=unspecified

Final scaling example:
  THREAD_LIST="1 2 4 8 16 32 64 128 256" REPEATS=3 ...

The BGZF and BCF source paths must have their normal TBI/CSI sidecars.
The script creates temporary symlink paths without sidecars for the no-index
cases; it never removes or changes the source indices.
EOF
}

if [[ ${1:-} == "--help" || ${1:-} == "-h" ]]; then
    usage
    exit 0
fi

: "${NG:?Set NG to the vcftools-ng executable}"
: "${ORIGINAL:?Set ORIGINAL to the VCFtools 0.1.17 executable}"
: "${SUBSET_BGZF:?Set SUBSET_BGZF to the real-subset BGZF VCF}"
: "${SUBSET_PLAIN:?Set SUBSET_PLAIN to the real-subset plain VCF}"
: "${SUBSET_BCF:?Set SUBSET_BCF to the real-subset BCF}"

THREAD_LIST=${THREAD_LIST:-"8 16"}
REPEATS=${REPEATS:-1}
WORKLOAD=${WORKLOAD:-counts}
EXPECTED_RECORDS=${EXPECTED_RECORDS:-2300000}
CACHE_STATE=${CACHE_STATE:-unspecified}
RESULT_ROOT=${1:-benchmarks/results/input-backends}

for executable in "$NG" "$ORIGINAL"; do
    if [[ ! -x "$executable" ]]; then
        printf 'Executable is missing or not executable: %s\n' "$executable" >&2
        exit 2
    fi
done
for input in "$SUBSET_BGZF" "$SUBSET_PLAIN" "$SUBSET_BCF"; do
    if [[ ! -f "$input" ]]; then
        printf 'Input is missing: %s\n' "$input" >&2
        exit 2
    fi
done
if [[ ! -f "${SUBSET_BGZF}.tbi" && ! -f "${SUBSET_BGZF}.csi" ]]; then
    printf 'BGZF indexed case requires %s.tbi or %s.csi\n' \
        "$SUBSET_BGZF" "$SUBSET_BGZF" >&2
    exit 2
fi
if [[ ! -f "${SUBSET_BCF}.csi" ]]; then
    printf 'BCF indexed case requires %s.csi\n' "$SUBSET_BCF" >&2
    exit 2
fi
if [[ ! "$REPEATS" =~ ^[1-9][0-9]*$ ||
      ! "$EXPECTED_RECORDS" =~ ^[1-9][0-9]*$ ]]; then
    printf 'REPEATS and EXPECTED_RECORDS must be positive integers\n' >&2
    exit 2
fi

case "$WORKLOAD" in
    counts)
        original_args=(--counts)
        ng_args=(--counts)
        required_suffixes=(.frq.count)
        ;;
    five-stats)
        original_args=(--freq --counts --missing-site --site-depth --site-mean-depth)
        ng_args=(--freq --counts --missing-site --site-depth --site-mean-depth)
        required_suffixes=(.frq .frq.count .lmiss .ldepth .ldepth.mean)
        ;;
    *)
        printf 'Unsupported WORKLOAD: %s (use counts or five-stats)\n' \
            "$WORKLOAD" >&2
        exit 2
        ;;
esac

NG=$(realpath "$NG")
ORIGINAL=$(realpath "$ORIGINAL")
SUBSET_BGZF=$(realpath "$SUBSET_BGZF")
SUBSET_PLAIN=$(realpath "$SUBSET_PLAIN")
SUBSET_BCF=$(realpath "$SUBSET_BCF")
mkdir -p "$RESULT_ROOT"
RESULT_ROOT=$(realpath "$RESULT_ROOT")
mkdir -p "$RESULT_ROOT/logs"
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/vcftools-ng-input-matrix.XXXXXX")
cleanup() {
    if [[ -n ${work_dir:-} && -d "$work_dir" ]]; then
        rm -rf -- "$work_dir"
    fi
}
trap cleanup EXIT INT TERM

ln -s "$SUBSET_BGZF" "$work_dir/no-index.vcf.gz"
ln -s "$SUBSET_BCF" "$work_dir/no-index.bcf"

manifest="$RESULT_ROOT/manifest.tsv"
timings="$RESULT_ROOT/timings.tsv"
comparisons="$RESULT_ROOT/comparisons.tsv"
printf 'case\tkind\tindexed\tpath\tbytes\texpected_records\n' >"$manifest"
printf 'case\tengine\tthreads\trepeat\twall_s\tuser_s\tsystem_s\tcpu_pct\tmax_rss_kb\texit\tstdout\tstderr\n' >"$timings"
printf 'case\tthreads\trepeat\tartifact\tstatus\tgolden_sha256\tcandidate_sha256\n' >"$comparisons"

declare -a case_names=(
    bgzf_tbi
    bgzf_no_index
    plain_vcf
    bcf_csi
    bcf_no_index
)
declare -a case_kinds=(gzvcf gzvcf vcf bcf bcf)
declare -a case_indexed=(yes no no yes no)
declare -a case_paths=(
    "$SUBSET_BGZF"
    "$work_dir/no-index.vcf.gz"
    "$SUBSET_PLAIN"
    "$SUBSET_BCF"
    "$work_dir/no-index.bcf"
)

run_timed() {
    local case_name=$1
    local engine=$2
    local threads=$3
    local repeat=$4
    local output_prefix=$5
    shift 5
    local timing_file="$work_dir/time-${case_name}-${engine}-${threads}-${repeat}.txt"
    local log_base="$RESULT_ROOT/logs/${case_name}-${engine}-t${threads}-r${repeat}"
    local status=0
    /usr/bin/time \
        -f '%e\t%U\t%S\t%P\t%M' \
        -o "$timing_file" \
        "$@" >"${log_base}.stdout.txt" 2>"${log_base}.stderr.txt" ||
        status=$?
    local values
    values=$(<"$timing_file")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$case_name" "$engine" "$threads" "$repeat" "$values" "$status" \
        "${log_base}.stdout.txt" "${log_base}.stderr.txt" \
        >>"$timings"
    if ((status != 0)); then
        printf 'Benchmark failed: case=%s engine=%s threads=%s repeat=%s\n' \
            "$case_name" "$engine" "$threads" "$repeat" >&2
        exit "$status"
    fi
    for suffix in "${required_suffixes[@]}"; do
        if [[ ! -f "${output_prefix}${suffix}" ]]; then
            printf 'Expected output is missing: %s\n' \
                "${output_prefix}${suffix}" >&2
            exit 2
        fi
    done
}

for index in "${!case_names[@]}"; do
    case_name=${case_names[$index]}
    kind=${case_kinds[$index]}
    indexed=${case_indexed[$index]}
    input=${case_paths[$index]}
    bytes=$(stat -Lc '%s' "$input")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$case_name" "$kind" "$indexed" "$input" "$bytes" \
        "$EXPECTED_RECORDS" >>"$manifest"

    original_input_option="--${kind}"
    ng_index_args=()
    if [[ "$case_name" == bgzf_no_index ||
          "$case_name" == bcf_no_index ]]; then
        ng_index_args=(--no-auto-index)
    fi
    golden_prefix="$work_dir/${case_name}-original"
    run_timed \
        "$case_name" original 1 1 "$golden_prefix" \
        "$ORIGINAL" "$original_input_option" "$input" \
        "${original_args[@]}" --out "$golden_prefix"

    for threads in $THREAD_LIST; do
        if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
            printf 'Invalid thread count in THREAD_LIST: %s\n' "$threads" >&2
            exit 2
        fi
        for repeat in $(seq 1 "$REPEATS"); do
            candidate_prefix="$work_dir/${case_name}-ng-t${threads}-r${repeat}"
            run_timed \
                "$case_name" vcftools-ng "$threads" "$repeat" \
                "$candidate_prefix" \
                "$NG" "$original_input_option" "$input" \
                --threads "$threads" "${ng_index_args[@]}" \
                "${ng_args[@]}" \
                --out "$candidate_prefix"
            for suffix in "${required_suffixes[@]}"; do
                golden="${golden_prefix}${suffix}"
                candidate="${candidate_prefix}${suffix}"
                if cmp -s "$golden" "$candidate"; then
                    comparison=PASS
                else
                    comparison=FAIL
                fi
                golden_sha256=$(sha256sum "$golden" | cut -d' ' -f1)
                candidate_sha256=$(sha256sum "$candidate" | cut -d' ' -f1)
                printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                    "$case_name" "$threads" "$repeat" "$suffix" \
                    "$comparison" "$golden_sha256" "$candidate_sha256" \
                    >>"$comparisons"
                if [[ "$comparison" != PASS ]]; then
                    printf 'Byte comparison failed: case=%s threads=%s repeat=%s artifact=%s\n' \
                        "$case_name" "$threads" "$repeat" "$suffix" >&2
                    exit 1
                fi
            done
        done
    done
done

awk -F '\t' '
    BEGIN {
        OFS = "\t"
        print "case", "threads", "repeat", "original_wall_s", \
              "ng_wall_s", "speedup", "exact"
    }
    NR == 1 { next }
    $2 == "original" { baseline[$1] = $5; next }
    $2 == "vcftools-ng" {
        speedup = ($5 > 0 ? baseline[$1] / $5 : 0)
        printf "%s\t%s\t%s\t%s\t%s\t%.4f\tPASS\n", \
               $1, $3, $4, baseline[$1], $5, speedup
    }
' "$timings" >"$RESULT_ROOT/summary.tsv"

{
    printf 'date\t%s\n' "$(date --iso-8601=seconds)"
    printf 'host\t%s\n' "$(hostname)"
    printf 'kernel\t%s\n' "$(uname -sr)"
    printf 'cpu_available\t%s\n' "$(nproc)"
    printf 'thread_list\t%s\n' "$THREAD_LIST"
    printf 'repeats\t%s\n' "$REPEATS"
    printf 'workload\t%s\n' "$WORKLOAD"
    printf 'expected_records\t%s\n' "$EXPECTED_RECORDS"
    printf 'cache_state\t%s\n' "$CACHE_STATE"
    if command -v lscpu >/dev/null 2>&1; then
        printf 'cpu_model\t%s\n' \
            "$(lscpu | awk -F: '/Model name/ {sub(/^[[:space:]]+/, "", $2); print $2; exit}')"
    fi
    if command -v lsblk >/dev/null 2>&1; then
        printf 'block_devices\t%s\n' \
            "$(lsblk -dn -o NAME,ROTA,TYPE | paste -sd ';' -)"
    fi
} >"$RESULT_ROOT/environment.tsv"

printf 'Input-backend benchmark completed: %s\n' "$RESULT_ROOT"
