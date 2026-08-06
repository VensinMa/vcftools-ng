#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 /path/to/vcftools-ng /path/to/profile.sh" >&2
    exit 2
fi

ng=$1
profile=$2
if [[ ! -x $ng ]]; then
    echo "vcftools-ng executable is not executable: $ng" >&2
    exit 2
fi
if [[ ! -f $profile ]]; then
    echo "Workload profile does not exist: $profile" >&2
    exit 2
fi

# shellcheck source=/dev/null
source "$profile"

required=(
    INPUT_OPTION INPUT_VCF ORACLE_ROOT RESULT_ROOT
    POSITION_SPARSE_SORTED POSITION_SPARSE_UNSORTED
    POSITION_DENSE POSITION_EXCLUDE_DENSE KEEP_25 KEEP_50
    POP_SMALL_1 POP_SMALL_2 POP_LARGE_1 POP_LARGE_2
    PI_WINDOW_SIZE PI_WINDOW_STEP TAJIMA_WINDOW_SIZE
    FST_WINDOW_SIZE FST_WINDOW_STEP
)
for name in "${required[@]}"; do
    if [[ -z ${!name:-} ]]; then
        echo "Profile variable is required: $name" >&2
        exit 2
    fi
done

case $INPUT_OPTION in
    --vcf|--gzvcf|--bcf|--input) ;;
    *)
        echo "Unsupported INPUT_OPTION: $INPUT_OPTION" >&2
        exit 2
        ;;
esac

inputs=(
    "$INPUT_VCF" "$POSITION_SPARSE_SORTED"
    "$POSITION_SPARSE_UNSORTED" "$POSITION_DENSE"
    "$POSITION_EXCLUDE_DENSE" "$KEEP_25" "$KEEP_50"
    "$POP_SMALL_1" "$POP_SMALL_2" "$POP_LARGE_1" "$POP_LARGE_2"
)
for path in "${inputs[@]}"; do
    if [[ ! -f $path ]]; then
        echo "Locked input/fixture does not exist: $path" >&2
        exit 2
    fi
done
if [[ ! -d $ORACLE_ROOT ]]; then
    echo "Hash-locked Original oracle directory is missing: $ORACLE_ROOT" >&2
    echo "The runner never regenerates Original oracles automatically." >&2
    exit 2
fi
if [[ ! -f $ORACLE_ROOT/SHA256SUMS ]]; then
    echo "Original oracle hash manifest is missing: $ORACLE_ROOT/SHA256SUMS" >&2
    exit 2
fi
(
    cd "$ORACLE_ROOT"
    sha256sum -c SHA256SUMS
)

BCFTOOLS=${BCFTOOLS:-bcftools}

if ! declare -p THREADS >/dev/null 2>&1; then
    THREADS=(1 4 8 16 32)
fi
REPEATS=${REPEATS:-3}
KEEP_REPEAT_OUTPUTS=${KEEP_REPEAT_OUTPUTS:-0}
if ((REPEATS < 1)); then
    echo "REPEATS must be at least one" >&2
    exit 2
fi

all_cases=(
    W01_counts
    W02_seven_filter_counts
    W03_positions_sparse_sorted
    W03_positions_sparse_unsorted
    W03_positions_dense
    W04_exclude_positions_dense
    W05_keep_25_counts
    W05_keep_50_counts
    W06_keep_50_bgzf_recode
    W07_window_pi_overlap
    W07_window_pi_nonoverlap
    W08_tajima
    W09_site_fst_small
    W09_site_fst_large
    W10_window_fst
    W10_window_fst_multiallelic
    W11_plain_recode
    W11_bgzf_recode
    W11_bcf_recode
    W12_multiallelic_counts
)
if [[ -n ${RUN_CASES:-} ]]; then
    read -r -a cases <<<"$RUN_CASES"
else
    cases=("${all_cases[@]}")
fi

seven_filters=(
    --min-alleles 2 --max-alleles 2
    --minGQ 10 --minQ 30 --min-meanDP 7
    --max-missing 0.9 --maf 0.1
)

mkdir -p "$RESULT_ROOT"/{logs,outputs,times}
sha256sum "${inputs[@]}" >"$RESULT_ROOT/input-fixtures.sha256"
summary="$RESULT_ROOT/runs.tsv"
if [[ ! -e $summary ]]; then
    printf '%s\n' \
        $'case\tthreads\trepeat\tgate\twall_s\tuser_s\tsystem_s\tcpu\tpeak_rss_kib\tkernel\tcomponents\tbackend\tkept_sites\ttotal_sites\toutput_rows\toutput_bytes\tsha256' \
        >"$summary"
fi

configure_case() {
    local case_id=$1
    args=()
    suffix=
    oracle_suffix=
    comparison=cmp
    case "$case_id" in
        W01_counts)
            args=(--counts); suffix=.frq.count ;;
        W02_seven_filter_counts)
            args=("${seven_filters[@]}" --counts); suffix=.frq.count ;;
        W03_positions_sparse_sorted)
            args=(--positions "$POSITION_SPARSE_SORTED" --counts)
            suffix=.frq.count ;;
        W03_positions_sparse_unsorted)
            args=(--positions "$POSITION_SPARSE_UNSORTED" --counts)
            suffix=.frq.count ;;
        W03_positions_dense)
            args=(--positions "$POSITION_DENSE" --counts)
            suffix=.frq.count ;;
        W04_exclude_positions_dense)
            args=(--exclude-positions "$POSITION_EXCLUDE_DENSE" --counts)
            suffix=.frq.count ;;
        W05_keep_25_counts)
            args=(--keep "$KEEP_25" --counts); suffix=.frq.count ;;
        W05_keep_50_counts)
            args=(--keep "$KEEP_50" --counts); suffix=.frq.count ;;
        W06_keep_50_bgzf_recode)
            args=(--keep "$KEEP_50" "${seven_filters[@]}"
                  --recode-vcf-gz --recode-INFO-all)
            suffix=.recode.vcf.gz; oracle_suffix=.recode.vcf
            comparison=gzip-cmp ;;
        W07_window_pi_overlap)
            args=(--window-pi "$PI_WINDOW_SIZE"
                  --window-pi-step "$PI_WINDOW_STEP")
            suffix=.windowed.pi ;;
        W07_window_pi_nonoverlap)
            args=(--window-pi "$PI_WINDOW_SIZE"
                  --window-pi-step "$PI_WINDOW_SIZE")
            suffix=.windowed.pi ;;
        W08_tajima)
            args=(--TajimaD "$TAJIMA_WINDOW_SIZE")
            suffix=.Tajima.D ;;
        W09_site_fst_small)
            args=(--weir-fst-pop "$POP_SMALL_1"
                  --weir-fst-pop "$POP_SMALL_2")
            suffix=.weir.fst ;;
        W09_site_fst_large)
            args=(--weir-fst-pop "$POP_LARGE_1"
                  --weir-fst-pop "$POP_LARGE_2")
            suffix=.weir.fst ;;
        W10_window_fst)
            args=(--min-alleles 2 --max-alleles 2
                  --weir-fst-pop "$POP_LARGE_1"
                  --weir-fst-pop "$POP_LARGE_2"
                  --fst-window-size "$FST_WINDOW_SIZE"
                  --fst-window-step "$FST_WINDOW_STEP")
            suffix=.windowed.weir.fst ;;
        W10_window_fst_multiallelic)
            args=(--weir-fst-pop "$POP_LARGE_1"
                  --weir-fst-pop "$POP_LARGE_2"
                  --fst-window-size "$FST_WINDOW_SIZE"
                  --fst-window-step "$FST_WINDOW_STEP")
            suffix=.windowed.weir.fst ;;
        W11_plain_recode)
            args=("${seven_filters[@]}" --recode-vcf --recode-INFO-all)
            suffix=.recode.vcf ;;
        W11_bgzf_recode)
            args=("${seven_filters[@]}" --recode-vcf-gz --recode-INFO-all)
            suffix=.recode.vcf.gz; oracle_suffix=.recode.vcf
            comparison=gzip-cmp ;;
        W11_bcf_recode)
            args=(--min-alleles 2 --max-alleles 2 --minQ 30
                  --recode-bcf --recode-INFO-all)
            suffix=.recode.bcf ;;
        W12_multiallelic_counts)
            args=(--min-alleles 3 --counts); suffix=.frq.count ;;
        *)
            echo "Unknown workload case: $case_id" >&2
            exit 2
            ;;
    esac
    if [[ -z $oracle_suffix ]]; then
        oracle_suffix=$suffix
    fi
}

validate_output() {
    local candidate=$1 oracle=$2 mode=$3 scratch=$4
    case $mode in
        cmp)
            cmp "$oracle" "$candidate"
            ;;
        gzip-cmp)
            gzip -dc -- "$candidate" >"$scratch"
            cmp "$oracle" "$scratch"
            ;;
        *)
            echo "Unknown comparison mode: $mode" >&2
            exit 2
            ;;
    esac
}

for case_id in "${cases[@]}"; do
    configure_case "$case_id"
    oracle="$ORACLE_ROOT/$case_id$oracle_suffix"
    if [[ ! -f $oracle ]]; then
        echo "Missing locked Original oracle: $oracle" >&2
        exit 2
    fi
    for threads in "${THREADS[@]}"; do
        if [[ ! $threads =~ ^[1-9][0-9]*$ ]]; then
            echo "Invalid thread count: $threads" >&2
            exit 2
        fi
        for ((repeat = 1; repeat <= REPEATS; ++repeat)); do
            stem="$case_id-t$threads-r$repeat"
            prefix="$RESULT_ROOT/outputs/$stem"
            log="$RESULT_ROOT/logs/$stem.log"
            time_file="$RESULT_ROOT/times/$stem.time.tsv"
            candidate="$prefix$suffix"
            /usr/bin/time -f $'%e\t%U\t%S\t%P\t%M' -o "$time_file" \
                "$ng" "$INPUT_OPTION" "$INPUT_VCF" \
                --threads "$threads" "${args[@]}" \
                --log-file "$log" --out "$prefix" \
                >/dev/null 2>"$RESULT_ROOT/logs/$stem.stderr"

            gate=REUSED
            if ((repeat == 1)); then
                validate_output "$candidate" "$oracle" "$comparison" \
                    "$RESULT_ROOT/outputs/$stem.decompressed"
                gate=PASS
            fi
            read -r wall user system cpu rss <"$time_file"
            kernel=$(sed -n 's/^Execution kernel: //p' "$log" | head -1)
            components=$(
                sed -n 's/^Execution components: //p' "$log" | head -1
            )
            backend=$(
                sed -n 's/^Input backend: \([^ ]*\).*/\1/p' "$log" |
                    head -1
            )
            kept=$(sed -n 's/^Kept sites: //p' "$log" | tail -1)
            total=$(sed -n 's/^Total sites: //p' "$log" | tail -1)
            case $suffix in
                *.vcf.gz)
                    rows=$(gzip -dc -- "$candidate" | wc -l)
                    ;;
                *.bcf)
                    if ! command -v "$BCFTOOLS" >/dev/null 2>&1; then
                        echo "bcftools is required to count BCF output rows" >&2
                        exit 2
                    fi
                    rows=$("$BCFTOOLS" view -H "$candidate" | wc -l)
                    ;;
                *)
                    rows=$(wc -l <"$candidate")
                    ;;
            esac
            bytes=$(stat -c %s "$candidate")
            hash=$(sha256sum "$candidate" | cut -d' ' -f1)
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$case_id" "$threads" "$repeat" "$gate" \
                "$wall" "$user" "$system" "$cpu" "$rss" \
                "$kernel" "$components" "$backend" "$kept" "$total" \
                "$rows" "$bytes" "$hash" >>"$summary"

            rm -f -- "$RESULT_ROOT/outputs/$stem.decompressed"
            if ((repeat > 1 && KEEP_REPEAT_OUTPUTS == 0)); then
                rm -f -- "$candidate"
            fi
        done
    done
done

echo "Representative workload matrix complete: $summary"
