#!/usr/bin/env bash
set -euo pipefail

source_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

NG=${NG:-"$source_root/build/vcftools-ng"}
ORIGINAL=${ORIGINAL:-"/home/vensin/anaconda3/pkgs/vcftools-0.1.17-pl5321h077b44d_0/bin/vcftools"}
BCFTOOLS=${BCFTOOLS:-"/home/vensin/software/bcftools/bcftools"}
INPUT_BGZF=${INPUT_BGZF:-"$source_root/tests/fixtures/osmanthus412.23chr_100k.vcf.gz"}
INPUT_PLAIN=${INPUT_PLAIN:-"$source_root/data/osmanthus412.subset.vcf"}
INPUT_BCF=${INPUT_BCF:-"$source_root/tests/fixtures/osmanthus412.23chr_100k.bcf"}
RESULT_ROOT=${1:-"$source_root/benchmarks/results/v0114-real-filter-subset"}
THREAD_LIST=${THREAD_LIST:-"1 2 4 8 16 32"}
EXPECTED_RECORDS=${EXPECTED_RECORDS:-2300000}

readonly filter_args=(
    --min-alleles 2 --max-alleles 2
    --minGQ 10 --minQ 30
    --min-meanDP 7 --max-missing 0.9 --maf 0.1
)
readonly cases=(
    bgzf_tbi
    bgzf_auto_csi
    bgzf_no_auto_index
    plain_vcf
    bcf_csi
    bcf_auto_csi
    bcf_no_auto_index
)

for executable in "$NG" "$ORIGINAL" "$BCFTOOLS"; do
    [[ -x "$executable" ]] || {
        printf 'Missing executable: %s\n' "$executable" >&2
        exit 2
    }
done
for input in "$INPUT_BGZF" "$INPUT_PLAIN" "$INPUT_BCF"; do
    [[ -f "$input" ]] || {
        printf 'Missing input: %s\n' "$input" >&2
        exit 2
    }
done
[[ -f "$INPUT_BGZF.tbi" || -f "$INPUT_BGZF.csi" ]] || {
    printf 'Indexed BGZF scenario requires TBI or CSI\n' >&2
    exit 2
}
[[ -f "$INPUT_BCF.csi" ]] || {
    printf 'Indexed BCF scenario requires CSI\n' >&2
    exit 2
}

NG=$(realpath "$NG")
ORIGINAL=$(realpath "$ORIGINAL")
BCFTOOLS=$(realpath "$BCFTOOLS")
INPUT_BGZF=$(realpath "$INPUT_BGZF")
INPUT_PLAIN=$(realpath "$INPUT_PLAIN")
INPUT_BCF=$(realpath "$INPUT_BCF")
mkdir -p "$RESULT_ROOT"/{golden,logs,runs,scratch}
RESULT_ROOT=$(realpath "$RESULT_ROOT")

case_parameters() {
    case "$1" in
        bgzf_tbi)
            kind=gzvcf
            source_input=$INPUT_BGZF
            baseline_kind=gzvcf
            unique_input=no
            extra=()
            ;;
        bgzf_auto_csi)
            kind=gzvcf
            source_input=$INPUT_BGZF
            baseline_kind=gzvcf
            unique_input=yes
            extra=(--bcftools "$BCFTOOLS")
            ;;
        bgzf_no_auto_index)
            kind=gzvcf
            source_input=$INPUT_BGZF
            baseline_kind=gzvcf
            unique_input=yes
            extra=(--no-auto-index --bcftools "$BCFTOOLS")
            ;;
        plain_vcf)
            kind=vcf
            source_input=$INPUT_PLAIN
            baseline_kind=vcf
            unique_input=no
            extra=()
            ;;
        bcf_csi)
            kind=bcf
            source_input=$INPUT_BCF
            baseline_kind=bcf
            unique_input=no
            extra=()
            ;;
        bcf_auto_csi)
            kind=bcf
            source_input=$INPUT_BCF
            baseline_kind=bcf
            unique_input=yes
            extra=(--bcftools "$BCFTOOLS")
            ;;
        bcf_no_auto_index)
            kind=bcf
            source_input=$INPUT_BCF
            baseline_kind=bcf
            unique_input=yes
            extra=(--no-auto-index --bcftools "$BCFTOOLS")
            ;;
        *)
            printf 'Unknown case: %s\n' "$1" >&2
            exit 2
            ;;
    esac
}

write_run_metadata() {
    local metadata=$1
    local case_name=$2
    local engine=$3
    local threads=$4
    local timing=$5
    local output=$6
    local exact=$7
    local hash
    hash=$(sha256sum "$output" | cut -d' ' -f1)
    local wall user system cpu rss
    IFS=$'\t' read -r wall user system cpu rss <"$timing"
    cpu=${cpu%\%}
    printf 'case\tengine\tthreads\twall_s\tuser_s\tsystem_s\tcpu_pct\tmax_rss_kb\toutput_bytes\toutput_sha256\texact\n' \
        >"$metadata"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$case_name" "$engine" "$threads" \
        "$wall" "$user" "$system" "$cpu" "$rss" \
        "$(stat -c '%s' "$output")" "$hash" "$exact" \
        >>"$metadata"
}

run_original() {
    local baseline_kind=$1
    local input=$2
    local golden="$RESULT_ROOT/golden/$baseline_kind.vcf"
    local metadata="$RESULT_ROOT/runs/original-$baseline_kind.tsv"
    [[ -s "$golden" && -s "$metadata" ]] && return
    local timing="$RESULT_ROOT/runs/original-$baseline_kind.time.txt"
    printf 'START original-%s %s\n' \
        "$baseline_kind" "$(date --iso-8601=seconds)"
    /usr/bin/time -f '%e\t%U\t%S\t%P\t%M' -o "$timing" \
        "$ORIGINAL" "--$baseline_kind" "$input" \
        "${filter_args[@]}" \
        --recode --recode-INFO-all --stdout \
        >"$golden" \
        2>"$RESULT_ROOT/logs/original-$baseline_kind.stderr.txt"
    write_run_metadata \
        "$metadata" "$baseline_kind" original 1 \
        "$timing" "$golden" ORACLE
    printf 'DONE original-%s\n' "$baseline_kind"
}

run_candidate() {
    local case_name=$1
    local threads=$2
    case_parameters "$case_name"
    local run_id="$case_name-t$threads"
    local metadata="$RESULT_ROOT/runs/$run_id.tsv"
    [[ -s "$metadata" ]] && return
    local scratch="$RESULT_ROOT/scratch/$run_id"
    rm -rf -- "$scratch"
    mkdir -p "$scratch"
    local input=$source_input
    if [[ "$unique_input" == yes ]]; then
        if [[ "$kind" == gzvcf ]]; then
            input="$scratch/input.vcf.gz"
        else
            input="$scratch/input.bcf"
        fi
        ln -s "$source_input" "$input"
    fi
    local output="$scratch/output.vcf"
    local timing="$RESULT_ROOT/runs/$run_id.time.txt"
    printf 'START %s %s\n' "$run_id" "$(date --iso-8601=seconds)"
    /usr/bin/time -f '%e\t%U\t%S\t%P\t%M' -o "$timing" \
        "$NG" "--$kind" "$input" --threads "$threads" \
        "${extra[@]}" "${filter_args[@]}" \
        --recode --recode-INFO-all --stdout \
        >"$output" \
        2>"$RESULT_ROOT/logs/$run_id.stderr.txt"
    local golden="$RESULT_ROOT/golden/$baseline_kind.vcf"
    cmp "$golden" "$output"
    write_run_metadata \
        "$metadata" "$case_name" vcftools-ng "$threads" \
        "$timing" "$output" PASS
    rm -rf -- "$scratch"
    printf 'DONE %s exact=PASS\n' "$run_id"
}

run_original gzvcf "$INPUT_BGZF"
run_original vcf "$INPUT_PLAIN"
run_original bcf "$INPUT_BCF"

for case_name in "${cases[@]}"; do
    for threads in $THREAD_LIST; do
        [[ "$threads" =~ ^[1-9][0-9]*$ ]] || {
            printf 'Invalid thread count: %s\n' "$threads" >&2
            exit 2
        }
        run_candidate "$case_name" "$threads"
    done
done

all_runs="$RESULT_ROOT/all-runs.tsv"
head -1 "$RESULT_ROOT/runs/original-gzvcf.tsv" >"$all_runs"
find "$RESULT_ROOT/runs" -maxdepth 1 -name '*.tsv' -type f -print0 |
    sort -z | xargs -0 -r -n1 tail -1 >>"$all_runs"

awk -F '\t' '
    BEGIN {
        OFS = "\t"
        print "case", "threads", "original_wall_s", "ng_wall_s", \
              "speedup", "cpu_pct", "max_rss_kb", "exact"
    }
    FNR == NR {
        if (FNR > 1 && $2 == "original") {
            baseline[$1] = $4
        }
        next
    }
    FNR == 1 { next }
    $2 == "original" {
        baseline[$1] = $4
        next
    }
    $2 == "vcftools-ng" {
        kind = $1 ~ /^bgzf_/ ? "gzvcf" :
               ($1 ~ /^bcf_/ ? "bcf" : "vcf")
        printf "%s\t%s\t%.2f\t%.2f\t%.4f\t%.1f\t%d\t%s\n", \
               $1, $3, baseline[kind], $4, baseline[kind] / $4, \
               $7, $8, $11
    }
' "$all_runs" "$all_runs" | {
    IFS= read -r header
    printf '%s\n' "$header"
    sort -t $'\t' -k1,1 -k2,2n
} >"$RESULT_ROOT/summary.tsv"

{
    printf 'key\tvalue\n'
    printf 'date\t%s\n' "$(date --iso-8601=seconds)"
    printf 'host\t%s\n' "$(hostname)"
    printf 'os\t%s\n' "$(
        . /etc/os-release
        printf '%s' "$PRETTY_NAME"
    )"
    printf 'kernel\t%s\n' "$(uname -srmo)"
    printf 'cpu_model\t%s\n' "$(
        lscpu | awk -F: '/Model name/ {
            sub(/^[[:space:]]+/, "", $2); print $2; exit
        }'
    )"
    printf 'logical_cpus\t%s\n' "$(nproc)"
    printf 'memory\t%s\n' "$(free -h | awk '/^Mem:/ {print $2}')"
    printf 'gpu\t%s\n' "$(
        if command -v nvidia-smi >/dev/null 2>&1; then
            nvidia-smi --query-gpu=name,memory.total \
                --format=csv,noheader | paste -sd ';' -
        else
            lspci 2>/dev/null |
                awk '/VGA compatible controller|3D controller/ {
                    sub(/^.*: /, ""); print
                }' | paste -sd ';' -
        fi
    )"
    printf 'block_devices\t%s\n' "$(
        lsblk -dn -o NAME,SIZE,ROTA,TYPE,MODEL | paste -sd ';' -
    )"
    printf 'ng_version\t%s\n' "$("$NG" --version)"
    printf 'git_commit\t%s\n' "$(git -C "$source_root" rev-parse HEAD)"
    printf 'original\t%s\n' "$ORIGINAL"
    printf 'bcftools\t%s\n' "$("$BCFTOOLS" --version | head -1)"
    printf 'threads\t%s\n' "$THREAD_LIST"
    printf 'expected_records\t%s\n' "$EXPECTED_RECORDS"
    printf 'filter_args\t%s\n' "${filter_args[*]}"
    printf 'output\t--recode --recode-INFO-all --stdout\n'
    printf 'cache\toperating_system_cache_not_flushed\n'
    printf 'input_bgzf\t%s\n' "$INPUT_BGZF"
    printf 'input_bgzf_sha256\t%s\n' \
        "$(sha256sum "$INPUT_BGZF" | cut -d' ' -f1)"
    printf 'input_plain\t%s\n' "$INPUT_PLAIN"
    printf 'input_plain_sha256\t%s\n' \
        "$(sha256sum "$INPUT_PLAIN" | cut -d' ' -f1)"
    printf 'input_bcf\t%s\n' "$INPUT_BCF"
    printf 'input_bcf_sha256\t%s\n' \
        "$(sha256sum "$INPUT_BCF" | cut -d' ' -f1)"
} >"$RESULT_ROOT/environment.tsv"

printf 'v0.11.4 real-filter matrix complete: %s\n' "$RESULT_ROOT"
