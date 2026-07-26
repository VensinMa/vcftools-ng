#!/usr/bin/env bash
set -euo pipefail

source_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

NG=${NG:-"$source_root/build/vcftools-ng"}
ORIGINAL=${ORIGINAL:-"/home/vensin/anaconda3/envs/vcftools/bin/vcftools"}
BCFTOOLS=${BCFTOOLS:-"/home/vensin/software/bcftools/bcftools"}
FULL_BGZF=${FULL_BGZF:-"$source_root/data/osmanthus412.snps.vcf.gz"}
FULL_PLAIN=${FULL_PLAIN:-"$source_root/data/osmanthus412.snps.vcf"}
FULL_BCF=${FULL_BCF:-"$source_root/data/osmanthus412.snps.bcf"}
THREAD_LIST=${THREAD_LIST:-"1 2 4 8 16 32"}
REPEATS=${REPEATS:-5}
EXPECTED_RECORDS=${EXPECTED_RECORDS:-11230392}
RESULT_ROOT=${1:-"$source_root/benchmarks/results/final-full-v0121"}

readonly required_threads="1 2 4 8 16 32"
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
readonly formats=(gzvcf vcf bcf)

[[ "$THREAD_LIST" == "$required_threads" ]] || {
    printf 'Release matrix requires THREAD_LIST="%s"\n' \
        "$required_threads" >&2
    exit 2
}
[[ "$REPEATS" == 5 ]] || {
    printf 'Release matrix requires REPEATS=5\n' >&2
    exit 2
}
for executable in "$NG" "$ORIGINAL" "$BCFTOOLS"; do
    [[ -x "$executable" ]] || {
        printf 'Missing executable: %s\n' "$executable" >&2
        exit 2
    }
done
for input in "$FULL_BGZF" "$FULL_PLAIN" "$FULL_BCF"; do
    [[ -f "$input" ]] || {
        printf 'Missing full input: %s\n' "$input" >&2
        exit 2
    }
done
[[ -f "$FULL_BGZF.tbi" || -f "$FULL_BGZF.csi" ]] || {
    printf 'BGZF indexed scenario requires TBI or CSI\n' >&2
    exit 2
}
[[ -f "$FULL_BCF.csi" ]] || {
    printf 'BCF indexed scenario requires CSI\n' >&2
    exit 2
}

NG=$(realpath "$NG")
ORIGINAL=$(realpath "$ORIGINAL")
BCFTOOLS=$(realpath "$BCFTOOLS")
FULL_BGZF=$(realpath "$FULL_BGZF")
FULL_PLAIN=$(realpath "$FULL_PLAIN")
FULL_BCF=$(realpath "$FULL_BCF")
mkdir -p "$RESULT_ROOT"/{gates,golden,logs,runs,scratch}
RESULT_ROOT=$(realpath "$RESULT_ROOT")

format_parameters() {
    case "$1" in
        gzvcf)
            source_input=$FULL_BGZF
            ;;
        vcf)
            source_input=$FULL_PLAIN
            ;;
        bcf)
            source_input=$FULL_BCF
            ;;
        *)
            printf 'Unknown format: %s\n' "$1" >&2
            exit 2
            ;;
    esac
}

case_parameters() {
    case "$1" in
        bgzf_tbi)
            kind=gzvcf
            source_input=$FULL_BGZF
            unique_input=no
            extra=()
            ;;
        bgzf_auto_csi)
            kind=gzvcf
            source_input=$FULL_BGZF
            unique_input=yes
            extra=(--bcftools "$BCFTOOLS")
            ;;
        bgzf_no_auto_index)
            kind=gzvcf
            source_input=$FULL_BGZF
            unique_input=yes
            extra=(--no-auto-index --bcftools "$BCFTOOLS")
            ;;
        plain_vcf)
            kind=vcf
            source_input=$FULL_PLAIN
            unique_input=no
            extra=()
            ;;
        bcf_csi)
            kind=bcf
            source_input=$FULL_BCF
            unique_input=no
            extra=()
            ;;
        bcf_auto_csi)
            kind=bcf
            source_input=$FULL_BCF
            unique_input=yes
            extra=(--bcftools "$BCFTOOLS")
            ;;
        bcf_no_auto_index)
            kind=bcf
            source_input=$FULL_BCF
            unique_input=yes
            extra=(--no-auto-index --bcftools "$BCFTOOLS")
            ;;
        *)
            printf 'Unknown scenario: %s\n' "$1" >&2
            exit 2
            ;;
    esac
}

completed_run() {
    local metadata=$1
    [[ -s "$metadata" ]] &&
        awk -F '\t' \
            'NR == 2 && $12 == 0 && $13 ~ /^(PASS|ORACLE)$/ {
                 ok=1
             }
             END {exit !ok}' \
            "$metadata"
}

write_metadata() {
    local metadata=$1
    local scenario=$2
    local format=$3
    local engine=$4
    local threads=$5
    local repeat=$6
    local timing=$7
    local output=$8
    local exit_status=$9
    local exact=${10}
    local wall=NA user=NA system=NA cpu=NA rss=NA
    if [[ -s "$timing" ]]; then
        IFS=$'\t' read -r wall user system cpu rss <"$timing"
        cpu=${cpu%\%}
    fi
    local bytes=0 hash=missing
    if [[ -s "$output" ]]; then
        bytes=$(stat -c '%s' "$output")
        hash=$(sha256sum "$output" | cut -d' ' -f1)
    fi
    printf 'scenario\tformat\tengine\tthreads\trepeat\twall_s\tuser_s\tsystem_s\tcpu_pct\tmax_rss_kb\toutput_bytes\texit\texact\toutput_sha256\n' \
        >"$metadata"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$scenario" "$format" "$engine" "$threads" "$repeat" \
        "$wall" "$user" "$system" "$cpu" "$rss" "$bytes" \
        "$exit_status" "$exact" "$hash" >>"$metadata"
}

run_original() {
    local format=$1
    local repeat=$2
    format_parameters "$format"
    local run_id="original-$format-r$repeat"
    local metadata="$RESULT_ROOT/runs/$run_id.tsv"
    local golden="$RESULT_ROOT/golden/$format.vcf"
    if completed_run "$metadata" &&
       { [[ "$repeat" != 1 ]] || [[ -s "$golden" ]]; }; then
        printf 'SKIP completed %s\n' "$run_id"
        return
    fi
    local scratch="$RESULT_ROOT/scratch/$run_id"
    rm -rf -- "$scratch"
    mkdir -p "$scratch"
    local output="$scratch/output.vcf"
    if [[ "$repeat" == 1 ]]; then
        output=$golden
        rm -f -- "$output"
    fi
    local timing="$RESULT_ROOT/runs/$run_id.time.txt"
    local status=0 exact=FAIL
    printf 'START %s %s\n' "$run_id" "$(date --iso-8601=seconds)"
    /usr/bin/time -f '%e\t%U\t%S\t%P\t%M' -o "$timing" \
        "$ORIGINAL" "--$format" "$source_input" \
        "${filter_args[@]}" \
        --recode --recode-INFO-all --stdout \
        >"$output" \
        2>"$RESULT_ROOT/logs/$run_id.stderr.txt" ||
        status=$?
    if ((status == 0)) && [[ -s "$output" ]]; then
        if [[ "$repeat" == 1 ]]; then
            exact=ORACLE
        elif cmp -s "$golden" "$output"; then
            exact=PASS
        fi
    fi
    write_metadata \
        "$metadata" "$format" "$format" original 1 "$repeat" \
        "$timing" "$output" "$status" "$exact"
    if [[ "$repeat" != 1 ]]; then
        rm -rf -- "$scratch"
    fi
    if ((status != 0)) || [[ "$exact" == FAIL ]]; then
        printf 'FAILED %s exit=%s exact=%s\n' \
            "$run_id" "$status" "$exact" >&2
        exit 1
    fi
    printf 'DONE %s exact=%s %s\n' \
        "$run_id" "$exact" "$(date --iso-8601=seconds)"
}

run_candidate() {
    local scenario=$1
    local threads=$2
    local repeat=$3
    case_parameters "$scenario"
    local run_id="$scenario-vcftools-ng-t$threads-r$repeat"
    local metadata="$RESULT_ROOT/runs/$run_id.tsv"
    if completed_run "$metadata"; then
        printf 'SKIP completed %s\n' "$run_id"
        return
    fi
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
    local status=0 exact=FAIL
    printf 'START %s %s\n' "$run_id" "$(date --iso-8601=seconds)"
    /usr/bin/time -f '%e\t%U\t%S\t%P\t%M' -o "$timing" \
        "$NG" "--$kind" "$input" --threads "$threads" \
        "${extra[@]}" "${filter_args[@]}" \
        --recode --recode-INFO-all --stdout \
        >"$output" \
        2>"$RESULT_ROOT/logs/$run_id.stderr.txt" ||
        status=$?
    local golden="$RESULT_ROOT/golden/$kind.vcf"
    if ((status == 0)) && [[ -s "$output" ]] &&
       cmp -s "$golden" "$output"; then
        exact=PASS
    fi
    write_metadata \
        "$metadata" "$scenario" "$kind" vcftools-ng \
        "$threads" "$repeat" "$timing" "$output" "$status" "$exact"
    rm -rf -- "$scratch"
    if ((status != 0)) || [[ "$exact" != PASS ]]; then
        printf 'FAILED %s exit=%s exact=%s\n' \
            "$run_id" "$status" "$exact" >&2
        exit 1
    fi
    printf 'DONE %s exact=PASS %s\n' \
        "$run_id" "$(date --iso-8601=seconds)"
}

manifest="$RESULT_ROOT/manifest.tsv"
if [[ ! -s "$manifest" ]]; then
    bgzf_index=$FULL_BGZF.tbi
    [[ -f "$bgzf_index" ]] || bgzf_index=$FULL_BGZF.csi
    {
        printf 'key\tvalue\n'
        printf 'started\t%s\n' "$(date --iso-8601=seconds)"
        printf 'release\tv0.12.1\n'
        printf 'title\tFused Site Statistics and Scalable Exact Recode\n'
        printf 'git_commit\t%s\n' "$(git -C "$source_root" rev-parse HEAD)"
        printf 'ng_version\t%s\n' "$("$NG" --version)"
        printf 'original\t%s\n' "$("$ORIGINAL" --version 2>&1 | head -1)"
        printf 'bcftools\t%s\n' "$("$BCFTOOLS" --version | head -1)"
        printf 'expected_records\t%s\n' "$EXPECTED_RECORDS"
        printf 'threads\t%s\n' "$THREAD_LIST"
        printf 'repeats\t%s\n' "$REPEATS"
        printf 'execution\tstrictly_serial\n'
        printf 'cache\toperating_system_cache_not_flushed\n'
        printf 'workload\t%s --recode --recode-INFO-all --stdout\n' \
            "${filter_args[*]}"
        for item in \
            "bgzf:$FULL_BGZF" \
            "bgzf_index:$bgzf_index" \
            "plain:$FULL_PLAIN" \
            "bcf:$FULL_BCF" \
            "bcf_index:$FULL_BCF.csi"
        do
            key=${item%%:*}
            path=${item#*:}
            printf '%s\t%s\n' "$key" "$path"
            printf '%s_bytes\t%s\n' "$key" "$(stat -c '%s' "$path")"
            printf '%s_sha256\t%s\n' \
                "$key" "$(sha256sum "$path" | cut -d' ' -f1)"
        done
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
    } >"$manifest"
fi

# Phase 1: establish one golden per actual format, then require every
# scenario/thread first repeat to pass before any timing repeats begin.
for format in "${formats[@]}"; do
    run_original "$format" 1
done
for scenario in "${cases[@]}"; do
    gate="$RESULT_ROOT/gates/$scenario.pass"
    if [[ ! -s "$gate" ]]; then
        for threads in $THREAD_LIST; do
            run_candidate "$scenario" "$threads" 1
        done
        printf 'scenario\tstatus\tvalidated\n%s\tPASS\t%s\n' \
            "$scenario" "$(date --iso-8601=seconds)" >"$gate"
        printf 'GATE PASS %s\n' "$scenario"
    fi
done

printf 'ALL FIRST-REPEAT GATES PASSED; STARTING REPEATS 2-%s\n' \
    "$REPEATS"
for repeat in $(seq 2 "$REPEATS"); do
    for format in "${formats[@]}"; do
        run_original "$format" "$repeat"
    done
done
for scenario in "${cases[@]}"; do
    for threads in $THREAD_LIST; do
        for repeat in $(seq 2 "$REPEATS"); do
            run_candidate "$scenario" "$threads" "$repeat"
        done
    done
done

all_runs="$RESULT_ROOT/all-runs.tsv"
head -1 "$RESULT_ROOT/runs/original-gzvcf-r1.tsv" >"$all_runs"
find "$RESULT_ROOT/runs" -maxdepth 1 -type f -name '*.tsv' -print0 |
    sort -z | xargs -0 -r -n1 tail -1 >>"$all_runs"

summary="$RESULT_ROOT/summary.tsv"
awk -F '\t' '
    BEGIN {
        OFS = "\t"
        print "scenario", "threads", "runs", "original_mean_wall_s", \
              "ng_mean_wall_s", "speedup", "mean_cpu_pct", \
              "max_rss_kb", "exact"
    }
    NR == 1 { next }
    $3 == "original" {
        original_wall[$2] += $6
        original_n[$2]++
        next
    }
    $3 == "vcftools-ng" {
        key = $1 SUBSEP $4
        format[key] = $2
        n[key]++
        wall[key] += $6
        cpu[key] += $9
        if (!(key in rss) || $10 > rss[key]) rss[key] = $10
        if ($13 != "PASS") exact[key] = "FAIL"
        else if (!(key in exact)) exact[key] = "PASS"
    }
    END {
        for (key in n) {
            split(key, fields, SUBSEP)
            baseline = original_wall[format[key]] /
                       original_n[format[key]]
            mean = wall[key] / n[key]
            printf "%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.1f\t%d\t%s\n", \
                   fields[1], fields[2], n[key], baseline, mean, \
                   baseline / mean, cpu[key] / n[key], rss[key], \
                   exact[key]
        }
    }
' "$all_runs" | {
    IFS= read -r header
    printf '%s\n' "$header"
    sort -t $'\t' -k1,1 -k2,2n
} >"$summary"

{
    printf 'finished\t%s\n' "$(date --iso-8601=seconds)"
    for format in "${formats[@]}"; do
        printf 'golden_%s_bytes\t%s\n' \
            "$format" "$(stat -c '%s' "$RESULT_ROOT/golden/$format.vcf")"
        printf 'golden_%s_sha256\t%s\n' \
            "$format" "$(
                sha256sum "$RESULT_ROOT/golden/$format.vcf" |
                    cut -d' ' -f1
            )"
    done
} >>"$manifest"

printf 'FINAL v0.12.1 RELEASE MATRIX PASS %s\n' "$RESULT_ROOT"
