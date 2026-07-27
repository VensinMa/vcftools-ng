#!/usr/bin/env bash
set -euo pipefail

source_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

NG=${NG:-"$source_root/build/vcftools-ng"}
BCFTOOLS=${BCFTOOLS:-"/home/vensin/software/bcftools/bcftools"}
FULL_BGZF=${FULL_BGZF:-"$source_root/data/osmanthus412.snps.vcf.gz"}
FULL_PLAIN=${FULL_PLAIN:-"$source_root/data/osmanthus412.snps.vcf"}
FULL_BCF=${FULL_BCF:-"$source_root/data/osmanthus412.snps.bcf"}
ORACLE_ROOT=${ORACLE_ROOT:-"$source_root/benchmarks/results/final-full-v0121/golden"}
ORACLE_RUN_ROOT=${ORACLE_RUN_ROOT:-"$source_root/benchmarks/results/final-full-v0121/runs"}
THREAD_LIST=${THREAD_LIST:-"1 2 4 8 16 32"}
REPEATS=${REPEATS:-5}
EXPECTED_RECORDS=${EXPECTED_RECORDS:-11230392}
GATE_ONLY=${GATE_ONLY:-0}
RESULT_ROOT=${1:-"$source_root/benchmarks/results/final-full-v0122"}

readonly required_threads="1 2 4 8 16 32"
readonly filter_args=(
    --min-alleles 2 --max-alleles 2
    --minGQ 10 --minQ 30
    --min-meanDP 7 --max-missing 0.9 --maf 0.1
)
readonly cases=(
    bgzf_tbi
    bgzf_auto_csi
    plain_vcf
    bcf_adaptive
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
[[ "$GATE_ONLY" == 0 || "$GATE_ONLY" == 1 ]] || {
    printf 'GATE_ONLY must be 0 or 1\n' >&2
    exit 2
}
for executable in "$NG" "$BCFTOOLS"; do
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
[[ -f "$FULL_BGZF.tbi" ]] || {
    printf 'BGZF indexed scenario requires TBI: %s.tbi\n' \
        "$FULL_BGZF" >&2
    exit 2
}
for format in "${formats[@]}"; do
    [[ -s "$ORACLE_ROOT/$format.vcf" ]] || {
        printf 'Locked Original oracle is missing: %s/%s.vcf\n' \
            "$ORACLE_ROOT" "$format" >&2
        exit 2
    }
    [[ -s "$ORACLE_RUN_ROOT/original-$format-r1.tsv" ]] || {
        printf 'Locked Original timing is missing: %s/original-%s-r1.tsv\n' \
            "$ORACLE_RUN_ROOT" "$format" >&2
        exit 2
    }
done

NG=$(realpath "$NG")
BCFTOOLS=$(realpath "$BCFTOOLS")
FULL_BGZF=$(realpath "$FULL_BGZF")
FULL_PLAIN=$(realpath "$FULL_PLAIN")
FULL_BCF=$(realpath "$FULL_BCF")
ORACLE_ROOT=$(realpath "$ORACLE_ROOT")
ORACLE_RUN_ROOT=$(realpath "$ORACLE_RUN_ROOT")
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
    extra=()
    unique_input=no
    case "$1" in
        bgzf_tbi)
            kind=gzvcf
            source_input=$FULL_BGZF
            ;;
        bgzf_auto_csi)
            kind=gzvcf
            source_input=$FULL_BGZF
            unique_input=yes
            extra=(--bcftools "$BCFTOOLS")
            ;;
        plain_vcf)
            kind=vcf
            source_input=$FULL_PLAIN
            ;;
        bcf_adaptive)
            kind=bcf
            source_input=$FULL_BCF
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
    local backend=${11}
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
    printf 'scenario\tformat\tengine\tthreads\trepeat\twall_s\tuser_s\tsystem_s\tcpu_pct\tmax_rss_kb\toutput_bytes\texit\texact\toutput_sha256\tbackend\n' \
        >"$metadata"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$scenario" "$format" "$engine" "$threads" "$repeat" \
        "$wall" "$user" "$system" "$cpu" "$rss" "$bytes" \
        "$exit_status" "$exact" "$hash" "$backend" >>"$metadata"
}

link_locked_golden() {
    local format=$1
    local target="$RESULT_ROOT/golden/$format.vcf"
    if [[ ! -e "$target" ]]; then
        ln "$ORACLE_ROOT/$format.vcf" "$target" 2>/dev/null ||
            ln -s "$ORACLE_ROOT/$format.vcf" "$target"
    fi
    cmp -s "$ORACLE_ROOT/$format.vcf" "$target" || {
        printf 'Linked golden changed: %s\n' "$target" >&2
        exit 1
    }
}

lock_original_baseline() {
    local format=$1
    link_locked_golden "$format"
    local run_id="original-locked-$format"
    local metadata="$RESULT_ROOT/runs/$run_id.tsv"
    if completed_run "$metadata"; then
        printf 'SKIP completed %s\n' "$run_id"
        return
    fi
    local locked="$ORACLE_RUN_ROOT/original-$format-r1.tsv"
    if ! awk -F '\t' \
        'NR == 2 && $3 == "original" && $12 == 0 &&
         $13 == "ORACLE" {ok=1} END {exit !ok}' "$locked"; then
        printf 'Locked Original metadata failed validation: %s\n' \
            "$locked" >&2
        exit 1
    fi
    IFS=$'\t' read -r _ _ _ _ _ wall user system cpu rss \
        _ _ _ _ < <(sed -n '2p' "$locked")
    local output="$RESULT_ROOT/golden/$format.vcf"
    local timing="$RESULT_ROOT/runs/$run_id.time.txt"
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$wall" "$user" "$system" "${cpu}%" "$rss" >"$timing"
    write_metadata \
        "$metadata" "$format" "$format" original 1 1 \
        "$timing" "$output" 0 ORACLE locked-v0.12.1
    printf 'LOCK PASS %s original_wall=%s oracle=%s\n' \
        "$format" "$wall" "$(sha256sum "$output" | cut -d' ' -f1)"
}

assert_adaptive_backend() {
    local scenario=$1
    local threads=$2
    local input=$3
    local backend=$4
    case "$scenario" in
        bgzf_tbi)
            if ((threads == 1)); then
                [[ "$backend" == stream ]]
            else
                [[ "$backend" == indexed-regions ]]
            fi
            ;;
        bgzf_auto_csi)
            if ((threads == 1)); then
                [[ "$backend" == stream && ! -e "$input.csi" ]]
            else
                [[ "$backend" == indexed-regions && -s "$input.csi" ]]
            fi
            ;;
        plain_vcf)
            if ((threads <= 2)); then
                [[ "$backend" == stream ]]
            else
                [[ "$backend" == plain-ranges ]]
            fi
            ;;
        bcf_adaptive)
            [[ "$backend" == stream ]]
            ;;
    esac
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
        input="$scratch/input.vcf.gz"
        ln -s "$source_input" "$input"
    fi
    local output="$scratch/output.vcf"
    local timing="$RESULT_ROOT/runs/$run_id.time.txt"
    local log="$RESULT_ROOT/logs/$run_id.stderr.txt"
    local status=0 exact=FAIL backend=missing
    printf 'START %s %s\n' "$run_id" "$(date --iso-8601=seconds)"
    /usr/bin/time -f '%e\t%U\t%S\t%P\t%M' -o "$timing" \
        "$NG" "--$kind" "$input" --threads "$threads" \
        "${extra[@]}" "${filter_args[@]}" \
        --log-file "$log" \
        --recode --recode-INFO-all --stdout \
        >"$output" 2>"$log" ||
        status=$?
    backend=$(
        sed -n 's/^Input backend: \([^ ]*\).*/\1/p' "$log" |
            head -1
    )
    if ((status == 0)) &&
       cmp -s "$RESULT_ROOT/golden/$kind.vcf" "$output" &&
       assert_adaptive_backend \
           "$scenario" "$threads" "$input" "$backend"; then
        exact=PASS
    fi
    write_metadata \
        "$metadata" "$scenario" "$kind" vcftools-ng \
        "$threads" "$repeat" "$timing" "$output" \
        "$status" "$exact" "$backend"
    rm -rf -- "$scratch"
    if ((status != 0)) || [[ "$exact" != PASS ]]; then
        printf 'FAILED %s exit=%s exact=%s backend=%s\n' \
            "$run_id" "$status" "$exact" "$backend" >&2
        exit 1
    fi
    printf 'DONE %s exact=PASS backend=%s %s\n' \
        "$run_id" "$backend" "$(date --iso-8601=seconds)"
}

refresh_reports() {
    local all_runs="$RESULT_ROOT/all-runs.tsv"
    head -1 "$RESULT_ROOT/runs/original-locked-gzvcf.tsv" >"$all_runs"
    find "$RESULT_ROOT/runs" -maxdepth 1 -type f -name '*.tsv' -print0 |
        sort -z | xargs -0 -r -n1 tail -1 >>"$all_runs"

    awk -F '\t' '
        BEGIN {
            OFS = "\t"
            print "scenario", "threads", "runs", \
                  "original_single_wall_s", "ng_mean_wall_s", \
                  "speedup", "mean_cpu_pct", "max_rss_kb", \
                  "backend", "exact"
        }
        NR == 1 { next }
        $3 == "original" {
            original_wall[$2] = $6
            next
        }
        $3 == "vcftools-ng" {
            key = $1 SUBSEP $4
            format[key] = $2
            n[key]++
            wall[key] += $6
            cpu[key] += $9
            if (!(key in rss) || $10 > rss[key]) rss[key] = $10
            if (!(key in selected_backend)) {
                selected_backend[key] = $15
            } else if (selected_backend[key] != $15) {
                selected_backend[key] = selected_backend[key] "," $15
            }
            if ($13 != "PASS") exact[key] = "FAIL"
            else if (!(key in exact)) exact[key] = "PASS"
        }
        END {
            for (key in n) {
                split(key, fields, SUBSEP)
                baseline = original_wall[format[key]]
                mean = wall[key] / n[key]
                printf "%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.1f\t%d\t%s\t%s\n", \
                       fields[1], fields[2], n[key], baseline, mean, \
                       baseline / mean, cpu[key] / n[key], rss[key], \
                       selected_backend[key], exact[key]
            }
        }
    ' "$all_runs" | {
        IFS= read -r header
        printf '%s\n' "$header"
        sort -t $'\t' -k1,1 -k2,2n
    } >"$RESULT_ROOT/summary.tsv"
}

manifest="$RESULT_ROOT/manifest.tsv"
if [[ ! -s "$manifest" ]]; then
    {
        printf 'key\tvalue\n'
        printf 'started\t%s\n' "$(date --iso-8601=seconds)"
        printf 'release\tv0.12.2\n'
        printf 'title\tAdaptive Indexing for Fast Full Scans\n'
        printf 'git_commit\t%s\n' "$(git -C "$source_root" rev-parse HEAD)"
        printf 'working_tree_sha256\t%s\n' "$(
            git -C "$source_root" diff --binary |
                sha256sum | cut -d' ' -f1
        )"
        printf 'ng_version\t%s\n' "$("$NG" --version)"
        printf 'original\tVCFtools (0.1.17), locked v0.12.1 baseline\n'
        printf 'original_runs_v0122\t0\n'
        printf 'original_policy\thash-validated retained oracle and single-run timing; not rerun\n'
        printf 'vcftools_ng_repeats\t%s\n' "$REPEATS"
        printf 'expected_records\t%s\n' "$EXPECTED_RECORDS"
        printf 'threads\t%s\n' "$THREAD_LIST"
        printf 'scenarios\tBGZF VCF + TBI;BGZF VCF + automatic CSI;Plain VCF;BCF adaptive\n'
        printf 'execution\tstrictly_serial\n'
        printf 'cache\toperating_system_cache_not_flushed\n'
        printf 'workload\t%s --recode --recode-INFO-all --stdout\n' \
            "${filter_args[*]}"
        for item in \
            "bgzf:$FULL_BGZF" \
            "bgzf_index:$FULL_BGZF.tbi" \
            "plain:$FULL_PLAIN" \
            "bcf:$FULL_BCF"
        do
            key=${item%%:*}
            path=${item#*:}
            printf '%s\t%s\n' "$key" "$path"
            printf '%s_bytes\t%s\n' "$key" "$(stat -c '%s' "$path")"
            printf '%s_sha256\t%s\n' \
                "$key" "$(sha256sum "$path" | cut -d' ' -f1)"
        done
        printf 'locked_oracle_root\t%s\n' "$ORACLE_ROOT"
        printf 'locked_oracle_run_root\t%s\n' "$ORACLE_RUN_ROOT"
        for format in "${formats[@]}"; do
            printf 'locked_oracle_%s_bytes\t%s\n' \
                "$format" "$(stat -c '%s' "$ORACLE_ROOT/$format.vcf")"
            printf 'locked_oracle_%s_sha256\t%s\n' \
                "$format" "$(
                    sha256sum "$ORACLE_ROOT/$format.vcf" |
                        cut -d' ' -f1
                )"
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

# Reuse the retained Original run because input bytes, Original version, and
# workload are unchanged. Any missing or invalid artifact stops the matrix.
for format in "${formats[@]}"; do
    lock_original_baseline "$format"
done

# All first repeats must pass exactness and backend-policy gates before timing
# repeats are allowed to begin.
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

refresh_reports
if ! grep -q '^first_repeat_gates_passed	' "$manifest"; then
    printf 'first_repeat_gates_passed\t%s\n' \
        "$(date --iso-8601=seconds)" >>"$manifest"
fi
if [[ "$GATE_ONLY" == 1 ]]; then
    printf 'FIRST-REPEAT v0.12.2 RELEASE GATES PASS %s\n' "$RESULT_ROOT"
    printf 'Repeats 2-%s remain pending; resume with GATE_ONLY=0.\n' \
        "$REPEATS"
    exit 0
fi

for repeat in $(seq 2 "$REPEATS"); do
    for scenario in "${cases[@]}"; do
        for threads in $THREAD_LIST; do
            run_candidate "$scenario" "$threads" "$repeat"
        done
    done
done

refresh_reports
if ! grep -q '^finished	' "$manifest"; then
    printf 'finished\t%s\n' "$(date --iso-8601=seconds)" >>"$manifest"
fi
printf 'FINAL v0.12.2 FOUR-SCENARIO RELEASE MATRIX PASS %s\n' "$RESULT_ROOT"
