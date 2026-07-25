#!/usr/bin/env bash
set -euo pipefail

source_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

NG=${NG:-"$source_root/build/vcftools-ng"}
ORIGINAL=${ORIGINAL:-"/home/vensin/anaconda3/pkgs/vcftools-0.1.17-pl5321h077b44d_0/bin/vcftools"}
BCFTOOLS=${BCFTOOLS:-"/home/vensin/software/bcftools/bcftools"}
FULL_BGZF=${FULL_BGZF:-"/home/vensin/workspace/Sweet_Osmanthus/05.variant_filter/02.vcftools_filter_snp_indel/412samples.SNP.biallelic.minGQ10.minQ30.meanDP6.maxmiss0.8.maf0.05.vcf.gz"}
FULL_PLAIN=${FULL_PLAIN:-"$source_root/data/osmanthus412.snps.vcf"}
FULL_BCF=${FULL_BCF:-"$source_root/data/osmanthus412.snps.bcf"}
RESULT_ROOT=${1:-"$source_root/benchmarks/results/final-full-v0112"}
THREAD_LIST=${THREAD_LIST:-"1 2 4 8 16 32"}
REPEATS=${REPEATS:-5}
EXPECTED_RECORDS=${EXPECTED_RECORDS:-11230392}
STRICT_FINAL_MATRIX=${STRICT_FINAL_MATRIX:-1}

case_names=(
    bgzf_tbi
    bgzf_auto_csi
    bgzf_no_auto_index
    plain_vcf
    bcf_csi
    bcf_auto_csi
    bcf_no_auto_index
)

for executable in "$NG" "$ORIGINAL" "$BCFTOOLS"; do
    if [[ ! -x "$executable" ]]; then
        printf 'Missing executable: %s\n' "$executable" >&2
        exit 2
    fi
done
for input in "$FULL_BGZF" "$FULL_PLAIN" "$FULL_BCF"; do
    if [[ ! -f "$input" ]]; then
        printf 'Missing full input: %s\n' "$input" >&2
        exit 2
    fi
done
if [[ ! -f "$FULL_BGZF.csi" && ! -f "$FULL_BGZF.tbi" ]]; then
    printf 'Indexed BGZF case requires a CSI or TBI sidecar\n' >&2
    exit 2
fi
if [[ ! -f "$FULL_BCF.csi" ]]; then
    printf 'Indexed BCF case requires a CSI sidecar\n' >&2
    exit 2
fi
if [[ "$STRICT_FINAL_MATRIX" == 1 && "$REPEATS" != 5 ]]; then
    printf 'Final benchmark requires REPEATS=5, got %s\n' "$REPEATS" >&2
    exit 2
fi
if [[ "$STRICT_FINAL_MATRIX" == 1 &&
      "$THREAD_LIST" != "1 2 4 8 16 32" ]]; then
    printf 'Final benchmark requires THREAD_LIST=\"1 2 4 8 16 32\"\n' >&2
    exit 2
fi
if [[ ! "$REPEATS" =~ ^[1-9][0-9]*$ ]]; then
    printf 'REPEATS must be a positive integer\n' >&2
    exit 2
fi

NG=$(realpath "$NG")
ORIGINAL=$(realpath "$ORIGINAL")
BCFTOOLS=$(realpath "$BCFTOOLS")
FULL_BGZF=$(realpath "$FULL_BGZF")
FULL_PLAIN=$(realpath "$FULL_PLAIN")
FULL_BCF=$(realpath "$FULL_BCF")
mkdir -p \
    "$RESULT_ROOT/runs" \
    "$RESULT_ROOT/logs" \
    "$RESULT_ROOT/golden" \
    "$RESULT_ROOT/gates" \
    "$RESULT_ROOT/scratch"
RESULT_ROOT=$(realpath "$RESULT_ROOT")

manifest="$RESULT_ROOT/manifest.tsv"
if [[ ! -s "$manifest" ]]; then
    printf 'key\tvalue\n' >"$manifest"
    {
        printf 'started\t%s\n' "$(date --iso-8601=seconds)"
        printf 'host\t%s\n' "$(hostname)"
        printf 'kernel\t%s\n' "$(uname -sr)"
        printf 'cpu_available\t%s\n' "$(nproc)"
        printf 'memory_bytes\t%s\n' "$(awk '/MemTotal/ {print $2 * 1024}' /proc/meminfo)"
        printf 'ng\t%s\n' "$NG"
        printf 'ng_version\t%s\n' "$("$NG" --version)"
        printf 'original\t%s\n' "$ORIGINAL"
        printf 'bcftools\t%s\n' "$BCFTOOLS"
        printf 'bgzf\t%s\n' "$FULL_BGZF"
        printf 'bgzf_bytes\t%s\n' "$(stat -c '%s' "$FULL_BGZF")"
        printf 'bgzf_sha256\t%s\n' "$(sha256sum "$FULL_BGZF" | cut -d' ' -f1)"
        printf 'plain\t%s\n' "$FULL_PLAIN"
        printf 'plain_bytes\t%s\n' "$(stat -c '%s' "$FULL_PLAIN")"
        printf 'bcf\t%s\n' "$FULL_BCF"
        printf 'bcf_bytes\t%s\n' "$(stat -c '%s' "$FULL_BCF")"
        printf 'bcf_sha256\t%s\n' "$(sha256sum "$FULL_BCF" | cut -d' ' -f1)"
        printf 'expected_records\t%s\n' "$EXPECTED_RECORDS"
        printf 'threads\t%s\n' "$THREAD_LIST"
        printf 'repeats\t%s\n' "$REPEATS"
        printf 'execution\tstrictly_serial\n'
        printf 'cache\toperating_system_cache_not_flushed\n'
        if [[ -f "$FULL_BGZF.csi" ]]; then
            printf 'bgzf_csi_sha256_before\t%s\n' \
                "$(sha256sum "$FULL_BGZF.csi" | cut -d' ' -f1)"
        fi
        if [[ -f "$FULL_BGZF.tbi" ]]; then
            printf 'bgzf_tbi_sha256_before\t%s\n' \
                "$(sha256sum "$FULL_BGZF.tbi" | cut -d' ' -f1)"
        fi
        printf 'bcf_csi_sha256_before\t%s\n' \
            "$(sha256sum "$FULL_BCF.csi" | cut -d' ' -f1)"
    } >>"$manifest"
fi

case_parameters() {
    local case_name=$1
    case "$case_name" in
        bgzf_tbi)
            kind=gzvcf
            source_input=$FULL_BGZF
            unique_input=no
            ng_extra=()
            ;;
        bgzf_auto_csi)
            kind=gzvcf
            source_input=$FULL_BGZF
            unique_input=yes
            ng_extra=(--bcftools "$BCFTOOLS")
            ;;
        bgzf_no_auto_index)
            kind=gzvcf
            source_input=$FULL_BGZF
            unique_input=yes
            ng_extra=(--no-auto-index --bcftools "$BCFTOOLS")
            ;;
        plain_vcf)
            kind=vcf
            source_input=$FULL_PLAIN
            unique_input=no
            ng_extra=()
            ;;
        bcf_csi)
            kind=bcf
            source_input=$FULL_BCF
            unique_input=no
            ng_extra=()
            ;;
        bcf_auto_csi)
            kind=bcf
            source_input=$FULL_BCF
            unique_input=yes
            ng_extra=(--bcftools "$BCFTOOLS")
            ;;
        bcf_no_auto_index)
            kind=bcf
            source_input=$FULL_BCF
            unique_input=yes
            ng_extra=(--no-auto-index --bcftools "$BCFTOOLS")
            ;;
        *)
            printf 'Unknown benchmark case: %s\n' "$case_name" >&2
            exit 2
            ;;
    esac
}

completed_run() {
    local metadata=$1
    [[ -s "$metadata" ]] &&
        awk -F '\t' 'NR == 2 && $11 == 0 && $12 ~ /^(PASS|ORACLE)$/ {ok=1} END {exit !ok}' \
            "$metadata"
}

run_one() {
    local case_name=$1
    local engine=$2
    local threads=$3
    local repeat=$4
    case_parameters "$case_name"

    local run_id="${case_name}-${engine}-t${threads}-r${repeat}"
    local metadata="$RESULT_ROOT/runs/$run_id.tsv"
    local golden="$RESULT_ROOT/golden/$case_name.frq.count"
    if completed_run "$metadata"; then
        if [[ "$engine" != original || "$repeat" != 1 || -s "$golden" ]]; then
            printf 'SKIP completed %s\n' "$run_id"
            return
        fi
    fi

    local scratch="$RESULT_ROOT/scratch/$run_id"
    if [[ -d "$scratch" ]]; then
        rm -rf -- "$scratch"
    fi
    mkdir -p "$scratch"

    local input=$source_input
    if [[ "$unique_input" == yes ]]; then
        input="$scratch/input"
        if [[ "$kind" == gzvcf ]]; then
            input+=".vcf.gz"
        else
            input+=".bcf"
        fi
        ln -s "$source_input" "$input"
    fi

    local output_prefix="$scratch/output"
    if [[ "$engine" == original && "$repeat" == 1 ]]; then
        output_prefix="$RESULT_ROOT/golden/$case_name"
        rm -f -- "$golden"
    fi

    local stdout_log="$RESULT_ROOT/logs/$run_id.stdout.txt"
    local stderr_log="$RESULT_ROOT/logs/$run_id.stderr.txt"
    local timing_file="$RESULT_ROOT/runs/$run_id.time.txt"
    local status=0
    printf 'START %s %s\n' "$run_id" "$(date --iso-8601=seconds)"
    if [[ "$engine" == original ]]; then
        /usr/bin/time \
            -f '%e\t%U\t%S\t%P\t%M' \
            -o "$timing_file" \
            "$ORIGINAL" "--$kind" "$input" --counts \
            --out "$output_prefix" \
            >"$stdout_log" 2>"$stderr_log" ||
            status=$?
    else
        /usr/bin/time \
            -f '%e\t%U\t%S\t%P\t%M' \
            -o "$timing_file" \
            "$NG" "--$kind" "$input" --threads "$threads" \
            "${ng_extra[@]}" --counts --out "$output_prefix" \
            >"$stdout_log" 2>"$stderr_log" ||
            status=$?
    fi

    local output="$output_prefix.frq.count"
    local exact=FAIL
    local output_hash=missing
    if ((status == 0)) && [[ -s "$output" ]]; then
        output_hash=$(sha256sum "$output" | cut -d' ' -f1)
        if [[ "$engine" == original && "$repeat" == 1 ]]; then
            exact=ORACLE
        elif cmp -s "$golden" "$output"; then
            exact=PASS
        fi
    fi

    local wall=NA user=NA system=NA cpu=NA rss=NA
    if [[ -s "$timing_file" ]]; then
        IFS=$'\t' read -r wall user system cpu rss <"$timing_file"
        cpu=${cpu%\%}
    fi
    printf 'case\tengine\tthreads\trepeat\twall_s\tuser_s\tsystem_s\tcpu_pct\tmax_rss_kb\toutput_sha256\texit\texact\n' \
        >"$metadata"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$case_name" "$engine" "$threads" "$repeat" \
        "$wall" "$user" "$system" "$cpu" "$rss" \
        "$output_hash" "$status" "$exact" >>"$metadata"

    if [[ "$engine" != original || "$repeat" != 1 ]]; then
        rm -f -- "$output"
    fi
    if [[ -d "$scratch" ]]; then
        rm -rf -- "$scratch"
    fi

    if ((status != 0)) || [[ "$exact" == FAIL ]]; then
        printf 'FAILED %s exit=%s exact=%s\n' \
            "$run_id" "$status" "$exact" >&2
        exit 1
    fi
    printf 'DONE %s wall=%s exact=%s %s\n' \
        "$run_id" "$wall" "$exact" "$(date --iso-8601=seconds)"
}

for case_name in "${case_names[@]}"; do
    gate="$RESULT_ROOT/gates/$case_name.pass"
    if [[ ! -s "$gate" ]]; then
        run_one "$case_name" original 1 1
        for threads in $THREAD_LIST; do
            run_one "$case_name" vcftools-ng "$threads" 1
        done
        printf 'case\tstatus\tvalidated\n%s\tPASS\t%s\n' \
            "$case_name" "$(date --iso-8601=seconds)" >"$gate"
        printf 'GATE PASS %s\n' "$case_name"
    fi

    if ((REPEATS > 1)); then
        for repeat in $(seq 2 "$REPEATS"); do
            run_one "$case_name" original 1 "$repeat"
        done
    fi
    for threads in $THREAD_LIST; do
        if ((REPEATS > 1)); then
            for repeat in $(seq 2 "$REPEATS"); do
                run_one "$case_name" vcftools-ng "$threads" "$repeat"
            done
        fi
    done
done

all_runs="$RESULT_ROOT/all-runs.tsv"
head -1 "$RESULT_ROOT/runs/${case_names[0]}-original-t1-r1.tsv" >"$all_runs"
find "$RESULT_ROOT/runs" -maxdepth 1 -type f -name '*.tsv' -print0 |
    sort -z |
    xargs -0 -r -n1 tail -n 1 >>"$all_runs"

summary="$RESULT_ROOT/summary.tsv"
awk -F '\t' '
    BEGIN {
        OFS = "\t"
        print "case", "engine", "threads", "runs", "mean_wall_s", \
              "min_wall_s", "max_wall_s", "mean_cpu_pct", \
              "max_rss_kb", "speedup_vs_original", "exact"
    }
    NR == 1 { next }
    {
        key = $1 SUBSEP $2 SUBSEP $3
        n[key]++
        wall[key] += $5
        cpu[key] += $8
        if (!(key in min) || $5 < min[key]) min[key] = $5
        if (!(key in max) || $5 > max[key]) max[key] = $5
        if (!(key in rss) || $9 > rss[key]) rss[key] = $9
        exact[key] = exact[key] == "FAIL" || $12 == "FAIL" ? "FAIL" : "PASS"
        if ($2 == "original") original[$1] += $5
        if ($2 == "original") original_n[$1]++
    }
    END {
        for (key in n) {
            split(key, fields, SUBSEP)
            baseline = original[fields[1]] / original_n[fields[1]]
            mean = wall[key] / n[key]
            printf "%s\t%s\t%s\t%d\t%.4f\t%.2f\t%.2f\t%.1f\t%d\t%.4f\t%s\n", \
                   fields[1], fields[2], fields[3], n[key], mean, \
                   min[key], max[key], cpu[key] / n[key], rss[key], \
                   baseline / mean, exact[key]
        }
    }
' "$all_runs" | {
    IFS= read -r header
    printf '%s\n' "$header"
    sort -t $'\t' -k1,1 -k2,2 -k3,3n
} >"$summary"

{
    printf 'finished\t%s\n' "$(date --iso-8601=seconds)"
    if [[ -f "$FULL_BGZF.csi" ]]; then
        printf 'bgzf_csi_sha256_after\t%s\n' \
            "$(sha256sum "$FULL_BGZF.csi" | cut -d' ' -f1)"
    fi
    if [[ -f "$FULL_BGZF.tbi" ]]; then
        printf 'bgzf_tbi_sha256_after\t%s\n' \
            "$(sha256sum "$FULL_BGZF.tbi" | cut -d' ' -f1)"
    fi
    printf 'bcf_csi_sha256_after\t%s\n' \
        "$(sha256sum "$FULL_BCF.csi" | cut -d' ' -f1)"
} >>"$manifest"

printf 'FINAL BENCHMARK COMPLETE %s\n' "$RESULT_ROOT"
