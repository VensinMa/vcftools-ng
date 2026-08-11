#!/usr/bin/env bash
set -euo pipefail

source_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

NG=${NG:-"$source_root/build-pgo-v0141-final/vcftools-ng"}
BCFTOOLS=${BCFTOOLS:-"/home/vensin/software/vcftools-ng-v0.13.0-linux-x86_64/libexec/bcftools.bin"}
FULL_BGZF=${FULL_BGZF:-"/home/vensin/workspace/Sweet_Osmanthus/05.variant_filter/02.vcftools_filter_snp_indel/412samples.SNP.biallelic.minGQ10.minQ30.meanDP6.maxmiss0.8.maf0.05.vcf.gz"}
FULL_PLAIN=${FULL_PLAIN:-"$source_root/data/osmanthus412.snps.vcf"}
FULL_BCF=${FULL_BCF:-"$source_root/data/osmanthus412.snps.bcf"}
ORACLE_ROOT=${ORACLE_ROOT:-"$source_root/benchmarks/results/final-full-v0121/golden"}
ORACLE_RUN_ROOT=${ORACLE_RUN_ROOT:-"$source_root/benchmarks/results/final-full-v0121/runs"}
THREAD_LIST=${THREAD_LIST:-"1 2 4 8 12 16 24 28 32"}
REPEATS=${REPEATS:-3}
EXPECTED_RECORDS=${EXPECTED_RECORDS:-11230392}
GATE_ONLY=${GATE_ONLY:-0}
VALIDATE_ONLY=${VALIDATE_ONLY:-0}
RESULT_ROOT=${1:-"$source_root/benchmarks/results/final-full-v0141"}

readonly required_threads="1 2 4 8 12 16 24 28 32"
readonly filter_args=(
    --min-alleles 2 --max-alleles 2
    --minGQ 10 --minQ 30
    --min-meanDP 7 --max-missing 0.9 --maf 0.1
)
readonly cases=(bgzf_tbi bgzf_auto_csi plain_vcf bcf_adaptive)
readonly formats=(gzvcf vcf bcf)

declare -Ar expected_sizes=(
    [bgzf]=18940264903
    [bgzf_tbi]=539057
    [plain]=122911664549
    [bcf]=22074294040
    [bcf_csi]=481397
    [oracle_gzvcf]=59434159204
    [oracle_vcf]=59434159621
    [oracle_bcf]=57211771106
)
declare -Ar expected_hashes=(
    [bgzf]=24d9dceeb6cf174c5112523cdc1488c55f92f8f5a484c89ac2f61807384313b4
    [bgzf_tbi]=f2dd96154a2f75f1c8bd49acc005036785ea03d8ab1d373e2c0a655f89bff371
    [plain]=22e0589561425f700f691f458994710ebec6d615c17d790dd500ad1c131d981e
    [bcf]=dfba81b4ed44985eae27eb883a0bcb1df51486ea323074d654e128221047324d
    [bcf_csi]=8b82bbae9a6977f67d122abd936a89a1cc01ab0b3869f522b53490bd39312a31
    [oracle_gzvcf]=7548416e01d4a318b81c5d1feb9429f60c7995205d66169242c3792af4c4fc14
    [oracle_vcf]=d4f2a15e8c5ad0cc12abf4a3ab308bb48f22adf8ec13776d72ceeaa1f8d402b8
    [oracle_bcf]=dde9edd98d5d05aa885e0e2f78a9696cfe0802c472d4dfd81ccffb49339107f3
)

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

[[ "$THREAD_LIST" == "$required_threads" ]] ||
    fail "release matrix requires THREAD_LIST=\"$required_threads\""
[[ "$REPEATS" == 3 ]] || fail "release matrix requires REPEATS=3"
[[ "$GATE_ONLY" == 0 || "$GATE_ONLY" == 1 ]] ||
    fail "GATE_ONLY must be 0 or 1"
[[ "$VALIDATE_ONLY" == 0 || "$VALIDATE_ONLY" == 1 ]] ||
    fail "VALIDATE_ONLY must be 0 or 1"
for executable in "$NG" "$BCFTOOLS"; do
    [[ -x "$executable" ]] || fail "missing executable: $executable"
done
bcftools_version=$($BCFTOOLS --version | head -n 2)
grep -Fqx 'bcftools 1.24' <<<"$bcftools_version" ||
    fail "automatic-CSI release gate requires bcftools 1.24"
grep -Fqx 'Using htslib 1.24' <<<"$bcftools_version" ||
    fail "automatic-CSI release gate requires HTSlib 1.24"
for input in "$FULL_BGZF" "$FULL_BGZF.tbi" "$FULL_PLAIN" \
             "$FULL_BCF" "$FULL_BCF.csi"; do
    [[ -f "$input" ]] || fail "missing locked input or index: $input"
done
for format in "${formats[@]}"; do
    [[ -f "$ORACLE_ROOT/$format.vcf" ]] ||
        fail "missing locked Original oracle: $ORACLE_ROOT/$format.vcf"
    [[ -s "$ORACLE_RUN_ROOT/original-$format-r1.tsv" ]] ||
        fail "missing locked Original timing: $ORACLE_RUN_ROOT/original-$format-r1.tsv"
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

ng_version=$($NG --version)
[[ "$ng_version" == "vcftools-ng 0.14.1" ]] ||
    fail "candidate must report vcftools-ng 0.14.1, got: $ng_version"

asset_path() {
    case "$1" in
        bgzf) printf '%s\n' "$FULL_BGZF" ;;
        bgzf_tbi) printf '%s\n' "$FULL_BGZF.tbi" ;;
        plain) printf '%s\n' "$FULL_PLAIN" ;;
        bcf) printf '%s\n' "$FULL_BCF" ;;
        bcf_csi) printf '%s\n' "$FULL_BCF.csi" ;;
        oracle_gzvcf) printf '%s\n' "$ORACLE_ROOT/gzvcf.vcf" ;;
        oracle_vcf) printf '%s\n' "$ORACLE_ROOT/vcf.vcf" ;;
        oracle_bcf) printf '%s\n' "$ORACLE_ROOT/bcf.vcf" ;;
        *) fail "unknown locked asset: $1" ;;
    esac
}

validate_assets() {
    local report="$RESULT_ROOT/asset-validation.tsv"
    local temporary="$report.tmp.$$"
    printf 'asset\tpath\tbytes\tsha256\tstatus\n' >"$temporary"
    local key path bytes digest
    for key in bgzf bgzf_tbi plain bcf bcf_csi \
               oracle_gzvcf oracle_vcf oracle_bcf; do
        path=$(asset_path "$key")
        bytes=$(stat -c '%s' "$path")
        [[ "$bytes" == "${expected_sizes[$key]}" ]] ||
            fail "$key size changed: expected ${expected_sizes[$key]}, got $bytes"
        printf 'HASH %s (%s bytes)\n' "$key" "$bytes"
        digest=$(sha256sum "$path" | awk '{print $1}')
        [[ "$digest" == "${expected_hashes[$key]}" ]] ||
            fail "$key SHA-256 changed: expected ${expected_hashes[$key]}, got $digest"
        printf '%s\t%s\t%s\t%s\tPASS\n' \
            "$key" "$path" "$bytes" "$digest" >>"$temporary"
    done
    mv -f -- "$temporary" "$report"
}

completed_run() {
    local metadata=$1
    [[ -s "$metadata" ]] &&
        awk -F '\t' 'NR == 2 && $12 == 0 && $13 ~ /^(PASS|ORACLE)$/ {
                         ok=1
                     }
                     END {exit !ok}' "$metadata"
}

write_metadata() {
    local metadata=$1 scenario=$2 format=$3 engine=$4 threads=$5 repeat=$6
    local timing=$7 output=$8 exit_status=$9 exact=${10} backend=${11}
    local output_hash=${12} kept=${13} total=${14} index_state=${15}
    local wall=NA user=NA system=NA cpu=NA rss=NA bytes=0
    if [[ -s "$timing" ]]; then
        IFS=$'\t' read -r wall user system cpu rss <"$timing"
        cpu=${cpu%\%}
    fi
    [[ -e "$output" ]] && bytes=$(stat -c '%s' "$output")
    local temporary="$metadata.tmp.$$"
    printf 'scenario\tformat\tengine\tthreads\trepeat\twall_s\tuser_s\tsystem_s\tcpu_pct\tmax_rss_kb\toutput_bytes\texit\texact\toutput_sha256\tbackend\tkept_sites\ttotal_sites\tindex_state\n' >"$temporary"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$scenario" "$format" "$engine" "$threads" "$repeat" \
        "$wall" "$user" "$system" "$cpu" "$rss" "$bytes" \
        "$exit_status" "$exact" "$output_hash" "$backend" "$kept" \
        "$total" "$index_state" >>"$temporary"
    mv -f -- "$temporary" "$metadata"
}

link_locked_golden() {
    local format=$1 target="$RESULT_ROOT/golden/$1.vcf"
    if [[ ! -e "$target" ]]; then
        ln -s "$ORACLE_ROOT/$format.vcf" "$target"
    fi
    [[ $(readlink -f "$target") == "$ORACLE_ROOT/$format.vcf" ]] ||
        fail "result golden does not reference the locked oracle: $target"
}

lock_original_baseline() {
    local format=$1
    link_locked_golden "$format"
    local metadata="$RESULT_ROOT/runs/original-locked-$format.tsv"
    local locked="$ORACLE_RUN_ROOT/original-$format-r1.tsv"
    if ! awk -F '\t' 'NR == 2 && $3 == "original" && $12 == 0 &&
                             $13 == "ORACLE" {ok=1}
                         END {exit !ok}' "$locked"; then
        fail "invalid locked Original metadata: $locked"
    fi
    IFS=$'\t' read -r _ _ _ _ _ wall user system cpu rss bytes \
        _ _ digest < <(sed -n '2p' "$locked")
    local timing="$RESULT_ROOT/runs/original-locked-$format.time.txt"
    printf '%s\t%s\t%s\t%s%%\t%s\n' \
        "$wall" "$user" "$system" "$cpu" "$rss" >"$timing"
    write_metadata "$metadata" "$format" "$format" original 1 1 \
        "$timing" "$ORACLE_ROOT/$format.vcf" 0 ORACLE locked-v0.12.1 \
        "$digest" NA "$EXPECTED_RECORDS" retained
}

case_parameters() {
    extra=()
    unique_input=no
    case "$1" in
        bgzf_tbi)
            kind=gzvcf; source_input=$FULL_BGZF ;;
        bgzf_auto_csi)
            kind=gzvcf; source_input=$FULL_BGZF; unique_input=yes
            extra=(--bcftools "$BCFTOOLS") ;;
        plain_vcf)
            kind=vcf; source_input=$FULL_PLAIN ;;
        bcf_adaptive)
            kind=bcf; source_input=$FULL_BCF ;;
        *) fail "unknown release scenario: $1" ;;
    esac
}

assert_backend() {
    local scenario=$1 threads=$2 input=$3 backend=$4
    case "$scenario" in
        bgzf_tbi)
            if ((threads == 1)); then
                [[ "$backend" == fast-filter-recode-bgzf* ]]
            else
                [[ "$backend" == fast-filter-recode-indexed-bgzf* ]]
            fi ;;
        bgzf_auto_csi)
            if ((threads == 1)); then
                [[ "$backend" == fast-filter-recode-bgzf* && ! -e "$input.csi" ]]
            else
                [[ "$backend" == fast-filter-recode-indexed-bgzf* && -s "$input.csi" ]]
            fi ;;
        plain_vcf)
            [[ "$backend" == fast-filter-recode-plain* ]] ;;
        bcf_adaptive)
            [[ "$backend" == stream* ]] ;;
    esac
}

ensure_scratch_capacity() {
    local expected=$1 available
    available=$(df -PB1 "$RESULT_ROOT/scratch" | awk 'NR == 2 {print $4}')
    ((available >= expected + 5 * 1024 * 1024 * 1024)) ||
        fail "insufficient scratch space: $available bytes available, $expected-byte output expected"
}

run_candidate() {
    local scenario=$1 threads=$2 repeat=$3
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
    local program_log="$RESULT_ROOT/logs/$run_id.log"
    local stderr_log="$RESULT_ROOT/logs/$run_id.stderr.txt"
    local status=0 exact=FAIL backend=missing output_hash=missing
    local kept=NA total=NA index_state=none
    ensure_scratch_capacity "${expected_sizes[oracle_$kind]}"
    printf 'START %s %s\n' "$run_id" "$(date --iso-8601=seconds)"
    /usr/bin/time -f '%e\t%U\t%S\t%P\t%M' -o "$timing" \
        "$NG" "--$kind" "$input" --threads "$threads" \
        "${extra[@]}" "${filter_args[@]}" \
        --log-file "$program_log" \
        --recode --recode-INFO-all --stdout \
        >"$output" 2>"$stderr_log" || status=$?
    backend=$(sed -n 's/^Input backend: //p' "$program_log" | head -1)
    read -r kept total < <(
        sed -n 's/^After filtering, kept \([0-9][0-9]*\) out of \([0-9][0-9]*\) sites$/\1 \2/p' \
            "$program_log" | tail -1
    ) || true
    kept=${kept:-NA}
    total=${total:-NA}
    [[ -s "$input.csi" ]] && index_state=created-csi
    if ((status == 0)) && [[ "$total" == "$EXPECTED_RECORDS" ]] &&
       cmp -s "$ORACLE_ROOT/$kind.vcf" "$output" &&
       assert_backend "$scenario" "$threads" "$input" "$backend"; then
        exact=PASS
        output_hash=${expected_hashes[oracle_$kind]}
    fi
    write_metadata "$metadata" "$scenario" "$kind" vcftools-ng \
        "$threads" "$repeat" "$timing" "$output" "$status" "$exact" \
        "$backend" "$output_hash" "$kept" "$total" "$index_state"
    if ((status != 0)) || [[ "$exact" != PASS ]]; then
        printf 'FAILED %s exit=%s exact=%s backend=%s kept=%s total=%s\n' \
            "$run_id" "$status" "$exact" "$backend" "$kept" "$total" >&2
        printf 'Scratch output retained for diagnosis: %s\n' "$scratch" >&2
        exit 1
    fi
    rm -rf -- "$scratch"
    printf 'DONE %s exact=PASS backend=%s %s\n' \
        "$run_id" "${backend%% (*}" "$(date --iso-8601=seconds)"
}

refresh_reports() {
    local all_runs="$RESULT_ROOT/all-runs.tsv"
    head -1 "$RESULT_ROOT/runs/original-locked-gzvcf.tsv" >"$all_runs"
    find "$RESULT_ROOT/runs" -maxdepth 1 -type f -name '*.tsv' -print0 |
        sort -z | xargs -0 -r -n1 tail -1 >>"$all_runs"
    awk -F '\t' '
        BEGIN {
            OFS="\t"
            print "scenario", "threads", "runs", "original_single_wall_s", \
                  "ng_mean_wall_s", "speedup", "mean_cpu_pct", \
                  "max_rss_kb", "backend", "exact"
        }
        NR == 1 {next}
        $3 == "original" {original_wall[$2]=$6; next}
        $3 == "vcftools-ng" {
            key=$1 SUBSEP $4; format[key]=$2; n[key]++; wall[key]+=$6
            cpu[key]+=$9
            if (!(key in rss) || $10 > rss[key]) rss[key]=$10
            if (!(key in backend)) backend[key]=$15
            else if (backend[key] != $15) backend[key]="mixed"
            if ($13 != "PASS") exact[key]="FAIL"
            else if (!(key in exact)) exact[key]="PASS"
        }
        END {
            for (key in n) {
                split(key, f, SUBSEP); base=original_wall[format[key]]
                mean=wall[key]/n[key]
                printf "%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.1f\t%d\t%s\t%s\n", \
                    f[1], f[2], n[key], base, mean, base/mean, \
                    cpu[key]/n[key], rss[key], backend[key], exact[key]
            }
        }
    ' "$all_runs" | {
        IFS= read -r header
        printf '%s\n' "$header"
        sort -t $'\t' -k1,1 -k2,2n
    } >"$RESULT_ROOT/summary.tsv"
}

validate_assets
for format in "${formats[@]}"; do
    lock_original_baseline "$format"
done

manifest="$RESULT_ROOT/manifest.tsv"
if [[ ! -s "$manifest" ]]; then
    {
        printf 'key\tvalue\n'
        printf 'started\t%s\n' "$(date --iso-8601=seconds)"
        printf 'release\tv0.14.1\n'
        printf 'title\tCapability-Planned Exact Analytics\n'
        printf 'git_commit\t%s\n' "$(git -C "$source_root" rev-parse HEAD)"
        printf 'tracked_diff_sha256\t%s\n' "$(git -C "$source_root" diff --binary | sha256sum | awk '{print $1}')"
        printf 'source_tree_sha256\t%s\n' "$(
            cd "$source_root"
            find CMakeLists.txt VERSION src -type f -print0 |
                LC_ALL=C sort -z | xargs -0 sha256sum | sha256sum |
                awk '{print $1}'
        )"
        printf 'candidate\t%s\n' "$NG"
        printf 'candidate_sha256\t%s\n' "$(sha256sum "$NG" | awk '{print $1}')"
        printf 'ng_version\t%s\n' "$ng_version"
        printf 'bcftools\t%s\n' "$(head -1 <<<"$bcftools_version")"
        printf 'bcftools_sha256\t%s\n' "$(sha256sum "$BCFTOOLS" | awk '{print $1}')"
        printf 'original\tVCFtools 0.1.17, locked v0.12.1 baseline\n'
        printf 'original_runs_v0141\t0\n'
        printf 'original_policy\thash-validated retained oracle and single-run timing; not rerun\n'
        printf 'vcftools_ng_max_repeats\t%s\n' "$REPEATS"
        printf 'expected_records\t%s\n' "$EXPECTED_RECORDS"
        printf 'threads\t%s\n' "$THREAD_LIST"
        printf 'scenarios\tBGZF VCF + TBI;BGZF VCF + automatic CSI;Plain VCF;BCF adaptive stream\n'
        printf 'execution\tstrictly_serial\n'
        printf 'cache\toperating_system_cache_not_flushed\n'
        printf 'workload\t%s --recode --recode-INFO-all --stdout\n' "${filter_args[*]}"
        printf 'asset_validation\t%s\n' "$RESULT_ROOT/asset-validation.tsv"
        printf 'host\t%s\n' "$(hostname)"
        printf 'os\t%s\n' "$(. /etc/os-release; printf '%s' "$PRETTY_NAME")"
        printf 'kernel\t%s\n' "$(uname -srmo)"
        printf 'cpu_model\t%s\n' "$(lscpu | awk -F: '/Model name/ {sub(/^[[:space:]]+/, "", $2); print $2; exit}')"
        printf 'logical_cpus\t%s\n' "$(nproc)"
        printf 'memory\t%s\n' "$(free -h | awk '/^Mem:/ {print $2}')"
    } >"$manifest"
fi

if [[ "$VALIDATE_ONLY" == 1 ]]; then
    printf 'LOCKED ASSET VALIDATION PASS %s\n' "$RESULT_ROOT"
    exit 0
fi

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
if ! grep -q $'^first_repeat_gates_passed\t' "$manifest"; then
    printf 'first_repeat_gates_passed\t%s\n' \
        "$(date --iso-8601=seconds)" >>"$manifest"
fi
if [[ "$GATE_ONLY" == 1 ]]; then
    printf 'FIRST-REPEAT v0.14.1 RELEASE GATES PASS (36/36) %s\n' "$RESULT_ROOT"
    printf 'Adaptive repeats remain pending; resume with GATE_ONLY=0.\n'
    exit 0
fi

for scenario in "${cases[@]}"; do
    for threads in $THREAD_LIST; do
        run_id="$scenario-vcftools-ng-t$threads-r1"
        first_wall=$(awk -F '\t' 'NR == 2 {print $6}' "$RESULT_ROOT/runs/$run_id.tsv")
        awk -v wall="$first_wall" 'BEGIN {exit !(wall <= 1800)}' || continue
        run_candidate "$scenario" "$threads" 2
        second_wall=$(awk -F '\t' 'NR == 2 {print $6}' \
            "$RESULT_ROOT/runs/$scenario-vcftools-ng-t$threads-r2.tsv")
        if awk -v a="$first_wall" -v b="$second_wall" 'BEGIN {
                   d=(a>b ? a-b : b-a); mean=(a+b)/2
                   exit !((a > 600 || b > 600) && d/mean < 0.10)
               }'; then
            printf 'SKIP third repeat %s t%s: >600 s and first two differ <10%%\n' \
                "$scenario" "$threads"
            continue
        fi
        run_candidate "$scenario" "$threads" 3
    done
done

refresh_reports
if ! grep -q $'^finished\t' "$manifest"; then
    printf 'finished\t%s\n' "$(date --iso-8601=seconds)" >>"$manifest"
fi
printf 'FINAL v0.14.1 FOUR-SCENARIO RELEASE MATRIX PASS %s\n' "$RESULT_ROOT"
