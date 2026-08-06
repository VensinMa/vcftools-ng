#!/usr/bin/env bash
set -euo pipefail

source_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)

PACKAGE_ROOT=${PACKAGE_ROOT:-"/home/vensin/software/vcftools-ng-v0.13.0-linux-x86_64"}
NG=${NG:-"$source_root/build-v0130/vcftools-ng"}
NG_CORE=${NG_CORE:-"$NG"}
BCFTOOLS=${BCFTOOLS:-"$source_root/benchmarks/run-bundled-bcftools.sh"}
BCFTOOLS_PRIVATE=${BCFTOOLS_PRIVATE:-"$PACKAGE_ROOT/libexec/bcftools.bin"}
PRIVATE_LIBDIR=${PRIVATE_LIBDIR:-"$PACKAGE_ROOT/lib"}
SSD_INPUT=${SSD_INPUT:-"$source_root/data/osmanthus412.snps.vcf.gz"}
SSD_INPUT_INDEX=${SSD_INPUT_INDEX:-"$SSD_INPUT.tbi"}
SSD_PLAIN_INPUT=${SSD_PLAIN_INPUT:-"$source_root/data/osmanthus412.snps.vcf"}
SSD_BCF_INPUT=${SSD_BCF_INPUT:-"$source_root/data/osmanthus412.snps.bcf"}
SSD_BCF_INDEX=${SSD_BCF_INDEX:-"$SSD_BCF_INPUT.csi"}
SSD_OUTPUT_ROOT=${SSD_OUTPUT_ROOT:-"/home/vensin/workspace/vcftools-ng-benchmark-v0130/ssd"}
HDD_ROOT=${HDD_ROOT:-"/home/data/vcftools-ng-benchmark-v0130"}
HDD_INPUT=${HDD_INPUT:-"$HDD_ROOT/input/osmanthus412.snps.vcf.gz"}
HDD_INPUT_INDEX=${HDD_INPUT_INDEX:-"$HDD_INPUT.tbi"}
HDD_PLAIN_INPUT=${HDD_PLAIN_INPUT:-"$HDD_ROOT/input/osmanthus412.snps.vcf"}
HDD_BCF_INPUT=${HDD_BCF_INPUT:-"$HDD_ROOT/input/osmanthus412.snps.bcf"}
HDD_BCF_INDEX=${HDD_BCF_INDEX:-"$HDD_BCF_INPUT.csi"}
HDD_OUTPUT_ROOT=${HDD_OUTPUT_ROOT:-"$HDD_ROOT/output"}
ORACLE_RUN_ROOT=${ORACLE_RUN_ROOT:-"$source_root/benchmarks/results/final-full-v0121/runs"}
ORACLE_ROOT=${ORACLE_ROOT:-"$source_root/benchmarks/results/final-full-v0121/golden"}
RESULT_ROOT=${1:-"$source_root/benchmarks/results/v0130-input-output-storage"}
THREAD_LIST=${THREAD_LIST:-"1 2 4 8 12 16 24 28 32"}
REPEATS=${REPEATS:-3}
GATE_ONLY=${GATE_ONLY:-0}
RUN_SMOKE=${RUN_SMOKE:-1}
SMOKE_ONLY=${SMOKE_ONLY:-0}
LONG_RUN_SECONDS=${LONG_RUN_SECONDS:-600}
LONG_RUN_VARIATION_PCT=${LONG_RUN_VARIATION_PCT:-10}
SINGLE_RUN_SECONDS=${SINGLE_RUN_SECONDS:-1800}

readonly expected_input_bytes=18940264903
readonly expected_input_sha256=24d9dceeb6cf174c5112523cdc1488c55f92f8f5a484c89ac2f61807384313b4
readonly expected_index_bytes=539057
readonly expected_index_sha256=f2dd96154a2f75f1c8bd49acc005036785ea03d8ab1d373e2c0a655f89bff371
readonly expected_plain_bytes=122911664549
readonly expected_plain_sha256=22e0589561425f700f691f458994710ebec6d615c17d790dd500ad1c131d981e
readonly expected_bcf_bytes=22074294040
readonly expected_bcf_sha256=dfba81b4ed44985eae27eb883a0bcb1df51486ea323074d654e128221047324d
readonly expected_bcf_index_bytes=481397
readonly expected_bcf_index_sha256=8b82bbae9a6977f67d122abd936a89a1cc01ab0b3869f522b53490bd39312a31
readonly expected_total_sites=11230392
readonly expected_kept_sites=5425725
readonly required_threads="1 2 4 8 12 16 24 28 32"
readonly -a input_scenarios=(bgzf_tbi bgzf_auto_csi plain_vcf bcf_adaptive)
readonly -a output_scenarios=(ssd_plain ssd_bgzf hdd_bgzf)
readonly -a filter_args=(
    --min-alleles 2 --max-alleles 2
    --minGQ 10 --minQ 30
    --min-meanDP 7 --max-missing 0.9 --maf 0.1
    --recode-INFO-all
)

fail() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}

[[ "$THREAD_LIST" == "$required_threads" ]] ||
    fail "release matrix requires THREAD_LIST=\"$required_threads\""
[[ "$REPEATS" == 3 ]] || fail "release matrix requires REPEATS=3"
[[ "$LONG_RUN_SECONDS" == 600 ]] ||
    fail "release matrix requires LONG_RUN_SECONDS=600"
[[ "$LONG_RUN_VARIATION_PCT" == 10 ]] ||
    fail "release matrix requires LONG_RUN_VARIATION_PCT=10"
[[ "$SINGLE_RUN_SECONDS" == 1800 ]] ||
    fail "release matrix requires SINGLE_RUN_SECONDS=1800"
[[ "$GATE_ONLY" == 0 || "$GATE_ONLY" == 1 ]] ||
    fail "GATE_ONLY must be 0 or 1"
[[ "$RUN_SMOKE" == 0 || "$RUN_SMOKE" == 1 ]] ||
    fail "RUN_SMOKE must be 0 or 1"
[[ "$SMOKE_ONLY" == 0 || "$SMOKE_ONLY" == 1 ]] ||
    fail "SMOKE_ONLY must be 0 or 1"
[[ -x "$NG" ]] || fail "vcftools-ng executable is missing: $NG"
[[ -x "$NG_CORE" ]] || fail "vcftools-ng core executable is missing: $NG_CORE"
[[ -x "$BCFTOOLS" ]] || fail "bcftools 1.24 executable is missing: $BCFTOOLS"
[[ -x "$BCFTOOLS_PRIVATE" ]] ||
    fail "private bcftools 1.24 executable is missing: $BCFTOOLS_PRIVATE"
[[ -d "$PRIVATE_LIBDIR" ]] || fail "portable private library directory is missing"
"$BCFTOOLS" --version | sed -n '1p' | grep -Fqx 'bcftools 1.24' ||
    fail "automatic CSI benchmark requires bcftools 1.24"
for path in \
    "$SSD_INPUT" "$SSD_INPUT_INDEX" "$SSD_PLAIN_INPUT" \
    "$SSD_BCF_INPUT" "$SSD_BCF_INDEX" \
    "$ORACLE_RUN_ROOT/original-gzvcf-r1.tsv" \
    "$ORACLE_RUN_ROOT/original-vcf-r1.tsv" \
    "$ORACLE_RUN_ROOT/original-bcf-r1.tsv" \
    "$ORACLE_ROOT/gzvcf.vcf" "$ORACLE_ROOT/vcf.vcf" \
    "$ORACLE_ROOT/bcf.vcf"; do
    [[ -r "$path" ]] || fail "required input is missing: $path"
done

NG=$(realpath "$NG")
NG_CORE=$(realpath "$NG_CORE")
BCFTOOLS=$(realpath "$BCFTOOLS")
BCFTOOLS_PRIVATE=$(realpath "$BCFTOOLS_PRIVATE")
PRIVATE_LIBDIR=$(realpath "$PRIVATE_LIBDIR")
SSD_INPUT=$(realpath "$SSD_INPUT")
SSD_INPUT_INDEX=$(realpath "$SSD_INPUT_INDEX")
SSD_PLAIN_INPUT=$(realpath "$SSD_PLAIN_INPUT")
SSD_BCF_INPUT=$(realpath "$SSD_BCF_INPUT")
SSD_BCF_INDEX=$(realpath "$SSD_BCF_INDEX")
ORACLE_RUN_ROOT=$(realpath "$ORACLE_RUN_ROOT")
ORACLE_ROOT=$(realpath "$ORACLE_ROOT")
mkdir -p "$RESULT_ROOT" "$RESULT_ROOT"/{gates,golden,logs,runs,smoke}
mkdir -p "$SSD_OUTPUT_ROOT" "$HDD_ROOT/input" "$HDD_OUTPUT_ROOT"
RESULT_ROOT=$(realpath "$RESULT_ROOT")
SSD_OUTPUT_ROOT=$(realpath "$SSD_OUTPUT_ROOT")
HDD_ROOT=$(realpath "$HDD_ROOT")
HDD_OUTPUT_ROOT=$(realpath "$HDD_OUTPUT_ROOT")

declare -A oracle_wall=() oracle_bytes=() oracle_sha256=() oracle_vcf=()
load_oracle() {
    local key=$1
    local metadata="$ORACLE_RUN_ROOT/original-$key-r1.tsv"
    local path="$ORACLE_ROOT/$key.vcf"
    local exit_status status
    IFS=$'\t' read -r _ _ _ _ _ oracle_wall["$key"] _ _ _ _ \
        oracle_bytes["$key"] exit_status status oracle_sha256["$key"] \
        < <(sed -n '2p' "$metadata")
    [[ "$exit_status" == 0 && "$status" == ORACLE ]] ||
        fail "locked Original metadata is invalid: $metadata"
    [[ "${oracle_bytes[$key]}" == "$(stat -c '%s' "$path")" ]] ||
        fail "locked Original VCF size changed: $path"
    oracle_vcf[$key]=$path
}
for key in gzvcf vcf bcf; do
    load_oracle "$key"
done

same_device() {
    [[ "$(stat -c '%d' "$1")" == "$(stat -c '%d' "$2")" ]]
}

mount_description() {
    findmnt -no SOURCE,FSTYPE,TARGET -T "$1" | tr -s ' '
}

device_rotational() {
    local source
    source=$(findmnt -no SOURCE -T "$1")
    local parent
    parent=$(lsblk -no PKNAME "$source" 2>/dev/null | head -1)
    if [[ -n "$parent" ]]; then
        cat "/sys/class/block/$parent/queue/rotational" 2>/dev/null ||
            printf 'unknown\n'
    else
        printf 'unknown\n'
    fi
}

available_bytes() {
    df --output=avail -B1 "$1" | awk 'NR == 2 {print $1}'
}

import_single_input_results() {
    local old="$source_root/benchmarks/results/v0130-output-storage"
    [[ -d "$old" ]] || return
    if [[ ! -s "$RESULT_ROOT/hdd-input-bgzf.tsv" &&
          -s "$old/hdd-input-validation.tsv" ]]; then
        cp -- "$old/hdd-input-validation.tsv" \
            "$RESULT_ROOT/hdd-input-bgzf.tsv"
    fi
    if [[ ! -s "$RESULT_ROOT/gates/output-smoke.pass" &&
          -s "$old/gates/smoke.pass" ]]; then
        cp -a -- "$old/smoke/." "$RESULT_ROOT/smoke/"
        cp -- "$old/gates/smoke.pass" \
            "$RESULT_ROOT/gates/output-smoke.pass"
    fi
    local threads old_id new_id
    for threads in $THREAD_LIST; do
        old_id="ssd_plain-t$threads-r1"
        new_id="bgzf_tbi__ssd_plain-t$threads-r1"
        if [[ ! -s "$RESULT_ROOT/runs/$new_id.tsv" &&
              -s "$old/runs/$old_id.tsv" ]]; then
            awk -F '\t' -v OFS='\t' '
                NR == 1 {print; next}
                {$1 = "bgzf_tbi__ssd_plain"; print}
            ' "$old/runs/$old_id.tsv" >"$RESULT_ROOT/runs/$new_id.tsv"
            suffix=time.txt
            [[ ! -e "$old/runs/$old_id.$suffix" ]] ||
                cp -- "$old/runs/$old_id.$suffix" \
                    "$RESULT_ROOT/runs/$new_id.$suffix"
            for suffix in log stderr.txt stdout.txt; do
                [[ ! -e "$old/logs/$old_id.$suffix" ]] ||
                    cp -- "$old/logs/$old_id.$suffix" \
                        "$RESULT_ROOT/logs/$new_id.$suffix"
            done
            [[ ! -e "$old/gates/ssd_plain-t$threads.raw.sha256" ]] ||
                cp -- "$old/gates/ssd_plain-t$threads.raw.sha256" \
                    "$RESULT_ROOT/gates/bgzf_tbi__ssd_plain-t$threads.raw.sha256"
        fi
    done
    if [[ -s "$old/gates/ssd_plain.pass" ]]; then
        cp -- "$old/gates/ssd_plain.pass" \
            "$RESULT_ROOT/gates/bgzf_tbi__ssd_plain.pass"
    fi
}

prepare_retained_file() {
    local label=$1 source=$2 target=$3 expected_bytes=$4 expected_hash=$5
    local validation="$RESULT_ROOT/hdd-input-$label.tsv"
    if [[ -s "$validation" ]] &&
       awk -F '\t' '$1 == "status" && $2 == "PASS" {ok=1}
            END {exit !ok}' "$validation" &&
       [[ -s "$target" ]] &&
       [[ "$(stat -c '%s' "$target")" == "$expected_bytes" ]]; then
        printf 'REUSE hash-validated HDD file: %s\n' "$target"
        return
    fi
    [[ "$(stat -c '%s' "$source")" == "$expected_bytes" ]] ||
        fail "$label source size changed: $source"
    printf 'Validating locked %s source SHA-256 once...\n' "$label"
    local source_hash
    source_hash=$(sha256sum "$source" | awk '{print $1}')
    [[ "$source_hash" == "$expected_hash" ]] ||
        fail "$label source SHA-256 changed: $source_hash"
    if [[ ! -e "$target" ]]; then
        local partial="$target.vcftools-ng-copy-partial"
        [[ ! -e "$partial" ]] ||
            fail "incomplete prior HDD copy requires inspection: $partial"
        printf 'Copying retained %s to HDD once...\n' "$label"
        cp --reflink=never -- "$source" "$partial"
        sync "$partial"
        mv -- "$partial" "$target"
        sync "$target"
    fi
    [[ "$(stat -c '%s' "$target")" == "$expected_bytes" ]] ||
        fail "$label HDD copy has the wrong size"
    printf 'Validating retained %s HDD SHA-256 once...\n' "$label"
    local target_hash
    target_hash=$(sha256sum "$target" | awk '{print $1}')
    [[ "$target_hash" == "$expected_hash" ]] ||
        fail "$label HDD SHA-256 mismatch: $target_hash"
    {
        printf 'key\tvalue\n'
        printf 'status\tPASS\n'
        printf 'validated\t%s\n' "$(date --iso-8601=seconds)"
        printf 'label\t%s\n' "$label"
        printf 'source\t%s\n' "$source"
        printf 'source_bytes\t%s\n' "$expected_bytes"
        printf 'source_sha256\t%s\n' "$source_hash"
        printf 'retained_hdd_file\t%s\n' "$target"
        printf 'retained_hdd_sha256\t%s\n' "$target_hash"
        printf 'copy_policy\tone-time; retained; never overwritten automatically\n'
    } >"$validation"
    printf 'HDD FILE VALIDATION PASS: %s\n' "$target"
}

if [[ ${IMPORT_SINGLE_INPUT_RESULTS:-0} == 1 ]]; then
    import_single_input_results
fi
prepare_retained_file bgzf "$SSD_INPUT" "$HDD_INPUT" \
    "$expected_input_bytes" "$expected_input_sha256"
prepare_retained_file bgzf-tbi "$SSD_INPUT_INDEX" "$HDD_INPUT_INDEX" \
    "$expected_index_bytes" "$expected_index_sha256"
prepare_retained_file plain-vcf "$SSD_PLAIN_INPUT" "$HDD_PLAIN_INPUT" \
    "$expected_plain_bytes" "$expected_plain_sha256"
prepare_retained_file bcf "$SSD_BCF_INPUT" "$HDD_BCF_INPUT" \
    "$expected_bcf_bytes" "$expected_bcf_sha256"
prepare_retained_file bcf-csi "$SSD_BCF_INDEX" "$HDD_BCF_INDEX" \
    "$expected_bcf_index_bytes" "$expected_bcf_index_sha256"
HDD_INPUT=$(realpath "$HDD_INPUT")
HDD_INPUT_INDEX=$(realpath "$HDD_INPUT_INDEX")
HDD_PLAIN_INPUT=$(realpath "$HDD_PLAIN_INPUT")
HDD_BCF_INPUT=$(realpath "$HDD_BCF_INPUT")
HDD_BCF_INDEX=$(realpath "$HDD_BCF_INDEX")
for path in "$SSD_INPUT" "$SSD_PLAIN_INPUT" "$SSD_BCF_INPUT"; do
    same_device "$path" "$SSD_OUTPUT_ROOT" ||
        fail "SSD input and output are not on the same filesystem/device: $path"
done
for path in "$HDD_INPUT" "$HDD_PLAIN_INPUT" "$HDD_BCF_INPUT"; do
    same_device "$path" "$HDD_OUTPUT_ROOT" ||
        fail "HDD input and output are not on the same filesystem/device: $path"
done
[[ "$(device_rotational "$SSD_OUTPUT_ROOT")" == 0 ]] ||
    fail "SSD output target is not reported as non-rotational"
[[ "$(device_rotational "$HDD_OUTPUT_ROOT")" == 1 ]] ||
    fail "HDD output target is not reported as rotational"

# The largest exact output is about 59.4 GB. Keep a conservative margin for
# the staged transactional file and filesystem headroom.
(( $(available_bytes "$SSD_OUTPUT_ROOT") > 80000000000 )) ||
    fail "SSD output target has less than 80 GB available"
(( $(available_bytes "$HDD_OUTPUT_ROOT") > 80000000000 )) ||
    fail "HDD output target has less than 80 GB available"

output_parameters() {
    local output_scenario=$1
    case "$output_scenario" in
        ssd_plain)
            storage=SSD
            output_root=$SSD_OUTPUT_ROOT
            output_format=VCF
            recode_option=--recode-vcf
            output_suffix=.recode.vcf
            ;;
        ssd_bgzf)
            storage=SSD
            output_root=$SSD_OUTPUT_ROOT
            output_format='BGZF VCF'
            recode_option=--recode-vcf-gz
            output_suffix=.recode.vcf.gz
            ;;
        hdd_bgzf)
            storage=HDD
            output_root=$HDD_OUTPUT_ROOT
            output_format='BGZF VCF'
            recode_option=--recode-vcf-gz
            output_suffix=.recode.vcf.gz
            ;;
        *) fail "unknown output scenario: $output_scenario" ;;
    esac
    filesystem=$(findmnt -no FSTYPE -T "$output_root")
}

input_parameters() {
    local input_scenario=$1
    unique_auto_index=no
    input_extra=()
    case "$input_scenario" in
        bgzf_tbi)
            input_kind=gzvcf
            oracle_key=gzvcf
            if [[ "$storage" == SSD ]]; then
                base_input=$SSD_INPUT
                required_sidecar=$SSD_INPUT_INDEX
            else
                base_input=$HDD_INPUT
                required_sidecar=$HDD_INPUT_INDEX
            fi
            ;;
        bgzf_auto_csi)
            input_kind=gzvcf
            oracle_key=gzvcf
            unique_auto_index=yes
            input_extra=(--bcftools "$BCFTOOLS")
            required_sidecar=
            if [[ "$storage" == SSD ]]; then
                base_input=$SSD_INPUT
            else
                base_input=$HDD_INPUT
            fi
            ;;
        plain_vcf)
            input_kind=vcf
            oracle_key=vcf
            required_sidecar=
            if [[ "$storage" == SSD ]]; then
                base_input=$SSD_PLAIN_INPUT
            else
                base_input=$HDD_PLAIN_INPUT
            fi
            ;;
        bcf_adaptive)
            input_kind=bcf
            oracle_key=bcf
            if [[ "$storage" == SSD ]]; then
                base_input=$SSD_BCF_INPUT
                required_sidecar=$SSD_BCF_INDEX
            else
                base_input=$HDD_BCF_INPUT
                required_sidecar=$HDD_BCF_INDEX
            fi
            ;;
        *) fail "unknown input scenario: $input_scenario" ;;
    esac
    [[ -r "$base_input" ]] || fail "input is missing: $base_input"
    [[ -z "$required_sidecar" || -r "$required_sidecar" ]] ||
        fail "required input sidecar is missing: $required_sidecar"
}

seconds_between() {
    awk -v start="$1" -v finish="$2" \
        'BEGIN {printf "%.6f", (finish - start) / 1000000000}'
}

completed_run() {
    local metadata=$1
    [[ -s "$metadata" ]] &&
        awk -F '\t' 'NR == 2 && $21 == 0 && $16 == "PASS" {ok=1}
             END {exit !ok}' "$metadata"
}

run_smoke() {
    local gate="$RESULT_ROOT/gates/output-smoke.pass"
    [[ "$RUN_SMOKE" == 1 ]] || return 0
    if [[ -s "$gate" ]]; then
        printf 'REUSE completed 23,000-site smoke gate\n'
        return
    fi
    local smoke_ssd="$source_root/tests/fixtures/osmanthus205.gatk.23chr_1k.vcf.gz"
    local smoke_tbi="$smoke_ssd.tbi"
    local smoke_golden="$source_root/tests/golden/gatk205-seven-filter.recode.recode.vcf"
    for path in "$smoke_ssd" "$smoke_tbi" "$smoke_golden"; do
        [[ -r "$path" ]] || fail "smoke fixture is missing: $path"
    done
    local smoke_hdd="$HDD_ROOT/input/osmanthus205.gatk.23chr_1k.vcf.gz"
    if [[ ! -e "$smoke_hdd" ]]; then
        cp --reflink=never -- "$smoke_ssd" "$smoke_hdd"
        cp --reflink=never -- "$smoke_tbi" "$smoke_hdd.tbi"
        sync "$smoke_hdd" "$smoke_hdd.tbi"
    fi
    cmp -s "$smoke_ssd" "$smoke_hdd" || fail "HDD smoke input changed"
    cmp -s "$smoke_tbi" "$smoke_hdd.tbi" || fail "HDD smoke TBI changed"
    local expected
    expected=$(sha256sum "$smoke_golden" | awk '{print $1}')
    printf 'scenario\tthreads\tcontent_sha256\texact\n' \
        >"$RESULT_ROOT/smoke/summary.tsv"
    local scenario threads
    for scenario in "${output_scenarios[@]}"; do
        output_parameters "$scenario"
        local smoke_input=$smoke_ssd
        [[ "$scenario" == hdd_bgzf ]] && smoke_input=$smoke_hdd
        for threads in 1 8 32; do
            local prefix="$output_root/smoke-$scenario-t$threads"
            local output="$prefix$output_suffix"
            local log="$RESULT_ROOT/smoke/$scenario-t$threads.log"
            rm -f -- "$output"
            "$NG" --gzvcf "$smoke_input" --threads "$threads" \
                "${filter_args[@]}" "$recode_option" \
                --log-file "$log" --out "$prefix" \
                >"$RESULT_ROOT/smoke/$scenario-t$threads.stdout" \
                2>"$RESULT_ROOT/smoke/$scenario-t$threads.stderr"
            sync "$output"
            local actual
            if [[ "$output_format" == VCF ]]; then
                actual=$(sha256sum "$output" | awk '{print $1}')
            else
                gzip -t "$output"
                actual=$(gzip -dc "$output" | sha256sum | awk '{print $1}')
            fi
            [[ "$actual" == "$expected" ]] ||
                fail "smoke compatibility failed: $scenario t$threads"
            printf '%s\t%s\t%s\tPASS\n' \
                "$scenario" "$threads" "$actual" \
                >>"$RESULT_ROOT/smoke/summary.tsv"
            rm -f -- "$output"
        done
    done
    sync -f "$SSD_OUTPUT_ROOT"
    sync -f "$HDD_OUTPUT_ROOT"
    printf 'status\tvalidated\nPASS\t%s\n' \
        "$(date --iso-8601=seconds)" >"$gate"
    printf '23,000-SITE THREE-SCENARIO SMOKE PASS\n'
}

write_manifest() {
    local manifest="$RESULT_ROOT/manifest.tsv"
    [[ -s "$manifest" ]] && return
    {
        printf 'key\tvalue\n'
        printf 'started\t%s\n' "$(date --iso-8601=seconds)"
        printf 'candidate\tv0.13.0\n'
        printf 'title\tFour input backends by three same-device output scenarios\n'
        printf 'git_commit\t%s\n' "$(git -C "$source_root" rev-parse HEAD)"
        printf 'working_tree_sha256\t%s\n' "$(
            git -C "$source_root" diff --binary | sha256sum | awk '{print $1}'
        )"
        printf 'ng\t%s\n' "$NG"
        printf 'ng_version\t%s\n' "$("$NG" --version)"
        printf 'ng_launcher_sha256\t%s\n' "$(sha256sum "$NG" | awk '{print $1}')"
        printf 'ng_core\t%s\n' "$NG_CORE"
        printf 'ng_core_sha256\t%s\n' "$(sha256sum "$NG_CORE" | awk '{print $1}')"
        printf 'bcftools\t%s\n' "$BCFTOOLS"
        printf 'bcftools_version\t%s\n' "$("$BCFTOOLS" --version | sed -n '1p')"
        printf 'bcftools_sha256\t%s\n' "$(sha256sum "$BCFTOOLS" | awk '{print $1}')"
        printf 'bcftools_private\t%s\n' "$BCFTOOLS_PRIVATE"
        printf 'bcftools_private_sha256\t%s\n' \
            "$(sha256sum "$BCFTOOLS_PRIVATE" | awk '{print $1}')"
        printf 'threads\t%s\n' "$THREAD_LIST"
        printf 'repeats\t%s\n' "$REPEATS"
        printf 'repeat_policy\tup to 3 repeats; a first run over 1800 seconds runs once; otherwise skip repeat 3 when either of the first two application wall times exceeds 600 seconds and their symmetric difference is below 10 percent\n'
        printf 'long_run_seconds\t%s\n' "$LONG_RUN_SECONDS"
        printf 'long_run_variation_pct\t%s\n' "$LONG_RUN_VARIATION_PCT"
        printf 'single_run_seconds\t%s\n' "$SINGLE_RUN_SECONDS"
        printf 'input_scenarios\t%s\n' "${input_scenarios[*]}"
        printf 'output_scenarios\t%s\n' "${output_scenarios[*]}"
        printf 'execution\tstrictly serial\n'
        printf 'cache_policy\tOS cache not dropped; first repeat retained separately\n'
        printf 'durability\tapplication wall plus targeted sync of final output\n'
        printf 'output_index\tnot built\n'
        printf 'filters\t%s\n' "${filter_args[*]}"
        printf 'expected_sites\t%s\n' "$expected_total_sites"
        printf 'expected_kept_sites\t%s\n' "$expected_kept_sites"
        printf 'ssd_input\t%s\n' "$SSD_INPUT"
        printf 'ssd_input_sha256\t%s\n' "$expected_input_sha256"
        printf 'ssd_index\t%s\n' "$SSD_INPUT_INDEX"
        printf 'ssd_index_sha256\t%s\n' "$expected_index_sha256"
        printf 'ssd_plain_input\t%s\n' "$SSD_PLAIN_INPUT"
        printf 'ssd_plain_sha256\t%s\n' "$expected_plain_sha256"
        printf 'ssd_bcf_input\t%s\n' "$SSD_BCF_INPUT"
        printf 'ssd_bcf_sha256\t%s\n' "$expected_bcf_sha256"
        printf 'hdd_bgzf_input\t%s\n' "$HDD_INPUT"
        printf 'hdd_plain_input\t%s\n' "$HDD_PLAIN_INPUT"
        printf 'hdd_bcf_input\t%s\n' "$HDD_BCF_INPUT"
        printf 'hdd_input_policy\tone-time copied, hash-validated, retained\n'
        local key
        for key in gzvcf vcf bcf; do
            printf 'oracle_%s\t%s\n' "$key" "${oracle_vcf[$key]}"
            printf 'oracle_%s_sha256\t%s\n' "$key" "${oracle_sha256[$key]}"
            printf 'oracle_%s_bytes\t%s\n' "$key" "${oracle_bytes[$key]}"
            printf 'original_%s_wall_s_context_only\t%s\n' \
                "$key" "${oracle_wall[$key]}"
        done
        printf 'original_rerun\tno\n'
        printf 'ssd_mount\t%s\n' "$(mount_description "$SSD_OUTPUT_ROOT")"
        printf 'ssd_rotational\t%s\n' "$(device_rotational "$SSD_OUTPUT_ROOT")"
        printf 'hdd_mount\t%s\n' "$(mount_description "$HDD_OUTPUT_ROOT")"
        printf 'hdd_rotational\t%s\n' "$(device_rotational "$HDD_OUTPUT_ROOT")"
        printf 'ssd_same_device\tyes\n'
        printf 'hdd_same_device\tyes\n'
        printf 'host\t%s\n' "$(hostname)"
        printf 'kernel\t%s\n' "$(uname -srmo)"
        printf 'cpu_model\t%s\n' "$(
            lscpu | awk -F: '/Model name/ {
                sub(/^[[:space:]]+/, "", $2); print $2; exit
            }'
        )"
        printf 'logical_cpus\t%s\n' "$(nproc)"
        printf 'memory\t%s\n' "$(free -h | awk '/^Mem:/ {print $2}')"
    } >"$manifest"
}

first_run_is_single_only() {
    local input_scenario=$1
    local output_scenario=$2
    local threads=$3
    local scenario="${input_scenario}__${output_scenario}"
    local first="$RESULT_ROOT/runs/$scenario-t$threads-r1.tsv"
    completed_run "$first" || return 1
    awk -F '\t' -v limit="$SINGLE_RUN_SECONDS" \
        'NR == 2 {exit !($6 > limit)}' "$first"
}

record_single_run_skip() {
    local input_scenario=$1
    local output_scenario=$2
    local threads=$3
    local repeat=$4
    local scenario="${input_scenario}__${output_scenario}"
    local first="$RESULT_ROOT/runs/$scenario-t$threads-r1.tsv"
    local record="$RESULT_ROOT/gates/$scenario-t$threads-r$repeat.skipped.tsv"
    local first_wall
    first_wall=$(awk -F '\t' 'NR == 2 {print $6}' "$first")
    {
        printf 'scenario\tthreads\trepeat_1_application_wall_s\taction\treason\n'
        printf '%s\t%s\t%s\tSKIP_REPEAT_%s\tfirst repeat exceeds 1800 seconds; single-repeat policy\n' \
            "$scenario" "$threads" "$first_wall" "$repeat"
    } >"$record"
    printf 'SKIP %s-t%s-r%s first=%ss exceeds %ss single-repeat limit\n' \
        "$scenario" "$threads" "$repeat" "$first_wall" \
        "$SINGLE_RUN_SECONDS"
}

third_repeat_not_needed() {
    local input_scenario=$1
    local output_scenario=$2
    local threads=$3
    local scenario="${input_scenario}__${output_scenario}"
    local first="$RESULT_ROOT/runs/$scenario-t$threads-r1.tsv"
    local second="$RESULT_ROOT/runs/$scenario-t$threads-r2.tsv"
    completed_run "$first" && completed_run "$second" || return 1

    local first_wall second_wall
    first_wall=$(awk -F '\t' 'NR == 2 {print $6}' "$first")
    second_wall=$(awk -F '\t' 'NR == 2 {print $6}' "$second")
    awk -v a="$first_wall" -v b="$second_wall" \
        -v long="$LONG_RUN_SECONDS" -v limit="$LONG_RUN_VARIATION_PCT" '
        BEGIN {
            mean = (a + b) / 2
            delta = mean > 0 ? 100 * (a > b ? a - b : b - a) / mean : 0
            exit !((a > long || b > long) && delta < limit)
        }'
}

record_third_repeat_skip() {
    local input_scenario=$1
    local output_scenario=$2
    local threads=$3
    local scenario="${input_scenario}__${output_scenario}"
    local first="$RESULT_ROOT/runs/$scenario-t$threads-r1.tsv"
    local second="$RESULT_ROOT/runs/$scenario-t$threads-r2.tsv"
    local record="$RESULT_ROOT/gates/$scenario-t$threads-r3.skipped.tsv"
    local first_wall second_wall delta
    first_wall=$(awk -F '\t' 'NR == 2 {print $6}' "$first")
    second_wall=$(awk -F '\t' 'NR == 2 {print $6}' "$second")
    delta=$(awk -v a="$first_wall" -v b="$second_wall" '
        BEGIN {
            mean = (a + b) / 2
            printf "%.4f", mean > 0 ? 100 * (a > b ? a - b : b - a) / mean : 0
        }')
    {
        printf 'scenario\tthreads\trepeat_1_application_wall_s\trepeat_2_application_wall_s\tsymmetric_difference_pct\taction\treason\n'
        printf '%s\t%s\t%s\t%s\t%s\tSKIP_REPEAT_3\tlonger than 600 seconds and first two repeats differ by less than 10 percent\n' \
            "$scenario" "$threads" "$first_wall" "$second_wall" "$delta"
    } >"$record"
    printf 'SKIP %s-t%s-r3 first=%ss second=%ss difference=%s%%\n' \
        "$scenario" "$threads" "$first_wall" "$second_wall" "$delta"
}

run_case() {
    local input_scenario=$1
    local output_scenario=$2
    local threads=$3
    local repeat=$4
    output_parameters "$output_scenario"
    input_parameters "$input_scenario"
    local scenario="${input_scenario}__${output_scenario}"
    local run_id="$scenario-t$threads-r$repeat"
    local metadata="$RESULT_ROOT/runs/$run_id.tsv"
    if completed_run "$metadata"; then
        printf 'SKIP completed %s\n' "$run_id"
        return
    fi
    local run_dir
    run_dir=$(mktemp -d "$output_root/$run_id.XXXXXX")
    local input=$base_input
    if [[ "$unique_auto_index" == yes ]]; then
        input="$run_dir/input.vcf.gz"
        ln -- "$base_input" "$input"
        [[ ! -e "$input.tbi" && ! -e "$input.csi" ]] ||
            fail "automatic-CSI run unexpectedly inherited an index: $input"
    fi
    local prefix="$run_dir/output"
    local output="$prefix$output_suffix"
    local timing="$RESULT_ROOT/runs/$run_id.time.txt"
    local log="$RESULT_ROOT/logs/$run_id.log"
    local stderr="$RESULT_ROOT/logs/$run_id.stderr.txt"
    local stdout="$RESULT_ROOT/logs/$run_id.stdout.txt"
    local status=0

    # Complete outstanding writes on the target filesystem before measuring
    # this run. This wait is deliberately outside the recorded interval.
    sync -f "$output_root"
    local started app_finished durable_finished
    started=$(date +%s%N)
    set +e
    /usr/bin/time -f $'%e\t%U\t%S\t%P\t%M' -o "$timing" \
        "$NG" "--$input_kind" "$input" --threads "$threads" \
        "${input_extra[@]}" "${filter_args[@]}" "$recode_option" \
        --log-file "$log" --out "$prefix" \
        >"$stdout" 2>"$stderr"
    status=$?
    set -e
    app_finished=$(date +%s%N)
    if ((status == 0)) && [[ -s "$output" ]]; then
        sync "$output"
    fi
    durable_finished=$(date +%s%N)

    local app_wall durable_wall flush_wall user_time system_time cpu rss
    app_wall=$(seconds_between "$started" "$app_finished")
    durable_wall=$(seconds_between "$started" "$durable_finished")
    flush_wall=$(seconds_between "$app_finished" "$durable_finished")
    user_time=NA system_time=NA cpu=NA rss=NA
    if [[ -s "$timing" ]]; then
        IFS=$'\t' read -r _ user_time system_time cpu rss <"$timing"
        cpu=${cpu%%%}
    fi

    local output_bytes=0 raw_hash=missing content_hash=missing
    local exact=FAIL deterministic=NA validation=not-run
    local backend=missing kept=missing total=missing
    if [[ -s "$log" ]]; then
        backend=$(sed -n 's/^Input backend: \([^ ]*\).*/\1/p' "$log" | head -1)
        read -r kept total < <(
            sed -n 's/^After filtering, kept \([0-9]*\) out of \([0-9]*\) sites$/\1 \2/p' "$log" |
                tail -1
        ) || true
    fi
    if ((status == 0)) && [[ -s "$output" ]]; then
        output_bytes=$(stat -c '%s' "$output")
        raw_hash=$(sha256sum "$output" | awk '{print $1}')
        if [[ "$output_format" == VCF ]]; then
            content_hash=$raw_hash
            validation=plain-raw-sha256
        else
            local gate_hash_file="$RESULT_ROOT/gates/$scenario-t$threads.raw.sha256"
            local canonical_raw="$RESULT_ROOT/golden/v0130-$oracle_key-bgzf.raw.sha256"
            if ((repeat > 1)) && [[ -s "$gate_hash_file" ]] &&
               [[ "$raw_hash" == "$(cut -f1 "$gate_hash_file")" ]]; then
                content_hash=${oracle_sha256[$oracle_key]}
                validation=raw-matches-decompressed-gate
            elif [[ -s "$canonical_raw" ]] &&
                 [[ "$raw_hash" == "$(cut -f1 "$canonical_raw")" ]]; then
                content_hash=${oracle_sha256[$oracle_key]}
                validation=raw-matches-validated-canonical
            else
                gzip -t "$output"
                content_hash=$(gzip -dc "$output" | sha256sum | awk '{print $1}')
                validation=full-decompressed-sha256
            fi
        fi
        if [[ "$content_hash" == "${oracle_sha256[$oracle_key]}" &&
              "$kept" == "$expected_kept_sites" &&
              "$total" == "$expected_total_sites" ]]; then
            exact=PASS
        fi
        if ((repeat == 1)) && [[ "$exact" == PASS ]]; then
            printf '%s\t%s\n' "$raw_hash" "$run_id" \
                >"$RESULT_ROOT/gates/$scenario-t$threads.raw.sha256"
        fi
        if [[ "$output_format" == 'BGZF VCF' ]]; then
            local canonical_raw="$RESULT_ROOT/golden/v0130-$oracle_key-bgzf.raw.sha256"
            if [[ ! -s "$canonical_raw" && "$output_scenario" == ssd_bgzf &&
                  "$threads" == 1 && "$repeat" == 1 && "$exact" == PASS ]]; then
                printf '%s\t%s\n' "$raw_hash" "$run_id" >"$canonical_raw"
            fi
            if [[ -s "$canonical_raw" &&
                  "$raw_hash" == "$(cut -f1 "$canonical_raw")" ]]; then
                deterministic=YES
            else
                deterministic=NO
            fi
        else
            deterministic=YES
        fi
    fi

    printf 'scenario\tstorage\tfilesystem\tthreads\trepeat\tapplication_wall_s\tdurable_wall_s\tflush_s\tuser_s\tsystem_s\tcpu_pct\tmax_rss_kb\toutput_bytes\traw_sha256\tcontent_sha256\texact\traw_deterministic\tbackend\tkept_sites\ttotal_sites\texit\tvalidation\tinput\toutput_format\n' \
        >"$metadata"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$scenario" "$storage" "$filesystem" "$threads" "$repeat" \
        "$app_wall" "$durable_wall" "$flush_wall" "$user_time" \
        "$system_time" "$cpu" "$rss" "$output_bytes" "$raw_hash" \
        "$content_hash" "$exact" "$deterministic" "$backend" "$kept" \
        "$total" "$status" "$validation" "$input" "$output_format" \
        >>"$metadata"

    if ((status != 0)) || [[ "$exact" != PASS ]]; then
        printf 'FAILED %s exit=%s exact=%s retained=%s\n' \
            "$run_id" "$status" "$exact" "$run_dir" >&2
        exit 1
    fi

    # Preserve one complete compressed v0.13.0 artifact and its hashes. All
    # other large candidate outputs are summarized and then removed.
    if [[ "$output_scenario" == ssd_bgzf && "$threads" == 1 && "$repeat" == 1 ]]; then
        local retained="$RESULT_ROOT/golden/v0130.$input_scenario.ssd-bgzf-t1.recode.vcf.gz"
        if [[ ! -e "$retained" ]]; then
            mv -- "$output" "$retained"
            sync "$retained"
        else
            [[ "$(sha256sum "$retained" | awk '{print $1}')" == "$raw_hash" ]] ||
                fail "retained BGZF golden differs: $retained"
        fi
        {
            printf 'raw_sha256\t%s\n' "$raw_hash"
            printf 'decompressed_sha256\t%s\n' "$content_hash"
            printf 'bytes\t%s\n' "$output_bytes"
            printf 'source_run\t%s\n' "$run_id"
        } >"$RESULT_ROOT/golden/v0130.$input_scenario.ssd-bgzf-t1.manifest.tsv"
    fi
    rm -f -- "$output"
    if [[ "$unique_auto_index" == yes ]]; then
        rm -f -- "$input.csi" "$input.tbi" "$input"
    fi
    rmdir "$run_dir"
    printf 'DONE %s app=%ss durable=%ss exact=PASS backend=%s\n' \
        "$run_id" "$app_wall" "$durable_wall" "$backend"
}

refresh_reports() {
    local all="$RESULT_ROOT/all-runs.tsv"
    local first
    first=$(find "$RESULT_ROOT/runs" -maxdepth 1 -type f -name '*.tsv' |
        sort | head -1)
    [[ -n "$first" ]] || return
    head -1 "$first" >"$all"
    find "$RESULT_ROOT/runs" -maxdepth 1 -type f -name '*.tsv' -print0 |
        sort -z | xargs -0 -r -n1 tail -1 >>"$all"

    local summary="$RESULT_ROOT/summary.tsv"
    printf 'input_scenario\toutput_scenario\tscenario\tstorage\tfilesystem\tthreads\truns\tapplication_mean_s\tapplication_median_s\tdurable_mean_s\tdurable_median_s\tflush_mean_s\tmean_cpu_pct\tmax_rss_kb\toutput_bytes\tapplication_speedup\tdurable_speedup\texact\traw_deterministic\tbackend\n' \
        >"$summary"
    declare -A app_baseline=() durable_baseline=()
    local input_scenario output_scenario scenario threads
    for input_scenario in "${input_scenarios[@]}"; do
      for output_scenario in "${output_scenarios[@]}"; do
        scenario="${input_scenario}__${output_scenario}"
        for threads in $THREAD_LIST; do
            mapfile -t app_values < <(
                awk -F '\t' -v s="$scenario" -v t="$threads" \
                    'NR > 1 && $1 == s && $4 == t {print $6}' "$all" |
                    sort -n
            )
            ((${#app_values[@]} > 0)) || continue
            mapfile -t durable_values < <(
                awk -F '\t' -v s="$scenario" -v t="$threads" \
                    'NR > 1 && $1 == s && $4 == t {print $7}' "$all" |
                    sort -n
            )
            local n=${#app_values[@]}
            local middle=$((n / 2))
            local app_median durable_median
            if ((n % 2 == 1)); then
                app_median=${app_values[$middle]}
                durable_median=${durable_values[$middle]}
            else
                app_median=$(awk -v a="${app_values[$((middle - 1))]}" \
                    -v b="${app_values[$middle]}" 'BEGIN {printf "%.6f", (a+b)/2}')
                durable_median=$(awk -v a="${durable_values[$((middle - 1))]}" \
                    -v b="${durable_values[$middle]}" 'BEGIN {printf "%.6f", (a+b)/2}')
            fi
            local row
            row=$(awk -F '\t' -v s="$scenario" -v t="$threads" '
                NR > 1 && $1 == s && $4 == t {
                    n++; app += $6; durable += $7; flush += $8; cpu += $11
                    if ($12 > rss) rss = $12
                    bytes = $13; storage = $2; fs = $3; backend[$18] = 1
                    if ($16 != "PASS") exact = "FAIL"
                    if ($17 != "YES") deterministic = "NO"
                }
                END {
                    if (exact == "") exact = "PASS"
                    if (deterministic == "") deterministic = "YES"
                    backend_list = ""
                    for (value in backend) {
                        if (backend_list != "") backend_list = backend_list ","
                        backend_list = backend_list value
                    }
                    printf "%s\t%s\t%.6f\t%.6f\t%.6f\t%.1f\t%d\t%s\t%s\t%s\t%s\n", \
                        storage, fs, app/n, durable/n, flush/n, cpu/n, rss, \
                        bytes, exact, deterministic, backend_list
                }
            ' "$all")
            local storage_value fs_value app_mean durable_mean flush_mean
            local cpu_mean rss_max bytes exact_value deterministic_value backend_value
            IFS=$'\t' read -r storage_value fs_value app_mean durable_mean \
                flush_mean cpu_mean rss_max bytes exact_value \
                deterministic_value backend_value <<<"$row"
            if [[ "$threads" == 1 ]]; then
                app_baseline[$scenario]=$app_mean
                durable_baseline[$scenario]=$durable_mean
            fi
            local app_speedup durable_speedup
            app_speedup=$(awk -v base="${app_baseline[$scenario]}" -v x="$app_mean" \
                'BEGIN {printf "%.4f", base/x}')
            durable_speedup=$(awk -v base="${durable_baseline[$scenario]}" -v x="$durable_mean" \
                'BEGIN {printf "%.4f", base/x}')
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$input_scenario" "$output_scenario" "$scenario" \
                "$storage_value" "$fs_value" "$threads" "$n" \
                "$app_mean" "$app_median" "$durable_mean" "$durable_median" \
                "$flush_mean" "$cpu_mean" "$rss_max" "$bytes" \
                "$app_speedup" "$durable_speedup" "$exact_value" \
                "$deterministic_value" "$backend_value" >>"$summary"
        done
      done
    done
}

run_smoke
write_manifest
if [[ "$SMOKE_ONLY" == 1 ]]; then
    printf 'V0.13.0 OUTPUT/STORAGE SMOKE-ONLY PASS\n'
    exit 0
fi

for input_scenario in "${input_scenarios[@]}"; do
  for output_scenario in "${output_scenarios[@]}"; do
    scenario="${input_scenario}__${output_scenario}"
    gate="$RESULT_ROOT/gates/$scenario.pass"
    if [[ -s "$gate" ]]; then
        continue
    fi
    for threads in $THREAD_LIST; do
        run_case "$input_scenario" "$output_scenario" "$threads" 1
        refresh_reports
    done
    printf 'scenario\tstatus\tvalidated\n%s\tPASS\t%s\n' \
        "$scenario" "$(date --iso-8601=seconds)" >"$gate"
    printf 'FIRST-REPEAT GATE PASS %s\n' "$scenario"
  done
done

refresh_reports
if ! grep -q $'^first_repeat_gates_passed\t' "$RESULT_ROOT/manifest.tsv"; then
    printf 'first_repeat_gates_passed\t%s\n' \
        "$(date --iso-8601=seconds)" >>"$RESULT_ROOT/manifest.tsv"
fi
if [[ "$GATE_ONLY" == 1 ]]; then
    printf 'V0.13.0 OUTPUT/STORAGE FIRST-REPEAT GATES PASS\n'
    exit 0
fi

for repeat in $(seq 2 "$REPEATS"); do
  for input_scenario in "${input_scenarios[@]}"; do
    for output_scenario in "${output_scenarios[@]}"; do
        for threads in $THREAD_LIST; do
            scenario="${input_scenario}__${output_scenario}"
            current="$RESULT_ROOT/runs/$scenario-t$threads-r$repeat.tsv"
            if ! completed_run "$current" &&
               first_run_is_single_only \
                   "$input_scenario" "$output_scenario" "$threads"; then
                record_single_run_skip \
                    "$input_scenario" "$output_scenario" "$threads" "$repeat"
                continue
            fi
            if ((repeat == 3)); then
                third="$RESULT_ROOT/runs/$scenario-t$threads-r3.tsv"
                if ! completed_run "$third" &&
                   third_repeat_not_needed \
                       "$input_scenario" "$output_scenario" "$threads"; then
                    record_third_repeat_skip \
                        "$input_scenario" "$output_scenario" "$threads"
                    continue
                fi
            fi
            run_case "$input_scenario" "$output_scenario" \
                "$threads" "$repeat"
            refresh_reports
        done
    done
  done
done

refresh_reports
if ! grep -q $'^finished\t' "$RESULT_ROOT/manifest.tsv"; then
    printf 'finished\t%s\n' "$(date --iso-8601=seconds)" \
        >>"$RESULT_ROOT/manifest.tsv"
fi
printf 'FINAL v0.13.0 FOUR-INPUT / THREE-OUTPUT SAME-DEVICE MATRIX PASS %s\n' \
    "$RESULT_ROOT"
