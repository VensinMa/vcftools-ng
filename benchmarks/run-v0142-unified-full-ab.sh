#!/usr/bin/env bash
set -euo pipefail

source_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
result_root=${RESULT_ROOT:-"$source_root/benchmarks/results/full-unified-v0142-ab"}
scratch_root=${SCRATCH_ROOT:-/home/vensin/workspace/vcftools-ng-benchmark-v0142/ssd}
binary_root=${BINARY_ROOT:-/home/vensin/workspace/vcftools-ng-benchmark-v0142/binaries}
v13=${V13:-/home/vensin/software/vcftools-ng-v0.13.0-linux-x86_64/bin/vcftools-ng}
bgzf=${FULL_BGZF:-/home/vensin/workspace/Sweet_Osmanthus/05.variant_filter/02.vcftools_filter_snp_indel/412samples.SNP.biallelic.minGQ10.minQ30.meanDP6.maxmiss0.8.maf0.05.vcf.gz}
plain=${FULL_PLAIN:-"$source_root/data/osmanthus412.snps.vcf"}
bcf=${FULL_BCF:-"$source_root/data/osmanthus412.snps.bcf"}
thread_list=${THREAD_LIST:-"1 2 4 8 12 16 24 28 32"}

readonly filter_args=(
    --min-alleles 2 --max-alleles 2 --minGQ 10 --minQ 30
    --min-meanDP 7 --max-missing 0.9 --maf 0.1
)

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

mkdir -p "$result_root/runs" "$result_root/logs" "$scratch_root" \
    "$binary_root/v0141" "$binary_root/v0142"

if [[ ! -x "$binary_root/v0141/vcftools-ng-v0.14.1-linux-x86_64/bin/vcftools-ng" ]]; then
    tar -xzf "$source_root/dist/vcftools-ng-v0.14.1-linux-x86_64.tar.gz" \
        -C "$binary_root/v0141"
fi
if [[ ! -x "$binary_root/v0142/vcftools-ng-v0.14.2-linux-x86_64/bin/vcftools-ng" ]]; then
    tar -xzf "$source_root/dist/vcftools-ng-v0.14.2-linux-x86_64.tar.gz" \
        -C "$binary_root/v0142"
fi
v14="$binary_root/v0141/vcftools-ng-v0.14.1-linux-x86_64/bin/vcftools-ng"
v142="$binary_root/v0142/vcftools-ng-v0.14.2-linux-x86_64/bin/vcftools-ng"

for executable in "$v13" "$v14" "$v142"; do
    [[ -x "$executable" ]] || fail "missing executable: $executable"
done
for input in "$bgzf" "$bgzf.tbi" "$plain" "$bcf"; do
    [[ -s "$input" ]] || fail "missing input or index: $input"
done
[[ $($v13 --version) == "vcftools-ng 0.13.0" ]] || fail "invalid v0.13.0 binary"
[[ $($v14 --version) == "vcftools-ng 0.14.1" ]] || fail "invalid v0.14.1 binary"
[[ $($v142 --version) == "vcftools-ng 0.14.2" ]] || fail "invalid v0.14.2 binary"

summary="$result_root/all-runs.tsv"
if [[ ! -e "$summary" ]]; then
    printf 'scenario\tthreads\trepeat\tversion\twall_s\tdurable_s\tflush_s\tcpu_pct\tmax_rss_kb\toutput_bytes\tsha256\texact\tkept\ttotal\tbackend\n' \
        >"$summary"
fi
cat >"$result_root/manifest.tsv" <<EOF
key\tvalue
workload\t${filter_args[*]} --recode-vcf --recode-INFO-all
versions\tv0.13.0;v0.14.1;v0.14.2 portable Linux x86_64
scenarios\tBGZF VCF + TBI;Plain VCF;BCF adaptive stream
threads\t$thread_list
output\tPlain VCF on /dev/nvme3n1p2
repeat_policy\tone strictly serial run per row; differences within 5% are treated as effectively tied
execution\tstrictly serial, rotating version order
EOF

run_one() {
    local scenario=$1 threads=$2 version=$3 binary=$4
    local flag input expected_hash expected_bytes
    case "$scenario" in
        bgzf_tbi)
            flag=--gzvcf
            input=$bgzf
            expected_hash=7548416e01d4a318b81c5d1feb9429f60c7995205d66169242c3792af4c4fc14
            expected_bytes=59434159204
            ;;
        plain_vcf)
            flag=--vcf
            input=$plain
            expected_hash=d4f2a15e8c5ad0cc12abf4a3ab308bb48f22adf8ec13776d72ceeaa1f8d402b8
            expected_bytes=59434159621
            ;;
        bcf)
            flag=--bcf
            input=$bcf
            expected_hash=dde9edd98d5d05aa885e0e2f78a9696cfe0802c472d4dfd81ccffb49339107f3
            expected_bytes=57211771106
            ;;
        *) fail "unknown scenario: $scenario" ;;
    esac

    local id="${scenario}-t${threads}-${version}-r1"
    local metadata="$result_root/runs/$id.tsv"
    if [[ -s "$metadata" ]] &&
       awk -F '\t' 'NR == 2 && $12 == "PASS" {ok=1} END {exit !ok}' \
           "$metadata"; then
        printf 'SKIP\t%s\n' "$id"
        return
    fi

    local scratch="$scratch_root/$id"
    local prefix="$scratch/output"
    local output="$prefix.recode.vcf"
    local timing="$result_root/runs/$id.time.txt"
    local stderr_log="$result_root/logs/$id.stderr.txt"
    rm -rf -- "$scratch"
    mkdir -p "$scratch"

    local start_ns application_end_ns durable_end_ns
    local wall durable flush cpu rss bytes digest backend status=0 exact=FAIL
    local kept=NA total=NA
    start_ns=$(date +%s%N)
    /usr/bin/time -f '%e\t%U\t%S\t%P\t%M' -o "$timing" \
        "$binary" "$flag" "$input" --threads "$threads" \
        "${filter_args[@]}" --recode-vcf --recode-INFO-all \
        --no-log-file --out "$prefix" \
        >/dev/null 2>"$stderr_log" || status=$?
    application_end_ns=$(date +%s%N)
    if ((status == 0)) && [[ -s "$output" ]]; then
        sync -f "$output"
    fi
    durable_end_ns=$(date +%s%N)

    IFS=$'\t' read -r wall _ _ cpu rss <"$timing"
    cpu=${cpu%%%}
    durable=$(awk -v a="$start_ns" -v b="$durable_end_ns" \
        'BEGIN {printf "%.3f", (b-a)/1e9}')
    flush=$(awk -v a="$application_end_ns" -v b="$durable_end_ns" \
        'BEGIN {printf "%.3f", (b-a)/1e9}')
    bytes=$(stat -c '%s' "$output" 2>/dev/null || printf 0)
    digest=$(sha256sum "$output" 2>/dev/null | awk '{print $1}' || printf missing)
    read -r kept total < <(
        sed -n 's/^After filtering, kept \([0-9][0-9]*\) out of \([0-9][0-9]*\) sites$/\1 \2/p' \
            "$stderr_log" | tail -1
    ) || true
    kept=${kept:-NA}
    total=${total:-NA}
    backend=$(sed -n 's/^Input backend: //p' "$stderr_log" | head -1)
    backend=${backend:-missing}
    if ((status == 0)) && [[ "$bytes" == "$expected_bytes" ]] &&
       [[ "$digest" == "$expected_hash" ]] &&
       [[ "$kept" == 5425725 && "$total" == 11230392 ]]; then
        exact=PASS
    fi

    {
        printf 'scenario\tthreads\trepeat\tversion\twall_s\tdurable_s\tflush_s\tcpu_pct\tmax_rss_kb\toutput_bytes\tsha256\texact\tkept\ttotal\tbackend\n'
        printf '%s\t%s\t1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$scenario" "$threads" "$version" "$wall" "$durable" \
            "$flush" "$cpu" "$rss" "$bytes" "$digest" "$exact" \
            "$kept" "$total" "$backend"
    } >"$metadata.tmp"
    mv -f -- "$metadata.tmp" "$metadata"
    sed -n '2p' "$metadata" >>"$summary"
    printf 'DONE\t%s\twall=%s\tdurable=%s\tcpu=%s\trss=%s\texact=%s\n' \
        "$id" "$wall" "$durable" "$cpu" "$rss" "$exact"
    if [[ "$exact" != PASS ]]; then
        printf 'FAILED output retained: %s\n' "$output" >&2
        exit 1
    fi
    rm -f -- "$output"
    rmdir "$scratch"
}

for scenario in bgzf_tbi plain_vcf bcf; do
    index=0
    for threads in $thread_list; do
        case $((index % 3)) in
            0) order=(v0.13.0 v0.14.1 v0.14.2) ;;
            1) order=(v0.14.2 v0.14.1 v0.13.0) ;;
            2) order=(v0.14.1 v0.13.0 v0.14.2) ;;
        esac
        for version in "${order[@]}"; do
            case "$version" in
                v0.13.0) binary=$v13 ;;
                v0.14.1) binary=$v14 ;;
                v0.14.2) binary=$v142 ;;
            esac
            run_one "$scenario" "$threads" "$version" "$binary"
        done
        index=$((index + 1))
    done
done

printf 'FULL_UNIFIED_V0142_FIRST_REPEAT_PASS\n'
