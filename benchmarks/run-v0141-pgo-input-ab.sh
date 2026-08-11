#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 CONTROL_BINARY PGO_BINARY RESULT_DIR" >&2
    exit 2
fi

control=$1
pgo=$2
result=$3
project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
root=$project_root/benchmarks/results/workload-matrix-230k-v0130
indexed=$root/osmanthus412.flags.23chr_10k.vcf.gz
noindex=$root/osmanthus412.flags.23chr_10k.noindex.vcf.gz
expected_input_sha=8b482fd4560ddc0a057629553638f17a4fde9eebf17bd8535dd9ccf6310af55a
repeats=${REPEATS:-3}
read -r -a threads <<<"${THREADS:-1 8 16 32}"
read -r -a input_cases <<<"${INPUT_CASES:-indexed noindex}"
read -r -a workloads <<<"${WORKLOADS:-W01_counts W02_ten_filter_counts}"

[[ -x $control && -x $pgo ]] || {
    echo "Both A/B executables must exist and be executable" >&2
    exit 2
}
[[ -f $indexed && -f $indexed.tbi && -f $noindex ]] || {
    echo "Locked 230k BGZF inputs or TBI are missing" >&2
    exit 2
}
for path in "$indexed" "$noindex"; do
    actual=$(sha256sum -- "$path" | awk '{print $1}')
    [[ $actual == "$expected_input_sha" ]] || {
        echo "BGZF input SHA-256 mismatch: $path $actual" >&2
        exit 1
    }
done
if [[ -e $noindex.tbi || -e $noindex.csi ]]; then
    echo "No-index BGZF fixture unexpectedly has a sidecar" >&2
    exit 1
fi
[[ ! -e $result ]] || {
    echo "Refusing to overwrite PGO input A/B result: $result" >&2
    exit 2
}

mkdir -p "$result"/{outputs,times}
sha256sum -- "$control" "$pgo" >"$result/binaries.sha256"
runs=$result/runs.tsv
printf '%s\n' $'input\tcase\tversion\tthreads\trepeat\twall_s\tuser_s\tsystem_s\tpeak_rss_kib\tbackend\tgate' >"$runs"

configure_workload() {
    case $1 in
        W01_counts)
            args=(--counts)
            suffix=.frq.count
            expected=d570b1be69cfded0946c95fa2375230c9859c5be3f2d1b55a7d6e520ce294815 ;;
        W02_ten_filter_counts)
            args=(--min-alleles 2 --max-alleles 2
                  --minQ 30 --minGQ 10 --minDP 5 --maxDP 30
                  --min-meanDP 7 --max-missing 0.9 --maf 0.1 --mac 2
                  --counts)
            suffix=.frq.count
            expected=3699c4c996fbfca6c78d8c4865283b890a3a44f728a7b427143c24cb8cb49992 ;;
        W06_keep_50pct_filtered_recode)
            args=(--input-backend stream
                  --keep "$root/fixtures/keep-50pct.samples.txt"
                  --min-alleles 2 --max-alleles 2
                  --minGQ 10 --minQ 30 --min-meanDP 7
                  --max-missing 0.9 --maf 0.1
                  --recode-vcf --recode-INFO-all)
            suffix=.recode.vcf
            expected=bc50b61b8fa7edfd4c8dde12061769c10c84d311416c676f8fb765a0ba4d79f7 ;;
        *) echo "Unknown PGO input workload: $1" >&2; exit 2 ;;
    esac
}

for repeat in $(seq 1 "$repeats"); do
    versions=(control pgo)
    if (( repeat % 2 == 0 )); then
        versions=(pgo control)
    fi
    for input_case in "${input_cases[@]}"; do
        input=$indexed
        [[ $input_case == noindex ]] && input=$noindex
        for case_id in "${workloads[@]}"; do
            configure_workload "$case_id"
            for thread in "${threads[@]}"; do
                for version in "${versions[@]}"; do
                    binary=$control
                    [[ $version == pgo ]] && binary=$pgo
                    stem=$input_case-$case_id-$version-t$thread-r$repeat
                    prefix=$result/outputs/$stem
                    time_file=$result/times/$stem.tsv
                    stderr=$result/$stem.stderr
                    /usr/bin/time -f $'%e\t%U\t%S\t%M' -o "$time_file" \
                        "$binary" --gzvcf "$input" --threads "$thread" \
                        "${args[@]}" --no-log-file --out "$prefix" \
                        >/dev/null 2>"$stderr"
                    artifact=$prefix$suffix
                    actual=$(sha256sum -- "$artifact" | awk '{print $1}')
                    [[ $actual == "$expected" ]] || {
                        echo "$stem output SHA-256 mismatch: $actual" >&2
                        exit 1
                    }
                    read -r wall user system rss <"$time_file"
                    backend=$(sed -n 's/^Input backend: //p' "$stderr" | head -1)
                    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
                        "$input_case" "$case_id" "$version" "$thread" \
                        "$repeat" "$wall" "$user" "$system" "$rss" \
                        "${backend:-unknown}" >>"$runs"
                    rm -f -- "$artifact"
                done
            done
        done
    done
done

{
    printf '%s\n' $'input\tcase\tversion\tthreads\tmean_wall_s\tmean_user_s'
    awk -F'\t' '
        BEGIN {OFS="\t"}
        NR > 1 {key=$1 FS $2 FS $3 FS $4; wall[key]+=$6; user[key]+=$7; count[key]++}
        END {for (key in count) print key, wall[key]/count[key], user[key]/count[key]}
    ' "$runs" | sort -t$'\t' -k1,1 -k2,2 -k4,4n -k3,3
} >"$result/means.tsv"

echo "v0.14.1 PGO BGZF input A/B complete: $runs"
