#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 CONTROL_BINARY CANDIDATE_BINARY RESULT_DIR" >&2
    exit 2
fi

control=$1
candidate=$2
result=$3
project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
root=$project_root/benchmarks/results/workload-matrix-230k-v0130
input=$root/osmanthus412.flags.23chr_10k.vcf
fixtures=$root/fixtures
oracles=$root/oracles
repeats=${REPEATS:-10}
read -r -a threads <<<"${THREADS:-8 16 32}"
read -r -a cases <<<"${CASES:-W01_counts W02_ten_filter_counts W10_window_fst_biallelic}"

[[ -x $control && -x $candidate ]] || {
    echo "Both A/B executables must exist and be executable" >&2
    exit 2
}
[[ ! -e $result ]] || {
    echo "Refusing to overwrite A/B result: $result" >&2
    exit 2
}
mkdir -p "$result"/{outputs,times}
sha256sum -- "$control" "$candidate" >"$result/binaries.sha256"
runs=$result/runs.tsv
printf '%s\n' $'case\tversion\tthreads\trepeat\twall_s\tuser_s\tsystem_s\tpeak_rss_kib\tgate' >"$runs"

configure() {
    case $1 in
        W01_counts)
            args=(--counts)
            suffix=.frq.count
            expected_sha=d570b1be69cfded0946c95fa2375230c9859c5be3f2d1b55a7d6e520ce294815
            oracle= ;;
        W02_ten_filter_counts)
            args=(--min-alleles 2 --max-alleles 2
                  --minQ 30 --minGQ 10 --minDP 5 --maxDP 30
                  --min-meanDP 7 --max-missing 0.9 --maf 0.1 --mac 2
                  --counts)
            suffix=.frq.count
            expected_sha=3699c4c996fbfca6c78d8c4865283b890a3a44f728a7b427143c24cb8cb49992
            oracle= ;;
        W06_keep_50pct_filtered_recode)
            args=(--keep "$fixtures/keep-50pct.samples.txt"
                  --min-alleles 2 --max-alleles 2
                  --minGQ 10 --minQ 30 --min-meanDP 7
                  --max-missing 0.9 --maf 0.1
                  --recode-vcf --recode-INFO-all)
            suffix=.recode.vcf
            expected_sha=bc50b61b8fa7edfd4c8dde12061769c10c84d311416c676f8fb765a0ba4d79f7
            oracle= ;;
        W09_site_fst_large_pair)
            args=(--weir-fst-pop "$fixtures/pop-large-wild-166.txt"
                  --weir-fst-pop "$fixtures/pop-large-cultivated-243.txt")
            suffix=.weir.fst
            expected_sha=
            oracle=$oracles/W09_site_fst_large_pair/oracle.weir.fst ;;
        W10_window_fst_biallelic)
            args=(--min-alleles 2 --max-alleles 2
                  --weir-fst-pop "$fixtures/pop-large-wild-166.txt"
                  --weir-fst-pop "$fixtures/pop-large-cultivated-243.txt"
                  --fst-window-size 100000 --fst-window-step 10000)
            suffix=.windowed.weir.fst
            expected_sha=
            oracle=$oracles/W10_window_fst_biallelic/oracle.windowed.weir.fst ;;
        W13_genotype_ld)
            args=(--geno-r2 --ld-window 200
                  --ld-window-bp 1000000000 --min-r2 0.99)
            suffix=.geno.ld
            expected_sha=d15d43a67e2365bdb1d308cdb61050a93ab71f7b735ca54d749ed469f95de7c5
            oracle= ;;
        W14_pca_64_samples)
            args=(--keep "$project_root/benchmarks/fixtures/v0132-pca-keep-64.txt"
                  --max-missing 1 --pca)
            suffix=.pca
            expected_sha=7c615aad3c47a59ed3e672299ec42ae0bd85482a4e95ce4112ccf4212e594794
            oracle= ;;
        *) echo "Unknown field-requirement A/B case: $1" >&2; exit 2 ;;
    esac
}

dd if="$input" of=/dev/null bs=32M status=none
for repeat in $(seq 1 "$repeats"); do
    versions=(control candidate)
    if (( repeat % 2 == 0 )); then
        versions=(candidate control)
    fi
    for case_id in "${cases[@]}"; do
        configure "$case_id"
        for thread in "${threads[@]}"; do
            for version in "${versions[@]}"; do
                binary=$control
                [[ $version == candidate ]] && binary=$candidate
                stem=$case_id-$version-t$thread-r$repeat
                prefix=$result/outputs/$stem
                time_file=$result/times/$stem.tsv
                /usr/bin/time -f $'%e\t%U\t%S\t%M' -o "$time_file" \
                    "$binary" --vcf "$input" --threads "$thread" \
                    "${args[@]}" --no-log-file --out "$prefix" \
                    >/dev/null 2>"$result/$stem.stderr"
                artifact=$prefix$suffix
                if [[ -n $oracle ]]; then
                    cmp "$oracle" "$artifact"
                else
                    actual_sha=$(sha256sum -- "$artifact" | awk '{print $1}')
                    [[ $actual_sha == "$expected_sha" ]] || {
                        echo "$stem output SHA-256 mismatch: $actual_sha" >&2
                        exit 1
                    }
                fi
                read -r wall user system rss <"$time_file"
                printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
                    "$case_id" "$version" "$thread" "$repeat" \
                    "$wall" "$user" "$system" "$rss" >>"$runs"
                rm -f -- "$artifact"
            done
        done
    done
done

{
    printf '%s\n' $'case\tversion\tthreads\tmean_wall_s\tmean_user_s'
    awk -F'\t' '
        BEGIN {OFS="\t"}
        NR > 1 {key=$1 FS $2 FS $3; wall[key]+=$5; user[key]+=$6; count[key]++}
        END {for (key in count) print key, wall[key]/count[key], user[key]/count[key]}
    ' "$runs" | sort -t$'\t' -k1,1 -k3,3n -k2,2
} >"$result/means.tsv"

echo "v0.14.1 candidate A/B complete: $runs"
