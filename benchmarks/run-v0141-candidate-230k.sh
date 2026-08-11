#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 /path/to/vcftools-ng RESULT_DIR" >&2
    exit 2
fi

ng=$1
result=$2
project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
root=$project_root/benchmarks/results/workload-matrix-230k-v0130
input=$root/osmanthus412.flags.23chr_10k.vcf
fixtures=$root/fixtures
oracles=$root/oracles
repeats=${REPEATS:-3}
read -r -a threads <<<"${THREADS:-1 8 16 32}"
read -r -a cases <<<"${CASES:-W01_counts W02_ten_filter_counts W06_keep_50pct_filtered_recode W09_site_fst_small_pair W09_site_fst_large_pair W10_window_fst_biallelic W13_genotype_ld W14_pca_64_samples}"

expected_input_sha=e62f8f06617371229825bb4747777469ac173c56aa4bc03cec469b29bbee1be2
expected_w01_sha=d570b1be69cfded0946c95fa2375230c9859c5be3f2d1b55a7d6e520ce294815
expected_w02_sha=3699c4c996fbfca6c78d8c4865283b890a3a44f728a7b427143c24cb8cb49992

[[ -x $ng ]] || { echo "vcftools-ng missing: $ng" >&2; exit 2; }
[[ -f $input ]] || { echo "230k input missing: $input" >&2; exit 2; }
[[ -d $fixtures && -f $oracles/SHA256SUMS ]] || {
    echo "Locked 230k fixtures/oracles are missing" >&2
    exit 2
}
[[ ! -e $result ]] || {
    echo "Refusing to overwrite benchmark result: $result" >&2
    exit 2
}

actual_input_sha=$(sha256sum -- "$input" | awk '{print $1}')
[[ $actual_input_sha == "$expected_input_sha" ]] || {
    echo "230k input SHA-256 mismatch: $actual_input_sha" >&2
    exit 1
}
(
    cd "$oracles"
    sha256sum --check --strict SHA256SUMS
)

mkdir -p "$result"/{logs,outputs,times,reference}
sha256sum -- "$ng" >"$result/binary.sha256"
printf '%s\t%s\n' input_sha256 "$actual_input_sha" >"$result/manifest.tsv"
printf '%s\t%s\n' repeats "$repeats" >>"$result/manifest.tsv"
printf '%s\t%s\n' threads "${threads[*]}" >>"$result/manifest.tsv"
printf '%s\t%s\n' cases "${cases[*]}" >>"$result/manifest.tsv"

runs=$result/runs.tsv
printf '%s\n' $'case\tthreads\trepeat\twall_s\tuser_s\tsystem_s\tcpu\tpeak_rss_kib\tbackend\tgate' >"$runs"

configure() {
    case $1 in
        W01_counts)
            args=(--counts)
            suffix=.frq.count
            gate_kind=sha256
            expected_sha=$expected_w01_sha ;;
        W02_ten_filter_counts)
            args=(--min-alleles 2 --max-alleles 2
                  --minQ 30 --minGQ 10 --minDP 5 --maxDP 30
                  --min-meanDP 7 --max-missing 0.9 --maf 0.1 --mac 2
                  --counts)
            suffix=.frq.count
            gate_kind=sha256
            expected_sha=$expected_w02_sha ;;
        W06_keep_50pct_filtered_recode)
            args=(--keep "$fixtures/keep-50pct.samples.txt"
                  --min-alleles 2 --max-alleles 2
                  --minGQ 10 --minQ 30 --min-meanDP 7
                  --max-missing 0.9 --maf 0.1
                  --recode-vcf --recode-INFO-all)
            suffix=.recode.vcf
            gate_kind=candidate-reference ;;
        W09_site_fst_small_pair)
            args=(--weir-fst-pop "$fixtures/pop-small-ancient-12.txt"
                  --weir-fst-pop "$fixtures/pop-small-asiaticus-13.txt")
            suffix=.weir.fst
            gate_kind=oracle
            oracle=$oracles/W09_site_fst_small_pair/oracle.weir.fst ;;
        W09_site_fst_large_pair)
            args=(--weir-fst-pop "$fixtures/pop-large-wild-166.txt"
                  --weir-fst-pop "$fixtures/pop-large-cultivated-243.txt")
            suffix=.weir.fst
            gate_kind=oracle
            oracle=$oracles/W09_site_fst_large_pair/oracle.weir.fst ;;
        W10_window_fst_biallelic)
            args=(--min-alleles 2 --max-alleles 2
                  --weir-fst-pop "$fixtures/pop-large-wild-166.txt"
                  --weir-fst-pop "$fixtures/pop-large-cultivated-243.txt"
                  --fst-window-size 100000 --fst-window-step 10000)
            suffix=.windowed.weir.fst
            gate_kind=oracle
            oracle=$oracles/W10_window_fst_biallelic/oracle.windowed.weir.fst ;;
        W13_genotype_ld)
            args=(--geno-r2 --ld-window 200
                  --ld-window-bp 1000000000 --min-r2 0.99)
            suffix=.geno.ld
            gate_kind=sha256
            expected_sha=d15d43a67e2365bdb1d308cdb61050a93ab71f7b735ca54d749ed469f95de7c5 ;;
        W14_pca_64_samples)
            args=(--keep "$project_root/benchmarks/fixtures/v0132-pca-keep-64.txt"
                  --max-missing 1 --pca)
            suffix=.pca
            gate_kind=sha256
            expected_sha=7c615aad3c47a59ed3e672299ec42ae0bd85482a4e95ce4112ccf4212e594794 ;;
        *)
            echo "Unknown v0.14.1 candidate workload: $1" >&2
            exit 2 ;;
    esac
}

# Use one warm-page-cache regime for this small candidate A/B. Release
# benchmarks record cold/warm state separately.
dd if="$input" of=/dev/null bs=32M status=none

for case_id in "${cases[@]}"; do
    unset expected_sha oracle
    configure "$case_id"
    for thread in "${threads[@]}"; do
        for repeat in $(seq 1 "$repeats"); do
            stem=$case_id-t$thread-r$repeat
            prefix=$result/outputs/$stem
            log=$result/logs/$stem.log
            time_file=$result/times/$stem.tsv
            /usr/bin/time -f $'%e\t%U\t%S\t%P\t%M' -o "$time_file" \
                "$ng" --vcf "$input" --threads "$thread" \
                "${args[@]}" --log-file "$log" --out "$prefix" \
                >/dev/null 2>"$result/logs/$stem.stderr"
            artifact=$prefix$suffix
            case $gate_kind in
                sha256)
                    actual_sha=$(sha256sum -- "$artifact" | awk '{print $1}')
                    [[ $actual_sha == "$expected_sha" ]] || {
                        echo "$stem output SHA-256 mismatch: $actual_sha" >&2
                        exit 1
                    } ;;
                oracle)
                    cmp "$oracle" "$artifact" ;;
                candidate-reference)
                    reference=$result/reference/$case_id$suffix
                    if [[ ! -e $reference ]]; then
                        [[ $thread == "${threads[0]}" && $repeat == 1 ]] || {
                            echo "Candidate recode reference was not created by the first run" >&2
                            exit 1
                        }
                        mv -- "$artifact" "$reference"
                        artifact=$reference
                    else
                        cmp "$reference" "$artifact"
                    fi ;;
            esac
            read -r wall user system cpu rss <"$time_file"
            backend=$(sed -n 's/^Input backend: //p' "$log" | head -1)
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
                "$case_id" "$thread" "$repeat" "$wall" "$user" \
                "$system" "$cpu" "$rss" "${backend:-unknown}" >>"$runs"
            [[ $artifact == "$result/reference/"* ]] || rm -f -- "$artifact"
        done
    done
done

{
    printf '%s\n' $'case\tthreads\tmean_wall_s\tmean_rss_kib'
    awk -F'\t' '
        BEGIN {OFS="\t"}
        NR > 1 {key=$1 FS $2; wall[key]+=$4; rss[key]+=$8; count[key]++}
        END {for (key in count) print key, wall[key]/count[key], rss[key]/count[key]}
    ' "$runs" | sort -t$'\t' -k1,1 -k2,2n
} >"$result/means.tsv"

echo "v0.14.1 candidate 230k baseline complete: $runs"
