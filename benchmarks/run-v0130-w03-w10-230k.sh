#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 /path/to/vcftools-ng RESULT_DIR" >&2
    exit 2
fi

ng=$1
result=$2
root=${BENCH_ROOT:-/home/vensin/workspace/vcftools-ng/benchmarks/results/workload-matrix-230k-v0130}
input=$root/osmanthus412.flags.23chr_10k.vcf
fixtures=$root/fixtures
oracles=$root/oracles
original=${ORIGINAL:-/home/vensin/anaconda3/envs/vcftools/bin/vcftools}
repeats=${REPEATS:-3}
read -r -a threads <<<"${THREADS:-1 4 8 16 32}"

cases=(
    W03_positions_1pct_sorted
    W03_positions_1pct_shuffled
    W03_positions_50pct_sorted
    W03_positions_50pct_shuffled
    W04_exclude_1pct_duplicates_absent
    W04_exclude_50pct_duplicates_absent
    W05_keep_25pct_counts
    W05_keep_50pct_counts
    W05_keep_100pct_counts
    W07_window_pi_overlap
    W07_window_pi_nonoverlap
    W08_tajima_100kb
    W09_site_fst_small_pair
    W09_site_fst_large_pair
    W10_window_fst_biallelic
)
if [[ -n ${CASES:-} ]]; then
    read -r -a cases <<<"$CASES"
fi

[[ -x $ng ]] || { echo "vcftools-ng missing: $ng" >&2; exit 2; }
[[ -f $input ]] || { echo "230k input missing: $input" >&2; exit 2; }
[[ -d $fixtures ]] || { echo "Fixture directory missing: $fixtures" >&2; exit 2; }

configure() {
    local case_id=$1
    case $case_id in
        W03_positions_1pct_sorted|W03_positions_1pct_shuffled|W03_positions_50pct_sorted|W03_positions_50pct_shuffled)
            local name=${case_id#W03_positions_}
            name=${name/_/.}
            args=(--positions "$fixtures/positions-$name.txt" --counts)
            suffix=.frq.count ;;
        W04_exclude_1pct_duplicates_absent)
            args=(--exclude-positions "$fixtures/exclude-1pct.duplicates-and-absent.txt" --counts)
            suffix=.frq.count ;;
        W04_exclude_50pct_duplicates_absent)
            args=(--exclude-positions "$fixtures/exclude-50pct.duplicates-and-absent.txt" --counts)
            suffix=.frq.count ;;
        W05_keep_25pct_counts)
            args=(--keep "$fixtures/keep-25pct.samples.txt" --counts)
            suffix=.frq.count ;;
        W05_keep_50pct_counts)
            args=(--keep "$fixtures/keep-50pct.samples.txt" --counts)
            suffix=.frq.count ;;
        W05_keep_100pct_counts)
            args=(--keep "$fixtures/keep-100pct.samples.txt" --counts)
            suffix=.frq.count ;;
        W07_window_pi_overlap)
            args=(--window-pi 100000 --window-pi-step 10000)
            suffix=.windowed.pi ;;
        W07_window_pi_nonoverlap)
            args=(--window-pi 100000 --window-pi-step 100000)
            suffix=.windowed.pi ;;
        W08_tajima_100kb)
            args=(--TajimaD 100000)
            suffix=.Tajima.D ;;
        W09_site_fst_small_pair)
            args=(--weir-fst-pop "$fixtures/pop-small-ancient-12.txt"
                  --weir-fst-pop "$fixtures/pop-small-asiaticus-13.txt")
            suffix=.weir.fst ;;
        W09_site_fst_large_pair)
            args=(--weir-fst-pop "$fixtures/pop-large-wild-166.txt"
                  --weir-fst-pop "$fixtures/pop-large-cultivated-243.txt")
            suffix=.weir.fst ;;
        W10_window_fst_biallelic)
            args=(--min-alleles 2 --max-alleles 2
                  --weir-fst-pop "$fixtures/pop-large-wild-166.txt"
                  --weir-fst-pop "$fixtures/pop-large-cultivated-243.txt"
                  --fst-window-size 100000 --fst-window-step 10000)
            suffix=.windowed.weir.fst ;;
        *) echo "Unknown 230k workload: $case_id" >&2; exit 2 ;;
    esac
}

if [[ ${GENERATE_ORACLES:-0} == 1 ]]; then
    [[ -x $original ]] || {
        echo "Original VCFtools missing: $original" >&2
        exit 2
    }
    [[ ! -e $oracles ]] || {
        echo "Refusing to overwrite oracle directory: $oracles" >&2
        exit 2
    }
    mkdir -p "$oracles" "$root/original-logs"
    printf '%s\n' $'case\twall_s\tuser_s\tsystem_s\tcpu\tpeak_rss_kib' \
        >"$root/original-runs.tsv"
    for case_id in "${cases[@]}"; do
        configure "$case_id"
        mkdir -p "$oracles/$case_id"
        prefix="$oracles/$case_id/original"
        time_file="$root/original-logs/$case_id.time.tsv"
        /usr/bin/time -f $'%e\t%U\t%S\t%P\t%M' -o "$time_file" \
            "$original" --vcf "$input" "${args[@]}" --out "$prefix" \
            >"$root/original-logs/$case_id.stdout" \
            2>"$root/original-logs/$case_id.stderr"
        mv "$prefix$suffix" "$oracles/$case_id/oracle$suffix"
        read -r wall user system cpu rss <"$time_file"
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$case_id" "$wall" "$user" "$system" "$cpu" "$rss" \
            >>"$root/original-runs.tsv"
    done
    (
        cd "$oracles"
        find . -type f -name 'oracle.*' -print0 | sort -z |
            xargs -0 sha256sum >SHA256SUMS
    )
fi

[[ -f $oracles/SHA256SUMS ]] || {
    echo "230k Original oracle manifest missing; run with GENERATE_ORACLES=1" >&2
    exit 2
}
(
    cd "$oracles"
    sha256sum -c SHA256SUMS
)

[[ ! -e $result ]] || {
    echo "Refusing to overwrite benchmark result: $result" >&2
    exit 2
}
mkdir -p "$result"/{logs,outputs,times}
runs=$result/runs.tsv
printf '%s\n' $'case\tthreads\trepeat\twall_s\tuser_s\tsystem_s\tcpu\tpeak_rss_kib\tbackend\tgate' >"$runs"

dd if="$input" of=/dev/null bs=32M status=none
for case_id in "${cases[@]}"; do
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
            cmp "$oracles/$case_id/oracle$suffix" "$prefix$suffix"
            read -r wall user system cpu rss <"$time_file"
            backend=$(sed -n 's/^Input backend: //p' "$log" | head -1)
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
                "$case_id" "$thread" "$repeat" "$wall" "$user" \
                "$system" "$cpu" "$rss" "${backend:-unknown}" >>"$runs"
            rm -f "$prefix$suffix"
        done
    done
done

echo "230k W03-W10 benchmark complete: $runs"
