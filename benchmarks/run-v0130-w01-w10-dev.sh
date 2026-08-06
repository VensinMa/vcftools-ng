#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 /path/to/vcftools-ng [result-directory]" >&2
    exit 2
fi

ng=$1
root=/home/vensin/workspace/vcftools-ng/benchmarks/results/workload-matrix-23k-v0130
result=${2:-$root/dev-loop}
input=$root/osmanthus412.flags.23chr_1k.vcf
fixtures=$root/fixtures
oracles=$root/oracles
repeats=${REPEATS:-3}
read -r -a threads <<< "${THREADS:-1 4 8 16 32}"
default_cases=(
    W01_counts
    W02_seven_filter_counts
    W03_positions_1pct_sorted
    W03_positions_1pct_shuffled
    W03_positions_50pct_sorted
    W03_positions_50pct_shuffled
    W04_exclude_1pct_duplicates_absent
    W04_exclude_50pct_duplicates_absent
    W05_keep_25pct_counts
    W05_keep_50pct_counts
    W05_keep_100pct_counts
    W06_keep_50pct_filtered_bgzf_recode
    W07_window_pi_overlap
    W07_window_pi_nonoverlap
    W08_tajima_100kb
    W09_site_fst_small_pair
    W09_site_fst_large_pair
    W10_window_fst_biallelic
    W10_window_fst_multiallelic
)
if [[ -n ${CASES:-} ]]; then
    read -r -a cases <<< "$CASES"
else
    cases=("${default_cases[@]}")
fi

[[ -x $ng ]] || { echo "Executable not found: $ng" >&2; exit 2; }
for path in "$input" "$fixtures" "$oracles"; do
    [[ -e $path ]] || { echo "Locked benchmark input missing: $path" >&2; exit 2; }
done

mkdir -p "$result"/{logs,outputs}
summary=$result/runs.tsv
printf '%s\n' $'case\tthreads\trepeat\twall_s\tstage_pipeline_s\tstage_setup_s\tstage_final_s\tbackend\tgate' > "$summary"

configure() {
    case $1 in
        W01_counts)
            case_name=$1
            args=(--counts)
            suffix=.frq.count ;;
        W02_seven_filter_counts)
            case_name=W02_seven_filter_counts
            args=(--min-alleles 2 --max-alleles 2 --minGQ 10 --minQ 30 --min-meanDP 7 --max-missing 0.9 --maf 0.1 --counts)
            suffix=.frq.count ;;
        W03_positions_1pct_sorted|W03_positions_1pct_shuffled|W03_positions_50pct_sorted|W03_positions_50pct_shuffled)
            case_name=$1
            position_name=${case_name#W03_positions_}
            position_name=${position_name/_/.}
            args=(--positions "$fixtures/positions-${position_name}.txt" --counts)
            suffix=.frq.count ;;
        W04_exclude_1pct_duplicates_absent)
            case_name=$1
            args=(--exclude-positions "$fixtures/exclude-1pct.duplicates-and-absent.txt" --counts)
            suffix=.frq.count ;;
        W04_exclude_50pct_duplicates_absent)
            case_name=$1
            args=(--exclude-positions "$fixtures/exclude-50pct.duplicates-and-absent.txt" --counts)
            suffix=.frq.count ;;
        W05_keep_25pct_counts|W05_keep_50pct_counts|W05_keep_100pct_counts)
            case_name=$1
            keep_percent=${case_name#W05_keep_}
            keep_percent=${keep_percent%pct_counts}
            args=(--keep "$fixtures/keep-${keep_percent}pct.samples.txt" --counts)
            suffix=.frq.count ;;
        W06_keep_50pct_filtered_bgzf_recode)
            case_name=$1
            args=(--keep "$fixtures/keep-50pct.samples.txt"
                  --min-alleles 2 --max-alleles 2 --minGQ 10 --minQ 30
                  --min-meanDP 7 --max-missing 0.9 --maf 0.1
                  --recode-vcf-gz --recode-INFO-all)
            suffix=.recode.vcf.gz
            oracle_suffix=.recode.vcf
            comparison=gzip-cmp ;;
        W07_window_pi_overlap)
            case_name=$1
            args=(--window-pi 100000 --window-pi-step 10000)
            suffix=.windowed.pi ;;
        W07_window_pi_nonoverlap)
            case_name=$1
            args=(--window-pi 100000 --window-pi-step 100000)
            suffix=.windowed.pi ;;
        W08_tajima_100kb)
            case_name=$1
            args=(--TajimaD 100000)
            suffix=.Tajima.D ;;
        W09_site_fst_small_pair)
            case_name=$1
            args=(--weir-fst-pop "$fixtures/pop-small-ancient-12.txt" --weir-fst-pop "$fixtures/pop-small-asiaticus-13.txt")
            suffix=.weir.fst ;;
        W09_site_fst_large_pair)
            case_name=$1
            args=(--weir-fst-pop "$fixtures/pop-large-wild-166.txt" --weir-fst-pop "$fixtures/pop-large-cultivated-243.txt")
            suffix=.weir.fst ;;
        W10_window_fst_biallelic)
            case_name=$1
            args=(--min-alleles 2 --max-alleles 2
                  --weir-fst-pop "$fixtures/pop-large-wild-166.txt"
                  --weir-fst-pop "$fixtures/pop-large-cultivated-243.txt"
                  --fst-window-size 100000 --fst-window-step 10000)
            suffix=.windowed.weir.fst ;;
        W10_window_fst_multiallelic)
            case_name=$1
            input=$root/osmanthus205.gatk.23chr_1k.vcf
            args=(--weir-fst-pop "$fixtures/gatk205.pop-wild-141.txt"
                  --weir-fst-pop "$fixtures/gatk205.pop-cultivated-61.txt"
                  --fst-window-size 100000 --fst-window-step 10000)
            suffix=.windowed.weir.fst ;;
        *) echo "Unknown development case: $1" >&2; exit 2 ;;
    esac
    input=${input:-$root/osmanthus412.flags.23chr_1k.vcf}
    oracle_suffix=${oracle_suffix:-$suffix}
    comparison=${comparison:-cmp}
}

stage_value() {
    local name=$1 log=$2
    awk -v target="Stage time [$name]:" \
        '$0 ~ "^Stage time" && index($0, target) == 1 {value=$(NF-1)} END {print value}' \
        "$log"
}

dd if="$input" of=/dev/null bs=16M status=none
for short_case in "${cases[@]}"; do
    unset input oracle_suffix comparison
    configure "$short_case"
    for thread in "${threads[@]}"; do
        for repeat in $(seq 1 "$repeats"); do
            stem=$case_name-t$thread-r$repeat
            prefix=$result/outputs/$stem
            log=$result/logs/$stem.log
            stderr=$result/logs/$stem.stderr
            start=$(date +%s%N)
            "$ng" --vcf "$input" --threads "$thread" --input-backend auto \
                "${args[@]}" --log-file "$log" --out "$prefix" \
                >/dev/null 2>"$stderr"
            end=$(date +%s%N)
            wall=$(awk -v a="$start" -v b="$end" 'BEGIN {printf "%.6f", (b-a)/1000000000}')
            case $comparison in
                cmp)
                    cmp "$oracles/$case_name/oracle$oracle_suffix" "$prefix$suffix"
                    ;;
                gzip-cmp)
                    gzip -dc -- "$prefix$suffix" > "$prefix.decompressed.vcf"
                    cmp "$oracles/$case_name/oracle$oracle_suffix" "$prefix.decompressed.vcf"
                    rm -f "$prefix.decompressed.vcf"
                    ;;
            esac
            pipeline=$(stage_value 'ordered input/compute/commit' "$log")
            [[ -n $pipeline ]] || pipeline=$(stage_value 'fused scan/filter/output' "$log")
            setup=$(stage_value 'pipeline setup' "$log")
            final=$(stage_value 'output finalization' "$log")
            backend=$(sed -n 's/^Input backend: //p' "$log" | head -n 1)
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
                "$case_name" "$thread" "$repeat" "$wall" \
                "${pipeline:-NA}" "${setup:-NA}" "${final:-NA}" \
                "${backend:-unknown}" >> "$summary"
            rm -f "$prefix$suffix"
        done
    done
done

echo "W01-W10 development feedback loop: PASS ($summary)"
