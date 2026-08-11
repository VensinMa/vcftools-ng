#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 /path/to/vcftools-ng [result-directory]" >&2
    exit 2
fi

ng=$1
project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
result=${2:-$project_root/benchmarks/results/v0132-development-gate}
matrix=$project_root/benchmarks/results/workload-matrix-23k-v0130
flags_vcf=$matrix/osmanthus412.flags.23chr_1k.vcf
gatk_vcf=$matrix/osmanthus205.gatk.23chr_1k.vcf
flags_bcf=$project_root/tests/fixtures/osmanthus412.flags.23chr_1k.bcf
pca_keep=$project_root/benchmarks/fixtures/v0132-pca-keep-64.txt
original=${ORIGINAL_VCFTOOLS:-/home/vensin/anaconda3/envs/vcftools/bin/vcftools}
reference_ng=${REFERENCE_NG:-/home/vensin/software/vcftools-ng-v0.13.1-linux-x86_64/bin/vcftools-ng}
read -r -a threads <<< "${THREADS:-1 4 8 16 32}"

for path in "$ng" "$flags_vcf" "$gatk_vcf" "$flags_bcf" "$pca_keep"; do
    [[ -e $path ]] || { echo "Required v0.13.2 gate asset missing: $path" >&2; exit 2; }
done

mkdir -p "$result"/{inputs,logs,oracles,outputs}
noindex_vcf=$result/inputs/osmanthus412.flags.23chr_1k.noindex.vcf.gz
if [[ ! -f $noindex_vcf ]]; then
    bgzip -c -- "$flags_vcf" > "$noindex_vcf"
fi
if [[ -e $noindex_vcf.tbi || -e $noindex_vcf.csi ]]; then
    echo "No-index BGZF fixture unexpectedly has a sidecar: $noindex_vcf" >&2
    exit 1
fi

ten_filters=(
    --min-alleles 2 --max-alleles 2
    --minQ 30 --minGQ 10 --minDP 5 --maxDP 30
    --min-meanDP 7 --max-missing 0.9 --maf 0.1 --mac 2
)
site_info_filters=(
    --keep-filtered q10 --keep-filtered PASS --remove-filtered Cluster
    --keep-INFO Hotspot --remove-INFO Artifact
)
sample_ft_filters=(
    --keep "$project_root/tests/fixtures/samples.keep.txt"
    --indv W-DA-8
    --remove "$project_root/tests/fixtures/samples.remove.txt"
    --remove-indv W-CP-5
)

generate_oracles() {
    [[ -x $original ]] || {
        echo "Original VCFtools is required to generate v0.13.2 oracles: $original" >&2
        exit 2
    }
    [[ -x $reference_ng ]] || {
        echo "v0.13.1 reference is required for LD/PCA oracle locking: $reference_ng" >&2
        exit 2
    }
    rm -f -- "$result/oracles"/*
    "$original" --vcf "$flags_vcf" "${ten_filters[@]}" --counts \
        --out "$result/oracles/ten-filter" >/dev/null 2>"$result/logs/oracle-ten-filter.stderr"
    "$original" --vcf "$flags_vcf" --hwe 0.05 --counts \
        --out "$result/oracles/hwe" >/dev/null 2>"$result/logs/oracle-hwe.stderr"
    "$original" --vcf "$flags_vcf" "${site_info_filters[@]}" \
        --recode --recode-INFO-all --out "$result/oracles/site-info" \
        >/dev/null 2>"$result/logs/oracle-site-info.stderr"
    # Original allows only one output mode per invocation. Generate the three
    # FT artifacts once, then lock them under a common oracle prefix so the
    # fused vcftools-ng invocation can be checked against all three.
    "$original" --vcf "$flags_vcf" "${sample_ft_filters[@]}" \
        --remove-filtered-geno-all --counts \
        --out "$result/oracles/ft-counts" \
        >/dev/null 2>"$result/logs/oracle-ft-counts.stderr"
    "$original" --vcf "$flags_vcf" "${sample_ft_filters[@]}" \
        --remove-filtered-geno-all --missing-site \
        --out "$result/oracles/ft-missing" \
        >/dev/null 2>"$result/logs/oracle-ft-missing.stderr"
    "$original" --vcf "$flags_vcf" "${sample_ft_filters[@]}" \
        --remove-filtered-geno-all --recode --recode-INFO-all \
        --out "$result/oracles/ft-recode" \
        >/dev/null 2>"$result/logs/oracle-ft-recode.stderr"
    mv -- "$result/oracles/ft-counts.frq.count" \
        "$result/oracles/ft.frq.count"
    mv -- "$result/oracles/ft-missing.lmiss" \
        "$result/oracles/ft.lmiss"
    mv -- "$result/oracles/ft-recode.recode.vcf" \
        "$result/oracles/ft.recode.vcf"
    "$original" --vcf "$gatk_vcf" --site-pi \
        --out "$result/oracles/site-pi" >/dev/null 2>"$result/logs/oracle-site-pi.stderr"
    "$original" --vcf "$gatk_vcf" --counts \
        --out "$result/oracles/site-pi-counts" >/dev/null 2>"$result/logs/oracle-site-pi-counts.stderr"

    # Original 0.1.17 writes the LD file and then segfaults for this RNC
    # header (OVI-010). v0.13.1 was already byte-gated against the completed
    # Original output, so it is the stable, successful reference here.
    "$reference_ng" --vcf "$flags_vcf" --threads 1 --geno-r2 \
        --ld-window 10 --ld-window-bp 100000 \
        --out "$result/oracles/ld" >/dev/null 2>"$result/logs/oracle-ld.stderr"
    "$reference_ng" --vcf "$flags_vcf" --threads 1 --keep "$pca_keep" \
        --max-missing 1 --pca --out "$result/oracles/pca" \
        >/dev/null 2>"$result/logs/oracle-pca.stderr"
    "$reference_ng" --bcf "$flags_bcf" --diff-bcf "$flags_bcf" \
        --threads 1 --diff-site-discordance --out "$result/oracles/diff" \
        >/dev/null 2>"$result/logs/oracle-diff.stderr"
    (
        cd "$result/oracles"
        sha256sum -- ten-filter.frq.count hwe.frq.count site-info.recode.vcf \
            ft.frq.count ft.lmiss ft.recode.vcf \
            site-pi.sites.pi site-pi-counts.frq.count \
            ld.geno.ld pca.pca diff.diff.sites > SHA256SUMS
    )
}

if [[ ${GENERATE_ORACLES:-0} == 1 ]]; then
    generate_oracles
fi
[[ -s $result/oracles/SHA256SUMS ]] || {
    echo "Locked v0.13.2 oracles are absent; run once with GENERATE_ORACLES=1" >&2
    exit 2
}
(
    cd "$result/oracles"
    sha256sum --check --strict SHA256SUMS
)

summary=$result/runs.tsv
printf '%s\n' $'case\tthreads\twall_s\tbackend\tgate' > "$summary"

run_case() {
    local case_name=$1 thread=$2 input_flag=$3 input_path=$4
    shift 4
    local prefix=$result/outputs/$case_name-t$thread
    local log=$result/logs/$case_name-t$thread.log
    local stderr=$result/logs/$case_name-t$thread.stderr
    local start end wall backend
    start=$(date +%s%N)
    "$ng" "$input_flag" "$input_path" --threads "$thread" "$@" \
        --log-file "$log" --out "$prefix" >/dev/null 2>"$stderr"
    end=$(date +%s%N)
    wall=$(awk -v a="$start" -v b="$end" 'BEGIN {printf "%.6f", (b-a)/1000000000}')
    backend=$(sed -n 's/^Input backend: //p' "$log" | head -n 1)
    if [[ -z $backend ]]; then
        backend=$(sed -n 's/^Diff backend: //p' "$log" | head -n 1)
    fi
    printf '%s\t%s\t%s\t%s\tPASS\n' \
        "$case_name" "$thread" "$wall" "${backend:-unknown}" >> "$summary"
}

for thread in "${threads[@]}"; do
    run_case ten-filter "$thread" --vcf "$flags_vcf" \
        "${ten_filters[@]}" --counts
    cmp "$result/oracles/ten-filter.frq.count" \
        "$result/outputs/ten-filter-t$thread.frq.count"

    run_case hwe "$thread" --vcf "$flags_vcf" --hwe 0.05 --counts
    cmp "$result/oracles/hwe.frq.count" \
        "$result/outputs/hwe-t$thread.frq.count"

    run_case site-info "$thread" --vcf "$flags_vcf" \
        "${site_info_filters[@]}" --recode-vcf --recode-INFO-all
    cmp "$result/oracles/site-info.recode.vcf" \
        "$result/outputs/site-info-t$thread.recode.vcf"

    run_case ft "$thread" --vcf "$flags_vcf" \
        "${sample_ft_filters[@]}" --remove-filtered-geno-all \
        --counts --missing-site --recode-vcf --recode-INFO-all
    cmp "$result/oracles/ft.frq.count" "$result/outputs/ft-t$thread.frq.count"
    cmp "$result/oracles/ft.lmiss" "$result/outputs/ft-t$thread.lmiss"
    cmp "$result/oracles/ft.recode.vcf" "$result/outputs/ft-t$thread.recode.vcf"

    run_case site-pi-multi "$thread" --vcf "$gatk_vcf" --site-pi --counts
    cmp "$result/oracles/site-pi.sites.pi" \
        "$result/outputs/site-pi-multi-t$thread.sites.pi"
    cmp "$result/oracles/site-pi-counts.frq.count" \
        "$result/outputs/site-pi-multi-t$thread.frq.count"

    run_case noindex-bgzf "$thread" --gzvcf "$noindex_vcf" \
        "${ten_filters[@]}" --counts
    cmp "$result/oracles/ten-filter.frq.count" \
        "$result/outputs/noindex-bgzf-t$thread.frq.count"

    run_case ld "$thread" --vcf "$flags_vcf" --geno-r2 \
        --ld-window 10 --ld-window-bp 100000
    cmp "$result/oracles/ld.geno.ld" "$result/outputs/ld-t$thread.geno.ld"

    run_case pca "$thread" --vcf "$flags_vcf" --keep "$pca_keep" \
        --max-missing 1 --pca
    cmp "$result/oracles/pca.pca" "$result/outputs/pca-t$thread.pca"

    run_case indexed-diff "$thread" --bcf "$flags_bcf" \
        --diff-bcf "$flags_bcf" --diff-site-discordance
    cmp "$result/oracles/diff.diff.sites" \
        "$result/outputs/indexed-diff-t$thread.diff.sites"
done

echo "v0.13.2 23k exact development gate: PASS ($summary)"
