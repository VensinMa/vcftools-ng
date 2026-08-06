#!/usr/bin/env bash
set -euo pipefail

ng=${1:?vcftools-ng executable required}
source_root=${2:?source root required}
bgzf="$source_root/tests/fixtures/osmanthus412.flags.23chr_1k.vcf.gz"
bcf="$source_root/tests/fixtures/osmanthus412.flags.23chr_1k.bcf"
work=$(mktemp -d "${TMPDIR:-/tmp}/vcftools-ng-index-policy.XXXXXX")
cleanup() {
    rm -rf -- "$work"
}
trap cleanup EXIT INT TERM

readonly recode_args=(
    --min-alleles 2 --max-alleles 2
    --minGQ 10 --minQ 30
    --min-meanDP 7 --max-missing 0.9 --maf 0.1
    --recode --recode-INFO-all --stdout
)

run_recode() {
    local kind=$1
    local input=$2
    local threads=$3
    local name=$4
    shift 4
    "$ng" "--$kind" "$input" --threads "$threads" "$@" \
        "${recode_args[@]}" >"$work/$name.vcf" 2>"$work/$name.log"
}

# Removed compatibility escape hatch: automatic policy is now authoritative.
if "$ng" --gzvcf "$bgzf" --no-auto-index --counts \
    --out "$work/removed-option" \
    >/dev/null 2>"$work/removed-option.log"; then
    printf '%s\n' '--no-auto-index unexpectedly remains accepted' >&2
    exit 1
fi
grep -q 'Unsupported option in phase 1: --no-auto-index' \
    "$work/removed-option.log"

# An existing BCF CSI must be ignored for a full-file recode.
run_recode bcf "$bcf" 8 bcf-auto
run_recode bcf "$bcf" 8 bcf-stream --input-backend stream
cmp "$work/bcf-stream.vcf" "$work/bcf-auto.vcf"
grep -q '^Input backend: stream ' "$work/bcf-auto.log"
grep -q 'BCF full-file recode favors streaming' "$work/bcf-auto.log"

# Missing BCF CSI must not be built for a full-file recode.
cp "$bcf" "$work/full-scan.bcf"
run_recode bcf "$work/full-scan.bcf" 8 bcf-no-index
test ! -e "$work/full-scan.bcf.csi"
grep -q '^Input backend: stream ' "$work/bcf-no-index.log"

# A selective BCF region query should build CSI and use indexed regions.
cp "$bcf" "$work/region.bcf"
run_recode bcf "$work/region.bcf" 1 bcf-region-auto --chr chr1
run_recode bcf "$bcf" 1 bcf-region-stream \
    --input-backend stream --chr chr1
cmp "$work/bcf-region-stream.vcf" "$work/bcf-region-auto.vcf"
test -s "$work/region.bcf.csi"
grep -q '^Input backend: indexed-regions ' "$work/bcf-region-auto.log"

# Existing BGZF indexes are skipped at one thread and used from two threads
# for a full recode.
run_recode gzvcf "$bgzf" 1 bgzf-tbi-t1
run_recode gzvcf "$bgzf" 2 bgzf-tbi-t2
grep -q '^Input backend: fast-filter-recode-bgzf ' \
    "$work/bgzf-tbi-t1.log"
grep -q '^Input backend: fast-filter-recode-indexed-bgzf ' \
    "$work/bgzf-tbi-t2.log"
cmp "$work/bgzf-tbi-t1.vcf" "$work/bgzf-tbi-t2.vcf"

# Missing BGZF index is not built at one thread, but automatic CSI is
# worthwhile for a multi-thread full recode.
cp "$bgzf" "$work/auto-t1.vcf.gz"
run_recode gzvcf "$work/auto-t1.vcf.gz" 1 bgzf-auto-t1
test ! -e "$work/auto-t1.vcf.gz.csi"
grep -q '^Input backend: fast-filter-recode-bgzf ' \
    "$work/bgzf-auto-t1.log"

cp "$bgzf" "$work/auto-t2.vcf.gz"
run_recode gzvcf "$work/auto-t2.vcf.gz" 2 bgzf-auto-t2
test -s "$work/auto-t2.vcf.gz.csi"
grep -q '^Input backend: fast-filter-recode-indexed-bgzf ' \
    "$work/bgzf-auto-t2.log"
cmp "$work/bgzf-auto-t1.vcf" "$work/bgzf-auto-t2.vcf"

# Compact full-scan statistics reuse a valid index only at the measured
# crossover and never build a one-use index.
"$ng" --gzvcf "$bgzf" --threads 2 --counts \
    --out "$work/counts-t2" >/dev/null 2>"$work/counts-t2.log"
"$ng" --gzvcf "$bgzf" --threads 4 --counts \
    --out "$work/counts-t4" >/dev/null 2>"$work/counts-t4.log"
grep -q '^Input backend: fast-counts-bgzf ' "$work/counts-t2.log"
grep -q '^Input backend: fast-counts-indexed-bgzf ' "$work/counts-t4.log"
cmp "$work/counts-t2.frq.count" "$work/counts-t4.frq.count"

cp "$bgzf" "$work/counts-no-index.vcf.gz"
"$ng" --gzvcf "$work/counts-no-index.vcf.gz" --threads 8 --counts \
    --out "$work/counts-no-index" \
    >/dev/null 2>"$work/counts-no-index.log"
test ! -e "$work/counts-no-index.vcf.gz.csi"
grep -q '^Input backend: fast-counts-bgzf ' \
    "$work/counts-no-index.log"

printf 'Adaptive index policy test passed\n'
