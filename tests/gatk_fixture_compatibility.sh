#!/usr/bin/env bash
set -euo pipefail

ng=${1:?vcftools-ng executable required}
source_root=${2:?source root required}
fixture="$source_root/tests/fixtures/osmanthus205.gatk.23chr_1k.vcf.gz"
bcf_fixture="$source_root/tests/fixtures/osmanthus205.gatk.23chr_1k.bcf"
golden_prefix="$source_root/tests/golden/gatk205-seven-filter"

required=(
    "$fixture" "$fixture.tbi" "$bcf_fixture" "$bcf_fixture.csi"
    "$golden_prefix.counts.frq.count"
    "$golden_prefix.recode.recode.vcf"
)
for path in "${required[@]}"; do
    if [[ ! -s $path ]]; then
        printf 'Missing local GATK compatibility artifact: %s\n' "$path" >&2
        exit 1
    fi
done

work=$(mktemp -d "${TMPDIR:-/tmp}/vcftools-ng-gatk-gate.XXXXXX")
cleanup() {
    rm -rf -- "$work"
}
trap cleanup EXIT INT TERM

filters=(
    --min-alleles 2 --max-alleles 2
    --minGQ 10 --minQ 30 --min-meanDP 7
    --max-missing 0.9 --maf 0.1
)

for threads in 1 8 32; do
    "$ng" --gzvcf "$fixture" --threads "$threads" \
        "${filters[@]}" --counts --recode-vcf --recode-INFO-all \
        --no-log-file --out "$work/vcf-t$threads" \
        >/dev/null 2>"$work/vcf-t$threads.stderr"
    cmp "$golden_prefix.counts.frq.count" \
        "$work/vcf-t$threads.frq.count"
    cmp "$golden_prefix.recode.recode.vcf" \
        "$work/vcf-t$threads.recode.vcf"
    grep -q 'Input backend: fast-filter-recode-' \
        "$work/vcf-t$threads.stderr"
done

for parser in generic specialized; do
    VCFTOOLS_NG_TEST_PARSER=$parser \
        "$ng" --gzvcf "$fixture" --threads 8 \
        "${filters[@]}" --counts --no-log-file \
        --out "$work/parser-$parser" \
        >/dev/null 2>"$work/parser-$parser.stderr"
    cmp "$golden_prefix.counts.frq.count" \
        "$work/parser-$parser.frq.count"
done
cmp "$work/parser-generic.frq.count" \
    "$work/parser-specialized.frq.count"

"$ng" --gzvcf "$fixture" --threads 8 \
    "${filters[@]}" --recode --recode-INFO-all \
    --no-log-file --out "$work/bgzf" \
    >/dev/null 2>"$work/bgzf.stderr"
gzip -dc -- "$work/bgzf.recode.vcf.gz" >"$work/bgzf.recode.vcf"
cmp "$golden_prefix.recode.recode.vcf" "$work/bgzf.recode.vcf"

"$ng" --bcf "$bcf_fixture" --threads 8 \
    "${filters[@]}" --counts \
    --no-log-file --out "$work/bcf" \
    >/dev/null 2>"$work/bcf.stderr"
cmp "$golden_prefix.counts.frq.count" "$work/bcf.frq.count"
grep -Eq 'Input backend: (stream|indexed-regions)' "$work/bcf.stderr"
if grep -q 'fast-filter-recode-' "$work/bcf.stderr"; then
    printf 'BCF input unexpectedly entered the direct VCF text kernel\n' >&2
    exit 1
fi

printf 'GATK 23,000-record Original compatibility gate passed\n'
