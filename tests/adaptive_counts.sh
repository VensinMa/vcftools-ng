#!/usr/bin/env bash
set -euo pipefail

ng=${1:?vcftools-ng executable required}
source_root=${2:?source root required}
fixture="$source_root/tests/fixtures/osmanthus412.flags.23chr_1k.vcf.gz"
bcf_fixture="$source_root/tests/fixtures/osmanthus412.flags.23chr_1k.bcf"
synthetic="$source_root/tests/fixtures/fast-counts.vcf"
polyploid="$source_root/tests/fixtures/fast-counts-polyploid.vcf"
synthetic_golden="$source_root/tests/golden/fast-counts.frq.count"

work=$(mktemp -d "${TMPDIR:-/tmp}/vcftools-ng-adaptive-counts.XXXXXX")
cleanup() {
    rm -rf -- "$work"
}
trap cleanup EXIT INT TERM

"$ng" --gzvcf "$fixture" --input-backend stream --threads 4 \
    --counts --out "$work/reference" \
    >/dev/null 2>"$work/reference.log"

for threads in 1 2; do
    "$ng" --vcf "$synthetic" --threads "$threads" \
        --counts --out "$work/plain-t$threads" \
        >/dev/null 2>"$work/plain-t$threads.log"
    cmp "$synthetic_golden" "$work/plain-t$threads.frq.count"
    grep -q 'Input backend: fast-counts-plain' \
        "$work/plain-t$threads.log"

    "$ng" --gzvcf "$fixture" --threads "$threads" \
        --counts --out "$work/bgzf-t$threads" \
        >/dev/null 2>"$work/bgzf-t$threads.log"
    cmp "$work/reference.frq.count" "$work/bgzf-t$threads.frq.count"
    grep -q 'Input backend: fast-counts-bgzf' \
        "$work/bgzf-t$threads.log"
done

for threads in 4 8 16 32; do
    "$ng" --vcf "$synthetic" --threads "$threads" \
        --counts --out "$work/plain-t$threads" \
        >/dev/null 2>"$work/plain-t$threads.log"
    cmp "$synthetic_golden" "$work/plain-t$threads.frq.count"
    grep -q 'Input backend: fast-counts-plain' \
        "$work/plain-t$threads.log"
done

for threads in 4 8 16 32; do
    "$ng" --gzvcf "$fixture" --threads "$threads" \
        --counts --out "$work/indexed-t$threads" \
        >/dev/null 2>"$work/indexed-t$threads.log"
    cmp "$work/reference.frq.count" \
        "$work/indexed-t$threads.frq.count"
    grep -q 'Input backend: fast-counts-indexed-bgzf' \
        "$work/indexed-t$threads.log"
done

if "$ng" --vcf "$polyploid" --threads 1 \
    --counts --out "$work/polyploid" \
    >/dev/null 2>"$work/polyploid.log"; then
    printf 'Adaptive counts path unexpectedly accepted polyploid GT\n' >&2
    exit 1
fi
grep -q 'Polyploid genotype is not supported' "$work/polyploid.log"

ln -s "$fixture" "$work/no-index.vcf.gz"
"$ng" --gzvcf "$work/no-index.vcf.gz" --threads 2 \
    --counts --out "$work/no-index-vcf" \
    >/dev/null 2>"$work/no-index-vcf.log"
cmp "$work/reference.frq.count" "$work/no-index-vcf.frq.count"
test ! -e "$work/no-index.vcf.gz.csi"
grep -q 'Input backend: fast-counts-bgzf' \
    "$work/no-index-vcf.log"

ln -s "$fixture" "$work/auto-index.vcf.gz"
"$ng" --gzvcf "$work/auto-index.vcf.gz" --threads 8 \
    --counts --out "$work/auto-index-vcf" \
    >/dev/null 2>"$work/auto-index-vcf.log"
cmp "$work/reference.frq.count" "$work/auto-index-vcf.frq.count"
test -f "$work/auto-index.vcf.gz.csi"
grep -q 'Input backend: fast-counts-indexed-bgzf' \
    "$work/auto-index-vcf.log"

ln -s "$bcf_fixture" "$work/no-index.bcf"
"$ng" --bcf "$work/no-index.bcf" --threads 2 \
    --counts --out "$work/no-index-bcf" \
    >/dev/null 2>"$work/no-index-bcf.log"
cmp "$work/reference.frq.count" "$work/no-index-bcf.frq.count"
test ! -e "$work/no-index.bcf.csi"
grep -q 'Input backend: stream' "$work/no-index-bcf.log"

printf 'Adaptive counts regression test passed\n'
