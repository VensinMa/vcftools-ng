#!/usr/bin/env bash
set -euo pipefail

ng=${1:?vcftools-ng executable required}
source_root=${2:?source root required}
fixture="$source_root/tests/fixtures/osmanthus412.flags.23chr_1k.vcf.gz"
bcf_fixture="$source_root/tests/fixtures/osmanthus412.flags.23chr_1k.bcf"
synthetic="$source_root/tests/fixtures/fast-counts.vcf"
site_stats="$source_root/tests/fixtures/fast-site-stats.vcf"
polyploid="$source_root/tests/fixtures/fast-counts-polyploid.vcf"
synthetic_golden="$source_root/tests/golden/fast-counts.frq.count"
site_stats_golden="$source_root/tests/golden/fast-site-stats"
filtered_stdout_sha256=292684f4994507b09b1ac339dcb03211293e9dd68fcf9b309ed006e1968818c4

work=$(mktemp -d "${TMPDIR:-/tmp}/vcftools-ng-adaptive-counts.XXXXXX")
cleanup() {
    rm -rf -- "$work"
}
trap cleanup EXIT INT TERM

compare_site_stats() {
    local prefix=$1
    cmp "$site_stats_golden.frq" "$prefix.frq"
    cmp "$site_stats_golden.frq.count" "$prefix.frq.count"
    cmp "$site_stats_golden.lmiss" "$prefix.lmiss"
    cmp "$site_stats_golden.ldepth" "$prefix.ldepth"
    cmp "$site_stats_golden.ldepth.mean" "$prefix.ldepth.mean"
    cmp "$site_stats_golden.lqual" "$prefix.lqual"
}

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

for threads in 1 2 4 8 16 32; do
    "$ng" --vcf "$site_stats" --threads "$threads" \
        --freq --counts --missing-site --site-depth \
        --site-mean-depth --site-quality \
        --out "$work/site-stats-t$threads" \
        >/dev/null 2>"$work/site-stats-t$threads.log"
    compare_site_stats "$work/site-stats-t$threads"
    grep -q 'Input backend: fast-site-stats-plain' \
        "$work/site-stats-t$threads.log"
done

for threads in 1 32; do
    "$ng" --vcf "$site_stats" --threads "$threads" \
        --freq2 --out "$work/freq2-t$threads" \
        >/dev/null 2>"$work/freq2-t$threads.log"
    cmp "$site_stats_golden.freq2.frq" \
        "$work/freq2-t$threads.frq"
done

for threads in 1 32; do
    "$ng" --gzvcf "$fixture" --threads "$threads" \
        --min-alleles 2 --max-alleles 2 \
        --minGQ 10 --minQ 30 --min-meanDP 7 \
        --max-missing 0.9 --maf 0.1 \
        --recode --recode-INFO-all --stdout \
        >"$work/filtered-stdout-t$threads.vcf" \
        2>"$work/filtered-stdout-t$threads.log"
    test "$filtered_stdout_sha256" = "$(
        sha256sum "$work/filtered-stdout-t$threads.vcf" |
            cut -d' ' -f1
    )"
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
"$ng" --gzvcf "$work/auto-index.vcf.gz" --threads 2 \
    --counts --out "$work/auto-index-vcf" \
    >/dev/null 2>"$work/auto-index-vcf.log"
cmp "$work/reference.frq.count" "$work/auto-index-vcf.frq.count"
test ! -e "$work/auto-index.vcf.gz.csi"
grep -q 'Input backend: fast-counts-bgzf' \
    "$work/auto-index-vcf.log"

ln -s "$bcf_fixture" "$work/auto-index.bcf"
"$ng" --bcf "$work/auto-index.bcf" --threads 1 \
    --counts --out "$work/auto-index-bcf" \
    >/dev/null 2>"$work/auto-index-bcf.log"
cmp "$work/reference.frq.count" \
    "$work/auto-index-bcf.frq.count"
test ! -e "$work/auto-index.bcf.csi"
grep -q 'Input backend: stream' \
    "$work/auto-index-bcf.log"

ln -s "$bcf_fixture" "$work/no-index.bcf"
"$ng" --bcf "$work/no-index.bcf" --threads 2 \
    --counts --out "$work/no-index-bcf" \
    >/dev/null 2>"$work/no-index-bcf.log"
cmp "$work/reference.frq.count" "$work/no-index-bcf.frq.count"
test ! -e "$work/no-index.bcf.csi"
grep -q 'Input backend: stream' "$work/no-index-bcf.log"

printf 'Adaptive counts regression test passed\n'
