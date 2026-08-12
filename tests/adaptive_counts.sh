#!/usr/bin/env bash
set -euo pipefail

ng=${1:?vcftools-ng executable required}
source_root=${2:?source root required}
fixture="$source_root/tests/fixtures/osmanthus412.flags.23chr_1k.vcf.gz"
bcf_fixture="$source_root/tests/fixtures/osmanthus412.flags.23chr_1k.bcf"
synthetic="$source_root/tests/fixtures/fast-counts.vcf"
genotype_forms="$source_root/tests/fixtures/genotype-forms.vcf"
site_stats="$source_root/tests/fixtures/fast-site-stats.vcf"
fast_recode="$source_root/tests/fixtures/fast-recode.vcf"
polyploid="$source_root/tests/fixtures/fast-counts-polyploid.vcf"
polyploid_positions="$source_root/tests/fixtures/polyploid.positions.txt"
synthetic_golden="$source_root/tests/golden/fast-counts.frq.count"
genotype_forms_golden="$source_root/tests/golden/genotype-forms.frq.count"
genotype_forms_missing_golden="$source_root/tests/golden/genotype-forms.lmiss"
genotype_forms_indv_golden="$source_root/tests/golden/genotype-forms.imiss"
genotype_forms_max_missing_golden="$source_root/tests/golden/genotype-forms.max-missing.frq.count"
site_stats_golden="$source_root/tests/golden/fast-site-stats"
fast_recode_info_golden="$source_root/tests/golden/fast-recode.info.vcf"
fast_recode_clear_golden="$source_root/tests/golden/fast-recode.clear-info.vcf"
polyploid_site_only_golden="$source_root/tests/golden/polyploid-site-only.recode.vcf"
polyploid_positions_golden="$source_root/tests/golden/polyploid-positions.recode.vcf"
filtered_stdout_sha256=292684f4994507b09b1ac339dcb03211293e9dd68fcf9b309ed006e1968818c4
filtered_bgzf_sha256=d8455a9642c220045d95705a74beb794faa9b6eda58ac22c792939f419236e14
real_filters=(
    --min-alleles 2 --max-alleles 2
    --minGQ 10 --minQ 30 --min-meanDP 7
    --max-missing 0.9 --maf 0.1
)

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

for threads in 1 4 8 16 32; do
    "$ng" --vcf "$fast_recode" --threads "$threads" \
        --minGQ 10 --recode-vcf --recode-INFO-all \
        --out "$work/recode-edge-info-t$threads" \
        >/dev/null 2>"$work/recode-edge-info-t$threads.log"
    cmp "$fast_recode_info_golden" \
        "$work/recode-edge-info-t$threads.recode.vcf"
    grep -q 'Input backend: fast-filter-recode-plain' \
        "$work/recode-edge-info-t$threads.log"
done

for threads in 1 32; do
    "$ng" --vcf "$genotype_forms" --threads "$threads" \
        --missing-indv --out "$work/genotype-forms-indv-t$threads" \
        >/dev/null 2>"$work/genotype-forms-indv-t$threads.log"
    cmp "$genotype_forms_indv_golden" \
        "$work/genotype-forms-indv-t$threads.imiss"
done

# Original's .lmiss haploid quirk must not leak into --max-missing, which
# still counts a three-character phased partial call as diploid (OVI-011).
for threads in 1 32; do
    "$ng" --vcf "$genotype_forms" --threads "$threads" \
        --max-missing 0.7 --counts \
        --out "$work/genotype-forms-max-missing-t$threads" \
        >/dev/null 2>"$work/genotype-forms-max-missing-t$threads.log"
    cmp "$genotype_forms_max_missing_golden" \
        "$work/genotype-forms-max-missing-t$threads.frq.count"
done

"$ng" --vcf "$fast_recode" --threads 8 \
    --minGQ 10 --recode-vcf \
    --out "$work/recode-edge-clear" \
    >/dev/null 2>"$work/recode-edge-clear.log"
cmp "$fast_recode_clear_golden" \
    "$work/recode-edge-clear.recode.vcf"

"$ng" --vcf "$fast_recode" --threads 8 \
    --minGQ 10 --recode-vcf --recode-vcf-gz --recode-INFO-all \
    --out "$work/recode-edge-dual" \
    >/dev/null 2>"$work/recode-edge-dual.log"
cmp "$fast_recode_info_golden" \
    "$work/recode-edge-dual.recode.vcf"
gzip -dc -- "$work/recode-edge-dual.recode.vcf.gz" \
    >"$work/recode-edge-dual.uncompressed.vcf"
cmp "$fast_recode_info_golden" \
    "$work/recode-edge-dual.uncompressed.vcf"

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

# The large-file planner uses the same explicit pread path.  Force both
# access modes on a compact fixture so their scientific output remains a
# deterministic differential gate independent of host RAM and file size.
VCFTOOLS_NG_TEST_PLAIN_ACCESS=pread \
    "$ng" --vcf "$synthetic" --threads 4 --counts \
    --out "$work/plain-forced-pread" \
    >/dev/null 2>"$work/plain-forced-pread.log"
VCFTOOLS_NG_TEST_PLAIN_ACCESS=mmap \
    "$ng" --vcf "$synthetic" --threads 4 --counts \
    --out "$work/plain-forced-mmap" \
    >/dev/null 2>"$work/plain-forced-mmap.log"
cmp "$work/plain-forced-pread.frq.count" \
    "$work/plain-forced-mmap.frq.count"
grep -q 'fused aligned byte ranges' \
    "$work/plain-forced-pread.log"
grep -q 'zero-copy mapped aligned ranges' \
    "$work/plain-forced-mmap.log"

for threads in 4 8 16 32; do
    "$ng" --vcf "$synthetic" --threads "$threads" \
        --counts --out "$work/plain-t$threads" \
        >/dev/null 2>"$work/plain-t$threads.log"
    cmp "$synthetic_golden" "$work/plain-t$threads.frq.count"
    grep -q 'Input backend: fast-counts-plain' \
        "$work/plain-t$threads.log"
done

# Cover every GT family observed in the 23-chromosome DeepVariant/GATK
# fixtures, plus valid reverse, trailing-partial-missing, phased-partial,
# haploid-called, and multi-digit-allele spellings.  The golden was generated
# by Original VCFtools 0.1.17.
for threads in 1 4 32; do
    "$ng" --vcf "$genotype_forms" --threads "$threads" \
        --counts --missing-site --out "$work/genotype-forms-t$threads" \
        >/dev/null 2>"$work/genotype-forms-t$threads.log"
    cmp "$genotype_forms_golden" \
        "$work/genotype-forms-t$threads.frq.count"
    cmp "$genotype_forms_missing_golden" \
        "$work/genotype-forms-t$threads.lmiss"
    grep -q 'Input backend: fast-site-stats-plain' \
        "$work/genotype-forms-t$threads.log"
done

# Force both scalar parser policies over the same GT grammar fixture. This is
# the baseline differential that future AVX2/AVX-512 parsers must also pass.
for parser in generic specialized; do
    VCFTOOLS_NG_TEST_PARSER=$parser \
        "$ng" --vcf "$genotype_forms" --threads 4 \
        --counts --missing-site \
        --out "$work/genotype-parser-$parser" \
        >/dev/null 2>"$work/genotype-parser-$parser.log"
    cmp "$genotype_forms_golden" \
        "$work/genotype-parser-$parser.frq.count"
    cmp "$genotype_forms_missing_golden" \
        "$work/genotype-parser-$parser.lmiss"
done
cmp "$work/genotype-parser-generic.frq.count" \
    "$work/genotype-parser-specialized.frq.count"
cmp "$work/genotype-parser-generic.lmiss" \
    "$work/genotype-parser-specialized.lmiss"

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
        "${real_filters[@]}" \
        --recode --recode-INFO-all \
        --out "$work/filtered-bgzf-recode-t$threads" \
        >/dev/null 2>"$work/filtered-bgzf-recode-t$threads.log"
    test "$filtered_bgzf_sha256" = "$(
        sha256sum \
            "$work/filtered-bgzf-recode-t$threads.recode.vcf.gz" |
            cut -d' ' -f1
    )"
    grep -q 'Input backend: fast-filter-recode-' \
        "$work/filtered-bgzf-recode-t$threads.log"
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

"$ng" --bcf "$bcf_fixture" --input-backend stream --threads 4 \
    "${real_filters[@]}" --counts --out "$work/filtered-reference" \
    >/dev/null 2>"$work/filtered-reference.log"

for parser in generic specialized; do
    VCFTOOLS_NG_TEST_PARSER=$parser \
        "$ng" --gzvcf "$fixture" --threads 4 \
        "${real_filters[@]}" --counts \
        --out "$work/deepvariant-parser-$parser" \
        >/dev/null 2>"$work/deepvariant-parser-$parser.log"
    cmp "$work/filtered-reference.frq.count" \
        "$work/deepvariant-parser-$parser.frq.count"
done
cmp "$work/deepvariant-parser-generic.frq.count" \
    "$work/deepvariant-parser-specialized.frq.count"

if VCFTOOLS_NG_TEST_PARSER=invalid \
    "$ng" --vcf "$synthetic" --threads 1 --counts \
    --out "$work/parser-invalid" \
    >/dev/null 2>"$work/parser-invalid.log"; then
    printf 'Invalid parser test override unexpectedly succeeded\n' >&2
    exit 1
fi
grep -q 'VCFTOOLS_NG_TEST_PARSER must be' "$work/parser-invalid.log"
test ! -e "$work/parser-invalid.frq.count"

"$ng" --gzvcf "$fixture" --threads 8 \
    "${real_filters[@]}" --counts --recode --recode-INFO-all \
    --out "$work/filtered-combined" \
    >/dev/null 2>"$work/filtered-combined.log"
cmp "$work/filtered-reference.frq.count" \
    "$work/filtered-combined.frq.count"
test "$filtered_bgzf_sha256" = "$(
    sha256sum "$work/filtered-combined.recode.vcf.gz" |
        cut -d' ' -f1
)"
grep -q 'Input backend: fast-filter-recode-' \
    "$work/filtered-combined.log"
gzip -dc -- "$fixture" >"$work/filtered.vcf"
for threads in 1 4 32; do
    "$ng" --vcf "$work/filtered.vcf" --threads "$threads" \
        "${real_filters[@]}" --counts \
        --out "$work/filtered-plain-t$threads" \
        >/dev/null 2>"$work/filtered-plain-t$threads.log"
    cmp "$work/filtered-reference.frq.count" \
        "$work/filtered-plain-t$threads.frq.count"
    grep -q 'Input backend: fast-counts-plain' \
        "$work/filtered-plain-t$threads.log"

    "$ng" --gzvcf "$fixture" --threads "$threads" \
        "${real_filters[@]}" --counts \
        --out "$work/filtered-bgzf-t$threads" \
        >/dev/null 2>"$work/filtered-bgzf-t$threads.log"
    cmp "$work/filtered-reference.frq.count" \
        "$work/filtered-bgzf-t$threads.frq.count"
    grep -q 'Input backend: fast-counts-' \
        "$work/filtered-bgzf-t$threads.log"
done

# Site-only operations keep GT opaque. Polyploid sample text therefore
# survives raw recode even though Original 0.1.17 rejects the same safe
# operation before writing records.
for threads in 1 32; do
    "$ng" --vcf "$polyploid" --threads "$threads" \
        --minQ 30 --recode-vcf --recode-INFO-all \
        --out "$work/polyploid-site-only-t$threads" \
        >/dev/null 2>"$work/polyploid-site-only-t$threads.log"
    cmp "$polyploid_site_only_golden" \
        "$work/polyploid-site-only-t$threads.recode.vcf"

    "$ng" --vcf "$polyploid" --threads "$threads" \
        --positions "$polyploid_positions" \
        --recode-vcf --recode-INFO-all \
        --out "$work/polyploid-positions-t$threads" \
        >/dev/null 2>"$work/polyploid-positions-t$threads.log"
    cmp "$polyploid_positions_golden" \
        "$work/polyploid-positions-t$threads.recode.vcf"
done

# Every representative triploid spelling must fail before publication when
# GT semantics are requested. Coordinate restriction also exercises the
# generic HTSlib pipeline rather than only the fused text parser.
for position in 1 2 3 4 5; do
    prefix="$work/polyploid-semantic-$position"
    if "$ng" --vcf "$polyploid" --threads 1 \
        --chr chr1 --from-bp "$position" --to-bp "$position" --counts \
        --out "$prefix" >/dev/null 2>"$prefix.log"; then
        printf 'Semantic path accepted polyploid GT at position %s\n' \
            "$position" >&2
        exit 1
    fi
    grep -q 'Polyploid genotype is not supported' "$prefix.log"
    test ! -e "$prefix.frq.count"
done

if "$ng" --vcf "$polyploid" --threads 1 \
    --minGQ 30 --recode-vcf --out "$work/polyploid-masking" \
    >/dev/null 2>"$work/polyploid-masking.log"; then
    printf 'Genotype masking accepted polyploid GT\n' >&2
    exit 1
fi
grep -q 'Polyploid genotype is not supported' \
    "$work/polyploid-masking.log"
test ! -e "$work/polyploid-masking.recode.vcf"

ln -s "$fixture" "$work/no-index.vcf.gz"
"$ng" --gzvcf "$work/no-index.vcf.gz" --threads 2 \
    --counts --out "$work/no-index-vcf" \
    >/dev/null 2>"$work/no-index-vcf.log"
cmp "$work/reference.frq.count" "$work/no-index-vcf.frq.count"
test ! -e "$work/no-index.vcf.gz.csi"
grep -q 'Input backend: fast-counts-bgzf' \
    "$work/no-index-vcf.log"

# A serial compressed feeder must not wake more compute workers than its
# bounded in-flight queue can use. GT-only work caps at four workers; FORMAT
# and frequency-heavy filtering caps at six. The strict total budget still
# governs smaller thread requests.
"$ng" --gzvcf "$work/no-index.vcf.gz" --threads 16 --counts \
    --out "$work/no-index-cap-light" \
    >/dev/null 2>"$work/no-index-cap-light.log"
cmp "$work/reference.frq.count" \
    "$work/no-index-cap-light.frq.count"
grep -Fq 'Stage concurrency: input 1, HTSlib I/O 3, HTSlib coordinator 1, compute 4' \
    "$work/no-index-cap-light.log"

"$ng" --gzvcf "$work/no-index.vcf.gz" --threads 16 \
    "${real_filters[@]}" --counts \
    --out "$work/no-index-cap-heavy" \
    >/dev/null 2>"$work/no-index-cap-heavy.log"
cmp "$work/filtered-reference.frq.count" \
    "$work/no-index-cap-heavy.frq.count"
grep -Fq 'Stage concurrency: input 1, HTSlib I/O 3, HTSlib coordinator 1, compute 6' \
    "$work/no-index-cap-heavy.log"

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

# BCF record decoding needs a larger HTSlib share than a BGZF VCF text
# stream.  The generic ordered pipeline also reserves the calling thread as
# its committer, so input + HTSlib + compute must be no greater than the
# requested total and, from three threads upward, exactly one below it.  The
# quarter-share curve is evidence-selected on the locked 230k BCF fixture and
# grows beyond the local 32-CPU host without encoding a machine-specific
# ceiling.
for spec in 2:1:0:0 3:1:0:1 4:1:0:2 8:1:2:3 16:1:4:9 32:1:8:21; do
    threads=${spec%%:*}
    remainder=${spec#*:}
    input_threads=${remainder%%:*}
    remainder=${remainder#*:}
    hts_io=${remainder%%:*}
    compute=${remainder#*:}
    prefix="$work/bcf-resource-t$threads"
    "$ng" --bcf "$bcf_fixture" --input-backend stream \
        --threads "$threads" --counts --out "$prefix" \
        >/dev/null 2>"$prefix.log"
    cmp "$work/reference.frq.count" "$prefix.frq.count"
    hts_coordinator=0
    if ((hts_io > 0)); then
        hts_coordinator=1
    fi
    grep -Fq \
        "Stage concurrency: input $input_threads, HTSlib I/O $hts_io, HTSlib coordinator $hts_coordinator, compute $compute" \
        "$prefix.log"
    planned=$((input_threads + hts_io + hts_coordinator + compute))
    if ((planned > threads)); then
        printf 'BCF plan exceeds budget: %s > %s\n' \
            "$planned" "$threads" >&2
        exit 1
    fi
    if ((threads >= 3 && planned + 1 != threads)); then
        printf 'BCF plan does not reserve exactly one committer: %s/%s\n' \
            "$planned" "$threads" >&2
        exit 1
    fi
done

printf 'Adaptive counts regression test passed\n'
