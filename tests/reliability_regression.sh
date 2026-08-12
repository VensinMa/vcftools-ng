#!/usr/bin/env bash
set -euo pipefail

ng=${1:?vcftools-ng executable required}
source_root=${2:?source root required}
fixture="$source_root/tests/fixtures/fast-site-stats.vcf"
numeric_fixture="$source_root/tests/fixtures/numeric-edge.vcf"
numeric_golden="$source_root/tests/golden/numeric-edge"
no_contig_fixture="$source_root/tests/fixtures/no-contig.vcf"
numeric_pop1="$source_root/tests/fixtures/numeric-population-1.txt"
numeric_pop2="$source_root/tests/fixtures/numeric-population-2.txt"
multiallelic_fixture="$source_root/tests/fixtures/multiallelic-cutovers.vcf"
multiallelic_golden="$source_root/tests/golden/multiallelic-cutovers"
mixed_ploidy_fixture="$source_root/tests/fixtures/mixed-ploidy.vcf"
unsorted_fixture="$source_root/tests/fixtures/unsorted-sites.vcf"
work=$(mktemp -d "${TMPDIR:-/tmp}/vcftools-ng-reliability.XXXXXX")
cleanup() {
    rm -rf -- "$work"
}
trap cleanup EXIT INT TERM

# Generate the compressed failure-injection fixture from committed text
# instead of depending on the large local Osmanthus fixture.  The latter is
# intentionally excluded from Git, so referencing it here made this test pass
# in the development workspace but fail in every clean CI checkout.  Varying
# positions keep the BGZF output comfortably larger than the truncation size.
bgzf_source="$work/generated-bgzf-source.vcf"
awk 'BEGIN {
    print "##fileformat=VCFv4.2"
    print "##contig=<ID=chr1,length=50000>"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    for (i = 1; i <= 50000; ++i) {
        print "chr1\t" i "\t.\tA\tG\t30\tPASS\t."
    }
}' >"$bgzf_source"
"$ng" --vcf "$bgzf_source" --threads 1 --recode \
    --out "$work/generated-bgzf" \
    >"$work/generated-bgzf.stdout" 2>"$work/generated-bgzf.stderr"
bgzf_fixture="$work/generated-bgzf.recode.vcf.gz"
test -s "$bgzf_fixture"
test "$(stat -c %s "$bgzf_fixture")" -gt 8192

# A failed multi-output publication must preserve every old destination and
# remove all private staged files.
printf 'previous scientific result\n' >"$work/transaction.lqual"
printf 'previous depth result\n' >"$work/transaction.ldepth"
if VCFTOOLS_NG_TEST_FAIL_OUTPUT_COMMIT_AFTER=1 \
    "$ng" --vcf "$fixture" --threads 1 --site-quality --site-depth \
    --out "$work/transaction" \
    >"$work/transaction.stdout" 2>"$work/transaction.stderr"; then
    echo 'Injected output commit failure unexpectedly succeeded' >&2
    exit 1
fi
grep -Fqx 'previous scientific result' "$work/transaction.lqual"
grep -Fqx 'previous depth result' "$work/transaction.ldepth"
grep -Fq 'Injected scientific-output partial commit failure' \
    "$work/transaction.stderr"
grep -Fq 'Exit status: failed' "$work/transaction.stderr"
test -z "$(find "$work" -name '*.vcftools-ng.tmp.*' -print -quit)"
test -z "$(find "$work" -name '*.vcftools-ng.backup.*' -print -quit)"

# If rollback itself fails, the private backup must remain recoverable and the
# log must name both the preserved backup and intended destination.
printf 'recoverable old result\n' >"$work/restore.lqual"
printf 'recoverable old depth\n' >"$work/restore.ldepth"
if VCFTOOLS_NG_TEST_FAIL_OUTPUT_COMMIT_AFTER=1 \
   VCFTOOLS_NG_TEST_FAIL_OUTPUT_RESTORE=1 \
    "$ng" --vcf "$fixture" --threads 1 --site-quality --site-depth \
    --out "$work/restore" \
    >"$work/restore.stdout" 2>"$work/restore.stderr"; then
    echo 'Injected output restore failure unexpectedly succeeded' >&2
    exit 1
fi
grep -Fq 'output recovery is incomplete; restore' "$work/restore.stderr"
restore_backup=$(find "$work" -name 'restore.lqual.vcftools-ng.backup.*' \
    -print -quit)
test -n "$restore_backup"
grep -Fqx 'recoverable old result' "$restore_backup"
grep -Fq "$restore_backup" "$work/restore.stderr"
grep -Fq "$work/restore.lqual" "$work/restore.stderr"

# A derived output name must never be allowed to replace any input file.
cp "$fixture" "$work/protected.frq"
protected_hash=$(sha256sum "$work/protected.frq" | awk '{print $1}')
if "$ng" --vcf "$work/protected.frq" --threads 1 --freq \
    --out "$work/protected" \
    >"$work/protected.stdout" 2>"$work/protected.stderr"; then
    echo 'Input/output destination collision unexpectedly succeeded' >&2
    exit 1
fi
grep -Fq 'Scientific output must not overwrite an input file' \
    "$work/protected.stderr"
test "$protected_hash" = "$(
    sha256sum "$work/protected.frq" | awk '{print $1}'
)"

# A failed log mirror must never poison stderr or the scientific result.
"$ng" --vcf "$fixture" --threads 1 --counts \
    --log-file /dev/full --out "$work/log-full" \
    >"$work/log-full.stdout" 2>"$work/log-full.stderr"
test -s "$work/log-full.frq.count"
grep -Fq 'Warning: log file write failed' "$work/log-full.stderr"
grep -Fq 'Log status: incomplete' "$work/log-full.stderr"
grep -Fq 'Exit status: success' "$work/log-full.stderr"

# Every floating-point option rejects non-finite and overflowing spellings
# during argument parsing, before a scientific output is staged.
float_options=(
    --min-r2 --minQ --minGQ --min-meanDP --max-meanDP
    --max-missing --maf --max-maf --hwe --non-ref-af
    --max-non-ref-af --non-ref-af-any --max-non-ref-af-any
)
for option in "${float_options[@]}"; do
    for value in nan NaN inf -inf 1e9999; do
        prefix="$work/nonfinite-${option#--}-${value//[^A-Za-z0-9]/_}"
        if "$ng" --vcf "$fixture" --threads 1 "$option" "$value" \
            --counts --out "$prefix" \
            >"$prefix.stdout" 2>"$prefix.stderr"; then
            echo "Non-finite value unexpectedly accepted: $option $value" \
                >&2
            exit 1
        fi
        grep -Fq "Invalid finite value for $option" "$prefix.stderr"
        test ! -e "$prefix.frq.count"
    done
done

# Automatic selection is capped at 128, while an explicit request may exceed
# 128 when the detected allocation permits it.
VCFTOOLS_NG_TEST_AVAILABLE_THREADS=256 \
    "$ng" --vcf "$fixture" --counts --out "$work/auto-cap" \
    >"$work/auto-cap.stdout" 2>"$work/auto-cap.stderr"
grep -Fq 'Requested threads: 256' "$work/auto-cap.stderr"
grep -Fq 'Effective threads: 128' "$work/auto-cap.stderr"
grep -Fq 'capped at 128' "$work/auto-cap.stderr"

VCFTOOLS_NG_TEST_AVAILABLE_THREADS=256 \
    "$ng" --vcf "$fixture" --threads 200 --counts \
    --out "$work/explicit-200" \
    >"$work/explicit-200.stdout" 2>"$work/explicit-200.stderr"
grep -Fq 'Effective threads: 200' "$work/explicit-200.stderr"

# The parallel input and compute allocations share, rather than duplicate,
# the requested CPU budget.
VCFTOOLS_NG_TEST_AVAILABLE_THREADS=8 \
    "$ng" --vcf "$fixture" --threads 8 --depth \
    --out "$work/budget" \
    >"$work/budget.stdout" 2>"$work/budget.stderr"
input_threads=$(awk '/^Input threads:/{print $3}' "$work/budget.stderr")
compute_threads=$(awk '/^Compute threads:/{print $3}' "$work/budget.stderr")
test "$((input_threads + compute_threads))" -le 8

# The file-descriptor ceiling is intersected with, never substituted for, the
# requested input-worker budget.
(
    ulimit -n 66
    VCFTOOLS_NG_TEST_AVAILABLE_THREADS=8 \
        "$ng" --vcf "$fixture" --input-backend plain --threads 8 --depth \
        --out "$work/fd-cap" \
        >"$work/fd-cap.stdout" 2>"$work/fd-cap.stderr"
)
grep -Fq 'Input threads: 1' "$work/fd-cap.stderr"
fd_compute_threads=$(
    awk '/^Compute threads:/{print $3}' "$work/fd-cap.stderr"
)
test "$((1 + fd_compute_threads))" -le 8

# The compressed-stream worker-cap sweep hook rejects invalid values without
# publishing a partial scientific output.
ln -s "$bgzf_fixture" "$work/invalid-cap-noindex.vcf.gz"
if VCFTOOLS_NG_TEST_COMPRESSED_COMPUTE_CAP=0 \
    "$ng" --gzvcf "$work/invalid-cap-noindex.vcf.gz" \
    --threads 8 --counts \
    --out "$work/invalid-compressed-cap" \
    >"$work/invalid-compressed-cap.stdout" \
    2>"$work/invalid-compressed-cap.stderr"; then
    echo 'Invalid compressed compute cap unexpectedly succeeded' >&2
    exit 1
fi
grep -Fq \
    'VCFTOOLS_NG_TEST_COMPRESSED_COMPUTE_CAP must be a positive integer' \
    "$work/invalid-compressed-cap.stderr"
test ! -e "$work/invalid-compressed-cap.frq.count"

# Lock Original VCFtools 0.1.17's observable 32-bit site-depth arithmetic
# and degenerate NaN formatting in both the fused and generic implementations.
"$ng" --vcf "$numeric_fixture" --threads 1 \
    --site-depth --site-mean-depth --out "$work/numeric-fast" \
    >"$work/numeric-fast.stdout" 2>"$work/numeric-fast.stderr"
cmp "$numeric_golden.ldepth" "$work/numeric-fast.ldepth"
cmp "$numeric_golden.ldepth.mean" "$work/numeric-fast.ldepth.mean"

"$ng" --vcf "$numeric_fixture" --threads 3 \
    --site-depth --site-mean-depth --depth --out "$work/numeric-generic" \
    >"$work/numeric-generic.stdout" 2>"$work/numeric-generic.stderr"
cmp "$numeric_golden.ldepth" "$work/numeric-generic.ldepth"
cmp "$numeric_golden.ldepth.mean" "$work/numeric-generic.ldepth.mean"
grep -Fq 'Input backend: plain-ranges' "$work/numeric-generic.stderr"

# The corrected-depth extension is explicit and separately locked. It must be
# identical in fused and generic implementations without changing exact mode.
for mode in fast generic; do
    if [[ $mode == fast ]]; then
        numeric_args=(--threads 1)
    else
        numeric_args=(--threads 3 --depth)
    fi
    "$ng" --vcf "$numeric_fixture" "${numeric_args[@]}" \
        --site-depth --site-mean-depth --corrected-depth-arithmetic \
        --out "$work/numeric-corrected-$mode" \
        >"$work/numeric-corrected-$mode.stdout" \
        2>"$work/numeric-corrected-$mode.stderr"
    cmp "$numeric_golden.corrected.ldepth" \
        "$work/numeric-corrected-$mode.ldepth"
    cmp "$numeric_golden.corrected.ldepth.mean" \
        "$work/numeric-corrected-$mode.ldepth.mean"
    grep -Fq -- '--corrected-depth-arithmetic' \
        "$work/numeric-corrected-$mode.stderr"
done

if "$ng" --vcf "$numeric_fixture" --threads 1 \
    --corrected-depth-arithmetic --counts \
    --out "$work/numeric-corrected-invalid" \
    >"$work/numeric-corrected-invalid.stdout" \
    2>"$work/numeric-corrected-invalid.stderr"; then
    echo 'Corrected depth without a site-depth output unexpectedly succeeded' \
        >&2
    exit 1
fi
grep -Fq 'requires --site-depth or --site-mean-depth' \
    "$work/numeric-corrected-invalid.stderr"

# Lock zero-called, all-missing, single-population-observation, empty-window,
# and NaN behavior across several outputs sharing one decoded scan.
for threads in 1 3; do
    prefix="$work/numeric-analyses-t$threads"
    "$ng" --vcf "$numeric_fixture" --threads "$threads" \
        --counts --missing-site --site-pi \
        --window-pi 10 --window-pi-step 5 --TajimaD 10 \
        --weir-fst-pop "$numeric_pop1" \
        --weir-fst-pop "$numeric_pop2" \
        --out "$prefix" >"$prefix.stdout" 2>"$prefix.stderr"
    cmp "$numeric_golden.frq.count" "$prefix.frq.count"
    cmp "$numeric_golden.lmiss" "$prefix.lmiss"
    cmp "$numeric_golden.sites.pi" "$prefix.sites.pi"
    cmp "$numeric_golden.windowed.pi" "$prefix.windowed.pi"
    cmp "$numeric_golden.Tajima.D" "$prefix.Tajima.D"
    cmp "$numeric_golden.weir.fst" "$prefix.weir.fst"
done

# Cross both inline-to-heap allele cutovers, including GT indices 15-17.
for threads in 1 3; do
    prefix="$work/multiallelic-t$threads"
    "$ng" --vcf "$multiallelic_fixture" --threads "$threads" \
        --counts --missing-site --site-pi --out "$prefix" \
        >"$prefix.stdout" 2>"$prefix.stderr"
    cmp "$multiallelic_golden.frq.count" "$prefix.frq.count"
    cmp "$multiallelic_golden.lmiss" "$prefix.lmiss"
    cmp "$multiallelic_golden.sites.pi" "$prefix.sites.pi"
done

# Site-only recode keeps a mid-stream triploid record opaque, while a semantic
# command fails transactionally after earlier diploid/haploid work.
for threads in 1 3; do
    prefix="$work/mixed-ploidy-recode-t$threads"
    "$ng" --vcf "$mixed_ploidy_fixture" --threads "$threads" \
        --minQ 0 --recode-vcf --recode-INFO-all --out "$prefix" \
        >"$prefix.stdout" 2>"$prefix.stderr"
    cmp "$mixed_ploidy_fixture" "$prefix.recode.vcf"
done
for threads in 1 3; do
    prefix="$work/mixed-ploidy-counts-t$threads"
    if "$ng" --vcf "$mixed_ploidy_fixture" --threads "$threads" \
        --counts --out "$prefix" \
        >"$prefix.stdout" 2>"$prefix.stderr"; then
        echo 'Mid-stream triploid semantic command unexpectedly succeeded' \
            >&2
        exit 1
    fi
    grep -Fq 'Polyploid genotype is not supported' "$prefix.stderr"
    test ! -e "$prefix.frq.count"
done

# Ordered window analyses reject decreasing positions/reappearing chromosome
# segments instead of emitting implementation-dependent plausible results.
# Site-local commands continue to accept physical record order.
"$ng" --vcf "$unsorted_fixture" --threads 3 --counts \
    --out "$work/unsorted-site-local" \
    >"$work/unsorted-site-local.stdout" \
    2>"$work/unsorted-site-local.stderr"
test -s "$work/unsorted-site-local.frq.count"
for threads in 1 3; do
    prefix="$work/unsorted-ordered-t$threads"
    if "$ng" --vcf "$unsorted_fixture" --threads "$threads" \
        --window-pi 250 --window-pi-step 100 --TajimaD 250 \
        --out "$prefix" >"$prefix.stdout" 2>"$prefix.stderr"; then
        echo 'Unsorted ordered analysis unexpectedly succeeded' >&2
        exit 1
    fi
    grep -Fq 'Ordered analysis requires a position-sorted input' \
        "$prefix.stderr"
    test ! -e "$prefix.windowed.pi"
    test ! -e "$prefix.Tajima.D"
done

# Population membership edges retain Original behavior: duplicate/unknown
# names are ignored, overlap is allowed, and an empty effective population
# yields Original-compatible NaN rows.
printf 'S1\nS1\nUNKNOWN\n' >"$work/pop-dup-1.txt"
printf 'S2\n' >"$work/pop-dup-2.txt"
printf 'S1\nS2\n' >"$work/pop-overlap-1.txt"
printf 'S2\n' >"$work/pop-overlap-2.txt"
printf 'UNKNOWN\n' >"$work/pop-empty-1.txt"
printf 'S2\n' >"$work/pop-empty-2.txt"
for threads in 1 3; do
    for population_case in dup overlap empty; do
        prefix="$work/pop-$population_case-t$threads"
        "$ng" --vcf "$numeric_fixture" --threads "$threads" \
            --weir-fst-pop "$work/pop-$population_case-1.txt" \
            --weir-fst-pop "$work/pop-$population_case-2.txt" \
            --out "$prefix" >"$prefix.stdout" 2>"$prefix.stderr"
        if [[ $population_case == overlap ]]; then
            cmp "$numeric_golden.overlap.weir.fst" "$prefix.weir.fst"
        else
            cmp "$numeric_golden.weir.fst" "$prefix.weir.fst"
        fi
    done
done

# HTSlib's parallel VCF parser cannot assign stable RIDs when ##contig lines
# are absent. Automatic mode must use the ordered stream rather than creating
# records with RID=-1; an explicitly forced plain-range backend fails clearly.
"$ng" --vcf "$no_contig_fixture" --threads 3 --depth --site-depth \
    --out "$work/no-contig-auto" \
    >"$work/no-contig-auto.stdout" 2>"$work/no-contig-auto.stderr"
test -s "$work/no-contig-auto.ldepth"
grep -Fq 'Input backend: stream' "$work/no-contig-auto.stderr"
grep -Fq 'without declared contigs requires an ordered HTSlib stream' \
    "$work/no-contig-auto.stderr"

if "$ng" --vcf "$no_contig_fixture" --threads 3 --input-backend plain \
    --depth --site-depth --out "$work/no-contig-forced" \
    >"$work/no-contig-forced.stdout" 2>"$work/no-contig-forced.stderr"; then
    echo 'Forced plain ranges without ##contig unexpectedly succeeded' >&2
    exit 1
fi
grep -Fq 'requires ##contig declarations' \
    "$work/no-contig-forced.stderr"
test ! -e "$work/no-contig-forced.ldepth"

# A truncated BGZF stream is an input failure, not a successful short file.
# No staged scientific output may be published.
cp -- "$bgzf_fixture" "$work/truncated.vcf.gz"
bgzf_size=$(stat -c %s "$work/truncated.vcf.gz")
truncate -s "$((bgzf_size - 4096))" "$work/truncated.vcf.gz"
if "$ng" --gzvcf "$work/truncated.vcf.gz" --threads 4 --counts \
    --out "$work/truncated" \
    >"$work/truncated.stdout" 2>"$work/truncated.stderr"; then
    echo 'Truncated BGZF input unexpectedly succeeded' >&2
    exit 1
fi
grep -Fq 'HTSlib failed while reading compressed VCF records' \
    "$work/truncated.stderr"
test ! -e "$work/truncated.frq.count"
test -z "$(find "$work" -name '*.vcftools-ng.tmp.*' -print -quit)"

# Reader and worker exceptions must cancel every wait, join all threads, and
# leave no published or staged scientific artifact. timeout turns a deadlock
# into a deterministic regression failure.
for injection in \
    VCFTOOLS_NG_TEST_FAIL_READER_AFTER_BATCHES=1 \
    VCFTOOLS_NG_TEST_FAIL_WORKER_AFTER_RECORDS=1; do
    prefix="$work/${injection%%=*}"
    set +e
    timeout 10s env "$injection" \
        "$ng" --vcf "$fixture" --threads 3 --depth --site-depth \
        --batch-size 2 --out "$prefix" \
        >"$prefix.stdout" 2>"$prefix.stderr"
    status=$?
    set -e
    test "$status" -eq 1
    grep -Fq 'Injected ordered-pipeline' "$prefix.stderr"
    grep -Fq 'Exit status: failed' "$prefix.stderr"
    test ! -e "$prefix.ldepth"
    test ! -e "$prefix.idepth"
done
test -z "$(find "$work" -name '*.vcftools-ng.tmp.*' -print -quit)"

# Compression and ordered-writer failures must release every BGZF waiter and
# must never publish a partial compressed VCF.
for injection in \
    VCFTOOLS_NG_TEST_FAIL_COMPRESSOR_AFTER_JOBS=1 \
    VCFTOOLS_NG_TEST_FAIL_WRITER_AFTER_BLOCKS=1; do
    prefix="$work/${injection%%=*}"
    set +e
    timeout 10s env "$injection" \
        "$ng" --vcf "$fixture" --threads 3 --recode \
        --out "$prefix" >"$prefix.stdout" 2>"$prefix.stderr"
    status=$?
    set -e
    test "$status" -eq 1
    grep -Eq 'Injected BGZF (compressor|writer) failure' "$prefix.stderr"
    grep -Fq 'Exit status: failed' "$prefix.stderr"
    test ! -e "$prefix.recode.vcf.gz"
done
test -z "$(find "$work" -name '*.vcftools-ng.tmp.*' -print -quit)"

# A record larger than both the plain reader block and the historical 8 MiB
# compressed batch target must remain one record. The final input line
# intentionally has no newline. Compare serial and aligned-range text kernels.
long_vcf="$work/long-record.vcf"
{
    printf '%s\n' '##fileformat=VCFv4.2'
    printf '%s\n' '##contig=<ID=chr1,length=1000>'
    printf '%s\n' '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
    printf 'chr1\t1\t.\tA\t'
    head -c 9000000 /dev/zero | tr '\0' A
    printf '\t30\tPASS\t.\tGT\t0/1'
} >"$long_vcf"
for threads in 1 3; do
    "$ng" --vcf "$long_vcf" --threads "$threads" --counts \
        --out "$work/long-t$threads" \
        >"$work/long-t$threads.stdout" 2>"$work/long-t$threads.stderr"
done
cmp "$work/long-t1.frq.count" "$work/long-t3.frq.count"
grep -Fq 'Input backend: fast-counts-plain' "$work/long-t3.stderr"

# CRLF is accepted consistently by serial and parallel text implementations.
printf '##fileformat=VCFv4.2\r\n##contig=<ID=chr1,length=10>\r\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\r\nchr1\t1\t.\tA\tG\t30\tPASS\t.\r\n' \
    >"$work/crlf.vcf"
for threads in 1 3; do
    "$ng" --vcf "$work/crlf.vcf" --threads "$threads" --site-quality \
        --out "$work/crlf-t$threads" \
        >"$work/crlf-t$threads.stdout" 2>"$work/crlf-t$threads.stderr"
done
cmp "$work/crlf-t1.lqual" "$work/crlf-t3.lqual"

# Exercise non-power-of-two thread budgets and ordered publication over enough
# Plain VCF data to create multiple aligned shards. This is a determinism gate,
# not a timing benchmark.
awk 'BEGIN {
    print "##fileformat=VCFv4.2"
    print "##contig=<ID=chr1,length=100000>"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    for (i = 1; i <= 50000; ++i) {
        print "chr1\t" i "\t.\tA\tG\t30\tPASS\t."
    }
}' >"$work/odd-threads.vcf"
"$ng" --vcf "$work/odd-threads.vcf" --threads 1 \
    --counts --site-quality --out "$work/odd-reference" \
    >"$work/odd-reference.stdout" 2>"$work/odd-reference.stderr"
for threads in 2 3 5 7 11 31; do
    prefix="$work/odd-t$threads"
    VCFTOOLS_NG_TEST_AVAILABLE_THREADS="$threads" \
        "$ng" --vcf "$work/odd-threads.vcf" --threads "$threads" \
        --counts --site-quality --out "$prefix" \
        >"$prefix.stdout" 2>"$prefix.stderr"
    cmp "$work/odd-reference.frq.count" "$prefix.frq.count"
    cmp "$work/odd-reference.lqual" "$prefix.lqual"
done

echo RELIABILITY_REGRESSION_PASS
