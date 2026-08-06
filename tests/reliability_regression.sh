#!/usr/bin/env bash
set -euo pipefail

ng=${1:?vcftools-ng executable required}
source_root=${2:?source root required}
fixture="$source_root/tests/fixtures/fast-site-stats.vcf"
work=$(mktemp -d "${TMPDIR:-/tmp}/vcftools-ng-reliability.XXXXXX")
cleanup() {
    rm -rf -- "$work"
}
trap cleanup EXIT INT TERM

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

echo RELIABILITY_REGRESSION_PASS
