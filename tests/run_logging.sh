#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 /path/to/vcftools-ng /path/to/repository" >&2
    exit 2
fi

ng=$1
repository=$2
fixture="$repository/tests/fixtures/fast-site-stats.vcf"
work=$(mktemp -d)
trap 'rm -rf -- "$work"' EXIT
bcftools_bin=${BCFTOOLS:-bcftools}

"$bcftools_bin" view -Oz -o "$work/indexed.vcf.gz" "$fixture"
"$bcftools_bin" index --tbi "$work/indexed.vcf.gz"
"$bcftools_bin" view -Ob -o "$work/indexed.bcf" "$fixture"
"$bcftools_bin" index --csi "$work/indexed.bcf"

assert_standard_log() {
    local log=$1
    grep -Eq '^vcftools-ng [0-9]+\.[0-9]+\.[0-9]+$' "$log"
    grep -Fqx 'Log format: vcftools-ng-text-v1' "$log"
    grep -Eq '^Start time: .+[+-][0-9]{2}:[0-9]{2}$' "$log"
    grep -Fq 'Command: ' "$log"
    grep -Fq 'Working directory: ' "$log"
    grep -Fqx 'Input format: Plain VCF' "$log"
    grep -Eq '^Input size: [0-9]+ bytes$' "$log"
    grep -Fqx 'Index policy: automatic' "$log"
    grep -Fqx 'Existing sidecar: none' "$log"
    grep -Fqx 'Index validation: not applicable' "$log"
    grep -Fqx 'Execution kernel: fused-text-filter-recode' "$log"
    grep -Fqx 'Execution components: vcf-recode' "$log"
    grep -Fqx 'Fused text fast path: selected' "$log"
    grep -Fqx 'Selected backend: fast-filter-recode-plain' "$log"
    grep -Fqx 'Index used: no' "$log"
    grep -Fqx 'Requested threads: 1' "$log"
    grep -Fqx 'Effective threads: 1' "$log"
    grep -Fqx 'Outputs:' "$log"
    grep -Fqx 'Filters: none' "$log"
    grep -Fqx 'Exit status: success' "$log"
    grep -Eq '^End time: .+[+-][0-9]{2}:[0-9]{2}$' "$log"
    grep -Eq '^Wall time: [0-9.]+ seconds$' "$log"
    grep -Eq '^User CPU time: [0-9.]+ seconds$' "$log"
    grep -Eq '^System CPU time: [0-9.]+ seconds$' "$log"
    grep -Eq '^Average CPU: [0-9.]+%$' "$log"
    grep -Eq '^Peak RSS: [0-9]+ KiB$' "$log"
    grep -Eq \
        '^Stage time \[fused scan/filter/output\]: [0-9.]+ seconds$' \
        "$log"
    grep -Eq \
        '^Stage time \[output finalization\]: [0-9.]+ seconds$' \
        "$log"
}

"$ng" \
    --vcf "$fixture" --threads 1 --recode \
    --out "$work/default" \
    >"$work/default.stdout" 2>"$work/default.stderr"
test -s "$work/default.log"
test -s "$work/default.recode.vcf.gz"
test ! -e "$work/default.recode.vcf"
assert_standard_log "$work/default.log"
cmp "$work/default.stderr" "$work/default.log"
grep -Fqx "  Recode BGZF VCF: $work/default.recode.vcf.gz" \
    "$work/default.log"
grep -Eq '^  Output bytes: [0-9]+$' "$work/default.log"

printf 'stale log contents\n' >"$work/custom.log"
"$ng" \
    --vcf "$fixture" --threads 1 --counts \
    --log-file "$work/custom.log" \
    --out "$work/custom-prefix" \
    >"$work/custom.stdout" 2>"$work/custom.stderr"
test -s "$work/custom.log"
test ! -e "$work/custom-prefix.log"
if grep -Fq 'stale log contents' "$work/custom.log"; then
    echo "Custom log was not overwritten" >&2
    exit 1
fi
cmp "$work/custom.stderr" "$work/custom.log"
grep -Fqx "Log file: $work/custom.log" "$work/custom.log"
grep -Fqx 'Decision: CSI/TBI is not applicable to Plain VCF' \
    "$work/custom.log"
grep -Fqx 'Selected backend: fast-counts-plain' "$work/custom.log"
grep -Fqx 'Execution kernel: fused-text-site-statistics' \
    "$work/custom.log"
grep -Fqx 'Execution components: site-statistics' "$work/custom.log"

"$ng" \
    --gzvcf "$work/indexed.vcf.gz" --threads 2 --recode \
    --out "$work/existing-tbi" \
    >"$work/existing-tbi.stdout" 2>"$work/existing-tbi.stderr"
grep -Fqx "Existing sidecar: $work/indexed.vcf.gz.tbi" \
    "$work/existing-tbi.log"
grep -Fqx 'Index type: TBI' "$work/existing-tbi.log"
grep -Fqx 'Decision: use existing index' "$work/existing-tbi.log"
grep -Fqx 'Selected backend: fast-filter-recode-indexed-bgzf' \
    "$work/existing-tbi.log"
grep -Fqx 'Index used: yes' "$work/existing-tbi.log"

"$ng" \
    --bcf "$work/indexed.bcf" --threads 2 --recode \
    --out "$work/bcf-stream" \
    >"$work/bcf-stream.stdout" 2>"$work/bcf-stream.stderr"
grep -Fqx "Existing sidecar: $work/indexed.bcf.csi" \
    "$work/bcf-stream.log"
grep -Fqx 'Index type: CSI' "$work/bcf-stream.log"
grep -Fqx \
    'Decision: preserve existing sidecars; do not build or use an index' \
    "$work/bcf-stream.log"
grep -Fqx 'Selected backend: stream' "$work/bcf-stream.log"
grep -Fqx 'Index used: no' "$work/bcf-stream.log"
grep -Fqx 'Execution kernel: generic-ordered-pipeline' \
    "$work/bcf-stream.log"
grep -Fqx 'Execution components: vcf-recode' "$work/bcf-stream.log"
grep -Fqx 'Fused text fast path: not selected' "$work/bcf-stream.log"
grep -Eq '^Stage time \[input/index planning\]: [0-9.]+ seconds$' \
    "$work/bcf-stream.log"
grep -Eq \
    '^Stage time \[ordered input/compute/commit\]: [0-9.]+ seconds$' \
    "$work/bcf-stream.log"
test -s "$work/indexed.bcf.csi"

cp "$work/indexed.vcf.gz" "$work/auto.vcf.gz"
"$ng" \
    --gzvcf "$work/auto.vcf.gz" --threads 2 --recode \
    --out "$work/auto-csi" \
    >"$work/auto-csi.stdout" 2>"$work/auto-csi.stderr"
test -s "$work/auto.vcf.gz.csi"
grep -Fqx 'Existing sidecar: none' "$work/auto-csi.log"
grep -Fqx 'Decision: build CSI' "$work/auto-csi.log"
grep -Fqx 'Index build threads: 2' "$work/auto-csi.log"
grep -Eq '^Index build time: [0-9.]+ seconds$' \
    "$work/auto-csi.log"
grep -Fqx 'Index build result: PASS' "$work/auto-csi.log"
grep -Fqx 'Index type: CSI' "$work/auto-csi.log"
grep -Fqx 'Selected backend: fast-filter-recode-indexed-bgzf' \
    "$work/auto-csi.log"
grep -Fqx 'Index used: yes' "$work/auto-csi.log"

"$ng" \
    --vcf "$fixture" --threads 1 --counts \
    --out "$work/same-path" \
    >"$work/same-path.stdout" 2>"$work/same-path.log"
grep -Fqx 'Exit status: success' "$work/same-path.log"
test "$(
    grep -Ec '^vcftools-ng [0-9]+\.[0-9]+\.[0-9]+$' \
        "$work/same-path.log"
)" -eq 1

"$ng" \
    --vcf "$fixture" --threads 1 --counts \
    --no-log-file --out "$work/no-log" \
    >"$work/no-log.stdout" 2>"$work/no-log.stderr"
test ! -e "$work/no-log.log"
grep -Fqx 'Log file: disabled (--no-log-file)' \
    "$work/no-log.stderr"
grep -Fqx 'Exit status: success' "$work/no-log.stderr"

"$ng" \
    --vcf "$fixture" --threads 1 --recode --stdout \
    --out "$work/stdout" \
    >"$work/stdout.vcf" 2>"$work/stdout.stderr"
test -s "$work/stdout.log"
cmp "$work/stdout.stderr" "$work/stdout.log"
grep -Fqx '  Recode VCF: stdout' "$work/stdout.log"
grep -Fq $'#CHROM\tPOS\tID\tREF\tALT' "$work/stdout.vcf"
if grep -Fq 'Log format:' "$work/stdout.vcf"; then
    echo "Log text contaminated VCF stdout" >&2
    exit 1
fi

if "$ng" \
    --vcf "$fixture" --threads 1 --recode \
    --input-backend indexed --out "$work/failed" \
    >"$work/failed.stdout" 2>"$work/failed.stderr"; then
    echo "Expected forced indexed Plain VCF run to fail" >&2
    exit 1
fi
test -s "$work/failed.log"
cmp "$work/failed.stderr" "$work/failed.log"
grep -Fq 'Error: --input-backend indexed requires' "$work/failed.log"
grep -Fqx 'Exit status: failed' "$work/failed.log"

if "$ng" \
    --vcf "$fixture" --threads 1 --counts \
    --unsupported-log-test --out "$work/parse-failed" \
    >"$work/parse-failed.stdout" 2>"$work/parse-failed.stderr"; then
    echo "Expected parse failure" >&2
    exit 1
fi
test -s "$work/parse-failed.log"
cmp "$work/parse-failed.stderr" "$work/parse-failed.log"
grep -Fq 'Error: Unsupported option' "$work/parse-failed.log"
grep -Fqx 'Exit status: failed' "$work/parse-failed.log"

if "$ng" \
    --vcf "$fixture" --threads 1 --counts \
    --out "$work/conflict" \
    --log-file "$work/conflict.custom.log" --no-log-file \
    >"$work/conflict.stdout" 2>"$work/conflict.stderr"; then
    echo "Expected conflicting log options to fail" >&2
    exit 1
fi
grep -Fq -- '--log-file and --no-log-file cannot be combined' \
    "$work/conflict.stderr"

mkdir "$work/log-is-directory"
if "$ng" \
    --vcf "$fixture" --threads 1 --counts \
    --log-file "$work/log-is-directory" \
    --out "$work/unwritable" \
    >"$work/unwritable.stdout" 2>"$work/unwritable.stderr"; then
    echo "Expected log creation failure" >&2
    exit 1
fi
grep -Fq 'Could not create log file' "$work/unwritable.stderr"

input_hash_before=$(sha256sum "$fixture" | cut -d' ' -f1)
if "$ng" \
    --vcf "$fixture" --threads 1 --counts \
    --log-file "$fixture" --out "$work/input-conflict" \
    >"$work/input-conflict.stdout" \
    2>"$work/input-conflict.stderr"; then
    echo "Expected input/log path conflict to fail" >&2
    exit 1
fi
test "$(sha256sum "$fixture" | cut -d' ' -f1)" = "$input_hash_before"
grep -Fq 'Log file must not overwrite an input file' \
    "$work/input-conflict.stderr"

printf 'protected output sentinel\n' >"$work/output-conflict.recode.vcf.gz"
if "$ng" \
    --vcf "$fixture" --threads 1 --recode \
    --log-file "$work/output-conflict.recode.vcf.gz" \
    --out "$work/output-conflict" \
    >"$work/output-conflict.stdout" \
    2>"$work/output-conflict.stderr"; then
    echo "Expected output/log path conflict to fail" >&2
    exit 1
fi
grep -Fqx 'protected output sentinel' \
    "$work/output-conflict.recode.vcf.gz"
grep -Fq 'Log file conflicts with output file' \
    "$work/output-conflict.stderr"

(
    cd "$work"
    "$ng" --help >help.txt 2>help.stderr
    "$ng" --version >version.txt 2>version.stderr
    test ! -e out.log
    test ! -s help.stderr
    test ! -s version.stderr
)

echo "RUN_LOGGING_PASS"
