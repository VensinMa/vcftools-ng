#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    printf 'Usage: %s /path/to/vcftools-ng /path/to/source-root\n' "$0" >&2
    exit 2
fi

ng=$1
source_root=$2
fixture="$source_root/tests/fixtures/osmanthus412.flags.23chr_1k.vcf.gz"
bcf_fixture="$source_root/tests/fixtures/osmanthus412.flags.23chr_1k.bcf"
bcftools_bin=${BCFTOOLS:-}
if [[ -z "$bcftools_bin" ]]; then
    bcftools_bin=$(command -v bcftools || true)
fi
if [[ -z "$bcftools_bin" ]]; then
    printf 'bcftools is required for the auto-index test\n' >&2
    exit 1
fi
work=$(mktemp -d "${TMPDIR:-/tmp}/vcftools-ng-input-test.XXXXXX")
fixture_tbi_sha256=$(sha256sum "$fixture.tbi" | cut -d' ' -f1)
fixture_csi_sha256=$(sha256sum "$bcf_fixture.csi" | cut -d' ' -f1)
cleanup() {
    rm -rf -- "$work"
}
trap cleanup EXIT INT TERM

gzip -dc "$fixture" >"$work/input.vcf"

"$ng" --gzvcf "$fixture" --input-backend stream --threads 4 \
    --counts --out "$work/stream" >/dev/null 2>"$work/stream.log"
"$ng" --gzvcf "$fixture" --input-backend indexed --threads 4 \
    --counts --out "$work/indexed-vcf" \
    >/dev/null 2>"$work/indexed-vcf.log"
"$ng" --bcf "$bcf_fixture" --input-backend indexed --threads 4 \
    --counts --out "$work/indexed-bcf" \
    >/dev/null 2>"$work/indexed-bcf.log"
"$ng" --vcf "$work/input.vcf" --input-backend plain --threads 4 \
    --counts --out "$work/plain" >/dev/null 2>"$work/plain.log"

cp "$fixture" "$work/auto-index.vcf.gz"
"$ng" --gzvcf "$work/auto-index.vcf.gz" --threads 4 \
    --bcftools "$bcftools_bin" \
    --freq --counts --out "$work/auto-index" \
    >/dev/null 2>"$work/auto-index.log"

cp "$bcf_fixture" "$work/auto-index.bcf"
"$ng" --bcf "$work/auto-index.bcf" --threads 4 \
    --bcftools "$bcftools_bin" \
    --freq --counts --out "$work/auto-index-bcf" \
    >/dev/null 2>"$work/auto-index-bcf.log"

cp "$fixture" "$work/no-auto-index.vcf.gz"
"$ng" --gzvcf "$work/no-auto-index.vcf.gz" --threads 4 \
    --no-auto-index --bcftools "$bcftools_bin" \
    --freq --counts --out "$work/no-auto-index" \
    >/dev/null 2>"$work/no-auto-index.log"

cp "$fixture" "$work/auto-threads.vcf.gz"
SLURM_CPUS_PER_TASK=3 \
    "$ng" --gzvcf "$work/auto-threads.vcf.gz" \
    --bcftools "$bcftools_bin" \
    --freq --counts --out "$work/auto-threads" \
    >/dev/null 2>"$work/auto-threads.log"

cp "$fixture" "$work/concurrent.vcf.gz"
"$ng" --gzvcf "$work/concurrent.vcf.gz" --threads 4 \
    --bcftools "$bcftools_bin" \
    --freq --counts --out "$work/concurrent-a" \
    >/dev/null 2>"$work/concurrent-a.log" &
first_pid=$!
"$ng" --gzvcf "$work/concurrent.vcf.gz" --threads 4 \
    --bcftools "$bcftools_bin" \
    --freq --counts --out "$work/concurrent-b" \
    >/dev/null 2>"$work/concurrent-b.log" &
second_pid=$!
wait "$first_pid"
wait "$second_pid"

cp "$fixture" "$work/missing-bcftools.vcf.gz"
"$ng" --gzvcf "$work/missing-bcftools.vcf.gz" --threads 4 \
    --bcftools "$work/does-not-exist" \
    --freq --counts --out "$work/missing-bcftools" \
    >/dev/null 2>"$work/missing-bcftools.log"
if "$ng" --gzvcf "$work/missing-bcftools.vcf.gz" \
    --input-backend indexed --threads 4 \
    --bcftools "$work/does-not-exist" \
    --freq --counts --out "$work/missing-bcftools-explicit" \
    >/dev/null 2>"$work/missing-bcftools-explicit.log"; then
    printf 'Explicit indexed backend unexpectedly accepted no index\n' >&2
    exit 1
fi

cp "$fixture" "$work/corrupt-index.vcf.gz"
cp "$fixture" "$work/corrupt-index.vcf.gz.csi"
corrupt_index_sha256=$(
    sha256sum "$work/corrupt-index.vcf.gz.csi" | cut -d' ' -f1
)
"$ng" --gzvcf "$work/corrupt-index.vcf.gz" --threads 4 \
    --bcftools "$bcftools_bin" \
    --freq --counts --out "$work/corrupt-index" \
    >/dev/null 2>"$work/corrupt-index.log"

cp "$fixture" "$work/valid-tbi-corrupt-csi.vcf.gz"
cp "$fixture.tbi" "$work/valid-tbi-corrupt-csi.vcf.gz.tbi"
cp "$fixture" "$work/valid-tbi-corrupt-csi.vcf.gz.csi"
dual_bad_csi_sha256=$(
    sha256sum "$work/valid-tbi-corrupt-csi.vcf.gz.csi" |
        cut -d' ' -f1
)
"$ng" --gzvcf "$work/valid-tbi-corrupt-csi.vcf.gz" --threads 4 \
    --bcftools "$bcftools_bin" \
    --freq --counts --out "$work/valid-tbi-corrupt-csi" \
    >/dev/null 2>"$work/valid-tbi-corrupt-csi.log"

cp "$fixture" "$work/stale-index.vcf.gz"
cp "$fixture.tbi" "$work/stale-index.vcf.gz.tbi"
touch -d '@2000000000' "$work/stale-index.vcf.gz"
stale_index_sha256=$(
    sha256sum "$work/stale-index.vcf.gz.tbi" | cut -d' ' -f1
)
"$ng" --gzvcf "$work/stale-index.vcf.gz" --threads 4 \
    --bcftools "$bcftools_bin" \
    --freq --counts --out "$work/stale-index" \
    >/dev/null 2>"$work/stale-index.log"

for candidate in indexed-vcf indexed-bcf plain; do
    cmp "$work/stream.frq.count" "$work/$candidate.frq.count"
done
for candidate in auto-index auto-index-bcf no-auto-index auto-threads \
                 concurrent-a concurrent-b missing-bcftools \
                 corrupt-index valid-tbi-corrupt-csi stale-index; do
    cmp "$work/stream.frq.count" "$work/$candidate.frq.count"
done

grep -q 'Input backend: stream' "$work/stream.log"
grep -q 'Input backend: indexed-regions' "$work/indexed-vcf.log"
grep -q 'Input backend: indexed-regions' "$work/indexed-bcf.log"
grep -q 'Input backend: plain-ranges' "$work/plain.log"
test ! -e "$work/input.vcf.csi"
test ! -e "$work/input.vcf.tbi"
if grep -q 'Auto-index:' "$work/plain.log"; then
    printf 'Plain VCF unexpectedly invoked automatic indexing\n' >&2
    exit 1
fi
grep -q 'Auto-index: no CSI/TBI sidecar found' "$work/auto-index.log"
grep -q 'index --csi --threads 4' "$work/auto-index.log"
grep -q 'Input backend: indexed-regions' "$work/auto-index.log"
test -s "$work/auto-index.vcf.gz.csi"

grep -q 'index --csi --threads 4' "$work/auto-index-bcf.log"
grep -q 'Input backend: indexed-regions' "$work/auto-index-bcf.log"
test -s "$work/auto-index.bcf.csi"

grep -q 'Input backend: stream' "$work/no-auto-index.log"
test ! -e "$work/no-auto-index.vcf.gz.csi"

grep -q 'index --csi --threads 3' "$work/auto-threads.log"
grep -q 'Threads: 3 (auto from SLURM_CPUS_PER_TASK)' \
    "$work/auto-threads.log"
grep -q 'Input backend: indexed-regions' "$work/auto-threads.log"
test -s "$work/auto-threads.vcf.gz.csi"

test -s "$work/concurrent.vcf.gz.csi"
grep -q 'Input backend: indexed-regions' "$work/concurrent-a.log"
grep -q 'Input backend: indexed-regions' "$work/concurrent-b.log"
test -z "$(find "$work" -name '*.vcftools-ng.tmp.*' -print -quit)"

grep -q 'Auto-index warning: automatic CSI construction failed' \
    "$work/missing-bcftools.log"
grep -q 'Input backend: stream' "$work/missing-bcftools.log"
grep -q 'Error: automatic CSI construction failed' \
    "$work/missing-bcftools-explicit.log"
test ! -e "$work/missing-bcftools.vcf.gz.csi"

grep -q 'protected sidecar is unusable' "$work/corrupt-index.log"
grep -q 'refusing to overwrite it' "$work/corrupt-index.log"
grep -q 'Input backend: stream' "$work/corrupt-index.log"
test "$corrupt_index_sha256" = "$(
    sha256sum "$work/corrupt-index.vcf.gz.csi" | cut -d' ' -f1
)"

grep -q 'protected sidecar is unusable' \
    "$work/valid-tbi-corrupt-csi.log"
grep -q 'Input backend: indexed-regions' \
    "$work/valid-tbi-corrupt-csi.log"
grep -q "via $work/valid-tbi-corrupt-csi.vcf.gz.tbi" \
    "$work/valid-tbi-corrupt-csi.log"
test "$dual_bad_csi_sha256" = "$(
    sha256sum "$work/valid-tbi-corrupt-csi.vcf.gz.csi" |
        cut -d' ' -f1
)"
test "$fixture_tbi_sha256" = "$(
    sha256sum "$work/valid-tbi-corrupt-csi.vcf.gz.tbi" |
        cut -d' ' -f1
)"

grep -q 'index is older than the data file' "$work/stale-index.log"
grep -q 'refusing to overwrite it' "$work/stale-index.log"
grep -q 'Input backend: stream' "$work/stale-index.log"
test "$stale_index_sha256" = "$(
    sha256sum "$work/stale-index.vcf.gz.tbi" | cut -d' ' -f1
)"

test "$fixture_tbi_sha256" = "$(
    sha256sum "$fixture.tbi" | cut -d' ' -f1
)"
test "$fixture_csi_sha256" = "$(
    sha256sum "$bcf_fixture.csi" | cut -d' ' -f1
)"

printf 'Input backend differential test passed\n'
