#!/usr/bin/env bash
set -euo pipefail

repository_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
archive=${1:-"$repository_root/dist/vcftools-ng-v0.12.3-linux-x86_64.tar.gz"}
checksum="${archive}.sha256"
archive_name=$(basename "$archive")
archive_directory=$(cd -- "$(dirname -- "$archive")" && pwd)

command -v docker >/dev/null 2>&1 || {
    echo "error: docker is required to test the portable archive" >&2
    exit 1
}
[[ -s "$archive" && -s "$checksum" ]] || {
    echo "error: archive or checksum is missing" >&2
    exit 1
}

(
    cd "$archive_directory"
    sha256sum -c "$(basename "$checksum")"
)

for image in centos:7 ubuntu:20.04; do
    docker run --rm -i \
        -e "ARCHIVE_NAME=$archive_name" \
        -v "$archive_directory:/dist:ro" \
        -v "$repository_root/tests/fixtures:/fixtures:ro" \
        -v "$repository_root/tests/golden:/golden:ro" \
        "$image" \
        bash -s <<'CONTAINER'
set -euo pipefail

mkdir -p /work
cd /work
tar -xzf "/dist/$ARCHIVE_NAME"
package_root="/work/${ARCHIVE_NAME%.tar.gz}"

"$package_root/bin/vcftools-ng" --version
"$package_root/bin/bcftools" --version | head -n 2
env -u NO_COLOR CLICOLOR_FORCE=1 \
    "$package_root/bin/vcftools-ng" --help > /work/help.color
NO_COLOR=1 CLICOLOR_FORCE=1 \
    "$package_root/bin/vcftools-ng" --help > /work/help.plain
grep -Fq $'\033[1;36m' /work/help.color
grep -Fq 'QUICK EXAMPLES:' /work/help.plain
if grep -q $'\033\\[' /work/help.plain; then
    echo "error: NO_COLOR help contains ANSI escapes" >&2
    exit 1
fi
if LD_LIBRARY_PATH="$package_root/lib" \
    ldd "$package_root/libexec/vcftools-ng.bin" |
    grep -q "not found"; then
    exit 1
fi

cp /fixtures/osmanthus412.flags.23chr_1k.vcf.gz /work/input.vcf.gz
"$package_root/bin/vcftools-ng" \
    --gzvcf /work/input.vcf.gz \
    --threads 4 \
    --keep-filtered q10 \
    --keep-filtered PASS \
    --remove-filtered Cluster \
    --keep-INFO Hotspot \
    --remove-INFO Artifact \
    --recode \
    --recode-INFO-all \
    --out /work/out

[[ -s /work/input.vcf.gz.csi ]]
cmp /golden/flags-site-info.recode.vcf /work/out.recode.vcf
printf 'PORTABLE_PASS %s\n' "$(
    . /etc/os-release
    printf '%s' "$PRETTY_NAME"
)"
CONTAINER
done
