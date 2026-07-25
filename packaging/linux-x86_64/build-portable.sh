#!/usr/bin/env bash
set -euo pipefail

repository_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
container_image=${VCFTOOLS_NG_BUILD_IMAGE:-quay.io/pypa/manylinux2014_x86_64:latest}
output_directory=${VCFTOOLS_NG_OUTPUT_DIR:-"${repository_root}/dist"}

command -v docker >/dev/null 2>&1 || {
    echo "error: docker is required to build the portable archive" >&2
    exit 1
}

mkdir -p "${output_directory}"

docker run --rm \
    -e "HOST_UID=$(id -u)" \
    -e "HOST_GID=$(id -g)" \
    -v "${repository_root}:/source:ro" \
    -v "${output_directory}:/output" \
    "${container_image}" \
    /source/packaging/linux-x86_64/build-in-container.sh

