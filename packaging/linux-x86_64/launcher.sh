#!/bin/sh
set -eu

package_root=$(CDPATH='' cd -- "$(dirname -- "$0")/.." && pwd)
LD_LIBRARY_PATH="${package_root}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
VCFTOOLS_NG_BCFTOOLS="${package_root}/libexec/bcftools.bin"
export LD_LIBRARY_PATH VCFTOOLS_NG_BCFTOOLS

exec "${package_root}/libexec/vcftools-ng.bin" "$@"
