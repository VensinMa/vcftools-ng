#!/bin/sh
set -eu

package_root=${VCFTOOLS_NG_BENCHMARK_PACKAGE_ROOT:-/home/vensin/software/vcftools-ng-v0.13.0-linux-x86_64}
private_bcftools=${package_root}/libexec/bcftools.bin
private_libdir=${package_root}/lib

if [ ! -x "$private_bcftools" ] || [ ! -d "$private_libdir" ]; then
    echo "error: validated portable bcftools package is unavailable: $package_root" >&2
    exit 127
fi

LD_LIBRARY_PATH="${private_libdir}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export LD_LIBRARY_PATH
exec "$private_bcftools" "$@"
