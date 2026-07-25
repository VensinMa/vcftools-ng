#!/bin/sh
set -eu

package_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
PATH="${package_root}/bin${PATH:+:${PATH}}"
LD_LIBRARY_PATH="${package_root}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export PATH LD_LIBRARY_PATH

exec "${package_root}/libexec/vcftools-ng.bin" "$@"

