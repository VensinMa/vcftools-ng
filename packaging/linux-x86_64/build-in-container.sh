#!/usr/bin/env bash
set -euo pipefail

readonly version=0.12.3
readonly htslib_version=1.23
readonly bcftools_version=1.23
readonly htslib_sha256=63927199ef9cea03096345b95d96cb600ae10385248b2ef670b0496c2ab7e4cd
readonly bcftools_sha256=5acde0ac38f7981da1b89d8851a1a425d1c275e1eb76581925c04ca4252c0778
readonly archive_name="vcftools-ng-v${version}-linux-x86_64"
readonly work_directory=/tmp/vcftools-ng-portable
readonly dependency_prefix="${work_directory}/prefix"
readonly package_root="${work_directory}/${archive_name}"

rm -rf "${work_directory}"
mkdir -p "${work_directory}/downloads" "${dependency_prefix}" \
    "${package_root}/bin" "${package_root}/lib" \
    "${package_root}/libexec" "${package_root}/LICENSES"

yum install -y \
    blas-devel \
    devtoolset-11-gcc-c++ \
    lapack-devel \
    zlib-devel \
    >/dev/null

source /opt/rh/devtoolset-11/enable
export CC=gcc
export CXX=g++

curl -fsSL \
    "https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2" \
    -o "${work_directory}/downloads/htslib.tar.bz2"
curl -fsSL \
    "https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2" \
    -o "${work_directory}/downloads/bcftools.tar.bz2"

echo "${htslib_sha256}  ${work_directory}/downloads/htslib.tar.bz2" \
    | sha256sum --check --status
echo "${bcftools_sha256}  ${work_directory}/downloads/bcftools.tar.bz2" \
    | sha256sum --check --status

tar -xjf "${work_directory}/downloads/htslib.tar.bz2" -C "${work_directory}"
tar -xjf "${work_directory}/downloads/bcftools.tar.bz2" -C "${work_directory}"

pushd "${work_directory}/htslib-${htslib_version}" >/dev/null
./configure \
    --prefix="${dependency_prefix}" \
    --disable-bz2 \
    --disable-lzma \
    >/dev/null
make -j"$(nproc)" >/dev/null
make install >/dev/null
popd >/dev/null

pushd "${work_directory}/bcftools-${bcftools_version}" >/dev/null
./configure \
    --prefix="${dependency_prefix}" \
    --with-htslib="${dependency_prefix}" \
    --disable-bcftools-plugins \
    >/dev/null
make -j"$(nproc)" >/dev/null
make install >/dev/null
popd >/dev/null

cmake \
    -S /source \
    -B "${work_directory}/build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=OFF \
    -DHTSLIB_ROOT="${dependency_prefix}" \
    >/dev/null
cmake --build "${work_directory}/build" --parallel "$(nproc)" >/dev/null

cp "${work_directory}/build/vcftools-ng" \
    "${package_root}/libexec/vcftools-ng.bin"
cp "${dependency_prefix}/bin/bcftools" \
    "${package_root}/libexec/bcftools.bin"
cp /source/packaging/linux-x86_64/launcher.sh \
    "${package_root}/bin/vcftools-ng"
cp /source/packaging/linux-x86_64/bcftools-launcher.sh \
    "${package_root}/bin/bcftools"
cp /source/packaging/linux-x86_64/README.portable.md \
    "${package_root}/README.md"
cp "${work_directory}/htslib-${htslib_version}/LICENSE" \
    "${package_root}/LICENSES/HTSlib-LICENSE"
cp "${work_directory}/bcftools-${bcftools_version}/LICENSE" \
    "${package_root}/LICENSES/BCFtools-LICENSE"

sed -i "s/@VERSION@/${version}/g" "${package_root}/README.md"
chmod 0755 "${package_root}/bin/vcftools-ng" \
    "${package_root}/bin/bcftools" \
    "${package_root}/libexec/bcftools.bin" \
    "${package_root}/libexec/vcftools-ng.bin"

export LD_LIBRARY_PATH="${dependency_prefix}/lib:${dependency_prefix}/lib64"

declare -a dependency_queue=(
    "${package_root}/libexec/vcftools-ng.bin"
    "${package_root}/libexec/bcftools.bin"
)
declare -a dependency_sources=()
queue_index=0

while ((queue_index < ${#dependency_queue[@]})); do
    dependency_file=${dependency_queue[queue_index]}
    queue_index=$((queue_index + 1))

    while IFS='|' read -r soname source_path; do
        case "${soname}" in
            libc.so.*|libm.so.*|libpthread.so.*|libdl.so.*|librt.so.*|\
            libresolv.so.*|libutil.so.*|ld-linux-x86-64.so.*)
                continue
                ;;
        esac

        target_path="${package_root}/lib/${soname}"
        if [[ ! -e "${target_path}" ]]; then
            cp -L "${source_path}" "${target_path}"
            chmod 0755 "${target_path}"
            dependency_queue+=("${target_path}")
            dependency_sources+=("${source_path}")
        fi
    done < <(
        ldd "${dependency_file}" \
            | awk '$2 == "=>" && $3 ~ /^\// { print $1 "|" $3 }'
    )
done

unset LD_LIBRARY_PATH
export LD_LIBRARY_PATH="${package_root}/lib"

for executable in \
    "${package_root}/libexec/vcftools-ng.bin" \
    "${package_root}/libexec/bcftools.bin"; do
    if ldd "${executable}" | grep -q "not found"; then
        ldd "${executable}" >&2
        echo "error: unresolved runtime dependency in ${executable}" >&2
        exit 1
    fi
done

for source_path in "${dependency_sources[@]}"; do
    if ! rpm -qf "${source_path}" >/dev/null 2>&1; then
        continue
    fi
    rpm_package=$(
        rpm -qf --queryformat '%{NAME}' "${source_path}" 2>/dev/null
    )
    license_target="${package_root}/LICENSES/system-${rpm_package}"
    mkdir -p "${license_target}"
    while IFS= read -r license_path; do
        [[ -f "${license_path}" ]] || continue
        cp -L "${license_path}" \
            "${license_target}/$(basename "${license_path}")"
    done < <(rpm -qL "${rpm_package}" | grep '^/usr/share/licenses/' || true)
    rmdir "${license_target}" 2>/dev/null || true
done

strip --strip-unneeded \
    "${package_root}/libexec/vcftools-ng.bin" \
    "${package_root}/libexec/bcftools.bin"
find "${package_root}/lib" -type f -exec strip --strip-unneeded {} + \
    2>/dev/null || true

"${package_root}/bin/vcftools-ng" --version
"${package_root}/bin/bcftools" --version | head -n 2

rm -f "/output/${archive_name}.tar.gz" \
    "/output/${archive_name}.tar.gz.sha256"
pushd "${work_directory}" >/dev/null
find "${archive_name}" -print0 \
    | LC_ALL=C sort -z \
    | tar \
        --null \
        --no-recursion \
        --files-from=- \
        --mtime='UTC 2026-07-25' \
        --owner=0 \
        --group=0 \
        --numeric-owner \
        -czf "/output/${archive_name}.tar.gz"
popd >/dev/null
pushd /output >/dev/null
sha256sum "${archive_name}.tar.gz" \
    >"${archive_name}.tar.gz.sha256"
popd >/dev/null

chown "${HOST_UID}:${HOST_GID}" \
    "/output/${archive_name}.tar.gz" \
    "/output/${archive_name}.tar.gz.sha256"
