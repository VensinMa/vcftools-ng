#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
htslib_root=${HTSLIB_ROOT:-/home/vensin/anaconda3/pkgs/htslib-1.23.1-h633afcb_0}
build_dir=${BUILD_DIR:-$project_root/build-pgo-v0141}
profile_dir=${PROFILE_DIR:-$build_dir/profile-data}
train_root=${TRAIN_ROOT:-$project_root/benchmarks/results/workload-matrix-23k-v0130}
input=$train_root/osmanthus412.flags.23chr_1k.vcf
fixtures=$train_root/fixtures

if [[ -e $build_dir ]]; then
    echo "Refusing to overwrite PGO build directory: $build_dir" >&2
    echo "Choose a fresh BUILD_DIR or remove the previous PGO build explicitly." >&2
    exit 2
fi
for path in "$input" "$fixtures/keep-50pct.samples.txt" \
            "$fixtures/pop-large-wild-166.txt" \
            "$fixtures/pop-large-cultivated-243.txt"; do
    [[ -e $path ]] || {
        echo "v0.14.1 PGO training asset missing: $path" >&2
        exit 2
    }
done

cmake -S "$project_root" -B "$build_dir" \
    -DCMAKE_BUILD_TYPE=Release \
    -DHTSLIB_ROOT="$htslib_root" \
    -DVCFTOOLS_NG_PGO=GENERATE \
    -DVCFTOOLS_NG_PGO_DIR="$profile_dir"
cmake --build "$build_dir" -j "${BUILD_JOBS:-$(nproc)}"

training_output=$(mktemp -d /tmp/vcftools-ng-pgo-training.XXXXXX)
trap 'rm -rf -- "$training_output"' EXIT
ng=$build_dir/vcftools-ng

run_training() {
    local name=$1
    shift
    "$ng" --vcf "$input" --threads "${TRAIN_THREADS:-8}" \
        "$@" --no-log-file --out "$training_output/$name" \
        >/dev/null 2>"$training_output/$name.stderr"
}

run_training counts --counts
run_training ten-filter \
    --min-alleles 2 --max-alleles 2 --minQ 30 --minGQ 10 \
    --minDP 5 --maxDP 30 --min-meanDP 7 --max-missing 0.9 \
    --maf 0.1 --mac 2 --counts
run_training sample-recode \
    --keep "$fixtures/keep-50pct.samples.txt" \
    --min-alleles 2 --max-alleles 2 --minQ 30 --minGQ 10 \
    --minDP 5 --maxDP 30 --min-meanDP 7 --max-missing 0.9 \
    --maf 0.1 --mac 2 --recode-vcf-gz --recode-INFO-all
run_training window-pi --window-pi 100000 --window-pi-step 10000
run_training window-fst \
    --min-alleles 2 --max-alleles 2 \
    --weir-fst-pop "$fixtures/pop-large-wild-166.txt" \
    --weir-fst-pop "$fixtures/pop-large-cultivated-243.txt" \
    --fst-window-size 100000 --fst-window-step 10000
run_training ld \
    --min-alleles 2 --max-alleles 2 --geno-r2 \
    --ld-window 200 --ld-window-bp 1000000000 --min-r2 0.99

find "$profile_dir" -type f -print -quit | grep -q . || {
    echo "Compiler produced no PGO profile data in $profile_dir" >&2
    exit 1
}

compiler_description=$(find "$build_dir/CMakeFiles" \
    -path '*/CMakeCXXCompiler.cmake' -print -quit)
[[ -n $compiler_description ]] || {
    echo "Could not locate CMake's C++ compiler description" >&2
    exit 1
}
compiler_id=$(sed -n \
    's/^set(CMAKE_CXX_COMPILER_ID "\([^"]*\)")/\1/p' \
    "$compiler_description" | head -1)
[[ -n $compiler_id ]] || {
    echo "Could not determine the C++ compiler family" >&2
    exit 1
}
if [[ $compiler_id == Clang || $compiler_id == AppleClang ]]; then
    llvm_profdata=${LLVM_PROFDATA:-$(command -v llvm-profdata || true)}
    [[ -n $llvm_profdata ]] || {
        echo "Clang PGO requires llvm-profdata" >&2
        exit 1
    }
    "$llvm_profdata" merge -output="$profile_dir/merged.profdata" \
        "$profile_dir"/*.profraw
fi

cmake -S "$project_root" -B "$build_dir" \
    -DCMAKE_BUILD_TYPE=Release \
    -DHTSLIB_ROOT="$htslib_root" \
    -DVCFTOOLS_NG_PGO=USE \
    -DVCFTOOLS_NG_PGO_DIR="$profile_dir"
cmake --build "$build_dir" -j "${BUILD_JOBS:-$(nproc)}"

echo "PGO build complete: $build_dir/vcftools-ng"
