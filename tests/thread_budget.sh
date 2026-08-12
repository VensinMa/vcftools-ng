#!/usr/bin/env bash
set -euo pipefail

ng=${1:?vcftools-ng executable required}
bcftools=${2:?bcftools executable required}
work=$(mktemp -d "${TMPDIR:-/tmp}/vcftools-ng-thread-budget.XXXXXX")
cleanup() {
    rm -rf -- "$work"
}
trap cleanup EXIT INT TERM

# Keep the process alive long enough for /proc sampling while remaining a
# small, generated CI fixture. Distinct positions provide many indexed shards.
plain="$work/input.vcf"
awk 'BEGIN {
    print "##fileformat=VCFv4.2"
    print "##contig=<ID=chr1,length=300000>"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"
    for (i = 1; i <= 300000; ++i) {
        print "chr1\t" i "\t.\tA\tG\t30\tPASS\t.\tGT\t0/1"
    }
}' >"$plain"

"$ng" --vcf "$plain" --threads 1 --recode --out "$work/source" \
    >"$work/source.stdout" 2>"$work/source.stderr"
input="$work/source.recode.vcf.gz"
"$bcftools" index --tbi "$input"
gzip -c "$plain" >"$work/input.vcf.gz"
"$bcftools" view --no-version -Ob -o "$work/input.bcf" "$plain"
"$bcftools" index --csi "$work/input.bcf"

gzip -dc "$work/source.recode.vcf.gz" >"$work/source.vcf"

# The strict budget is the number of CPUs on which the complete process tree
# may run. I/O and ordered-queue workers may overlap while blocked, but every
# thread and child inherits the same N-CPU affinity. Cover the three text-input
# shapes because their caller/worker sharing differs.
shopt -s nullglob
affinity_cpu_count() {
    awk '/^Cpus_allowed_list:/{print $2}' "$1" | awk -F, '{
        total=0
        for (i=1; i<=NF; ++i) {
            split($i, part, "-")
            total += part[2] == "" ? 1 : part[2]-part[1]+1
        }
        print total
    }'
}

tree_max_affinity_count() {
    local root=$1
    local maximum=0
    local -a pending=("$root")
    while ((${#pending[@]} > 0)); do
        local index=$((${#pending[@]} - 1))
        local pid=${pending[$index]}
        unset 'pending[index]'
        [[ -r "/proc/$pid/status" ]] || continue
        local count
        count=$(affinity_cpu_count "/proc/$pid/status")
        if ((count > maximum)); then
            maximum=$count
        fi
        local children_file="/proc/$pid/task/$pid/children"
        if [[ -r "$children_file" ]]; then
            local children
            children=$(cat "$children_file" 2>/dev/null || true)
            for child in $children; do
                pending+=("$child")
            done
        fi
    done
    printf '%s\n' "$maximum"
}

run_and_sample() {
    local budget=$1
    local scenario=$2
    shift 2
    local prefix="$work/${scenario}-${budget}"
    "$ng" "$@" --threads "$budget" --recode --out "$prefix" \
        >"$prefix.stdout" 2>"$prefix.stderr" &
    local pid=$!
    local peak=0
    : >"$prefix.affinity"
    while [[ -d "/proc/$pid/task" ]]; do
        local tasks
        tasks=$(tree_thread_count "$pid")
        if ((tasks > peak)); then
            peak=$tasks
        fi
        tree_max_affinity_count "$pid" \
            >"$prefix.affinity" 2>/dev/null || true
    done
    if ! wait "$pid"; then
        sed -n '1,240p' "$prefix.stderr" >&2
        return 1
    fi
    # Worker pools may overlap while blocked on I/O or ordered queues. The
    # process-wide affinity is the strict CPU budget: every thread and child
    # inherits at most N runnable CPUs even when more waitable pthreads exist.
    local affinity_count
    affinity_count=$(<"$prefix.affinity")
    if ((affinity_count > budget)); then
        echo "CPU budget exceeded: scenario=$scenario requested=$budget affinity=$affinity_count" >&2
        return 1
    fi
    grep -Fq "Threads: $budget (user specified)" "$prefix.stderr"
    if ! grep -Fq 'Output compression threads:' "$prefix.stderr"; then
        grep -Fq 'VCF compression ' "$prefix.stderr"
    fi
    if ! grep -Fq 'Coordinator threads:' "$prefix.stderr"; then
        grep -Fq 'BCF serialization ' "$prefix.stderr"
    fi
    gzip -dc "$prefix.recode.vcf.gz" >"$prefix.vcf"
    if [[ "$scenario" == general-bgzf ]]; then
        if ((budget == 1)); then
            cp "$prefix.vcf" "$work/general-reference.vcf"
        else
            cmp "$work/general-reference.vcf" "$prefix.vcf"
        fi
    else
        cmp "$work/source.vcf" "$prefix.vcf"
    fi
}

# Count the root process and every live descendant.  Automatic index builds
# execute bcftools as a child, so sampling only /proc/$pid/task would miss a
# duplicated subprocess thread pool.
tree_thread_count() {
    local root=$1
    local total=0
    local -a pending=("$root")
    while ((${#pending[@]} > 0)); do
        local index=$((${#pending[@]} - 1))
        local pid=${pending[$index]}
        unset 'pending[index]'
        local task_paths=("/proc/$pid/task"/*)
        if [[ ! -d "/proc/$pid/task" ]]; then
            continue
        fi
        total=$((total + ${#task_paths[@]}))
        local children_file="/proc/$pid/task/$pid/children"
        if [[ -r "$children_file" ]]; then
            local children
            children=$(cat "$children_file" 2>/dev/null || true)
            for child in $children; do
                pending+=("$child")
            done
        fi
    done
    printf '%s\n' "$total"
}

run_and_sample_output() {
    local budget=$1
    local scenario=$2
    local suffix=$3
    shift 3
    local prefix="$work/${scenario}-${budget}"
    "$ng" "$@" --threads "$budget" --out "$prefix" \
        >"$prefix.stdout" 2>"$prefix.stderr" &
    local pid=$!
    local peak=0
    : >"$prefix.affinity"
    while [[ -d "/proc/$pid/task" ]]; do
        local tasks
        tasks=$(tree_thread_count "$pid")
        if ((tasks > peak)); then
            peak=$tasks
        fi
        tree_max_affinity_count "$pid" \
            >"$prefix.affinity" 2>/dev/null || true
    done
    if ! wait "$pid"; then
        sed -n '1,240p' "$prefix.stderr" >&2
        return 1
    fi
    local affinity_count
    affinity_count=$(<"$prefix.affinity")
    if ((affinity_count > budget)); then
        echo "CPU budget exceeded: scenario=$scenario requested=$budget affinity=$affinity_count" >&2
        return 1
    fi
    if ((budget == 1)); then
        cp "$prefix$suffix" "$work/$scenario.reference"
    else
        cmp "$work/$scenario.reference" "$prefix$suffix"
    fi
}

for budget in 1 2 4 8; do
    # Automatic BGZF is a stream at one thread and indexed at two or more,
    # covering the real adaptive transitions without forcing the general
    # compatibility pipeline.
    run_and_sample "$budget" adaptive-bgzf --gzvcf "$input"
    # Ordinary gzip can only stream and therefore exercises the fused
    # HTSlib-reader/compute/output allocation at every budget.
    run_and_sample "$budget" gzip-stream --gzvcf "$work/input.vcf.gz"
    run_and_sample "$budget" plain --vcf "$plain"
done

# An unindexed BGZF recode invokes bcftools to construct CSI at two or more
# threads. Include the child process in the same strict whole-workflow budget.
for budget in 1 2 4 8; do
    cp "$input" "$work/auto-index-$budget.vcf.gz"
    run_and_sample_output "$budget" auto-index \
        .recode.vcf.gz \
        --gzvcf "$work/auto-index-$budget.vcf.gz" --recode
    if ((budget == 1)); then
        test ! -e "$work/auto-index-$budget.vcf.gz.csi"
    else
        test -s "$work/auto-index-$budget.vcf.gz.csi"
    fi
done

# Force the general HTSlib pipeline and its output allocator. Sample
# selection deliberately makes compressed text recode ineligible for the
# fused parser, exercising reader/input/compute/compression coordination.
printf 'S1\n' >"$work/keep.txt"
for budget in 1 2 4 8; do
    run_and_sample "$budget" general-bgzf \
        --gzvcf "$input" --keep "$work/keep.txt"
done


# Exercise output serialization/compression, two-input comparison, and the
# post-scan LD/PCA pools. Their one-thread outputs are independent goldens;
# every higher budget must preserve exact bytes while never exceeding N live
# threads across the process tree.
for budget in 1 2 4 8; do
    run_and_sample_output "$budget" bcf-output .recode.bcf \
        --bcf "$work/input.bcf" --recode-bcf --recode-INFO-all
    run_and_sample_output "$budget" bcf-diff .diff.sites \
        --bcf "$work/input.bcf" --diff-bcf "$work/input.bcf" \
        --diff-site-discordance
    run_and_sample_output "$budget" genotype-ld .geno.ld \
        --bcf "$work/input.bcf" --geno-r2 --ld-window 2
    run_and_sample_output "$budget" pca .pca \
        --bcf "$work/input.bcf" --pca
done

echo THREAD_BUDGET_PASS
