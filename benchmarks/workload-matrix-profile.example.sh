# Copy this file outside the repository, replace every placeholder, and pass
# the copy to run-workload-matrix.sh. Paths may contain spaces.

INPUT_OPTION=--gzvcf
INPUT_VCF=/absolute/path/to/locked-development-subset.vcf.gz
ORACLE_ROOT=/absolute/path/to/hash-locked-original-oracles
RESULT_ROOT=/absolute/path/to/workload-matrix-results

# Deterministic fixtures produced once and retained with SHA256SUMS.
POSITION_SPARSE_SORTED=/absolute/path/to/positions-1pct.sorted.txt
POSITION_SPARSE_UNSORTED=/absolute/path/to/positions-1pct.unsorted.txt
POSITION_DENSE=/absolute/path/to/positions-50pct.txt
POSITION_EXCLUDE_DENSE=/absolute/path/to/exclude-positions-50pct.txt
KEEP_25=/absolute/path/to/keep-25pct.txt
KEEP_50=/absolute/path/to/keep-50pct.txt

# Real population pairs. Use the smallest and largest selected-sample pairs.
POP_SMALL_1=/absolute/path/to/small-population-1.txt
POP_SMALL_2=/absolute/path/to/small-population-2.txt
POP_LARGE_1=/absolute/path/to/large-population-1.txt
POP_LARGE_2=/absolute/path/to/large-population-2.txt

# Current Osmanthus production profile.
PI_WINDOW_SIZE=100000
PI_WINDOW_STEP=10000
TAJIMA_WINDOW_SIZE=100000
FST_WINDOW_SIZE=100000
FST_WINDOW_STEP=10000

# Development defaults. The runner accepts a shell array here.
THREADS=(1 4 8 16 32)
REPEATS=3

# A space-separated subset may be supplied for targeted profiling.
# RUN_CASES="W01_counts W02_seven_filter_counts"

# Keep repeat-one outputs and remove later candidate artifacts after their
# timing/hash records have been written. Set to 1 to retain every repeat.
KEEP_REPEAT_OUTPUTS=0
