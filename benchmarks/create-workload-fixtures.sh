#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 INPUT.vcf RESULT_DIR REFERENCE_FIXTURE_DIR" >&2
    exit 2
fi

input=$1
result=$2
reference=$3
fixtures=$result/fixtures

[[ -f $input ]] || { echo "Input VCF missing: $input" >&2; exit 2; }
[[ -d $reference ]] || {
    echo "Reference fixture directory missing: $reference" >&2
    exit 2
}
[[ ! -e $fixtures ]] || {
    echo "Refusing to overwrite fixture directory: $fixtures" >&2
    exit 2
}

mkdir -p "$fixtures"
all=$fixtures/all.positions.txt
awk -F '\t' '!/^#/ {print $1 "\t" $2}' "$input" >"$all"

awk 'NR % 100 == 0' "$all" >"$fixtures/positions-1pct.sorted.txt"
tac "$fixtures/positions-1pct.sorted.txt" \
    >"$fixtures/positions-1pct.shuffled.txt"
awk 'NR % 2 == 0' "$all" >"$fixtures/positions-50pct.sorted.txt"
tac "$fixtures/positions-50pct.sorted.txt" \
    >"$fixtures/positions-50pct.shuffled.txt"

awk '
    { print }
    NR % 101 == 0 { print }
    END {
        for (item = 1; item <= 64; ++item)
            print "__absent_contig__\t" item
    }
' "$fixtures/positions-1pct.sorted.txt" \
    >"$fixtures/exclude-1pct.duplicates-and-absent.txt"
awk '
    { print }
    NR % 1009 == 0 { print }
    END {
        for (item = 1; item <= 64; ++item)
            print "__absent_contig__\t" item
    }
' "$fixtures/positions-50pct.sorted.txt" \
    >"$fixtures/exclude-50pct.duplicates-and-absent.txt"

for name in \
    keep-25pct.samples.txt keep-50pct.samples.txt \
    keep-100pct.samples.txt pop-small-ancient-12.txt \
    pop-small-asiaticus-13.txt pop-large-wild-166.txt \
    pop-large-cultivated-243.txt; do
    cp -- "$reference/$name" "$fixtures/$name"
done

awk -F '\t' '/^#CHROM/ {
    for (column = 10; column <= NF; ++column)
        print $column
}' "$input" >"$fixtures/all.samples.txt"

sha256sum "$input" >"$result/input.sha256"
(
    cd "$fixtures"
    sha256sum ./* >SHA256SUMS
)

records=$(awk '!/^#/ {count++} END {print count+0}' "$input")
samples=$(wc -l <"$fixtures/all.samples.txt")
printf 'Input: %s\nRecords: %s\nSamples: %s\n' \
    "$input" "$records" "$samples" >"$result/FIXTURE_INFO.txt"

echo "Created deterministic workload fixtures: $fixtures"
