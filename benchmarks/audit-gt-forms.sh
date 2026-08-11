#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 INPUT.vcf[.gz]" >&2
    exit 2
fi

input=$1
[[ -f $input ]] || { echo "GT audit input not found: $input" >&2; exit 2; }

scan_gt() {
    awk '
        BEGIN { FS = "\t"; OFS = "\t" }
        /^#/ { next }
        {
            ++records
            gt_index = 0
            n_format = split($9, format, ":")
            for (field = 1; field <= n_format; ++field) {
                if (format[field] == "GT") {
                    gt_index = field
                    break
                }
            }
            if (gt_index == 0) {
                ++records_without_gt
                next
            }
            for (sample = 10; sample <= NF; ++sample) {
                n_sample = split($sample, values, ":")
                gt = gt_index <= n_sample ? values[gt_index] : "<ABSENT>"
                ++count[gt]
                if (!(gt in example))
                    example[gt] = $1 ":" $2
            }
        }
        END {
            print "#records", records
            print "#records_without_GT", records_without_gt + 0
            print "GT", "cells", "first_site"
            for (gt in count)
                print gt, count[gt], example[gt]
        }
    '
}

sort_audit() {
    IFS= read -r records
    IFS= read -r without_gt
    IFS= read -r header
    printf '%s\n%s\n%s\n' "$records" "$without_gt" "$header"
    LC_ALL=C sort -t $'\t' -k2,2nr -k1,1
}

if [[ $input == *.gz ]]; then
    gzip -dc -- "$input" | scan_gt | sort_audit
else
    scan_gt < "$input" | sort_audit
fi
