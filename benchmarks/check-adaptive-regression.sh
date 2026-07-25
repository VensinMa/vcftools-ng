#!/usr/bin/env bash
set -euo pipefail

if (( $# != 1 )); then
    printf 'Usage: %s SUMMARY.tsv\n' "$0" >&2
    exit 2
fi

summary=$1
if [[ ! -s "$summary" ]]; then
    printf 'Missing benchmark summary: %s\n' "$summary" >&2
    exit 2
fi

awk -F '\t' '
    NR == 1 {
        for (column = 1; column <= NF; ++column) {
            field[$column] = column
        }
        detailed = ("engine" in field)
        required["bgzf_auto_csi"] = 1
        required["bgzf_no_auto_index"] = 1
        required["bgzf_tbi"] = 1
        required["bcf_auto_csi"] = 1
        required["bcf_csi"] = 1
        required["bcf_no_auto_index"] = 1
        required["plain_vcf"] = 1
        next
    }
    detailed && $field["engine"] != "vcftools-ng" { next }
    $field["threads"] != 1 && $field["threads"] != 2 &&
    $field["threads"] != 4 && $field["threads"] != 8 &&
    $field["threads"] != 16 && $field["threads"] != 32 { next }
    !($field["case"] in required) { next }
    {
        key = $field["case"] SUBSEP $field["threads"]
        seen[key] = 1
        speedup = $field["speedup_vs_original"]
        if (!detailed) {
            speedup = $field["speedup"]
        }
        exact = $field["exact"]
        if (exact != "PASS" || speedup + 0.0 < 1.0) {
            printf "FAIL\t%s\tthreads=%s\tspeedup=%s\texact=%s\n", \
                $field["case"], $field["threads"], speedup, exact
            failed = 1
        } else {
            printf "PASS\t%s\tthreads=%s\tspeedup=%s\texact=%s\n", \
                $field["case"], $field["threads"], speedup, exact
        }
    }
    END {
        for (case_name in required) {
            split("1 2 4 8 16 32", thread_values, " ")
            for (thread_index = 1;
                 thread_index <= 6; ++thread_index) {
                threads = thread_values[thread_index]
                key = case_name SUBSEP threads
                if (!(key in seen)) {
                    printf "MISSING\t%s\tthreads=%d\n", \
                        case_name, threads
                    failed = 1
                }
            }
        }
        exit failed
    }
' "$summary"
