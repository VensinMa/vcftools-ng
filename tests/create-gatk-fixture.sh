#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source_vcf=${1:?source BGZF VCF required}
output_prefix=${2:-"$project_root/tests/fixtures/osmanthus205.gatk.23chr_1k"}
threads=${3:-32}
bcftools_bin=${BCFTOOLS:-bcftools}
bgzip_bin=${BGZIP:-bgzip}
tabix_bin=${TABIX:-tabix}

if [[ ! -f $source_vcf ]]; then
    printf 'Source VCF does not exist: %s\n' "$source_vcf" >&2
    exit 1
fi
for path in \
    "$output_prefix.vcf.gz" "$output_prefix.vcf.gz.tbi" \
    "$output_prefix.bcf" "$output_prefix.bcf.csi" \
    "$output_prefix.provenance.txt"; do
    if [[ -e $path ]]; then
        printf 'Refusing to replace existing fixture artifact: %s\n' \
            "$path" >&2
        exit 1
    fi
done

mkdir -p "$project_root/tests/output" "$(dirname "$output_prefix")"
work=$(mktemp -d "$project_root/tests/output/gatk-fixture.XXXXXX")
completed=false
cleanup() {
    if $completed; then
        rm -rf -- "$work"
    else
        printf 'Preserved failed fixture workspace for recovery: %s\n' \
            "$work" >&2
    fi
}
trap cleanup EXIT INT TERM

plain_vcf="$work/subset.vcf"
source_sha="$work/source.sha256"

# Read the 121-GiB source only once. tee keeps hashing the compressed source
# after awk has obtained 1,000 records from Chr23 and closes the extraction
# branch. A SIGPIPE from the now-unneeded decompressor is therefore expected.
set +e
set +o pipefail
tee -p >(sha256sum >"$source_sha") <"$source_vcf" |
    "$bgzip_bin" -dc |
    awk -F '\t' '
        /^#/ {
            print
            next
        }
        $1 ~ /^Chr(0[1-9]|1[0-9]|2[0-3])$/ && count[$1] < 1000 {
            print
            count[$1]++
            if (count[$1] == 1000) {
                print "Extracted 1,000 records from " $1 > "/dev/stderr"
                if ($1 == "Chr23")
                    exit
            }
        }
    ' >"$plain_vcf"
pipeline_status=("${PIPESTATUS[@]}")
set -o pipefail
set -e
if (( pipeline_status[0] != 0 || pipeline_status[2] != 0 )); then
    printf 'Source read or record extraction failed (tee=%s awk=%s)\n' \
        "${pipeline_status[0]}" "${pipeline_status[2]}" >&2
    exit 1
fi
if [[ ! -s $source_sha ]]; then
    printf 'Source SHA-256 was not completed\n' >&2
    exit 1
fi

counts="$work/counts.tsv"
awk -F '\t' '!/^#/ {count[$1]++} END {
    for (chromosome = 1; chromosome <= 23; ++chromosome) {
        name = sprintf("Chr%02d", chromosome)
        print name "\t" count[name]
        if (count[name] != 1000)
            failed = 1
    }
    exit failed
}' "$plain_vcf" >"$counts"
record_count=$(awk '!/^#/ {count++} END {print count+0}' "$plain_vcf")
sample_count=$(
    awk -F '\t' '/^#CHROM/ {print (NF > 9 ? NF - 9 : 0)}' "$plain_vcf"
)
if [[ $record_count != 23000 || $sample_count != 205 ]]; then
    printf 'Unexpected fixture shape: records=%s samples=%s\n' \
        "$record_count" "$sample_count" >&2
    exit 1
fi

vcf_gz="$work/subset.vcf.gz"
"$bgzip_bin" --threads "$threads" --output "$vcf_gz" "$plain_vcf"
"$tabix_bin" -p vcf "$vcf_gz"

bcf="$work/subset.bcf"
"$bcftools_bin" view --no-version --threads "$threads" \
    -Ob -o "$bcf" "$vcf_gz"
"$bcftools_bin" index --csi --threads "$threads" "$bcf"

vcf_records=$("$bcftools_bin" index -n "$vcf_gz")
bcf_records=$("$bcftools_bin" index -n "$bcf" 2>/dev/null)
if [[ $vcf_records != 23000 || $bcf_records != 23000 ]]; then
    printf 'Index record-count mismatch: VCF=%s BCF=%s\n' \
        "$vcf_records" "$bcf_records" >&2
    exit 1
fi

provenance="$work/provenance.txt"
fixture_name=$(basename "$output_prefix")
{
    printf 'Fixture: %s\n' "$fixture_name"
    printf 'Source: %s\n' "$source_vcf"
    printf 'Source bytes: %s\n' "$(stat -c %s "$source_vcf")"
    printf 'Source mtime: %s\n' "$(stat -c %y "$source_vcf")"
    printf 'Source SHA-256: %s\n' "$(awk '{print $1}' "$source_sha")"
    printf 'Selection: first 1,000 records from Chr01 through Chr23\n'
    printf 'Samples: %s\nRecords: %s\n' "$sample_count" "$record_count"
    printf 'Extraction preserves the source VCF header and selected text records.\n'
    printf '\nPer-contig records:\n'
    cat "$counts"
    printf '\nArtifact SHA-256:\n'
    printf '%s  %s.vcf.gz\n' \
        "$(sha256sum "$vcf_gz" | awk '{print $1}')" "$fixture_name"
    printf '%s  %s.vcf.gz.tbi\n' \
        "$(sha256sum "$vcf_gz.tbi" | awk '{print $1}')" "$fixture_name"
    printf '%s  %s.bcf\n' \
        "$(sha256sum "$bcf" | awk '{print $1}')" "$fixture_name"
    printf '%s  %s.bcf.csi\n' \
        "$(sha256sum "$bcf.csi" | awk '{print $1}')" "$fixture_name"
    printf '\nTool versions:\n'
    "$bcftools_bin" --version | sed -n '1,2p'
    "$bgzip_bin" --help 2>&1 | sed -n '2p'
} >"$provenance"

mv "$vcf_gz" "$output_prefix.vcf.gz"
mv "$vcf_gz.tbi" "$output_prefix.vcf.gz.tbi"
mv "$bcf" "$output_prefix.bcf"
mv "$bcf.csi" "$output_prefix.bcf.csi"
mv "$provenance" "$output_prefix.provenance.txt"

completed=true
printf 'Created GATK compatibility fixture: %s\n' "$output_prefix"
