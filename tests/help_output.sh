#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 /path/to/vcftools-ng" >&2
    exit 2
fi

ng=$1
work=$(mktemp -d)
trap 'rm -rf -- "$work"' EXIT

"$ng" -h >"$work/short.txt" 2>"$work/short.err"
"$ng" --help >"$work/long.txt" 2>"$work/long.err"
env -u NO_COLOR CLICOLOR_FORCE=1 \
    "$ng" -h >"$work/color.txt" 2>"$work/color.err"
NO_COLOR=1 CLICOLOR_FORCE=1 \
    "$ng" -h >"$work/no-color.txt" 2>"$work/no-color.err"

test ! -s "$work/short.err"
test ! -s "$work/long.err"
test ! -s "$work/color.err"
test ! -s "$work/no-color.err"
cmp "$work/short.txt" "$work/long.txt"
cmp "$work/short.txt" "$work/no-color.txt"

grep -Fq $'\033[1;36m' "$work/color.txt"
grep -Fq $'\033[1;34m' "$work/color.txt"
grep -Fq $'\033[32m' "$work/color.txt"
grep -Fq $'\033[1;33m' "$work/color.txt"
sed $'s/\033\\[[0-9;]*m//g' \
    "$work/color.txt" >"$work/color-stripped.txt"
cmp "$work/short.txt" "$work/color-stripped.txt"

grep -Eq '^vcftools-ng [0-9]+\.[0-9]+\.[0-9]+$' "$work/short.txt"
grep -Fq 'Usage: vcftools-ng INPUT [FILTERS] OUTPUT [OPTIONS]' \
    "$work/short.txt"

for heading in \
    'QUICK EXAMPLES:' \
    'GENERAL OPTIONS:' \
    'RUN LOGGING:' \
    'INPUT OPTIONS (choose one):' \
    'SITE STATISTICS OUTPUT:' \
    'INDIVIDUAL STATISTICS OUTPUT:' \
    'DIVERSITY AND POPULATION STATISTICS:' \
    'LD AND PCA OUTPUT:' \
    'RECODE AND FORMAT OUTPUT:' \
    'TWO-FILE COMPARISON:' \
    'CHROMOSOME, POSITION, AND BED FILTERS:' \
    'SAMPLE FILTERS:' \
    'ALLELE, QUALITY, DEPTH, AND MISSINGNESS FILTERS:' \
    'FREQUENCY, COUNT, AND HWE FILTERS:' \
    'VCF FILTER, INFO, AND GENOTYPE FILTERS:' \
    'ADAPTIVE INPUT AND INDEX POLICY (default --input-backend auto):' \
    'IMPORTANT COMBINATION RULES:' \
    'COMPATIBILITY:' \
    'TERMINAL COLORS:'
do
    grep -Fqx "$heading" "$work/short.txt"
done

options=(
    --help --version --vcf --gzvcf --bcf --input --out
    --log-file --no-log-file --threads
    --batch-size --compat --input-backend --bcftools
    --freq --freq2 --counts --missing-site --site-depth --site-mean-depth
    --corrected-depth-arithmetic
    --depth --missing-indv --het --hardy --site-quality
    --site-pi --window-pi --window-pi-step --TajimaD
    --weir-fst-pop --fst-window-size --fst-window-step
    --geno-r2 --ld-window --ld-window-min --ld-window-bp
    --ld-window-bp-min --min-r2 --pca --pca-no-norm
    --recode --recode-vcf --recode-bcf --recode-vcf-gz
    --recode-INFO-all --stdout
    --diff --gzdiff --diff-bcf --diff-site --diff-indv
    --diff-site-discordance --diff-indv-discordance
    --chr --not-chr --from-bp --to-bp --positions --exclude-positions
    --bed --exclude-bed --thin
    --indv --remove-indv --keep --remove
    --min-alleles --max-alleles --remove-indels --keep-only-indels --minQ
    --minGQ --minDP --maxDP --min-meanDP --max-meanDP
    --max-missing --max-missing-count
    --maf --max-maf --mac --max-mac --hwe
    --non-ref-af --max-non-ref-af --non-ref-af-any --max-non-ref-af-any
    --non-ref-ac --max-non-ref-ac --non-ref-ac-any --max-non-ref-ac-any
    --keep-filtered --remove-filtered --remove-filtered-all
    --keep-INFO --remove-INFO
    --remove-filtered-geno --remove-filtered-geno-all
)

for option in "${options[@]}"; do
    grep -Fq -- "$option" "$work/short.txt" || {
        echo "Help is missing supported option: $option" >&2
        exit 1
    }
done

awk '
    length($0) > 100 {
        print "Help line exceeds 100 columns: " $0 > "/dev/stderr"
        failed = 1
    }
    END {
        if (NR < 150) {
            print "Help output is unexpectedly short" > "/dev/stderr"
            failed = 1
        }
        exit failed
    }
' "$work/short.txt"

echo "HELP_OUTPUT_PASS"
