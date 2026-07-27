#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/vcf.h>
#include <zlib.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

#include <unistd.h>

#include "fast_counts.h"
#include "input_source.h"

namespace {

constexpr const char* kVersion = "0.12.3";

struct Options {
    std::string input;
    bool input_bcf = false;
    std::string output_prefix = "out";
    unsigned threads =
        vcftools_ng::input::detect_available_threads().count;
    unsigned requested_threads = threads;
    bool threads_explicit = false;
    std::size_t batch_size = 2048;
    vcftools_ng::input::Backend input_backend =
        vcftools_ng::input::Backend::automatic;
    std::string bcftools_path = "bcftools";

    bool output_freq = false;
    bool output_freq2 = false;
    bool output_counts = false;
    bool output_missing_site = false;
    bool output_site_depth = false;
    bool output_site_mean_depth = false;
    bool output_individual_depth = false;
    bool output_individual_missingness = false;
    bool output_heterozygosity = false;
    bool output_hardy_weinberg = false;
    bool output_site_quality = false;
    bool output_site_pi = false;
    int pi_window_size = 0;
    int pi_window_step = 0;
    int tajima_window_size = 0;
    std::vector<std::string> fst_population_files;
    int fst_window_size = 0;
    int fst_window_step = 0;
    bool output_genotype_r2 = false;
    int ld_snp_window_size = std::numeric_limits<int>::max();
    int ld_snp_window_min = -1;
    int ld_bp_window_size = std::numeric_limits<int>::max();
    int ld_bp_window_min = -1;
    double min_r2 = -1.0;
    bool output_pca = false;
    bool pca_normalise = true;
    bool output_recode = false;
    bool output_recode_bcf = false;
    bool output_recode_vcf_gz = false;
    bool recode_info_all = false;
    bool output_stdout = false;
    std::string diff_input;
    bool output_diff_sites_in_files = false;
    bool output_diff_individuals_in_files = false;
    bool output_diff_site_discordance = false;
    bool output_diff_individual_discordance = false;

    std::set<std::string> chromosomes_to_keep;
    std::set<std::string> chromosomes_to_exclude;
    int start_position = -1;
    int end_position = std::numeric_limits<int>::max();
    std::string positions_file;
    std::string exclude_positions_file;
    std::string bed_file;
    bool bed_exclude = false;
    std::set<std::string> site_filters_to_keep;
    std::set<std::string> site_filters_to_remove;
    bool remove_all_filtered_sites = false;
    std::set<std::string> info_flags_to_keep;
    std::set<std::string> info_flags_to_remove;
    std::set<std::string> genotype_filters_to_remove;
    bool remove_all_filtered_genotypes = false;
    std::set<std::string> samples_to_keep;
    std::set<std::string> samples_to_exclude;
    std::vector<std::string> sample_keep_files;
    std::vector<std::string> sample_exclude_files;

    int min_alleles = -1;
    int max_alleles = std::numeric_limits<int>::max();
    bool remove_indels = false;
    bool keep_only_indels = false;
    double min_qual = -1.0;
    double min_gq = -1.0;
    int min_dp = -1;
    int max_dp = std::numeric_limits<int>::max();
    double min_mean_dp = -1.0;
    double max_mean_dp = std::numeric_limits<double>::max();
    double min_call_rate = 0.0;
    int max_missing_count = std::numeric_limits<int>::max();
    double min_maf = -1.0;
    double max_maf = std::numeric_limits<double>::max();
    int min_mac = -1;
    int max_mac = std::numeric_limits<int>::max();
    double min_hwe = -1.0;
    int thin_distance = -1;
    double min_non_ref_af = -1.0;
    double max_non_ref_af = std::numeric_limits<double>::max();
    double min_non_ref_af_any = -1.0;
    double max_non_ref_af_any = std::numeric_limits<double>::max();
    int min_non_ref_ac = -1;
    int max_non_ref_ac = std::numeric_limits<int>::max();
    int min_non_ref_ac_any = -1;
    int max_non_ref_ac_any = std::numeric_limits<int>::max();
};

bool can_use_fused_site_stats(const Options& options) {
    const bool any_site_stat =
        options.output_freq ||
        options.output_freq2 ||
        options.output_counts ||
        options.output_missing_site ||
        options.output_site_depth ||
        options.output_site_mean_depth ||
        options.output_site_quality;
    const bool site_stats_only =
        any_site_stat &&
        !options.output_individual_depth &&
        !options.output_individual_missingness &&
        !options.output_heterozygosity &&
        !options.output_hardy_weinberg &&
        !options.output_site_pi &&
        options.pi_window_size == 0 &&
        options.pi_window_step == 0 &&
        options.tajima_window_size == 0 &&
        options.fst_population_files.empty() &&
        options.fst_window_size == 0 &&
        options.fst_window_step == 0 &&
        !options.output_genotype_r2 &&
        options.ld_snp_window_size ==
            std::numeric_limits<int>::max() &&
        options.ld_snp_window_min == -1 &&
        options.ld_bp_window_size ==
            std::numeric_limits<int>::max() &&
        options.ld_bp_window_min == -1 &&
        options.min_r2 == -1.0 &&
        !options.output_pca &&
        !options.output_recode &&
        !options.output_recode_bcf &&
        !options.output_recode_vcf_gz &&
        !options.recode_info_all &&
        options.diff_input.empty() &&
        !options.output_diff_sites_in_files &&
        !options.output_diff_individuals_in_files &&
        !options.output_diff_site_discordance &&
        !options.output_diff_individual_discordance;
    const bool no_selection =
        options.chromosomes_to_keep.empty() &&
        options.chromosomes_to_exclude.empty() &&
        options.start_position == -1 &&
        options.end_position == std::numeric_limits<int>::max() &&
        options.positions_file.empty() &&
        options.exclude_positions_file.empty() &&
        options.bed_file.empty() &&
        options.site_filters_to_keep.empty() &&
        options.site_filters_to_remove.empty() &&
        !options.remove_all_filtered_sites &&
        options.info_flags_to_keep.empty() &&
        options.info_flags_to_remove.empty() &&
        options.genotype_filters_to_remove.empty() &&
        !options.remove_all_filtered_genotypes &&
        options.samples_to_keep.empty() &&
        options.samples_to_exclude.empty() &&
        options.sample_keep_files.empty() &&
        options.sample_exclude_files.empty();
    const bool no_numeric_filters =
        options.min_alleles == -1 &&
        options.max_alleles == std::numeric_limits<int>::max() &&
        !options.remove_indels &&
        !options.keep_only_indels &&
        options.min_qual == -1.0 &&
        options.min_gq == -1.0 &&
        options.min_dp == -1 &&
        options.max_dp == std::numeric_limits<int>::max() &&
        options.min_mean_dp == -1.0 &&
        options.max_mean_dp == std::numeric_limits<double>::max() &&
        options.min_call_rate == 0.0 &&
        options.max_missing_count == std::numeric_limits<int>::max() &&
        options.min_maf == -1.0 &&
        options.max_maf == std::numeric_limits<double>::max() &&
        options.min_mac == -1 &&
        options.max_mac == std::numeric_limits<int>::max() &&
        options.min_hwe == -1.0 &&
        options.thin_distance == -1 &&
        options.min_non_ref_af == -1.0 &&
        options.max_non_ref_af == std::numeric_limits<double>::max() &&
        options.min_non_ref_af_any == -1.0 &&
        options.max_non_ref_af_any ==
            std::numeric_limits<double>::max() &&
        options.min_non_ref_ac == -1 &&
        options.max_non_ref_ac == std::numeric_limits<int>::max() &&
        options.min_non_ref_ac_any == -1 &&
        options.max_non_ref_ac_any == std::numeric_limits<int>::max();
    const bool backend_allows_fused_stream =
        options.input_backend ==
            vcftools_ng::input::Backend::automatic ||
        (options.threads <= 2 &&
         options.input_backend ==
             vcftools_ng::input::Backend::stream);
    return site_stats_only &&
           no_selection &&
           no_numeric_filters &&
           backend_allows_fused_stream;
}

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

bool input_is_bcf(const std::string& path);

std::string require_value(int argc, char** argv, int& i) {
    if (i + 1 >= argc) {
        fail("Missing value for " + std::string(argv[i]));
    }
    return argv[++i];
}

unsigned parse_unsigned(const std::string& value, const std::string& option) {
    std::size_t used = 0;
    const unsigned long parsed = std::stoul(value, &used);
    if (used != value.size() || parsed == 0 ||
        parsed > std::numeric_limits<unsigned>::max()) {
        fail("Invalid value for " + option + ": " + value);
    }
    return static_cast<unsigned>(parsed);
}

int parse_int(const std::string& value, const std::string& option) {
    std::size_t used = 0;
    const long parsed = std::stol(value, &used);
    if (used != value.size() || parsed < std::numeric_limits<int>::min() ||
        parsed > std::numeric_limits<int>::max()) {
        fail("Invalid value for " + option + ": " + value);
    }
    return static_cast<int>(parsed);
}

double parse_double(const std::string& value, const std::string& option) {
    std::size_t used = 0;
    const double parsed = std::stod(value, &used);
    if (used != value.size()) {
        fail("Invalid value for " + option + ": " + value);
    }
    return parsed;
}

bool help_color_enabled() {
    if (std::getenv("NO_COLOR") != nullptr) {
        return false;
    }
    const char* force = std::getenv("CLICOLOR_FORCE");
    if (force != nullptr && std::strcmp(force, "0") != 0) {
        return true;
    }
    const char* term = std::getenv("TERM");
    if (term != nullptr && std::strcmp(term, "dumb") == 0) {
        return false;
    }
    return isatty(STDOUT_FILENO) == 1;
}

void emit_help(const std::string& plain) {
    if (!help_color_enabled()) {
        std::cout << plain;
        return;
    }

    constexpr const char* reset = "\033[0m";
    constexpr const char* title = "\033[1;36m";
    constexpr const char* heading = "\033[1;34m";
    constexpr const char* usage = "\033[1;33m";
    constexpr const char* option = "\033[32m";
    constexpr const char* example = "\033[36m";
    constexpr const char* warning = "\033[33m";
    constexpr const char* link = "\033[4;36m";

    std::istringstream input(plain);
    std::string line;
    bool first_line = true;
    while (std::getline(input, line)) {
        if (first_line) {
            std::cout << title << line << reset;
            first_line = false;
        } else if (
            !line.empty() && line.front() != ' ' &&
            line.back() == ':') {
            std::cout << heading << line << reset;
        } else if (line.rfind("Usage:", 0) == 0) {
            std::cout << usage << "Usage:" << reset
                      << line.substr(6);
        } else if (
            line.rfind("    vcftools-ng", 0) == 0 ||
            line.rfind("      --", 0) == 0 ||
            line.rfind("    bcftools ", 0) == 0) {
            std::cout << example << line << reset;
        } else if (line.rfind("  * ", 0) == 0) {
            std::cout << warning << line << reset;
        } else if (line.rfind("  http", 0) == 0) {
            std::cout << link << line << reset;
        } else if (line.rfind("  -", 0) == 0) {
            const std::size_t description = line.find("  ", 2);
            if (description != std::string::npos &&
                description <= 32) {
                std::cout
                    << line.substr(0, 2)
                    << option
                    << line.substr(2, description - 2)
                    << reset
                    << line.substr(description);
            } else {
                std::cout << line;
            }
        } else {
            std::cout << line;
        }
        std::cout << '\n';
    }
}

void print_help() {
    std::ostringstream help;
    help << "vcftools-ng " << kVersion << R"HELP(

High-performance, byte-compatible successor to VCFtools 0.1.17.
Adaptive VCF/BGZF/BCF input, shared statistics, exact recode, and filtering.

Usage: vcftools-ng INPUT [FILTERS] OUTPUT [OPTIONS]

Choose exactly one input and at least one output. Most single-file outputs can
be combined and are produced in one scan.

QUICK EXAMPLES:

  Filter and write an Original-compatible uncompressed VCF:
    vcftools-ng --gzvcf input.vcf.gz --threads 16 \
      --min-alleles 2 --max-alleles 2 --minQ 30 --minGQ 10 \
      --min-meanDP 7 --max-missing 0.9 --maf 0.1 \
      --recode --recode-INFO-all --out filtered

  Filter and write BGZF VCF directly (vcftools-ng extension):
    vcftools-ng --gzvcf input.vcf.gz --threads 16 \
      --minQ 30 --max-missing 0.9 --maf 0.1 \
      --recode-vcf-gz --recode-INFO-all --out filtered
    Output: filtered.recode.vcf.gz

  Generate several compatible statistics in one scan:
    vcftools-ng --bcf input.bcf --threads 16 \
      --freq --counts --missing-site --site-mean-depth --out stats

  Query one region and count alleles:
    vcftools-ng --gzvcf input.vcf.gz --chr chr1 \
      --from-bp 1 --to-bp 1000000 --counts --out chr1

GENERAL OPTIONS:
  -h, --help                   Print this help and exit
  --version                    Print the vcftools-ng version and exit
  --out PREFIX                 Prefix for output files (default: out)
  -t, --threads N              Total CPU-thread budget (default: auto-detect)
  --batch-size N               Records per pipeline batch (default: 2048)
  --compat exact               Exact VCFtools 0.1.17 compatibility mode
                               (default and currently the only mode)
  --input-backend MODE         auto|stream|plain|indexed (default: auto)
                               Advanced diagnostic/performance override
  --bcftools FILE              bcftools executable for profitable CSI builds
                               (default: bcftools)

INPUT OPTIONS (choose one):
  --vcf FILE                   Uncompressed VCF input
  --gzvcf FILE                 BGZF/gzip-compressed VCF input
  --bcf FILE                   BCF input
  --input FILE                 Auto-detect VCF, BGZF VCF, gzip VCF, or BCF
                               (vcftools-ng extension)

SITE STATISTICS OUTPUT:
  --freq                       Allele frequencies with allele labels (.frq)
  --freq2                      Frequencies in VCF allele order (.frq)
  --counts                     Allele counts (.frq.count)
  --missing-site               Per-site missingness (.lmiss)
  --site-depth                 Per-site total-depth statistics (.ldepth)
  --site-mean-depth            Per-site mean depth (.ldepth.mean)
  --site-quality               Per-site QUAL value (.lqual)
  --hardy                      Per-site Hardy-Weinberg statistics (.hwe)

INDIVIDUAL STATISTICS OUTPUT:
  --depth                      Mean depth per individual (.idepth)
  --missing-indv               Missingness per individual (.imiss)
  --het                        Heterozygosity/inbreeding per individual (.het)

DIVERSITY AND POPULATION STATISTICS:
  --site-pi                    Nucleotide diversity per site (.sites.pi)
  --window-pi N                Nucleotide diversity in N-bp windows
                               (.windowed.pi)
  --window-pi-step N           Step size for --window-pi
  --TajimaD N                  Tajima's D in N-bp windows (.Tajima.D)
  --weir-fst-pop FILE          Population sample file for Weir-Cockerham FST;
                               repeat once per population
  --fst-window-size N          FST window size in bp (.windowed.weir.fst)
  --fst-window-step N          Step size for windowed FST

LD AND PCA OUTPUT:
  --geno-r2                    Genotype-correlation LD for biallelic sites
                               (.geno.ld)
  --ld-window N                Maximum SNP count in an LD window
  --ld-window-min N            Minimum SNP separation for LD pairs
  --ld-window-bp N             Maximum physical LD distance in bp
  --ld-window-bp-min N         Minimum physical LD distance in bp
  --min-r2 FLOAT               Report LD pairs with r2 >= FLOAT
  --pca                        Normalized-genotype PCA (.pca)
  --pca-no-norm                PCA without genotype normalization (.pca)

RECODE AND FORMAT OUTPUT:
  --recode                     Write uncompressed VCF (PREFIX.recode.vcf)
  --recode-vcf-gz              Write deterministic parallel BGZF VCF
                               (PREFIX.recode.vcf.gz; vcftools-ng extension)
  --recode-bcf                 Write BCF (PREFIX.recode.bcf)
  --recode-INFO-all            Retain all INFO fields in recoded output
  --stdout                     Send plain --recode VCF to stdout

  --recode and --recode-vcf-gz may be combined to write both files in one
  scan. --stdout is valid only with plain --recode. v0.12.3 does not create
  an index for new output; run:
    bcftools index --tbi --threads N PREFIX.recode.vcf.gz

TWO-FILE COMPARISON:
  --diff FILE                  Compare with a second uncompressed VCF
  --gzdiff FILE                Compare with a second compressed VCF
  --diff-bcf FILE              Compare with a second BCF
  --diff-site                  Report site membership differences
  --diff-indv                  Report individual membership differences
  --diff-site-discordance      Report per-site genotype discordance
  --diff-indv-discordance      Report per-individual genotype discordance

CHROMOSOME, POSITION, AND BED FILTERS:
  --chr CHROM                  Keep CHROM; may be repeated
  --not-chr CHROM              Exclude CHROM; may be repeated
  --from-bp POS                Keep positions >= POS; requires exactly one --chr
  --to-bp POS                  Keep positions <= POS; requires exactly one --chr
  --positions FILE             Keep sites listed as CHROM and POS
  --exclude-positions FILE     Exclude sites listed as CHROM and POS
  --bed FILE                   Keep sites overlapping BED intervals
  --exclude-bed FILE           Exclude sites overlapping BED intervals
  --thin N                     Keep sites at least N bp apart

SAMPLE FILTERS:
  --indv SAMPLE                Keep one sample; may be repeated
  --remove-indv SAMPLE         Exclude one sample; may be repeated
  --keep FILE                  Keep samples listed in FILE; may be repeated
  --remove FILE                Exclude samples listed in FILE; may be repeated

ALLELE, QUALITY, DEPTH, AND MISSINGNESS FILTERS:
  --min-alleles N              Keep sites with at least N alleles
  --max-alleles N              Keep sites with at most N alleles
  --remove-indels              Exclude indels
  --keep-only-indels           Keep indels and exclude other variants
  --minQ FLOAT                 Keep sites with QUAL >= FLOAT
  --minGQ FLOAT                Mask genotypes with GQ < FLOAT
  --minDP N                    Mask genotypes with depth < N
  --maxDP N                    Mask genotypes with depth > N
  --min-meanDP FLOAT           Keep sites with mean depth >= FLOAT
  --max-meanDP FLOAT           Keep sites with mean depth <= FLOAT
  --max-missing FLOAT          Keep sites with call rate >= FLOAT
  --max-missing-count N        Keep sites with at most N missing genotypes

FREQUENCY, COUNT, AND HWE FILTERS:
  --maf FLOAT                  Keep sites with MAF >= FLOAT
  --max-maf FLOAT              Keep sites with MAF <= FLOAT
  --mac N                      Keep sites with minor-allele count >= N
  --max-mac N                  Keep sites with minor-allele count <= N
  --hwe FLOAT                  Keep biallelic sites with HWE p >= FLOAT
  --non-ref-af FLOAT           Minimum total non-reference allele frequency
  --max-non-ref-af FLOAT       Maximum total non-reference allele frequency
  --non-ref-af-any FLOAT       Minimum frequency of any non-reference allele
  --max-non-ref-af-any FLOAT   Maximum frequency of every non-reference allele
  --non-ref-ac N               Minimum total non-reference allele count
  --max-non-ref-ac N           Maximum total non-reference allele count
  --non-ref-ac-any N           Minimum count of any non-reference allele
  --max-non-ref-ac-any N       Maximum count of every non-reference allele

VCF FILTER, INFO, AND GENOTYPE FILTERS:
  --keep-filtered FLAG         Keep sites carrying FILTER FLAG; may be repeated
  --remove-filtered FLAG       Exclude sites carrying FILTER FLAG; repeatable
  --remove-filtered-all        Exclude every site whose FILTER is not PASS
  --keep-INFO FLAG             Keep sites carrying INFO flag FLAG; repeatable
  --remove-INFO FLAG           Exclude sites carrying INFO flag FLAG; repeatable
  --remove-filtered-geno FLAG  Mask genotypes carrying FORMAT/FT FLAG
  --remove-filtered-geno-all   Mask every genotype with non-PASS FORMAT/FT

ADAPTIVE INPUT AND INDEX POLICY (default --input-backend auto):
  Plain VCF       Stream at 1-2 threads; aligned byte ranges at 3+ threads.
                  Plain VCF cannot use CSI/TBI.
  BGZF recode     Stream at 1 thread; reuse/build TBI/CSI at 2+ threads.
  BCF recode      Stream for full-file recode even when CSI exists.
  Region query    Reuse/build a valid index for selective BGZF/BCF access.
  Full-scan stats Reuse an existing index from 4 threads; do not build a
                  one-use index.

  Existing .tbi/.csi sidecars are validated independently and never
  overwritten. If a protected sidecar is stale or invalid, automatic mode
  warns and falls back instead of replacing it.

IMPORTANT COMBINATION RULES:
  * --freq and --freq2 write the same artifact and cannot be combined.
  * Diff outputs cannot be combined with single-file outputs.
  * --chr and --not-chr cannot be used together.
  * --bed and --exclude-bed cannot be used together.
  * --remove-indels and --keep-only-indels cannot be used together.
  * Exact BCF recode rejects genotype masking where Original output is corrupt.
  * At least one output option is required.

COMPATIBILITY:
  VCFtools 0.1.17 is the exact-output oracle. Parameters described as
  compatible have real-data complete-file comparison gates. vcftools-ng-only
  extensions include --input, --threads, --input-backend, and
  --recode-vcf-gz. The project does not yet implement every VCFtools option.

TERMINAL COLORS:
  Colors are enabled automatically on an interactive terminal.
  Set NO_COLOR=1 to disable or CLICOLOR_FORCE=1 to force colored help.

Documentation:
  https://github.com/VensinMa/vcftools-ng
  https://github.com/VensinMa/vcftools-ng/blob/master/docs/parameter-compatibility.md
)HELP";
    emit_help(help.str());
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_help();
            std::exit(0);
        } else if (arg == "--version") {
            std::cout << "vcftools-ng " << kVersion << '\n';
            std::exit(0);
        } else if (arg == "--vcf" || arg == "--gzvcf" || arg == "--bcf" ||
                   arg == "--input") {
            if (!options.input.empty()) {
                fail("Only one input file may be specified");
            }
            options.input = require_value(argc, argv, i);
            options.input_bcf = arg == "--bcf";
        } else if (arg == "--out") {
            options.output_prefix = require_value(argc, argv, i);
        } else if (arg == "--threads" || arg == "-t") {
            options.requested_threads =
                parse_unsigned(require_value(argc, argv, i), arg);
            options.threads = options.requested_threads;
            options.threads_explicit = true;
        } else if (arg == "--batch-size") {
            options.batch_size =
                parse_unsigned(require_value(argc, argv, i), arg);
        } else if (arg == "--compat") {
            const auto mode = require_value(argc, argv, i);
            if (mode != "exact") {
                fail("Only --compat exact is implemented");
            }
        } else if (arg == "--input-backend") {
            options.input_backend =
                vcftools_ng::input::parse_backend(
                    require_value(argc, argv, i));
        } else if (arg == "--bcftools") {
            options.bcftools_path =
                require_value(argc, argv, i);
        } else if (arg == "--freq") {
            options.output_freq = true;
        } else if (arg == "--freq2") {
            options.output_freq2 = true;
        } else if (arg == "--counts") {
            options.output_counts = true;
        } else if (arg == "--missing-site") {
            options.output_missing_site = true;
        } else if (arg == "--site-depth") {
            options.output_site_depth = true;
        } else if (arg == "--site-mean-depth") {
            options.output_site_mean_depth = true;
        } else if (arg == "--depth") {
            options.output_individual_depth = true;
        } else if (arg == "--missing-indv") {
            options.output_individual_missingness = true;
        } else if (arg == "--het") {
            options.output_heterozygosity = true;
        } else if (arg == "--hardy") {
            options.output_hardy_weinberg = true;
        } else if (arg == "--site-quality") {
            options.output_site_quality = true;
        } else if (arg == "--site-pi") {
            options.output_site_pi = true;
        } else if (arg == "--window-pi") {
            options.pi_window_size =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--window-pi-step") {
            options.pi_window_step =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--TajimaD") {
            options.tajima_window_size =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--weir-fst-pop") {
            const std::string path = require_value(argc, argv, i);
            options.fst_population_files.push_back(path);
            options.sample_keep_files.push_back(path);
        } else if (arg == "--fst-window-size") {
            options.fst_window_size =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--fst-window-step") {
            options.fst_window_step =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--geno-r2") {
            options.output_genotype_r2 = true;
            options.min_alleles = 2;
            options.max_alleles = 2;
        } else if (arg == "--ld-window") {
            options.ld_snp_window_size =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--ld-window-min") {
            options.ld_snp_window_min =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--ld-window-bp") {
            options.ld_bp_window_size =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--ld-window-bp-min") {
            options.ld_bp_window_min =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--min-r2") {
            options.min_r2 =
                parse_double(require_value(argc, argv, i), arg);
        } else if (arg == "--pca" || arg == "--pca-no-norm") {
            options.output_pca = true;
            options.pca_normalise = arg == "--pca";
            options.min_alleles = 2;
            options.max_alleles = 2;
        } else if (arg == "--recode") {
            options.output_recode = true;
        } else if (arg == "--recode-bcf") {
            options.output_recode_bcf = true;
        } else if (arg == "--recode-vcf-gz") {
            options.output_recode_vcf_gz = true;
        } else if (arg == "--recode-INFO-all") {
            options.recode_info_all = true;
        } else if (arg == "--stdout") {
            options.output_stdout = true;
        } else if (
            arg == "--diff" || arg == "--gzdiff" ||
            arg == "--diff-bcf") {
            if (!options.diff_input.empty()) {
                fail("Only one diff input file may be specified");
            }
            options.diff_input = require_value(argc, argv, i);
        } else if (arg == "--diff-site") {
            options.output_diff_sites_in_files = true;
        } else if (arg == "--diff-indv") {
            options.output_diff_individuals_in_files = true;
        } else if (arg == "--diff-site-discordance") {
            options.output_diff_site_discordance = true;
        } else if (arg == "--diff-indv-discordance") {
            options.output_diff_individual_discordance = true;
        } else if (arg == "--chr") {
            options.chromosomes_to_keep.insert(
                require_value(argc, argv, i));
        } else if (arg == "--not-chr") {
            options.chromosomes_to_exclude.insert(
                require_value(argc, argv, i));
        } else if (arg == "--from-bp") {
            options.start_position =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--to-bp") {
            options.end_position =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--positions") {
            options.positions_file = require_value(argc, argv, i);
        } else if (arg == "--exclude-positions") {
            options.exclude_positions_file =
                require_value(argc, argv, i);
        } else if (arg == "--bed" || arg == "--exclude-bed") {
            if (!options.bed_file.empty()) {
                fail("Multiple --bed/--exclude-bed options cannot be used");
            }
            options.bed_file = require_value(argc, argv, i);
            options.bed_exclude = arg == "--exclude-bed";
        } else if (arg == "--keep-filtered") {
            options.site_filters_to_keep.insert(
                require_value(argc, argv, i));
        } else if (arg == "--remove-filtered") {
            options.site_filters_to_remove.insert(
                require_value(argc, argv, i));
        } else if (arg == "--remove-filtered-all") {
            options.remove_all_filtered_sites = true;
        } else if (arg == "--keep-INFO") {
            options.info_flags_to_keep.insert(
                require_value(argc, argv, i));
        } else if (arg == "--remove-INFO") {
            options.info_flags_to_remove.insert(
                require_value(argc, argv, i));
        } else if (arg == "--indv") {
            options.samples_to_keep.insert(
                require_value(argc, argv, i));
        } else if (arg == "--remove-indv") {
            options.samples_to_exclude.insert(
                require_value(argc, argv, i));
        } else if (arg == "--keep") {
            options.sample_keep_files.push_back(
                require_value(argc, argv, i));
        } else if (arg == "--remove") {
            options.sample_exclude_files.push_back(
                require_value(argc, argv, i));
        } else if (arg == "--min-alleles") {
            options.min_alleles =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--max-alleles") {
            options.max_alleles =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--remove-indels") {
            options.remove_indels = true;
        } else if (arg == "--keep-only-indels") {
            options.keep_only_indels = true;
        } else if (arg == "--minQ") {
            options.min_qual =
                parse_double(require_value(argc, argv, i), arg);
        } else if (arg == "--minGQ") {
            options.min_gq =
                parse_double(require_value(argc, argv, i), arg);
        } else if (arg == "--minDP") {
            options.min_dp = parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--maxDP") {
            options.max_dp = parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--remove-filtered-geno") {
            options.genotype_filters_to_remove.insert(
                require_value(argc, argv, i));
        } else if (arg == "--remove-filtered-geno-all") {
            options.remove_all_filtered_genotypes = true;
        } else if (arg == "--min-meanDP") {
            options.min_mean_dp =
                parse_double(require_value(argc, argv, i), arg);
        } else if (arg == "--max-meanDP") {
            options.max_mean_dp =
                parse_double(require_value(argc, argv, i), arg);
        } else if (arg == "--max-missing") {
            options.min_call_rate =
                parse_double(require_value(argc, argv, i), arg);
        } else if (arg == "--max-missing-count") {
            options.max_missing_count =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--maf") {
            options.min_maf =
                parse_double(require_value(argc, argv, i), arg);
        } else if (arg == "--max-maf") {
            options.max_maf =
                parse_double(require_value(argc, argv, i), arg);
        } else if (arg == "--mac") {
            options.min_mac =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--max-mac") {
            options.max_mac =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--hwe") {
            options.max_alleles = 2;
            options.min_hwe =
                parse_double(require_value(argc, argv, i), arg);
        } else if (arg == "--thin") {
            options.thin_distance =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--non-ref-af") {
            options.min_non_ref_af =
                parse_double(require_value(argc, argv, i), arg);
        } else if (arg == "--max-non-ref-af") {
            options.max_non_ref_af =
                parse_double(require_value(argc, argv, i), arg);
        } else if (arg == "--non-ref-af-any") {
            options.min_non_ref_af_any =
                parse_double(require_value(argc, argv, i), arg);
        } else if (arg == "--max-non-ref-af-any") {
            options.max_non_ref_af_any =
                parse_double(require_value(argc, argv, i), arg);
        } else if (arg == "--non-ref-ac") {
            options.min_non_ref_ac =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--max-non-ref-ac") {
            options.max_non_ref_ac =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--non-ref-ac-any") {
            options.min_non_ref_ac_any =
                parse_int(require_value(argc, argv, i), arg);
        } else if (arg == "--max-non-ref-ac-any") {
            options.max_non_ref_ac_any =
                parse_int(require_value(argc, argv, i), arg);
        } else {
            fail("Unsupported option in phase 1: " + arg);
        }
    }

    if (options.input.empty()) {
        fail("Input file required");
    }
    if (options.threads_explicit) {
        options.threads = std::min(
            options.requested_threads,
            vcftools_ng::input::detect_available_threads().count);
    }
    if (!(options.output_freq || options.output_freq2 ||
          options.output_counts ||
          options.output_missing_site || options.output_site_depth ||
          options.output_site_mean_depth ||
          options.output_individual_depth ||
          options.output_individual_missingness ||
          options.output_heterozygosity ||
          options.output_hardy_weinberg ||
          options.output_site_quality || options.output_site_pi ||
          options.pi_window_size > 0 ||
          options.tajima_window_size > 0 ||
          !options.fst_population_files.empty() ||
          options.output_genotype_r2 ||
          options.output_pca ||
          options.output_recode || options.output_recode_bcf ||
          options.output_recode_vcf_gz ||
          options.output_diff_sites_in_files ||
          options.output_diff_individuals_in_files ||
          options.output_diff_site_discordance ||
          options.output_diff_individual_discordance)) {
        fail("At least one output option is required");
    }
    const bool diff_output =
        options.output_diff_sites_in_files ||
        options.output_diff_individuals_in_files ||
        options.output_diff_site_discordance ||
        options.output_diff_individual_discordance;
    if (diff_output != !options.diff_input.empty()) {
        fail(
            "A diff input and at least one diff output option must be "
            "specified together");
    }
    const bool non_diff_output =
        options.output_freq || options.output_freq2 ||
        options.output_counts ||
        options.output_missing_site || options.output_site_depth ||
        options.output_site_mean_depth ||
        options.output_individual_depth ||
        options.output_individual_missingness ||
        options.output_heterozygosity ||
        options.output_hardy_weinberg ||
        options.output_site_quality || options.output_site_pi ||
        options.pi_window_size > 0 ||
        options.tajima_window_size > 0 ||
        !options.fst_population_files.empty() ||
        options.output_genotype_r2 || options.output_pca ||
        options.output_recode || options.output_recode_bcf ||
        options.output_recode_vcf_gz;
    if (diff_output && non_diff_output) {
        fail("Diff outputs cannot yet be combined with single-file outputs");
    }
    if (options.output_freq && options.output_freq2) {
        fail("--freq and --freq2 write the same artifact and cannot combine");
    }
    if (options.output_recode_bcf &&
        (options.min_gq > 0.0 || options.min_dp > 0 ||
         options.max_dp != std::numeric_limits<int>::max() ||
         options.remove_all_filtered_genotypes ||
         !options.genotype_filters_to_remove.empty())) {
        fail(
            "--recode-bcf exact mode does not support genotype masking: "
            "VCFtools 0.1.17 corrupts this BCF-input oracle path");
    }
    if (options.output_stdout &&
        (!options.output_recode ||
         options.output_recode_bcf ||
         options.output_recode_vcf_gz)) {
        fail(
            "--stdout is supported with plain --recode output only");
    }
    if (options.keep_only_indels && options.remove_indels) {
        fail("Cannot use --keep-only-indels and --remove-indels together");
    }
    if (!options.chromosomes_to_keep.empty() &&
        !options.chromosomes_to_exclude.empty()) {
        fail("Cannot specify chromosomes to keep and to exclude");
    }
    if (options.end_position < options.start_position) {
        fail("End position must be greater than start position");
    }
    if ((options.start_position != -1 ||
         options.end_position != std::numeric_limits<int>::max()) &&
        options.chromosomes_to_keep.size() != 1) {
        fail("A range requires exactly one --chr");
    }
    if (options.max_alleles < options.min_alleles) {
        fail("--max-alleles must not be less than --min-alleles");
    }
    if (options.max_dp < options.min_dp) {
        fail("--maxDP must not be less than --minDP");
    }
    if (options.pi_window_size < 0) {
        fail("Pi Window size must be > 0");
    }
    if (options.tajima_window_size < 0) {
        fail("Tajima D bin size must be > 0");
    }
    if (options.fst_window_size < 0) {
        fail("Fst window size must be > 0");
    }
    if (!options.fst_population_files.empty() &&
        options.fst_population_files.size() < 2) {
        fail("Require at least two populations to estimate Fst");
    }
    if (options.max_mean_dp < options.min_mean_dp) {
        fail("--max-meanDP must not be less than --min-meanDP");
    }
    if (options.max_maf < options.min_maf) {
        fail("--max-maf must not be less than --maf");
    }
    if (options.max_mac < options.min_mac) {
        fail("--max-mac must not be less than --mac");
    }
    if (options.max_non_ref_af < options.min_non_ref_af ||
        options.max_non_ref_af_any < options.min_non_ref_af_any) {
        fail("Maximum non-reference AF must not be less than its minimum");
    }
    if (options.max_non_ref_ac < options.min_non_ref_ac ||
        options.max_non_ref_ac_any < options.min_non_ref_ac_any) {
        fail("Maximum non-reference AC must not be less than its minimum");
    }
    if (options.min_call_rate < 0.0 || options.min_call_rate > 1.0) {
        fail("--max-missing must be between 0 and 1");
    }
    options.input_bcf = input_is_bcf(options.input);
    return options;
}

struct HtsFileDeleter {
    void operator()(htsFile* file) const {
        if (file != nullptr) {
            hts_close(file);
        }
    }
};

struct HeaderDeleter {
    void operator()(bcf_hdr_t* header) const {
        if (header != nullptr) {
            bcf_hdr_destroy(header);
        }
    }
};

struct RecordDeleter {
    void operator()(bcf1_t* record) const {
        if (record != nullptr) {
            bcf_destroy(record);
        }
    }
};

struct IndexDeleter {
    void operator()(hts_idx_t* index) const {
        if (index != nullptr) {
            hts_idx_destroy(index);
        }
    }
};

struct IteratorDeleter {
    void operator()(hts_itr_t* iterator) const {
        if (iterator != nullptr) {
            hts_itr_destroy(iterator);
        }
    }
};

using HtsFilePtr = std::unique_ptr<htsFile, HtsFileDeleter>;
using HeaderPtr = std::unique_ptr<bcf_hdr_t, HeaderDeleter>;
using RecordPtr = std::unique_ptr<bcf1_t, RecordDeleter>;
using IndexPtr = std::unique_ptr<hts_idx_t, IndexDeleter>;
using IteratorPtr = std::unique_ptr<hts_itr_t, IteratorDeleter>;

std::vector<uint8_t> read_raw_bcf_header(const std::string& path) {
    std::unique_ptr<BGZF, decltype(&bgzf_close)> input(
        bgzf_open(path.c_str(), "r"), &bgzf_close);
    if (!input) {
        fail("Could not open BCF input for exact header preservation: " +
             path);
    }
    std::array<uint8_t, 9> prefix{};
    if (bgzf_read(input.get(), prefix.data(), prefix.size()) !=
        static_cast<ssize_t>(prefix.size())) {
        fail("Could not read BCF header prefix: " + path);
    }
    const std::array<uint8_t, 5> bcf_magic{
        'B', 'C', 'F', 2, 2};
    if (!std::equal(
            bcf_magic.begin(), bcf_magic.end(), prefix.begin())) {
        fail(
            "--recode-bcf exact mode requires BCF input so the original "
            "header bytes can be preserved");
    }
    const uint32_t text_length =
        static_cast<uint32_t>(prefix[5]) |
        (static_cast<uint32_t>(prefix[6]) << 8) |
        (static_cast<uint32_t>(prefix[7]) << 16) |
        (static_cast<uint32_t>(prefix[8]) << 24);
    std::vector<uint8_t> header(prefix.begin(), prefix.end());
    header.resize(prefix.size() + text_length);
    std::size_t offset = prefix.size();
    while (offset < header.size()) {
        const ssize_t count = bgzf_read(
            input.get(), header.data() + offset, header.size() - offset);
        if (count <= 0) {
            fail("Truncated BCF header: " + path);
        }
        offset += static_cast<std::size_t>(count);
    }
    return header;
}

std::vector<std::string> split_original_header_field(
    const std::string& value, char delimiter) {
    std::vector<std::string> fields;
    std::istringstream input(value);
    std::string field;
    while (std::getline(input, field, delimiter)) {
        fields.push_back(field);
    }
    return fields;
}

std::optional<std::string> reprint_original_bcf_meta_line(
    const std::string& line) {
    if (line.starts_with("##fileformat=")) {
        return std::nullopt;
    }
    std::string kind;
    if (line.starts_with("##INFO=<")) {
        kind = "INFO";
    } else if (line.starts_with("##FILTER=<")) {
        kind = "FILTER";
    } else if (line.starts_with("##FORMAT=<")) {
        kind = "FORMAT";
    } else if (line.starts_with("##contig=<")) {
        kind = "contig";
    } else {
        return line;
    }

    const std::size_t begin = line.find('<');
    const std::size_t end = line.rfind('>');
    const std::string details =
        line.substr(
            begin + 1,
            end == std::string::npos
                ? std::string::npos
                : end - begin - 1);
    std::string id;
    std::string number;
    std::string type;
    std::string description;
    std::string source;
    std::string version;
    std::string length;
    std::string assembly;
    std::string other;
    for (const auto& token :
         split_original_header_field(details, ',')) {
        const auto pair =
            split_original_header_field(token, '=');
        if (pair.size() < 2) {
            continue;
        }
        const auto& key = pair[0];
        const auto& value = pair[1];
        if (key == "ID") {
            id = value;
        } else if (key == "Number") {
            number = value;
        } else if (key == "Type") {
            type = value;
        } else if (key == "Description") {
            description = value;
        } else if (key == "Source" && kind == "INFO") {
            source = value;
        } else if (key == "Version" && kind == "INFO") {
            version = value;
        } else if (key == "length" && kind == "contig") {
            length = value;
        } else if (key == "assembly" && kind == "contig") {
            assembly = value;
        } else if (key != "IDX" || kind == "contig") {
            if (!other.empty()) {
                other += ',';
            }
            other += token;
        }
    }
    if (kind == "FILTER" && id == "PASS") {
        return std::nullopt;
    }

    std::ostringstream output;
    output << "##" << kind << "=<";
    if (!id.empty()) {
        output << "ID=" << id;
    }
    if (!number.empty()) {
        output << ",Number=" << number;
    }
    if (!type.empty()) {
        output << ",Type=" << type;
    }
    if (!description.empty()) {
        output << ",Description=" << description;
    }
    if (!source.empty()) {
        output << ",Source=" << source;
    }
    if (!version.empty()) {
        output << ",Version=" << version;
    }
    if (!length.empty()) {
        output << ",Length=" << length;
    }
    if (!assembly.empty()) {
        output << ",Assembly=" << assembly;
    }
    if (!other.empty()) {
        output << ',' << other;
    }
    output << '>';
    return output.str();
}

bool input_is_bcf(const std::string& path) {
    HtsFilePtr input(hts_open(path.c_str(), "r"));
    if (!input) {
        fail("Could not inspect input format: " + path);
    }
    const htsFormat* format = hts_get_format(input.get());
    return format != nullptr && format->format == bcf;
}

std::string original_compatible_bcf_vcf_header(
    const std::string& path, bcf_hdr_t* output_header) {
    const auto raw = read_raw_bcf_header(path);
    std::string text(
        reinterpret_cast<const char*>(raw.data() + 9),
        raw.size() - 9);
    const std::size_t terminator = text.find('\0');
    if (terminator != std::string::npos) {
        text.resize(terminator);
    }
    const bool has_idx = text.find("IDX=") != std::string::npos;
    std::ostringstream output;
    std::istringstream lines(text);
    std::string line;
    while (std::getline(lines, line)) {
        if (!line.starts_with("##")) {
            continue;
        }
        if (has_idx) {
            const auto reprinted =
                reprint_original_bcf_meta_line(line);
            if (reprinted.has_value()) {
                output << *reprinted << '\n';
            }
        } else {
            output << line << '\n';
        }
    }
    output << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    const int samples = bcf_hdr_nsamples(output_header);
    if (samples > 0) {
        output << "\tFORMAT";
        for (int sample = 0; sample < samples; ++sample) {
            output << '\t' << output_header->samples[sample];
        }
    }
    output << '\n';
    return output.str();
}

class DeterministicBgzfWriter {
public:
    DeterministicBgzfWriter(
        const std::string& path, unsigned compression_threads)
        : output_(path, std::ios::binary),
          maximum_in_flight_(
              std::max<std::size_t>(4, compression_threads * 2)) {
        if (!output_) {
            fail("Could not open BGZF output: " + path);
        }
        const unsigned worker_count =
            std::max(1u, compression_threads);
        workers_.reserve(worker_count);
        for (unsigned worker = 0; worker < worker_count; ++worker) {
            workers_.emplace_back([this] { compression_worker(); });
        }
        writer_ = std::thread([this] { ordered_writer(); });
    }

    DeterministicBgzfWriter(const DeterministicBgzfWriter&) = delete;
    DeterministicBgzfWriter& operator=(
        const DeterministicBgzfWriter&) = delete;

    ~DeterministicBgzfWriter() {
        try {
            finish();
        } catch (...) {
        }
    }

    void append(const void* data, std::size_t length) {
        rethrow_failure();
        const auto* bytes = static_cast<const uint8_t*>(data);
        while (length > 0) {
            const std::size_t count =
                std::min(length, kBlockSize - current_.size());
            current_.insert(
                current_.end(), bytes, bytes + count);
            bytes += count;
            length -= count;
            if (current_.size() == kBlockSize) {
                enqueue_current();
            }
        }
    }

    void finish() {
        if (finished_) {
            rethrow_failure();
            return;
        }
        finished_ = true;
        if (!current_.empty()) {
            enqueue_current();
        }
        enqueue(std::vector<uint8_t>{});
        {
            std::lock_guard lock(mutex_);
            producer_done_ = true;
            total_jobs_ = next_job_id_;
        }
        work_available_.notify_all();
        completed_available_.notify_all();
        for (auto& worker : workers_) {
            worker.join();
        }
        writer_.join();
        output_.close();
        rethrow_failure();
        if (!output_) {
            fail("Could not finish BGZF output");
        }
    }

private:
    static constexpr std::size_t kBlockSize = 0xff00;
    static constexpr std::size_t kMaximumBlockSize = 0x10000;

    struct Job {
        std::size_t id = 0;
        std::vector<uint8_t> input;
    };

    static void put_u16(uint8_t* output, uint16_t value) {
        output[0] = static_cast<uint8_t>(value);
        output[1] = static_cast<uint8_t>(value >> 8);
    }

    static void put_u32(uint8_t* output, uint32_t value) {
        output[0] = static_cast<uint8_t>(value);
        output[1] = static_cast<uint8_t>(value >> 8);
        output[2] = static_cast<uint8_t>(value >> 16);
        output[3] = static_cast<uint8_t>(value >> 24);
    }

    static std::vector<uint8_t> compress_block(
        const std::vector<uint8_t>& input) {
        static constexpr std::array<uint8_t, 18> header{
            31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0,
            66, 67, 2, 0, 0, 0};
        std::vector<uint8_t> output(kMaximumBlockSize);
        std::copy(header.begin(), header.end(), output.begin());

        z_stream stream{};
        stream.next_in = const_cast<Bytef*>(
            reinterpret_cast<const Bytef*>(input.data()));
        stream.avail_in = static_cast<uInt>(input.size());
        stream.next_out = output.data() + header.size();
        stream.avail_out = static_cast<uInt>(
            output.size() - header.size() - 8);
        if (deflateInit2(
                &stream, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8,
                Z_DEFAULT_STRATEGY) != Z_OK) {
            fail("Could not initialise deterministic BGZF compression");
        }
        const int status = deflate(&stream, Z_FINISH);
        const int close_status = deflateEnd(&stream);
        if (status != Z_STREAM_END || close_status != Z_OK) {
            fail("Could not compress deterministic BGZF block");
        }
        const std::size_t compressed_size =
            header.size() + stream.total_out + 8;
        if (compressed_size > kMaximumBlockSize) {
            fail("Deterministic BGZF block exceeded 64 KiB");
        }
        put_u16(
            output.data() + 16,
            static_cast<uint16_t>(compressed_size - 1));
        const uint32_t crc = crc32(
            crc32(0L, Z_NULL, 0),
            reinterpret_cast<const Bytef*>(input.data()),
            static_cast<uInt>(input.size()));
        put_u32(output.data() + compressed_size - 8, crc);
        put_u32(
            output.data() + compressed_size - 4,
            static_cast<uint32_t>(input.size()));
        output.resize(compressed_size);
        return output;
    }

    void enqueue_current() {
        std::vector<uint8_t> full;
        full.swap(current_);
        current_.reserve(kBlockSize);
        enqueue(std::move(full));
    }

    void enqueue(std::vector<uint8_t> input) {
        std::unique_lock lock(mutex_);
        slot_available_.wait(lock, [&] {
            return failure_ != nullptr ||
                   in_flight_ < maximum_in_flight_;
        });
        if (failure_) {
            std::rethrow_exception(failure_);
        }
        pending_.push_back(Job{next_job_id_++, std::move(input)});
        ++in_flight_;
        lock.unlock();
        work_available_.notify_one();
    }

    void compression_worker() {
        try {
            while (true) {
                Job job;
                {
                    std::unique_lock lock(mutex_);
                    work_available_.wait(lock, [&] {
                        return failure_ != nullptr ||
                               !pending_.empty() || producer_done_;
                    });
                    if (failure_) {
                        return;
                    }
                    if (pending_.empty()) {
                        if (producer_done_) {
                            return;
                        }
                        continue;
                    }
                    job = std::move(pending_.front());
                    pending_.pop_front();
                }
                auto compressed = compress_block(job.input);
                {
                    std::lock_guard lock(mutex_);
                    completed_.emplace(
                        job.id, std::move(compressed));
                }
                completed_available_.notify_one();
            }
        } catch (...) {
            record_failure(std::current_exception());
        }
    }

    void ordered_writer() {
        try {
            std::size_t next = 0;
            while (true) {
                std::vector<uint8_t> block;
                {
                    std::unique_lock lock(mutex_);
                    completed_available_.wait(lock, [&] {
                        return failure_ != nullptr ||
                               completed_.contains(next) ||
                               (producer_done_ && next == total_jobs_);
                    });
                    if (failure_) {
                        return;
                    }
                    if (producer_done_ && next == total_jobs_) {
                        return;
                    }
                    auto found = completed_.find(next);
                    block = std::move(found->second);
                    completed_.erase(found);
                    --in_flight_;
                    ++next;
                }
                slot_available_.notify_one();
                output_.write(
                    reinterpret_cast<const char*>(block.data()),
                    static_cast<std::streamsize>(block.size()));
                if (!output_) {
                    fail("Could not write deterministic BGZF block");
                }
            }
        } catch (...) {
            record_failure(std::current_exception());
        }
    }

    void record_failure(std::exception_ptr failure) {
        {
            std::lock_guard lock(mutex_);
            if (!failure_) {
                failure_ = failure;
            }
        }
        work_available_.notify_all();
        completed_available_.notify_all();
        slot_available_.notify_all();
    }

    void rethrow_failure() {
        std::lock_guard lock(mutex_);
        if (failure_) {
            std::rethrow_exception(failure_);
        }
    }

    std::ofstream output_;
    std::vector<uint8_t> current_;
    std::size_t maximum_in_flight_ = 0;
    std::mutex mutex_;
    std::condition_variable work_available_;
    std::condition_variable completed_available_;
    std::condition_variable slot_available_;
    std::deque<Job> pending_;
    std::map<std::size_t, std::vector<uint8_t>> completed_;
    std::vector<std::thread> workers_;
    std::thread writer_;
    std::exception_ptr failure_;
    std::size_t next_job_id_ = 0;
    std::size_t total_jobs_ = 0;
    std::size_t in_flight_ = 0;
    bool producer_done_ = false;
    bool finished_ = false;
};

class ExactBcfOutput {
public:
    ExactBcfOutput(
        const std::string& input_path, const std::string& output_path,
        unsigned compression_threads)
        : compressor_(output_path, compression_threads) {
        const auto header = read_raw_bcf_header(input_path);
        compressor_.append(header.data(), header.size());

        int descriptors[2]{-1, -1};
        if (pipe(descriptors) != 0) {
            fail("Could not create BCF record serialization pipe");
        }
        read_descriptor_ = descriptors[0];
        hFILE* stream = hdopen(descriptors[1], "w");
        if (!stream) {
            close(descriptors[0]);
            close(descriptors[1]);
            fail("Could not open BCF record serialization stream");
        }
        record_output_.reset(
            hts_hopen(stream, "vcftools-ng-records.bcf", "wbu"));
        if (!record_output_) {
            const int close_status = hclose(stream);
            (void)close_status;
            close(descriptors[0]);
            fail("Could not initialise uncompressed BCF serializer");
        }
        reader_ = std::thread([this] { read_serialized_records(); });
    }

    ExactBcfOutput(const ExactBcfOutput&) = delete;
    ExactBcfOutput& operator=(const ExactBcfOutput&) = delete;

    ~ExactBcfOutput() {
        try {
            finish();
        } catch (...) {
        }
    }

    void write(bcf_hdr_t* header, bcf1_t* record) {
        if (bcf_write(record_output_.get(), header, record) != 0) {
            fail("Could not serialize BCF record");
        }
    }

    void finish() {
        if (finished_) {
            if (reader_failure_) {
                std::rethrow_exception(reader_failure_);
            }
            return;
        }
        finished_ = true;
        record_output_.reset();
        reader_.join();
        if (reader_failure_) {
            std::rethrow_exception(reader_failure_);
        }
        compressor_.finish();
    }

private:
    void read_serialized_records() {
        try {
            std::array<uint8_t, 1 << 16> buffer{};
            while (true) {
                const ssize_t count =
                    ::read(read_descriptor_, buffer.data(), buffer.size());
                if (count == 0) {
                    break;
                }
                if (count < 0) {
                    fail("Could not read serialized BCF records");
                }
                compressor_.append(
                    buffer.data(), static_cast<std::size_t>(count));
            }
            close(read_descriptor_);
            read_descriptor_ = -1;
        } catch (...) {
            reader_failure_ = std::current_exception();
            if (read_descriptor_ >= 0) {
                close(read_descriptor_);
                read_descriptor_ = -1;
            }
        }
    }

    DeterministicBgzfWriter compressor_;
    HtsFilePtr record_output_;
    int read_descriptor_ = -1;
    std::thread reader_;
    std::exception_ptr reader_failure_;
    bool finished_ = false;
};

class SampleSelection {
public:
    SampleSelection(const Options& options, bcf_hdr_t* input_header)
        : active_(
              !options.samples_to_keep.empty() ||
              !options.samples_to_exclude.empty() ||
              !options.sample_keep_files.empty() ||
              !options.sample_exclude_files.empty()) {
        population_count_ = options.fst_population_files.size();
        std::set<std::string> keep = options.samples_to_keep;
        std::set<std::string> exclude = options.samples_to_exclude;
        load_sample_files(options.sample_keep_files, keep);
        load_sample_files(options.sample_exclude_files, exclude);

        const bool has_keep_filter =
            !options.samples_to_keep.empty() ||
            !options.sample_keep_files.empty();
        const int sample_count = bcf_hdr_nsamples(input_header);
        population_memberships_.resize(sample_count);
        for (std::size_t population = 0;
             population < options.fst_population_files.size();
             ++population) {
            std::set<std::string> members;
            load_sample_files(
                {options.fst_population_files[population]}, members);
            for (int sample = 0; sample < sample_count; ++sample) {
                if (members.contains(input_header->samples[sample])) {
                    population_memberships_[sample].push_back(
                        static_cast<int>(population));
                }
            }
        }
        for (int sample = 0; sample < sample_count; ++sample) {
            const std::string name = input_header->samples[sample];
            if (has_keep_filter && !keep.contains(name)) {
                continue;
            }
            if (exclude.contains(name)) {
                continue;
            }
            indices_.push_back(sample);
        }
        if (indices_.empty()) {
            fail("Sample filters removed all individuals");
        }

        if (active_) {
            std::vector<char*> names;
            names.reserve(indices_.size());
            for (const int sample : indices_) {
                names.push_back(input_header->samples[sample]);
            }
            imap_.resize(indices_.size());
            output_header_.reset(bcf_hdr_subset(
                input_header, static_cast<int>(names.size()),
                names.data(), imap_.data()));
            if (!output_header_ ||
                bcf_hdr_nsamples(output_header_.get()) !=
                    static_cast<int>(indices_.size())) {
                fail("Could not construct sample-subset VCF header");
            }
        }
    }

    const std::vector<int>& indices() const {
        return indices_;
    }

    int count() const {
        return static_cast<int>(indices_.size());
    }

    bool active() const {
        return active_;
    }

    std::size_t population_count() const {
        return population_count_;
    }

    const std::vector<int>& populations_for_sample(int sample) const {
        return population_memberships_[sample];
    }

    bcf_hdr_t* output_header(bcf_hdr_t* input_header) const {
        return active_ ? output_header_.get() : input_header;
    }

    void subset_record(
        bcf1_t* record, bcf_hdr_t* worker_output_header) const {
        if (!active_) {
            return;
        }
        if (bcf_subset(
                worker_output_header, record,
                static_cast<int>(imap_.size()),
                const_cast<int*>(imap_.data())) != 0) {
            fail("Could not subset VCF record samples");
        }
    }

private:
    static void load_sample_files(
        const std::vector<std::string>& paths,
        std::set<std::string>& names) {
        for (const auto& path : paths) {
            std::ifstream input(path);
            if (!input) {
                fail("Could not open sample file: " + path);
            }
            std::string line;
            while (std::getline(input, line)) {
                std::istringstream fields(line);
                std::string name;
                if (fields >> name) {
                    names.insert(name);
                }
            }
        }
    }

    bool active_;
    std::vector<int> indices_;
    std::vector<std::vector<int>> population_memberships_;
    std::size_t population_count_ = 0;
    std::vector<int> imap_;
    HeaderPtr output_header_;
};

class SiteSelector {
public:
    SiteSelector(const Options& options, bcf_hdr_t* header)
        : options_(options),
          active_(
              !options.chromosomes_to_keep.empty() ||
              !options.chromosomes_to_exclude.empty() ||
              options.start_position != -1 ||
              options.end_position !=
                  std::numeric_limits<int>::max() ||
              !options.positions_file.empty() ||
              !options.exclude_positions_file.empty() ||
              !options.bed_file.empty()) {
        load_positions(
            options.positions_file, header, positions_to_keep_);
        load_positions(
            options.exclude_positions_file, header,
            positions_to_exclude_);
        load_bed(options.bed_file, header);
    }

    bool keep(bcf1_t* record, bcf_hdr_t* header) const {
        if (!active_) {
            return true;
        }
        if (record->rid < 0) {
            return false;
        }
        const std::string chromosome =
            bcf_hdr_id2name(header, record->rid);
        const int position = record->pos + 1;

        if (!options_.chromosomes_to_keep.empty() &&
            options_.chromosomes_to_keep.find(chromosome) ==
                options_.chromosomes_to_keep.end()) {
            return false;
        }
        if (!options_.chromosomes_to_exclude.empty() &&
            options_.chromosomes_to_exclude.find(chromosome) !=
                options_.chromosomes_to_exclude.end()) {
            return false;
        }
        if (options_.start_position != -1 &&
            position < options_.start_position) {
            return false;
        }
        if (options_.end_position != std::numeric_limits<int>::max() &&
            position > options_.end_position) {
            return false;
        }

        const int reference_length = static_cast<int>(
            std::max<hts_pos_t>(1, record->rlen));
        if (!options_.positions_file.empty()) {
            bool found = false;
            for (int offset = 0; offset < reference_length; ++offset) {
                if (positions_to_keep_.contains(
                        make_key(record->rid, position + offset))) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                return false;
            }
        }
        if (!options_.exclude_positions_file.empty()) {
            for (int offset = 0; offset < reference_length; ++offset) {
                if (positions_to_exclude_.contains(
                        make_key(record->rid, position + offset))) {
                    return false;
                }
            }
        }
        if (!options_.bed_file.empty()) {
            bcf_unpack(record, BCF_UN_STR);
            int variant_end = position;
            for (int allele = 0; allele < record->n_allele; ++allele) {
                variant_end = std::max(
                    variant_end,
                    position +
                        static_cast<int>(
                            std::string(record->d.allele[allele]).size()) -
                        1);
            }
            const bool overlaps = overlaps_bed(
                record->rid, position, variant_end);
            if (overlaps == options_.bed_exclude) {
                return false;
            }
        }
        return true;
    }

private:
    static uint64_t make_key(int rid, int position) {
        return (static_cast<uint64_t>(static_cast<uint32_t>(rid)) << 32) |
               static_cast<uint32_t>(position);
    }

    static void load_positions(
        const std::string& path, bcf_hdr_t* header,
        std::unordered_set<uint64_t>& output) {
        if (path.empty()) {
            return;
        }
        HtsFilePtr input(hts_open(path.c_str(), "r"));
        if (!input) {
            fail("Could not open positions file: " + path);
        }
        kstring_t line{0, 0, nullptr};
        while (hts_getline(input.get(), '\n', &line) >= 0) {
            if (line.l == 0 || line.s[0] == '#') {
                continue;
            }
            std::istringstream fields(
                std::string(line.s, line.l));
            std::string chromosome;
            int position = 0;
            if (!(fields >> chromosome >> position)) {
                std::free(line.s);
                fail("Invalid positions line in " + path);
            }
            const int rid = bcf_hdr_name2id(header, chromosome.c_str());
            if (rid >= 0) {
                output.insert(make_key(rid, position));
            }
        }
        std::free(line.s);
    }

    void load_bed(const std::string& path, bcf_hdr_t* header) {
        if (path.empty()) {
            return;
        }
        std::ifstream input(path);
        if (!input) {
            fail("Could not open BED file: " + path);
        }
        bed_intervals_.resize(header->n[BCF_DT_CTG]);
        bed_prefix_max_end_.resize(header->n[BCF_DT_CTG]);

        std::string line;
        std::getline(input, line);
        while (std::getline(input, line)) {
            std::istringstream fields(line);
            std::string chromosome;
            int start = 0;
            int end = 0;
            if (!(fields >> chromosome >> start >> end)) {
                continue;
            }
            const int rid = bcf_hdr_name2id(header, chromosome.c_str());
            if (rid >= 0) {
                bed_intervals_[rid].emplace_back(start, end);
            }
        }

        for (std::size_t rid = 0; rid < bed_intervals_.size(); ++rid) {
            auto& intervals = bed_intervals_[rid];
            std::sort(intervals.begin(), intervals.end());
            auto& prefix = bed_prefix_max_end_[rid];
            prefix.reserve(intervals.size());
            int maximum_end = std::numeric_limits<int>::min();
            for (const auto& interval : intervals) {
                maximum_end = std::max(maximum_end, interval.second);
                prefix.push_back(maximum_end);
            }
        }
    }

    bool overlaps_bed(int rid, int variant_start, int variant_end) const {
        if (rid < 0 ||
            static_cast<std::size_t>(rid) >= bed_intervals_.size()) {
            return false;
        }
        const auto& intervals = bed_intervals_[rid];
        const auto first_not_before_end = std::lower_bound(
            intervals.begin(), intervals.end(), variant_end,
            [](const std::pair<int, int>& interval, int end) {
                return interval.first < end;
            });
        if (first_not_before_end == intervals.begin()) {
            return false;
        }
        const std::size_t last_candidate =
            static_cast<std::size_t>(
                std::distance(intervals.begin(), first_not_before_end) - 1);
        return bed_prefix_max_end_[rid][last_candidate] >= variant_start;
    }

    const Options& options_;
    bool active_;
    std::unordered_set<uint64_t> positions_to_keep_;
    std::unordered_set<uint64_t> positions_to_exclude_;
    std::vector<std::vector<std::pair<int, int>>> bed_intervals_;
    std::vector<std::vector<int>> bed_prefix_max_end_;
};

struct Scratch {
    int32_t* gt = nullptr;
    int gt_capacity = 0;
    int32_t* dp = nullptr;
    int dp_capacity = 0;
    int32_t* gq = nullptr;
    int gq_capacity = 0;
    char** ft = nullptr;
    int ft_capacity = 0;
    std::vector<uint8_t> genotype_filtered;
    std::vector<std::vector<uint32_t>> population_homozygotes;
    std::vector<std::vector<uint32_t>> population_heterozygotes;
    kstring_t recode_buffer{0, 0, nullptr};

    ~Scratch() {
        std::free(gt);
        std::free(dp);
        std::free(gq);
        if (ft != nullptr) {
            std::free(ft[0]);
            std::free(ft);
        }
        std::free(recode_buffer.s);
    }
};

struct SiteResult {
    bool kept = false;
    std::string chrom;
    int64_t pos = 0;
    std::vector<std::string> alleles;
    std::vector<uint32_t> allele_counts;
    uint32_t non_missing_chromosomes = 0;
    uint32_t n_data = 0;
    uint32_t n_genotype_filtered = 0;
    uint32_t n_missing = 0;
    uint32_t sum_depth = 0;
    uint32_t sumsq_depth = 0;
    uint32_t depth_count = 0;
    double quality = -1.0;
    bool fully_diploid = true;
    bool heterozygosity_eligible = false;
    double expected_homozygosity = 0.0;
    uint32_t hwe_hom_ref = 0;
    uint32_t hwe_hets = 0;
    uint32_t hwe_hom_alt = 0;
    bool fst_eligible = false;
    double fst_sum_a = 0.0;
    double fst_sum_all = 0.0;
    double fst = std::numeric_limits<double>::quiet_NaN();
    bcf1_t* output_record = nullptr;
    std::string recode_line;
};

enum IndividualState : uint8_t {
    kGenotypeFiltered = 1U << 0,
    kFirstAlleleMissing = 1U << 1,
    kCompleteDiploidGenotype = 1U << 2,
    kObservedHomozygote = 1U << 3,
};

struct BatchAnalysisPayload {
    std::size_t sample_count = 0;
    std::vector<uint8_t> individual_state;
    std::vector<int32_t> individual_depth;
    std::vector<int8_t> genotype_dosage;

    void reset(std::size_t row_count, std::size_t selected_samples,
               const Options& options) {
        sample_count = selected_samples;
        const std::size_t cell_count = row_count * sample_count;
        if (options.output_individual_missingness ||
            options.output_heterozygosity) {
            individual_state.assign(cell_count, 0);
        }
        if (options.output_individual_depth) {
            individual_depth.assign(cell_count, -1);
        }
        if (options.output_genotype_r2 || options.output_pca) {
            genotype_dosage.assign(cell_count, -1);
        }
    }

    uint8_t* state_row(std::size_t row) {
        return individual_state.empty()
                   ? nullptr
                   : individual_state.data() + row * sample_count;
    }

    int32_t* depth_row(std::size_t row) {
        return individual_depth.empty()
                   ? nullptr
                   : individual_depth.data() + row * sample_count;
    }

    int8_t* dosage_row(std::size_t row) {
        return genotype_dosage.empty()
                   ? nullptr
                   : genotype_dosage.data() + row * sample_count;
    }
};

class ThinSelector {
public:
    explicit ThinSelector(int minimum_distance)
        : minimum_distance_(minimum_distance) {}

    bool keep(const SiteResult& result) {
        if (minimum_distance_ < 1) {
            return true;
        }
        if (result.chrom == chromosome_ &&
            result.pos - position_ < minimum_distance_) {
            return false;
        }
        chromosome_ = result.chrom;
        position_ = result.pos;
        return true;
    }

private:
    int minimum_distance_;
    std::string chromosome_;
    int64_t position_ = 0;
};

int legacy_scalar_value(const int32_t* values, int count, int sample,
                        int sample_count) {
    if (count <= 0) {
        return -1;
    }
    const int stride = count / sample_count;
    if (stride <= 0) {
        return -1;
    }
    const int32_t value = values[sample * stride];
    if (value == bcf_int32_missing || value == bcf_int32_vector_end) {
        return -1;
    }
    return value;
}

bool passes_site_filter_flags(const Options& options, bcf_hdr_t* header,
                               bcf1_t* record) {
    if (!options.remove_all_filtered_sites &&
        options.site_filters_to_keep.empty() &&
        options.site_filters_to_remove.empty()) {
        return true;
    }

    if (!options.site_filters_to_keep.empty()) {
        bool keep = false;
        for (int index = 0; index < record->d.n_flt; ++index) {
            const char* name = bcf_hdr_int2id(
                header, BCF_DT_ID, record->d.flt[index]);
            if (name != nullptr &&
                options.site_filters_to_keep.contains(name)) {
                keep = true;
                break;
            }
        }
        if (!keep) {
            return false;
        }
    }

    if (record->d.n_flt > 0) {
        const char* first = bcf_hdr_int2id(
            header, BCF_DT_ID, record->d.flt[0]);
        if (first != nullptr && std::string_view(first) == "PASS") {
            return true;
        }
    }
    if (options.remove_all_filtered_sites && record->d.n_flt > 0) {
        return false;
    }
    for (int index = 0; index < record->d.n_flt; ++index) {
        const char* name = bcf_hdr_int2id(
            header, BCF_DT_ID, record->d.flt[index]);
        if (name != nullptr &&
            options.site_filters_to_remove.contains(name)) {
            return false;
        }
    }
    return true;
}

bool info_flag_present(bcf_hdr_t* header, bcf1_t* record,
                       const std::string& name) {
    return bcf_get_info(header, record, name.c_str()) != nullptr;
}

bool passes_info_flags(const Options& options, bcf_hdr_t* header,
                       bcf1_t* record) {
    if (!options.info_flags_to_keep.empty()) {
        bool keep = false;
        for (const auto& name : options.info_flags_to_keep) {
            if (info_flag_present(header, record, name)) {
                keep = true;
                break;
            }
        }
        if (!keep) {
            return false;
        }
    }
    for (const auto& name : options.info_flags_to_remove) {
        if (info_flag_present(header, record, name)) {
            return false;
        }
    }
    return true;
}

void validate_info_flag_filters(const Options& options, bcf_hdr_t* header) {
    const auto validate = [header](const std::string& name) {
        const int id = bcf_hdr_id2int(
            header, BCF_DT_ID, name.c_str());
        if (!bcf_hdr_idinfo_exists(header, BCF_HL_INFO, id) ||
            bcf_hdr_id2type(header, BCF_HL_INFO, id) != BCF_HT_FLAG) {
            fail("Using INFO flag filtering on non-Flag type " + name +
                 " will not work correctly");
        }
    };
    for (const auto& name : options.info_flags_to_keep) {
        validate(name);
    }
    for (const auto& name : options.info_flags_to_remove) {
        validate(name);
    }
}

bool genotype_filter_matches(
    std::string_view value, const std::set<std::string>& filters,
    bool remove_all) {
    std::string_view first;
    bool have_first = false;
    std::size_t start = 0;
    while (start <= value.size()) {
        const std::size_t end = value.find(';', start);
        const std::string_view token = value.substr(
            start, end == std::string_view::npos
                       ? value.size() - start
                       : end - start);
        if (!token.empty() && token != ".") {
            if (!have_first) {
                first = token;
                have_first = true;
            }
            if (!remove_all && filters.contains(std::string(token))) {
                return true;
            }
        }
        if (end == std::string_view::npos) {
            break;
        }
        start = end + 1;
    }
    return remove_all && have_first && first != "PASS";
}

bool is_legacy_indel(const bcf1_t* record) {
    const std::size_t ref_length = std::string(record->d.allele[0]).size();
    if (ref_length != 1) {
        return true;
    }
    for (int allele = 1; allele < record->n_allele; ++allele) {
        if (std::string(record->d.allele[allele]).size() != ref_length) {
            return true;
        }
    }
    return false;
}

struct ExactHweProbabilities {
    double two_sided = 1.0;
    double heterozygote_deficit = 1.0;
    double heterozygote_excess = 1.0;
};

ExactHweProbabilities exact_hwe_probabilities(
    int observed_hets, int observed_hom_ref, int observed_hom_alt) {
    ExactHweProbabilities result;
    if (observed_hom_ref + observed_hom_alt + observed_hets == 0) {
        return result;
    }

    const int observed_hom_common =
        std::max(observed_hom_ref, observed_hom_alt);
    const int observed_hom_rare =
        std::min(observed_hom_ref, observed_hom_alt);
    const int rare_copies = 2 * observed_hom_rare + observed_hets;
    const int genotypes =
        observed_hets + observed_hom_common + observed_hom_rare;

    std::vector<double> heterozygote_probabilities(
        static_cast<std::size_t>(rare_copies + 1), 0.0);
    int midpoint =
        rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
    if ((rare_copies & 1) ^ (midpoint & 1)) {
        ++midpoint;
    }

    int current_hets = midpoint;
    int current_hom_rare = (rare_copies - midpoint) / 2;
    int current_hom_common =
        genotypes - current_hets - current_hom_rare;
    heterozygote_probabilities[midpoint] = 1.0;
    double sum = 1.0;

    for (current_hets = midpoint; current_hets > 1; current_hets -= 2) {
        heterozygote_probabilities[current_hets - 2] =
            heterozygote_probabilities[current_hets] * current_hets *
            (current_hets - 1.0) /
            (4.0 * (current_hom_rare + 1.0) *
             (current_hom_common + 1.0));
        sum += heterozygote_probabilities[current_hets - 2];
        ++current_hom_rare;
        ++current_hom_common;
    }

    current_hets = midpoint;
    current_hom_rare = (rare_copies - midpoint) / 2;
    current_hom_common =
        genotypes - current_hets - current_hom_rare;
    for (current_hets = midpoint; current_hets <= rare_copies - 2;
         current_hets += 2) {
        heterozygote_probabilities[current_hets + 2] =
            heterozygote_probabilities[current_hets] * 4.0 *
            current_hom_rare * current_hom_common /
            ((current_hets + 2.0) * (current_hets + 1.0));
        sum += heterozygote_probabilities[current_hets + 2];
        --current_hom_rare;
        --current_hom_common;
    }

    for (double& probability : heterozygote_probabilities) {
        probability /= sum;
    }

    result.heterozygote_excess =
        heterozygote_probabilities[observed_hets];
    for (int index = observed_hets + 1; index <= rare_copies; ++index) {
        result.heterozygote_excess +=
            heterozygote_probabilities[index];
    }
    result.heterozygote_deficit =
        heterozygote_probabilities[observed_hets];
    for (int index = observed_hets - 1; index >= 0; --index) {
        result.heterozygote_deficit +=
            heterozygote_probabilities[index];
    }

    result.two_sided = 0.0;
    for (const double probability : heterozygote_probabilities) {
        if (probability <=
            heterozygote_probabilities[observed_hets]) {
            result.two_sided += probability;
        }
    }
    result.two_sided = std::min(result.two_sided, 1.0);
    return result;
}

double exact_hwe_pvalue(int observed_hets, int observed_hom_ref,
                        int observed_hom_alt) {
    return exact_hwe_probabilities(
               observed_hets, observed_hom_ref, observed_hom_alt)
        .two_sided;
}

struct FstContribution {
    double sum_a = 0.0;
    double sum_all = 0.0;
    double fst = std::numeric_limits<double>::quiet_NaN();
};

FstContribution calculate_fst(
    const std::vector<std::vector<uint32_t>>& homozygotes,
    const std::vector<std::vector<uint32_t>>& heterozygotes,
    std::size_t allele_count) {
    const std::size_t population_count = homozygotes.size();
    std::vector<double> n(population_count, 0.0);
    std::vector<std::vector<double>> p(
        population_count, std::vector<double>(allele_count, 0.0));
    std::vector<double> pbar(allele_count, 0.0);
    std::vector<double> hbar(allele_count, 0.0);
    std::vector<double> ssqr(allele_count, 0.0);
    double sum_nsqr = 0.0;

    for (std::size_t population = 0;
         population < population_count; ++population) {
        for (std::size_t allele = 0; allele < allele_count; ++allele) {
            n[population] +=
                homozygotes[population][allele] +
                0.5 * heterozygotes[population][allele];
            p[population][allele] =
                heterozygotes[population][allele] +
                2 * homozygotes[population][allele];
            pbar[allele] += p[population][allele];
            hbar[allele] +=
                heterozygotes[population][allele];
        }
        for (std::size_t allele = 0; allele < allele_count; ++allele) {
            p[population][allele] /= (2.0 * n[population]);
        }
        sum_nsqr += n[population] * n[population];
    }

    const double n_sum =
        std::accumulate(n.begin(), n.end(), 0.0);
    const double nbar = n_sum / population_count;
    for (std::size_t allele = 0; allele < allele_count; ++allele) {
        pbar[allele] /= (n_sum * 2.0);
        hbar[allele] /= n_sum;
    }
    for (std::size_t allele = 0; allele < allele_count; ++allele) {
        for (std::size_t population = 0;
             population < population_count; ++population) {
            ssqr[allele] +=
                n[population] *
                (p[population][allele] - pbar[allele]) *
                (p[population][allele] - pbar[allele]);
        }
        ssqr[allele] /=
            ((population_count - 1.0) * nbar);
    }
    const double nc =
        (n_sum - (sum_nsqr / n_sum)) /
        (population_count - 1.0);
    const double r = static_cast<double>(population_count);

    FstContribution result;
    for (std::size_t allele = 0; allele < allele_count; ++allele) {
        const double a =
            (ssqr[allele] -
             (pbar[allele] * (1.0 - pbar[allele]) -
              (((r - 1.0) * ssqr[allele]) / r) -
              (hbar[allele] / 4.0)) /
                 (nbar - 1.0)) *
            nbar / nc;
        const double b =
            (pbar[allele] * (1.0 - pbar[allele]) -
             (ssqr[allele] * (r - 1.0) / r) -
             hbar[allele] *
                 (((2.0 * nbar) - 1.0) / (4.0 * nbar))) *
            nbar / (nbar - 1.0);
        const double c = hbar[allele] / 2.0;
        if (!std::isnan(a) && !std::isnan(b) &&
            !std::isnan(c)) {
            result.sum_a += a;
            result.sum_all += a + b + c;
        }
    }
    result.fst = result.sum_a / result.sum_all;
    return result;
}

extern "C" void dgeev_(
    char* jobvl, char* jobvr, int* n, double* a, int* lda,
    double* wr, double* wi, double* vl, int* ldvl, double* vr,
    int* ldvr, double* work, int* lwork, int* info);

struct EigenResult {
    std::vector<double> real;
    std::vector<double> imaginary;
    std::vector<double> vectors;
};

EigenResult legacy_dgeev(
    const std::vector<double>& row_major_matrix, int dimension) {
    EigenResult result{
        std::vector<double>(dimension),
        std::vector<double>(dimension),
        std::vector<double>(
            static_cast<std::size_t>(dimension) * dimension)};
    std::vector<double> column_major(
        static_cast<std::size_t>(dimension) * dimension);
    for (int row = 0; row < dimension; ++row) {
        for (int column = 0; column < dimension; ++column) {
            column_major[row + column * dimension] =
                row_major_matrix[
                    static_cast<std::size_t>(row) * dimension +
                    column];
        }
    }
    char jobvl = 'N';
    char jobvr = 'V';
    int leading_dimension = dimension;
    int left_leading_dimension = dimension;
    int right_leading_dimension = dimension;
    std::vector<double> left_vectors(
        static_cast<std::size_t>(dimension) * dimension);
    std::vector<double> work(
        static_cast<std::size_t>(4) * dimension);
    int work_size = 4 * dimension;
    int info = 0;
    dgeev_(
        &jobvl, &jobvr, &dimension, column_major.data(),
        &leading_dimension, result.real.data(),
        result.imaginary.data(), left_vectors.data(),
        &left_leading_dimension, result.vectors.data(),
        &right_leading_dimension, work.data(), &work_size, &info);
    if (info != 0) {
        fail("LAPACK dgeev failed with status " + std::to_string(info));
    }

    std::vector<double> magnitude(dimension);
    for (int index = 0; index < dimension; ++index) {
        magnitude[index] =
            result.real[index] * result.real[index] +
            result.imaginary[index] * result.imaginary[index];
    }
    for (int pass = 0; pass < dimension; ++pass) {
        for (int index = 0; index < dimension - 1; ++index) {
            if (std::fabs(magnitude[index]) <
                std::fabs(magnitude[index + 1])) {
                std::swap(magnitude[index], magnitude[index + 1]);
                std::swap(result.real[index], result.real[index + 1]);
                std::swap(
                    result.imaginary[index],
                    result.imaginary[index + 1]);
                for (int row = 0; row < dimension; ++row) {
                    std::swap(
                        result.vectors[
                            row + index * dimension],
                        result.vectors[
                            row + (index + 1) * dimension]);
                }
            }
        }
    }
    return result;
}

void mask_filtered_genotypes(bcf_hdr_t* header, bcf1_t* record,
                             const std::vector<uint8_t>& genotype_filtered,
                             int32_t* genotypes, int genotype_count,
                             int sample_count, int max_ploidy) {
    if (genotype_count <= 0 || max_ploidy <= 0) {
        return;
    }
    for (int sample = 0; sample < sample_count; ++sample) {
        if (!genotype_filtered[sample]) {
            continue;
        }
        const int base = sample * max_ploidy;
        for (int copy = 0; copy < max_ploidy; ++copy) {
            if (genotypes[base + copy] == bcf_int32_vector_end) {
                break;
            }
            genotypes[base + copy] = bcf_gt_missing;
        }
    }
    if (bcf_update_genotypes(
            header, record, genotypes, genotype_count) != 0) {
        fail("Could not mask filtered genotypes");
    }
}

void remove_info_column(std::string& line) {
    std::size_t field_start = 0;
    for (int field = 0; field < 7; ++field) {
        field_start = line.find('\t', field_start);
        if (field_start == std::string::npos) {
            fail("Could not locate INFO column in recoded VCF record");
        }
        ++field_start;
    }
    const std::size_t field_end = line.find('\t', field_start);
    if (field_end == std::string::npos) {
        const std::size_t newline = line.find('\n', field_start);
        line.replace(
            field_start,
            (newline == std::string::npos ? line.size() : newline) -
                field_start,
            ".");
    } else {
        line.replace(field_start, field_end - field_start, ".");
    }
}

void sort_filter_column(std::string& line) {
    std::size_t field_start = 0;
    for (int field = 0; field < 6; ++field) {
        field_start = line.find('\t', field_start);
        if (field_start == std::string::npos) {
            fail("Could not locate FILTER column in recoded VCF record");
        }
        ++field_start;
    }
    const std::size_t field_end = line.find('\t', field_start);
    if (field_end == std::string::npos) {
        fail("Could not locate end of FILTER column");
    }
    const std::string_view filter(
        line.data() + field_start, field_end - field_start);
    if (filter == "." || filter.find(';') == std::string_view::npos) {
        return;
    }

    std::vector<std::string> names;
    std::size_t start = 0;
    while (start <= filter.size()) {
        const std::size_t end = filter.find(';', start);
        names.emplace_back(filter.substr(
            start, end == std::string_view::npos
                       ? filter.size() - start
                       : end - start));
        if (end == std::string_view::npos) {
            break;
        }
        start = end + 1;
    }
    std::sort(names.begin(), names.end());
    std::string sorted = names.front();
    for (std::size_t index = 1; index < names.size(); ++index) {
        sorted.push_back(';');
        sorted += names[index];
    }
    line.replace(field_start, field_end - field_start, sorted);
}

void apply_original_bcf_string_format_quirk(
    bcf_hdr_t* header, bcf1_t* record, std::string& line) {
    bcf_unpack(record, BCF_UN_FMT);
    std::vector<bool> generic_string(record->n_fmt, false);
    bool any = false;
    for (int index = 0; index < record->n_fmt; ++index) {
        const int id = record->d.fmt[index].id;
        const char* name =
            bcf_hdr_int2id(header, BCF_DT_ID, id);
        if (name == nullptr ||
            std::strcmp(name, "GT") == 0 ||
            std::strcmp(name, "FT") == 0) {
            continue;
        }
        if (bcf_hdr_id2type(header, BCF_HL_FMT, id) ==
            BCF_HT_STR) {
            generic_string[index] = true;
            any = true;
        }
    }
    if (!any) {
        return;
    }

    std::size_t sample_start = 0;
    for (int column = 0; column < 9; ++column) {
        sample_start = line.find('\t', sample_start);
        if (sample_start == std::string::npos) {
            return;
        }
        ++sample_start;
    }
    struct Replacement {
        std::size_t begin;
        std::size_t length;
        char value;
    };
    std::vector<Replacement> replacements;
    replacements.reserve(
        static_cast<std::size_t>(bcf_hdr_nsamples(header)));
    std::size_t field_start = sample_start;
    int sample = 0;
    int format_index = 0;
    while (field_start < line.size() &&
           sample < bcf_hdr_nsamples(header)) {
        const std::size_t field_end =
            line.find_first_of(":\t\n", field_start);
        const std::size_t end =
            field_end == std::string::npos
                ? line.size()
                : field_end;
        if (format_index < record->n_fmt &&
            generic_string[format_index]) {
            const bcf_fmt_t& format =
                record->d.fmt[format_index];
            char last = '\0';
            if (format.p != nullptr && format.size > 0) {
                last = reinterpret_cast<const char*>(
                    format.p)[
                    static_cast<std::size_t>(sample) *
                        format.size +
                    format.size - 1];
            }
            replacements.push_back(
                {field_start, end - field_start, last});
        }
        if (field_end == std::string::npos ||
            line[field_end] == '\n') {
            break;
        }
        if (line[field_end] == ':') {
            ++format_index;
        } else {
            ++sample;
            format_index = 0;
        }
        field_start = field_end + 1;
    }
    for (auto replacement = replacements.rbegin();
         replacement != replacements.rend(); ++replacement) {
        line.replace(
            replacement->begin, replacement->length,
            1, replacement->value);
    }
}

void apply_original_bcf_missing_gt_quirk(
    bcf_hdr_t* header, bcf1_t* record,
    const std::vector<uint8_t>& genotype_filtered,
    const std::vector<int>& selected_samples,
    std::string& line) {
    bcf_unpack(record, BCF_UN_FMT);
    int genotype_index = -1;
    for (int index = 0; index < record->n_fmt; ++index) {
        const char* name = bcf_hdr_int2id(
            header, BCF_DT_ID, record->d.fmt[index].id);
        if (name != nullptr && std::strcmp(name, "GT") == 0) {
            genotype_index = index;
            break;
        }
    }
    if (genotype_index < 0) {
        return;
    }

    std::size_t sample_start = 0;
    for (int column = 0; column < 9; ++column) {
        sample_start = line.find('\t', sample_start);
        if (sample_start == std::string::npos) {
            return;
        }
        ++sample_start;
    }
    std::vector<std::size_t> missing_alleles;
    std::size_t field_start = sample_start;
    int sample = 0;
    int format_index = 0;
    while (field_start < line.size() &&
           sample < static_cast<int>(selected_samples.size())) {
        const std::size_t field_end =
            line.find_first_of(":\t\n", field_start);
        const std::size_t end =
            field_end == std::string::npos
                ? line.size()
                : field_end;
        if (format_index == genotype_index &&
            !genotype_filtered[selected_samples[sample]] &&
            end - field_start > 1) {
            for (std::size_t position = field_start;
                 position < end; ++position) {
                if (line[position] == '.') {
                    missing_alleles.push_back(position);
                }
            }
        }
        if (field_end == std::string::npos ||
            line[field_end] == '\n') {
            break;
        }
        if (line[field_end] == ':') {
            ++format_index;
        } else {
            ++sample;
            format_index = 0;
        }
        field_start = field_end + 1;
    }
    for (auto position = missing_alleles.rbegin();
         position != missing_alleles.rend(); ++position) {
        line.replace(*position, 1, "-1");
    }
}

std::string format_recode_line(const Options& options, bcf_hdr_t* header,
                               bcf1_t* record, Scratch& scratch,
                               const std::vector<uint8_t>& genotype_filtered,
                               const std::vector<int>& selected_samples) {
    scratch.recode_buffer.l = 0;
    if (vcf_format1(header, record, &scratch.recode_buffer) != 0) {
        fail("Could not format recoded VCF record");
    }
    std::string line(
        scratch.recode_buffer.s, scratch.recode_buffer.l);
    if (record->d.n_flt > 1) {
        sort_filter_column(line);
    }
    if (options.input_bcf) {
        apply_original_bcf_string_format_quirk(
            header, record, line);
        apply_original_bcf_missing_gt_quirk(
            header, record, genotype_filtered,
            selected_samples, line);
    }
    if (!options.recode_info_all) {
        remove_info_column(line);
    }
    if (line.empty() || line.back() != '\n') {
        line.push_back('\n');
    }
    return line;
}

void clear_info_fields(bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_INFO);
    std::vector<int> ids;
    ids.reserve(record->n_info);
    for (int index = 0; index < record->n_info; ++index) {
        ids.push_back(record->d.info[index].key);
    }
    for (const int id : ids) {
        const char* key = bcf_hdr_int2id(header, BCF_DT_ID, id);
        if (key == nullptr) {
            continue;
        }
        const int type =
            bcf_hdr_id2type(header, BCF_HL_INFO, id);
        if (bcf_update_info(
                header, record, key, nullptr, 0, type) != 0) {
            fail("Could not remove INFO field " + std::string(key));
        }
    }
}

SiteResult process_site(const Options& options,
                        const SiteSelector& selector,
                        const SampleSelection& samples,
                        bcf_hdr_t* header,
                        bcf_hdr_t* worker_output_header,
                        bcf1_t* record,
                        Scratch& scratch,
                        uint8_t* individual_state,
                        int32_t* individual_depth,
                        int8_t* genotype_dosage) {
    SiteResult result;
    if (!selector.keep(record, header)) {
        return result;
    }
    bcf_unpack(record, BCF_UN_STR | BCF_UN_FLT);

    if (!passes_site_filter_flags(options, header, record) ||
        !passes_info_flags(options, header, record)) {
        return result;
    }

    if (record->n_allele < options.min_alleles ||
        record->n_allele > options.max_alleles) {
        return result;
    }

    const bool indel = is_legacy_indel(record);
    if ((options.remove_indels && indel) ||
        (options.keep_only_indels && !indel)) {
        return result;
    }

    const double quality =
        bcf_float_is_missing(record->qual) ? -1.0 : record->qual;
    if (options.min_qual >= 0.0 && quality < options.min_qual) {
        return result;
    }

    const int sample_count = bcf_hdr_nsamples(header);
    const bool frequency_filter_active =
        options.min_maf > 0.0 || options.max_maf < 1.0 ||
        options.min_non_ref_af > 0.0 ||
        options.max_non_ref_af < 1.0 ||
        options.min_non_ref_af_any > 0.0 ||
        options.max_non_ref_af_any < 1.0 ||
        options.min_non_ref_ac > 0 ||
        options.max_non_ref_ac != std::numeric_limits<int>::max() ||
        options.min_non_ref_ac_any > 0 ||
        options.max_non_ref_ac_any !=
            std::numeric_limits<int>::max() ||
        options.min_mac > 0 ||
        options.max_mac != std::numeric_limits<int>::max() ||
        options.min_hwe > 0.0 ||
        options.min_call_rate > 0.0 ||
        options.max_missing_count != std::numeric_limits<int>::max();
    const bool genotype_filter_active =
        options.min_gq > 0.0 || options.min_dp > 0 ||
        options.max_dp != std::numeric_limits<int>::max() ||
        options.remove_all_filtered_genotypes ||
        !options.genotype_filters_to_remove.empty();
    const bool need_gt =
                         options.output_freq ||
                         options.output_freq2 ||
                         options.output_counts ||
                         options.output_missing_site ||
                         options.output_individual_missingness ||
                         options.output_heterozygosity ||
                         options.output_hardy_weinberg ||
                         options.output_site_pi ||
                         options.pi_window_size > 0 ||
                         options.tajima_window_size > 0 ||
                         !options.fst_population_files.empty() ||
                         options.output_genotype_r2 ||
                         options.output_pca ||
                         frequency_filter_active ||
                         ((options.output_recode ||
                           options.output_recode_bcf ||
                           options.output_recode_vcf_gz) &&
                          genotype_filter_active);
    const bool need_dp = options.output_site_depth ||
                         options.output_site_mean_depth ||
                         options.output_individual_depth ||
                         options.min_mean_dp > 0.0 ||
                         options.max_mean_dp !=
                             std::numeric_limits<double>::max() ||
                         options.min_dp > 0 ||
                         options.max_dp != std::numeric_limits<int>::max();

    int gt_count = -1;
    int max_ploidy = 0;
    if (need_gt) {
        gt_count = bcf_get_genotypes(
            header, record, &scratch.gt, &scratch.gt_capacity);
        if (gt_count <= 0 || sample_count == 0 ||
            gt_count % sample_count != 0) {
            fail("GT field required at " +
                 std::string(bcf_hdr_id2name(header, record->rid)) + ":" +
                 std::to_string(record->pos + 1));
        }
        max_ploidy = gt_count / sample_count;
    }

    int dp_count = -1;
    if (need_dp) {
        dp_count = bcf_get_format_int32(
            header, record, "DP", &scratch.dp, &scratch.dp_capacity);
    }

    int gq_count = -1;
    if (options.min_gq > 0.0) {
        gq_count = bcf_get_format_int32(
            header, record, "GQ", &scratch.gq, &scratch.gq_capacity);
    }

    int ft_count = -1;
    if (options.remove_all_filtered_genotypes ||
        !options.genotype_filters_to_remove.empty()) {
        ft_count = bcf_get_format_string(
            header, record, "FT", &scratch.ft, &scratch.ft_capacity);
    }

    if (options.min_mean_dp > 0.0 ||
        options.max_mean_dp != std::numeric_limits<double>::max()) {
        double raw_sum = 0.0;
        for (const int sample : samples.indices()) {
            const int depth =
                legacy_scalar_value(scratch.dp, dp_count, sample, sample_count);
            if (depth >= 0) {
                raw_sum += depth;
            }
        }
        const double raw_mean = raw_sum / samples.count();
        if (raw_mean < options.min_mean_dp ||
            raw_mean > options.max_mean_dp) {
            return result;
        }
    }

    result.allele_counts.assign(record->n_allele, 0);
    if (samples.population_count() > 0) {
        scratch.population_homozygotes.resize(
            samples.population_count());
        scratch.population_heterozygotes.resize(
            samples.population_count());
        for (std::size_t population = 0;
             population < samples.population_count(); ++population) {
            scratch.population_homozygotes[population].assign(
                record->n_allele, 0);
            scratch.population_heterozygotes[population].assign(
                record->n_allele, 0);
        }
    }
    scratch.genotype_filtered.assign(sample_count, 0);
    auto& genotype_filtered = scratch.genotype_filtered;

    for (const int sample : samples.indices()) {
        if (options.min_gq > 0.0 && gq_count > 0) {
            const int quality_value = legacy_scalar_value(
                scratch.gq, gq_count, sample, sample_count);
            if (quality_value < options.min_gq) {
                genotype_filtered[sample] = true;
            }
        }

        if ((options.min_dp > 0 ||
             options.max_dp != std::numeric_limits<int>::max()) &&
            dp_count > 0) {
            const int depth =
                legacy_scalar_value(scratch.dp, dp_count, sample, sample_count);
            if (depth < options.min_dp || depth > options.max_dp) {
                genotype_filtered[sample] = true;
            }
        }
        if (ft_count > sample && scratch.ft != nullptr &&
            genotype_filter_matches(
                scratch.ft[sample], options.genotype_filters_to_remove,
                options.remove_all_filtered_genotypes)) {
            genotype_filtered[sample] = true;
        }
    }

    uint32_t total_ploidy_for_filter = 0;
    int hwe_hom_ref = 0;
    int hwe_hets = 0;
    int hwe_hom_alt = 0;
    result.fully_diploid = true;
    if (need_gt) {
        for (std::size_t selected_index = 0;
             selected_index < samples.indices().size(); ++selected_index) {
            const int sample = samples.indices()[selected_index];
            const int base = sample * max_ploidy;
            int actual_ploidy = 0;
            for (int copy = 0; copy < max_ploidy; ++copy) {
                if (scratch.gt[base + copy] == bcf_int32_vector_end) {
                    break;
                }
                ++actual_ploidy;
            }
            total_ploidy_for_filter += actual_ploidy;

            if (genotype_filtered[sample]) {
                ++result.n_genotype_filtered;
                if (individual_state != nullptr) {
                    individual_state[selected_index] |=
                        kGenotypeFiltered;
                }
                continue;
            }
            if (actual_ploidy != 2) {
                result.fully_diploid = false;
            }

            const int legacy_copies = std::min(actual_ploidy, 2);
            result.n_data += legacy_copies;
            for (int copy = 0; copy < legacy_copies; ++copy) {
                const int32_t encoded = scratch.gt[base + copy];
                if (bcf_gt_is_missing(encoded)) {
                    ++result.n_missing;
                    continue;
                }
                const int allele = bcf_gt_allele(encoded);
                if (allele >= 0 && allele < record->n_allele) {
                    ++result.allele_counts[allele];
                    ++result.non_missing_chromosomes;
                }
            }

            if (individual_state != nullptr) {
                if (actual_ploidy == 0 ||
                    bcf_gt_is_missing(scratch.gt[base])) {
                    individual_state[selected_index] |=
                        kFirstAlleleMissing;
                }
            }

            if (actual_ploidy >= 2) {
                const int32_t first = scratch.gt[base];
                const int32_t second = scratch.gt[base + 1];
                if (!bcf_gt_is_missing(first) &&
                    !bcf_gt_is_missing(second)) {
                    const int first_allele = bcf_gt_allele(first);
                    const int second_allele = bcf_gt_allele(second);
                    if (genotype_dosage != nullptr) {
                        genotype_dosage[selected_index] =
                            static_cast<int8_t>(
                                (first_allele == 0) +
                                (second_allele == 0));
                    }
                    if (individual_state != nullptr) {
                        individual_state[selected_index] |=
                            kCompleteDiploidGenotype;
                        if (first_allele == second_allele) {
                            individual_state[selected_index] |=
                                kObservedHomozygote;
                        }
                    }
                    if (record->n_allele == 2 &&
                        (options.min_hwe > 0.0 ||
                         options.output_hardy_weinberg ||
                         options.output_heterozygosity)) {
                        if (first_allele != second_allele) {
                            ++hwe_hets;
                        } else if (first_allele == 0) {
                            ++hwe_hom_ref;
                        } else if (first_allele == 1) {
                            ++hwe_hom_alt;
                        } else {
                            fail("Unknown allele in biallelic genotype");
                        }
                    }
                    for (const int population :
                         samples.populations_for_sample(sample)) {
                        if (first_allele == second_allele) {
                            ++scratch.population_homozygotes
                                  [population][first_allele];
                        } else {
                            ++scratch.population_heterozygotes
                                  [population][first_allele];
                            ++scratch.population_heterozygotes
                                  [population][second_allele];
                        }
                    }
                }
            }
        }
    }

    if (frequency_filter_active) {
        if (options.min_call_rate > 0.0) {
            const double call_rate =
                result.non_missing_chromosomes /
                static_cast<double>(total_ploidy_for_filter);
            if (call_rate < options.min_call_rate) {
                return result;
            }
        }
        if (options.max_missing_count != std::numeric_limits<int>::max() &&
            total_ploidy_for_filter - result.non_missing_chromosomes >
                static_cast<uint32_t>(options.max_missing_count)) {
            return result;
        }
        double maf = std::numeric_limits<double>::max();
        int failed_non_ref_af_any = 0;
        for (std::size_t allele = 0;
             allele < result.allele_counts.size(); ++allele) {
            const double frequency =
                result.allele_counts[allele] /
                static_cast<double>(result.non_missing_chromosomes);
            maf = std::min(maf, std::min(frequency, 1.0 - frequency));
            if (allele > 0 &&
                (frequency < options.min_non_ref_af ||
                 frequency > options.max_non_ref_af)) {
                return result;
            }
            if (allele > 0 &&
                (frequency < options.min_non_ref_af_any ||
                 frequency > options.max_non_ref_af_any)) {
                ++failed_non_ref_af_any;
            }
        }
        if ((options.min_non_ref_af > 0.0 ||
             options.max_non_ref_af < 1.0) &&
            failed_non_ref_af_any ==
                static_cast<int>(result.allele_counts.size()) - 1) {
            return result;
        }
        if (options.min_maf > 0.0 || options.max_maf < 1.0) {
            if (maf < options.min_maf || maf > options.max_maf) {
                return result;
            }
        }
        int failed_non_ref_ac_any = 0;
        for (std::size_t allele = 1;
             allele < result.allele_counts.size(); ++allele) {
            const int64_t count = result.allele_counts[allele];
            if (count < options.min_non_ref_ac ||
                count > options.max_non_ref_ac) {
                return result;
            }
            if (count < options.min_non_ref_ac_any ||
                count > options.max_non_ref_ac_any) {
                ++failed_non_ref_ac_any;
            }
        }
        if ((options.min_non_ref_ac_any > 0 ||
             options.max_non_ref_ac_any !=
                 std::numeric_limits<int>::max()) &&
            failed_non_ref_ac_any ==
                static_cast<int>(result.allele_counts.size()) - 1) {
            return result;
        }
        if (options.min_mac > 0 ||
            options.max_mac != std::numeric_limits<int>::max()) {
            if (record->n_allele <= 1 && options.min_mac > 0) {
                return result;
            }
            const uint32_t mac = *std::min_element(
                result.allele_counts.begin(), result.allele_counts.end());
            if (static_cast<int64_t>(mac) < options.min_mac ||
                static_cast<int64_t>(mac) > options.max_mac) {
                return result;
            }
        }
        if (options.min_hwe > 0.0 &&
            exact_hwe_pvalue(
                hwe_hets, hwe_hom_ref, hwe_hom_alt) < options.min_hwe) {
            return result;
        }
    }

    if (need_dp) {
        for (std::size_t selected_index = 0;
             selected_index < samples.indices().size(); ++selected_index) {
            const int sample = samples.indices()[selected_index];
            if (genotype_filtered[sample]) {
                continue;
            }
            const int depth =
                legacy_scalar_value(scratch.dp, dp_count, sample, sample_count);
            if (depth >= 0) {
                if (individual_depth != nullptr) {
                    individual_depth[selected_index] = depth;
                }
                result.sum_depth += static_cast<uint32_t>(depth);
                result.sumsq_depth +=
                    static_cast<uint32_t>(depth * depth);
                ++result.depth_count;
            }
        }
    }

    result.chrom = bcf_hdr_id2name(header, record->rid);
    result.pos = record->pos + 1;
    result.quality = quality;
    result.hwe_hom_ref = static_cast<uint32_t>(hwe_hom_ref);
    result.hwe_hets = static_cast<uint32_t>(hwe_hets);
    result.hwe_hom_alt = static_cast<uint32_t>(hwe_hom_alt);
    if (samples.population_count() > 0 && result.fully_diploid) {
        const auto fst = calculate_fst(
            scratch.population_homozygotes,
            scratch.population_heterozygotes,
            static_cast<std::size_t>(record->n_allele));
        result.fst_eligible = true;
        result.fst_sum_a = fst.sum_a;
        result.fst_sum_all = fst.sum_all;
        result.fst = fst.fst;
    }
    if (record->n_allele == 2 && result.fully_diploid &&
        result.non_missing_chromosomes > 0) {
        const double alt_frequency =
            result.allele_counts[1] /
            static_cast<double>(result.non_missing_chromosomes);
        if (alt_frequency > std::numeric_limits<double>::epsilon() &&
            1.0 - alt_frequency >
                std::numeric_limits<double>::epsilon()) {
            result.heterozygosity_eligible = true;
            result.expected_homozygosity =
                1.0 -
                (2.0 * alt_frequency * (1.0 - alt_frequency) *
                 (result.non_missing_chromosomes /
                  (result.non_missing_chromosomes - 1.0)));
        }
    }
    result.alleles.reserve(record->n_allele);
    for (int allele = 0; allele < record->n_allele; ++allele) {
        result.alleles.emplace_back(record->d.allele[allele]);
    }
    if (options.output_recode || options.output_recode_bcf ||
        options.output_recode_vcf_gz) {
        if (genotype_filter_active) {
            mask_filtered_genotypes(
                header, record, genotype_filtered, scratch.gt, gt_count,
                sample_count, max_ploidy);
        }
        samples.subset_record(record, worker_output_header);
        if (!options.recode_info_all &&
            (options.output_recode_bcf ||
             options.output_recode_vcf_gz)) {
            clear_info_fields(
                worker_output_header, record);
        }
        if (options.output_recode ||
            options.output_recode_vcf_gz) {
            result.recode_line =
                format_recode_line(
                    options, worker_output_header, record,
                    scratch, genotype_filtered,
                    samples.indices());
        }
        if (options.output_recode_bcf) {
            result.output_record = record;
        }
    }
    result.kept = true;
    return result;
}

class Outputs {
public:
    Outputs(
        const Options& options, bcf_hdr_t* header,
        bool sample_selection_active)
        : options_(options),
          output_header_(header),
          selected_chromosome_count_(
              static_cast<uint64_t>(bcf_hdr_nsamples(header)) * 2) {
        const int selected_samples = bcf_hdr_nsamples(header);
        sample_names_.reserve(selected_samples);
        for (int sample = 0; sample < selected_samples; ++sample) {
            sample_names_.emplace_back(header->samples[sample]);
        }
        individual_depth_sum_.assign(selected_samples, 0.0);
        individual_depth_count_.assign(selected_samples, 0);
        individual_data_count_.assign(selected_samples, 0);
        individual_filtered_count_.assign(selected_samples, 0);
        individual_missing_count_.assign(selected_samples, 0);
        heterozygosity_sites_.assign(selected_samples, 0);
        observed_homozygotes_.assign(selected_samples, 0);
        expected_homozygotes_.assign(selected_samples, 0.0);

        if (options.output_freq || options.output_freq2) {
            freq_ = open(".frq");
            freq_
                << (options.output_freq
                        ? "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\n"
                        : "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\n");
        }
        if (options.output_counts) {
            counts_ = open(".frq.count");
            counts_ << "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\n";
        }
        if (options.output_missing_site) {
            missing_ = open(".lmiss");
            missing_
                << "CHR\tPOS\tN_DATA\tN_GENOTYPE_FILTERED\tN_MISS\tF_MISS\n";
        }
        if (options.output_site_depth) {
            depth_ = open(".ldepth");
            depth_ << "CHROM\tPOS\tSUM_DEPTH\tSUMSQ_DEPTH\n";
        }
        if (options.output_site_mean_depth) {
            mean_depth_ = open(".ldepth.mean");
            mean_depth_ << "CHROM\tPOS\tMEAN_DEPTH\tVAR_DEPTH\n";
        }
        if (options.output_individual_depth) {
            individual_depth_ = open(".idepth");
            individual_depth_ << "INDV\tN_SITES\tMEAN_DEPTH\n";
        }
        if (options.output_individual_missingness) {
            individual_missingness_ = open(".imiss");
            individual_missingness_
                << "INDV\tN_DATA\tN_GENOTYPES_FILTERED\tN_MISS\tF_MISS\n";
        }
        if (options.output_heterozygosity) {
            heterozygosity_ = open(".het");
            heterozygosity_
                << "INDV\tO(HOM)\tE(HOM)\tN_SITES\tF\n";
        }
        if (options.output_hardy_weinberg) {
            hardy_weinberg_ = open(".hwe");
            hardy_weinberg_
                << "CHR\tPOS\tOBS(HOM1/HET/HOM2)"
                << "\tE(HOM1/HET/HOM2)\tChiSq_HWE\tP_HWE"
                << "\tP_HET_DEFICIT\tP_HET_EXCESS\n";
        }
        if (options.output_site_quality) {
            site_quality_ = open(".lqual");
            site_quality_ << "CHROM\tPOS\tQUAL\n";
        }
        if (options.output_site_pi) {
            site_pi_ = open(".sites.pi");
            site_pi_ << "CHROM\tPOS\tPI\n";
        }
        if (options.pi_window_size > 0) {
            window_pi_ = open(".windowed.pi");
            window_pi_
                << "CHROM\tBIN_START\tBIN_END\tN_VARIANTS"
                << "\tN_MONOMORPHIC\tPI\n";
            pi_window_step_ =
                options.pi_window_step <= 0 ||
                        options.pi_window_step > options.pi_window_size
                    ? options.pi_window_size
                    : options.pi_window_step;
        }
        if (options.tajima_window_size > 0) {
            tajima_ = open(".Tajima.D");
            tajima_ << "CHROM\tBIN_START\tN_SNPS\tTajimaD\n";
            initialise_tajima_constants();
        }
        if (!options.fst_population_files.empty()) {
            if (options.fst_window_size > 0) {
                window_fst_ = open(".windowed.weir.fst");
                window_fst_
                    << "CHROM\tBIN_START\tBIN_END\tN_VARIANTS"
                    << "\tWEIGHTED_FST\tMEAN_FST\n";
                fst_window_step_ =
                    options.fst_window_step <= 0 ||
                            options.fst_window_step >
                                options.fst_window_size
                        ? options.fst_window_size
                        : options.fst_window_step;
            } else {
                site_fst_ = open(".weir.fst");
                site_fst_
                    << "CHROM\tPOS\tWEIR_AND_COCKERHAM_FST\n";
            }
        }
        if (options.output_genotype_r2) {
            genotype_ld_ = open(".geno.ld");
            genotype_ld_ << "CHR\tPOS1\tPOS2\tN_INDV\tR^2\n";
        }
        if (options.output_pca) {
            pca_ = open(".pca");
        }
        if (options.output_recode ||
            options.output_recode_vcf_gz) {
            std::string bcf_header_text;
            kstring_t formatted_header{0, 0, nullptr};
            const char* header_data = nullptr;
            std::size_t header_size = 0;
            if (input_is_bcf(options.input)) {
                bcf_header_text =
                    original_compatible_bcf_vcf_header(
                        options.input, header);
                header_data = bcf_header_text.data();
                header_size = bcf_header_text.size();
            } else {
                if (bcf_hdr_format(
                        header, 0, &formatted_header) != 0) {
                    std::free(formatted_header.s);
                    fail("Could not format VCF header");
                }
                header_data = formatted_header.s;
                header_size = formatted_header.l;
            }
            if (options.output_recode) {
                if (options.output_stdout) {
                    recode_stream_ = &std::cout;
                } else {
                    recode_ = open(".recode.vcf");
                    recode_stream_ = &recode_;
                }
                recode_stream_->write(
                    header_data,
                    static_cast<std::streamsize>(header_size));
            }
            if (options.output_recode_vcf_gz) {
                recode_vcf_gz_ =
                    std::make_unique<DeterministicBgzfWriter>(
                        options.output_prefix + ".recode.vcf.gz",
                        options.threads);
                recode_vcf_gz_->append(header_data, header_size);
            }
            std::free(formatted_header.s);
            if (options.output_recode && !*recode_stream_) {
                fail("Could not write recoded VCF header");
            }
        }
        if (options.output_recode_bcf) {
            if (sample_selection_active) {
                fail(
                    "--recode-bcf exact mode does not yet support sample "
                    "subsetting because it changes the raw BCF header");
            }
            recode_bcf_ = std::make_unique<ExactBcfOutput>(
                options.input,
                options.output_prefix + ".recode.bcf",
                options.threads);
        }
    }

    void write(const SiteResult& result,
               const uint8_t* individual_state,
               const int32_t* individual_depth,
               const int8_t* genotype_dosage) {
        if (!result.kept) {
            return;
        }
        if (options_.output_freq || options_.output_freq2) {
            freq_ << result.chrom << '\t' << result.pos << '\t'
                  << result.alleles.size() << '\t'
                  << result.non_missing_chromosomes;
            for (std::size_t i = 0; i < result.alleles.size(); ++i) {
                const double frequency =
                    result.allele_counts[i] /
                    static_cast<double>(result.non_missing_chromosomes);
                freq_ << '\t';
                if (options_.output_freq) {
                    freq_ << result.alleles[i] << ':';
                }
                freq_ << frequency;
            }
            freq_ << '\n';
        }
        if (options_.output_counts) {
            counts_ << result.chrom << '\t' << result.pos << '\t'
                    << result.alleles.size() << '\t'
                    << result.non_missing_chromosomes;
            for (std::size_t i = 0; i < result.alleles.size(); ++i) {
                counts_ << '\t' << result.alleles[i] << ':'
                        << result.allele_counts[i];
            }
            counts_ << '\n';
        }
        if (options_.output_missing_site) {
            missing_ << result.chrom << '\t' << result.pos << '\t'
                     << result.n_data << '\t'
                     << result.n_genotype_filtered << '\t'
                     << result.n_missing << '\t'
                     << result.n_missing /
                            static_cast<double>(result.n_data)
                     << '\n';
        }
        if (options_.output_site_depth) {
            depth_ << result.chrom << '\t' << result.pos << '\t'
                   << result.sum_depth << '\t' << result.sumsq_depth << '\n';
        }
        if (options_.output_site_mean_depth) {
            const double mean =
                result.sum_depth / static_cast<double>(result.depth_count);
            const double variance =
                ((result.sumsq_depth /
                      static_cast<double>(result.depth_count)) -
                 (mean * mean)) *
                result.depth_count /
                static_cast<double>(result.depth_count - 1);
            mean_depth_ << result.chrom << '\t' << result.pos << '\t' << mean
                        << '\t' << variance << '\n';
        }
        if (options_.output_site_quality) {
            site_quality_ << result.chrom << '\t' << result.pos << '\t'
                          << result.quality << '\n';
        }
        if ((options_.output_site_pi || options_.pi_window_size > 0) &&
            result.fully_diploid) {
            uint64_t mismatches = 0;
            for (const uint32_t allele_count : result.allele_counts) {
                mismatches +=
                    static_cast<uint64_t>(allele_count) *
                    (result.non_missing_chromosomes - allele_count);
            }
            if (options_.output_site_pi) {
                const uint64_t pairs =
                    static_cast<uint64_t>(
                        result.non_missing_chromosomes) *
                    (result.non_missing_chromosomes - 1);
                site_pi_ << result.chrom << '\t' << result.pos << '\t'
                         << mismatches / static_cast<double>(pairs)
                         << '\n';
            }
            if (options_.pi_window_size > 0 && mismatches > 0) {
                consume_window_pi(result, mismatches);
            }
        }
        if (options_.tajima_window_size > 0 &&
            result.alleles.size() == 2 && result.fully_diploid) {
            consume_tajima(result);
        }
        if (result.fst_eligible) {
            if (options_.fst_window_size > 0) {
                if (!std::isnan(result.fst)) {
                    consume_window_fst(result);
                }
            } else if (!options_.fst_population_files.empty()) {
                site_fst_ << result.chrom << '\t' << result.pos << '\t'
                          << result.fst << '\n';
            }
        }
        if (options_.output_genotype_r2 || options_.output_pca) {
            LdSite site;
            site.chrom = result.chrom;
            site.pos = result.pos;
            site.dosage.assign(
                genotype_dosage,
                genotype_dosage + sample_names_.size());
            if (options_.output_pca) {
                if (!result.fully_diploid) {
                    fail(
                        "PCA only works for fully diploid sites. "
                        "Non-diploid site at " +
                        result.chrom + ":" +
                        std::to_string(result.pos));
                }
                site.alt_frequency =
                    result.allele_counts[1] /
                    static_cast<double>(
                        result.non_missing_chromosomes);
                site.pca_eligible =
                    site.alt_frequency >
                        std::numeric_limits<double>::epsilon() &&
                    site.alt_frequency <
                        1.0 -
                            std::numeric_limits<double>::epsilon();
                if (site.pca_eligible &&
                    std::find(
                        site.dosage.begin(), site.dosage.end(),
                        static_cast<int8_t>(-1)) !=
                        site.dosage.end()) {
                    fail(
                        "PCA exact mode requires complete genotypes; "
                        "use --max-missing 1");
                }
            }
            ld_sites_.push_back(std::move(site));
        }
        if (options_.output_hardy_weinberg &&
            result.alleles.size() == 2 && result.fully_diploid) {
            const double total =
                result.hwe_hom_ref + result.hwe_hets +
                result.hwe_hom_alt;
            const double frequency =
                result.allele_counts[0] /
                static_cast<double>(result.non_missing_chromosomes);
            const double expected_hom_ref =
                frequency * frequency * total;
            const double expected_hets =
                2.0 * frequency * (1.0 - frequency) * total;
            const double expected_hom_alt =
                (1.0 - frequency) * (1.0 - frequency) * total;
            const double chi_square =
                ((result.hwe_hom_ref - expected_hom_ref) *
                 (result.hwe_hom_ref - expected_hom_ref)) /
                    expected_hom_ref +
                ((result.hwe_hets - expected_hets) *
                 (result.hwe_hets - expected_hets)) /
                    expected_hets +
                ((result.hwe_hom_alt - expected_hom_alt) *
                 (result.hwe_hom_alt - expected_hom_alt)) /
                    expected_hom_alt;
            const auto probabilities = exact_hwe_probabilities(
                result.hwe_hets, result.hwe_hom_ref,
                result.hwe_hom_alt);
            hardy_weinberg_ << result.chrom << '\t' << result.pos << '\t'
                            << result.hwe_hom_ref << '/'
                            << result.hwe_hets << '/'
                            << result.hwe_hom_alt;
            hardy_weinberg_.precision(2);
            hardy_weinberg_ << std::fixed << '\t' << expected_hom_ref
                            << '/' << expected_hets << '/'
                            << expected_hom_alt;
            hardy_weinberg_.precision(6);
            hardy_weinberg_ << std::scientific << '\t' << chi_square
                            << '\t' << probabilities.two_sided
                            << '\t'
                            << probabilities.heterozygote_deficit
                            << '\t'
                            << probabilities.heterozygote_excess
                            << '\n';
        }
        for (std::size_t sample = 0; sample < sample_names_.size();
             ++sample) {
            if (options_.output_individual_depth &&
                individual_depth[sample] >= 0) {
                individual_depth_sum_[sample] +=
                    individual_depth[sample];
                ++individual_depth_count_[sample];
            }
            if (options_.output_individual_missingness) {
                const uint8_t state = individual_state[sample];
                if ((state & kGenotypeFiltered) != 0) {
                    ++individual_filtered_count_[sample];
                } else {
                    ++individual_data_count_[sample];
                    if ((state & kFirstAlleleMissing) != 0) {
                        ++individual_missing_count_[sample];
                    }
                }
            }
            if (options_.output_heterozygosity &&
                result.heterozygosity_eligible &&
                (individual_state[sample] &
                 kCompleteDiploidGenotype) != 0) {
                ++heterozygosity_sites_[sample];
                if ((individual_state[sample] &
                     kObservedHomozygote) != 0) {
                    ++observed_homozygotes_[sample];
                }
                expected_homozygotes_[sample] +=
                    result.expected_homozygosity;
            }
        }
        if (options_.output_recode) {
            recode_stream_->write(
                result.recode_line.data(),
                static_cast<std::streamsize>(result.recode_line.size()));
        }
        if (options_.output_recode_bcf) {
            recode_bcf_->write(output_header_, result.output_record);
        }
        if (options_.output_recode_vcf_gz) {
            recode_vcf_gz_->append(
                result.recode_line.data(),
                result.recode_line.size());
        }
    }

    void finish() {
        if (finished_) {
            return;
        }
        finished_ = true;
        for (std::size_t sample = 0; sample < sample_names_.size();
             ++sample) {
            if (options_.output_individual_depth) {
                const double mean =
                    individual_depth_sum_[sample] /
                    individual_depth_count_[sample];
                individual_depth_
                    << sample_names_[sample] << '\t'
                    << individual_depth_count_[sample] << '\t'
                    << mean << '\n';
            }
            if (options_.output_individual_missingness) {
                individual_missingness_
                    << sample_names_[sample] << '\t'
                    << individual_data_count_[sample] << '\t'
                    << individual_filtered_count_[sample] << '\t'
                    << individual_missing_count_[sample] << '\t'
                    << individual_missing_count_[sample] /
                           static_cast<double>(
                               individual_data_count_[sample])
                    << '\n';
            }
            if (options_.output_heterozygosity &&
                heterozygosity_sites_[sample] > 0) {
                const double inbreeding =
                    (observed_homozygotes_[sample] -
                     expected_homozygotes_[sample]) /
                    (heterozygosity_sites_[sample] -
                     expected_homozygotes_[sample]);
                heterozygosity_ << std::fixed
                                << sample_names_[sample] << '\t'
                                << observed_homozygotes_[sample] << '\t';
                heterozygosity_.precision(1);
                heterozygosity_ << expected_homozygotes_[sample] << '\t';
                heterozygosity_.precision(5);
                heterozygosity_ << heterozygosity_sites_[sample] << '\t'
                                << inbreeding << '\n';
            }
        }
        finish_window_pi();
        finish_tajima();
        finish_window_fst();
        finish_genotype_ld();
        finish_pca();
        if (recode_bcf_) {
            recode_bcf_->finish();
        }
        if (recode_vcf_gz_) {
            recode_vcf_gz_->finish();
        }
        if (recode_stream_ != nullptr) {
            recode_stream_->flush();
            if (!*recode_stream_) {
                fail("Could not finish recoded VCF output");
            }
        }
    }

private:
    using PiWindow = std::array<uint64_t, 4>;

    struct LdSite {
        std::string chrom;
        int64_t pos = 0;
        std::vector<int8_t> dosage;
        double alt_frequency = 0.0;
        bool pca_eligible = false;
    };

    void consume_window_pi(
        const SiteResult& result, uint64_t mismatches) {
        const int first_unclamped = static_cast<int>(
            std::ceil(
                (result.pos - options_.pi_window_size) /
                static_cast<double>(pi_window_step_)));
        const int first = std::max(0, first_unclamped);
        const int last = static_cast<int>(
            std::ceil(
                result.pos / static_cast<double>(pi_window_step_)));
        if (result.chrom != previous_pi_chromosome_) {
            pi_chromosomes_.push_back(result.chrom);
            previous_pi_chromosome_ = result.chrom;
            pi_windows_[result.chrom].resize(1, PiWindow{0, 0, 0, 0});
        }
        auto& windows = pi_windows_[result.chrom];
        if (last >= static_cast<int>(windows.size())) {
            windows.resize(
                static_cast<std::size_t>(last + 1),
                PiWindow{0, 0, 0, 0});
        }
        const uint64_t comparisons =
            static_cast<uint64_t>(
                result.non_missing_chromosomes) *
            (result.non_missing_chromosomes - 1);
        for (int index = first; index < last; ++index) {
            auto& window = windows[static_cast<std::size_t>(index)];
            ++window[0];
            window[1] += comparisons;
            window[2] += mismatches;
            if (result.allele_counts[0] <
                result.non_missing_chromosomes) {
                ++window[3];
            }
        }
    }

    void finish_window_pi() {
        if (options_.pi_window_size <= 0) {
            return;
        }
        const uint64_t monomorphic_comparisons =
            selected_chromosome_count_ *
            (selected_chromosome_count_ - 1);
        for (const auto& chromosome : pi_chromosomes_) {
            const auto& windows = pi_windows_[chromosome];
            for (std::size_t index = 0; index < windows.size(); ++index) {
                const auto& window = windows[index];
                if (window[3] == 0 && window[2] == 0) {
                    continue;
                }
                const uint64_t monomorphic_sites =
                    static_cast<uint64_t>(options_.pi_window_size) -
                    window[0];
                const uint64_t pairs =
                    window[1] +
                    monomorphic_sites * monomorphic_comparisons;
                const double pi = window[2] / static_cast<double>(pairs);
                window_pi_
                    << chromosome << '\t'
                    << index * pi_window_step_ + 1 << '\t'
                    << index * pi_window_step_ +
                           options_.pi_window_size
                    << '\t' << window[3] << '\t'
                    << monomorphic_sites << '\t' << pi << '\n';
            }
        }
    }

    void initialise_tajima_constants() {
        const uint64_t n = selected_chromosome_count_;
        for (uint64_t index = 1; index < n; ++index) {
            tajima_a1_ += 1.0 / static_cast<double>(index);
            tajima_a2_ +=
                1.0 / static_cast<double>(index * index);
        }
        const double b1 =
            static_cast<double>(n + 1) / 3.0 /
            static_cast<double>(n - 1);
        const double b2 =
            2.0 * static_cast<double>(n * n + n + 3) / 9.0 /
            static_cast<double>(n) / static_cast<double>(n - 1);
        const double c1 = b1 - (1.0 / tajima_a1_);
        const double c2 =
            b2 -
            (static_cast<double>(n + 2) /
             static_cast<double>(tajima_a1_ * n)) +
            (tajima_a2_ / tajima_a1_ / tajima_a1_);
        tajima_e1_ = c1 / tajima_a1_;
        tajima_e2_ =
            c2 / ((tajima_a1_ * tajima_a1_) + tajima_a2_);
    }

    void consume_tajima(const SiteResult& result) {
        const double reciprocal_window =
            1.0 / static_cast<double>(options_.tajima_window_size);
        const std::size_t index = static_cast<std::size_t>(
            result.pos * reciprocal_window);
        auto& windows = tajima_windows_[result.chrom];
        if (index >= windows.size()) {
            windows.resize(index + 1, std::pair<int, double>{0, 0.0});
        }
        if (result.chrom != previous_tajima_chromosome_) {
            tajima_chromosomes_.push_back(result.chrom);
            previous_tajima_chromosome_ = result.chrom;
        }
        const double frequency =
            result.allele_counts[0] /
            static_cast<double>(result.non_missing_chromosomes);
        if (frequency > 0.0 && frequency < 1.0) {
            ++windows[index].first;
            windows[index].second +=
                frequency * (1.0 - frequency);
        }
    }

    void finish_tajima() {
        if (options_.tajima_window_size <= 0) {
            return;
        }
        const double n =
            static_cast<double>(selected_chromosome_count_);
        for (const auto& chromosome : tajima_chromosomes_) {
            bool output = false;
            const auto& windows = tajima_windows_[chromosome];
            for (std::size_t index = 0; index < windows.size(); ++index) {
                const int segregating_sites = windows[index].first;
                double tajima_d =
                    std::numeric_limits<double>::quiet_NaN();
                if (segregating_sites > 0) {
                    const double pi =
                        2.0 * windows[index].second * n / (n - 1.0);
                    const double theta =
                        segregating_sites / tajima_a1_;
                    const double variance =
                        tajima_e1_ * segregating_sites +
                        tajima_e2_ * segregating_sites *
                            (segregating_sites - 1);
                    tajima_d = (pi - theta) / std::sqrt(variance);
                    output = true;
                }
                if (output) {
                    tajima_
                        << chromosome << '\t'
                        << index * options_.tajima_window_size << '\t'
                        << segregating_sites << '\t'
                        << tajima_d << '\n';
                }
            }
        }
    }

    void consume_window_fst(const SiteResult& result) {
        if (result.chrom != previous_fst_chromosome_) {
            fst_chromosomes_.push_back(result.chrom);
            previous_fst_chromosome_ = result.chrom;
        }
        const int first_unclamped = static_cast<int>(
            std::ceil(
                (result.pos - options_.fst_window_size) /
                static_cast<double>(fst_window_step_)));
        const int first = std::max(0, first_unclamped);
        const int last = static_cast<int>(
            std::ceil(
                result.pos / static_cast<double>(fst_window_step_)));
        auto& windows = fst_windows_[result.chrom];
        for (int index = first; index < last; ++index) {
            if (index >= static_cast<int>(windows.size())) {
                windows.resize(
                    static_cast<std::size_t>(index + 1),
                    std::array<double, 4>{0.0, 0.0, 0.0, 0.0});
            }
            auto& window = windows[static_cast<std::size_t>(index)];
            window[0] += result.fst_sum_a;
            window[1] += result.fst_sum_all;
            window[2] += result.fst;
            ++window[3];
        }
    }

    void finish_window_fst() {
        if (options_.fst_window_size <= 0) {
            return;
        }
        for (const auto& chromosome : fst_chromosomes_) {
            const auto& windows = fst_windows_[chromosome];
            for (std::size_t index = 0; index < windows.size(); ++index) {
                const auto& window = windows[index];
                if (window[1] != 0.0 &&
                    !std::isnan(window[0]) &&
                    !std::isnan(window[1]) && window[3] > 0.0) {
                    const double weighted_fst =
                        window[0] / window[1];
                    const double mean_fst =
                        window[2] / window[3];
                    window_fst_
                        << chromosome << '\t'
                        << index * fst_window_step_ + 1 << '\t'
                        << index * fst_window_step_ +
                               options_.fst_window_size
                        << '\t' << window[3] << '\t'
                        << weighted_fst << '\t' << mean_fst << '\n';
                }
            }
        }
    }

    static std::pair<double, int> genotype_r2(
        const LdSite& first, const LdSite& second) {
        double x = 0.0;
        double x2 = 0.0;
        double y = 0.0;
        double y2 = 0.0;
        double xy = 0.0;
        int individual_count = 0;
        for (std::size_t sample = 0;
             sample < first.dosage.size(); ++sample) {
            if (first.dosage[sample] < 0 ||
                second.dosage[sample] < 0) {
                continue;
            }
            double sx = first.dosage[sample];
            double sy = second.dosage[sample];
            x += sx;
            y += sy;
            xy += sx * sy;
            sx *= sx;
            sy *= sy;
            x2 += sx;
            y2 += sy;
            ++individual_count;
        }
        x /= individual_count;
        x2 /= individual_count;
        y /= individual_count;
        y2 /= individual_count;
        xy /= individual_count;
        const double variance_x = x2 - x * x;
        const double variance_y = y2 - y * y;
        const double covariance = xy - x * y;
        return {
            covariance * covariance /
                (variance_x * variance_y),
            individual_count};
    }

    void finish_genotype_ld() {
        if (!options_.output_genotype_r2 || ld_sites_.size() < 2) {
            return;
        }
        std::vector<std::string> ordered_lines(ld_sites_.size() - 1);
        std::atomic<std::size_t> next_left{0};
        const std::size_t worker_count = std::min<std::size_t>(
            options_.threads, ld_sites_.size() - 1);
        std::vector<std::thread> workers;
        workers.reserve(worker_count);
        const int minimum_snp_distance =
            std::max(1, options_.ld_snp_window_min);
        for (std::size_t worker = 0; worker < worker_count; ++worker) {
            workers.emplace_back([&] {
                while (true) {
                    const std::size_t left =
                        next_left.fetch_add(1);
                    if (left + 1 >= ld_sites_.size()) {
                        return;
                    }
                    const auto& first = ld_sites_[left];
                    std::ostringstream lines;
                    for (std::size_t right = left + 1;
                         right < ld_sites_.size(); ++right) {
                        if (static_cast<int>(right - left) >
                            options_.ld_snp_window_size) {
                            break;
                        }
                        if (right <
                            left +
                                static_cast<std::size_t>(
                                    minimum_snp_distance)) {
                            continue;
                        }
                        const auto& second = ld_sites_[right];
                        if (first.chrom != second.chrom) {
                            continue;
                        }
                        const int64_t distance =
                            second.pos - first.pos;
                        if (distance < options_.ld_bp_window_min) {
                            continue;
                        }
                        if (distance > options_.ld_bp_window_size) {
                            break;
                        }
                        const auto [r2, individual_count] =
                            genotype_r2(first, second);
                        if (options_.min_r2 > 0.0 &&
                            (r2 < options_.min_r2 ||
                             std::isnan(r2))) {
                            continue;
                        }
                        lines << first.chrom << '\t' << first.pos
                              << '\t' << second.pos << '\t'
                              << individual_count << '\t' << r2
                              << '\n';
                    }
                    ordered_lines[left] = std::move(lines).str();
                }
            });
        }
        for (auto& worker : workers) {
            worker.join();
        }
        for (const auto& lines : ordered_lines) {
            genotype_ld_ << lines;
        }
    }

    void finish_pca() {
        if (!options_.output_pca) {
            return;
        }
        const int individual_count =
            static_cast<int>(sample_names_.size());
        const std::size_t site_count = static_cast<std::size_t>(
            std::count_if(
                ld_sites_.begin(), ld_sites_.end(),
                [](const LdSite& site) {
                    return site.pca_eligible;
                }));
        if (static_cast<std::size_t>(individual_count) >= site_count) {
            fail(
                "PCA computation requires that there are more sites "
                "than individuals");
        }

        std::vector<double> matrix(
            static_cast<std::size_t>(individual_count) *
                individual_count,
            0.0);
        std::atomic<std::size_t> next_cell{0};
        const std::size_t cell_count =
            static_cast<std::size_t>(individual_count) *
            individual_count;
        const std::size_t worker_count = std::min<std::size_t>(
            options_.threads, cell_count);
        std::vector<std::thread> workers;
        workers.reserve(worker_count);
        for (std::size_t worker = 0; worker < worker_count; ++worker) {
            workers.emplace_back([&] {
                while (true) {
                    const std::size_t cell =
                        next_cell.fetch_add(1);
                    if (cell >= cell_count) {
                        return;
                    }
                    const int row = static_cast<int>(
                        cell / individual_count);
                    const int column = static_cast<int>(
                        cell % individual_count);
                    if (column < row) {
                        continue;
                    }
                    double sum = 0.0;
                    for (const auto& site : ld_sites_) {
                        if (!site.pca_eligible) {
                            continue;
                        }
                        const double mean =
                            site.alt_frequency * 2.0;
                        const double divisor =
                            1.0 /
                            std::sqrt(
                                site.alt_frequency *
                                (1.0 - site.alt_frequency));
                        const double row_genotype =
                            2.0 - site.dosage[row];
                        const double column_genotype =
                            2.0 - site.dosage[column];
                        const double row_value =
                            options_.pca_normalise
                                ? (row_genotype - mean) * divisor
                                : row_genotype - mean;
                        const double column_value =
                            options_.pca_normalise
                                ? (column_genotype - mean) * divisor
                                : column_genotype - mean;
                        sum += row_value * column_value;
                    }
                    matrix[cell] = sum;
                }
            });
        }
        for (auto& worker : workers) {
            worker.join();
        }
        for (int row = 0; row < individual_count; ++row) {
            for (int column = 0; column < row; ++column) {
                matrix[
                    static_cast<std::size_t>(row) *
                        individual_count +
                    column] =
                    matrix[
                        static_cast<std::size_t>(column) *
                            individual_count +
                        row];
            }
        }
        for (double& value : matrix) {
            value /= site_count;
        }

        const auto eigen = legacy_dgeev(matrix, individual_count);
        for (const double imaginary : eigen.imaginary) {
            if (imaginary != 0.0) {
                fail("Complex eigenvalue");
            }
        }
        pca_ << "INDV";
        for (int index = 0; index < individual_count; ++index) {
            pca_ << "\tEIG_" << index;
        }
        pca_ << '\n' << "EIGENVALUE";
        for (const double value : eigen.real) {
            pca_ << '\t' << value;
        }
        pca_ << '\n';
        for (int row = 0; row < individual_count; ++row) {
            pca_ << sample_names_[row];
            for (int column = 0; column < individual_count; ++column) {
                pca_ << '\t'
                     << eigen.vectors[
                            row + column * individual_count];
            }
            pca_ << '\n';
        }
    }

    std::ofstream open(const std::string& suffix) {
        const std::string path = options_.output_prefix + suffix;
        std::ofstream stream(path);
        if (!stream) {
            fail("Could not open output file: " + path);
        }
        return stream;
    }

    const Options& options_;
    bcf_hdr_t* output_header_;
    std::ofstream freq_;
    std::ofstream counts_;
    std::ofstream missing_;
    std::ofstream depth_;
    std::ofstream mean_depth_;
    std::ofstream individual_depth_;
    std::ofstream individual_missingness_;
    std::ofstream heterozygosity_;
    std::ofstream hardy_weinberg_;
    std::ofstream site_quality_;
    std::ofstream site_pi_;
    std::ofstream window_pi_;
    std::ofstream tajima_;
    std::ofstream site_fst_;
    std::ofstream window_fst_;
    std::ofstream genotype_ld_;
    std::ofstream pca_;
    std::ofstream recode_;
    std::ostream* recode_stream_ = nullptr;
    std::unique_ptr<ExactBcfOutput> recode_bcf_;
    std::unique_ptr<DeterministicBgzfWriter> recode_vcf_gz_;
    std::vector<std::string> sample_names_;
    std::vector<double> individual_depth_sum_;
    std::vector<uint64_t> individual_depth_count_;
    std::vector<uint64_t> individual_data_count_;
    std::vector<uint64_t> individual_filtered_count_;
    std::vector<uint64_t> individual_missing_count_;
    std::vector<uint64_t> heterozygosity_sites_;
    std::vector<uint64_t> observed_homozygotes_;
    std::vector<double> expected_homozygotes_;
    uint64_t selected_chromosome_count_ = 0;
    int pi_window_step_ = 0;
    std::map<std::string, std::vector<PiWindow>> pi_windows_;
    std::vector<std::string> pi_chromosomes_;
    std::string previous_pi_chromosome_;
    std::map<
        std::string, std::vector<std::pair<int, double>>>
        tajima_windows_;
    std::vector<std::string> tajima_chromosomes_;
    std::string previous_tajima_chromosome_;
    double tajima_a1_ = 0.0;
    double tajima_a2_ = 0.0;
    double tajima_e1_ = 0.0;
    double tajima_e2_ = 0.0;
    int fst_window_step_ = 0;
    std::map<
        std::string, std::vector<std::array<double, 4>>>
        fst_windows_;
    std::vector<std::string> fst_chromosomes_;
    std::string previous_fst_chromosome_;
    std::vector<LdSite> ld_sites_;
    bool finished_ = false;
};

struct PipelineSummary {
    uint64_t total = 0;
    uint64_t kept = 0;
};

class OrderedCommitter {
public:
    OrderedCommitter(
        const Options& options, bcf_hdr_t* output_header,
        bool sample_selection_active)
        : outputs_(
              options, output_header, sample_selection_active),
          thin_selector_(options.thin_distance) {}

    void commit(std::vector<SiteResult>& results,
                BatchAnalysisPayload& analysis) {
        for (std::size_t row = 0; row < results.size(); ++row) {
            auto& result = results[row];
            if (result.kept && !thin_selector_.keep(result)) {
                result.kept = false;
            }
            outputs_.write(
                result, analysis.state_row(row),
                analysis.depth_row(row),
                analysis.dosage_row(row));
            ++summary_.total;
            summary_.kept += result.kept;
        }
    }

    void finish() {
        outputs_.finish();
    }

    const PipelineSummary& summary() const {
        return summary_;
    }

private:
    Outputs outputs_;
    ThinSelector thin_selector_;
    PipelineSummary summary_;
};

struct PipelineBatch {
    std::size_t id = 0;
    std::vector<vcftools_ng::input::RecordPtr> records;
    std::vector<SiteResult> results;
    BatchAnalysisPayload analysis;
    std::atomic<std::size_t> remaining_slices{0};
};

struct PipelineSlice {
    std::shared_ptr<PipelineBatch> batch;
    std::size_t begin = 0;
    std::size_t end = 0;
};

struct PipelineState {
    std::mutex mutex;
    std::condition_variable work_available;
    std::condition_variable completed_available;
    std::condition_variable batch_slot_available;
    std::deque<PipelineSlice> ready_slices;
    std::map<std::size_t, std::shared_ptr<PipelineBatch>>
        completed_batches;
    std::size_t in_flight_batches = 0;
    std::size_t produced_batches = 0;
    bool reader_done = false;
    bool cancelled = false;
    std::exception_ptr error;
};

PipelineSummary run_ordered_pipeline(
    const Options& options, const SiteSelector& selector,
    const SampleSelection& samples, bcf_hdr_t* header,
    vcftools_ng::input::OrderedShardSource& source,
    unsigned compute_threads, OrderedCommitter& committer) {
    constexpr std::size_t max_in_flight_batches = 3;
    const std::size_t slices_per_worker = 4;
    const std::size_t target_slice =
        options.batch_size /
        std::max<std::size_t>(
            1, static_cast<std::size_t>(compute_threads) *
                   slices_per_worker);
    const std::size_t slice_size =
        std::clamp<std::size_t>(target_slice, 64, 256);
    PipelineState state;
    std::vector<HeaderPtr> worker_output_headers;
    worker_output_headers.reserve(compute_threads);
    for (unsigned worker = 0; worker < compute_threads; ++worker) {
        HeaderPtr duplicate(
            bcf_hdr_dup(samples.output_header(header)));
        if (!duplicate) {
            fail("Could not duplicate worker output header");
        }
        worker_output_headers.push_back(std::move(duplicate));
    }

    const auto record_failure = [&](std::exception_ptr error) {
        {
            std::lock_guard lock(state.mutex);
            if (!state.error) {
                state.error = error;
            }
            state.cancelled = true;
        }
        state.work_available.notify_all();
        state.completed_available.notify_all();
        state.batch_slot_available.notify_all();
    };

    std::thread reader([&] {
        try {
            bool eof = false;
            while (!eof) {
                {
                    std::unique_lock lock(state.mutex);
                    state.batch_slot_available.wait(lock, [&] {
                        return state.cancelled ||
                               state.in_flight_batches <
                                   max_in_flight_batches;
                    });
                    if (state.cancelled) {
                        return;
                    }
                }

                auto batch = std::make_shared<PipelineBatch>();
                batch->records =
                    source.next_batch(options.batch_size);
                eof = batch->records.empty();

                if (batch->records.empty()) {
                    std::lock_guard lock(state.mutex);
                    state.reader_done = true;
                    state.work_available.notify_all();
                    state.completed_available.notify_all();
                    return;
                }

                batch->results.resize(batch->records.size());
                batch->analysis.reset(
                    batch->records.size(),
                    static_cast<std::size_t>(samples.count()), options);
                const std::size_t slice_count =
                    (batch->records.size() + slice_size - 1) /
                    slice_size;
                batch->remaining_slices.store(
                    slice_count, std::memory_order_relaxed);
                {
                    std::lock_guard lock(state.mutex);
                    batch->id = state.produced_batches++;
                    ++state.in_flight_batches;
                    for (std::size_t begin = 0;
                         begin < batch->records.size();
                         begin += slice_size) {
                        state.ready_slices.push_back(PipelineSlice{
                            batch, begin,
                            std::min(
                                begin + slice_size,
                                batch->records.size())});
                    }
                    if (eof) {
                        state.reader_done = true;
                    }
                }
                state.work_available.notify_all();
                if (eof) {
                    state.completed_available.notify_all();
                    return;
                }
            }
        } catch (...) {
            record_failure(std::current_exception());
        }
    });

    std::vector<std::thread> workers;
    workers.reserve(compute_threads);
    for (unsigned worker_index = 0; worker_index < compute_threads;
         ++worker_index) {
        workers.emplace_back([&, worker_index] {
            Scratch scratch;
            bcf_hdr_t* worker_output_header =
                worker_output_headers[worker_index].get();
            try {
                while (true) {
                    PipelineSlice slice;
                    {
                        std::unique_lock lock(state.mutex);
                        state.work_available.wait(lock, [&] {
                            return state.cancelled ||
                                   !state.ready_slices.empty() ||
                                   state.reader_done;
                        });
                        if (state.cancelled) {
                            return;
                        }
                        if (state.ready_slices.empty()) {
                            if (state.reader_done) {
                                return;
                            }
                            continue;
                        }
                        slice = std::move(state.ready_slices.front());
                        state.ready_slices.pop_front();
                    }

                    for (std::size_t index = slice.begin;
                         index < slice.end; ++index) {
                        slice.batch->results[index] = process_site(
                            options, selector, samples, header,
                            worker_output_header,
                            slice.batch->records[index].get(), scratch,
                            slice.batch->analysis.state_row(index),
                            slice.batch->analysis.depth_row(index),
                            slice.batch->analysis.dosage_row(index));
                    }
                    if (slice.batch->remaining_slices.fetch_sub(
                            1, std::memory_order_acq_rel) == 1) {
                        {
                            std::lock_guard lock(state.mutex);
                            state.completed_batches.emplace(
                                slice.batch->id, slice.batch);
                        }
                        state.completed_available.notify_one();
                    }
                }
            } catch (...) {
                record_failure(std::current_exception());
            }
        });
    }

    std::size_t next_batch = 0;
    try {
        while (true) {
            std::shared_ptr<PipelineBatch> batch;
            {
                std::unique_lock lock(state.mutex);
                state.completed_available.wait(lock, [&] {
                    return state.cancelled ||
                           state.completed_batches.contains(next_batch) ||
                           (state.reader_done &&
                            next_batch == state.produced_batches);
                });
                if (state.cancelled) {
                    break;
                }
                if (state.reader_done &&
                    next_batch == state.produced_batches) {
                    break;
                }
                auto completed =
                    state.completed_batches.find(next_batch);
                batch = std::move(completed->second);
                state.completed_batches.erase(completed);
            }

            committer.commit(batch->results, batch->analysis);
            ++next_batch;
            {
                std::lock_guard lock(state.mutex);
                --state.in_flight_batches;
            }
            state.batch_slot_available.notify_one();
        }
    } catch (...) {
        record_failure(std::current_exception());
    }

    reader.join();
    for (auto& worker : workers) {
        worker.join();
    }
    if (state.error) {
        std::rethrow_exception(state.error);
    }
    return committer.summary();
}

bool has_unsupported_diff_filters(const Options& options) {
    return
        !options.site_filters_to_keep.empty() ||
        !options.site_filters_to_remove.empty() ||
        options.remove_all_filtered_sites ||
        !options.info_flags_to_keep.empty() ||
        !options.info_flags_to_remove.empty() ||
        !options.genotype_filters_to_remove.empty() ||
        options.remove_all_filtered_genotypes ||
        options.min_alleles != -1 ||
        options.max_alleles != std::numeric_limits<int>::max() ||
        options.remove_indels || options.keep_only_indels ||
        options.min_qual >= 0.0 || options.min_gq >= 0.0 ||
        options.min_dp >= 0 ||
        options.max_dp != std::numeric_limits<int>::max() ||
        options.min_mean_dp >= 0.0 ||
        options.max_mean_dp != std::numeric_limits<double>::max() ||
        options.min_call_rate > 0.0 ||
        options.max_missing_count != std::numeric_limits<int>::max() ||
        options.min_maf >= 0.0 ||
        options.max_maf != std::numeric_limits<double>::max() ||
        options.min_mac >= 0 ||
        options.max_mac != std::numeric_limits<int>::max() ||
        options.min_hwe >= 0.0 || options.thin_distance >= 0 ||
        options.min_non_ref_af >= 0.0 ||
        options.max_non_ref_af != std::numeric_limits<double>::max() ||
        options.min_non_ref_af_any >= 0.0 ||
        options.max_non_ref_af_any !=
            std::numeric_limits<double>::max() ||
        options.min_non_ref_ac >= 0 ||
        options.max_non_ref_ac != std::numeric_limits<int>::max() ||
        options.min_non_ref_ac_any >= 0 ||
        options.max_non_ref_ac_any != std::numeric_limits<int>::max();
}

class DiffRecordStream {
public:
    DiffRecordStream(
        const std::string& path, const Options& options,
        unsigned io_threads)
        : path_(path),
          input_(hts_open(path.c_str(), "r")) {
        if (!input_) {
            fail("Could not open diff input: " + path);
        }
        if (io_threads > 0 &&
            hts_set_threads(
                input_.get(), static_cast<int>(io_threads)) != 0) {
            fail("Could not enable diff input threads: " + path);
        }
        header_.reset(bcf_hdr_read(input_.get()));
        if (!header_) {
            fail("Could not read diff input header: " + path);
        }
        selector_ =
            std::make_unique<SiteSelector>(options, header_.get());
        record_.reset(bcf_init());
        advance();
    }

    bool present() const {
        return present_;
    }

    bcf1_t* record() const {
        return record_.get();
    }

    bcf_hdr_t* header() const {
        return header_.get();
    }

    std::string chromosome() const {
        return bcf_hdr_id2name(header_.get(), record_->rid);
    }

    int position() const {
        return record_->pos + 1;
    }

    void advance() {
        present_ = false;
        while (bcf_read(input_.get(), header_.get(), record_.get()) == 0) {
            if (!selector_->keep(record_.get(), header_.get())) {
                continue;
            }
            bcf_unpack(record_.get(), BCF_UN_ALL);
            present_ = true;
            return;
        }
    }

private:
    std::string path_;
    HtsFilePtr input_;
    HeaderPtr header_;
    std::unique_ptr<SiteSelector> selector_;
    RecordPtr record_;
    bool present_ = false;
};

struct DiffGenotypeScratch {
    int32_t* first = nullptr;
    int first_capacity = 0;
    int32_t* second = nullptr;
    int second_capacity = 0;

    ~DiffGenotypeScratch() {
        std::free(first);
        std::free(second);
    }
};

using DiffSampleMap =
    std::map<std::string, std::pair<int, int>>;

DiffSampleMap make_diff_sample_map(
    const SampleSelection& first_selection, bcf_hdr_t* first_header,
    const SampleSelection& second_selection, bcf_hdr_t* second_header) {
    DiffSampleMap samples;
    for (const int sample : first_selection.indices()) {
        samples[first_header->samples[sample]].first = sample;
        samples[first_header->samples[sample]].second = -1;
    }
    for (const int sample : second_selection.indices()) {
        const std::string name = second_header->samples[sample];
        const auto found = samples.find(name);
        if (found == samples.end()) {
            samples[name] = {-1, sample};
        } else {
            found->second.second = sample;
        }
    }
    return samples;
}

std::string reference_allele(bcf1_t* record) {
    return record->n_allele > 0 ? record->d.allele[0] : "";
}

std::string alternate_alleles(bcf1_t* record) {
    std::ostringstream output;
    for (int allele = 1; allele < record->n_allele; ++allele) {
        if (allele > 1) {
            output << ',';
        }
        output << record->d.allele[allele];
    }
    return output.str();
}

struct DiploidGenotype {
    int first = -1;
    int second = -1;

    bool fully_missing() const {
        return first < 0 && second < 0;
    }
};

DiploidGenotype get_diploid_genotype(
    const int32_t* genotypes, int count, int sample_count, int sample) {
    if (count <= 0 || sample_count <= 0) {
        return {};
    }
    const int ploidy = count / sample_count;
    if (ploidy < 2) {
        return {};
    }
    const int32_t first = genotypes[sample * ploidy];
    const int32_t second = genotypes[sample * ploidy + 1];
    const auto decode = [](int32_t value) {
        return value == bcf_int32_vector_end ||
                       bcf_gt_is_missing(value)
                   ? -1
                   : bcf_gt_allele(value);
    };
    return {decode(first), decode(second)};
}

bool same_unphased_genotype(
    const DiploidGenotype& first, const DiploidGenotype& second) {
    return
        (first.first == second.first &&
         first.second == second.second) ||
        (first.first == second.second &&
         first.second == second.first);
}

bool same_unphased_allele_strings(
    const DiploidGenotype& first, bcf1_t* first_record,
    const DiploidGenotype& second, bcf1_t* second_record) {
    const auto allele_string = [](bcf1_t* record, int allele) {
        return allele < 0
            ? std::string(".")
            : std::string(record->d.allele[allele]);
    };
    const std::string first_a =
        allele_string(first_record, first.first);
    const std::string first_b =
        allele_string(first_record, first.second);
    const std::string second_a =
        allele_string(second_record, second.first);
    const std::string second_b =
        allele_string(second_record, second.second);
    return
        (first_a == second_a && first_b == second_b) ||
        (first_a == second_b && first_b == second_a);
}

struct DiffSiteCounts {
    uint64_t common_called = 0;
    uint64_t discordant = 0;
};

DiffSiteCounts compare_site_genotypes(
    bcf_hdr_t* first_header, bcf1_t* first_record,
    bcf_hdr_t* second_header, bcf1_t* second_record,
    const DiffSampleMap& samples, DiffGenotypeScratch& scratch,
    std::map<std::string, DiffSiteCounts>* individual_counts) {
    const int first_count = bcf_get_genotypes(
        first_header, first_record, &scratch.first,
        &scratch.first_capacity);
    const int second_count = bcf_get_genotypes(
        second_header, second_record, &scratch.second,
        &scratch.second_capacity);
    const bool matching_alleles =
        reference_allele(first_record) ==
            reference_allele(second_record) &&
        alternate_alleles(first_record) ==
            alternate_alleles(second_record);

    DiffSiteCounts site;
    for (const auto& [name, indices] : samples) {
        if (indices.first < 0 || indices.second < 0) {
            continue;
        }
        const auto first = get_diploid_genotype(
            scratch.first, first_count,
            bcf_hdr_nsamples(first_header), indices.first);
        const auto second = get_diploid_genotype(
            scratch.second, second_count,
            bcf_hdr_nsamples(second_header), indices.second);
        if (first.fully_missing() || second.fully_missing()) {
            continue;
        }
        ++site.common_called;
        bool match = false;
        if (matching_alleles) {
            match = same_unphased_genotype(first, second);
        } else {
            match = same_unphased_allele_strings(
                first, first_record, second, second_record);
        }
        if (!match) {
            ++site.discordant;
        }
        if (individual_counts != nullptr) {
            auto& individual = (*individual_counts)[name];
            ++individual.common_called;
            individual.discordant += !match;
        }
    }
    return site;
}

void write_diff_individuals(
    const Options& options, const DiffSampleMap& samples) {
    std::ofstream output(
        options.output_prefix + ".diff.indv_in_files");
    if (!output) {
        fail("Could not open diff individual output");
    }
    output << "INDV\tFILES\n";
    for (const auto& [name, indices] : samples) {
        output << name << '\t';
        if (indices.first >= 0 && indices.second >= 0) {
            output << 'B';
        } else if (indices.first >= 0) {
            output << '1';
        } else {
            output << '2';
        }
        output << '\n';
    }
}

int run_diff(const Options& options) {
    if (has_unsupported_diff_filters(options)) {
        fail(
            "This v0.11 diff engine currently supports chromosome, "
            "position, BED and sample selection filters; other filters "
            "must be applied before comparison");
    }
    if (options.output_diff_sites_in_files &&
        (options.output_diff_site_discordance ||
         options.output_diff_individual_discordance)) {
        fail(
            "--diff-site cannot yet share a scan with discordance outputs "
            "when overlapping indels are present");
    }

    const unsigned io_threads =
        options.threads > 2
            ? std::min(4u, std::max(1u, options.threads / 2))
            : 0;
    DiffRecordStream first(options.input, options, io_threads);
    DiffRecordStream second(options.diff_input, options, io_threads);
    SampleSelection first_selection(options, first.header());
    SampleSelection second_selection(options, second.header());
    const auto samples = make_diff_sample_map(
        first_selection, first.header(),
        second_selection, second.header());

    if (options.output_diff_individuals_in_files) {
        write_diff_individuals(options, samples);
    }
    const bool needs_site_scan =
        options.output_diff_sites_in_files ||
        options.output_diff_site_discordance ||
        options.output_diff_individual_discordance;
    if (!needs_site_scan) {
        return 0;
    }

    std::ofstream sites_in_files;
    if (options.output_diff_sites_in_files) {
        sites_in_files.open(
            options.output_prefix + ".diff.sites_in_files");
        if (!sites_in_files) {
            fail("Could not open diff site-membership output");
        }
        sites_in_files
            << "CHROM\tPOS1\tPOS2\tIN_FILE"
            << "\tREF1\tREF2\tALT1\tALT2\n";
    }
    std::ofstream site_discordance;
    if (options.output_diff_site_discordance) {
        site_discordance.open(
            options.output_prefix + ".diff.sites");
        if (!site_discordance) {
            fail("Could not open site discordance output");
        }
        site_discordance
            << "CHROM\tPOS\tFILES\tMATCHING_ALLELES"
            << "\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n";
    }

    std::map<std::string, int> chromosome_order;
    for (int rid = 0;
         rid < first.header()->n[BCF_DT_CTG]; ++rid) {
        chromosome_order.emplace(
            bcf_hdr_id2name(first.header(), rid), rid);
    }
    std::map<std::string, DiffSiteCounts> individual_counts;
    DiffGenotypeScratch scratch;
    uint64_t common_sites = 0;
    uint64_t first_only = 0;
    uint64_t second_only = 0;

    while (first.present() || second.present()) {
        int comparison = 0;
        if (!second.present()) {
            comparison = -1;
        } else if (!first.present()) {
            comparison = 1;
        } else if (first.chromosome() == second.chromosome()) {
            comparison =
                first.position() < second.position()
                    ? -1
                    : (first.position() > second.position() ? 1 : 0);
        } else {
            const auto first_rank =
                chromosome_order.find(first.chromosome());
            const auto second_rank =
                chromosome_order.find(second.chromosome());
            if (first_rank == chromosome_order.end() ||
                second_rank == chromosome_order.end()) {
                fail(
                    "Both diff inputs must use the same chromosome "
                    "ordering");
            }
            comparison =
                first_rank->second < second_rank->second ? -1 : 1;
        }

        if (comparison < 0 && options.output_diff_sites_in_files &&
            second.present() &&
            first.chromosome() == second.chromosome() &&
            second.position() <
                first.position() +
                    static_cast<int>(
                        reference_allele(first.record()).size())) {
            sites_in_files
                << first.chromosome() << '\t' << first.position()
                << '\t' << second.position() << "\tO\t"
                << reference_allele(first.record()) << '\t'
                << reference_allele(second.record()) << '\t'
                << alternate_alleles(first.record()) << '\t'
                << alternate_alleles(second.record()) << '\n';
            first.advance();
            second.advance();
            continue;
        }
        if (comparison > 0 && options.output_diff_sites_in_files &&
            first.present() &&
            first.chromosome() == second.chromosome() &&
            first.position() <
                second.position() +
                    static_cast<int>(
                        reference_allele(second.record()).size())) {
            sites_in_files
                << first.chromosome() << '\t' << first.position()
                << '\t' << second.position() << "\tO\t"
                << reference_allele(first.record()) << '\t'
                << reference_allele(second.record()) << '\t'
                << alternate_alleles(first.record()) << '\t'
                << alternate_alleles(second.record()) << '\n';
            first.advance();
            second.advance();
            continue;
        }

        if (comparison < 0) {
            ++first_only;
            if (options.output_diff_sites_in_files) {
                sites_in_files
                    << first.chromosome() << '\t' << first.position()
                    << "\t.\t1\t"
                    << reference_allele(first.record()) << "\t.\t"
                    << alternate_alleles(first.record()) << "\t.\n";
            }
            if (options.output_diff_site_discordance) {
                site_discordance
                    << first.chromosome() << '\t' << first.position()
                    << "\t1\t0\t0\t0\t"
                    << std::numeric_limits<double>::quiet_NaN()
                    << '\n';
            }
            first.advance();
            continue;
        }
        if (comparison > 0) {
            ++second_only;
            if (options.output_diff_sites_in_files) {
                sites_in_files
                    << second.chromosome() << "\t.\t"
                    << second.position() << "\t2\t.\t"
                    << reference_allele(second.record()) << "\t.\t"
                    << alternate_alleles(second.record()) << '\n';
            }
            if (options.output_diff_site_discordance) {
                site_discordance
                    << second.chromosome() << '\t' << second.position()
                    << "\t2\t0\t0\t0\t"
                    << std::numeric_limits<double>::quiet_NaN()
                    << '\n';
            }
            second.advance();
            continue;
        }

        std::string first_ref = reference_allele(first.record());
        std::string second_ref = reference_allele(second.record());
        if (first_ref == "N" || first_ref == "." || first_ref.empty()) {
            first_ref = second_ref;
        }
        if (second_ref == "N" || second_ref == "." ||
            second_ref.empty()) {
            second_ref = first_ref;
        }
        if (first_ref != second_ref &&
            first_ref != "N" && second_ref != "N" &&
            first_ref != "." && second_ref != "." &&
            !first_ref.empty() && !second_ref.empty()) {
            if (options.output_diff_sites_in_files) {
                sites_in_files
                    << first.chromosome() << '\t' << first.position()
                    << '\t' << second.position() << "\tO\t"
                    << first_ref << '\t' << second_ref << '\t'
                    << alternate_alleles(first.record()) << '\t'
                    << alternate_alleles(second.record()) << '\n';
            }
            first.advance();
            second.advance();
            continue;
        }

        ++common_sites;
        const bool matching_alleles =
            first_ref == second_ref &&
            alternate_alleles(first.record()) ==
                alternate_alleles(second.record());
        if (options.output_diff_sites_in_files) {
            sites_in_files
                << first.chromosome() << '\t' << first.position()
                << '\t' << second.position() << "\tB\t"
                << first_ref << '\t' << second_ref << '\t'
                << alternate_alleles(first.record()) << '\t'
                << alternate_alleles(second.record()) << '\n';
        }
        if (options.output_diff_site_discordance ||
            options.output_diff_individual_discordance) {
            const auto counts = compare_site_genotypes(
                first.header(), first.record(),
                second.header(), second.record(), samples, scratch,
                options.output_diff_individual_discordance
                    ? &individual_counts
                    : nullptr);
            if (options.output_diff_site_discordance) {
                site_discordance
                    << first.chromosome() << '\t' << first.position()
                    << "\tB\t" << matching_alleles << '\t'
                    << counts.common_called << '\t'
                    << counts.discordant << '\t'
                    << counts.discordant /
                           static_cast<double>(counts.common_called)
                    << '\n';
            }
        }
        first.advance();
        second.advance();
    }

    if (options.output_diff_individual_discordance) {
        std::ofstream output(options.output_prefix + ".diff.indv");
        if (!output) {
            fail("Could not open individual discordance output");
        }
        output
            << "INDV\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n";
        for (const auto& [name, indices] : samples) {
            (void)indices;
            const auto counts = individual_counts[name];
            output << name << '\t' << counts.common_called << '\t'
                   << counts.discordant << '\t'
                   << counts.discordant /
                          static_cast<double>(counts.common_called)
                   << '\n';
        }
    }

    std::cerr << "Diff common sites: " << common_sites << "\n"
              << "Diff sites only in first: " << first_only << "\n"
              << "Diff sites only in second: " << second_only << "\n";
    return 0;
}

int run(const Options& options) {
    if (!options.diff_input.empty()) {
        return run_diff(options);
    }
    if (can_use_fused_site_stats(options)) {
        vcftools_ng::input::SourceOptions fast_options{
            .path = options.input,
            .requested_backend = options.input_backend,
            .total_threads = options.threads,
            .target_batch_records = options.batch_size,
            .parallel_safe = true,
            .workload =
                vcftools_ng::input::WorkloadProfile::
                    compact_site_statistics,
            .bcftools_path = options.bcftools_path,
            .index_path = {},
            .selected_contigs = {},
            .start_position = -1,
            .end_position = std::numeric_limits<int>::max(),
        };
        const vcftools_ng::FastSiteStatPlan fast_plan{
            .freq = options.output_freq,
            .freq2 = options.output_freq2,
            .counts = options.output_counts,
            .missing_site = options.output_missing_site,
            .site_depth = options.output_site_depth,
            .site_mean_depth = options.output_site_mean_depth,
            .site_quality = options.output_site_quality,
        };
        const auto fast = vcftools_ng::run_fast_text_site_stats(
            options.output_prefix, fast_options, fast_plan);
        if (fast.has_value()) {
            const auto detected_threads =
                vcftools_ng::input::detect_available_threads();
            std::cerr
                << "vcftools-ng " << kVersion << "\n"
                << "Input: " << options.input << "\n"
                << "Input backend: " << fast->backend
                << " (" << fast->description << ")\n"
                << "Samples: " << fast->samples << "\n"
                << "Threads: " << options.threads
                << (options.threads_explicit &&
                            options.requested_threads != options.threads
                        ? " (capped from user request " +
                              std::to_string(
                                  options.requested_threads) +
                              " by " + detected_threads.source + ")"
                        : options.threads_explicit
                        ? " (user specified)"
                        : " (auto from " +
                              detected_threads.source + ")")
                << "\n"
                << "Stage concurrency: input "
                << fast->input_threads
                << ", HTSlib I/O " << fast->hts_io_threads
                << ", compute fused\n"
                << "Planned input shards: "
                << fast->planned_shards << "\n"
                << "Selected samples: " << fast->samples << "\n"
                << "Scheduler: adaptive fused text site statistics\n"
                << "After filtering, kept " << fast->kept
                << " out of " << fast->total << " sites\n";
            return 0;
        }
    }
    vcftools_ng::input::SourceOptions source_options{
        .path = options.input,
        .requested_backend = options.input_backend,
        .total_threads = options.threads,
        .target_batch_records = options.batch_size,
        .parallel_safe = true,
        .workload =
            options.output_recode ||
                    options.output_recode_bcf ||
                    options.output_recode_vcf_gz
                ? vcftools_ng::input::WorkloadProfile::full_recode
                : vcftools_ng::input::WorkloadProfile::general,
        .bcftools_path = options.bcftools_path,
        .index_path = {},
        .selected_contigs = options.chromosomes_to_keep,
        .start_position = options.start_position,
        .end_position = options.end_position,
    };
    auto source =
        vcftools_ng::input::make_ordered_source(source_options);
    bcf_hdr_t* header = source->header();
    validate_info_flag_filters(options, header);
    const auto detected_threads =
        vcftools_ng::input::detect_available_threads();
    const auto& resources = source->resources();

    std::cerr << "vcftools-ng " << kVersion << "\n"
              << "Input: " << options.input << "\n"
              << "Input backend: " << source->backend_name()
              << " (" << source->description() << ")\n"
              << "Samples: " << bcf_hdr_nsamples(header) << "\n"
              << "Threads: " << options.threads
              << (options.threads_explicit &&
                          options.requested_threads != options.threads
                      ? " (capped from user request " +
                            std::to_string(options.requested_threads) +
                            " by " + detected_threads.source + ")"
                      : options.threads_explicit
                      ? " (user specified)"
                      : " (auto from " +
                            detected_threads.source + ")")
              << "\n"
              << "Stage concurrency: input "
              << resources.input_threads
              << ", HTSlib I/O " << resources.hts_io_threads
              << ", compute " << resources.compute_threads << "\n"
              << "Planned input shards: "
              << source->planned_shards() << "\n"
              << "Batch size: " << options.batch_size << "\n";

    SiteSelector selector(options, header);
    SampleSelection samples(options, header);
    std::cerr << "Selected samples: " << samples.count() << "\n";
    OrderedCommitter committer(
        options, samples.output_header(header), samples.active());
    std::cerr << "Scheduler: bounded ordered pipeline (3 batches)\n";
    const PipelineSummary summary = run_ordered_pipeline(
        options, selector, samples, header, *source,
        resources.compute_threads, committer);
    committer.finish();

    std::cerr << "After filtering, kept " << summary.kept << " out of "
              << summary.total
              << " sites\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::locale::global(std::locale::classic());
        return run(parse_options(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
