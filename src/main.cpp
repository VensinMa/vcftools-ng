#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
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
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr const char* kVersion = "0.9.0";

struct Options {
    std::string input;
    std::string output_prefix = "out";
    unsigned threads = std::max(1u, std::thread::hardware_concurrency());
    std::size_t batch_size = 2048;

    bool output_freq = false;
    bool output_counts = false;
    bool output_missing_site = false;
    bool output_site_depth = false;
    bool output_site_mean_depth = false;
    bool output_recode = false;
    bool recode_info_all = false;

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

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

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

void print_help() {
    std::cout
        << "vcftools-ng " << kVersion << "\n\n"
        << "Batched exact-compatibility statistics and recode engine.\n\n"
        << "Input:\n"
        << "  --vcf FILE | --gzvcf FILE | --bcf FILE | --input FILE\n"
        << "Execution:\n"
        << "  --out PREFIX\n"
        << "  -t, --threads N\n"
        << "  --batch-size N\n"
        << "  --compat exact\n"
        << "Outputs (may be combined in one scan):\n"
        << "  --freq --counts --missing-site --site-depth --site-mean-depth\n"
        << "  --recode [--recode-INFO-all]\n"
        << "Filters:\n"
        << "  --chr CHROM --not-chr CHROM\n"
        << "  --from-bp POS --to-bp POS\n"
        << "  --positions FILE --exclude-positions FILE\n"
        << "  --bed FILE | --exclude-bed FILE\n"
        << "  --keep-filtered FLAG --remove-filtered FLAG\n"
        << "  --remove-filtered-all\n"
        << "  --keep-INFO FLAG --remove-INFO FLAG\n"
        << "  --indv SAMPLE --remove-indv SAMPLE\n"
        << "  --keep FILE --remove FILE\n"
        << "  --min-alleles N --max-alleles N\n"
        << "  --remove-indels --keep-only-indels --minQ FLOAT\n"
        << "  --minGQ FLOAT --minDP N --maxDP N\n"
        << "  --remove-filtered-geno FLAG --remove-filtered-geno-all\n"
        << "  --min-meanDP FLOAT --max-meanDP FLOAT\n"
        << "  --max-missing FLOAT --max-missing-count N\n"
        << "  --maf FLOAT --max-maf FLOAT\n"
        << "  --mac N --max-mac N --hwe FLOAT\n"
        << "  --non-ref-af FLOAT --max-non-ref-af FLOAT\n"
        << "  --non-ref-af-any FLOAT --max-non-ref-af-any FLOAT\n"
        << "  --non-ref-ac N --max-non-ref-ac N\n"
        << "  --non-ref-ac-any N --max-non-ref-ac-any N\n"
        << "  --thin N\n";
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
        } else if (arg == "--out") {
            options.output_prefix = require_value(argc, argv, i);
        } else if (arg == "--threads" || arg == "-t") {
            options.threads =
                parse_unsigned(require_value(argc, argv, i), arg);
        } else if (arg == "--batch-size") {
            options.batch_size =
                parse_unsigned(require_value(argc, argv, i), arg);
        } else if (arg == "--compat") {
            const auto mode = require_value(argc, argv, i);
            if (mode != "exact") {
                fail("Only --compat exact is implemented");
            }
        } else if (arg == "--input-backend") {
            const auto backend = require_value(argc, argv, i);
            if (backend != "auto") {
                fail("Only --input-backend auto is implemented");
            }
        } else if (arg == "--freq") {
            options.output_freq = true;
        } else if (arg == "--counts") {
            options.output_counts = true;
        } else if (arg == "--missing-site") {
            options.output_missing_site = true;
        } else if (arg == "--site-depth") {
            options.output_site_depth = true;
        } else if (arg == "--site-mean-depth") {
            options.output_site_mean_depth = true;
        } else if (arg == "--recode") {
            options.output_recode = true;
        } else if (arg == "--recode-INFO-all") {
            options.recode_info_all = true;
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
    if (!(options.output_freq || options.output_counts ||
          options.output_missing_site || options.output_site_depth ||
          options.output_site_mean_depth || options.output_recode)) {
        fail("At least one output option is required");
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

class SampleSelection {
public:
    SampleSelection(const Options& options, bcf_hdr_t* input_header)
        : active_(
              !options.samples_to_keep.empty() ||
              !options.samples_to_exclude.empty() ||
              !options.sample_keep_files.empty() ||
              !options.sample_exclude_files.empty()) {
        std::set<std::string> keep = options.samples_to_keep;
        std::set<std::string> exclude = options.samples_to_exclude;
        load_sample_files(options.sample_keep_files, keep);
        load_sample_files(options.sample_exclude_files, exclude);

        const bool has_keep_filter =
            !options.samples_to_keep.empty() ||
            !options.sample_keep_files.empty();
        const int sample_count = bcf_hdr_nsamples(input_header);
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

    bcf_hdr_t* output_header(bcf_hdr_t* input_header) const {
        return active_ ? output_header_.get() : input_header;
    }

    void subset_record(bcf1_t* record) const {
        if (!active_) {
            return;
        }
        if (bcf_subset(
                output_header_.get(), record,
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
    std::string recode_line;
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

double exact_hwe_pvalue(int observed_hets, int observed_hom_ref,
                        int observed_hom_alt) {
    if (observed_hom_ref + observed_hom_alt + observed_hets == 0) {
        return 1.0;
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

    double p_hwe = 0.0;
    for (const double probability : heterozygote_probabilities) {
        if (probability <=
            heterozygote_probabilities[observed_hets]) {
            p_hwe += probability;
        }
    }
    return std::min(p_hwe, 1.0);
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

std::string format_recode_line(const Options& options, bcf_hdr_t* header,
                               bcf1_t* record, Scratch& scratch) {
    scratch.recode_buffer.l = 0;
    if (vcf_format1(header, record, &scratch.recode_buffer) != 0) {
        fail("Could not format recoded VCF record");
    }
    std::string line(
        scratch.recode_buffer.s, scratch.recode_buffer.l);
    if (record->d.n_flt > 1) {
        sort_filter_column(line);
    }
    if (!options.recode_info_all) {
        remove_info_column(line);
    }
    if (line.empty() || line.back() != '\n') {
        line.push_back('\n');
    }
    return line;
}

SiteResult process_site(const Options& options,
                        const SiteSelector& selector,
                        const SampleSelection& samples,
                        bcf_hdr_t* header, bcf1_t* record,
                        Scratch& scratch) {
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
    const bool need_gt = options.output_freq || options.output_counts ||
                         options.output_missing_site ||
                         frequency_filter_active ||
                         (options.output_recode && genotype_filter_active);
    const bool need_dp = options.output_site_depth ||
                         options.output_site_mean_depth ||
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
    if (need_gt) {
        for (const int sample : samples.indices()) {
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
                continue;
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

            if (options.min_hwe > 0.0 && actual_ploidy >= 2) {
                const int32_t first = scratch.gt[base];
                const int32_t second = scratch.gt[base + 1];
                if (!bcf_gt_is_missing(first) &&
                    !bcf_gt_is_missing(second)) {
                    const int first_allele = bcf_gt_allele(first);
                    const int second_allele = bcf_gt_allele(second);
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
        for (const int sample : samples.indices()) {
            if (genotype_filtered[sample]) {
                continue;
            }
            const int depth =
                legacy_scalar_value(scratch.dp, dp_count, sample, sample_count);
            if (depth >= 0) {
                result.sum_depth += static_cast<uint32_t>(depth);
                result.sumsq_depth +=
                    static_cast<uint32_t>(depth * depth);
                ++result.depth_count;
            }
        }
    }

    result.chrom = bcf_hdr_id2name(header, record->rid);
    result.pos = record->pos + 1;
    result.alleles.reserve(record->n_allele);
    for (int allele = 0; allele < record->n_allele; ++allele) {
        result.alleles.emplace_back(record->d.allele[allele]);
    }
    if (options.output_recode) {
        if (genotype_filter_active) {
            mask_filtered_genotypes(
                header, record, genotype_filtered, scratch.gt, gt_count,
                sample_count, max_ploidy);
        }
        samples.subset_record(record);
        result.recode_line =
            format_recode_line(
                options, samples.output_header(header), record, scratch);
    }
    result.kept = true;
    return result;
}

class Outputs {
public:
    Outputs(const Options& options, bcf_hdr_t* header) : options_(options) {
        if (options.output_freq) {
            freq_ = open(".frq");
            freq_ << "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\n";
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
        if (options.output_recode) {
            recode_ = open(".recode.vcf");
            kstring_t header_text{0, 0, nullptr};
            if (bcf_hdr_format(header, 0, &header_text) != 0) {
                std::free(header_text.s);
                fail("Could not format VCF header");
            }
            recode_.write(
                header_text.s,
                static_cast<std::streamsize>(header_text.l));
            std::free(header_text.s);
            if (!recode_) {
                fail("Could not write recoded VCF header");
            }
        }
    }

    void write(const SiteResult& result) {
        if (!result.kept) {
            return;
        }
        if (options_.output_freq) {
            freq_ << result.chrom << '\t' << result.pos << '\t'
                  << result.alleles.size() << '\t'
                  << result.non_missing_chromosomes;
            for (std::size_t i = 0; i < result.alleles.size(); ++i) {
                const double frequency =
                    result.allele_counts[i] /
                    static_cast<double>(result.non_missing_chromosomes);
                freq_ << '\t' << result.alleles[i] << ':' << frequency;
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
        if (options_.output_recode) {
            recode_.write(
                result.recode_line.data(),
                static_cast<std::streamsize>(result.recode_line.size()));
        }
    }

private:
    std::ofstream open(const std::string& suffix) {
        const std::string path = options_.output_prefix + suffix;
        std::ofstream stream(path);
        if (!stream) {
            fail("Could not open output file: " + path);
        }
        return stream;
    }

    const Options& options_;
    std::ofstream freq_;
    std::ofstream counts_;
    std::ofstream missing_;
    std::ofstream depth_;
    std::ofstream mean_depth_;
    std::ofstream recode_;
};

struct PipelineSummary {
    uint64_t total = 0;
    uint64_t kept = 0;
};

class OrderedCommitter {
public:
    OrderedCommitter(const Options& options, bcf_hdr_t* output_header)
        : outputs_(options, output_header),
          thin_selector_(options.thin_distance) {}

    void commit(std::vector<SiteResult>& results) {
        for (auto& result : results) {
            if (result.kept && !thin_selector_.keep(result)) {
                result.kept = false;
            }
            outputs_.write(result);
            ++summary_.total;
            summary_.kept += result.kept;
        }
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
    std::vector<RecordPtr> records;
    std::vector<SiteResult> results;
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
    const SampleSelection& samples, bcf_hdr_t* header, htsFile* input,
    hts_itr_t* iterator, OrderedCommitter& committer) {
    constexpr std::size_t max_in_flight_batches = 3;
    const std::size_t slices_per_worker = 4;
    const std::size_t target_slice =
        options.batch_size /
        std::max<std::size_t>(
            1, static_cast<std::size_t>(options.threads) *
                   slices_per_worker);
    const std::size_t slice_size =
        std::clamp<std::size_t>(target_slice, 64, 256);
    PipelineState state;

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
            RecordPtr reusable(bcf_init());
            if (!reusable) {
                fail("Could not allocate HTSlib reader record");
            }
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
                batch->records.reserve(options.batch_size);
                for (std::size_t index = 0;
                     index < options.batch_size; ++index) {
                    const int status =
                        iterator != nullptr
                            ? bcf_itr_next(input, iterator, reusable.get())
                            : bcf_read(input, header, reusable.get());
                    if (status != 0) {
                        eof = true;
                        break;
                    }
                    bcf1_t* duplicate = bcf_dup(reusable.get());
                    if (!duplicate) {
                        fail("Could not duplicate HTSlib record");
                    }
                    batch->records.emplace_back(duplicate);
                }

                if (batch->records.empty()) {
                    std::lock_guard lock(state.mutex);
                    state.reader_done = true;
                    state.work_available.notify_all();
                    state.completed_available.notify_all();
                    return;
                }

                batch->results.resize(batch->records.size());
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
    workers.reserve(options.threads);
    for (unsigned worker_index = 0; worker_index < options.threads;
         ++worker_index) {
        workers.emplace_back([&] {
            Scratch scratch;
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
                            slice.batch->records[index].get(), scratch);
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

            committer.commit(batch->results);
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

int run(const Options& options) {
    HtsFilePtr input(hts_open(options.input.c_str(), "r"));
    if (!input) {
        fail("Could not open input file: " + options.input);
    }

    const int io_threads =
        options.threads > 1 ? std::min(4u, options.threads - 1) : 0;
    if (io_threads > 0 && hts_set_threads(input.get(), io_threads) != 0) {
        fail("Could not enable HTSlib input threads");
    }

    HeaderPtr header(bcf_hdr_read(input.get()));
    if (!header) {
        fail("Could not read VCF/BCF header");
    }
    validate_info_flag_filters(options, header.get());

    std::cerr << "vcftools-ng " << kVersion << "\n"
              << "Input: " << options.input << "\n"
              << "Samples: " << bcf_hdr_nsamples(header.get()) << "\n"
              << "Threads: " << options.threads
              << " (HTSlib I/O: " << io_threads << ")\n"
              << "Batch size: " << options.batch_size << "\n";

    SiteSelector selector(options, header.get());
    SampleSelection samples(options, header.get());
    std::cerr << "Selected samples: " << samples.count() << "\n";
    IndexPtr index;
    IteratorPtr iterator;
    if (options.chromosomes_to_keep.size() == 1 &&
        hts_get_format(input.get())->format == bcf) {
        index.reset(bcf_index_load(options.input.c_str()));
        if (index) {
            std::string region = *options.chromosomes_to_keep.begin();
            if (options.start_position != -1 ||
                options.end_position !=
                    std::numeric_limits<int>::max()) {
                const int start =
                    options.start_position == -1
                        ? 1
                        : options.start_position;
                region += ":" + std::to_string(start);
                if (options.end_position !=
                    std::numeric_limits<int>::max()) {
                    region += "-" +
                              std::to_string(options.end_position);
                }
            }
            iterator.reset(
                bcf_itr_querys(index.get(), header.get(), region.c_str()));
            if (!iterator) {
                fail("Could not query BCF region: " + region);
            }
            std::cerr << "Indexed region: " << region << "\n";
        }
    }
    OrderedCommitter committer(
        options, samples.output_header(header.get()));
    std::cerr << "Scheduler: bounded ordered pipeline (3 batches)\n";
    const PipelineSummary summary = run_ordered_pipeline(
        options, selector, samples, header.get(), input.get(),
        iterator.get(), committer);

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
