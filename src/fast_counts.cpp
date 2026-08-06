#include "fast_counts.h"
#include "output_transaction.h"

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <unistd.h>

namespace vcftools_ng {
namespace {

constexpr std::uint64_t kMib = 1024ULL * 1024ULL;
constexpr std::uint64_t kMinimumShardBytes = 1ULL * kMib;
constexpr std::uint64_t kMaximumShardBytes = 64ULL * kMib;
constexpr std::uint64_t kTargetActiveInputBytes = 256ULL * kMib;
constexpr std::size_t kMaximumShards = 65536;
constexpr std::size_t kOutputFlushBytes = 8ULL * 1024ULL * 1024ULL;
constexpr std::size_t kInlineAlleles = 16;
constexpr std::size_t kAdvancedInlineAlleles = 4;
constexpr std::size_t kArtifactCount = 6;
constexpr std::string_view kFreqHeader =
    "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\n";
constexpr std::string_view kFreq2Header =
    "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\n";
constexpr std::string_view kCountsHeader =
    "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\n";
constexpr std::string_view kMissingHeader =
    "CHR\tPOS\tN_DATA\tN_GENOTYPE_FILTERED\tN_MISS\tF_MISS\n";
constexpr std::string_view kDepthHeader =
    "CHROM\tPOS\tSUM_DEPTH\tSUMSQ_DEPTH\n";
constexpr std::string_view kMeanDepthHeader =
    "CHROM\tPOS\tMEAN_DEPTH\tVAR_DEPTH\n";
constexpr std::string_view kQualityHeader =
    "CHROM\tPOS\tQUAL\n";

enum class Artifact : std::size_t {
    freq = 0,
    counts = 1,
    missing = 2,
    depth = 3,
    mean_depth = 4,
    quality = 5,
};

constexpr std::size_t artifact_index(Artifact artifact) {
    return static_cast<std::size_t>(artifact);
}

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

struct HtsFileDeleter {
    void operator()(htsFile* input) const noexcept {
        if (input != nullptr) {
            hts_close(input);
        }
    }
};

using HtsFilePtr = std::unique_ptr<htsFile, HtsFileDeleter>;

struct FileDescriptor {
    explicit FileDescriptor(int value) : value(value) {}
    ~FileDescriptor() {
        if (value >= 0) {
            ::close(value);
        }
    }
    FileDescriptor(const FileDescriptor&) = delete;
    FileDescriptor& operator=(const FileDescriptor&) = delete;
    int value = -1;
};

struct ReadOnlyMapping {
    ReadOnlyMapping(
        const std::string& path, std::uint64_t size, bool enabled) {
        if (!enabled || size == 0 ||
            size > static_cast<std::uint64_t>(
                       std::numeric_limits<std::size_t>::max())) {
            return;
        }
        FileDescriptor descriptor(::open(path.c_str(), O_RDONLY));
        if (descriptor.value < 0) {
            return;
        }
        void* mapped = ::mmap(
            nullptr, static_cast<std::size_t>(size), PROT_READ,
            MAP_PRIVATE, descriptor.value, 0);
        if (mapped != MAP_FAILED) {
            data = static_cast<const char*>(mapped);
            length = static_cast<std::size_t>(size);
        }
    }

    ~ReadOnlyMapping() {
        if (data != nullptr) {
            ::munmap(const_cast<char*>(data), length);
        }
    }

    ReadOnlyMapping(const ReadOnlyMapping&) = delete;
    ReadOnlyMapping& operator=(const ReadOnlyMapping&) = delete;

    const char* data = nullptr;
    std::size_t length = 0;
};

struct HeaderDeleter {
    void operator()(bcf_hdr_t* header) const noexcept {
        if (header != nullptr) {
            bcf_hdr_destroy(header);
        }
    }
};

struct TbxDeleter {
    void operator()(tbx_t* index) const noexcept {
        if (index != nullptr) {
            tbx_destroy(index);
        }
    }
};

struct IteratorDeleter {
    void operator()(hts_itr_t* iterator) const noexcept {
        if (iterator != nullptr) {
            hts_itr_destroy(iterator);
        }
    }
};

using HeaderPtr = std::unique_ptr<bcf_hdr_t, HeaderDeleter>;
using TbxPtr = std::unique_ptr<tbx_t, TbxDeleter>;
using IteratorPtr = std::unique_ptr<hts_itr_t, IteratorDeleter>;

struct KStringBuffer {
    kstring_t value{0, 0, nullptr};

    ~KStringBuffer() {
        std::free(value.s);
    }
};

struct HeaderLayout {
    std::uint64_t data_start = 0;
    std::size_t samples = 0;
    std::vector<std::string> sample_names;
};

std::string read_text_vcf_header(
    htsFile* input, const std::string& path) {
    KStringBuffer line;
    std::string header;
    while (hts_getline(input, '\n', &line.value) >= 0) {
        std::string_view text(line.value.s, line.value.l);
        if (!text.empty() && text.back() == '\r') {
            text.remove_suffix(1);
        }
        if (text.empty() || text.front() != '#') {
            fail("VCF data appeared before #CHROM header: " + path);
        }
        header.append(text);
        header.push_back('\n');
        if (text.rfind("#CHROM\t", 0) == 0) {
            return header;
        }
    }
    fail("Could not find #CHROM header: " + path);
}

struct PlainShard {
    std::uint64_t begin = 0;
    std::uint64_t end = 0;
};

struct IndexedShard {
    int tid = -1;
    hts_pos_t begin = 0;
    hts_pos_t end = HTS_POS_MAX;
};

struct ShardOutput {
    std::array<std::string, kArtifactCount> text;
    std::string recode_text;
    std::vector<std::uint8_t> genotype_filtered;
    struct AdvancedSite {
        std::string chromosome;
        int position = 0;
        std::uint32_t allele_count = 0;
        std::array<std::uint32_t, kAdvancedInlineAlleles>
            inline_allele_counts{};
        std::vector<std::uint32_t> heap_allele_counts;
        std::uint32_t non_missing_chromosomes = 0;
        bool fully_diploid = true;
        bool fst_eligible = false;
        double fst_sum_a = 0.0;
        double fst_sum_all = 0.0;
        double fst = std::numeric_limits<double>::quiet_NaN();

        [[nodiscard]] std::uint32_t allele_count_at(
            std::size_t index) const {
            return heap_allele_counts.empty()
                       ? inline_allele_counts[index]
                       : heap_allele_counts[index];
        }
    };
    std::vector<AdvancedSite> advanced_sites;
    std::uint64_t records = 0;
    std::uint64_t kept = 0;
};

class AdvancedOutputs {
public:
    AdvancedOutputs(
        const std::string& prefix, const FastSiteStatPlan& plan,
        std::size_t selected_samples)
        : plan_(plan),
          selected_chromosomes_(
              static_cast<std::uint64_t>(selected_samples) * 2) {
        if (plan_.pi_window_size > 0) {
            pi_ = open(prefix + ".windowed.pi");
            pi_ << "CHROM\tBIN_START\tBIN_END\tN_VARIANTS"
                   "\tN_MONOMORPHIC\tPI\n";
            pi_step_ =
                plan_.pi_window_step <= 0 ||
                        plan_.pi_window_step > plan_.pi_window_size
                    ? plan_.pi_window_size
                    : plan_.pi_window_step;
        }
        if (plan_.tajima_window_size > 0) {
            tajima_ = open(prefix + ".Tajima.D");
            tajima_ << "CHROM\tBIN_START\tN_SNPS\tTajimaD\n";
            initialise_tajima_constants();
        }
        if (!plan_.fst_population_files.empty()) {
            if (plan_.fst_window_size > 0) {
                window_fst_ = open(prefix + ".windowed.weir.fst");
                window_fst_
                    << "CHROM\tBIN_START\tBIN_END\tN_VARIANTS"
                       "\tWEIGHTED_FST\tMEAN_FST\n";
                fst_step_ =
                    plan_.fst_window_step <= 0 ||
                            plan_.fst_window_step > plan_.fst_window_size
                        ? plan_.fst_window_size
                        : plan_.fst_window_step;
            } else {
                site_fst_ = open(prefix + ".weir.fst");
                site_fst_
                    << "CHROM\tPOS\tWEIR_AND_COCKERHAM_FST\n";
            }
        }
    }

    void consume(const std::vector<ShardOutput::AdvancedSite>& sites) {
        for (const auto& site : sites) {
            consume(site);
        }
    }

    void finish(const std::string& prefix) {
        finish_pi();
        finish_tajima();
        finish_window_fst();
        finish_stream(pi_, prefix + ".windowed.pi",
                      plan_.pi_window_size > 0);
        finish_stream(tajima_, prefix + ".Tajima.D",
                      plan_.tajima_window_size > 0);
        finish_stream(site_fst_, prefix + ".weir.fst",
                      !plan_.fst_population_files.empty() &&
                          plan_.fst_window_size <= 0);
        finish_stream(window_fst_, prefix + ".windowed.weir.fst",
                      plan_.fst_window_size > 0);
    }

private:
    using PiWindow = std::array<std::uint64_t, 4>;

    static std::ofstream open(const std::string& path) {
        return vcftools_ng::output::open_stream(
            path, std::ios::binary | std::ios::trunc);
    }

    static void finish_stream(
        std::ofstream& stream, const std::string& path, bool enabled) {
        if (enabled) {
            vcftools_ng::output::finish_stream(stream, path);
        }
    }

    void consume(const ShardOutput::AdvancedSite& site) {
        if (plan_.pi_window_size > 0 && site.fully_diploid) {
            std::uint64_t mismatches = 0;
            for (std::size_t allele = 0;
                 allele < site.allele_count; ++allele) {
                const std::uint32_t count =
                    site.allele_count_at(allele);
                mismatches += static_cast<std::uint64_t>(count) *
                              (site.non_missing_chromosomes - count);
            }
            if (mismatches > 0) {
                consume_pi(site, mismatches);
            }
        }
        if (plan_.tajima_window_size > 0 &&
            site.allele_count == 2 && site.fully_diploid) {
            consume_tajima(site);
        }
        if (site.fst_eligible) {
            if (plan_.fst_window_size > 0) {
                if (!std::isnan(site.fst)) {
                    consume_window_fst(site);
                }
            } else if (!plan_.fst_population_files.empty()) {
                site_fst_ << site.chromosome << '\t' << site.position
                          << '\t' << site.fst << '\n';
            }
        }
    }

    void consume_pi(
        const ShardOutput::AdvancedSite& site,
        std::uint64_t mismatches) {
        const int first = std::max(
            0, static_cast<int>(std::ceil(
                   (site.position - plan_.pi_window_size) /
                   static_cast<double>(pi_step_))));
        const int last = static_cast<int>(std::ceil(
            site.position / static_cast<double>(pi_step_)));
        if (site.chromosome != previous_pi_chromosome_) {
            pi_chromosomes_.push_back(site.chromosome);
            previous_pi_chromosome_ = site.chromosome;
            pi_windows_[site.chromosome].resize(1, PiWindow{0, 0, 0, 0});
        }
        auto& windows = pi_windows_[site.chromosome];
        if (last >= static_cast<int>(windows.size())) {
            windows.resize(
                static_cast<std::size_t>(last + 1), PiWindow{0, 0, 0, 0});
        }
        const std::uint64_t comparisons =
            static_cast<std::uint64_t>(site.non_missing_chromosomes) *
            (site.non_missing_chromosomes - 1);
        for (int index = first; index < last; ++index) {
            auto& window = windows[static_cast<std::size_t>(index)];
            ++window[0];
            window[1] += comparisons;
            window[2] += mismatches;
            if (site.allele_count_at(0) <
                site.non_missing_chromosomes) {
                ++window[3];
            }
        }
    }

    void finish_pi() {
        if (plan_.pi_window_size <= 0) {
            return;
        }
        const std::uint64_t monomorphic_comparisons =
            selected_chromosomes_ * (selected_chromosomes_ - 1);
        for (const auto& chromosome : pi_chromosomes_) {
            const auto& windows = pi_windows_[chromosome];
            for (std::size_t index = 0; index < windows.size(); ++index) {
                const auto& window = windows[index];
                if (window[3] == 0 && window[2] == 0) {
                    continue;
                }
                const std::uint64_t monomorphic_sites =
                    static_cast<std::uint64_t>(plan_.pi_window_size) -
                    window[0];
                const std::uint64_t pairs =
                    window[1] +
                    monomorphic_sites * monomorphic_comparisons;
                pi_ << chromosome << '\t' << index * pi_step_ + 1 << '\t'
                    << index * pi_step_ + plan_.pi_window_size << '\t'
                    << window[3] << '\t' << monomorphic_sites << '\t'
                    << window[2] / static_cast<double>(pairs) << '\n';
            }
        }
    }

    void initialise_tajima_constants() {
        const std::uint64_t n = selected_chromosomes_;
        for (std::uint64_t index = 1; index < n; ++index) {
            tajima_a1_ += 1.0 / static_cast<double>(index);
            tajima_a2_ += 1.0 / static_cast<double>(index * index);
        }
        const double b1 = static_cast<double>(n + 1) / 3.0 /
                          static_cast<double>(n - 1);
        const double b2 =
            2.0 * static_cast<double>(n * n + n + 3) / 9.0 /
            static_cast<double>(n) / static_cast<double>(n - 1);
        const double c1 = b1 - (1.0 / tajima_a1_);
        const double c2 =
            b2 - (static_cast<double>(n + 2) /
                  static_cast<double>(tajima_a1_ * n)) +
            (tajima_a2_ / tajima_a1_ / tajima_a1_);
        tajima_e1_ = c1 / tajima_a1_;
        tajima_e2_ =
            c2 / ((tajima_a1_ * tajima_a1_) + tajima_a2_);
    }

    void consume_tajima(const ShardOutput::AdvancedSite& site) {
        const std::size_t index = static_cast<std::size_t>(
            site.position /
            static_cast<double>(plan_.tajima_window_size));
        if (site.chromosome != previous_tajima_chromosome_) {
            tajima_chromosomes_.push_back(site.chromosome);
            previous_tajima_chromosome_ = site.chromosome;
        }
        auto& windows = tajima_windows_[site.chromosome];
        if (index >= windows.size()) {
            windows.resize(index + 1, std::pair<int, double>{0, 0.0});
        }
        const double frequency =
            site.allele_count_at(0) /
            static_cast<double>(site.non_missing_chromosomes);
        if (frequency > 0.0 && frequency < 1.0) {
            ++windows[index].first;
            windows[index].second += frequency * (1.0 - frequency);
        }
    }

    void finish_tajima() {
        if (plan_.tajima_window_size <= 0) {
            return;
        }
        const double n = static_cast<double>(selected_chromosomes_);
        for (const auto& chromosome : tajima_chromosomes_) {
            bool output = false;
            const auto& windows = tajima_windows_[chromosome];
            for (std::size_t index = 0; index < windows.size(); ++index) {
                const int segregating_sites = windows[index].first;
                double tajima_d =
                    std::numeric_limits<double>::quiet_NaN();
                if (segregating_sites > 0) {
                    const double pi = 2.0 * windows[index].second * n /
                                      (n - 1.0);
                    const double theta = segregating_sites / tajima_a1_;
                    const double variance =
                        tajima_e1_ * segregating_sites +
                        tajima_e2_ * segregating_sites *
                            (segregating_sites - 1);
                    tajima_d = (pi - theta) / std::sqrt(variance);
                    output = true;
                }
                if (output) {
                    tajima_ << chromosome << '\t'
                            << index * plan_.tajima_window_size << '\t'
                            << segregating_sites << '\t' << tajima_d
                            << '\n';
                }
            }
        }
    }

    void consume_window_fst(const ShardOutput::AdvancedSite& site) {
        if (site.chromosome != previous_fst_chromosome_) {
            fst_chromosomes_.push_back(site.chromosome);
            previous_fst_chromosome_ = site.chromosome;
        }
        const int first = std::max(
            0, static_cast<int>(std::ceil(
                   (site.position - plan_.fst_window_size) /
                   static_cast<double>(fst_step_))));
        const int last = static_cast<int>(std::ceil(
            site.position / static_cast<double>(fst_step_)));
        auto& windows = fst_windows_[site.chromosome];
        for (int index = first; index < last; ++index) {
            if (index >= static_cast<int>(windows.size())) {
                windows.resize(
                    static_cast<std::size_t>(index + 1),
                    std::array<double, 4>{0.0, 0.0, 0.0, 0.0});
            }
            auto& window = windows[static_cast<std::size_t>(index)];
            window[0] += site.fst_sum_a;
            window[1] += site.fst_sum_all;
            window[2] += site.fst;
            ++window[3];
        }
    }

    void finish_window_fst() {
        if (plan_.fst_window_size <= 0) {
            return;
        }
        for (const auto& chromosome : fst_chromosomes_) {
            const auto& windows = fst_windows_[chromosome];
            for (std::size_t index = 0; index < windows.size(); ++index) {
                const auto& window = windows[index];
                if (window[1] != 0.0 && !std::isnan(window[0]) &&
                    !std::isnan(window[1]) && window[3] > 0.0) {
                    window_fst_
                        << chromosome << '\t' << index * fst_step_ + 1
                        << '\t' << index * fst_step_ + plan_.fst_window_size
                        << '\t' << window[3] << '\t'
                        << window[0] / window[1] << '\t'
                        << window[2] / window[3] << '\n';
                }
            }
        }
    }

    const FastSiteStatPlan& plan_;
    std::uint64_t selected_chromosomes_ = 0;
    std::ofstream pi_;
    std::ofstream tajima_;
    std::ofstream site_fst_;
    std::ofstream window_fst_;
    int pi_step_ = 0;
    int fst_step_ = 0;
    std::map<std::string, std::vector<PiWindow>, std::less<>> pi_windows_;
    std::vector<std::string> pi_chromosomes_;
    std::string previous_pi_chromosome_;
    std::map<std::string, std::vector<std::pair<int, double>>, std::less<>>
        tajima_windows_;
    std::vector<std::string> tajima_chromosomes_;
    std::string previous_tajima_chromosome_;
    double tajima_a1_ = 0.0;
    double tajima_a2_ = 0.0;
    double tajima_e1_ = 0.0;
    double tajima_e2_ = 0.0;
    std::map<std::string, std::vector<std::array<double, 4>>, std::less<>>
        fst_windows_;
    std::vector<std::string> fst_chromosomes_;
    std::string previous_fst_chromosome_;
};

class OrderedArtifactSet {
public:
    OrderedArtifactSet(
        const std::string& prefix,
        const FastSiteStatPlan& plan,
        std::size_t selected_samples)
        : prefix_(prefix), plan_(plan),
          advanced_(prefix, plan, selected_samples) {
        if (plan.recode && !plan.recode_sink) {
            fail("Fast recode output has no destination");
        }
        if (plan.freq || plan.freq2) {
            open(
                Artifact::freq, prefix + ".frq",
                plan.freq ? kFreqHeader : kFreq2Header);
        }
        if (plan.counts) {
            open(
                Artifact::counts, prefix + ".frq.count",
                kCountsHeader);
        }
        if (plan.missing_site) {
            open(
                Artifact::missing, prefix + ".lmiss",
                kMissingHeader);
        }
        if (plan.site_depth) {
            open(
                Artifact::depth, prefix + ".ldepth",
                kDepthHeader);
        }
        if (plan.site_mean_depth) {
            open(
                Artifact::mean_depth, prefix + ".ldepth.mean",
                kMeanDepthHeader);
        }
        if (plan.site_quality) {
            open(
                Artifact::quality, prefix + ".lqual",
                kQualityHeader);
        }
    }

    void append(const ShardOutput& shard) {
        for (std::size_t index = 0; index < outputs_.size(); ++index) {
            if (!enabled_[index]) {
                continue;
            }
            outputs_[index].write(
                shard.text[index].data(),
                static_cast<std::streamsize>(
                    shard.text[index].size()));
        }
        if (plan_.recode && !shard.recode_text.empty()) {
            plan_.recode_sink(shard.recode_text);
        }
        advanced_.consume(shard.advanced_sites);
    }

    void validate() {
        for (std::size_t index = 0; index < outputs_.size(); ++index) {
            if (enabled_[index]) {
                vcftools_ng::output::finish_stream(
                    outputs_[index], final_paths_[index]);
            }
        }
        advanced_.finish(prefix_);
    }

private:
    void open(
        Artifact artifact, const std::string& path,
        std::string_view header) {
        const std::size_t index = artifact_index(artifact);
        final_paths_[index] = path;
        outputs_[index] = vcftools_ng::output::open_stream(
            path, std::ios::binary | std::ios::trunc);
        enabled_[index] = true;
        outputs_[index].write(
            header.data(),
            static_cast<std::streamsize>(header.size()));
    }

    std::array<std::ofstream, kArtifactCount> outputs_;
    std::array<std::string, kArtifactCount> final_paths_;
    std::array<bool, kArtifactCount> enabled_{};
    std::string prefix_;
    const FastSiteStatPlan& plan_;
    AdvancedOutputs advanced_;
};

unsigned descriptor_limited_workers(
    unsigned requested, std::size_t shards);

std::size_t count_tabs(std::string_view text) {
    return static_cast<std::size_t>(
        std::count(text.begin(), text.end(), '\t'));
}

HeaderLayout read_plain_header(
    const std::string& path, std::uint64_t file_size) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        fail("Could not open plain VCF: " + path);
    }
    std::string line;
    while (std::getline(input, line)) {
        if (line.rfind("#CHROM\t", 0) == 0) {
            const auto position = input.tellg();
            HeaderLayout layout;
            layout.data_start =
                position < 0
                    ? file_size
                    : static_cast<std::uint64_t>(position);
            std::istringstream fields(line);
            std::string field;
            for (std::size_t column = 0; column < 9; ++column) {
                if (!(fields >> field)) {
                    break;
                }
            }
            while (fields >> field) {
                layout.sample_names.push_back(field);
            }
            layout.samples = layout.sample_names.size();
            return layout;
        }
        if (line.empty() || line.front() != '#') {
            fail("VCF data appeared before #CHROM header: " + path);
        }
    }
    fail("Could not find #CHROM header: " + path);
}

void load_position_file(
    const std::string& path,
    std::map<std::string, std::vector<int>, std::less<>>& positions) {
    if (path.empty()) {
        return;
    }
    std::ifstream input(path);
    if (!input) {
        fail("Could not open positions file: " + path);
    }
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') {
            continue;
        }
        std::istringstream fields(line);
        std::string chromosome;
        int position = 0;
        if (!(fields >> chromosome >> position)) {
            fail("Invalid positions line in " + path);
        }
        positions[std::move(chromosome)].push_back(position);
    }
    for (auto& [chromosome, values] : positions) {
        (void)chromosome;
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end());
    }
}

void load_sample_files(
    const std::vector<std::string>& paths,
    std::set<std::string>& samples) {
    for (const auto& path : paths) {
        std::ifstream input(path);
        if (!input) {
            fail("Could not open sample file: " + path);
        }
        std::string line;
        while (std::getline(input, line)) {
            std::istringstream fields(line);
            std::string sample;
            if (fields >> sample) {
                samples.insert(std::move(sample));
            }
        }
    }
}

std::vector<std::size_t> selected_sample_indices(
    const FastSiteStatPlan& plan, const HeaderLayout& header) {
    if (!plan.sample_selection_active()) {
        return {};
    }
    std::set<std::string> keep = plan.samples_to_keep;
    std::set<std::string> exclude = plan.samples_to_exclude;
    load_sample_files(plan.sample_keep_files, keep);
    load_sample_files(plan.sample_exclude_files, exclude);
    const bool has_keep =
        !plan.samples_to_keep.empty() ||
        !plan.sample_keep_files.empty();
    std::vector<std::size_t> selected;
    selected.reserve(header.samples);
    for (std::size_t index = 0;
         index < header.sample_names.size(); ++index) {
        const auto& sample = header.sample_names[index];
        if (has_keep && !keep.contains(sample)) {
            continue;
        }
        if (exclude.contains(sample)) {
            continue;
        }
        selected.push_back(index);
    }
    if (selected.empty()) {
        fail("Sample filters removed all individuals");
    }
    if (selected.size() == header.samples) {
        return {};
    }
    return selected;
}

std::string subset_text_vcf_header(
    std::string_view header,
    const std::vector<std::size_t>& selected_samples) {
    if (selected_samples.empty()) {
        return std::string(header);
    }
    const std::size_t chrom_begin = header.rfind("#CHROM\t");
    if (chrom_begin == std::string_view::npos) {
        fail("Could not find #CHROM line while subsetting VCF header");
    }
    const std::size_t chrom_end = header.find('\n', chrom_begin);
    const std::string_view line = header.substr(
        chrom_begin,
        chrom_end == std::string_view::npos
            ? std::string_view::npos
            : chrom_end - chrom_begin);
    std::array<std::string_view, 9> fixed{};
    std::size_t begin = 0;
    for (std::size_t column = 0; column < fixed.size(); ++column) {
        const std::size_t end = line.find('\t', begin);
        if (end == std::string_view::npos) {
            fail("VCF #CHROM header has fewer than nine fixed columns");
        }
        fixed[column] = line.substr(begin, end - begin);
        begin = end + 1;
    }

    std::string output;
    output.reserve(header.size());
    output.append(header.substr(0, chrom_begin));
    for (std::size_t column = 0; column < fixed.size(); ++column) {
        if (column != 0) {
            output.push_back('\t');
        }
        output.append(fixed[column]);
    }
    std::size_t selected_cursor = 0;
    for (std::size_t sample = 0;
         selected_cursor < selected_samples.size(); ++sample) {
        const std::size_t end = line.find('\t', begin);
        const std::string_view value = line.substr(
            begin,
            end == std::string_view::npos
                ? std::string_view::npos
                : end - begin);
        if (sample == selected_samples[selected_cursor]) {
            output.push_back('\t');
            output.append(value);
            ++selected_cursor;
        }
        if (end == std::string_view::npos) {
            break;
        }
        begin = end + 1;
    }
    if (selected_cursor != selected_samples.size()) {
        fail("Selected sample index is absent from VCF #CHROM header");
    }
    output.push_back('\n');
    if (chrom_end != std::string_view::npos &&
        chrom_end + 1 < header.size()) {
        output.append(header.substr(chrom_end + 1));
    }
    return output;
}

struct PositionSelectionCursor {
    std::string chromosome;
    const std::vector<int>* keep = nullptr;
    const std::vector<int>* exclude = nullptr;
};

void update_position_cursor(
    PositionSelectionCursor& cursor, std::string_view chromosome,
    const FastSiteStatPlan& plan) {
    if (cursor.chromosome == chromosome) {
        return;
    }
    cursor.chromosome.assign(chromosome);
    const auto keep = plan.positions_to_keep.find(chromosome);
    cursor.keep = keep == plan.positions_to_keep.end()
                      ? nullptr
                      : &keep->second;
    const auto exclude = plan.positions_to_exclude.find(chromosome);
    cursor.exclude = exclude == plan.positions_to_exclude.end()
                         ? nullptr
                         : &exclude->second;
}

bool contains_position(const std::vector<int>* positions, int position) {
    return positions != nullptr &&
           std::binary_search(
               positions->begin(), positions->end(), position);
}

std::uint64_t align_plain_offset(
    std::ifstream& input, std::uint64_t offset,
    std::uint64_t data_start, std::uint64_t file_size) {
    if (offset <= data_start) {
        return data_start;
    }
    if (offset >= file_size) {
        return file_size;
    }
    input.clear();
    input.seekg(static_cast<std::streamoff>(offset - 1));
    char previous = '\0';
    input.get(previous);
    if (previous == '\n') {
        return offset;
    }
    input.seekg(static_cast<std::streamoff>(offset));
    std::string remainder;
    std::getline(input, remainder);
    const auto position = input.tellg();
    return position < 0
               ? file_size
               : static_cast<std::uint64_t>(position);
}

std::vector<PlainShard> build_plain_shards(
    const std::string& path, std::uint64_t data_start,
    std::uint64_t file_size, unsigned requested_threads) {
    if (file_size <= data_start) {
        return {};
    }
    const std::uint64_t bytes = file_size - data_start;
    const std::uint64_t shard_bytes = std::clamp<std::uint64_t>(
        kTargetActiveInputBytes /
            std::max<std::uint64_t>(1, requested_threads),
        kMinimumShardBytes, kMaximumShardBytes);
    const std::size_t count = std::min<std::size_t>(
        kMaximumShards,
        static_cast<std::size_t>(
            (bytes + shard_bytes - 1) / shard_bytes));
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        fail("Could not align plain VCF shards: " + path);
    }
    std::vector<std::uint64_t> boundaries(count + 1);
    boundaries.front() = data_start;
    boundaries.back() = file_size;
    for (std::size_t index = 1; index < count; ++index) {
        boundaries[index] = align_plain_offset(
            input,
            data_start +
                static_cast<std::uint64_t>(
                    (static_cast<long double>(bytes) * index) /
                    count),
            data_start, file_size);
    }
    std::vector<PlainShard> shards;
    shards.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        if (boundaries[index] < boundaries[index + 1]) {
            shards.push_back(
                PlainShard{boundaries[index], boundaries[index + 1]});
        }
    }
    return shards;
}

void append_unsigned(std::string& output, std::uint64_t value) {
    std::array<char, 32> buffer{};
    const auto converted = std::to_chars(
        buffer.data(), buffer.data() + buffer.size(), value);
    if (converted.ec != std::errc{}) {
        fail("Could not format allele count");
    }
    output.append(buffer.data(), converted.ptr);
}

void append_floating(std::string& output, double value) {
    std::array<char, 64> buffer{};
    const auto converted = std::to_chars(
        buffer.data(), buffer.data() + buffer.size(), value,
        std::chars_format::general, 6);
    if (converted.ec != std::errc{}) {
        fail("Could not format floating site statistic");
    }
    output.append(buffer.data(), converted.ptr);
}

struct FormatIndices {
    std::optional<std::size_t> gt;
    std::optional<std::size_t> dp;
    std::optional<std::size_t> gq;
};

FormatIndices format_indices(std::string_view format) {
    FormatIndices result;
    std::size_t index = 0;
    std::size_t begin = 0;
    while (begin <= format.size()) {
        const std::size_t end = format.find(':', begin);
        const std::string_view key = format.substr(
            begin,
            end == std::string_view::npos
                ? std::string_view::npos
                : end - begin);
        if (key == "GT") {
            result.gt = index;
        } else if (key == "DP") {
            result.dp = index;
        } else if (key == "GQ") {
            result.gq = index;
        }
        if (end == std::string_view::npos) {
            break;
        }
        begin = end + 1;
        ++index;
    }
    return result;
}

struct SampleFields {
    std::string_view gt;
    std::string_view dp;
    std::string_view gq;
};

SampleFields select_sample_fields(
    std::string_view sample,
    const FormatIndices& indices) {
    SampleFields result;
    const std::size_t required =
        static_cast<std::size_t>(indices.gt.has_value()) +
        static_cast<std::size_t>(indices.dp.has_value()) +
        static_cast<std::size_t>(indices.gq.has_value());
    if (required == 0) {
        return result;
    }
    if (indices.gt == 0 && !indices.dp.has_value() &&
        !indices.gq.has_value()) {
        if (sample.size() >= 4 && sample[3] == ':' &&
            (sample[1] == '/' || sample[1] == '|')) {
            result.gt = sample.substr(0, 3);
        } else {
            const std::size_t end = sample.find(':');
            result.gt = sample.substr(0, end);
        }
        return result;
    }
    std::size_t index = 0;
    std::size_t begin = 0;
    std::size_t found = 0;
    while (begin <= sample.size()) {
        const std::size_t end = sample.find(':', begin);
        const std::string_view value = sample.substr(
            begin,
            end == std::string_view::npos
                ? std::string_view::npos
                : end - begin);
        if (indices.gt.has_value() && index == *indices.gt) {
            result.gt = value;
            ++found;
        }
        if (indices.dp.has_value() && index == *indices.dp) {
            result.dp = value;
            ++found;
        }
        if (indices.gq.has_value() && index == *indices.gq) {
            result.gq = value;
            ++found;
        }
        if (found == required) {
            break;
        }
        if (end == std::string_view::npos) {
            break;
        }
        begin = end + 1;
        ++index;
    }
    return result;
}

void append_masked_genotype(
    std::string_view /* genotype */, std::string& output) {
    // VCFtools 0.1.17 recodes every filtered GT as an unphased
    // diploid missing genotype, including haploid and phased inputs.
    output.append("./.");
}

void append_recode_sample(
    std::string_view sample, std::size_t gt_index,
    bool filtered, std::string& output) {
    if (!filtered) {
        output.append(sample);
        return;
    }
    std::size_t index = 0;
    std::size_t begin = 0;
    while (begin <= sample.size()) {
        const std::size_t end = sample.find(':', begin);
        if (index > 0) {
            output.push_back(':');
        }
        const std::string_view value = sample.substr(
            begin,
            end == std::string_view::npos
                ? std::string_view::npos
                : end - begin);
        if (index == gt_index) {
            append_masked_genotype(value, output);
        } else {
            output.append(value);
        }
        if (end == std::string_view::npos) {
            return;
        }
        begin = end + 1;
        ++index;
    }
}

void append_sorted_filter_column(
    std::string_view filter, std::string& output) {
    if (filter == "." ||
        filter.find(';') == std::string_view::npos) {
        output.append(filter);
        return;
    }
    std::vector<std::string_view> names;
    std::size_t begin = 0;
    while (begin <= filter.size()) {
        const std::size_t end = filter.find(';', begin);
        names.push_back(filter.substr(
            begin,
            end == std::string_view::npos
                ? std::string_view::npos
                : end - begin));
        if (end == std::string_view::npos) {
            break;
        }
        begin = end + 1;
    }
    std::sort(names.begin(), names.end());
    for (std::size_t index = 0; index < names.size(); ++index) {
        if (index > 0) {
            output.push_back(';');
        }
        output.append(names[index]);
    }
}

void append_recode_record(
    std::string_view line,
    const std::array<std::string_view, 9>& columns,
    std::size_t samples_begin, std::size_t sample_count,
    const FormatIndices& indices,
    const std::vector<std::uint8_t>& genotype_filtered,
    const std::vector<std::size_t>& selected_samples,
    bool info_all, std::string& output) {
    const bool any_filtered = std::any_of(
        genotype_filtered.begin(), genotype_filtered.end(),
        [](std::uint8_t value) { return value != 0; });
    const bool filter_already_canonical =
        columns[6].find(';') == std::string_view::npos;
    if (info_all && !any_filtered && filter_already_canonical &&
        selected_samples.empty()) {
        output.append(line);
        output.push_back('\n');
        return;
    }

    for (std::size_t column = 0; column < 7; ++column) {
        if (column > 0) {
            output.push_back('\t');
        }
        if (column == 6) {
            append_sorted_filter_column(columns[column], output);
        } else {
            output.append(columns[column]);
        }
    }
    output.push_back('\t');
    output.append(info_all ? columns[7] : std::string_view("."));
    if (sample_count > 0) {
        output.push_back('\t');
        output.append(columns[8]);
    }

    std::size_t begin = samples_begin;
    std::size_t selected_cursor = 0;
    for (std::size_t sample_index = 0;
         sample_index < sample_count; ++sample_index) {
        const std::size_t end = line.find('\t', begin);
        const std::string_view sample = line.substr(
            begin,
            end == std::string_view::npos
                ? std::string_view::npos
                : end - begin);
        const bool selected =
            selected_samples.empty() ||
            (selected_cursor < selected_samples.size() &&
             selected_samples[selected_cursor] == sample_index);
        if (selected) {
            output.push_back('\t');
            append_recode_sample(
                sample, indices.gt.value_or(0),
                indices.gt.has_value() &&
                    genotype_filtered[sample_index] != 0,
                output);
            if (!selected_samples.empty()) {
                ++selected_cursor;
            }
        }
        begin =
            end == std::string_view::npos
                ? line.size()
                : end + 1;
    }
    output.push_back('\n');
}

struct GenotypeSummary {
    std::uint64_t called = 0;
    std::uint32_t ploidy = 0;
    std::uint32_t missing = 0;
    std::array<int, 2> alleles{-1, -1};
};

int parse_depth(std::string_view value) {
    if (value.empty() || value == ".") {
        return -1;
    }
    std::int64_t depth = 0;
    const auto parsed = std::from_chars(
        value.data(), value.data() + value.size(), depth);
    if (parsed.ec != std::errc{} ||
        parsed.ptr != value.data() + value.size() ||
        depth < std::numeric_limits<std::int32_t>::min() ||
        depth > std::numeric_limits<std::int32_t>::max()) {
        fail("Invalid DP value in VCF record");
    }
    return depth < 0 ? -1 : static_cast<int>(depth);
}

double parse_quality(std::string_view value) {
    if (value.empty() || value == ".") {
        return -1.0;
    }
    float quality = 0.0F;
    const auto parsed = std::from_chars(
        value.data(), value.data() + value.size(), quality,
        std::chars_format::general);
    if (parsed.ec != std::errc{} ||
        parsed.ptr != value.data() + value.size()) {
        fail("Invalid QUAL value in VCF record");
    }
    return static_cast<double>(quality);
}

template <bool CaptureAlleles, typename Counts>
GenotypeSummary count_genotype(
    std::string_view genotype, Counts& counts,
    std::size_t valid_alleles, bool include) {
    GenotypeSummary result;
    if (genotype.size() == 3 &&
        (genotype[1] == '/' || genotype[1] == '|') &&
        (genotype[0] == '.' ||
         (genotype[0] >= '0' && genotype[0] <= '9')) &&
        (genotype[2] == '.' ||
         (genotype[2] >= '0' && genotype[2] <= '9'))) {
        result.ploidy = 2;
        constexpr std::array<std::size_t, 2> offsets{0, 2};
        for (std::size_t copy = 0; copy < offsets.size(); ++copy) {
            const char encoded = genotype[offsets[copy]];
            if (encoded == '.') {
                ++result.missing;
                continue;
            }
            const std::size_t allele =
                static_cast<std::size_t>(encoded - '0');
            if (allele >= valid_alleles) {
                fail("Invalid GT allele in VCF record");
            }
            if (include) {
                ++counts[allele];
                ++result.called;
            }
            if constexpr (CaptureAlleles) {
                result.alleles[copy] = static_cast<int>(allele);
            }
        }
        return result;
    }
    std::size_t begin = 0;
    while (begin < genotype.size()) {
        std::size_t end = begin;
        while (end < genotype.size() &&
               genotype[end] != '/' && genotype[end] != '|') {
            ++end;
        }
        const std::string_view allele =
            genotype.substr(begin, end - begin);
        ++result.ploidy;
        if (result.ploidy > 2) {
            fail(
                "Polyploid genotype is not supported by "
                "VCFtools compatibility mode");
        }
        if (allele.empty() || allele == ".") {
            ++result.missing;
        } else {
            std::uint64_t index = 0;
            bool valid = false;
            if (allele.size() == 1 && allele.front() >= '0' &&
                allele.front() <= '9') {
                index = static_cast<std::uint64_t>(
                    allele.front() - '0');
                valid = true;
            } else {
                const auto parsed = std::from_chars(
                    allele.data(), allele.data() + allele.size(), index);
                valid = parsed.ec == std::errc{} &&
                        parsed.ptr == allele.data() + allele.size();
            }
            if (!valid || index >= valid_alleles) {
                fail("Invalid GT allele in VCF record");
            }
            if (include) {
                ++counts[static_cast<std::size_t>(index)];
                ++result.called;
            }
            if constexpr (CaptureAlleles) {
                result.alleles[result.ploidy - 1] =
                    static_cast<int>(index);
            }
        }
        begin = end + 1;
    }
    return result;
}

struct FastFstContribution {
    double sum_a = 0.0;
    double sum_all = 0.0;
    double fst = std::numeric_limits<double>::quiet_NaN();
};

FastFstContribution calculate_fast_fst(
    const std::vector<std::vector<std::uint32_t>>& homozygotes,
    const std::vector<std::vector<std::uint32_t>>& heterozygotes,
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
            n[population] += homozygotes[population][allele] +
                             0.5 * heterozygotes[population][allele];
            p[population][allele] =
                heterozygotes[population][allele] +
                2 * homozygotes[population][allele];
            pbar[allele] += p[population][allele];
            hbar[allele] += heterozygotes[population][allele];
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
    FastFstContribution result;
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
        if (!std::isnan(a) && !std::isnan(b) && !std::isnan(c)) {
            result.sum_a += a;
            result.sum_all += a + b + c;
        }
    }
    result.fst = result.sum_a / result.sum_all;
    return result;
}

FastFstContribution calculate_two_population_biallelic_fst(
    const std::array<std::array<std::uint32_t, 2>, 2>& homozygotes,
    const std::array<std::array<std::uint32_t, 2>, 2>& heterozygotes) {
    constexpr std::size_t population_count = 2;
    constexpr std::size_t allele_count = 2;
    std::array<double, population_count> n{};
    std::array<std::array<double, allele_count>, population_count> p{};
    std::array<double, allele_count> pbar{};
    std::array<double, allele_count> hbar{};
    std::array<double, allele_count> ssqr{};
    double sum_nsqr = 0.0;
    for (std::size_t population = 0;
         population < population_count; ++population) {
        for (std::size_t allele = 0; allele < allele_count; ++allele) {
            n[population] += homozygotes[population][allele] +
                             0.5 * heterozygotes[population][allele];
            p[population][allele] =
                heterozygotes[population][allele] +
                2 * homozygotes[population][allele];
            pbar[allele] += p[population][allele];
            hbar[allele] += heterozygotes[population][allele];
        }
        for (std::size_t allele = 0; allele < allele_count; ++allele) {
            p[population][allele] /= (2.0 * n[population]);
        }
        sum_nsqr += n[population] * n[population];
    }
    const double n_sum = std::accumulate(n.begin(), n.end(), 0.0);
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
        ssqr[allele] /= ((population_count - 1.0) * nbar);
    }
    const double nc =
        (n_sum - (sum_nsqr / n_sum)) /
        (population_count - 1.0);
    const double r = static_cast<double>(population_count);
    FastFstContribution result;
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
        if (!std::isnan(a) && !std::isnan(b) && !std::isnan(c)) {
            result.sum_a += a;
            result.sum_all += a + b + c;
        }
    }
    result.fst = result.sum_a / result.sum_all;
    return result;
}

bool append_site_stat_record(
    std::string_view line, std::size_t sample_count,
    const std::vector<std::size_t>& selected_samples,
    const FastSiteStatPlan& plan, PositionSelectionCursor& position_cursor,
    ShardOutput& output) {
    if (line.empty() || line.front() == '#') {
        return false;
    }
    if (line.back() == '\r') {
        line.remove_suffix(1);
    }

    std::array<std::string_view, 9> columns{};
    std::size_t begin = 0;
    for (std::size_t column = 0; column < columns.size(); ++column) {
        const std::size_t end = line.find('\t', begin);
        if (end == std::string_view::npos) {
            if (column < 7) {
                fail("VCF record has fewer than eight fixed columns");
            }
            columns[column] = line.substr(begin);
            begin = line.size();
            break;
        }
        columns[column] = line.substr(begin, end - begin);
        begin = end + 1;
    }
    const std::size_t samples_begin = begin;

    int position = 0;
    const auto parsed_position = std::from_chars(
        columns[1].data(),
        columns[1].data() + columns[1].size(), position);
    if (parsed_position.ec != std::errc{} ||
        parsed_position.ptr != columns[1].data() + columns[1].size()) {
        fail("Invalid POS value in VCF record");
    }
    if (plan.position_selection_active()) {
        update_position_cursor(position_cursor, columns[0], plan);
    }
    if (!plan.positions_file.empty() &&
        !contains_position(position_cursor.keep, position)) {
        return false;
    }
    if (!plan.exclude_positions_file.empty() &&
        contains_position(position_cursor.exclude, position)) {
        return false;
    }

    std::array<std::string_view, kInlineAlleles> inline_alleles{};
    std::size_t allele_count = 1;
    inline_alleles[0] = columns[3];
    std::vector<std::string_view> heap_alleles;
    std::string_view alternate = columns[4];
    if (!alternate.empty() && alternate != ".") {
        while (true) {
            const std::size_t separator = alternate.find(',');
            const std::string_view allele =
                alternate.substr(0, separator);
            if (allele_count < inline_alleles.size()) {
                inline_alleles[allele_count] = allele;
            } else {
                if (heap_alleles.empty()) {
                    heap_alleles.assign(
                        inline_alleles.begin(),
                        inline_alleles.end());
                }
                heap_alleles.push_back(allele);
            }
            ++allele_count;
            if (separator == std::string_view::npos) {
                break;
            }
            alternate.remove_prefix(separator + 1);
        }
    }
    if (static_cast<int>(allele_count) < plan.min_alleles ||
        static_cast<int>(allele_count) > plan.max_alleles) {
        return false;
    }
    const double quality = parse_quality(columns[5]);
    if (plan.min_quality >= 0.0 &&
        quality < plan.min_quality) {
        return false;
    }

    std::array<std::uint64_t, kInlineAlleles> inline_counts{};
    std::vector<std::uint64_t> heap_counts;
    if (allele_count > inline_counts.size()) {
        heap_counts.assign(allele_count, 0);
    }
    const FormatIndices indices = format_indices(columns[8]);
    std::uint64_t called = 0;
    std::uint32_t n_data = 0;
    std::uint32_t n_missing = 0;
    std::uint32_t sum_depth = 0;
    std::uint32_t sumsq_depth = 0;
    std::uint32_t depth_count = 0;
    double raw_depth_sum = 0.0;
    std::uint64_t total_ploidy = 0;
    const bool need_gt =
        plan.freq || plan.freq2 ||
        plan.counts || plan.missing_site ||
        plan.min_call_rate > 0.0 || plan.min_maf > 0.0 ||
        plan.advanced_statistics_active();
    const bool need_dp =
        plan.site_depth || plan.site_mean_depth ||
        plan.min_mean_depth > 0.0;
    if (plan.recode) {
        output.genotype_filtered.assign(sample_count, 0);
    }
    const bool subset_samples = !selected_samples.empty();
    const std::size_t analysis_sample_count =
        subset_samples ? selected_samples.size() : sample_count;
    std::size_t selected_cursor = 0;
    bool fully_diploid = true;
    std::vector<std::vector<std::uint32_t>> population_homozygotes;
    std::vector<std::vector<std::uint32_t>> population_heterozygotes;
    const bool scalar_fst =
        plan.fst_population_files.size() == 2 && allele_count == 2;
    std::array<std::array<std::uint32_t, 2>, 2>
        scalar_population_homozygotes{};
    std::array<std::array<std::uint32_t, 2>, 2>
        scalar_population_heterozygotes{};
    if (!plan.fst_population_files.empty() && !scalar_fst) {
        population_homozygotes.assign(
            plan.fst_population_files.size(),
            std::vector<std::uint32_t>(allele_count, 0));
        population_heterozygotes.assign(
            plan.fst_population_files.size(),
            std::vector<std::uint32_t>(allele_count, 0));
    }
    for (std::size_t sample_index = 0;
         sample_index < sample_count; ++sample_index) {
        std::string_view sample;
        if (begin < line.size()) {
            const std::size_t end = line.find('\t', begin);
            sample = line.substr(
                begin,
                end == std::string_view::npos
                    ? std::string_view::npos
                    : end - begin);
            begin =
                end == std::string_view::npos
                    ? line.size()
                    : end + 1;
        }
        if (subset_samples &&
            (selected_cursor >= selected_samples.size() ||
             selected_samples[selected_cursor] != sample_index)) {
            continue;
        }
        if (subset_samples) {
            ++selected_cursor;
        }
        const SampleFields fields =
            select_sample_fields(sample, indices);
        const bool genotype_filtered =
            plan.min_genotype_quality > 0.0 &&
            indices.gq.has_value() &&
            parse_depth(fields.gq) < plan.min_genotype_quality;
        if (plan.recode) {
            output.genotype_filtered[sample_index] =
                genotype_filtered;
        }
        if (need_gt && indices.gt.has_value()) {
            GenotypeSummary genotype;
            if (heap_counts.empty()) {
                genotype = !plan.fst_population_files.empty()
                               ? count_genotype<true>(
                                     fields.gt, inline_counts,
                                     allele_count, !genotype_filtered)
                               : count_genotype<false>(
                                     fields.gt, inline_counts,
                                     allele_count, !genotype_filtered);
            } else {
                genotype = !plan.fst_population_files.empty()
                               ? count_genotype<true>(
                                     fields.gt, heap_counts,
                                     allele_count, !genotype_filtered)
                               : count_genotype<false>(
                                     fields.gt, heap_counts,
                                     allele_count, !genotype_filtered);
            }
            total_ploidy += genotype.ploidy;
            called += genotype.called;
            if (!genotype_filtered) {
                if (plan.advanced_statistics_active() &&
                    genotype.ploidy != 2) {
                    fully_diploid = false;
                }
                n_data += genotype.ploidy;
                n_missing += genotype.missing;
                if (genotype.ploidy >= 2 &&
                    genotype.alleles[0] >= 0 &&
                    genotype.alleles[1] >= 0 &&
                    sample_index < plan.population_memberships.size()) {
                    for (const std::size_t population :
                         plan.population_memberships[sample_index]) {
                        const auto first = static_cast<std::size_t>(
                            genotype.alleles[0]);
                        const auto second = static_cast<std::size_t>(
                            genotype.alleles[1]);
                        if (first == second) {
                            if (scalar_fst) {
                                ++scalar_population_homozygotes
                                      [population][first];
                            } else {
                                ++population_homozygotes
                                      [population][first];
                            }
                        } else {
                            if (scalar_fst) {
                                ++scalar_population_heterozygotes
                                      [population][first];
                                ++scalar_population_heterozygotes
                                      [population][second];
                            } else {
                                ++population_heterozygotes
                                      [population][first];
                                ++population_heterozygotes
                                      [population][second];
                            }
                        }
                    }
                }
            }
        }
        if (need_dp &&
            indices.dp.has_value()) {
            const int depth = parse_depth(fields.dp);
            if (depth >= 0) {
                raw_depth_sum += depth;
            }
            if (!genotype_filtered && depth >= 0) {
                const auto unsigned_depth =
                    static_cast<std::uint32_t>(depth);
                sum_depth += unsigned_depth;
                sumsq_depth += unsigned_depth * unsigned_depth;
                ++depth_count;
            }
        }
    }
    if (!indices.gt.has_value() && plan.missing_site) {
        const std::uint64_t absent =
            static_cast<std::uint64_t>(analysis_sample_count) * 2;
        if (absent > std::numeric_limits<std::uint32_t>::max()) {
            fail("Too many samples for site missingness output");
        }
        n_data = static_cast<std::uint32_t>(absent);
        n_missing = n_data;
    }

    if (plan.min_mean_depth > 0.0 &&
        raw_depth_sum / static_cast<double>(analysis_sample_count) <
            plan.min_mean_depth) {
        return false;
    }
    if (plan.min_call_rate > 0.0 &&
        called / static_cast<double>(total_ploidy) <
            plan.min_call_rate) {
        return false;
    }
    if (plan.min_maf > 0.0) {
        double maf = std::numeric_limits<double>::max();
        if (called > 0) {
            for (std::size_t allele = 0;
                 allele < allele_count; ++allele) {
                const double frequency =
                    (heap_counts.empty()
                         ? inline_counts[allele]
                         : heap_counts[allele]) /
                    static_cast<double>(called);
                maf = std::min(
                    maf, std::min(frequency, 1.0 - frequency));
            }
        }
        if (maf < plan.min_maf) {
            return false;
        }
    }

    const auto allele_at =
        [&](std::size_t index) -> std::string_view {
        return heap_alleles.empty()
                   ? inline_alleles[index]
                   : heap_alleles[index];
    };
    const auto allele_count_at =
        [&](std::size_t index) -> std::uint64_t {
        return heap_counts.empty()
                   ? inline_counts[index]
                   : heap_counts[index];
    };
    const auto append_position =
        [&](std::string& text) {
        text.append(columns[0]);
        text.push_back('\t');
        text.append(columns[1]);
    };

    if (plan.freq || plan.freq2) {
        std::string& text =
            output.text[artifact_index(Artifact::freq)];
        append_position(text);
        text.push_back('\t');
        append_unsigned(text, allele_count);
        text.push_back('\t');
        append_unsigned(text, called);
        for (std::size_t allele = 0;
             allele < allele_count; ++allele) {
            text.push_back('\t');
            if (plan.freq) {
                text.append(allele_at(allele));
                text.push_back(':');
            }
            append_floating(
                text,
                allele_count_at(allele) /
                    static_cast<double>(called));
        }
        text.push_back('\n');
    }
    if (plan.counts) {
        std::string& text =
            output.text[artifact_index(Artifact::counts)];
        append_position(text);
        text.push_back('\t');
        append_unsigned(text, allele_count);
        text.push_back('\t');
        append_unsigned(text, called);
        for (std::size_t allele = 0;
             allele < allele_count; ++allele) {
            text.push_back('\t');
            text.append(allele_at(allele));
            text.push_back(':');
            append_unsigned(text, allele_count_at(allele));
        }
        text.push_back('\n');
    }
    if (plan.missing_site) {
        std::string& text =
            output.text[artifact_index(Artifact::missing)];
        append_position(text);
        text.push_back('\t');
        append_unsigned(text, n_data);
        text.append("\t0\t");
        append_unsigned(text, n_missing);
        text.push_back('\t');
        append_floating(
            text,
            n_missing / static_cast<double>(n_data));
        text.push_back('\n');
    }
    if (plan.site_depth) {
        std::string& text =
            output.text[artifact_index(Artifact::depth)];
        append_position(text);
        text.push_back('\t');
        append_unsigned(text, sum_depth);
        text.push_back('\t');
        append_unsigned(text, sumsq_depth);
        text.push_back('\n');
    }
    if (plan.site_mean_depth) {
        std::string& text =
            output.text[artifact_index(Artifact::mean_depth)];
        const double mean =
            sum_depth / static_cast<double>(depth_count);
        const double variance =
            ((sumsq_depth /
                  static_cast<double>(depth_count)) -
             (mean * mean)) *
            depth_count /
            static_cast<double>(depth_count - 1);
        append_position(text);
        text.push_back('\t');
        append_floating(text, mean);
        text.push_back('\t');
        append_floating(text, variance);
        text.push_back('\n');
    }
    if (plan.site_quality) {
        std::string& text =
            output.text[artifact_index(Artifact::quality)];
        append_position(text);
        text.push_back('\t');
        append_floating(text, quality);
        text.push_back('\n');
    }
    if (plan.advanced_statistics_active()) {
        ShardOutput::AdvancedSite site;
        site.chromosome.assign(columns[0]);
        site.position = position;
        site.allele_count = static_cast<std::uint32_t>(allele_count);
        if (allele_count > site.inline_allele_counts.size()) {
            site.heap_allele_counts.resize(allele_count);
        }
        for (std::size_t allele = 0; allele < allele_count; ++allele) {
            if (site.heap_allele_counts.empty()) {
                site.inline_allele_counts[allele] =
                    static_cast<std::uint32_t>(allele_count_at(allele));
            } else {
                site.heap_allele_counts[allele] =
                    static_cast<std::uint32_t>(allele_count_at(allele));
            }
        }
        site.non_missing_chromosomes =
            static_cast<std::uint32_t>(called);
        site.fully_diploid = fully_diploid;
        if (!plan.fst_population_files.empty() && fully_diploid) {
            const auto fst =
                scalar_fst
                    ? calculate_two_population_biallelic_fst(
                          scalar_population_homozygotes,
                          scalar_population_heterozygotes)
                    : calculate_fast_fst(
                          population_homozygotes,
                          population_heterozygotes, allele_count);
            site.fst_eligible = true;
            site.fst_sum_a = fst.sum_a;
            site.fst_sum_all = fst.sum_all;
            site.fst = fst.fst;
        }
        output.advanced_sites.push_back(std::move(site));
    }
    if (plan.recode) {
        append_recode_record(
            line, columns, samples_begin, sample_count,
            indices, output.genotype_filtered,
            selected_samples,
            plan.recode_info_all, output.recode_text);
    }
    return true;
}

ShardOutput process_plain_shard(
    int descriptor, const char* mapped_input, const PlainShard& shard,
    std::string& input_bytes, std::size_t samples,
    const std::vector<std::size_t>& selected_samples,
    const FastSiteStatPlan& plan) {
    const std::uint64_t length = shard.end - shard.begin;
    if (length >
        static_cast<std::uint64_t>(
            std::numeric_limits<std::size_t>::max())) {
        fail("Plain VCF shard is too large for this platform");
    }
    std::string_view shard_bytes;
    if (mapped_input != nullptr) {
        shard_bytes = std::string_view(
            mapped_input + static_cast<std::size_t>(shard.begin),
            static_cast<std::size_t>(length));
    } else {
        input_bytes.resize(static_cast<std::size_t>(length));
        std::size_t read_total = 0;
        while (read_total < input_bytes.size()) {
            const ssize_t count = ::pread(
                descriptor, input_bytes.data() + read_total,
                input_bytes.size() - read_total,
                static_cast<off_t>(shard.begin + read_total));
            if (count < 0) {
                if (errno == EINTR) {
                    continue;
                }
                const std::string error = std::strerror(errno);
                fail("Could not read plain VCF shard: " + error);
            }
            if (count == 0) {
                fail("Plain VCF ended inside an aligned shard");
            }
            read_total += static_cast<std::size_t>(count);
        }
        shard_bytes = input_bytes;
    }

    ShardOutput result;
    const std::size_t reserve = shard_bytes.size() / 48;
    if (plan.freq || plan.freq2) {
        result.text[artifact_index(Artifact::freq)].reserve(reserve);
    }
    if (plan.counts) {
        result.text[artifact_index(Artifact::counts)].reserve(reserve);
    }
    if (plan.missing_site) {
        result.text[artifact_index(Artifact::missing)].reserve(reserve);
    }
    if (plan.site_depth) {
        result.text[artifact_index(Artifact::depth)].reserve(reserve / 2);
    }
    if (plan.site_mean_depth) {
        result.text[artifact_index(Artifact::mean_depth)].reserve(reserve);
    }
    if (plan.site_quality) {
        result.text[artifact_index(Artifact::quality)].reserve(reserve / 2);
    }
    if (plan.recode) {
        result.recode_text.reserve(shard_bytes.size());
    }
    std::size_t begin = 0;
    PositionSelectionCursor position_cursor;
    while (begin < shard_bytes.size()) {
        std::size_t end = shard_bytes.find('\n', begin);
        if (end == std::string::npos) {
            end = shard_bytes.size();
        }
        if (end > begin) {
            result.kept += append_site_stat_record(
                shard_bytes.substr(begin, end - begin),
                samples, selected_samples, plan, position_cursor,
                result);
            ++result.records;
        }
        begin = end == shard_bytes.size() ? end : end + 1;
    }
    return result;
}

FastSiteStatsSummary run_plain_site_stats(
    const std::string& input_path,
    const std::string& output_prefix,
    unsigned requested_threads,
    const FastSiteStatPlan& plan) {
    const std::uint64_t file_size =
        std::filesystem::file_size(input_path);
    const HeaderLayout header =
        read_plain_header(input_path, file_size);
    FastSiteStatPlan effective_plan = plan;
    if (!effective_plan.fst_population_files.empty()) {
        effective_plan.population_memberships.assign(
            header.samples, {});
        for (std::size_t population = 0;
             population < effective_plan.fst_population_files.size();
             ++population) {
            std::set<std::string> members;
            load_sample_files(
                {effective_plan.fst_population_files[population]},
                members);
            for (std::size_t sample = 0;
                 sample < header.sample_names.size(); ++sample) {
                if (members.contains(header.sample_names[sample])) {
                    effective_plan.population_memberships[sample].push_back(
                        population);
                }
            }
        }
    }
    const auto selected_samples =
        selected_sample_indices(effective_plan, header);
    if (effective_plan.recode) {
        HtsFilePtr header_input(hts_open(input_path.c_str(), "r"));
        if (!header_input) {
            fail("Could not open plain VCF header: " + input_path);
        }
        const std::string raw_header =
            read_text_vcf_header(header_input.get(), input_path);
        effective_plan.recode_sink(subset_text_vcf_header(
            raw_header, selected_samples));
    }
    const auto shards = build_plain_shards(
        input_path, header.data_start, file_size,
        requested_threads);
    const bool prefer_mapping =
        requested_threads > 1 ||
        plan.position_selection_active() ||
        plan.sample_selection_active();
    const ReadOnlyMapping mapping(
        input_path, file_size, prefer_mapping);
    const unsigned worker_count = descriptor_limited_workers(
        requested_threads, shards.size());

    OrderedArtifactSet outputs(
        output_prefix, effective_plan,
        selected_samples.empty() ? header.samples
                                 : selected_samples.size());

    std::vector<std::optional<ShardOutput>> completed(shards.size());
    std::mutex mutex;
    std::condition_variable available;
    std::condition_variable capacity_available;
    std::size_t next_commit = 0;
    const std::size_t maximum_ahead =
        std::max<std::size_t>(1, worker_count);
    std::atomic<bool> cancelled{false};
    std::exception_ptr error;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (unsigned worker = 0; worker < worker_count; ++worker) {
        workers.emplace_back([&, worker] {
            try {
                FileDescriptor descriptor(
                    mapping.data == nullptr
                        ? ::open(input_path.c_str(), O_RDONLY)
                        : -1);
                if (mapping.data == nullptr && descriptor.value < 0) {
                    fail(
                        "Could not open plain VCF worker: " +
                        std::string(std::strerror(errno)));
                }
                std::string input_bytes;
                for (std::size_t index = worker;
                     index < shards.size();
                     index += worker_count) {
                    if (cancelled.load(std::memory_order_relaxed)) {
                        return;
                    }
                    auto result = process_plain_shard(
                        descriptor.value, mapping.data, shards[index],
                        input_bytes, header.samples,
                        selected_samples, effective_plan);
                    {
                        std::unique_lock lock(mutex);
                        capacity_available.wait(lock, [&] {
                            return cancelled.load(
                                       std::memory_order_relaxed) ||
                                   index <
                                       next_commit + maximum_ahead;
                        });
                        if (cancelled.load(
                                std::memory_order_relaxed)) {
                            return;
                        }
                        completed[index] = std::move(result);
                    }
                    available.notify_all();
                }
            } catch (...) {
                {
                    std::lock_guard lock(mutex);
                    if (!error) {
                        error = std::current_exception();
                    }
                    cancelled.store(true, std::memory_order_relaxed);
                }
                available.notify_all();
                capacity_available.notify_all();
            }
        });
    }

    std::uint64_t records = 0;
    std::uint64_t kept = 0;
    try {
        for (std::size_t index = 0; index < shards.size(); ++index) {
            std::optional<ShardOutput> result;
            {
                std::unique_lock lock(mutex);
                available.wait(lock, [&] {
                    return completed[index].has_value() ||
                           cancelled.load(std::memory_order_relaxed);
                });
                if (cancelled.load(std::memory_order_relaxed)) {
                    break;
                }
                result = std::move(completed[index]);
                completed[index].reset();
                next_commit = index + 1;
            }
            capacity_available.notify_all();
            outputs.append(*result);
            records += result->records;
            kept += result->kept;
        }
    } catch (...) {
        std::lock_guard lock(mutex);
        if (!error) {
            error = std::current_exception();
        }
    }
    cancelled.store(true, std::memory_order_relaxed);
    available.notify_all();
    capacity_available.notify_all();
    for (auto& worker : workers) {
        worker.join();
    }
    if (error) {
        std::rethrow_exception(error);
    }
    outputs.validate();
    return FastSiteStatsSummary{
        .total = records,
        .kept = kept,
        .samples = selected_samples.empty()
                       ? header.samples
                       : selected_samples.size(),
        .input_threads = worker_count,
        .hts_io_threads = 0,
        .planned_shards = shards.size(),
        .backend =
            effective_plan.recode
                ? "fast-filter-recode-plain"
                : plan.counts_only()
                ? "fast-counts-plain"
                : "fast-site-stats-plain",
        .description =
            effective_plan.recode
                ? mapping.data != nullptr
                      ? "zero-copy mapped aligned ranges with direct filtering/recode"
                      : "fused aligned byte ranges with direct filtering/recode"
                : mapping.data != nullptr
                ? "zero-copy mapped aligned ranges with direct site statistics"
                : "fused aligned byte ranges with direct site statistics",
    };
}

FastSiteStatsSummary run_compressed_site_stats(
    HtsFilePtr input, const std::string& input_path,
    const std::string& output_prefix, unsigned requested_threads,
    bool bgzf, const FastSiteStatPlan& plan) {
    const unsigned io_threads =
        bgzf && requested_threads > 1
            ? std::min(3u, requested_threads - 1)
            : 0;
    if (io_threads > 0 &&
        hts_set_threads(input.get(), static_cast<int>(io_threads)) != 0) {
        fail("Could not enable BGZF decompression threads: " + input_path);
    }
    OrderedArtifactSet outputs(output_prefix, plan, 0);

    KStringBuffer line;
    std::size_t samples = 0;
    bool saw_header = false;
    std::uint64_t records = 0;
    std::uint64_t kept = 0;
    ShardOutput buffer;
    PositionSelectionCursor position_cursor;
    const std::vector<std::size_t> all_samples;
    std::size_t buffered_bytes = 0;
    while (hts_getline(
               input.get(), '\n', &line.value) >= 0) {
        const std::string_view text(line.value.s, line.value.l);
        if (!saw_header) {
            if (text.rfind("#CHROM\t", 0) == 0) {
                const std::size_t columns = count_tabs(text) + 1;
                samples = columns > 9 ? columns - 9 : 0;
                saw_header = true;
            } else if (text.empty() || text.front() != '#') {
                fail(
                    "VCF data appeared before #CHROM header: " +
                    input_path);
            }
            continue;
        }
        kept += append_site_stat_record(
            text, samples, all_samples, plan, position_cursor, buffer);
        ++records;
        buffered_bytes = 0;
        for (const auto& artifact : buffer.text) {
            buffered_bytes += artifact.size();
        }
        buffered_bytes += buffer.recode_text.size();
        if (buffered_bytes >= kOutputFlushBytes) {
            outputs.append(buffer);
            for (auto& artifact : buffer.text) {
                artifact.clear();
            }
            buffer.recode_text.clear();
        }
    }
    if (!saw_header) {
        fail("Could not find #CHROM header: " + input_path);
    }
    outputs.append(buffer);
    outputs.validate();
    return FastSiteStatsSummary{
        .total = records,
        .kept = kept,
        .samples = samples,
        .input_threads = 1,
        .hts_io_threads = io_threads,
        .planned_shards = 1,
        .backend =
            plan.recode
                ? (bgzf
                       ? "fast-filter-recode-bgzf"
                       : "fast-filter-recode-gzip")
                : plan.counts_only()
                ? (bgzf
                       ? "fast-counts-bgzf"
                       : "fast-counts-gzip")
                : (bgzf
                       ? "fast-site-stats-bgzf"
                       : "fast-site-stats-gzip"),
        .description =
            plan.recode
                ? (bgzf
                       ? "direct filtering/recode with BGZF decompression overlap"
                       : "direct filtering/recode from compressed VCF")
                : bgzf
                ? "direct site statistics with BGZF decompression overlap"
                : "direct site statistics from compressed VCF",
    };
}

hts_pos_t record_position(std::string_view line) {
    const std::size_t first = line.find('\t');
    const std::size_t second =
        first == std::string_view::npos
            ? std::string_view::npos
            : line.find('\t', first + 1);
    if (first == std::string_view::npos ||
        second == std::string_view::npos) {
        fail("VCF record has no POS column");
    }
    std::uint64_t position = 0;
    const std::string_view value =
        line.substr(first + 1, second - first - 1);
    const auto parsed = std::from_chars(
        value.data(), value.data() + value.size(), position);
    if (parsed.ec != std::errc{} ||
        parsed.ptr != value.data() + value.size() ||
        position == 0 ||
        position - 1 >
            static_cast<std::uint64_t>(HTS_POS_MAX)) {
        fail("Invalid VCF POS value");
    }
    return static_cast<hts_pos_t>(position - 1);
}

hts_pos_t header_contig_length(
    bcf_hdr_t* header, const char* name) {
    const int rid = bcf_hdr_name2id(header, name);
    if (rid < 0 || rid >= header->n[BCF_DT_CTG] ||
        header->id[BCF_DT_CTG][rid].val == nullptr) {
        return HTS_POS_MAX;
    }
    const hts_pos_t length = static_cast<hts_pos_t>(
        header->id[BCF_DT_CTG][rid].val->info[0]);
    return length > 0 ? length : HTS_POS_MAX;
}

std::vector<IndexedShard> build_indexed_shards(
    tbx_t* index, bcf_hdr_t* header, unsigned threads) {
    int count = 0;
    const char** names = tbx_seqnames(index, &count);
    if (count > 0 && names == nullptr) {
        fail("Could not read indexed VCF contig names");
    }
    struct Contig {
        int tid = -1;
        hts_pos_t length = HTS_POS_MAX;
        std::uint64_t records = 0;
    };
    std::vector<Contig> contigs;
    contigs.reserve(std::max(0, count));
    std::uint64_t total_records = 0;
    for (int tid = 0; tid < count; ++tid) {
        if (bcf_hdr_name2id(header, names[tid]) < 0) {
            std::free(const_cast<char**>(names));
            fail("Indexed contig is absent from VCF header");
        }
        std::uint64_t mapped = 0;
        std::uint64_t unmapped = 0;
        if (hts_idx_get_stat(
                index->idx, tid, &mapped, &unmapped) != 0) {
            std::free(const_cast<char**>(names));
            fail("Could not read VCF index statistics");
        }
        contigs.push_back(Contig{
            .tid = tid,
            .length = header_contig_length(header, names[tid]),
            .records = mapped,
        });
        total_records += mapped;
    }
    std::free(const_cast<char**>(names));

    const std::size_t target_shards = std::min<std::size_t>(
        kMaximumShards,
        std::max<std::size_t>(
            contigs.size(),
            static_cast<std::size_t>(std::max(1u, threads)) * 8));
    std::vector<IndexedShard> shards;
    shards.reserve(target_shards + contigs.size());
    for (const auto& contig : contigs) {
        if (contig.length == HTS_POS_MAX ||
            contig.length <= 1 ||
            total_records == 0 ||
            contig.records == 0) {
            shards.push_back(IndexedShard{
                .tid = contig.tid,
                .begin = 0,
                .end = HTS_POS_MAX,
            });
            continue;
        }
        const std::size_t proportional =
            static_cast<std::size_t>(std::ceil(
                (static_cast<long double>(target_shards) *
                 static_cast<long double>(contig.records)) /
                static_cast<long double>(total_records)));
        const std::size_t windows = std::max<std::size_t>(
            1,
            std::min<std::size_t>(
                proportional,
                static_cast<std::size_t>(contig.length)));
        const hts_pos_t width = std::max<hts_pos_t>(
            1,
            (contig.length + static_cast<hts_pos_t>(windows) - 1) /
                static_cast<hts_pos_t>(windows));
        for (hts_pos_t begin = 0; begin < contig.length;) {
            const hts_pos_t end =
                std::min(contig.length, begin + width);
            shards.push_back(IndexedShard{
                .tid = contig.tid,
                .begin = begin,
                .end = end,
            });
            begin = end;
        }
    }
    return shards;
}

unsigned descriptor_limited_workers(
    unsigned requested, std::size_t shards) {
    unsigned limit = std::max(1u, requested);
    struct rlimit descriptors {};
    if (getrlimit(RLIMIT_NOFILE, &descriptors) == 0 &&
        descriptors.rlim_cur != RLIM_INFINITY) {
        constexpr rlim_t reserve = 64;
        const rlim_t available =
            descriptors.rlim_cur > reserve
                ? descriptors.rlim_cur - reserve
                : 1;
        limit = std::min<unsigned>(
            limit,
            static_cast<unsigned>(std::min<rlim_t>(
                available,
                std::numeric_limits<unsigned>::max())));
    }
    return std::min<unsigned>(
        limit,
        std::max<std::size_t>(1, shards));
}

FastSiteStatsSummary run_indexed_site_stats(
    const std::string& input_path,
    const std::string& index_path,
    const std::string& output_prefix,
    unsigned requested_threads,
    const FastSiteStatPlan& plan) {
    HtsFilePtr probe(hts_open(input_path.c_str(), "r"));
    if (!probe) {
        fail("Could not open indexed VCF: " + input_path);
    }
    HeaderPtr header(bcf_hdr_read(probe.get()));
    if (!header) {
        fail("Could not read indexed VCF header: " + input_path);
    }
    const std::size_t samples =
        static_cast<std::size_t>(bcf_hdr_nsamples(header.get()));
    probe.reset();

    TbxPtr index(tbx_index_load2(
        input_path.c_str(), index_path.c_str()));
    if (!index) {
        fail("Could not load selected VCF index: " + index_path);
    }
    const auto shards = build_indexed_shards(
        index.get(), header.get(), requested_threads);
    const unsigned worker_count = descriptor_limited_workers(
        requested_threads, shards.size());
    const std::vector<std::size_t> all_samples;

    OrderedArtifactSet outputs(output_prefix, plan, samples);

    std::vector<std::optional<ShardOutput>> completed(shards.size());
    std::mutex mutex;
    std::condition_variable available;
    std::condition_variable capacity_available;
    std::size_t next_commit = 0;
    const std::size_t maximum_ahead =
        std::max<std::size_t>(1, worker_count);
    std::atomic<std::size_t> next_shard{0};
    std::atomic<bool> cancelled{false};
    std::exception_ptr error;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (unsigned worker = 0; worker < worker_count; ++worker) {
        workers.emplace_back([&] {
            try {
                HtsFilePtr input(hts_open(input_path.c_str(), "r"));
                if (!input) {
                    fail("Could not open indexed VCF worker");
                }
                KStringBuffer line;
                while (!cancelled.load(std::memory_order_relaxed)) {
                    const std::size_t ordinal =
                        next_shard.fetch_add(
                            1, std::memory_order_relaxed);
                    if (ordinal >= shards.size()) {
                        return;
                    }
                    const auto& shard = shards[ordinal];
                    IteratorPtr iterator(tbx_itr_queryi(
                        index.get(), shard.tid,
                        shard.begin, shard.end));
                    if (!iterator) {
                        fail("Could not create tabix region iterator");
                    }
                    ShardOutput result;
                    PositionSelectionCursor position_cursor;
                    while (tbx_itr_next(
                               input.get(), index.get(),
                               iterator.get(), &line.value) >= 0) {
                        const std::string_view text(
                            line.value.s, line.value.l);
                        const hts_pos_t position =
                            record_position(text);
                        if (position < shard.begin ||
                            (shard.end != HTS_POS_MAX &&
                             position >= shard.end)) {
                            continue;
                        }
                        result.kept += append_site_stat_record(
                            text, samples, all_samples, plan,
                            position_cursor, result);
                        ++result.records;
                    }
                    {
                        std::unique_lock lock(mutex);
                        capacity_available.wait(lock, [&] {
                            return cancelled.load(
                                       std::memory_order_relaxed) ||
                                   ordinal <
                                       next_commit + maximum_ahead;
                        });
                        if (cancelled.load(
                                std::memory_order_relaxed)) {
                            return;
                        }
                        completed[ordinal] = std::move(result);
                    }
                    available.notify_all();
                }
            } catch (...) {
                {
                    std::lock_guard lock(mutex);
                    if (!error) {
                        error = std::current_exception();
                    }
                    cancelled.store(true, std::memory_order_relaxed);
                }
                available.notify_all();
                capacity_available.notify_all();
            }
        });
    }

    std::uint64_t records = 0;
    std::uint64_t kept = 0;
    try {
        for (std::size_t ordinal = 0;
             ordinal < shards.size(); ++ordinal) {
            std::optional<ShardOutput> result;
            {
                std::unique_lock lock(mutex);
                available.wait(lock, [&] {
                    return completed[ordinal].has_value() ||
                           cancelled.load(std::memory_order_relaxed);
                });
                if (cancelled.load(std::memory_order_relaxed)) {
                    break;
                }
                result = std::move(completed[ordinal]);
                completed[ordinal].reset();
                next_commit = ordinal + 1;
            }
            capacity_available.notify_all();
            outputs.append(*result);
            records += result->records;
            kept += result->kept;
        }
    } catch (...) {
        std::lock_guard lock(mutex);
        if (!error) {
            error = std::current_exception();
        }
    }
    cancelled.store(true, std::memory_order_relaxed);
    available.notify_all();
    capacity_available.notify_all();
    for (auto& worker : workers) {
        worker.join();
    }
    if (error) {
        std::rethrow_exception(error);
    }
    outputs.validate();
    return FastSiteStatsSummary{
        .total = records,
        .kept = kept,
        .samples = samples,
        .input_threads = worker_count,
        .hts_io_threads = 0,
        .planned_shards = shards.size(),
        .backend =
            plan.recode
                ? "fast-filter-recode-indexed-bgzf"
                : plan.counts_only()
                ? "fast-counts-indexed-bgzf"
                : "fast-site-stats-indexed-bgzf",
        .description =
            (plan.recode
                 ? "ordered tabix regions with direct filtering/recode via "
                 : "ordered tabix regions with direct site statistics via ") +
            index_path,
    };
}

}  // namespace

std::optional<FastSiteStatsSummary> run_fast_text_site_stats(
    const std::string& output_prefix,
    const input::SourceOptions& options,
    const FastSiteStatPlan& plan) {
    FastSiteStatPlan effective_plan = plan;
    load_position_file(
        effective_plan.positions_file,
        effective_plan.positions_to_keep);
    load_position_file(
        effective_plan.exclude_positions_file,
        effective_plan.positions_to_exclude);
    const std::string& input_path = options.path;
    const unsigned threads = std::min(
        std::max(1u, options.total_threads),
        input::detect_available_threads().count);
    HtsFilePtr probe(hts_open(input_path.c_str(), "r"));
    if (!probe) {
        fail("Could not open input file: " + input_path);
    }
    const htsFormat* format = hts_get_format(probe.get());
    if (format == nullptr || format->format != vcf) {
        return std::nullopt;
    }
    const auto compression = format->compression;
    if (compression != no_compression &&
        (effective_plan.position_selection_active() ||
         effective_plan.sample_selection_active())) {
        return std::nullopt;
    }
    if (compression == no_compression) {
        probe.reset();
        return run_plain_site_stats(
            input_path, output_prefix, threads, effective_plan);
    }
    const bool is_bgzf = compression == bgzf;
    if (effective_plan.recode) {
        const std::string header =
            read_text_vcf_header(probe.get(), input_path);
        effective_plan.recode_sink(header);
    }
    if (is_bgzf) {
        probe.reset();
        input::SourceOptions effective_options = options;
        effective_options.total_threads = threads;
        const std::string index_path =
            input::prepare_variant_index(effective_options);
        if (!index_path.empty()) {
            return run_indexed_site_stats(
                input_path, index_path, output_prefix,
                threads, effective_plan);
        }
        probe.reset(hts_open(input_path.c_str(), "r"));
        if (!probe) {
            fail("Could not reopen compressed VCF: " + input_path);
        }
    }
    if (!is_bgzf && effective_plan.recode) {
        probe.reset(hts_open(input_path.c_str(), "r"));
        if (!probe) {
            fail("Could not reopen compressed VCF: " + input_path);
        }
    }
    return run_compressed_site_stats(
        std::move(probe), input_path, output_prefix,
        threads, is_bgzf, effective_plan);
}

}  // namespace vcftools_ng
