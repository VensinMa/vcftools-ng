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
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include <fcntl.h>
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
    std::uint64_t records = 0;
    std::uint64_t kept = 0;
};

class OrderedArtifactSet {
public:
    OrderedArtifactSet(
        const std::string& prefix,
        const FastSiteStatPlan& plan)
        : plan_(plan) {
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
    }

    void validate() {
        for (std::size_t index = 0; index < outputs_.size(); ++index) {
            if (enabled_[index]) {
                vcftools_ng::output::finish_stream(
                    outputs_[index], final_paths_[index]);
            }
        }
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
    const FastSiteStatPlan& plan_;
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
            const std::size_t columns = count_tabs(line) + 1;
            layout.samples = columns > 9 ? columns - 9 : 0;
            return layout;
        }
        if (line.empty() || line.front() != '#') {
            fail("VCF data appeared before #CHROM header: " + path);
        }
    }
    fail("Could not find #CHROM header: " + path);
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
    std::size_t index = 0;
    std::size_t begin = 0;
    while (begin <= sample.size()) {
        const std::size_t end = sample.find(':', begin);
        const std::string_view value = sample.substr(
            begin,
            end == std::string_view::npos
                ? std::string_view::npos
                : end - begin);
        if (indices.gt.has_value() && index == *indices.gt) {
            result.gt = value;
        }
        if (indices.dp.has_value() && index == *indices.dp) {
            result.dp = value;
        }
        if (indices.gq.has_value() && index == *indices.gq) {
            result.gq = value;
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
    bool info_all, std::string& output) {
    const bool any_filtered = std::any_of(
        genotype_filtered.begin(), genotype_filtered.end(),
        [](std::uint8_t value) { return value != 0; });
    const bool filter_already_canonical =
        columns[6].find(';') == std::string_view::npos;
    if (info_all && !any_filtered && filter_already_canonical) {
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
    for (std::size_t sample_index = 0;
         sample_index < sample_count; ++sample_index) {
        const std::size_t end = line.find('\t', begin);
        const std::string_view sample = line.substr(
            begin,
            end == std::string_view::npos
                ? std::string_view::npos
                : end - begin);
        output.push_back('\t');
        append_recode_sample(
            sample, indices.gt.value_or(0),
            indices.gt.has_value() &&
                genotype_filtered[sample_index] != 0,
            output);
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

template <typename Counts>
GenotypeSummary count_genotype(
    std::string_view genotype, Counts& counts,
    std::size_t valid_alleles, bool include) {
    GenotypeSummary result;
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
            const auto parsed = std::from_chars(
                allele.data(), allele.data() + allele.size(), index);
            if (parsed.ec != std::errc{} ||
                parsed.ptr != allele.data() + allele.size() ||
                index >= valid_alleles) {
                fail("Invalid GT allele in VCF record");
            }
            if (include) {
                ++counts[static_cast<std::size_t>(index)];
                ++result.called;
            }
        }
        begin = end + 1;
    }
    return result;
}

bool append_site_stat_record(
    std::string_view line, std::size_t sample_count,
    const FastSiteStatPlan& plan, ShardOutput& output) {
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
        plan.min_call_rate > 0.0 || plan.min_maf > 0.0;
    const bool need_dp =
        plan.site_depth || plan.site_mean_depth ||
        plan.min_mean_depth > 0.0;
    if (plan.recode) {
        output.genotype_filtered.assign(sample_count, 0);
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
            const GenotypeSummary genotype =
                heap_counts.empty()
                    ? count_genotype(
                          fields.gt, inline_counts, allele_count,
                          !genotype_filtered)
                    : count_genotype(
                          fields.gt, heap_counts, allele_count,
                          !genotype_filtered);
            total_ploidy += genotype.ploidy;
            called += genotype.called;
            if (!genotype_filtered) {
                n_data += genotype.ploidy;
                n_missing += genotype.missing;
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
            static_cast<std::uint64_t>(sample_count) * 2;
        if (absent > std::numeric_limits<std::uint32_t>::max()) {
            fail("Too many samples for site missingness output");
        }
        n_data = static_cast<std::uint32_t>(absent);
        n_missing = n_data;
    }

    if (plan.min_mean_depth > 0.0 &&
        raw_depth_sum / static_cast<double>(sample_count) <
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
    if (plan.recode) {
        append_recode_record(
            line, columns, samples_begin, sample_count,
            indices, output.genotype_filtered,
            plan.recode_info_all, output.recode_text);
    }
    return true;
}

ShardOutput process_plain_shard(
    int descriptor, const PlainShard& shard,
    std::string& input_bytes, std::size_t samples,
    const FastSiteStatPlan& plan) {
    const std::uint64_t length = shard.end - shard.begin;
    if (length >
        static_cast<std::uint64_t>(
            std::numeric_limits<std::size_t>::max())) {
        fail("Plain VCF shard is too large for this platform");
    }
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

    ShardOutput result;
    const std::size_t reserve = input_bytes.size() / 48;
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
        result.recode_text.reserve(input_bytes.size());
    }
    std::size_t begin = 0;
    while (begin < input_bytes.size()) {
        std::size_t end = input_bytes.find('\n', begin);
        if (end == std::string::npos) {
            end = input_bytes.size();
        }
        if (end > begin) {
            result.kept += append_site_stat_record(
                std::string_view(
                    input_bytes.data() + begin, end - begin),
                samples, plan, result);
            ++result.records;
        }
        begin = end == input_bytes.size() ? end : end + 1;
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
    const auto shards = build_plain_shards(
        input_path, header.data_start, file_size,
        requested_threads);
    const unsigned worker_count = descriptor_limited_workers(
        requested_threads, shards.size());

    OrderedArtifactSet outputs(output_prefix, plan);

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
                    ::open(input_path.c_str(), O_RDONLY));
                if (descriptor.value < 0) {
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
                        descriptor.value, shards[index],
                        input_bytes, header.samples, plan);
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
        .samples = header.samples,
        .input_threads = worker_count,
        .hts_io_threads = 0,
        .planned_shards = shards.size(),
        .backend =
            plan.recode
                ? "fast-filter-recode-plain"
                : plan.counts_only()
                ? "fast-counts-plain"
                : "fast-site-stats-plain",
        .description =
            plan.recode
                ? "fused aligned byte ranges with direct filtering/recode"
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
    OrderedArtifactSet outputs(output_prefix, plan);

    KStringBuffer line;
    std::size_t samples = 0;
    bool saw_header = false;
    std::uint64_t records = 0;
    std::uint64_t kept = 0;
    ShardOutput buffer;
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
            text, samples, plan, buffer);
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

    OrderedArtifactSet outputs(output_prefix, plan);

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
                            text, samples, plan, result);
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
    if (plan.recode) {
        const std::string header =
            read_text_vcf_header(probe.get(), input_path);
        plan.recode_sink(header);
    }
    if (compression == no_compression) {
        probe.reset();
        return run_plain_site_stats(
            input_path, output_prefix, threads, plan);
    }
    const bool is_bgzf = compression == bgzf;
    if (is_bgzf) {
        probe.reset();
        input::SourceOptions effective_options = options;
        effective_options.total_threads = threads;
        const std::string index_path =
            input::prepare_variant_index(effective_options);
        if (!index_path.empty()) {
            return run_indexed_site_stats(
                input_path, index_path, output_prefix,
                threads, plan);
        }
        probe.reset(hts_open(input_path.c_str(), "r"));
        if (!probe) {
            fail("Could not reopen compressed VCF: " + input_path);
        }
    }
    if (!is_bgzf && plan.recode) {
        probe.reset(hts_open(input_path.c_str(), "r"));
        if (!probe) {
            fail("Could not reopen compressed VCF: " + input_path);
        }
    }
    return run_compressed_site_stats(
        std::move(probe), input_path, output_prefix,
        threads, is_bgzf, plan);
}

}  // namespace vcftools_ng
