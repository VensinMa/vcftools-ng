#include "fast_counts.h"

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
constexpr std::string_view kCountsHeader =
    "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\n";

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
    std::string text;
    std::uint64_t records = 0;
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

std::string_view subfield(
    std::string_view sample, std::size_t index) {
    std::size_t begin = 0;
    for (std::size_t current = 0; current < index; ++current) {
        const std::size_t separator = sample.find(':', begin);
        if (separator == std::string_view::npos) {
            return {};
        }
        begin = separator + 1;
    }
    const std::size_t end = sample.find(':', begin);
    return sample.substr(
        begin,
        end == std::string_view::npos
            ? std::string_view::npos
            : end - begin);
}

std::optional<std::size_t> genotype_index(
    std::string_view format) {
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
            return index;
        }
        if (end == std::string_view::npos) {
            break;
        }
        begin = end + 1;
        ++index;
    }
    return std::nullopt;
}

template <typename Counts>
std::uint64_t count_genotype(
    std::string_view genotype, Counts& counts) {
    std::uint64_t called = 0;
    std::size_t ploidy = 0;
    std::size_t begin = 0;
    while (begin < genotype.size()) {
        std::size_t end = begin;
        while (end < genotype.size() &&
               genotype[end] != '/' && genotype[end] != '|') {
            ++end;
        }
        const std::string_view allele =
            genotype.substr(begin, end - begin);
        ++ploidy;
        if (ploidy > 2) {
            fail(
                "Polyploid genotype is not supported by "
                "VCFtools compatibility mode");
        }
        if (!allele.empty() && allele != ".") {
            std::uint64_t index = 0;
            const auto parsed = std::from_chars(
                allele.data(), allele.data() + allele.size(), index);
            if (parsed.ec != std::errc{} ||
                parsed.ptr != allele.data() + allele.size() ||
                index >= counts.size()) {
                fail("Invalid GT allele in VCF record");
            }
            ++counts[static_cast<std::size_t>(index)];
            ++called;
        }
        begin = end + 1;
    }
    return called;
}

void append_count_record(
    std::string_view line, std::string& output) {
    if (line.empty() || line.front() == '#') {
        return;
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

    std::array<std::uint64_t, kInlineAlleles> inline_counts{};
    std::vector<std::uint64_t> heap_counts;
    if (allele_count > inline_counts.size()) {
        heap_counts.assign(allele_count, 0);
    }
    const auto gt_index = genotype_index(columns[8]);
    std::uint64_t called = 0;
    while (begin < line.size() && gt_index.has_value()) {
        const std::size_t end = line.find('\t', begin);
        const std::string_view sample = line.substr(
            begin,
            end == std::string_view::npos
                ? std::string_view::npos
                : end - begin);
        const std::string_view genotype =
            subfield(sample, *gt_index);
        called += heap_counts.empty()
                      ? count_genotype(genotype, inline_counts)
                      : count_genotype(genotype, heap_counts);
        if (end == std::string_view::npos) {
            break;
        }
        begin = end + 1;
    }

    output.append(columns[0]);
    output.push_back('\t');
    output.append(columns[1]);
    output.push_back('\t');
    append_unsigned(output, allele_count);
    output.push_back('\t');
    append_unsigned(output, called);
    for (std::size_t allele = 0; allele < allele_count; ++allele) {
        output.push_back('\t');
        output.append(
            heap_alleles.empty()
                ? inline_alleles[allele]
                : heap_alleles[allele]);
        output.push_back(':');
        append_unsigned(
            output,
            heap_counts.empty()
                ? inline_counts[allele]
                : heap_counts[allele]);
    }
    output.push_back('\n');
}

ShardOutput process_plain_shard(
    int descriptor, const PlainShard& shard,
    std::string& input_bytes) {
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
    result.text.reserve(input_bytes.size() / 32);
    std::size_t begin = 0;
    while (begin < input_bytes.size()) {
        std::size_t end = input_bytes.find('\n', begin);
        if (end == std::string::npos) {
            end = input_bytes.size();
        }
        if (end > begin) {
            append_count_record(
                std::string_view(
                    input_bytes.data() + begin, end - begin),
                result.text);
            ++result.records;
        }
        begin = end == input_bytes.size() ? end : end + 1;
    }
    return result;
}

FastCountsSummary run_plain_counts(
    const std::string& input_path,
    const std::string& output_prefix,
    unsigned requested_threads) {
    const std::uint64_t file_size =
        std::filesystem::file_size(input_path);
    const HeaderLayout header =
        read_plain_header(input_path, file_size);
    const auto shards = build_plain_shards(
        input_path, header.data_start, file_size,
        requested_threads);
    const unsigned worker_count = descriptor_limited_workers(
        requested_threads, shards.size());

    std::ofstream output(
        output_prefix + ".frq.count",
        std::ios::binary | std::ios::trunc);
    if (!output) {
        fail("Could not open counts output");
    }
    output.write(kCountsHeader.data(), kCountsHeader.size());

    std::vector<std::optional<ShardOutput>> completed(shards.size());
    std::mutex mutex;
    std::condition_variable available;
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
                        input_bytes);
                    {
                        std::lock_guard lock(mutex);
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
            }
        });
    }

    std::uint64_t records = 0;
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
        }
        output.write(result->text.data(), result->text.size());
        records += result->records;
    }
    cancelled.store(true, std::memory_order_relaxed);
    for (auto& worker : workers) {
        worker.join();
    }
    if (error) {
        std::rethrow_exception(error);
    }
    if (!output) {
        fail("Could not write counts output");
    }
    return FastCountsSummary{
        .total = records,
        .kept = records,
        .samples = header.samples,
        .input_threads = worker_count,
        .hts_io_threads = 0,
        .planned_shards = shards.size(),
        .backend = "fast-counts-plain",
        .description =
            "fused aligned byte ranges with direct GT counting",
    };
}

FastCountsSummary run_compressed_counts(
    HtsFilePtr input, const std::string& input_path,
    const std::string& output_prefix, unsigned requested_threads,
    bool bgzf) {
    const unsigned io_threads =
        bgzf && requested_threads > 1
            ? std::min(3u, requested_threads - 1)
            : 0;
    if (io_threads > 0 &&
        hts_set_threads(input.get(), static_cast<int>(io_threads)) != 0) {
        fail("Could not enable BGZF decompression threads: " + input_path);
    }
    std::ofstream output(
        output_prefix + ".frq.count",
        std::ios::binary | std::ios::trunc);
    if (!output) {
        fail("Could not open counts output");
    }
    output.write(kCountsHeader.data(), kCountsHeader.size());

    KStringBuffer line;
    std::size_t samples = 0;
    bool saw_header = false;
    std::uint64_t records = 0;
    std::string buffer;
    buffer.reserve(kOutputFlushBytes + 4096);
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
        append_count_record(text, buffer);
        ++records;
        if (buffer.size() >= kOutputFlushBytes) {
            output.write(buffer.data(), buffer.size());
            buffer.clear();
        }
    }
    if (!saw_header) {
        fail("Could not find #CHROM header: " + input_path);
    }
    output.write(buffer.data(), buffer.size());
    if (!output) {
        fail("Could not write counts output");
    }
    return FastCountsSummary{
        .total = records,
        .kept = records,
        .samples = samples,
        .input_threads = 1,
        .hts_io_threads = io_threads,
        .planned_shards = 1,
        .backend = bgzf
                       ? "fast-counts-bgzf"
                       : "fast-counts-gzip",
        .description =
            bgzf
                ? "direct GT counting with BGZF decompression overlap"
                : "direct GT counting from compressed VCF",
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

FastCountsSummary run_indexed_counts(
    const std::string& input_path,
    const std::string& index_path,
    const std::string& output_prefix,
    unsigned requested_threads) {
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

    std::ofstream output(
        output_prefix + ".frq.count",
        std::ios::binary | std::ios::trunc);
    if (!output) {
        fail("Could not open counts output");
    }
    output.write(kCountsHeader.data(), kCountsHeader.size());

    std::vector<std::optional<ShardOutput>> completed(shards.size());
    std::mutex mutex;
    std::condition_variable available;
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
                        append_count_record(text, result.text);
                        ++result.records;
                    }
                    {
                        std::lock_guard lock(mutex);
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
            }
        });
    }

    std::uint64_t records = 0;
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
        }
        output.write(result->text.data(), result->text.size());
        records += result->records;
    }
    cancelled.store(true, std::memory_order_relaxed);
    for (auto& worker : workers) {
        worker.join();
    }
    if (error) {
        std::rethrow_exception(error);
    }
    if (!output) {
        fail("Could not write counts output");
    }
    return FastCountsSummary{
        .total = records,
        .kept = records,
        .samples = samples,
        .input_threads = worker_count,
        .hts_io_threads = 0,
        .planned_shards = shards.size(),
        .backend = "fast-counts-indexed-bgzf",
        .description =
            "ordered tabix regions with direct GT counting via " +
            index_path,
    };
}

}  // namespace

std::optional<FastCountsSummary> run_fast_text_counts(
    const std::string& output_prefix,
    const input::SourceOptions& options) {
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
    if (format->compression == no_compression) {
        probe.reset();
        return run_plain_counts(
            input_path, output_prefix, threads);
    }
    const bool is_bgzf = format->compression == bgzf;
    if (is_bgzf && threads >= 3) {
        probe.reset();
        input::SourceOptions effective_options = options;
        effective_options.total_threads = threads;
        if (threads <= 4) {
            effective_options.auto_index = false;
        }
        const std::string index_path =
            input::prepare_variant_index(effective_options);
        if (!index_path.empty()) {
            return run_indexed_counts(
                input_path, index_path, output_prefix, threads);
        }
        probe.reset(hts_open(input_path.c_str(), "r"));
        if (!probe) {
            fail("Could not reopen compressed VCF: " + input_path);
        }
    }
    return run_compressed_counts(
        std::move(probe), input_path, output_prefix,
        threads, is_bgzf);
}

}  // namespace vcftools_ng
