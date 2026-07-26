#include "input_source.h"

#include <htslib/hts.h>
#include <htslib/tbx.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>
#include <utility>

#include <fcntl.h>
#include <sched.h>
#include <sys/file.h>
#include <sys/resource.h>
#include <sys/stat.h>
#if defined(__linux__)
#include <sys/sysmacros.h>
#endif
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace vcftools_ng::input {
namespace {

constexpr std::uint64_t kMib = 1024ULL * 1024ULL;
constexpr std::uint64_t kMinimumPlainShardBytes = 4ULL * kMib;
constexpr std::uint64_t kTargetPlainShardBytes = 64ULL * kMib;
constexpr std::size_t kMaximumShards = 65536;

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

std::optional<bool> detect_rotational_storage(
    const std::string& path) {
#if defined(__linux__)
    struct stat status {};
    if (::stat(path.c_str(), &status) != 0) {
        return std::nullopt;
    }
    const std::filesystem::path device_link =
        std::filesystem::path("/sys/dev/block") /
        (std::to_string(major(status.st_dev)) + ":" +
         std::to_string(minor(status.st_dev)));
    std::error_code error;
    std::filesystem::path device =
        std::filesystem::canonical(device_link, error);
    if (error) {
        return std::nullopt;
    }
    while (!device.empty() &&
           device != device.root_path()) {
        std::ifstream rotational(
            device / "queue" / "rotational");
        int value = -1;
        if (rotational >> value) {
            return value != 0;
        }
        device = device.parent_path();
    }
#else
    (void)path;
#endif
    return std::nullopt;
}

struct HtsFileDeleter {
    void operator()(htsFile* file) const noexcept {
        if (file != nullptr) {
            hts_close(file);
        }
    }
};

struct HeaderDeleter {
    void operator()(bcf_hdr_t* header) const noexcept {
        if (header != nullptr) {
            bcf_hdr_destroy(header);
        }
    }
};

struct IndexDeleter {
    void operator()(hts_idx_t* index) const noexcept {
        if (index != nullptr) {
            hts_idx_destroy(index);
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

struct TbxDeleter {
    void operator()(tbx_t* index) const noexcept {
        if (index != nullptr) {
            tbx_destroy(index);
        }
    }
};

using HtsFilePtr = std::unique_ptr<htsFile, HtsFileDeleter>;
using HeaderPtr = std::unique_ptr<bcf_hdr_t, HeaderDeleter>;
using IndexPtr = std::unique_ptr<hts_idx_t, IndexDeleter>;
using IteratorPtr = std::unique_ptr<hts_itr_t, IteratorDeleter>;
using TbxPtr = std::unique_ptr<tbx_t, TbxDeleter>;

HeaderPtr read_header(htsFile* input, const std::string& path) {
    HeaderPtr header(bcf_hdr_read(input));
    if (!header) {
        fail("Could not read VCF/BCF header: " + path);
    }
    return header;
}

HeaderPtr duplicate_header(bcf_hdr_t* source) {
    HeaderPtr duplicate(bcf_hdr_dup(source));
    if (!duplicate) {
        fail("Could not duplicate VCF/BCF header");
    }
    return duplicate;
}

HtsFilePtr open_input(const std::string& path) {
    HtsFilePtr input(hts_open(path.c_str(), "r"));
    if (!input) {
        fail("Could not open input file: " + path);
    }
    return input;
}

bool path_entry_exists(const std::string& path) {
    std::error_code error;
    const auto status = std::filesystem::symlink_status(path, error);
    return !error &&
           status.type() != std::filesystem::file_type::not_found;
}

struct IndexValidation {
    bool valid = false;
    std::string reason;
};

struct KStringBuffer {
    kstring_t value{0, 0, nullptr};

    ~KStringBuffer() {
        std::free(value.s);
    }
};

struct NameArrayDeleter {
    void operator()(const char** names) const noexcept {
        std::free(const_cast<char**>(names));
    }
};

IndexValidation validate_index_file(
    const std::string& data_path,
    const std::string& index_path, bool is_bcf,
    bcf_hdr_t* header) {
    if (!std::filesystem::is_regular_file(index_path)) {
        return {
            false,
            "not a readable regular file"};
    }

    std::error_code data_time_error;
    std::error_code index_time_error;
    const auto data_time =
        std::filesystem::last_write_time(
            data_path, data_time_error);
    const auto index_time =
        std::filesystem::last_write_time(
            index_path, index_time_error);
    if (!data_time_error && !index_time_error &&
        index_time < data_time) {
        return {
            false,
            "index is older than the data file"};
    }

    HtsFilePtr probe_input = open_input(data_path);
    if (is_bcf) {
        IndexPtr index(bcf_index_load2(
            data_path.c_str(), index_path.c_str()));
        if (!index) {
            return {
                false,
                "HTSlib could not load it as a BCF CSI index"};
        }
        const int indexed_sequences = hts_idx_nseq(index.get());
        if (indexed_sequences < 0 ||
            indexed_sequences > header->n[BCF_DT_CTG]) {
            return {
                false,
                "indexed sequence count is incompatible with "
                "the BCF header"};
        }
        for (int tid = 0; tid < indexed_sequences; ++tid) {
            std::uint64_t mapped = 0;
            std::uint64_t unmapped = 0;
            if (hts_idx_get_stat(
                    index.get(), tid, &mapped, &unmapped) != 0) {
                return {
                    false,
                    "could not read BCF index statistics"};
            }
            if (mapped == 0) {
                continue;
            }
            IteratorPtr iterator(bcf_itr_queryi(
                index.get(), tid, 0, HTS_POS_MAX));
            RecordPtr record(bcf_init());
            if (!iterator || !record ||
                bcf_itr_next(
                    probe_input.get(), iterator.get(),
                    record.get()) < 0 ||
                record->rid != tid) {
                return {
                    false,
                    "BCF index probe did not return the expected "
                    "contig"};
            }
        }
        return {true, {}};
    }

    TbxPtr index(tbx_index_load2(
        data_path.c_str(), index_path.c_str()));
    if (!index) {
        return {
            false,
            "HTSlib could not load it as a VCF TBI/CSI index"};
    }
    if ((index->conf.preset & 0xffff) != TBX_VCF) {
        return {
            false,
            "tabix preset is not VCF"};
    }
    int indexed_sequences = 0;
    std::unique_ptr<const char*, NameArrayDeleter> names(
        tbx_seqnames(index.get(), &indexed_sequences));
    if (indexed_sequences < 0 ||
        (indexed_sequences > 0 && !names)) {
        return {
            false,
            "could not read indexed VCF contig names"};
    }
    HeaderPtr probe_header = duplicate_header(header);
    KStringBuffer line;
    for (int tid = 0; tid < indexed_sequences; ++tid) {
        const char* name = names.get()[tid];
        const int expected_rid =
            name == nullptr
                ? -1
                : bcf_hdr_name2id(header, name);
        if (expected_rid < 0) {
            return {
                false,
                "indexed contig is absent from the VCF header"};
        }
        std::uint64_t mapped = 0;
        std::uint64_t unmapped = 0;
        if (hts_idx_get_stat(
                index->idx, tid, &mapped, &unmapped) != 0) {
            return {
                false,
                "could not read VCF index statistics"};
        }
        if (mapped == 0) {
            continue;
        }
        IteratorPtr iterator(tbx_itr_queryi(
            index.get(), tid, 0, HTS_POS_MAX));
        RecordPtr record(bcf_init());
        line.value.l = 0;
        if (!iterator || !record ||
            tbx_itr_next(
                probe_input.get(), index.get(), iterator.get(),
                &line.value) < 0 ||
            vcf_parse1(
                &line.value, probe_header.get(),
                record.get()) != 0 ||
            record->rid != expected_rid) {
            return {
                false,
                "VCF index probe did not return the expected "
                "contig"};
        }
    }
    return {true, {}};
}

struct SidecarInspection {
    bool any_present = false;
    std::string selected_path;
    std::vector<std::string> invalid;
};

SidecarInspection inspect_sidecars(
    const std::string& data_path, bool is_bcf,
    bcf_hdr_t* header) {
    SidecarInspection result;
    for (const std::string extension : {".csi", ".tbi"}) {
        const std::string candidate = data_path + extension;
        if (!path_entry_exists(candidate)) {
            continue;
        }
        result.any_present = true;
        const IndexValidation validation =
            validate_index_file(
                data_path, candidate, is_bcf, header);
        if (validation.valid && result.selected_path.empty()) {
            result.selected_path = candidate;
        } else if (!validation.valid) {
            result.invalid.push_back(
                candidate + ": " + validation.reason);
        }
    }
    return result;
}

std::string invalid_sidecar_summary(
    const SidecarInspection& inspection) {
    std::string result;
    for (const auto& diagnostic : inspection.invalid) {
        if (!result.empty()) {
            result += "; ";
        }
        result += diagnostic;
    }
    return result;
}

class FileDescriptor {
public:
    explicit FileDescriptor(int descriptor = -1)
        : descriptor_(descriptor) {}

    ~FileDescriptor() {
        if (descriptor_ >= 0) {
            close(descriptor_);
        }
    }

    FileDescriptor(const FileDescriptor&) = delete;
    FileDescriptor& operator=(const FileDescriptor&) = delete;

    int get() const noexcept {
        return descriptor_;
    }

private:
    int descriptor_;
};

std::uint64_t available_memory_bytes() {
    std::ifstream memory("/proc/meminfo");
    std::string key;
    std::uint64_t value_kib = 0;
    std::string unit;
    while (memory >> key >> value_kib >> unit) {
        if (key == "MemAvailable:") {
            return value_kib * 1024ULL;
        }
    }
    return 0;
}

bool can_prefetch_to_page_cache(
    const std::string& path, unsigned threads) {
    if (threads < 4 ||
        !std::filesystem::is_regular_file(path)) {
        return false;
    }
    constexpr std::uint64_t maximum_prefetch =
        64ULL * 1024ULL * 1024ULL * 1024ULL;
    const std::uint64_t available = available_memory_bytes();
    if (available == 0) {
        return false;
    }
    const std::uint64_t bytes =
        std::filesystem::file_size(path);
    return bytes <= maximum_prefetch &&
           bytes <= available / 2;
}

bool prefetch_to_page_cache(const std::string& path) {
    FileDescriptor input(::open(path.c_str(), O_RDONLY));
    if (input.get() < 0) {
        std::cerr
            << "Storage prefetch warning: could not open "
            << path << ": " << std::strerror(errno) << "\n";
        return false;
    }
#if defined(POSIX_FADV_SEQUENTIAL)
    (void)::posix_fadvise(
        input.get(), 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    constexpr std::size_t buffer_bytes =
        16ULL * 1024ULL * 1024ULL;
    std::vector<char> buffer(buffer_bytes);
    const auto started = std::chrono::steady_clock::now();
    std::uint64_t total = 0;
    while (true) {
        const ssize_t count =
            ::read(input.get(), buffer.data(), buffer.size());
        if (count == 0) {
            break;
        }
        if (count < 0) {
            if (errno == EINTR) {
                continue;
            }
            std::cerr
                << "Storage prefetch warning: read failed: "
                << std::strerror(errno) << "\n";
            return false;
        }
        total += static_cast<std::uint64_t>(count);
    }
    const double seconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started)
            .count();
    std::cerr
        << "Storage prefetch: cached " << total
        << " bytes from rotational storage in "
        << seconds << " seconds\n";
    return true;
}

struct IndexBuildResult {
    bool success = false;
    std::string detail;
    std::string index_path;

    IndexBuildResult(
        bool succeeded, std::string message,
        std::string selected_index = {})
        : success(succeeded),
          detail(std::move(message)),
          index_path(std::move(selected_index)) {}
};

IndexBuildResult build_csi_index(
    const SourceOptions& options, bool is_bcf,
    bcf_hdr_t* header) {
    const int input_descriptor =
        open(options.path.c_str(), O_RDONLY | O_CLOEXEC);
    if (input_descriptor < 0) {
        return {
            false,
            "could not open input for automatic indexing: " +
                std::string(std::strerror(errno))};
    }
    FileDescriptor input_lock(input_descriptor);
    if (flock(input_lock.get(), LOCK_EX) != 0) {
        return {
            false,
            "could not lock input for automatic indexing: " +
                std::string(std::strerror(errno))};
    }

    SidecarInspection current =
        inspect_sidecars(options.path, is_bcf, header);
    if (!current.selected_path.empty()) {
        return {
            true,
            "a valid index was created by another process while waiting",
            current.selected_path};
    }
    if (current.any_present) {
        return {
            false,
            "an unusable sidecar appeared while waiting; refusing "
            "to overwrite it: " +
                invalid_sidecar_summary(current)};
    }

    const std::string final_path = options.path + ".csi";
    const std::string temporary_path =
        final_path + ".vcftools-ng.tmp." +
        std::to_string(static_cast<long long>(getpid())) + "." +
        std::to_string(
            std::chrono::steady_clock::now()
                .time_since_epoch()
                .count());
    const std::string thread_count =
        std::to_string(std::max(1u, options.total_threads));

    std::cerr
        << "Auto-index: no CSI/TBI sidecar found; running "
        << options.bcftools_path << " index --csi --threads "
        << thread_count << " for " << options.path << "\n";

    const pid_t child = fork();
    if (child < 0) {
        return {
            false,
            "could not start bcftools: " +
                std::string(std::strerror(errno))};
    }
    if (child == 0) {
        close(input_lock.get());
        execlp(
            options.bcftools_path.c_str(),
            options.bcftools_path.c_str(),
            "index", "--csi", "--threads",
            thread_count.c_str(), "--output",
            temporary_path.c_str(), options.path.c_str(),
            static_cast<char*>(nullptr));
        dprintf(
            STDERR_FILENO,
            "vcftools-ng: could not execute %s: %s\n",
            options.bcftools_path.c_str(), std::strerror(errno));
        _exit(127);
    }

    int status = 0;
    pid_t waited = -1;
    do {
        waited = waitpid(child, &status, 0);
    } while (waited < 0 && errno == EINTR);
    if (waited < 0) {
        std::error_code ignored;
        std::filesystem::remove(temporary_path, ignored);
        return {
            false,
            "could not wait for bcftools: " +
                std::string(std::strerror(errno))};
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        std::error_code ignored;
        std::filesystem::remove(temporary_path, ignored);
        const std::string outcome =
            WIFSIGNALED(status)
                ? "terminated by signal " +
                      std::to_string(WTERMSIG(status))
                : "exited with status " +
                      std::to_string(WEXITSTATUS(status));
        return {
            false,
            "bcftools index " + outcome};
    }

    const IndexValidation temporary_validation =
        validate_index_file(
            options.path, temporary_path, is_bcf, header);
    if (!temporary_validation.valid) {
        std::error_code ignored;
        std::filesystem::remove(temporary_path, ignored);
        return {
            false,
            "bcftools produced an unusable CSI: " +
                temporary_validation.reason};
    }

    current = inspect_sidecars(options.path, is_bcf, header);
    if (!current.selected_path.empty()) {
        std::error_code ignored;
        std::filesystem::remove(temporary_path, ignored);
        return {
            true,
            "a valid index was created by another process",
            current.selected_path};
    }
    if (current.any_present) {
        std::error_code ignored;
        std::filesystem::remove(temporary_path, ignored);
        return {
            false,
            "an unusable sidecar appeared during indexing; "
            "refusing to overwrite it: " +
                invalid_sidecar_summary(current)};
    }

    if (link(temporary_path.c_str(), final_path.c_str()) != 0) {
        const int publish_errno = errno;
        std::error_code ignored;
        std::filesystem::remove(temporary_path, ignored);
        if (publish_errno == EEXIST) {
            current =
                inspect_sidecars(options.path, is_bcf, header);
            if (!current.selected_path.empty()) {
                return {
                    true,
                    "a valid index was created by another process",
                    current.selected_path};
            }
        }
        return {
            false,
            "could not publish CSI index " + final_path +
                ": " + std::string(std::strerror(publish_errno))};
    }
    std::error_code cleanup_error;
    std::filesystem::remove(temporary_path, cleanup_error);
    if (cleanup_error) {
        std::cerr
            << "Auto-index warning: CSI was published but temporary "
            << "link could not be removed: "
            << cleanup_error.message() << "\n";
    }
    if (!std::filesystem::is_regular_file(final_path)) {
        return {
            false,
            "bcftools completed but CSI index was not created"};
    }

    std::cerr
        << "Auto-index: created " << final_path << "\n";
    return {
        true,
        "created CSI with " + thread_count +
            " requested bcftools threads",
        final_path};
}

unsigned descriptor_limited_input_threads(unsigned requested);

ResourcePlan plan_resources(
    unsigned requested_threads, bool parallel_input,
    bool compressed_stream, bool text_parsing = false,
    std::optional<bool> rotational = std::nullopt,
    bool page_cache_prefetched = false) {
    ResourcePlan plan;
    plan.total_threads = std::max(1u, requested_threads);
    plan.storage_profile_known = rotational.has_value();
    plan.rotational_storage = rotational.value_or(false);
    plan.page_cache_prefetched = page_cache_prefetched;
    if (parallel_input) {
        if (plan.rotational_storage &&
            !plan.page_cache_prefetched) {
            plan.input_threads = 1;
            plan.compute_threads = plan.total_threads;
        } else {
            plan.input_threads = std::min(
                plan.total_threads,
                descriptor_limited_input_threads(
                    plan.total_threads));
            plan.compute_threads =
                text_parsing
                    ? std::max(
                          1u, (plan.total_threads + 1u) / 2u)
                    : std::max(
                          1u, (plan.total_threads + 2u) / 3u);
        }
    } else {
        plan.compute_threads = plan.total_threads;
    }
    if (!parallel_input && compressed_stream &&
        plan.total_threads > 2) {
        plan.hts_io_threads = std::min(
            8u, std::max(1u, plan.total_threads / 4u));
        plan.compute_threads = std::max(
            1u, plan.total_threads - plan.hts_io_threads);
    }
    return plan;
}

std::string storage_note(const ResourcePlan& plan) {
    if (!plan.storage_profile_known) {
        return "storage profile unknown";
    }
    if (plan.rotational_storage) {
        return plan.page_cache_prefetched
                   ? "rotational storage prefetched into page cache"
                   : "rotational storage, low-seek input";
    }
    return "non-rotational storage";
}

unsigned descriptor_limited_input_threads(unsigned requested) {
    struct rlimit descriptor_limit {};
    if (::getrlimit(RLIMIT_NOFILE, &descriptor_limit) != 0 ||
        descriptor_limit.rlim_cur == RLIM_INFINITY) {
        return requested;
    }
    constexpr rlim_t reserved_descriptors = 64;
    constexpr rlim_t descriptors_per_input_worker = 2;
    if (descriptor_limit.rlim_cur <= reserved_descriptors) {
        return 1;
    }
    const rlim_t available_workers =
        (descriptor_limit.rlim_cur - reserved_descriptors) /
        descriptors_per_input_worker;
    return std::max(
        1u,
        static_cast<unsigned>(std::min<rlim_t>(
            available_workers,
            std::numeric_limits<unsigned>::max())));
}

class StreamSource final : public OrderedShardSource {
public:
    StreamSource(
        const SourceOptions& options, std::string fallback_reason = {})
        : path_(options.path),
          input_(open_input(path_)) {
        const auto* format = hts_get_format(input_.get());
        const bool compressed =
            format != nullptr &&
            format->compression != no_compression;
        resources_ = plan_resources(
            options.total_threads, false, compressed, false,
            detect_rotational_storage(path_));
        if (resources_.hts_io_threads > 0 &&
            hts_set_threads(
                input_.get(),
                static_cast<int>(resources_.hts_io_threads)) != 0) {
            fail("Could not enable HTSlib input threads: " + path_);
        }
        header_ = read_header(input_.get(), path_);
        fallback_reason_ = std::move(fallback_reason);
    }

    bcf_hdr_t* header() const noexcept override {
        return header_.get();
    }

    std::vector<RecordPtr> next_batch(
        std::size_t maximum_records) override {
        std::vector<RecordPtr> records;
        records.reserve(maximum_records);
        while (records.size() < maximum_records) {
            RecordPtr record(bcf_init());
            if (!record) {
                fail("Could not allocate HTSlib input record");
            }
            if (bcf_read(input_.get(), header_.get(), record.get()) != 0) {
                break;
            }
            records.push_back(std::move(record));
        }
        return records;
    }

    const ResourcePlan& resources() const noexcept override {
        return resources_;
    }

    std::size_t planned_shards() const noexcept override {
        return 0;
    }

    std::string backend_name() const override {
        return "stream";
    }

    std::string description() const override {
        std::string result =
            "HTSlib ordered stream, " +
            storage_note(resources_);
        if (!fallback_reason_.empty()) {
            result += " (fallback: " + fallback_reason_ + ")";
        }
        return result;
    }

private:
    std::string path_;
    HtsFilePtr input_;
    HeaderPtr header_;
    ResourcePlan resources_;
    std::string fallback_reason_;
};

struct ShardSpec {
    std::size_t ordinal = 0;
    std::uint64_t byte_begin = 0;
    std::uint64_t byte_end = 0;
    int rid = -1;
    int index_tid = -1;
    hts_pos_t position_begin = 0;
    hts_pos_t position_end = HTS_POS_MAX;
};

struct RecordShard {
    std::size_t ordinal = 0;
    std::vector<RecordPtr> records;
};

class ParallelShardSource : public OrderedShardSource {
public:
    ParallelShardSource(
        SourceOptions options, HeaderPtr header,
        std::vector<ShardSpec> shards, ResourcePlan resources)
        : options_(std::move(options)),
          header_(std::move(header)),
          shards_(std::move(shards)),
          resources_(resources),
          maximum_ahead_(
              std::max<std::size_t>(
                  2, static_cast<std::size_t>(
                         std::max(1u, resources_.input_threads)))) {}

    ~ParallelShardSource() override {
        stop_workers();
    }

    bcf_hdr_t* header() const noexcept override {
        return header_.get();
    }

    std::vector<RecordPtr> next_batch(
        std::size_t maximum_records) override {
        std::vector<RecordPtr> batch;
        batch.reserve(maximum_records);
        while (batch.size() < maximum_records) {
            if (current_record_ < current_shard_.records.size()) {
                batch.push_back(std::move(
                    current_shard_.records[current_record_++]));
                continue;
            }
            current_shard_.records.clear();
            current_record_ = 0;

            std::unique_lock lock(mutex_);
            completed_available_.wait(lock, [&] {
                return cancelled_ || failure_ ||
                       completed_.contains(next_to_consume_) ||
                       next_to_consume_ == shards_.size();
            });
            if (failure_) {
                std::rethrow_exception(failure_);
            }
            if (cancelled_ || next_to_consume_ == shards_.size()) {
                break;
            }
            auto found = completed_.find(next_to_consume_);
            current_shard_ = std::move(found->second);
            completed_.erase(found);
            ++next_to_consume_;
            lock.unlock();
            assignment_available_.notify_all();
        }
        return batch;
    }

    const ResourcePlan& resources() const noexcept override {
        return resources_;
    }

    std::size_t planned_shards() const noexcept override {
        return shards_.size();
    }

protected:
    const SourceOptions& options() const noexcept {
        return options_;
    }

    void start_workers() {
        const unsigned count =
            std::max(1u, resources_.input_threads);
        workers_.reserve(count);
        for (unsigned worker = 0; worker < count; ++worker) {
            workers_.emplace_back([this, worker] {
                try {
                    worker_loop(worker);
                } catch (...) {
                    {
                        std::lock_guard lock(mutex_);
                        if (!failure_) {
                            failure_ = std::current_exception();
                        }
                        cancelled_ = true;
                    }
                    completed_available_.notify_all();
                    assignment_available_.notify_all();
                }
            });
        }
    }

    void stop_workers() noexcept {
        {
            std::lock_guard lock(mutex_);
            cancelled_ = true;
        }
        completed_available_.notify_all();
        assignment_available_.notify_all();
        for (auto& worker : workers_) {
            if (worker.joinable()) {
                worker.join();
            }
        }
        workers_.clear();
    }

    virtual RecordShard read_shard(
        const ShardSpec& shard, unsigned worker) = 0;

private:
    void worker_loop(unsigned worker) {
        while (true) {
            ShardSpec shard;
            {
                std::unique_lock lock(mutex_);
                assignment_available_.wait(lock, [&] {
                    return cancelled_ ||
                           (next_to_assign_ < shards_.size() &&
                            next_to_assign_ <
                                next_to_consume_ + maximum_ahead_);
                });
                if (cancelled_) {
                    return;
                }
                shard = shards_[next_to_assign_++];
            }

            RecordShard result = read_shard(shard, worker);
            {
                std::lock_guard lock(mutex_);
                if (cancelled_) {
                    return;
                }
                completed_.emplace(
                    result.ordinal, std::move(result));
            }
            completed_available_.notify_all();
        }
    }

    SourceOptions options_;
    HeaderPtr header_;
    std::vector<ShardSpec> shards_;
    ResourcePlan resources_;
    std::size_t maximum_ahead_;
    std::vector<std::thread> workers_;
    std::mutex mutex_;
    std::condition_variable completed_available_;
    std::condition_variable assignment_available_;
    std::map<std::size_t, RecordShard> completed_;
    std::size_t next_to_assign_ = 0;
    std::size_t next_to_consume_ = 0;
    bool cancelled_ = false;
    std::exception_ptr failure_;
    RecordShard current_shard_;
    std::size_t current_record_ = 0;
};

std::uint64_t find_plain_data_start(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        fail("Could not open plain VCF: " + path);
    }
    std::string line;
    while (std::getline(input, line)) {
        if (line.rfind("#CHROM", 0) == 0) {
            const auto position = input.tellg();
            if (position < 0) {
                fail("Could not find plain VCF data offset: " + path);
            }
            return static_cast<std::uint64_t>(position);
        }
        if (line.empty() || line[0] != '#') {
            fail("VCF data appeared before #CHROM header: " + path);
        }
    }
    return static_cast<std::uint64_t>(
        std::filesystem::file_size(path));
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

std::vector<ShardSpec> build_plain_shards(
    const SourceOptions& options, unsigned input_threads) {
    const std::uint64_t data_start =
        find_plain_data_start(options.path);
    const std::uint64_t file_size =
        std::filesystem::file_size(options.path);
    if (file_size <= data_start) {
        return {};
    }
    const std::uint64_t data_bytes = file_size - data_start;
    const std::uint64_t target_shard_bytes = std::clamp<std::uint64_t>(
        (256ULL * kMib) /
            std::max<unsigned>(1, input_threads),
        kMinimumPlainShardBytes,
        kTargetPlainShardBytes);
    const std::size_t by_threads =
        static_cast<std::size_t>(
            std::max(1u, input_threads)) *
        8;
    const std::size_t by_size = static_cast<std::size_t>(
        (data_bytes + target_shard_bytes - 1) /
        target_shard_bytes);
    const std::size_t maximum_useful =
        std::max<std::size_t>(
            1, static_cast<std::size_t>(
                   (data_bytes + kMinimumPlainShardBytes - 1) /
                   kMinimumPlainShardBytes));
    const std::size_t count = std::clamp<std::size_t>(
        std::max(by_threads, by_size), 1,
        std::min(kMaximumShards, maximum_useful));

    std::ifstream input(options.path, std::ios::binary);
    if (!input) {
        fail("Could not align plain VCF shards: " + options.path);
    }
    std::vector<std::uint64_t> boundaries(count + 1);
    boundaries.front() = data_start;
    boundaries.back() = file_size;
    for (std::size_t index = 1; index < count; ++index) {
        const std::uint64_t nominal =
            data_start + (data_bytes * index) / count;
        boundaries[index] = align_plain_offset(
            input, nominal, data_start, file_size);
    }

    std::vector<ShardSpec> shards;
    shards.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        if (boundaries[index] == boundaries[index + 1]) {
            continue;
        }
        shards.push_back(ShardSpec{
            .ordinal = shards.size(),
            .byte_begin = boundaries[index],
            .byte_end = boundaries[index + 1],
        });
    }
    return shards;
}

class PlainRangeSource final : public ParallelShardSource {
public:
    PlainRangeSource(
        const SourceOptions& options, HeaderPtr header,
        ResourcePlan resources)
        : ParallelShardSource(
              options, std::move(header),
              build_plain_shards(
                  options, resources.input_threads),
              resources) {
        worker_contexts_.resize(resources.input_threads);
        for (auto& context : worker_contexts_) {
            context.header = duplicate_header(this->header());
        }
        start_workers();
    }

    ~PlainRangeSource() override {
        stop_workers();
    }

    std::string backend_name() const override {
        return "plain-ranges";
    }

    std::string description() const override {
        return "plain VCF aligned byte ranges, " +
               storage_note(resources());
    }

private:
    struct WorkerContext {
        HeaderPtr header;
        kstring_t line{0, 0, nullptr};

        WorkerContext() = default;
        WorkerContext(WorkerContext&& other) noexcept
            : header(std::move(other.header)),
              line(other.line) {
            other.line = {0, 0, nullptr};
        }
        WorkerContext& operator=(WorkerContext&& other) noexcept {
            if (this != &other) {
                std::free(line.s);
                header = std::move(other.header);
                line = other.line;
                other.line = {0, 0, nullptr};
            }
            return *this;
        }
        WorkerContext(const WorkerContext&) = delete;
        WorkerContext& operator=(const WorkerContext&) = delete;
        ~WorkerContext() {
            std::free(line.s);
        }
    };

    RecordShard read_shard(
        const ShardSpec& shard, unsigned worker) override {
        auto& context = worker_contexts_[worker];
        const std::uint64_t length =
            shard.byte_end - shard.byte_begin;
        if (length >
            static_cast<std::uint64_t>(
                std::numeric_limits<std::size_t>::max())) {
            fail("Plain VCF shard is too large for this platform");
        }
        std::string bytes(static_cast<std::size_t>(length), '\0');
        const int descriptor =
            ::open(options().path.c_str(), O_RDONLY);
        if (descriptor < 0) {
            fail(
                "Could not open plain VCF shard: " +
                std::string(std::strerror(errno)));
        }
#if defined(POSIX_FADV_SEQUENTIAL)
        (void)::posix_fadvise(
            descriptor, static_cast<off_t>(shard.byte_begin),
            static_cast<off_t>(length), POSIX_FADV_SEQUENTIAL);
#endif
        std::size_t read_total = 0;
        while (read_total < bytes.size()) {
            const ssize_t count = ::pread(
                descriptor, bytes.data() + read_total,
                bytes.size() - read_total,
                static_cast<off_t>(
                    shard.byte_begin + read_total));
            if (count < 0) {
                const std::string error = std::strerror(errno);
                ::close(descriptor);
                fail("Could not read plain VCF shard: " + error);
            }
            if (count == 0) {
                ::close(descriptor);
                fail("Plain VCF ended inside an aligned shard");
            }
            read_total += static_cast<std::size_t>(count);
        }
        ::close(descriptor);

        RecordShard result;
        result.ordinal = shard.ordinal;
        std::size_t begin = 0;
        while (begin < bytes.size()) {
            std::size_t end = bytes.find('\n', begin);
            if (end == std::string::npos) {
                end = bytes.size();
            }
            std::size_t length_without_eol = end - begin;
            if (length_without_eol > 0 &&
                bytes[begin + length_without_eol - 1] == '\r') {
                --length_without_eol;
            }
            if (length_without_eol > 0) {
                context.line.l = 0;
                if (kputsn(
                        bytes.data() + begin,
                        length_without_eol,
                        &context.line) < 0) {
                    fail("Could not allocate plain VCF parse buffer");
                }
                RecordPtr record(bcf_init());
                if (!record ||
                    vcf_parse1(
                        &context.line, context.header.get(),
                        record.get()) != 0) {
                    fail(
                        "Could not parse plain VCF record in shard " +
                        std::to_string(shard.ordinal));
                }
                result.records.push_back(std::move(record));
            }
            begin = end == bytes.size() ? end : end + 1;
        }
        return result;
    }

    std::vector<WorkerContext> worker_contexts_;
};

struct IndexedLayout {
    bool bcf = false;
    std::vector<ShardSpec> shards;
};

hts_pos_t contig_length(bcf_hdr_t* header, int rid) {
    if (rid < 0 || rid >= header->n[BCF_DT_CTG] ||
        header->id[BCF_DT_CTG][rid].val == nullptr) {
        return HTS_POS_MAX;
    }
    const hts_pos_t length = static_cast<hts_pos_t>(
        header->id[BCF_DT_CTG][rid].val->info[0]);
    return length > 0 ? length : HTS_POS_MAX;
}

struct IndexedContig {
    std::string name;
    int rid = -1;
    int index_tid = -1;
    hts_pos_t length = HTS_POS_MAX;
    std::uint64_t mapped_records = 0;
};

std::vector<IndexedContig> indexed_contigs(
    const SourceOptions& options, bcf_hdr_t* header,
    bool is_bcf) {
    std::vector<IndexedContig> result;
    if (is_bcf) {
        IndexPtr index(bcf_index_load2(
            options.path.c_str(),
            options.index_path.c_str()));
        if (!index) {
            fail(
                "Could not load selected BCF CSI index: " +
                options.index_path);
        }
        for (int rid = 0; rid < header->n[BCF_DT_CTG]; ++rid) {
            const std::string name = bcf_hdr_id2name(header, rid);
            if (!options.selected_contigs.empty() &&
                !options.selected_contigs.contains(name)) {
                continue;
            }
            std::uint64_t mapped = 0;
            std::uint64_t unmapped = 0;
            (void)hts_idx_get_stat(
                index.get(), rid, &mapped, &unmapped);
            result.push_back(IndexedContig{
                .name = name,
                .rid = rid,
                .index_tid = rid,
                .length = contig_length(header, rid),
                .mapped_records = mapped,
            });
        }
        return result;
    }

    TbxPtr index(tbx_index_load2(
        options.path.c_str(), options.index_path.c_str()));
    if (!index) {
        fail(
            "Could not load selected VCF TBI/CSI index: " +
            options.index_path);
    }
    int count = 0;
    const char** names = tbx_seqnames(index.get(), &count);
    if (count > 0 && names == nullptr) {
        fail("Could not read indexed VCF contig names");
    }
    for (int tid = 0; tid < count; ++tid) {
        const std::string name = names[tid];
        if (!options.selected_contigs.empty() &&
            !options.selected_contigs.contains(name)) {
            continue;
        }
        const int rid = bcf_hdr_name2id(header, name.c_str());
        if (rid < 0) {
            std::free(const_cast<char**>(names));
            fail("Indexed contig is absent from VCF header: " + name);
        }
        std::uint64_t mapped = 0;
        std::uint64_t unmapped = 0;
        (void)hts_idx_get_stat(
            index->idx, tid, &mapped, &unmapped);
        result.push_back(IndexedContig{
            .name = name,
            .rid = rid,
            .index_tid = tid,
            .length = contig_length(header, rid),
            .mapped_records = mapped,
        });
    }
    std::free(const_cast<char**>(names));
    return result;
}

IndexedLayout build_indexed_layout(
    const SourceOptions& options, bcf_hdr_t* header,
    const htsFormat& format, unsigned input_threads) {
    IndexedLayout layout;
    layout.bcf = format.format == bcf;
    if (!layout.bcf && format.format != vcf) {
        fail("Indexed regions require VCF or BCF input");
    }
    const auto contigs =
        indexed_contigs(options, header, layout.bcf);
    if (contigs.empty()) {
        return layout;
    }

    hts_pos_t known_total_length = 0;
    std::uint64_t known_total_records = 0;
    for (const auto& contig : contigs) {
        if (contig.length != HTS_POS_MAX &&
            known_total_length <= HTS_POS_MAX - contig.length) {
            known_total_length += contig.length;
        }
        known_total_records += contig.mapped_records;
    }
    const std::size_t record_target =
        std::max<std::size_t>(
            64,
            std::min<std::size_t>(
                std::max<std::size_t>(
                    64,
                    8192 /
                        std::max<unsigned>(
                            1, input_threads)),
                std::max<std::size_t>(
                    256,
                    options.target_batch_records / 2)));
    const std::size_t by_records =
        known_total_records == 0
            ? 0
            : static_cast<std::size_t>(
                  (known_total_records + record_target - 1) /
                  record_target);
    const std::size_t target_shards = std::min<std::size_t>(
        kMaximumShards,
        std::max<std::size_t>(
            std::max(contigs.size(), by_records),
            static_cast<std::size_t>(
                std::max(1u, input_threads)) *
                8));

    for (const auto& contig : contigs) {
        hts_pos_t begin =
            options.start_position >= 1
                ? options.start_position - 1
                : 0;
        const bool open_ended =
            options.end_position ==
            std::numeric_limits<int>::max();
        hts_pos_t end =
            open_ended
                ? contig.length
                : options.end_position;
        if (end <= begin) {
            continue;
        }
        if (end == HTS_POS_MAX) {
            layout.shards.push_back(ShardSpec{
                .ordinal = layout.shards.size(),
                .rid = contig.rid,
                .index_tid = contig.index_tid,
                .position_begin = begin,
                .position_end = end,
            });
            continue;
        }
        const hts_pos_t span = end - begin;
        std::size_t windows = 1;
        if (known_total_records > 0 &&
            contig.mapped_records > 0) {
            std::uint64_t records_in_span =
                contig.mapped_records;
            if (contig.length != HTS_POS_MAX &&
                span < contig.length) {
                records_in_span = std::max<std::uint64_t>(
                    1,
                    static_cast<std::uint64_t>(
                        (static_cast<long double>(
                             contig.mapped_records) *
                         static_cast<long double>(span)) /
                        static_cast<long double>(
                            contig.length)));
            }
            const std::size_t record_windows =
                static_cast<std::size_t>(
                    (records_in_span + record_target - 1) /
                    record_target);
            const std::size_t balancing_windows =
                static_cast<std::size_t>(
                    (static_cast<long double>(target_shards) *
                     static_cast<long double>(records_in_span)) /
                    static_cast<long double>(
                        known_total_records));
            windows = std::max<std::size_t>(
                1, std::max(
                       record_windows, balancing_windows));
        } else if (known_total_length > 0) {
            windows = std::max<std::size_t>(
                1, static_cast<std::size_t>(
                       (static_cast<long double>(
                            target_shards) *
                        static_cast<long double>(span)) /
                       static_cast<long double>(
                           known_total_length)));
        }
        windows = std::min<std::size_t>(
            windows, static_cast<std::size_t>(span));
        const hts_pos_t window =
            std::max<hts_pos_t>(
                1, (span + static_cast<hts_pos_t>(windows) - 1) /
                       static_cast<hts_pos_t>(windows));
        for (hts_pos_t start = begin; start < end;) {
            const hts_pos_t stop =
                std::min(end, start + window);
            const hts_pos_t owned_end =
                open_ended && stop == end
                    ? HTS_POS_MAX
                    : stop;
            layout.shards.push_back(ShardSpec{
                .ordinal = layout.shards.size(),
                .rid = contig.rid,
                .index_tid = contig.index_tid,
                .position_begin = start,
                .position_end = owned_end,
            });
            start = stop;
        }
    }
    return layout;
}

class IndexedRegionSource final : public ParallelShardSource {
public:
    IndexedRegionSource(
        const SourceOptions& options, HeaderPtr header,
        const htsFormat& format, ResourcePlan resources)
        : IndexedRegionSource(
              options,
              duplicate_header(header.get()),
              resources,
              build_indexed_layout(
                  options, header.get(), format,
                  resources.input_threads)) {}

    ~IndexedRegionSource() override {
        stop_workers();
    }

    std::string backend_name() const override {
        return "indexed-regions";
    }

    std::string description() const override {
        const std::string adapter =
            is_bcf_
                ? "BCF CSI ordered regions"
                : "BGZF VCF TBI/CSI ordered regions";
        return adapter + " via " + options().index_path +
               ", " + storage_note(resources());
    }

private:
    IndexedRegionSource(
        const SourceOptions& options, HeaderPtr header,
        ResourcePlan resources,
        IndexedLayout layout)
        : ParallelShardSource(
              options, std::move(header),
              std::move(layout.shards), resources),
          is_bcf_(layout.bcf) {
        contexts_.resize(resources.input_threads);
        for (auto& context : contexts_) {
            context.input = open_input(options.path);
            context.header =
                read_header(context.input.get(), options.path);
            if (is_bcf_) {
                context.bcf_index.reset(
                    bcf_index_load2(
                        options.path.c_str(),
                        options.index_path.c_str()));
                if (!context.bcf_index) {
                    fail(
                        "Could not load selected BCF CSI index: " +
                        options.index_path);
                }
            } else {
                context.tbx_index.reset(
                    tbx_index_load2(
                        options.path.c_str(),
                        options.index_path.c_str()));
                if (!context.tbx_index) {
                    fail(
                        "Could not load selected VCF TBI/CSI index: " +
                        options.index_path);
                }
            }
        }
        start_workers();
    }

    struct WorkerContext {
        HtsFilePtr input;
        HeaderPtr header;
        IndexPtr bcf_index;
        TbxPtr tbx_index;
        kstring_t line{0, 0, nullptr};

        WorkerContext() = default;
        WorkerContext(WorkerContext&& other) noexcept
            : input(std::move(other.input)),
              header(std::move(other.header)),
              bcf_index(std::move(other.bcf_index)),
              tbx_index(std::move(other.tbx_index)),
              line(other.line) {
            other.line = {0, 0, nullptr};
        }
        WorkerContext& operator=(WorkerContext&& other) noexcept {
            if (this != &other) {
                std::free(line.s);
                input = std::move(other.input);
                header = std::move(other.header);
                bcf_index = std::move(other.bcf_index);
                tbx_index = std::move(other.tbx_index);
                line = other.line;
                other.line = {0, 0, nullptr};
            }
            return *this;
        }
        WorkerContext(const WorkerContext&) = delete;
        WorkerContext& operator=(const WorkerContext&) = delete;
        ~WorkerContext() {
            std::free(line.s);
        }
    };

    RecordShard read_shard(
        const ShardSpec& shard, unsigned worker) override {
        auto& context = contexts_[worker];
        RecordShard result;
        result.ordinal = shard.ordinal;
        if (is_bcf_) {
            IteratorPtr iterator(bcf_itr_queryi(
                context.bcf_index.get(), shard.index_tid,
                shard.position_begin, shard.position_end));
            if (!iterator) {
                fail(
                    "Could not query BCF indexed shard " +
                    std::to_string(shard.ordinal));
            }
            while (true) {
                RecordPtr record(bcf_init());
                if (!record) {
                    fail("Could not allocate indexed BCF record");
                }
                const int status = bcf_itr_next(
                    context.input.get(), iterator.get(),
                    record.get());
                if (status < 0) {
                    break;
                }
                if (record->rid == shard.rid &&
                    record->pos >= shard.position_begin &&
                    record->pos < shard.position_end) {
                    result.records.push_back(std::move(record));
                }
            }
            return result;
        }

        IteratorPtr iterator(tbx_itr_queryi(
            context.tbx_index.get(), shard.index_tid,
            shard.position_begin, shard.position_end));
        if (!iterator) {
            fail(
                "Could not query VCF indexed shard " +
                std::to_string(shard.ordinal));
        }
        while (tbx_itr_next(
                   context.input.get(), context.tbx_index.get(),
                   iterator.get(), &context.line) >= 0) {
            RecordPtr record(bcf_init());
            if (!record ||
                vcf_parse1(
                    &context.line, context.header.get(),
                    record.get()) != 0) {
                fail(
                    "Could not parse indexed VCF record in shard " +
                    std::to_string(shard.ordinal));
            }
            if (record->rid == shard.rid &&
                record->pos >= shard.position_begin &&
                record->pos < shard.position_end) {
                result.records.push_back(std::move(record));
            }
        }
        return result;
    }

    bool is_bcf_;
    std::vector<WorkerContext> contexts_;
};

struct AdaptiveIndexPolicy {
    bool inspect_sidecars = false;
    bool build_missing = false;
    std::string reason;
};

bool has_selective_region(const SourceOptions& options) {
    return !options.selected_contigs.empty() ||
           options.start_position >= 1 ||
           options.end_position != 0x7fffffff;
}

AdaptiveIndexPolicy choose_index_policy(
    const SourceOptions& options, bool is_bcf) {
    if (options.requested_backend == Backend::indexed_regions) {
        return {
            true, true,
            "forced indexed backend"};
    }
    if (options.requested_backend != Backend::automatic) {
        return {
            false, false,
            options.requested_backend == Backend::stream
                ? "forced stream backend"
                : "non-indexed backend requested"};
    }
    if (!options.parallel_safe) {
        return {
            false, false,
            "adaptive policy: output requires ordered streaming"};
    }
    if (has_selective_region(options)) {
        return {
            true, true,
            "adaptive policy: selective region query favors indexed access"};
    }

    if (is_bcf) {
        if (options.workload == WorkloadProfile::full_recode) {
            return {
                false, false,
                "adaptive policy: BCF full-file recode favors streaming"};
        }
        if (options.total_threads >= 4) {
            return {
                true, false,
                "adaptive policy: BCF full-scan statistics reuse an "
                "existing index from four threads"};
        }
        return {
            false, false,
            "adaptive policy: low-thread BCF full scan favors streaming"};
    }

    if (options.workload == WorkloadProfile::full_recode) {
        if (options.total_threads >= 2) {
            return {
                true, true,
                "adaptive policy: multi-thread BGZF recode favors indexed "
                "access"};
        }
        return {
            false, false,
            "adaptive policy: one-thread BGZF recode avoids index overhead"};
    }
    if (options.total_threads >= 4) {
        return {
            true, false,
            "adaptive policy: BGZF full-scan statistics reuse an existing "
            "index from four threads"};
    }
    return {
        false, false,
        "adaptive policy: low-thread BGZF full scan favors streaming"};
}

}  // namespace

void RecordDeleter::operator()(bcf1_t* record) const noexcept {
    if (record != nullptr) {
        bcf_destroy(record);
    }
}

AvailableThreads detect_available_threads() {
    constexpr const char* scheduler_variables[] = {
        "SLURM_CPUS_PER_TASK",
        "PBS_NP",
        "NSLOTS",
        "LSB_DJOB_NUMPROC",
    };
    for (const char* variable : scheduler_variables) {
        const char* value = std::getenv(variable);
        if (value == nullptr || *value == '\0') {
            continue;
        }
        char* end = nullptr;
        const unsigned long parsed = std::strtoul(value, &end, 10);
        if (end != value && *end == '\0' && parsed > 0 &&
            parsed <= std::numeric_limits<unsigned>::max()) {
            return {
                static_cast<unsigned>(parsed), variable};
        }
    }

#if defined(__linux__)
    cpu_set_t affinity;
    CPU_ZERO(&affinity);
    if (sched_getaffinity(
            0, sizeof(affinity), &affinity) == 0) {
        const int count = CPU_COUNT(&affinity);
        if (count > 0) {
            return {
                static_cast<unsigned>(count), "CPU affinity"};
        }
    }
#endif

    return {
        std::max(1u, std::thread::hardware_concurrency()),
        "hardware_concurrency"};
}

Backend parse_backend(const std::string& value) {
    if (value == "auto") {
        return Backend::automatic;
    }
    if (value == "stream") {
        return Backend::stream;
    }
    if (value == "plain") {
        return Backend::plain_ranges;
    }
    if (value == "indexed") {
        return Backend::indexed_regions;
    }
    fail(
        "Unsupported --input-backend: " + value +
        " (use auto, stream, plain, or indexed)");
}

std::string prepare_variant_index(const SourceOptions& options) {
    auto inspection = open_input(options.path);
    const htsFormat* format_pointer =
        hts_get_format(inspection.get());
    if (format_pointer == nullptr) {
        fail("Could not inspect input format: " + options.path);
    }
    const htsFormat format = *format_pointer;
    HeaderPtr header =
        read_header(inspection.get(), options.path);
    inspection.reset();

    const bool is_bgzf_variant =
        (format.format == vcf || format.format == bcf) &&
        format.compression == bgzf;
    if (!is_bgzf_variant) {
        return {};
    }

    const AdaptiveIndexPolicy policy =
        choose_index_policy(options, format.format == bcf);
    if (!policy.inspect_sidecars) {
        std::cerr << "Adaptive index: " << policy.reason << "\n";
        return {};
    }

    SidecarInspection indexes = inspect_sidecars(
        options.path, format.format == bcf, header.get());
    for (const auto& diagnostic : indexes.invalid) {
        std::cerr
            << "Index warning: protected sidecar is unusable: "
            << diagnostic << "\n";
    }
    if (!indexes.selected_path.empty()) {
        return indexes.selected_path;
    }
    if (indexes.any_present) {
        const std::string message =
            "existing CSI/TBI sidecar is unusable; refusing to "
            "overwrite it: " +
            invalid_sidecar_summary(indexes);
        if (options.requested_backend == Backend::indexed_regions) {
            fail(message);
        }
        std::cerr << "Adaptive-index warning: " << message << "\n";
        return {};
    }
    if (!policy.build_missing ||
        !std::filesystem::is_regular_file(options.path)) {
        std::cerr << "Adaptive index: " << policy.reason
                  << "; no valid index available, using stream\n";
        return {};
    }

    const IndexBuildResult build = build_csi_index(
        options, format.format == bcf, header.get());
    if (!build.success) {
        if (options.requested_backend == Backend::indexed_regions) {
            fail(
                "automatic CSI construction failed: " +
                build.detail);
        }
        std::cerr
            << "Adaptive-index warning: automatic CSI construction "
               "failed: "
            << build.detail << "\n";
        return {};
    }
    return build.index_path;
}

std::unique_ptr<OrderedShardSource> make_ordered_source(
    const SourceOptions& options) {
    SourceOptions prepared = options;
    auto inspection = open_input(options.path);
    const htsFormat* format_pointer =
        hts_get_format(inspection.get());
    if (format_pointer == nullptr) {
        fail("Could not inspect input format: " + options.path);
    }
    const htsFormat format = *format_pointer;
    HeaderPtr header =
        read_header(inspection.get(), options.path);
    inspection.reset();

    const bool is_plain_vcf =
        format.format == vcf &&
        format.compression == no_compression;
    const bool is_bgzf_variant =
        (format.format == vcf || format.format == bcf) &&
        format.compression == bgzf;
    const std::optional<bool> rotational =
        detect_rotational_storage(options.path);
    const bool automatic_plain_ranges_worthwhile =
        options.parallel_safe && options.total_threads >= 3;
    const AdaptiveIndexPolicy index_policy =
        choose_index_policy(options, format.format == bcf);
    std::string auto_index_failure;
    SidecarInspection indexes;
    if (is_bgzf_variant && index_policy.inspect_sidecars) {
        indexes = inspect_sidecars(
            options.path, format.format == bcf, header.get());
        prepared.index_path = indexes.selected_path;
        for (const auto& diagnostic : indexes.invalid) {
            std::cerr
                << "Index warning: protected sidecar is unusable: "
                << diagnostic << "\n";
        }
        if (prepared.index_path.empty() &&
            indexes.any_present) {
            auto_index_failure =
                "existing CSI/TBI sidecar is unusable; refusing "
                "to overwrite it: " +
                invalid_sidecar_summary(indexes);
        } else if (
            prepared.index_path.empty() &&
            index_policy.build_missing &&
            std::filesystem::is_regular_file(options.path)) {
            const IndexBuildResult build =
                build_csi_index(
                    options, format.format == bcf,
                    header.get());
            if (build.success) {
                prepared.index_path = build.index_path;
            } else {
                auto_index_failure =
                    "automatic CSI construction failed: " +
                    build.detail;
            }
        }
        if (!auto_index_failure.empty()) {
            std::cerr
                << "Adaptive-index warning: "
                << auto_index_failure << "\n";
        }
    } else if (
        is_bgzf_variant &&
        options.requested_backend == Backend::automatic) {
        std::cerr << "Adaptive index: "
                  << index_policy.reason << "\n";
    }
    const bool is_indexable =
        is_bgzf_variant && !prepared.index_path.empty();
    bool page_cache_prefetched = false;
    if (options.requested_backend == Backend::automatic &&
        options.parallel_safe &&
        rotational.value_or(false) &&
        (is_plain_vcf || is_indexable) &&
        can_prefetch_to_page_cache(
            options.path, options.total_threads)) {
        page_cache_prefetched =
            prefetch_to_page_cache(options.path);
    }

    Backend selected = options.requested_backend;
    if (selected == Backend::automatic) {
        if (automatic_plain_ranges_worthwhile && is_plain_vcf) {
            selected = Backend::plain_ranges;
        } else if (
            index_policy.inspect_sidecars &&
            is_indexable &&
            (!rotational.value_or(false) ||
             page_cache_prefetched)) {
            selected = Backend::indexed_regions;
        } else {
            selected = Backend::stream;
        }
    }
    if (selected == Backend::plain_ranges && !is_plain_vcf) {
        fail("--input-backend plain requires uncompressed VCF input");
    }
    if (selected == Backend::indexed_regions && !is_indexable) {
        if (!auto_index_failure.empty()) {
            fail(auto_index_failure);
        }
        fail(
            "--input-backend indexed requires BGZF VCF/BCF "
            "with a usable TBI/CSI sidecar; adaptive index "
            "construction was unavailable");
    }
    if ((selected == Backend::plain_ranges ||
         selected == Backend::indexed_regions) &&
        !options.parallel_safe) {
        fail(
            "The requested output currently requires the ordered "
            "stream backend");
    }

    if (selected == Backend::plain_ranges) {
        ResourcePlan resources =
            plan_resources(
                options.total_threads, true, false, true,
                rotational, page_cache_prefetched);
        return std::make_unique<PlainRangeSource>(
            prepared, std::move(header), resources);
    }
    if (selected == Backend::indexed_regions) {
        ResourcePlan resources =
            plan_resources(
                options.total_threads, true, false,
                format.format == vcf, rotational,
                page_cache_prefetched);
        try {
            return std::make_unique<IndexedRegionSource>(
                prepared, std::move(header), format, resources);
        } catch (const std::exception& error) {
            if (options.requested_backend != Backend::automatic) {
                throw;
            }
            return std::make_unique<StreamSource>(
                prepared, error.what());
        }
    }
    std::string stream_reason = auto_index_failure;
    if (stream_reason.empty() &&
        options.requested_backend == Backend::automatic) {
        stream_reason = index_policy.reason;
        if (index_policy.inspect_sidecars && !is_indexable) {
            stream_reason += "; no valid index available";
        } else if (
            index_policy.inspect_sidecars && is_indexable &&
            rotational.value_or(false) &&
            !page_cache_prefetched) {
            stream_reason =
                "adaptive policy: rotational storage favors streaming";
        }
    }
    return std::make_unique<StreamSource>(
        prepared, std::move(stream_reason));
}

}  // namespace vcftools_ng::input
