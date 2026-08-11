#pragma once

#include <htslib/vcf.h>

#include <cstddef>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace vcftools_ng::input {

struct RecordDeleter {
    void operator()(bcf1_t* record) const noexcept;
};

using RecordPtr = std::unique_ptr<bcf1_t, RecordDeleter>;

enum class Backend {
    automatic,
    stream,
    plain_ranges,
    indexed_regions,
};

enum class WorkloadProfile {
    general,
    compact_site_statistics,
    full_recode,
};

struct AvailableThreads {
    unsigned count = 1;
    std::string source;
};

struct ResourcePlan {
    unsigned total_threads = 1;
    unsigned input_threads = 0;
    unsigned compute_threads = 1;
    unsigned hts_io_threads = 0;
    bool storage_profile_known = false;
    bool rotational_storage = false;
    bool page_cache_prefetched = false;
};

// Shared resource planner for ordered Plain/gzip/BGZF stream readers. The
// input reader and HTSlib workers are included in the total CPU budget.
ResourcePlan plan_stream_resources(
    unsigned total_threads, bool compressed);

struct SourceOptions {
    std::string path;
    Backend requested_backend = Backend::automatic;
    unsigned total_threads = 1;
    std::size_t target_batch_records = 2048;
    bool parallel_safe = true;
    WorkloadProfile workload = WorkloadProfile::general;
    std::string bcftools_path = "bcftools";
    std::string index_path;
    std::set<std::string> selected_contigs;
    int start_position = -1;
    int end_position = 0x7fffffff;
};

class OrderedShardSource {
public:
    virtual ~OrderedShardSource() = default;

    virtual bcf_hdr_t* header() const noexcept = 0;
    virtual std::vector<RecordPtr> next_batch(
        std::size_t maximum_records) = 0;
    virtual const ResourcePlan& resources() const noexcept = 0;
    virtual std::size_t planned_shards() const noexcept = 0;
    virtual std::string backend_name() const = 0;
    virtual std::string description() const = 0;
};

AvailableThreads detect_available_threads();
Backend parse_backend(const std::string& value);
std::string describe_input_format(const std::string& path);
std::string describe_storage(const std::string& path);
std::string prepare_variant_index(const SourceOptions& options);
std::unique_ptr<OrderedShardSource> make_ordered_source(
    const SourceOptions& options);

}  // namespace vcftools_ng::input
