#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <vector>

#include "input_source.h"

namespace vcftools_ng {

enum class FastCapability : std::uint16_t {
    need_gt = 1U << 0U,
    need_dp = 1U << 1U,
    need_gq = 1U << 2U,
    need_ft = 1U << 3U,
    frequency_filter = 1U << 4U,
    mean_depth_filter = 1U << 5U,
};

constexpr std::uint16_t capability_mask(FastCapability capability) {
    return static_cast<std::uint16_t>(capability);
}

struct FastSiteStatPlan {
    bool freq = false;
    bool freq2 = false;
    bool counts = false;
    bool missing_site = false;
    bool site_depth = false;
    bool site_mean_depth = false;
    bool corrected_depth_arithmetic = false;
    bool site_quality = false;
    bool site_pi = false;
    int pi_window_size = 0;
    int pi_window_step = 0;
    int tajima_window_size = 0;
    std::vector<std::string> fst_population_files;
    int fst_window_size = 0;
    int fst_window_step = 0;
    bool recode = false;
    bool recode_info_all = false;
    std::function<void(std::string_view)> recode_sink;
    // Optional strict-budget allocations for fused BGZF recode. Index policy
    // still uses SourceOptions::total_threads; only active worker creation is
    // constrained by these values.
    std::optional<unsigned> input_worker_budget;
    std::optional<unsigned> stream_thread_budget;
    std::string positions_file;
    std::string exclude_positions_file;
    std::map<std::string, std::vector<int>, std::less<>> positions_to_keep;
    std::map<std::string, std::vector<int>, std::less<>> positions_to_exclude;
    std::set<std::string> samples_to_keep;
    std::set<std::string> samples_to_exclude;
    std::vector<std::string> sample_keep_files;
    std::vector<std::string> sample_exclude_files;
    std::vector<std::vector<std::size_t>> population_memberships;
    std::vector<std::uint8_t> population_roles;
    std::set<std::string, std::less<>> site_filters_to_keep;
    std::set<std::string, std::less<>> site_filters_to_remove;
    bool remove_all_filtered_sites = false;
    std::set<std::string, std::less<>> info_flags_to_keep;
    std::set<std::string, std::less<>> info_flags_to_remove;
    std::set<std::string, std::less<>> genotype_filters_to_remove;
    bool remove_all_filtered_genotypes = false;
    int min_alleles = -1;
    int max_alleles = std::numeric_limits<int>::max();
    bool remove_indels = false;
    bool keep_only_indels = false;
    double min_quality = -1.0;
    double min_genotype_quality = -1.0;
    int min_genotype_depth = -1;
    int max_genotype_depth = std::numeric_limits<int>::max();
    double min_mean_depth = -1.0;
    double max_mean_depth = std::numeric_limits<double>::max();
    double min_call_rate = 0.0;
    int max_missing_count = std::numeric_limits<int>::max();
    double min_maf = -1.0;
    double max_maf = std::numeric_limits<double>::max();
    int min_mac = -1;
    int max_mac = std::numeric_limits<int>::max();
    double min_hwe = -1.0;
    double min_non_ref_af = -1.0;
    double max_non_ref_af = std::numeric_limits<double>::max();
    double min_non_ref_af_any = -1.0;
    double max_non_ref_af_any = std::numeric_limits<double>::max();
    int min_non_ref_ac = -1;
    int max_non_ref_ac = std::numeric_limits<int>::max();
    int min_non_ref_ac_any = -1;
    int max_non_ref_ac_any = std::numeric_limits<int>::max();

    [[nodiscard]] bool counts_only() const {
        return counts && !freq && !freq2 && !missing_site &&
               !site_depth && !site_mean_depth && !site_quality &&
               !site_pi &&
               pi_window_size == 0 && tajima_window_size == 0 &&
               fst_population_files.empty() &&
               !recode;
    }

    [[nodiscard]] bool advanced_statistics_active() const {
        return pi_window_size > 0 || tajima_window_size > 0 ||
               !fst_population_files.empty();
    }

    [[nodiscard]] bool sample_selection_active() const {
        return !samples_to_keep.empty() || !samples_to_exclude.empty() ||
               !sample_keep_files.empty() ||
               !sample_exclude_files.empty();
    }

    [[nodiscard]] bool position_selection_active() const {
        return !positions_file.empty() ||
               !exclude_positions_file.empty();
    }
};

struct FastSiteStatsSummary {
    std::uint64_t total = 0;
    std::uint64_t kept = 0;
    std::size_t samples = 0;
    unsigned input_threads = 1;
    unsigned hts_io_threads = 0;
    unsigned hts_coordinator_threads = 0;
    unsigned compute_threads = 0;
    std::size_t planned_shards = 1;
    std::string backend;
    std::string description;
};

std::optional<FastSiteStatsSummary> run_fast_text_site_stats(
    const std::string& output_prefix,
    const input::SourceOptions& options,
    const FastSiteStatPlan& plan,
    std::uint16_t capabilities);

}  // namespace vcftools_ng
