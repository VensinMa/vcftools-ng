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

struct FastSiteStatPlan {
    bool freq = false;
    bool freq2 = false;
    bool counts = false;
    bool missing_site = false;
    bool site_depth = false;
    bool site_mean_depth = false;
    bool site_quality = false;
    int pi_window_size = 0;
    int pi_window_step = 0;
    int tajima_window_size = 0;
    std::vector<std::string> fst_population_files;
    int fst_window_size = 0;
    int fst_window_step = 0;
    bool recode = false;
    bool recode_info_all = false;
    std::function<void(std::string_view)> recode_sink;
    std::string positions_file;
    std::string exclude_positions_file;
    std::map<std::string, std::vector<int>, std::less<>> positions_to_keep;
    std::map<std::string, std::vector<int>, std::less<>> positions_to_exclude;
    std::set<std::string> samples_to_keep;
    std::set<std::string> samples_to_exclude;
    std::vector<std::string> sample_keep_files;
    std::vector<std::string> sample_exclude_files;
    std::vector<std::vector<std::size_t>> population_memberships;
    int min_alleles = -1;
    int max_alleles = std::numeric_limits<int>::max();
    double min_quality = -1.0;
    double min_genotype_quality = -1.0;
    double min_mean_depth = -1.0;
    double min_call_rate = 0.0;
    double min_maf = -1.0;

    [[nodiscard]] bool counts_only() const {
        return counts && !freq && !freq2 && !missing_site &&
               !site_depth && !site_mean_depth && !site_quality &&
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
    std::size_t planned_shards = 1;
    std::string backend;
    std::string description;
};

std::optional<FastSiteStatsSummary> run_fast_text_site_stats(
    const std::string& output_prefix,
    const input::SourceOptions& options,
    const FastSiteStatPlan& plan);

}  // namespace vcftools_ng
