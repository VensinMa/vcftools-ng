#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>

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

    [[nodiscard]] bool counts_only() const {
        return counts && !freq && !freq2 && !missing_site &&
               !site_depth && !site_mean_depth && !site_quality;
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
