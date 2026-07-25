#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>

#include "input_source.h"

namespace vcftools_ng {

struct FastCountsSummary {
    std::uint64_t total = 0;
    std::uint64_t kept = 0;
    std::size_t samples = 0;
    unsigned input_threads = 1;
    unsigned hts_io_threads = 0;
    std::size_t planned_shards = 1;
    std::string backend;
    std::string description;
};

std::optional<FastCountsSummary> run_fast_text_counts(
    const std::string& output_prefix,
    const input::SourceOptions& options);

}  // namespace vcftools_ng
