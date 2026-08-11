#pragma once

#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>

namespace vcftools_ng {

// Windowed population statistics require nondecreasing positions and
// contiguous chromosome segments. Rejecting an unsorted stream prevents a
// plausible-looking result whose meaning depends on implementation details.
class OrderedSiteValidator {
public:
    void observe(std::string_view chromosome, std::int64_t position) {
        if (current_chromosome_.empty()) {
            current_chromosome_.assign(chromosome);
            last_position_ = position;
            return;
        }
        if (chromosome == current_chromosome_) {
            if (position < last_position_) {
                throw std::runtime_error(
                    "Ordered analysis requires a position-sorted input; "
                    "position decreased within chromosome " +
                    current_chromosome_);
            }
            last_position_ = position;
            return;
        }
        finished_chromosomes_.insert(current_chromosome_);
        if (finished_chromosomes_.contains(chromosome)) {
            throw std::runtime_error(
                "Ordered analysis requires contiguous chromosome segments; "
                "chromosome reappeared after a later segment: " +
                std::string(chromosome));
        }
        current_chromosome_.assign(chromosome);
        last_position_ = position;
    }

private:
    std::set<std::string, std::less<>> finished_chromosomes_;
    std::string current_chromosome_;
    std::int64_t last_position_ = 0;
};

}  // namespace vcftools_ng
