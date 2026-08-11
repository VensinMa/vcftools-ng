#pragma once

#include <cstdint>
#include <limits>
#include <stdexcept>

namespace vcftools_ng {

// VCFtools 0.1.17 exposes unsigned-32-bit wrapping in --site-depth and
// --site-mean-depth. Keep those exact values while also maintaining checked
// 64-bit sums for the explicit corrected-depth extension.
class DepthAccumulator {
public:
    void add(std::uint32_t depth) {
        legacy_sum_ += depth;
        const std::uint64_t square =
            static_cast<std::uint64_t>(depth) * depth;
        legacy_sumsq_ += static_cast<std::uint32_t>(square);

        if (corrected_sum_ >
            std::numeric_limits<std::uint64_t>::max() - depth) {
            throw std::overflow_error(
                "Corrected site-depth sum exceeds 64-bit range");
        }
        if (corrected_sumsq_ >
            std::numeric_limits<std::uint64_t>::max() - square) {
            throw std::overflow_error(
                "Corrected site-depth sum of squares exceeds 64-bit range");
        }
        corrected_sum_ += depth;
        corrected_sumsq_ += square;
        if (count_ == std::numeric_limits<std::uint64_t>::max()) {
            throw std::overflow_error(
                "Corrected site-depth observation count exceeds 64-bit range");
        }
        ++count_;
    }

    [[nodiscard]] std::uint32_t legacy_sum() const noexcept {
        return legacy_sum_;
    }
    [[nodiscard]] std::uint32_t legacy_sumsq() const noexcept {
        return legacy_sumsq_;
    }
    [[nodiscard]] std::uint64_t corrected_sum() const noexcept {
        return corrected_sum_;
    }
    [[nodiscard]] std::uint64_t corrected_sumsq() const noexcept {
        return corrected_sumsq_;
    }
    [[nodiscard]] std::uint64_t count() const noexcept {
        return count_;
    }

private:
    std::uint32_t legacy_sum_ = 0;
    std::uint32_t legacy_sumsq_ = 0;
    std::uint64_t corrected_sum_ = 0;
    std::uint64_t corrected_sumsq_ = 0;
    std::uint64_t count_ = 0;
};

}  // namespace vcftools_ng
