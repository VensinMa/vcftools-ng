#include "exact_hwe.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace vcftools_ng {

ExactHweProbabilities exact_hwe_probabilities(
    int observed_hets, int observed_hom_ref, int observed_hom_alt) {
    ExactHweProbabilities result;
    if (observed_hom_ref + observed_hom_alt + observed_hets == 0) {
        return result;
    }

    const int observed_hom_common =
        std::max(observed_hom_ref, observed_hom_alt);
    const int observed_hom_rare =
        std::min(observed_hom_ref, observed_hom_alt);
    const int rare_copies = 2 * observed_hom_rare + observed_hets;
    const int genotypes =
        observed_hets + observed_hom_common + observed_hom_rare;

    std::vector<double> heterozygote_probabilities(
        static_cast<std::size_t>(rare_copies + 1), 0.0);
    const auto midpoint_numerator =
        static_cast<std::int64_t>(rare_copies) *
        (2LL * genotypes - rare_copies);
    int midpoint = static_cast<int>(
        midpoint_numerator / (2LL * genotypes));
    if ((rare_copies & 1) ^ (midpoint & 1)) {
        ++midpoint;
    }

    int current_hets = midpoint;
    int current_hom_rare = (rare_copies - midpoint) / 2;
    int current_hom_common =
        genotypes - current_hets - current_hom_rare;
    heterozygote_probabilities[midpoint] = 1.0;
    double sum = 1.0;

    for (current_hets = midpoint; current_hets > 1; current_hets -= 2) {
        heterozygote_probabilities[current_hets - 2] =
            heterozygote_probabilities[current_hets] * current_hets *
            (current_hets - 1.0) /
            (4.0 * (current_hom_rare + 1.0) *
             (current_hom_common + 1.0));
        sum += heterozygote_probabilities[current_hets - 2];
        ++current_hom_rare;
        ++current_hom_common;
    }

    current_hets = midpoint;
    current_hom_rare = (rare_copies - midpoint) / 2;
    current_hom_common =
        genotypes - current_hets - current_hom_rare;
    for (current_hets = midpoint; current_hets <= rare_copies - 2;
         current_hets += 2) {
        heterozygote_probabilities[current_hets + 2] =
            heterozygote_probabilities[current_hets] * 4.0 *
            current_hom_rare * current_hom_common /
            ((current_hets + 2.0) * (current_hets + 1.0));
        sum += heterozygote_probabilities[current_hets + 2];
        --current_hom_rare;
        --current_hom_common;
    }

    for (double& probability : heterozygote_probabilities) {
        probability /= sum;
    }

    result.heterozygote_excess =
        heterozygote_probabilities[observed_hets];
    for (int index = observed_hets + 1; index <= rare_copies; ++index) {
        result.heterozygote_excess +=
            heterozygote_probabilities[index];
    }
    result.heterozygote_deficit =
        heterozygote_probabilities[observed_hets];
    for (int index = observed_hets - 1; index >= 0; --index) {
        result.heterozygote_deficit +=
            heterozygote_probabilities[index];
    }

    result.two_sided = 0.0;
    for (const double probability : heterozygote_probabilities) {
        if (probability <= heterozygote_probabilities[observed_hets]) {
            result.two_sided += probability;
        }
    }
    result.two_sided = std::min(result.two_sided, 1.0);
    return result;
}

double exact_hwe_pvalue(
    int observed_hets, int observed_hom_ref, int observed_hom_alt) {
    return exact_hwe_probabilities(
               observed_hets, observed_hom_ref, observed_hom_alt)
        .two_sided;
}

}  // namespace vcftools_ng
