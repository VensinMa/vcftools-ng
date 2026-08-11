#pragma once

namespace vcftools_ng {

struct ExactHweProbabilities {
    double two_sided = 1.0;
    double heterozygote_deficit = 1.0;
    double heterozygote_excess = 1.0;
};

ExactHweProbabilities exact_hwe_probabilities(
    int observed_hets, int observed_hom_ref, int observed_hom_alt);

double exact_hwe_pvalue(
    int observed_hets, int observed_hom_ref, int observed_hom_alt);

}  // namespace vcftools_ng
