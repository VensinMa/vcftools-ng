#include "resource_budget.h"
#include "input_source.h"

#include <array>
#include <climits>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

}  // namespace

int main() {
    using vcftools_ng::resources::plan_fused_bgzf_recode;
    using vcftools_ng::resources::plan_generic_outputs;
    using vcftools_ng::resources::valid;

    constexpr std::array<unsigned, 18> budgets{
        0, 1, 2, 3, 4, 7, 8, 12, 16, 24, 28, 32,
        64, 128, 256, 512, 1024, UINT_MAX};
    unsigned previous_compression = 0;
    unsigned previous_input = 0;
    for (const unsigned requested : budgets) {
        const auto plan = plan_fused_bgzf_recode(requested);
        require(valid(plan), "invalid representative resource plan");
        require(
            plan.total_threads == (requested == 0 ? 1u : requested),
            "requested total was not preserved");
        require(
            plan.compression_threads >= previous_compression,
            "compression allocation is not monotonic");
        if (plan.total_threads >= 3) {
            require(
                plan.input_worker_threads >= previous_input,
                "input allocation is not monotonic");
        }
        previous_compression = plan.compression_threads;
        previous_input = plan.input_worker_threads;
    }

    for (unsigned budget = 1; budget <= 65536; ++budget) {
        const auto plan = plan_fused_bgzf_recode(budget);
        require(valid(plan), "invalid exhaustive resource plan");
        require(
            plan.stream_thread_budget == plan.total_threads,
            "fused stream did not receive the CPU budget");
        for (const bool vcf : {false, true}) {
            for (const bool bcf : {false, true}) {
                const auto generic =
                    plan_generic_outputs(budget, vcf, bcf);
                require(
                    valid(generic),
                    "invalid generic output resource plan");
            }
        }
        for (const bool compressed : {false, true}) {
            for (const bool bcf : {false, true}) {
                const auto stream =
                    vcftools_ng::input::plan_stream_resources(
                        budget, compressed, bcf);
                require(
                    stream.total_threads == budget,
                    "stream planner changed total budget");
                require(
                    stream.input_threads + stream.compute_threads +
                            stream.hts_io_threads +
                            stream.hts_coordinator_threads ==
                        budget,
                    "fused stream plan exceeds total budget");
            }
        }
    }
    std::cout << "RESOURCE_BUDGET_PASS\n";
}
