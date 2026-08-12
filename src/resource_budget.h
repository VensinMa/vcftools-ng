#pragma once

#include <algorithm>
#include <limits>

namespace vcftools_ng::resources {

struct FusedBgzfRecodePlan {
    unsigned total_threads = 1;
    // Logical input concurrency shown to users.  At one thread this is the
    // calling thread; input_worker_threads counts only background workers.
    unsigned input_threads = 1;
    unsigned input_worker_threads = 0;
    // Budget passed to a sequential compressed-input pipeline.  That path
    // shares the caller between reading, committing, and ordered output.
    unsigned stream_thread_budget = 1;
    unsigned compression_threads = 0;
    // Background coordinators.  The calling thread is the coordinator for
    // multi-thread plans but is not an additional thread at total_threads=1.
    unsigned coordinator_threads = 0;
};

struct GenericOutputPlan {
    unsigned total_threads = 1;
    unsigned pipeline_threads = 1;
    unsigned vcf_compression_threads = 0;
    unsigned bcf_compression_threads = 0;
    unsigned bcf_serialization_threads = 0;
};

inline unsigned integer_sqrt(unsigned value) noexcept {
    unsigned root = 0;
    while (root < value &&
           root + 1u <= value / (root + 1u)) {
        ++root;
    }
    return root;
}

inline FusedBgzfRecodePlan plan_fused_bgzf_recode(
    unsigned requested_threads) noexcept {
    const unsigned total_threads = std::max(1u, requested_threads);
    if (total_threads == 1) {
        return FusedBgzfRecodePlan{
            .total_threads = 1,
            .input_threads = 1,
            .input_worker_threads = 0,
            .stream_thread_budget = 1,
            .compression_threads = 0,
            .coordinator_threads = 0,
        };
    }

    // Indexed input workers frequently block on storage while BGZF output
    // workers block on ordered handoff. Keeping both pools independently
    // ready, under the process-wide N-CPU affinity enforced by main(), avoids
    // leaving cores idle without allowing the workflow to consume more than
    // N CPU cores. A 4-32 CPU real-data sweep selected approximately
    // 2*sqrt(N) compression workers; this grows conservatively on large hosts.
    const unsigned compression_threads =
        std::min(total_threads,
                 std::max(1u, 2u * integer_sqrt(total_threads)));
    return FusedBgzfRecodePlan{
        .total_threads = total_threads,
        .input_threads = total_threads,
        .input_worker_threads = total_threads,
        .stream_thread_budget = total_threads,
        .compression_threads = compression_threads,
        .coordinator_threads = 1,
    };
}

inline bool valid(const FusedBgzfRecodePlan& plan) noexcept {
    if (plan.total_threads == 0 || plan.input_threads == 0 ||
        plan.stream_thread_budget == 0) {
        return false;
    }
    if (plan.total_threads == 1) {
        return plan.input_worker_threads == 0 &&
               plan.compression_threads == 0 &&
               plan.coordinator_threads == 0 &&
               plan.stream_thread_budget == 1;
    }
    return plan.input_worker_threads == plan.total_threads &&
           plan.input_threads == plan.total_threads &&
           plan.stream_thread_budget == plan.total_threads &&
           plan.compression_threads > 0 &&
           plan.compression_threads <= plan.total_threads &&
           plan.coordinator_threads == 1;
}

inline GenericOutputPlan plan_generic_outputs(
    unsigned requested_threads, bool vcf_bgzf, bool bcf) noexcept {
    const unsigned total_threads = std::max(1u, requested_threads);
    GenericOutputPlan plan{.total_threads = total_threads};
    if (!vcf_bgzf && !bcf) {
        plan.pipeline_threads = total_threads;
        return plan;
    }

    // BCF serialization feeds deterministic compression through a pipe when
    // the budget permits.  At one thread it uses a synchronous temporary
    // byte stream instead of creating a hidden reader thread.
    if (bcf && total_threads >= 2) {
        plan.bcf_serialization_threads = 1;
    }
    const unsigned reserved = 1 + plan.bcf_serialization_threads;
    const unsigned available =
        total_threads > reserved ? total_threads - reserved : 0;
    const unsigned sinks = static_cast<unsigned>(vcf_bgzf) +
                           static_cast<unsigned>(bcf);
    unsigned compression = sinks == 0
                               ? 0
                               : std::min(
                                     integer_sqrt(total_threads - 1),
                                     available);
    while (compression > 0) {
        if (vcf_bgzf &&
            (!bcf || plan.vcf_compression_threads <=
                         plan.bcf_compression_threads)) {
            ++plan.vcf_compression_threads;
        } else {
            ++plan.bcf_compression_threads;
        }
        --compression;
    }
    plan.pipeline_threads =
        total_threads - plan.bcf_serialization_threads -
        plan.vcf_compression_threads -
        plan.bcf_compression_threads;
    return plan;
}

inline bool valid(const GenericOutputPlan& plan) noexcept {
    return plan.total_threads > 0 && plan.pipeline_threads > 0 &&
           plan.pipeline_threads + plan.vcf_compression_threads +
                   plan.bcf_compression_threads +
                   plan.bcf_serialization_threads ==
               plan.total_threads;
}

}  // namespace vcftools_ng::resources
