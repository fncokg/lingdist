#include "expected_counts.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

// All counts are in log space to avoid underflow. The number of valid optimal paths can be very large, while the weights (probabilities) are, obviously, very small.

namespace lingdist
{
    void prepare_expected_count_workspace(const TokenSequence &source, const TokenSequence &target,
                                          const IndexedCosts &cost, ExpectedCountWorkspace &workspace)
    {
        ed_dp_with_dag(source, target, cost, workspace.dag);
        const auto &dag = workspace.dag;
        auto &forward = workspace.log_forward;
        auto &backward = workspace.log_backward;
        forward.assign(dag.distance.size(), negative_infinity);
        backward.assign(dag.distance.size(), negative_infinity);
        forward.front() = 0.0; // log(1): one empty prefix.
        for (std::size_t cell = 1; cell < forward.size(); ++cell)
            for_each_predecessor(dag, cell, [&](std::size_t predecessor, std::uint8_t)
                                 { forward[cell] = logaddexp(forward[cell], forward[predecessor]); });

        backward.back() = 0.0; // One empty suffix at the terminal cell.
        for (std::size_t cell = backward.size(); cell-- > 1;)
            for_each_predecessor(dag, cell, [&](std::size_t predecessor, std::uint8_t)
                                 { backward[predecessor] = logaddexp(backward[predecessor], backward[cell]); });
    }

    double alignment_edge_probability(const ExpectedCountWorkspace &workspace,
                                      std::size_t predecessor, std::size_t cell)
    {
        // Every complete path through u->v is one prefix to u and one suffix from v. Its probability is F(u)*B(v)/F(t), without materializing K.
        const double prefix = workspace.log_forward[predecessor];
        const double suffix = workspace.log_backward[cell];
        const double total = workspace.log_forward.back();
        if (!std::isfinite(prefix) || !std::isfinite(total) ||
            std::isnan(suffix) || suffix == std::numeric_limits<double>::infinity())
            throw std::runtime_error("Invalid path count in the optimal alignment DAG.");
        if (suffix == negative_infinity)
            return 0.0;
        const double log_probability = (prefix + suffix) - total;
        const double roundoff = 64.0 * std::numeric_limits<double>::epsilon() *
                                std::max({1.0, std::abs(prefix), std::abs(suffix), std::abs(total)});
        if (std::isnan(log_probability) || log_probability > roundoff)
            throw std::runtime_error("Invalid probability in the optimal alignment DAG.");
        // Only rounding-sized positive log probabilities are corrected. Tiny
        // probabilities may underflow naturally; there is no cutoff threshold.
        return std::exp(std::min(0.0, log_probability));
    }
}
