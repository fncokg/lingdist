#pragma once

#include "dp.hpp"
#include <cmath>
#include <vector>

namespace lingdist
{
    struct ExpectedCountWorkspace
    {
        AlignmentDAG dag;
        std::vector<double> log_forward;
        std::vector<double> log_backward;
    };

    // Recompute values for this word pair and cost snapshot, reusing capacity.
    void prepare_expected_count_workspace(const TokenSequence &source, const TokenSequence &target,
                                          const IndexedCosts &cost, ExpectedCountWorkspace &workspace);

    double alignment_edge_probability(const ExpectedCountWorkspace &workspace,
                                      std::size_t predecessor, std::size_t cell);

    // The sink receives directed, unsmoothed expected occurrences. Task weights, unordered PMI aggregation and smoothing belong to the caller.
    template <class Emit>
    void accumulate_expected_alignment_counts(const TokenSequence &source, const TokenSequence &target,
                                              const IndexedCosts &cost, ExpectedCountWorkspace &workspace,
                                              Emit &&emit)
    {
        prepare_expected_count_workspace(source, target, cost, workspace);
        const auto &dag = workspace.dag;
        for (std::size_t cell = 1; cell < dag.distance.size(); ++cell)
        {
            // A locally optimal prefix may not lead to the final cell.
            if (!std::isfinite(workspace.log_backward[cell]))
                continue;
            const auto i = cell / dag.ncol, j = cell % dag.ncol;
            for_each_predecessor(dag, cell, [&](std::size_t predecessor, std::uint8_t direction)
                                 {
                const auto a = direction == FROM_INS ? GAP_TOKEN_ID : source[j - 1];
                const auto b = direction == FROM_DEL ? GAP_TOKEN_ID : target[i - 1];
                // Matches remain in F/B: removing them would change the path
                // distribution. Only the emitted PMI observations exclude them.
                if (a == b) return;
                const double probability = alignment_edge_probability(workspace, predecessor, cell);
                if (probability > 0.0) emit(a, b, probability); });
        }
    }
}
