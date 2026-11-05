#pragma once

#include "indexed_costs.hpp"
#include "numeric.hpp"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace lingdist
{

    // Basic DP distance
    double ed_dp(const TokenSequence &vec1, const TokenSequence &vec2, const IndexedCosts &cost, const std::string &normalize_method);

    // DP distance with predecessor tracking for path enumeration.

    // Here we use a bitmask to represent the predecessor directions,
    // so that all possible 8 directions can be represented in a single byte.
    // For example, 5(binary 0101) means that the current cell has two predecessors: one from the left (deletion) and one from the top-left (substitution).
    inline constexpr std::uint8_t FROM_SUB = 1;
    inline constexpr std::uint8_t FROM_DEL = 2;
    inline constexpr std::uint8_t FROM_INS = 4;

    // Rows for target, columns for source.
    // distance and predecessor_mask are stored in row-major order, with the first row and column representing the empty prefix.
    struct AlignmentDAG
    {
        std::size_t nrow = 0, ncol = 0;
        std::vector<double> distance;
        std::vector<std::uint8_t> predecessor_mask;
    };

    void ed_dp_with_dag(const TokenSequence &source, const TokenSequence &target,
                        const IndexedCosts &cost, AlignmentDAG &out);

    // Fixed order also preserves the existing display DFS traversal order.
    template <class Emit>
    void for_each_predecessor(const AlignmentDAG &dag, std::size_t cell, Emit &&emit)
    {
        const auto mask = dag.predecessor_mask[cell];
        if (mask & FROM_SUB)
            emit(cell - dag.ncol - 1, FROM_SUB);
        if (mask & FROM_DEL)
            emit(cell - 1, FROM_DEL);
        if (mask & FROM_INS)
            emit(cell - dag.ncol, FROM_INS);
    }

    double normalize_edit_dist(double dist, std::size_t len1, std::size_t len2, const std::string &method);

} // namespace lingdist
