#pragma once

#include "helpers.hpp"
#include "indexed_costs.hpp"

namespace lingdist
{

    // Basic DP distance
    double edit_dist_core_dp(const TokenSequence &vec1, const TokenSequence &vec2, const IndexedCosts &cost, const std::string &normalize_method);

    // DP with backpointers
    std::vector<std::vector<DistPts>> edit_dist_core_dp_record(const TokenSequence &vec1, const TokenSequence &vec2, const IndexedCosts &cost);

    double normalize_edit_dist(double dist, std::size_t len1, std::size_t len2, const std::string &method);

} // namespace lingdist
