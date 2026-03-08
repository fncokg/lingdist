#pragma once

#include "helpers.hpp"
#include "cost_table.hpp"

namespace lingdist
{

    // Basic DP distance
    double edit_dist_core_dp(const StrVec &vec1, const StrVec &vec2, const CostTable &cost, const String &normalize_method);

    // DP with backpointers
    std::vector<std::vector<DistPts>> edit_dist_core_dp_record(const StrVec &vec1, const StrVec &vec2, const CostTable &cost);

    double normalize_edit_dist(double dist, const StrVec &vec1, const StrVec &vec2, const String &method);

} // namespace lingdist
