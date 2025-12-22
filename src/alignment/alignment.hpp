#pragma once

#include <Rcpp.h>
#include "helpers.hpp"
#include "dp.hpp"
#include "cost_table.hpp"

namespace lingdist
{

    struct AlignmentResult
    {
        double distance;
        std::vector<std::vector<std::pair<std::string, std::string>>> alignments;
    };

    // Enumerate all matched optimal paths from recorded DP table
    std::vector<std::vector<Point>> get_matched_paths(std::vector<std::vector<DistPts>> dist, int max_paths = 10);

    AlignmentResult get_string_alignment_result(const StrVec &chars1, const StrVec &chars2, const CostTable &cost, int max_paths = 10);

    // Build alignment report list for two token sequences
    Rcpp::List get_string_alignment(const StrVec &chars1, const StrVec &chars2, const CostTable &cost, int max_paths = 10);

} // namespace lingdist
