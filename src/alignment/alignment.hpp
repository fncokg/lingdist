#pragma once

#include <Rcpp.h>
#include "helpers.hpp"
#include "dp.hpp"
#include "indexed_costs.hpp"

namespace lingdist
{

    struct AlignmentResult
    {
        double distance;
        std::vector<std::vector<std::pair<TokenId, TokenId>>> alignments;
    };

    // Paths contain flat DAG cell indices. Enumeration is only for display;
    // PMI training uses expected_counts and never constructs complete paths.
    using AlignmentPath = std::vector<std::size_t>;
    std::vector<AlignmentPath> get_matched_paths(const AlignmentDAG &dag, int max_paths = 10);

    AlignmentResult get_alignment_result(const TokenSequence &chars1, const TokenSequence &chars2, const IndexedCosts &cost, int max_paths = 10, const std::string &normalize_method = "longest");

    // Build alignment report list for two token sequences
    Rcpp::List get_string_alignment(const TokenSequence &chars1, const TokenSequence &chars2, const IndexedCosts &cost, int max_paths = 10, const std::string &normalize_method = "longest");

} // namespace lingdist
