#include "alignment.hpp"
#include "indexed_costs.hpp"
#include <algorithm>
#include <stdexcept>
#include <utility>

using namespace Rcpp;

namespace lingdist
{
    std::vector<AlignmentPath> get_matched_paths(const AlignmentDAG &dag, int max_paths)
    {
        if (max_paths < 1)
            throw std::invalid_argument("Path limit must be positive.");
        std::vector<AlignmentPath> result;
        std::vector<std::pair<std::size_t, AlignmentPath>> stack;
        stack.emplace_back(dag.distance.size() - 1, AlignmentPath{});
        while (!stack.empty())
        {
            const auto cell = stack.back().first;
            auto path = std::move(stack.back().second);
            stack.pop_back();
            path.push_back(cell);
            if (cell == 0)
            {
                std::reverse(path.begin(), path.end());
                result.push_back(std::move(path));
                if (result.size() >= static_cast<std::size_t>(max_paths))
                    break;
            }
            else
            {
                for_each_predecessor(dag, cell, [&](std::size_t predecessor, std::uint8_t)
                                     { stack.emplace_back(predecessor, path); });
            }
        }
        return result;
    }

    // Decode token IDs only when formatting the R display report.
    AlignmentResult get_alignment_result(const TokenSequence &chars1, const TokenSequence &chars2, const IndexedCosts &cost, int max_paths, const std::string &normalize_method)
    {
        AlignmentResult report;
        TokenId char1, char2;
        AlignmentDAG dag;
        ed_dp_with_dag(chars1, chars2, cost, dag);
        auto paths = get_matched_paths(dag, max_paths);
        double result = normalize_edit_dist(dag.distance.back(), chars1.size(), chars2.size(), normalize_method);
        report.distance = result;

        for (auto &path : paths)
        {
            std::vector<std::pair<TokenId, TokenId>> alignments;
            for (std::size_t j = 1; j < path.size(); j++)
            {
                const auto curr_x = path[j] / dag.ncol, curr_y = path[j] % dag.ncol;
                const auto prev_x = path[j - 1] / dag.ncol, prev_y = path[j - 1] % dag.ncol;
                char1 = curr_y == prev_y ? GAP_TOKEN_ID : chars1[curr_y - 1];
                char2 = curr_x == prev_x ? GAP_TOKEN_ID : chars2[curr_x - 1];
                alignments.emplace_back(std::make_pair(char1, char2));
            }
            report.alignments.emplace_back(alignments);
        }
        return report;
    }

    // format alignment result as R List
    List get_string_alignment(const TokenSequence &chars1, const TokenSequence &chars2, const IndexedCosts &cost, int max_paths, const std::string &normalize_method)
    {
        List report;
        List path_dfs;
        AlignmentResult result = get_alignment_result(chars1, chars2, cost, max_paths, normalize_method);
        TokenId char1, char2;

        for (auto &path : result.alignments)
        {
            StringVector chars1_col, chars2_col, operation_col;
            NumericVector cost_col;
            for (std::size_t j = 0; j < path.size(); j++)
            {
                char1 = path[j].first;
                char2 = path[j].second;
                const std::string opt = char1 == GAP_TOKEN_ID ? "insert" : (char2 == GAP_TOKEN_ID ? "delete" : (char1 == char2 ? "same" : "substitute"));
                chars1_col.push_back(Rcpp::String(cost.symbol_table().name(char1), CE_UTF8));
                chars2_col.push_back(Rcpp::String(cost.symbol_table().name(char2), CE_UTF8));
                operation_col.push_back(opt);
                cost_col.push_back(cost.get(char1, char2));
            }
            DataFrame path_df = DataFrame::create(Named("Chars1") = chars1_col, Named("Chars2") = chars2_col, Named("Operation") = operation_col, Named("Cost") = cost_col);
            path_dfs.push_back(path_df);
        }
        report["alignments"] = path_dfs;
        report["distance"] = result.distance;
        return report;
    }

} // namespace lingdist
