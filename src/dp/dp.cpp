#include "dp.hpp"
#include <algorithm>

#include "helpers.hpp"

namespace lingdist
{
    // Token IDs and GAP_TOKEN_ID all refer to cost.symbol_table(); string
    // resolution is complete before entering this kernel.
    double edit_dist_core_dp(const TokenSequence &vec1, const TokenSequence &vec2, const IndexedCosts &cost, const std::string &normalize_method)
    {
        std::size_t len1 = vec1.size(), len2 = vec2.size();
        std::size_t nrow = len2 + 1, ncol = len1 + 1;

        std::vector<double> dist(nrow * ncol, 0.0);

        // Use 1d `dist` as 2d array
        auto at = [&](std::size_t r, std::size_t c) -> double &
        { return dist[r * ncol + c]; };

        // Initialize DP table: the first row and first column
        at(0, 0) = 0.0;
        for (std::size_t j = 1; j < ncol; j++)
        {
            // delete the first j symbols of vec1 to match an empty vec2
            at(0, j) = at(0, j - 1) + cost.get(vec1[j - 1], GAP_TOKEN_ID);
        }
        for (std::size_t i = 1; i < nrow; i++)
        {
            // insert the first i symbols of vec2 to match an empty vec1
            at(i, 0) = at(i - 1, 0) + cost.get(GAP_TOKEN_ID, vec2[i - 1]);
        }

        // Fill DP table
        for (std::size_t i = 1; i < nrow; i++)
        {
            double cost_insert = cost.get(GAP_TOKEN_ID, vec2[i - 1]);
            for (std::size_t j = 1; j < ncol; j++)
            {
                double cost_delete = cost.get(vec1[j - 1], GAP_TOKEN_ID);
                double cost_sub = cost.get(vec1[j - 1], vec2[i - 1]);
                at(i, j) = std::min({
                    at(i - 1, j - 1) + cost_sub,
                    at(i, j - 1) + cost_delete,
                    at(i - 1, j) + cost_insert,
                });
            }
        }
        return normalize_edit_dist(at(len2, len1), len1, len2, normalize_method);
    }

    // DP with backpointers, return full DP table with paths
    std::vector<std::vector<DistPts>> edit_dist_core_dp_record(const TokenSequence &vec1, const TokenSequence &vec2, const IndexedCosts &cost)
    {
        std::size_t len1 = vec1.size(), len2 = vec2.size();
        std::size_t nrow = len2 + 1, ncol = len1 + 1;

        // The only difference: dist table stores (distance, source points) pairs
        std::vector<std::vector<DistPts>> dist(nrow, std::vector<DistPts>(ncol, std::make_pair(0.0, Points{})));

        for (std::size_t j = 1; j < ncol; j++)
        {
            dist[0][j] = std::make_pair(
                dist[0][j - 1].first + cost.get(vec1[j - 1], GAP_TOKEN_ID),
                Points({std::make_pair(0, static_cast<int>(j - 1))}));
        }
        for (std::size_t i = 1; i < nrow; i++)
        {
            dist[i][0] = std::make_pair(
                dist[i - 1][0].first + cost.get(GAP_TOKEN_ID, vec2[i - 1]),
                Points({std::make_pair(static_cast<int>(i - 1), 0)}));
        }

        for (std::size_t i = 1; i < nrow; i++)
        {
            double cost_insert = cost.get(GAP_TOKEN_ID, vec2[i - 1]);
            for (std::size_t j = 1; j < ncol; j++)
            {
                double cost_delete = cost.get(vec1[j - 1], GAP_TOKEN_ID);
                double cost_sub = cost.get(vec1[j - 1], vec2[i - 1]);
                std::vector<double> possible_values({
                    dist[i - 1][j - 1].first + cost_sub,
                    dist[i][j - 1].first + cost_delete,
                    dist[i - 1][j].first + cost_insert,
                });
                double min_value = *std::min_element(possible_values.begin(), possible_values.end());
                Points source_points;
                if (double_equal(possible_values[0], min_value))
                    source_points.emplace_back(static_cast<int>(i - 1), static_cast<int>(j - 1));
                if (double_equal(possible_values[1], min_value))
                    source_points.emplace_back(static_cast<int>(i), static_cast<int>(j - 1));
                if (double_equal(possible_values[2], min_value))
                    source_points.emplace_back(static_cast<int>(i - 1), static_cast<int>(j));
                dist[i][j] = std::make_pair(min_value, source_points);
            }
        }
        return dist;
    }

    // For performance reasons, we do not check the validity of `method` here. Checking should be done at the higher level before calling this function.
    double normalize_edit_dist(double dist, std::size_t len1, std::size_t len2, const std::string &method)
    {
        // fastest case: no normalization
        if (method == "none")
        {
            return dist; // no normalization
        }

        if (method == "longest")
        {
            size_t max_len = std::max(len1, len2);
            if (max_len == 0)
                return 0.0; // both strings are empty
            return dist / static_cast<double>(max_len);
        }
        else
        {
            return dist; // no normalization
        }
    }

} // namespace lingdist
