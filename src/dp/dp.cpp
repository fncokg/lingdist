#include "dp.hpp"
#include <algorithm>

namespace lingdist
{
    // Basic DP distance, return only distance value
    double edit_dist_core_dp(const StrVec &vec1, const StrVec &vec2, const CostTable &cost)
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
            at(0, j) = at(0, j - 1) + cost.get_cost(EMPTY, vec1[j - 1]);
        }
        for (std::size_t i = 1; i < nrow; i++)
        {
            at(i, 0) = at(i - 1, 0) + cost.get_cost(EMPTY, vec2[i - 1]);
        }

        // Fill DP table
        for (std::size_t i = 1; i < nrow; i++)
        {
            for (std::size_t j = 1; j < ncol; j++)
            {
                at(i, j) = std::min({
                    at(i - 1, j - 1) + cost.get_cost(vec1[j - 1], vec2[i - 1]),
                    at(i, j - 1) + cost.get_cost(vec1[j - 1], EMPTY),
                    at(i - 1, j) + cost.get_cost(EMPTY, vec2[i - 1]),
                });
            }
        }
        return at(len2, len1);
    }

    // DP with backpointers, return full DP table with paths
    std::vector<std::vector<DistPts>> edit_dist_core_dp_record(const StrVec &vec1, const StrVec &vec2, const CostTable &cost)
    {
        std::size_t len1 = vec1.size(), len2 = vec2.size();
        std::size_t nrow = len2 + 1, ncol = len1 + 1;

        // The only difference: dist table stores (distance, source points) pairs
        std::vector<std::vector<DistPts>> dist(nrow, std::vector<DistPts>(ncol, std::make_pair(0.0, Points{})));

        for (std::size_t j = 1; j < ncol; j++)
        {
            dist[0][j] = std::make_pair(
                dist[0][j - 1].first + cost.get_cost(EMPTY, vec1[j - 1]),
                Points({std::make_pair(0, static_cast<int>(j - 1))}));
        }
        for (std::size_t i = 1; i < nrow; i++)
        {
            dist[i][0] = std::make_pair(
                dist[i - 1][0].first + cost.get_cost(EMPTY, vec2[i - 1]),
                Points({std::make_pair(static_cast<int>(i - 1), 0)}));
        }

        for (std::size_t i = 1; i < nrow; i++)
        {
            for (std::size_t j = 1; j < ncol; j++)
            {
                std::vector<double> possible_values({
                    dist[i - 1][j - 1].first + cost.get_cost(vec1[j - 1], vec2[i - 1]),
                    dist[i][j - 1].first + cost.get_cost(vec1[j - 1], EMPTY),
                    dist[i - 1][j].first + cost.get_cost(EMPTY, vec2[i - 1]),
                });
                double min_value = *std::min_element(possible_values.begin(), possible_values.end());
                Points source_points;
                if (possible_values[0] == min_value)
                    source_points.emplace_back(static_cast<int>(i - 1), static_cast<int>(j - 1));
                if (possible_values[1] == min_value)
                    source_points.emplace_back(static_cast<int>(i), static_cast<int>(j - 1));
                if (possible_values[2] == min_value)
                    source_points.emplace_back(static_cast<int>(i - 1), static_cast<int>(j));
                dist[i][j] = std::make_pair(min_value, source_points);
            }
        }
        return dist;
    }

} // namespace lingdist
