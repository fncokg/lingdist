#include "dp.hpp"
#include <algorithm>

#include "helpers.hpp"

namespace lingdist
{
    // Token IDs and GAP_TOKEN_ID all refer to cost.symbol_table(); string resolution is complete before entering this kernel.
    double ed_dp(const TokenSequence &vec1, const TokenSequence &vec2, const IndexedCosts &cost, const std::string &normalize_method)
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
            double cost_ins = cost.get(GAP_TOKEN_ID, vec2[i - 1]);
            for (std::size_t j = 1; j < ncol; j++)
            {
                double cost_del = cost.get(vec1[j - 1], GAP_TOKEN_ID);
                double cost_sub = cost.get(vec1[j - 1], vec2[i - 1]);
                at(i, j) = std::min({
                    at(i - 1, j - 1) + cost_sub,
                    at(i, j - 1) + cost_del,
                    at(i - 1, j) + cost_ins,
                });
            }
        }
        return normalize_edit_dist(at(len2, len1), len1, len2, normalize_method);
    }

    void ed_dp_with_dag(const TokenSequence &vec1, const TokenSequence &vec2,
                        const IndexedCosts &cost, AlignmentDAG &out)
    {
        const auto cells = (vec1.size() + 1) * (vec2.size() + 1);
        out.nrow = vec2.size() + 1;
        out.ncol = vec1.size() + 1;
        // Reset the complete active area when reusing a larger workspace.
        out.distance.assign(cells, 0.0);
        out.predecessor_mask.assign(cells, 0);
        auto &dist = out.distance;
        auto &masks = out.predecessor_mask;
        const auto ncol = out.ncol;
        for (std::size_t j = 1; j < ncol; ++j)
        {
            dist[j] = dist[j - 1] + cost.get(vec1[j - 1], GAP_TOKEN_ID);
            masks[j] = FROM_DEL;
        }
        for (std::size_t i = 1; i < out.nrow; ++i)
        {
            const auto cell = i * ncol;
            dist[cell] = dist[cell - ncol] + cost.get(GAP_TOKEN_ID, vec2[i - 1]);
            masks[cell] = FROM_INS;
        }
        for (std::size_t i = 1; i < out.nrow; ++i)
        {
            const auto insertion = cost.get(GAP_TOKEN_ID, vec2[i - 1]);
            for (std::size_t j = 1; j < ncol; ++j)
            {
                const auto cell = i * ncol + j;
                const double cost_sub = dist[cell - ncol - 1] + cost.get(vec1[j - 1], vec2[i - 1]);
                const double cost_del = dist[cell - 1] + cost.get(vec1[j - 1], GAP_TOKEN_ID);
                const double cost_ins = dist[cell - ncol] + insertion;
                const auto best = std::min({cost_sub, cost_del, cost_ins});
                dist[cell] = best;
                // Preserve the existing absolute local tolerance. Paths in this
                // graph need not lie within one EPS of the final minimum cost.
                if (double_equal(cost_sub, best))
                    masks[cell] |= FROM_SUB;
                if (double_equal(cost_del, best))
                    masks[cell] |= FROM_DEL;
                if (double_equal(cost_ins, best))
                    masks[cell] |= FROM_INS;
            }
        }
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
