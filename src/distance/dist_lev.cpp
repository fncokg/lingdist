#include <Rcpp.h>
#include <RcppThread.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>
#include <string>
#include <utility>

#include "helpers.hpp"
#include "cost_table.hpp"
#include "dp.hpp"
#include "alignment.hpp"
#include "dist.hpp"

using namespace Rcpp;

namespace
{

    double edit_dist_row(const std::vector<lingdist::StrVec> &row1,
                         const std::vector<lingdist::StrVec> &row2,
                         const lingdist::CostTable &cost)
    {
        std::size_t ncols = row1.size();
        double ttlDist = 0.0;
        int nword = 0;
        for (std::size_t coli = 0; coli < ncols; coli++)
        {
            const auto &chars1 = row1[coli];
            const auto &chars2 = row2[coli];
            if (!chars1.empty() && !chars2.empty())
            {
                ttlDist += lingdist::edit_dist_core_dp(chars1, chars2, cost);
                nword += 1;
            }
        }
        if (nword <= 0)
            return -1.0;
        return ttlDist / static_cast<double>(nword);
    }

}

namespace lingdist
{

    DataFrame edit_dist_df(const DataFrame &data, Nullable<DataFrame> cost_mat, const String &delim, bool squareform, bool symmetric, bool parallel, int n_threads)
    {
        lingdist::CostTable cost;
        if (cost_mat.isNotNull())
        {
            DataFrame cost_mat_ = DataFrame(cost_mat);
            cost = lingdist::build_cost_table(cost_mat_);
        }
        else
        {
            lingdist::StrVec chars = lingdist::get_all_unique_chars(data, delim);
            chars.push_back(lingdist::EMPTY);
            cost = lingdist::build_default_cost_table(chars);
        }

        auto rows_vector = lingdist::split_df(data, delim);
        auto [lab1Col, lab2Col, row_pairs] = lingdist::get_row_pairs(data, symmetric);

        int n_row_pairs = static_cast<int>(row_pairs.size());

        std::vector<double> dists(row_pairs.size());
        RcppThread::ProgressBar bar(row_pairs.size(), 1);
        std::function<void(std::int32_t)> loop_body = [&](std::int32_t idx)
        {
            auto [rowi, rowj] = row_pairs[idx];
            dists[idx] = edit_dist_row(rows_vector[rowi], rows_vector[rowj], cost);
            bar++;
        };
        if (parallel)
        {
            RcppThread::parallelFor(0, static_cast<std::int32_t>(n_row_pairs), loop_body, n_threads);
        }
        else
        {
            lingdist::singleFor(0, n_row_pairs, loop_body);
        }
        DataFrame result = DataFrame::create(Named("lab1") = lab1Col, Named("lab2") = lab2Col, Named("dist") = NumericVector::import(dists.begin(), dists.end()));
        if (squareform)
        {
            result = lingdist::long2squareform(result, symmetric);
        }
        return result;
    }

}