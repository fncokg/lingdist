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

    DataFrame edit_dist_df(const DataFrame &data, CostTable cost, const String &delim, bool squareform, bool symmetric, bool parallel, int n_threads, bool check_missing_cost)
    {
        if (check_missing_cost)
        {
            lingdist::StrVec missing_symbols = cost.check_missing_symbols(data, delim);
            if (!missing_symbols.empty())
            {
                std::string missing_symbols_str;
                for (const auto &s : missing_symbols)
                {
                    missing_symbols_str += s + ", ";
                }
                Rprintf("Warning: The cost table is missing the following symbols found in the data: %s. Cost values for these symbols will be treated as 1.0 (deletion/insertion/substitution) or 0.0 (matching) by default.\n",
                        missing_symbols_str.c_str());
            }
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