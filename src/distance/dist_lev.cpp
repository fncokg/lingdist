#include <Rcpp.h>
#include <RcppThread.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>
#include <string>
#include <utility>
#include <memory>

#include "helpers.hpp"
#include "cost_table.hpp"
#include "dp.hpp"
#include "alignment.hpp"
#include "dist.hpp"

using namespace Rcpp;

namespace
{

    std::vector<double> edit_dist_row_vec(const std::vector<lingdist::StrVec> &row1, const std::vector<lingdist::StrVec> &row2, const lingdist::CostTable &cost, const String &normalize_method)
    {
        std::size_t ncols = row1.size();
        std::vector<double> dists(ncols, NA_REAL);
        for (std::size_t coli = 0; coli < ncols; coli++)
        {
            const auto &chars1 = row1[coli];
            const auto &chars2 = row2[coli];
            if (!chars1.empty() && !chars2.empty())
            {
                dists[coli] = lingdist::edit_dist_core_dp(chars1, chars2, cost, normalize_method);
            }
        }
        return dists;
    }

    double edit_dist_row(const std::vector<lingdist::StrVec> &row1,
                         const std::vector<lingdist::StrVec> &row2,
                         const lingdist::CostTable &cost,
                         const String &normalize_method)
    {
        std::vector<double> dists = edit_dist_row_vec(row1, row2, cost, normalize_method);
        double sum_dist = 0.0, nwords = 0.0;
        for (auto &dist : dists)
        {
            if (!Rcpp::NumericVector::is_na(dist))
            {
                sum_dist += dist;
                nwords += 1.0;
            }
        }
        if (nwords <= 0.0)
            return NA_REAL;
        return sum_dist / nwords;
    }

}

namespace lingdist
{

    DataFrame edit_dist_df(const DataFrame &data, const CostTable &cost, const String &delim, bool detailed, const String &normalize_method, bool squareform, bool symmetric, bool parallel, int n_threads, bool check_missing_cost, bool quiet)
    {
        if (normalize_method != "longest" && normalize_method != "none")
        {
            stop("Unsupported normalize_method. Supported values are 'longest' and 'none'.");
        }

        if (check_missing_cost && !cost.is_fast)
        {
            lingdist::StrVec missing_symbols = cost.check_missing_symbols(data, delim);
            if (!missing_symbols.empty())
            {
                std::string missing_symbols_str;
                for (const auto &s : missing_symbols)
                {
                    missing_symbols_str += s + ", ";
                }
                warning("Warning: The cost table is missing the following symbols found in the data: %s. Cost values for these symbols will be treated as 1.0 (deletion/insertion/substitution) or 0.0 (matching) by default.\n",
                        missing_symbols_str.c_str());
            }
        }
        auto rows_vector = lingdist::split_df(data, delim);
        auto [lab1Col, lab2Col, row_pairs] = lingdist::get_row_pairs(data, symmetric);

        int n_row_pairs = static_cast<int>(row_pairs.size());

        std::unique_ptr<lingdist::SafeProgressBar> bar;
        if (!quiet)
            bar = std::make_unique<lingdist::SafeProgressBar>(row_pairs.size(), 1);
        std::function<void(std::int32_t)> loop_body;
        std::vector<std::vector<double>> dists_vec(row_pairs.size());
        std::vector<double> dists(row_pairs.size());
        if (detailed)
        {

            loop_body = [&](std::int32_t idx)
            {
                auto [rowi, rowj] = row_pairs[idx];
                dists_vec[idx] = edit_dist_row_vec(rows_vector[rowi], rows_vector[rowj], cost, normalize_method);
                if (bar)
                    (*bar)++;
            };
        }
        else
        {

            loop_body = [&](std::int32_t idx)
            {
                auto [rowi, rowj] = row_pairs[idx];
                dists[idx] = edit_dist_row(rows_vector[rowi], rows_vector[rowj], cost, normalize_method);
                if (bar)
                    (*bar)++;
            };
        }

        if (parallel)
        {
            RcppThread::parallelFor(0, static_cast<std::int32_t>(n_row_pairs), loop_body, n_threads);
        }
        else
        {
            lingdist::singleFor(0, n_row_pairs, loop_body);
        }
        if (detailed)
        {
            return gen_dist_df_detailed(dists_vec, lab1Col, lab2Col, as<StrVec>(data.names()));
        }
        else
        {
            return gen_dist_df(dists, lab1Col, lab2Col, squareform, symmetric);
        }
    }

}