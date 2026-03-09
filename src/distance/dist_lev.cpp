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
#include "dispatcher.hpp"

using namespace Rcpp;

namespace lingdist
{

    DataFrame edit_dist_df(const DataFrame &data, const CostTable &cost, const String &delim, bool detailed, const String &normalize_method, const String &form_strategy, const String &form_delim, std::vector<double> form_weights, bool squareform, bool symmetric, int n_threads, bool check_missing_cost, bool quiet)
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
        auto rows_vector = lingdist::split_df2(data, form_delim, delim, form_strategy == "off");
        auto [lab1Col, lab2Col, row_pairs] = lingdist::get_row_pairs(data, symmetric);

        int n_row_pairs = static_cast<int>(row_pairs.size());

        std::unique_ptr<lingdist::SafeProgressBar> bar;
        if (!quiet)
            bar = std::make_unique<lingdist::SafeProgressBar>(row_pairs.size(), 1);
        std::vector<std::vector<double>> dists_vec(row_pairs.size());

        lingdist::DistFunc dist_func = [&](const lingdist::StrVec &chars1, const lingdist::StrVec &chars2)
        { return lingdist::edit_dist_core_dp(chars1, chars2, cost, normalize_method); };

        lingdist::FormsDispatcher dispatcher = make_dispatcher(form_strategy, dist_func, form_weights);

        std::function<void(std::int32_t)> loop_body = [&](std::int32_t idx)
        {
            auto [rowi, rowj] = row_pairs[idx];
            dists_vec[idx] = dist_row_pair(rows_vector[rowi], rows_vector[rowj], dispatcher);
            if (bar)
                (*bar)++;
        };
        RcppThread::parallelFor(0, static_cast<std::int32_t>(n_row_pairs), loop_body, n_threads);

        return gen_dist_df(dists_vec, lab1Col, lab2Col, detailed, as<StrVec>(data.names()), squareform, symmetric);
    }

}