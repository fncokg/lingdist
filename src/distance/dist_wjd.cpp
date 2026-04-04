#include <Rcpp.h>
#include <RcppThread.h>

#include "helpers.hpp"
#include "dispatcher.hpp"

using namespace Rcpp;

namespace
{

    double wjd_form(const lingdist::StrVec &form1,
                    const lingdist::StrVec &form2,
                    const std::vector<double> &cate_level_weights)
    {
        size_t len1 = form1.size(), len2 = form2.size();
        if (len1 == 0 || len2 == 0)
            return 0.0;
        size_t min_len = std::min(len1, len2);
        if (min_len > cate_level_weights.size())
        {
            Rcpp::stop("Number of categories exceeds length of cate_level_weights.");
            return NA_REAL;
        }
        for (size_t i = 0; i < min_len; i++)
        {
            if (form1[i] != form2[i])
            {
                return cate_level_weights[i];
            }
        }
        if (len1 == len2)
            return 0.0;
        // ["A","2"] vs ["A","2","a"]
        return cate_level_weights[min_len];
    }
}

namespace lingdist
{
    DataFrame wjd_df(const DataFrame &data, const String &cate_delim, const std::vector<double> &cate_weights, const String &form_strategy, const String &form_delim, const std::vector<double> &form_weights, bool detailed, bool squareform, int n_threads, bool quiet)
    {
        auto rows_vector = lingdist::split_df2(data, form_delim, cate_delim, form_strategy == "off");
        auto [lab1Col, lab2Col, row_pairs] = lingdist::get_row_pairs(data, true);

        int n_row_pairs = static_cast<int>(row_pairs.size());

        std::unique_ptr<lingdist::SafeProgressBar> bar;
        if (!quiet)
            bar = std::make_unique<lingdist::SafeProgressBar>(row_pairs.size(), 1);

        std::vector<std::vector<double>> dists_vec(row_pairs.size());

        lingdist::DistFunc dist_func = [&](const lingdist::StrVec &chars1, const lingdist::StrVec &chars2)
        { return wjd_form(chars1, chars2, cate_weights); };

        lingdist::FormsDispatcher dispatcher = make_dispatcher(form_strategy, dist_func, form_weights);

        std::function<void(std::int32_t)> loop_body = [&](std::int32_t idx)
        {
            auto [rowi, rowj] = row_pairs[idx];
            dists_vec[idx] = dist_row_pair(rows_vector[rowi], rows_vector[rowj], dispatcher);
            if (bar)
                (*bar)++;
        };
        RcppThread::parallelFor(0, static_cast<std::int32_t>(n_row_pairs), loop_body, n_threads);

        return gen_dist_df(dists_vec, lab1Col, lab2Col, detailed, as<StrVec>(data.names()), squareform, true);
    }
}