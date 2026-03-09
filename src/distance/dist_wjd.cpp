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
        return cate_level_weights[min_len - 1];
    }

    std::vector<double> wjd_row_vec(const std::vector<std::vector<lingdist::StrVec>> &row1,
                                    const std::vector<std::vector<lingdist::StrVec>> &row2,
                                    const std::vector<double> &cate_level_weights,
                                    const std::vector<double> &multi_form_weights)
    {
        std::vector<double> dists(row1.size(), NA_REAL);
        double sum_dist = 0.0, nitems = 0.0;
        if (row1.size() != row2.size() || row1.size() == 0)
            return dists;
        std::function<double(const lingdist::StrVec &, const lingdist::StrVec &)> dist_func = [&](const lingdist::StrVec &cells1, const lingdist::StrVec &cells2)
        { return wjd_form(cells1, cells2, cate_level_weights); };
        for (size_t icol = 0; icol < row1.size(); icol++)
        {
            const auto &cells1 = row1[icol];
            const auto &cells2 = row2[icol];
            dists[icol] = dispatcher_weighting(cells1, cells2, dist_func, multi_form_weights);
        }
        return dists;
    }

    double wjd_row(const std::vector<std::vector<lingdist::StrVec>> &row1,
                   const std::vector<std::vector<lingdist::StrVec>> &row2,
                   const std::vector<double> &cate_level_weights,
                   const std::vector<double> &multi_form_weights)
    {
        std::vector<double> dists = wjd_row_vec(row1, row2, cate_level_weights, multi_form_weights);
        return lingdist::nan_mean(dists);
    }
}

namespace lingdist
{
    DataFrame wjd_df(const DataFrame &data, const std::vector<double> &cate_level_weights,
                     const std::vector<double> &multi_form_weights, const String &form_delim, const String &cate_delim, bool detailed, bool squareform, int n_threads, bool quiet)
    {

        auto rows_vector = lingdist::split_df2(data, form_delim, cate_delim);
        auto [lab1Col, lab2Col, row_pairs] = lingdist::get_row_pairs(data, true);

        int n_row_pairs = static_cast<int>(row_pairs.size());

        // RcppThread::ProgressBar *bar = nullptr;
        std::unique_ptr<lingdist::SafeProgressBar> bar;
        if (!quiet)
            bar = std::make_unique<lingdist::SafeProgressBar>(row_pairs.size(), 1);
        std::function<void(std::int32_t)> loop_body;
        std::vector<double> dists(row_pairs.size());
        std::vector<std::vector<double>> dists_vec(row_pairs.size());
        if (detailed)
        {
            loop_body = [&](std::int32_t idx)
            {
                auto [rowi, rowj] = row_pairs[idx];
                dists_vec[idx] = wjd_row_vec(rows_vector[rowi], rows_vector[rowj], cate_level_weights, multi_form_weights);
                if (bar)
                    (*bar)++;
            };
            RcppThread::parallelFor(0, static_cast<std::int32_t>(n_row_pairs), loop_body, n_threads);
            return gen_dist_df_detailed(dists_vec, lab1Col, lab2Col, as<StrVec>(data.names()));
        }
        else
        {
            loop_body = [&](std::int32_t idx)
            {
                auto [rowi, rowj] = row_pairs[idx];
                dists[idx] = wjd_row(rows_vector[rowi], rows_vector[rowj], cate_level_weights, multi_form_weights);
                if (bar)
                    (*bar)++;
            };
            RcppThread::parallelFor(0, static_cast<std::int32_t>(n_row_pairs), loop_body, n_threads);
            return gen_dist_df(dists, lab1Col, lab2Col, squareform, true);
        }
    }
}