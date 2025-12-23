#include <Rcpp.h>
#include <RcppThread.h>

#include "helpers.hpp"

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
        for (size_t i = 0; i < min_len; i++)
        {
            if (form1[i] != form2[i])
            {
                return cate_level_weights[i];
            }
        }
        // ["A2"] vs ["A2","a"]
        return cate_level_weights[min_len - 1];
    }

    double wjd_row(const std::vector<std::vector<lingdist::StrVec>> &row1,
                   const std::vector<std::vector<lingdist::StrVec>> &row2,
                   const std::vector<double> &cate_level_weights,
                   const std::vector<double> &multi_form_weights)
    {
        double sum_dist = 0.0, nitems = 0.0;
        if (row1.size() != row2.size() || row1.size() == 0)
            return 0.0;

        for (size_t icol = 0; icol < row1.size(); icol++)
        {
            const auto &cells1 = row1[icol];
            const auto &cells2 = row2[icol];
            size_t nforms1 = cells1.size(), nforms2 = cells2.size();
            if (nforms1 != 0 && nforms2 != 0)
            {
                double this_dist = 0.0, sum_weights = 0.0;
                size_t max_len = std::max(nforms1, nforms2);
                for (size_t i = 0; i < max_len; i++)
                {
                    const auto &form1 = i < nforms1 ? cells1[i] : cells1.back();
                    const auto &form2 = i < nforms2 ? cells2[i] : cells2.back();
                    this_dist += wjd_form(form1, form2, cate_level_weights) * multi_form_weights[i];
                    // only the first max_len weights are counted
                    sum_weights += multi_form_weights[i];
                }
                if (sum_weights > 0.0)
                {
                    this_dist /= sum_weights;
                }
                sum_dist += this_dist;
                nitems += 1.0;
            }
        }
        return nitems > 0.0 ? sum_dist / nitems : 0.0;
    }
}

namespace lingdist
{
    DataFrame wjd_df(const DataFrame &data, const std::vector<double> &cate_level_weights,
                     const std::vector<double> &multi_form_weights, const String &form_delim, const String &cate_delim, bool squareform, bool parallel, int n_threads)
    {

        auto rows_vector = lingdist::split_df2(data, form_delim, cate_delim);
        auto [lab1Col, lab2Col, row_pairs] = lingdist::get_row_pairs(data, true);

        int n_row_pairs = static_cast<int>(row_pairs.size());

        std::vector<double> dists(row_pairs.size());
        RcppThread::ProgressBar bar(row_pairs.size(), 1);
        std::function<void(std::int32_t)> loop_body = [&](std::int32_t idx)
        {
            auto [rowi, rowj] = row_pairs[idx];
            dists[idx] = wjd_row(rows_vector[rowi], rows_vector[rowj], cate_level_weights, multi_form_weights);
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
            result = lingdist::long2squareform(result, true);
        }
        return result;
    }
}