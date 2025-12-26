#pragma once

#include <Rcpp.h>

using namespace Rcpp;

namespace lingdist
{
    List pmi_df(const DataFrame &data, const String &delim = "", bool squareform = false, bool parallel = false, int n_threads = 4, int max_epochs = 20, double tol = 1e-4, int alignment_max_paths = 3, bool verbose = true);

    DataFrame edit_dist_df(const DataFrame &data, CostTable cost, const String &delim = "", bool squareform = false, bool symmetric = true, bool parallel = false, int n_threads = 2, bool check_missing_cost = true);

    DataFrame wjd_df(const DataFrame &data, const std::vector<double> &cate_level_weights,
                     const std::vector<double> &multi_form_weights, const String &form_delim = "#", const String &cate_delim = "_", bool squareform = false, bool parallel = false, int n_threads = 2);

}