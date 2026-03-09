#pragma once

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

namespace lingdist
{
    List pmi_df(const DataFrame &data, const String &delim = "", bool detailed = false, const String &normalize_method = "longest", bool squareform = false, int n_threads = 2, int max_epochs = 20, double tol = 1e-4, int alignment_max_paths = 3, bool quiet = false);

    DataFrame edit_dist_df(const DataFrame &data, const CostTable &cost, const String &delim = "", bool detailed = false, const String &normalize_method = "longest", const String &form_strategy = "off", const String &form_delim = "#", std::vector<double> form_weights = std::vector<double>(), bool squareform = false, bool symmetric = true, int n_threads = 2, bool check_missing_cost = true, bool quiet = false);

    DataFrame wjd_df(const DataFrame &data, const String &cate_delim = "_", const std::vector<double> &cate_weights = std::vector<double>(), const String &form_strategy = "weighting", const String &form_delim = "#", const std::vector<double> &form_weights = std::vector<double>(), bool detailed = false, bool squareform = false, int n_threads = 2, bool quiet = false);

}