#pragma once

#include <Rcpp.h>

using namespace Rcpp;

namespace lingdist
{
    List pmi_df(const DataFrame &data, const String &delim = "", bool squareform = false, bool parallel = false, int n_threads = 4, int max_epochs = 20, double tol = 1e-4, int alignment_max_paths = 3, bool verbose = true);

    DataFrame edit_dist_df(const DataFrame &data, Nullable<DataFrame> cost_mat = R_NilValue, const String &delim = "", bool squareform = false, bool symmetric = true, bool parallel = false, int n_threads = 2);

}