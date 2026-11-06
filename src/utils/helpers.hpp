#pragma once

#include <Rcpp.h>
#include <RcppThread.h>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <functional>
#include <numeric>
#include "numeric.hpp"

using namespace Rcpp;

namespace lingdist
{
    // --- Common Types & Constants ---

    // Public constants
    extern const std::string EMPTY;
    extern const double LOG2;

    // Safe progressbar
    // To suppress warnings about inheriting from class with virtual functions
    // ProgressBar destructor is not virtual
    struct SafeProgressBar final : public RcppThread::ProgressBar
    {
        using RcppThread::ProgressBar::ProgressBar;
    };

    // Common type aliases
    using StrVec = std::vector<std::string>;
    using FormsDispatcher = std::function<double(const std::vector<lingdist::StrVec> &, const std::vector<lingdist::StrVec> &)>;
    using DistFunc = std::function<double(const lingdist::StrVec &, const lingdist::StrVec &)>;

    // loop utility
    void singleFor(int begin, int end, std::function<void(int)> f);

    std::vector<double> dist_row_pair(const std::vector<std::vector<lingdist::StrVec>> &row1, const std::vector<std::vector<lingdist::StrVec>> &row2, FormsDispatcher dispatcher);

    // --- String Utils ---

    // Convert R text to UTF-8 and split by a literal delimiter; empty means code points.
    StrVec split(const Rcpp::String &target, const Rcpp::String &delim);

    std::vector<std::vector<std::vector<StrVec>>> split_df2(const DataFrame &data, const String &delim1, const String &delim2, bool ignore_delim1);

    // --- DataFrame Utils ---

    // Convert long-form (lab1, lab2, dist) to square DataFrame
    Rcpp::DataFrame long2squareform(const Rcpp::DataFrame &data, bool symmetric, double default_diag = 0.0);

    std::tuple<StringVector, StringVector, std::vector<std::pair<int, int>>> get_row_pairs(const DataFrame &data, bool symmetric);

    // Generate distance DataFrame from distance vectors and labels
    DataFrame gen_dist_df(std::vector<std::vector<double>> dists_vec, const StringVector &lab1Col, const StringVector &lab2Col, bool detailed, const lingdist::StrVec &col_names, bool squareform, bool symmetric);

    // Validate/generate weights before dispatching work to parallel threads.
    std::vector<double> prepare_weights(const std::vector<std::vector<std::vector<StrVec>>> &rows,
                                        std::vector<double> weights, bool categories, const char *name);

    // --- Vector Utils ---
    double nan_mean(const std::vector<double> &vec);

} // namespace lingdist
