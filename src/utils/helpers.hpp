#pragma once

#include <Rcpp.h>
#include <RcppThread.h>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>

using namespace Rcpp;

namespace lingdist
{
    // --- Common Types & Constants ---

    // Public constants
    extern const std::string EMPTY;
    extern const double LOG2;
    extern const double EPS;

    // Check double equality within EPS
    inline bool double_equal(double a, double b)
    {
        return std::abs(a - b) < EPS;
    }

    // Safe progressbar
    // To suppress warnings about inheriting from class with virtual functions
    // ProgressBar destructor is not virtual
    struct SafeProgressBar final : public RcppThread::ProgressBar
    {
        using RcppThread::ProgressBar::ProgressBar;
    };

    // Single-threaded for loop utility
    void singleFor(int begin, int end, std::function<void(int)> f);

    // Common type aliases
    using Point = std::pair<int, int>;
    using Points = std::vector<Point>;
    using DistPts = std::pair<double, Points>;
    using StrVec = std::vector<std::string>;

    // --- String Utils ---

    // Split a single Rcpp String by a delimiter String into StrVec
    StrVec split(const Rcpp::String &target, const Rcpp::String &delim);

    // Split an entire DataFrame into rows of token vectors
    std::vector<std::vector<StrVec>> split_df(const Rcpp::DataFrame &data, const Rcpp::String &delim);

    std::vector<std::vector<std::vector<StrVec>>> split_df2(const DataFrame &data, const String &delim1, const String &delim2);

    // Extract all unique tokens from a DataFrame
    StrVec get_all_unique_syms(const Rcpp::DataFrame &data, const Rcpp::String &delim, bool include_empty = true);

    // --- DataFrame Utils ---

    // Convert long-form (lab1, lab2, dist) to square DataFrame
    Rcpp::DataFrame long2squareform(const Rcpp::DataFrame &data, bool symmetric, double default_diag = 0.0);

    std::tuple<StringVector, StringVector, std::vector<std::pair<int, int>>> get_row_pairs(const DataFrame &data, bool symmetric);

} // namespace lingdist
