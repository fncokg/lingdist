#pragma once

#include <Rcpp.h>
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

    // Extract all unique tokens from a DataFrame
    StrVec get_all_unique_chars(const Rcpp::DataFrame &data, const Rcpp::String &delim);

    // --- DataFrame Utils ---

    // Convert long-form (lab1, lab2, dist) to square DataFrame
    Rcpp::DataFrame long2squareform(const Rcpp::DataFrame &data, bool symmetric);

    std::tuple<StringVector, StringVector, std::vector<std::pair<int, int>>> get_row_pairs(const DataFrame &data, bool symmetric);

} // namespace lingdist
