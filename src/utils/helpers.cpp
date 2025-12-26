#include "helpers.hpp"
#include <set>
#include <unordered_map>
#include <vector>
#include <string>

using namespace Rcpp;

namespace lingdist
{
    const std::string EMPTY = "_NULL_";
    const double LOG2 = std::log(2.0);

    void singleFor(int begin,
                   int end,
                   std::function<void(int)> f)
    {
        for (int i = begin; i < end; i++)
        {
            f(i);
        }
    }

    StrVec split(const String &target, const String &delim)
    {
        Function r_split("strsplit");
        List res = r_split(target, delim);
        return res[0];
    }

    std::vector<std::vector<StrVec>> split_df(const DataFrame &data, const String &delim)
    {
        int nrows = data.nrow(), ncols = data.ncol();
        std::vector<std::vector<StrVec>> rows_vector;
        rows_vector.reserve(nrows);
        for (int i = 0; i < nrows; i++)
        {
            std::vector<StrVec> row_vector;
            row_vector.reserve(ncols);
            for (int j = 0; j < ncols; j++)
            {
                const StrVec &col = data[j];
                row_vector.push_back(split(col[i], delim));
            }
            rows_vector.push_back(std::move(row_vector));
        }
        return rows_vector;
    }

    std::vector<std::vector<std::vector<StrVec>>> split_df2(const DataFrame &data, const String &delim1, const String &delim2)
    {
        int nrows = data.nrow(), ncols = data.ncol();
        std::vector<std::vector<std::vector<StrVec>>> rows_vector;
        rows_vector.reserve(nrows);
        for (int i = 0; i < nrows; i++)
        {
            std::vector<std::vector<StrVec>> row_vector;
            row_vector.reserve(ncols);
            for (int j = 0; j < ncols; j++)
            {
                const StrVec &col = data[j];
                std::vector<StrVec> cell_vector;
                for (auto &part : split(col[i], delim1))
                {
                    cell_vector.push_back(split(part, delim2));
                }
                row_vector.push_back(std::move(cell_vector));
            }
            rows_vector.push_back(std::move(row_vector));
        }
        return rows_vector;
    }

    std::tuple<StringVector, StringVector, std::vector<std::pair<int, int>>> get_row_pairs(const DataFrame &data, bool symmetric)
    {
        int nrows = data.nrow();
        CharacterVector rownames = data.attr("row.names");
        StringVector lab1Col, lab2Col;

        std::vector<std::pair<int, int>> row_pairs;
        int inner_loop_start;
        for (int rowi = 0; rowi < nrows - 1; rowi++)
        {
            inner_loop_start = symmetric ? rowi + 1 : 0;
            for (int rowj = inner_loop_start; rowj < nrows; rowj++)
            {
                String label1 = rownames[rowi];
                String label2 = rownames[rowj];
                row_pairs.emplace_back(rowi, rowj);
                lab1Col.push_back(label1);
                lab2Col.push_back(label2);
            }
        }
        return std::make_tuple(lab1Col, lab2Col, row_pairs);
    }

    StrVec get_all_unique_chars(const DataFrame &data, const String &delim)
    {
        auto rows_vector = split_df(data, delim);
        StrVec all_chars;
        for (auto &row : rows_vector)
        {
            for (auto &chars : row)
            {
                all_chars.insert(all_chars.end(), chars.begin(), chars.end());
            }
        }
        std::set<std::string> unique_chars(all_chars.begin(), all_chars.end());
        return StrVec(unique_chars.begin(), unique_chars.end());
    }

    DataFrame long2squareform(const DataFrame &data, bool symmetric, double default_diag)
    {
        StringVector chars1_col = data[0];
        StringVector chars2_col = data[1];
        NumericVector dist_col = data[2];

        String char1, char2;
        double dist;
        std::vector<String> all_chars(chars1_col.begin(), chars1_col.end());
        all_chars.insert(all_chars.end(), chars2_col.begin(), chars2_col.end());
        std::set<String> unique_chars_set(all_chars.begin(), all_chars.end());
        std::vector<String> unique_chars(unique_chars_set.begin(), unique_chars_set.end());
        StringVector names;
        std::unordered_map<String, int> char2idx;
        DataFrame result;
        for (std::size_t i = 0; i < unique_chars.size(); i++)
        {
            char2idx[unique_chars[i]] = static_cast<int>(i);
            names.push_back(unique_chars[i]);
            result.push_back(NumericVector(unique_chars.size(), NA_REAL));
        }
        int nrow = data.nrow();
        for (int i = 0; i < nrow; i++)
        {
            char1 = chars1_col[i];
            char2 = chars2_col[i];
            dist = dist_col[i];
            NumericVector col1 = result[char2idx[char1]];
            col1[char2idx[char2]] = dist;
            if (symmetric)
            {
                NumericVector col2 = result[char2idx[char2]];
                col2[char2idx[char1]] = dist;
            }
        }
        for (size_t i = 0; i < unique_chars.size(); i++)
        {
            NumericVector col = result[i];
            if (NumericVector::is_na(col[i]))
            {
                col[i] = default_diag;
            }
        }
        result.attr("row.names") = names;
        result.attr("names") = names;
        return result;
    }

} // namespace lingdist
