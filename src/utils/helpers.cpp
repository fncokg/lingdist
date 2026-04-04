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
    const double EPS = 1e-9;

    void singleFor(int begin,
                   int end,
                   std::function<void(int)> f)
    {
        for (int i = begin; i < end; i++)
        {
            f(i);
        }
    }

    std::vector<double> dist_row_pair(const std::vector<std::vector<lingdist::StrVec>> &row1, const std::vector<std::vector<lingdist::StrVec>> &row2, FormsDispatcher dispatcher)
    {
        size_t ncols = row1.size();
        std::vector<double> dists(ncols, NA_REAL);
        for (size_t coli = 0; coli < ncols; coli++)
        {
            const auto &cells1 = row1[coli];
            const auto &cells2 = row2[coli];
            if (!cells1.empty() && !cells2.empty())
            {
                dists[coli] = dispatcher(cells1, cells2);
            }
        }
        return dists;
    }

    // TODO: This is heavy, optimize later
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

    std::vector<std::vector<std::vector<StrVec>>> split_df2(const DataFrame &data, const String &delim1, const String &delim2, bool ignore_delim1)
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
                if (ignore_delim1)
                {
                    cell_vector.push_back(split(col[i], delim2));
                }
                else
                {
                    for (auto &part : split(col[i], delim1))
                    {
                        cell_vector.push_back(split(part, delim2));
                    }
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
        int npairs = symmetric ? nrows * (nrows - 1) / 2 : nrows * nrows;
        CharacterVector rownames = data.attr("row.names");
        StringVector lab1Col(npairs), lab2Col(npairs);

        std::vector<std::pair<int, int>> row_pairs(npairs);
        int outer_loop_end = symmetric ? nrows - 1 : nrows;

        int pair_idx = 0;

        for (int rowi = 0; rowi < outer_loop_end; rowi++)
        {
            int inner_loop_start = symmetric ? rowi + 1 : 0;
            String label1 = rownames[rowi];
            for (int rowj = inner_loop_start; rowj < nrows; rowj++, pair_idx++)
            {

                String label2 = rownames[rowj];
                // row_pairs.emplace_back(rowi, rowj);
                row_pairs[pair_idx] = std::make_pair(rowi, rowj);
                lab1Col[pair_idx] = label1;
                lab2Col[pair_idx] = label2;
            }
        }
        return std::make_tuple(lab1Col, lab2Col, row_pairs);
    }

    StrVec get_all_unique_syms(const DataFrame &data, const String &delim, bool include_empty)
    {
        auto rows_vector = split_df(data, delim);
        StrVec all_syms;
        for (auto &row : rows_vector)
        {
            for (auto &syms : row)
            {
                all_syms.insert(all_syms.end(), syms.begin(), syms.end());
            }
        }
        std::set<std::string> unique_syms(all_syms.begin(), all_syms.end());
        StrVec unique_syms_vec(unique_syms.begin(), unique_syms.end());
        if (include_empty)
        {
            unique_syms_vec.push_back(lingdist::EMPTY);
        }
        return unique_syms_vec;
    }

    DataFrame long2squareform(const DataFrame &data, bool symmetric, double default_diag)
    {
        StringVector lab1_col = data[0];
        StringVector lab2_col = data[1];
        NumericVector dist_col = data[2];

        // get unique labels
        std::vector<String> all_labs(lab1_col.begin(), lab1_col.end());
        all_labs.insert(all_labs.end(), lab2_col.begin(), lab2_col.end());
        std::set<String> unique_labs_set(all_labs.begin(), all_labs.end());
        std::vector<String> unique_labs(unique_labs_set.begin(), unique_labs_set.end());

        // initialize result DataFrame
        StringVector names;
        std::unordered_map<String, int> lab2index;
        DataFrame result;
        for (std::size_t i = 0; i < unique_labs.size(); i++)
        {
            lab2index[unique_labs[i]] = static_cast<int>(i);
            names.push_back(unique_labs[i]);
            result.push_back(NumericVector(unique_labs.size(), NA_REAL));
        }

        // fill in distances
        String char1, char2;
        double dist;
        int nrow = data.nrow();
        for (int i = 0; i < nrow; i++)
        {
            char1 = lab1_col[i];
            char2 = lab2_col[i];
            dist = dist_col[i];
            NumericVector col1 = result[lab2index[char1]];
            col1[lab2index[char2]] = dist;
            if (symmetric)
            {
                NumericVector col2 = result[lab2index[char2]];
                col2[lab2index[char1]] = dist;
            }
        }
        for (size_t i = 0; i < unique_labs.size(); i++)
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

    DataFrame gen_dist_df(std::vector<std::vector<double>> dists_vec, const StringVector &lab1Col, const StringVector &lab2Col, bool detailed, const lingdist::StrVec &col_names, bool squareform, bool symmetric)
    {
        DataFrame result = DataFrame::create(Named("lab1") = lab1Col, Named("lab2") = lab2Col);
        if (detailed)
        {
            size_t n_vars = col_names.size();
            int n_pairs = lab1Col.size();
            for (size_t i = 0; i < n_vars; i++)
            {
                NumericVector dist_vec(n_pairs);
                for (int j = 0; j < n_pairs; j++)
                {
                    dist_vec[j] = dists_vec[j][i];
                }
                result.push_back(dist_vec, col_names[i]);
            }
        }
        else
        {
            std::vector<double> dists(dists_vec.size());
            for (size_t i = 0; i < dists_vec.size(); i++)
            {
                dists[i] = lingdist::nan_mean(dists_vec[i]);
            }
            NumericVector dist_col = NumericVector::import(dists.begin(), dists.end());
            result.push_back(dist_col, "dist");
            if (squareform)
            {
                result = lingdist::long2squareform(result, symmetric);
            }
        }
        return result;
    }

    double nan_mean(const std::vector<double> &vec)
    {
        double sum = 0.0;
        int count = 0;
        for (double val : vec)
        {
            if (!Rcpp::NumericVector::is_na(val))
            {
                sum += val;
                count++;
            }
        }
        if (count == 0)
            return NA_REAL;
        return sum / count;
    }

} // namespace lingdist
