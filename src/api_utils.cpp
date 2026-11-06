#include <Rcpp.h>

#include "helpers.hpp"
#include "cost_io.hpp"
#include "prepare_input.hpp"
#include "text.hpp"

using namespace Rcpp;

//' Generate a default cost matrix
//'
//' Generate a default cost matrix containing all possible characters in the raw data with all diagonal values set to 0 and others set to 1. This avoids constructing the matrix from scratch.
//'
//' @param data DataFrame to be computed.
//' @param delim Literal delimiter separating atomic symbols. Empty splits UTF-8 text into Unicode code points.
//' @param form_delim Literal delimiter separating forms within a cell, default "#". Empty forms are discarded.
//' @return Cost matrix containing all possible characters in the raw data with all diagonal values set to 0 and others set to 1.
//' @examples
//' df <- as.data.frame(rbind(a=c("ph_l_i_z","k_o_l"),b=c("b_l_i_s", "k_a:_l")))
//' default.cost <- generate_default_cost_matrix(df, "_")
//[[Rcpp::export]]
DataFrame generate_default_cost_matrix(const DataFrame &data, const String &delim = "", const String &form_delim = "#")
{
    const auto input = lingdist::prepare_input(data, delim, form_delim);
    return lingdist::write_costs(lingdist::IndexedCosts::uniform(input.symbol_table()));
}

//' Convert long table to square form
//'
//' Convert a distance dataframe in long table form to a square matrix form. Values in the diagonal positions are automatically filled with `default_diag`. If you want self-defined diagonal values, define them in the long table.
//'
//' @param data Dataframe in long table form. The first and second columns are labels and the third column stores the distance values.
//' @param symmetric Whether the distance matrix is symmetric (if cost matrix is not, then the distance matrix is also not).
//' @param default_diag The default value to fill in the diagonal positions. By default it is 0.0.
//' @return Dataframe in square matrix form, rownames and colnames are labels. If the long table only contains \eqn{C_n^2} rows and `symmetric` is set to FALSE, then only lower triangle positions in the result are filled.
//' @examples
//' data <- as.data.frame(list(chars1=c("a","a","b"),chars2=c("b","c","c"),dist=c(1,2,3)))
//' mat <- long2squareform(data)
//[[Rcpp::export]]
DataFrame long2squareform(const DataFrame &data, bool symmetric = true, double default_diag = 0.0)
{
    return lingdist::long2squareform(data, symmetric, default_diag);
}

//' Check which symbols in `data` are missing from a cost matrix
//'
//' Given a cost matrix (square-form DataFrame) and a data table of tokenized strings,
//' return a character vector of symbols that appear in `data` but are not present in the
//' cost matrix row/column names. This is useful to validate a user-supplied cost matrix
//' before running distance computations.
//'
//' @param cost_mat DataFrame representing a cost matrix in square form. Row names and
//'   column names should be the symbols included in the matrix.
//' @param data DataFrame containing the tokenized strings; used to extract symbols to check.
//' @param delim Literal delimiter used in `data` to split atomic symbols
//'   (empty string means split UTF-8 text into Unicode code points).
//' @param form_delim Literal delimiter separating forms within a cell, default "#". Empty forms are discarded.
//' @return A character vector of symbols that are present in `data` but missing from `cost_mat`.
//' @examples
//' cost.mat <- data.frame()
//' df <- as.data.frame(rbind(a=c("ph_l_i_z","k_o_l"), b=c("b_l_i_s","k_a:_l")))
//' missing <- check_cost_mat_symbols(cost.mat, df, "_")
//[[Rcpp::export]]
StringVector check_cost_mat_symbols(const DataFrame &cost_mat, const DataFrame &data, const String &delim, const String &form_delim = "#")
{
    const auto input = lingdist::prepare_input(data, delim, form_delim);
    const auto missing_symbols = lingdist::missing_cost_symbols(input.symbol_table(), cost_mat);
    return lingdist::wrap_utf8(missing_symbols);
}
