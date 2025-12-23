#include <Rcpp.h>

#include "helpers.hpp"
#include "cost_table.hpp"
#include "alignment.hpp"

using namespace Rcpp;

//' Generate a default cost matrix
//'
//' Generate a default cost matrix containing all possible characters in the raw data with all diagonal values set to 0 and others set to 1. This avoids constructing the matrix from scratch.
//'
//' @param data DataFrame to be computed.
//' @param delim The delimiter separating atomic symbols.
//' @return Cost matrix containing all possible characters in the raw data with all diagonal values set to 0 and others set to 1.
//' @examples
//' df <- as.data.frame(rbind(a=c("ph_l_i_z","k_o_l"),b=c("b_l_i_s", "k_a:_l")))
//' default.cost <- generate_default_cost_matrix(df, "_")
//[[Rcpp::export]]
DataFrame generate_default_cost_matrix(const DataFrame &data, const String &delim = "")
{
    lingdist::StrVec chars = lingdist::get_all_unique_chars(data, delim);
    return lingdist::build_default_cost_table(chars).to_dataframe();
}

//' Convert long table to square form
//'
//' Convert a distance dataframe in long table form to a square matrix form.
//'
//' @param data Dataframe in long table form. The first and second columns are labels and the third column stores the distance values.
//' @param symmetric Whether the distance matrix is symmetric (if cost matrix is not, then the distance matrix is also not).
//' @return Dataframe in square matrix form, rownames and colnames are labels. If the long table only contains \eqn{C_n^2} rows and `symmetric` is set to FALSE, then only lower triangle positions in the result are filled.
//' @examples
//' data <- as.data.frame(list(chars1=c("a","a","b"),chars2=c("b","c","c"),dist=c(1,2,3)))
//' mat <- long2squareform(data)
//[[Rcpp::export]]
DataFrame long2squareform(const DataFrame &data, bool symmetric = true)
{
    return lingdist::long2squareform(data, symmetric);
}