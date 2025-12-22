#include <Rcpp.h>
#include <math.h>

#include "alignment.hpp"
#include "dist.hpp"

using namespace Rcpp;

namespace
{
    inline std::vector<double> make_default_weights(int n)
    {
        std::vector<double> weights;
        for (int i = 0; i < n; ++i)
        {
            weights.push_back(std::pow(0.5, i));
        }
        return weights;
    }
}

//' Compute edit distance between two strings
//'
//' Compute edit distance between two strings and get all possible alignment scenarios. Custom cost matrix is supported. Symbols separated by custom delimiters are supported.
//'
//' @param str1 String to be compared.
//' @param str2 String to be compared.
//' @param cost_mat Dataframe in squareform indicating the cost values when one symbol is deleted, inserted or substituted by another. Rownames and colnames are symbols. `cost_mat[char1,"_NULL_"]` indicates the cost value of deleting char1 and `cost_mat["_NULL_",char1]` is the cost value of inserting it. When an operation is not defined in the cost_mat, it is set 0 when the two symbols are the same, otherwise 1.
//' @param delim The delimiter in `str1` and `str2` separating atomic symbols.
//' @param return_alignments Whether to return alignment details
//' @return A list containing a `distance` element storing the distance result. If `return_alignments` is TRUE, then an `alignments` element is present which is a list of dataframes with each storing a possible best alignment scenario.
//' @examples
//' cost.mat <- data.frame()
//' dist <- edit_dist_string("leaf","leaves")$distance
//' dist <- edit_dist_string("ph_l_i_z","p_l_i_s",cost_mat=cost.mat,delim="_")$distance
//' alignments <- edit_dist_string("ph_l_i_z","p_l_i_s",delim="_",return_alignments=TRUE)$alignments
//[[Rcpp::export]]
List string_edit_dist(const String &str1, const String &str2, Nullable<DataFrame> cost_mat = R_NilValue, const String &delim = "", bool return_alignments = false)
{
    lingdist::CostTable cost;
    if (cost_mat.isNotNull())
    {
        DataFrame cost_mat_ = DataFrame(cost_mat);
        cost = lingdist::build_cost_table(cost_mat_);
    }
    lingdist::StrVec chars1 = lingdist::split(str1, delim);
    lingdist::StrVec chars2 = lingdist::split(str2, delim);
    if (return_alignments)
    {
        return lingdist::get_string_alignment(chars1, chars2, cost);
    }
    else
    {
        double result = lingdist::edit_dist_core_dp(chars1, chars2, cost);
        List report;
        report["distance"] = result;
        return report;
    }
}

//' Compute edit distance between all row pairs of a dataframe
//'
//' Compute average edit distance between all row pairs of a dataframe, empty or NA cells are ignored. If all values in a row are not valid strings, all average distances involving this row is set to -1.
//'
//' @param data DataFrame with n rows and m columns indicating there are n languages or dialects to involve in the calculation and there are at most m words to base on, in which the rownames are the language ids.
//' @param cost_mat Dataframe in squareform indicating the cost values when one symbol is deleted, inserted or substituted by another. Rownames and colnames are symbols. `cost_mat[char1,"_NULL_"]` indicates the cost value of deleting char1 and `cost_mat["_NULL_",char1]` is the cost value of inserting it. When an operation is not defined in the cost_mat, it is set 0 when the two symbols are the same, otherwise 1.
//' @param delim The delimiter separating atomic symbols.
//' @param squareform Whether to return a dataframe in squareform.
//' @param symmetric Whether the result matrix is symmetric. This depends on whether the `cost_mat` is symmetric.
//' @param parallel Whether to parallelize the computation.
//' @param n_threads The number of threads is used to parallelize the computation. Only meaningful if `parallel` is TRUE.
//' @return A dataframe in long table form if `squareform` is FALSE, otherwise in squareform. If `symmetric` is TRUE, the long table form has \eqn{C_n^2} rows, otherwise \eqn{n^2} rows.
//' @examples
//' df <- as.data.frame(rbind(a=c("a_bc_d","d_bc_a"),b=c("b_bc_d","d_bc_a")))
//' cost.mat <- data.frame()
//' result <- edit_dist_df(df, cost_mat=cost.mat, delim="_")
//' result <- edit_dist_df(df, cost_mat=cost.mat, delim="_", squareform=TRUE)
//' result <- edit_dist_df(df, cost_mat=cost.mat, delim="_", parallel=TRUE, n_threads=4)
//[[Rcpp::export]]
DataFrame pw_edit_dist(const DataFrame &data, Nullable<DataFrame> cost_mat = R_NilValue, const String &delim = "", bool squareform = false, bool symmetric = true, bool parallel = false, int n_threads = 2)
{
    return lingdist::edit_dist_df(data, cost_mat, delim, squareform, symmetric, parallel, n_threads);
}

//' Compute PMI distance between all row pairs of a dataframe
//'
//' Compute PMI (Pointwise Mutual Information) distance between all row pairs of a dataframe.
//'
//' @param data DataFrame with n rows and m columns indicating there are n languages or dialects to involve in the calculation and there are at most m words to base on, in which the rownames are the language ids.
//' @param delim The delimiter separating atomic symbols.
//' @param squareform Whether to return a dataframe in squareform.
//' @param parallel Whether to parallelize the computation.
//' @param n_threads The number of threads is used to parallelize the computation. Only meaningful if `parallel` is TRUE.
//' @param max_epochs Maximum number of epochs for EM algorithm.
//' @param tol Tolerance for convergence.
//' @param alignment_max_paths Maximum number of paths to consider in alignment.
//' @param verbose Whether to print progress.
//' @return A list containing the following components:
//' \item{result}{A dataframe of distances, either in long table form or square form.}
//' \item{cost}{The final cost matrix used for distance calculation.}
//' \item{prev_cost}{The cost matrix from the previous iteration.}
//' @examples
//' df <- as.data.frame(rbind(a=c("a_bc_d","d_bc_a"),b=c("b_bc_d","d_bc_a")))
//' result <- pw_pmi_dist(df, delim="_")
//' result <- pw_pmi_dist(df, delim="_", squareform=TRUE)
//' result <- pw_pmi_dist(df, delim="_", parallel=TRUE, n_threads=4)
//[[Rcpp::export]]
List pw_pmi_dist(const DataFrame &data, const String &delim = "", bool squareform = false, bool parallel = false, int n_threads = 4, int max_epochs = 20, double tol = 1e-4, int alignment_max_paths = 3, bool verbose = true)
{
    return lingdist::pmi_df(data, delim, squareform, parallel, n_threads, max_epochs, tol, alignment_max_paths, verbose);
}

//' Compute Weighted Jaccard Distance between all row pairs of a dataframe
//'
//' Compute Weighted Jaccard Distance (WJD) between all row pairs of a dataframe. This metric is suitable for categorical data with hierarchical structure and multiple forms.
//' For example, a cell value "A_2_a#A_3#B_5_c" indicates there are 3 word forms for this lexical item. The first form "A_2_a" belongs to category A, subcategory 2, and sub-subcategory a. The second form "A_3" belongs to category A and subcategory 3. The third form "B_5_c" belongs to category B, subcategory 5, and sub-subcategory c. The delimiters can be customized via `form_delim` and `cate_delim`.
//'
//' @param data DataFrame with n rows and m columns indicating there are n languages or dialects to involve in the calculation and there are at most m words to base on, in which the rownames are the language ids.
//' @param cate_level_weights Numeric vector of weights for different levels of categories for a single word form. For example, if a word form is "A_2_a", the first weight applies to "A", the second to "2", and the third to "a". If not provided, default weights are used.
//' @param multi_form_weights Numeric vector of weights for different word forms of a single lexical item (ordered). For example, if a lexical item is "A_2_a#A_3#B_5_c", the first weight applies to "A_2_a", the second to "A_3", etc. If not provided, default weights are used.
//' @param form_delim The delimiter separating different word forms of a single lexical item. Default is "#".
//' @param cate_delim The delimiter separating different levels of categories within a word form. Default is "_".
//' @param squareform Whether to return a dataframe in squareform.
//' @param parallel Whether to parallelize the computation.
//' @param n_threads The number of threads is used to parallelize the computation. Only meaningful if `parallel` is TRUE.
//' @return A dataframe in long table form if `squareform` is FALSE, otherwise in squareform.
//' @examples
//' df <- as.data.frame(rbind(a=c("A_1#B_2","C_3"),b=c("A_1","C_3_x")))
//' result <- pw_wjd(df, form_delim="#", cate_delim="_")
//[[Rcpp::export]]
DataFrame pw_wjd(const DataFrame &data, Nullable<NumericVector> cate_level_weights = R_NilValue, Nullable<NumericVector> multi_form_weights = R_NilValue, const String &form_delim = "#", const String &cate_delim = "_", bool squareform = false, bool parallel = false, int n_threads = 2)
{
    std::vector<double> cate_weights, form_weights;
    cate_weights = cate_level_weights.isNotNull() ? as<std::vector<double>>(NumericVector(cate_level_weights)) : make_default_weights(10);
    form_weights = multi_form_weights.isNotNull() ? as<std::vector<double>>(NumericVector(multi_form_weights)) : make_default_weights(10);
    return lingdist::wjd_df(data, cate_weights, form_weights, form_delim, cate_delim, squareform, parallel, n_threads);
}