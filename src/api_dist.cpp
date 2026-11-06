#include <Rcpp.h>
#include <math.h>

#include "alignment.hpp"
#include "cost_io.hpp"
#include "dist.hpp"

using namespace Rcpp;

namespace
{
    std::vector<double> read_weights(Nullable<NumericVector> weights, const char *name)
    {
        // NULL requests adaptive defaults; an explicitly empty vector is invalid (B01/B02).
        if (weights.isNull())
            return {};
        auto values = as<std::vector<double>>(NumericVector(weights));
        if (values.empty())
            stop("%s must not be empty; use NULL for default weights.", name);
        return values;
    }
}

//' Compute edit distance between two strings
//'
//' Compute edit distance between two strings and optionally display up to 10 optimal alignments. Custom cost matrix is supported. Symbols separated by custom delimiters are supported.
//'
//' @param str1 String to be compared.
//' @param str2 String to be compared.
//' @param cost_mat Dataframe in squareform indicating the cost values when one symbol is deleted, inserted or substituted by another. Rownames and colnames are symbols. `cost_mat[char1,"_NULL_"]` indicates the cost value of deleting char1 and `cost_mat["_NULL_",char1]` is the cost value of inserting it. When an operation is not defined in the cost_mat, it is set 0 when the two symbols are the same, otherwise 1. When `cost_mat` is NULL, general cost rules are used: substitution cost is `default_sub_cost` if the two symbols are different and 0 if they are the same; insertion and deletion cost is `default_ins_del_cost`.
//' @param delim Literal delimiter in `str1` and `str2` separating atomic symbols. Empty splits UTF-8 text into Unicode code points, not grapheme clusters.
//' @param normalize_method Method to normalize the distance. Can be "longest" (default) to normalize by the length of the longer string, or "none" to return the raw distance without normalization.
//' @param return_alignments Whether to return alignment details for up to 10 paths. This display limit does not apply to PMI training.
//' @param default_sub_cost Default substitution cost when `cost_mat` is NULL.
//' @param default_ins_del_cost Default insertion and deletion cost when `cost_mat` is NULL.
//' @param quiet Whether to suppress all output messages.
//' @return A list containing a `distance` element storing the distance result. If `return_alignments` is TRUE, then an `alignments` element is present which is a list of dataframes with each storing a possible best alignment scenario.
//' @seealso [generate_default_cost_matrix()]
//' @examples
//' cost.mat <- data.frame()
//' dist <- string_edit_dist("leaf","leaves")$distance
//' res <- string_edit_dist("ph_l_i_z","p_l_i_s",cost_mat=cost.mat,delim="_")
//' alignments <- res$alignments
//[[Rcpp::export]]
List string_edit_dist(const String &str1, const String &str2, Nullable<DataFrame> cost_mat = R_NilValue, const String &delim = "", const String &normalize_method = "longest", bool return_alignments = false, double default_sub_cost = 1.0, double default_ins_del_cost = 1.0, bool quiet = false)
{
    const auto chars1 = lingdist::split(str1, delim);
    const auto chars2 = lingdist::split(str2, delim);
    auto alphabet = chars1;
    alphabet.insert(alphabet.end(), chars2.begin(), chars2.end());
    const lingdist::SymbolTable symbols(alphabet);
    const auto cost = cost_mat.isNull()
                          ? lingdist::IndexedCosts::uniform(symbols, default_sub_cost, default_ins_del_cost)
                          : lingdist::read_costs(symbols, DataFrame(cost_mat));
    const auto source = symbols.encode(chars1), target = symbols.encode(chars2);
    const std::string method = normalize_method;
    if (return_alignments)
    {
        return lingdist::get_string_alignment(source, target, cost, 10, method);
    }
    else
    {
        double result = lingdist::ed_dp(source, target, cost, method);
        List report;
        report["distance"] = result;
        return report;
    }
}

//' Compute edit distance between all row pairs of a dataframe
//'
//' Compute average edit distance between all row pairs of a dataframe, empty or NA cells are ignored. If a row pair has no jointly observed cells, its average distance is NA.
//'
//' @param data DataFrame with n rows and m columns indicating there are n languages or dialects to involve in the calculation and there are at most m words to base on, in which the rownames are the language ids.
//' @param cost_mat Dataframe in squareform indicating the cost values when one symbol is deleted, inserted or substituted by another. Rownames and colnames are symbols. `cost_mat[char1,"_NULL_"]` indicates the cost value of deleting char1 and `cost_mat["_NULL_",char1]` is the cost value of inserting it. When an operation is not defined in the cost_mat, it is set 0 when the two symbols are the same, otherwise 1. When `cost_mat` is NULL, general cost rules are used: substitution cost is `default_sub_cost` if the two symbols are different and 0 if they are the same; insertion and deletion cost is `default_ins_del_cost`.
//' @param delim Literal delimiter separating atomic symbols. Empty splits UTF-8 text into Unicode code points.
//' @param detailed Whether to return detailed information. When FALSE (default), returns the average distance across all columns (words) for each row pair. When TRUE, returns column-wise distances for each row pair, showing how distance varies across words.
//' @param normalize_method Method to normalize the distance. Can be "longest" (default) to normalize by the length of the longer string, or "none" to return the raw distance without normalization.
//' @param squareform Whether to return a dataframe in squareform.
//' @param symmetric Whether the result matrix is symmetric. This depends on whether the `cost_mat` is symmetric.
//' @param n_threads Number of threads used for computation. Use 1 for sequential computation.
//' @param form_strategy How to combine word forms: off (default), weighting, mean_all, min_all, or first.
//' @param form_delim Literal delimiter separating word forms when form_strategy is not off.
//' @param form_weights Finite, non-negative weights covering all observed forms for weighting. NULL generates weights 1, 0.5, 0.25, ... to the required length. An all-zero effective weight sum returns NA.
//' @param check_missing_cost Whether to check if all symbols in `data` are defined in `cost_mat`. If TRUE, an warning message is printed when there are missing symbols. You are recommended to set this to `TRUE` unless you are sure all symbols are defined in `cost_mat` and you want to skip the check for performance consideration.
//' @param default_sub_cost Default substitution cost when `cost_mat` is NULL.
//' @param default_ins_del_cost Default insertion and deletion cost when `cost_mat` is NULL.
//' @param quiet Whether to suppress all output messages.
//' @param cache Whether to precompute and reuse distances for unique directed word pairs within this call. Default FALSE skips caching. TRUE can speed up data with repeated word forms, at the cost of task planning and extra memory; small or low-repetition data may be faster with FALSE.
//' @return A dataframe in long table form if `squareform` is FALSE, otherwise in squareform. If `symmetric` is TRUE, the long table form has \eqn{C_n^2} rows, otherwise \eqn{n^2} rows.
//' @examples
//' df <- as.data.frame(rbind(a=c("ph_l_i_z","k_o_l"),b=c("b_l_i_s", "k_a:_l")))
//' df2 <- as.data.frame(rbind(a=c("ph_l_i_z#b_l_i_s","k_o_l"),b=c("b_l_i_s", "k_a:_l#k_o")))
//' cost.mat <- data.frame()
//' result <- pw_edit_dist(df, delim="_")
//' result <- pw_edit_dist(df, delim="_", cache=TRUE)
//' result <- pw_edit_dist(df, delim="_", default_sub_cost=2.0, default_ins_del_cost=1.5)
//' result <- pw_edit_dist(df, cost_mat=cost.mat, delim="_")
//' result <- pw_edit_dist(df, cost_mat=cost.mat, delim="_", squareform=TRUE)
//' result <- pw_edit_dist(df2, cost_mat=cost.mat, delim="_", form_strategy="weighting", form_delim="#", form_weights=c(0.7, 0.3))
//[[Rcpp::export]]
DataFrame pw_edit_dist(const DataFrame &data, Nullable<DataFrame> cost_mat = R_NilValue, const String &delim = "", bool detailed = false, const String &normalize_method = "longest", const String &form_strategy = "off", const String &form_delim = "#", Nullable<NumericVector> form_weights = R_NilValue, bool squareform = false, bool symmetric = true, int n_threads = 2, bool check_missing_cost = true, double default_sub_cost = 1.0, double default_ins_del_cost = 1.0, bool quiet = false, bool cache = false)
{
    const auto form_weights_vec = read_weights(form_weights, "form_weights");
    return lingdist::edit_dist_df(data, cost_mat, delim, detailed, normalize_method, form_strategy,
                                  form_delim, form_weights_vec, squareform, symmetric, n_threads,
                                  check_missing_cost, quiet, default_sub_cost, default_ins_del_cost, cache);
}

//' Compute PMI distance between all row pairs of a dataframe
//'
//' Compute PMI (Pointwise Mutual Information) distance between all row pairs of a dataframe.
//'
//' Training uses the expected nonmatching token-pair occurrences over all paths
//' in the optimal alignment DAG, with equal probability for each complete path.
//' A DP predecessor is retained when its candidate cost differs from the local
//' minimum by less than 1e-9. Matching edges affect path probabilities but are
//' excluded from PMI observations. Form weights are applied separately; there
//' is no training path limit or length weighting.
//'
//' @param data DataFrame with n rows and m columns indicating there are n languages or dialects to involve in the calculation and there are at most m words to base on, in which the rownames are the language ids.
//' @param delim Literal delimiter separating atomic symbols. Empty splits UTF-8 text into Unicode code points.
//' @param detailed Whether to return detailed information. When FALSE (default), returns the average distance across all columns (words) for each row pair. When TRUE, returns column-wise distances for each row pair, showing how distance varies across words.
//' @param normalize_method Method to normalize the distance. Can be "longest" (default) to normalize by the length of the longer string, or "none" to return the raw distance without normalization.
//' @param squareform Whether to return a dataframe in squareform.
//' @param n_threads Positive number of threads used for computation. Use 1 for sequential computation.
//' @param max_epochs Maximum number of epochs; must be at least 1.
//' @param tol Finite positive tolerance for convergence.
//' @param damping_alpha Finite update weight in [0, 1]. Each new cost is (1 - damping_alpha) times the previous cost plus damping_alpha times the mapped PMI target. Default 1 applies the full update. Convergence uses the damped change; 0 keeps the initial costs and converges after one iteration.
//' @param smoothing_lambda Finite, strictly positive pseudocount for each unordered nonmatching token pair, including character-gap pairs. Applied once after merging all observations. Diagonal pairs receive no prior; with V symbols including the gap, each marginal receives (V - 1) * smoothing_lambda. Default 0.001.
//' @param min_edit_cost Finite lower bound in (0, 1] for nonmatching costs, including insertions and deletions. Nonconstant negative PMI scores are mapped to [min_edit_cost, 1] before damping; equal scores map to 1. Matching costs remain 0. Default 0.01.
//' @param form_strategy How to pair forms in both training and prediction: "weighting" (default) pads by repeating the last form and normalizes position weights within each cell pair; "mean_all" averages all combinations; "first" selects the first form; "off" disables form splitting. "min_all" is not supported for PMI training.
//' @param form_delim Literal delimiter separating forms within a cell, default "#". Empty forms are discarded before pairing.
//' @param form_weights Finite, non-negative position weights for "weighting". NULL generates 1, 0.5, 0.25, ... . The used prefix is normalized separately for each cell pair. Zero effective weight makes the cell pair invalid for both training and prediction.
//' @param quiet Whether to suppress all output messages.
//' @return A list containing the following components:
//' \item{result}{A dataframe of distances, either in long table form or square form.}
//' \item{cost}{The final cost matrix used for distance calculation.}
//' \item{sum_diff}{The sum of absolute differences between the cost matrices of the last two iterations.}
//' \item{mean_diff}{The mean of absolute differences between the cost matrices of the last two iterations.}
//' \item{converged}{Whether a completed cost update met the convergence tolerance.}
//' \item{n_epochs}{Number of cost updates actually completed.}
//' @examples
//' df <- as.data.frame(rbind(a=c("ph_l_i_z","k_o_l"),b=c("b_l_i_s", "k_a:_l")))
//' result <- pw_pmi_dist(df, delim="_")
//' result <- pw_pmi_dist(df, delim="_", squareform=TRUE)
//' result <- pw_pmi_dist(df, delim="_", n_threads=4)
//' forms <- data.frame(word=c("a_b#b", "b_a#a#b"))
//' result <- pw_pmi_dist(forms, delim="_", form_strategy="weighting")
//' result <- pw_pmi_dist(forms, delim="_", form_strategy="mean_all")
//[[Rcpp::export]]
List pw_pmi_dist(const DataFrame &data, const String &delim = "", bool detailed = false, const String &normalize_method = "longest", bool squareform = false, int n_threads = 4, int max_epochs = 20, double tol = 1e-4, double damping_alpha = 1.0, double smoothing_lambda = 0.001, double min_edit_cost = 0.01, const String &form_strategy = "weighting", const String &form_delim = "#", Nullable<NumericVector> form_weights = R_NilValue, bool quiet = false)
{
    return lingdist::pmi_df(data, delim, detailed, normalize_method, squareform, n_threads, max_epochs, tol, damping_alpha, smoothing_lambda, min_edit_cost,
                            form_strategy, form_delim, read_weights(form_weights, "form_weights"), quiet);
}

//' Compute Weighted Jaccard Distance between all row pairs of a dataframe
//'
//' Compute Weighted Jaccard Distance (WJD) between all row pairs of a dataframe. This metric is suitable for categorical data with hierarchical structure and multiple forms.
//' For example, a cell value "A_2_a#A_3#B_5_c" indicates there are 3 word forms for this lexical item. The first form "A_2_a" belongs to category A, subcategory 2, and sub-subcategory a. The second form "A_3" belongs to category A and subcategory 3. The third form "B_5_c" belongs to category B, subcategory 5, and sub-subcategory c. The delimiters can be customized via `form_delim` and `cate_delim`.
//'
//' @param data DataFrame with n rows and m columns indicating there are n languages or dialects to involve in the calculation and there are at most m words to base on, in which the rownames are the language ids.
//' @param cate_weights Numeric vector of weights for different levels of categories for a single word form. For example, if a word form is "A_2_a", the first weight applies to "A", the second to "2", and the third to "a". Weights must be finite, non-negative, and cover all observed levels. NULL generates weights 1, 0.5, 0.25, ... to the required depth.
//' @param form_weights Numeric vector of weights for different word forms of a single lexical item (ordered). For example, if a lexical item is "A_2_a#A_3#B_5_c", the first weight applies to "A_2_a", the second to "A_3", etc. Weights must be finite, non-negative, and cover all observed forms. NULL generates weights 1, 0.5, 0.25, ... to the required length.
//' @param form_strategy How to combine word forms: weighting, mean_all, min_all, first, or off.
//' @param form_delim Literal delimiter separating different word forms of a single lexical item. Default is "#".
//' @param cate_delim Literal delimiter separating different levels of categories within a word form. Default is "_".
//' @param detailed Whether to return detailed information. When FALSE (default), returns the average distance across all columns (lexical items) for each row pair. When TRUE, returns column-wise distances for each row pair, showing how distance varies across lexical items.
//' @param squareform Whether to return a dataframe in squareform.
//' @param n_threads Number of threads used for computation. Use 1 for sequential computation.
//' @param quiet Whether to suppress all output messages.
//' @return A dataframe in long table form if `squareform` is FALSE, otherwise in squareform.
//' @examples
//' df <- as.data.frame(rbind(a=c("A_1#B_2","C_3"),b=c("A_1","C_3_x")))
//' result <- pw_wjd(df, form_delim="#", cate_delim="_")
//[[Rcpp::export]]
DataFrame pw_wjd(const DataFrame &data, const String &cate_delim = "_", Nullable<NumericVector> cate_weights = R_NilValue, const String &form_strategy = "weighting", const String &form_delim = "#", Nullable<NumericVector> form_weights = R_NilValue, bool detailed = false, bool squareform = false, int n_threads = 2, bool quiet = false)
{
    std::vector<double> cate_weights_vec = read_weights(cate_weights, "cate_weights");
    std::vector<double> form_weights_vec = read_weights(form_weights, "form_weights");
    return lingdist::wjd_df(data, cate_delim, cate_weights_vec, form_strategy, form_delim, form_weights_vec, detailed, squareform, n_threads, quiet);
}
