#define point std::pair<int, int>
#define points std::vector<point>
#define double_points std::pair<double, points>
#define cost_umap std::unordered_map<std::string, double>
#define str_vec std::vector<std::string>
#include <Rcpp.h>
#include <RcppThread.h>
#include <unordered_map>
#include <vector>
#include <set>
#include <string>
#include <utility>
using namespace Rcpp;

std::string EMPTY("_NULL_");
std::string DELIM("_#_");

////R functions////

/*
    Wrapper of R `strsplit` function.
*/
str_vec
split(const String &str, const String &delim)
{
    Function r_strsplit("strsplit");
    List result = r_strsplit(str, delim);
    return as<str_vec>(result[0]);
}

std::vector<std::vector<str_vec>> split_df(const DataFrame &data, const String &delim)
{
    int nrows = data.nrow(), ncols = data.ncol();
    std::vector<std::vector<str_vec>> rows_vector;
    for (int i = 0; i < nrows; i++)
    {
        std::vector<str_vec> row_vector;
        for (int j = 0; j < ncols; j++)
        {
            const str_vec &col = data[j];
            str_vec chars = split(col[i], delim);
            row_vector.push_back(chars);
        }
        rows_vector.push_back(row_vector);
    }
    return rows_vector;
}

cost_umap mat2umap(const DataFrame &cost_mat)
{
    cost_umap umap;
    String key;

    if (cost_mat.nrow() == 0 || cost_mat.ncol() == 0)
    {
        return umap;
    }
    StringVector _row_names = cost_mat.attr("row.names");
    StringVector _col_names = cost_mat.names();
    str_vec row_names = as<str_vec>(_row_names);
    str_vec col_names = as<str_vec>(_col_names);
    for (std::size_t col = 0; col < col_names.size(); col++)
    {
        NumericVector this_col = cost_mat[col];
        for (std::size_t row = 0; row < row_names.size(); row++)
        {
            std::string key = col_names[col] + DELIM + row_names[row];
            umap[key] = this_col[row];
        }
    }
    return umap;
}

void free_2d_arr_memory(double **arr, int nrow)
{
    for (int i = 0; i < nrow; i++)
    {
        delete[] arr[i];
    }
    delete[] arr;
}

double get_cost(const std::string &str1, const std::string &str2, const cost_umap &cost)
{
    std::string key = str1 + DELIM + str2;
    if (cost.find(key) != cost.end())
    {
        return cost.at(key);
    }
    else
    {
        return str1 == str2 ? 0 : 1;
    }
}

////Dist functions////

double edit_dist_core_dp(const str_vec &vec1, const str_vec &vec2, const cost_umap &cost)
{
    std::size_t len1 = vec1.size(), len2 = vec2.size(), i, j, nrow = len2 + 1, ncol = len1 + 1;

    double **dist = new double *[nrow];

    for (i = 0; i < nrow; i++)
    {
        dist[i] = new double[ncol];
    }

    dist[0][0] = 0;

    for (j = 1; j < ncol; j++)
    {
        dist[0][j] = dist[0][j - 1] + get_cost(EMPTY, vec1[j - 1], cost);
    }
    for (i = 1; i < nrow; i++)
    {
        dist[i][0] = dist[i - 1][0] + get_cost(EMPTY, vec2[i - 1], cost);
    }

    for (i = 1; i < nrow; i++)
    {
        for (j = 1; j < ncol; j++)
        {

            dist[i][j] = std::min({
                dist[i - 1][j - 1] + get_cost(vec1[j - 1], vec2[i - 1], cost), // sub
                dist[i][j - 1] + get_cost(vec1[j - 1], EMPTY, cost),           // vec1 insert
                dist[i - 1][j] + get_cost(EMPTY, vec2[i - 1], cost),           // vec2 insert
            });
        }
    }

    double result = dist[vec2.size()][vec1.size()];
    free_2d_arr_memory(dist, vec2.size() + 1);
    return result;
}

std::vector<std::vector<double_points>> edit_dist_core_dp_record(const str_vec &vec1, const str_vec &vec2, const cost_umap &cost)
{
    std::size_t len1 = vec1.size(), len2 = vec2.size(), i, j, nrow = len2 + 1, ncol = len1 + 1;

    std::vector<std::vector<double_points>> dist;
    for (i = 0; i < nrow; i++)
    {
        std::vector<double_points> this_row;
        for (j = 0; j < ncol; j++)
        {
            this_row.push_back(std::make_pair(0, points()));
        }
        dist.push_back(this_row);
    }

    for (j = 1; j < ncol; j++)
    {
        dist[0][j] = std::make_pair(
            dist[0][j - 1].first + get_cost(EMPTY, vec1[j - 1], cost),
            points({std::make_pair(0, j - 1)}));
    }
    for (i = 1; i < nrow; i++)
    {
        dist[i][0] = std::make_pair(
            dist[i - 1][0].first + get_cost(EMPTY, vec2[i - 1], cost),
            points({std::make_pair(i - 1, 0)}));
    }

    for (i = 1; i < nrow; i++)
    {
        for (j = 1; j < ncol; j++)
        {

            std::vector<double> possible_values({
                dist[i - 1][j - 1].first + get_cost(vec1[j - 1], vec2[i - 1], cost), // sub
                dist[i][j - 1].first + get_cost(vec1[j - 1], EMPTY, cost),           // vec1 insert
                dist[i - 1][j].first + get_cost(EMPTY, vec2[i - 1], cost),           // vec2 insert
            });

            double min_value = *std::min_element(possible_values.begin(), possible_values.end());
            points source_points;
            if (possible_values[0] == min_value)
            {
                source_points.push_back(std::make_pair(i - 1, j - 1));
            }
            if (possible_values[1] == min_value)
            {
                source_points.push_back(std::make_pair(i, j - 1));
            }
            if (possible_values[2] == min_value)
            {
                source_points.push_back(std::make_pair(i - 1, j));
            }
            dist[i][j] = std::make_pair(
                min_value,
                source_points);
        }
    }

    return dist;
}

std::vector<std::vector<point>> get_matched_paths(std::vector<std::vector<double_points>> dist)
{
    int nrow = dist.size(), ncol = dist[0].size(), this_x = nrow - 1, this_y = ncol - 1;
    std::vector<std::vector<point>> result;
    std::vector<std::pair<point, std::vector<point>>> stack;
    std::vector<point> this_path;
    while (true)
    {
        if (this_x == 0 && this_y == 0)
        {
            this_path.insert(this_path.begin(), std::make_pair(this_x, this_y));
            result.push_back(this_path);
            if (stack.size() == 0)
            {
                break;
            }
            else
            {
                auto stack_back = stack.back();
                stack.pop_back();
                this_x = stack_back.first.first;
                this_y = stack_back.first.second;
                this_path = stack_back.second;
            }
        }
        this_path.insert(this_path.begin(), std::make_pair(this_x, this_y));
        points next_points = dist[this_x][this_y].second;
        this_x = next_points[0].first;
        this_y = next_points[0].second;
        if (next_points.size() > 1)
        {
            for (std::size_t i = 1; i < next_points.size(); i++)
            {
                stack.push_back(std::make_pair(next_points[i], this_path));
            }
        }
    }
    return result;
}

double edit_dist_row(const std::vector<str_vec> &row1, const std::vector<str_vec> &row2, const cost_umap &cost)
{
    std::size_t ncols = row1.size();
    double currentDist, ttlDist = 0., meanDist;
    int nword = 0;
    for (std::size_t coli = 0; coli < ncols; coli++)
    {
        const str_vec &chars1 = row1[coli], &chars2 = row2[coli];

        if (chars1.size() != 0 && chars2.size() != 0)
        {
            currentDist = edit_dist_core_dp(chars1, chars2, cost);
            ttlDist += currentDist;
            nword += 1;
        }
    }
    if (nword <= 0)
    {
        meanDist = -1;
    }
    else
    {
        meanDist = ttlDist / nword;
    }

    return meanDist;
}
//' Convert long table to square form
//'
//' Convert a distance dataframe in long table form to a square matrix form.
//'
//' @param data Dataframe in long table form. The first and second columns are labels and the third column stores the distance values.
//' @param symmetric Whether the distance matrix are symmetric (if cost matrix is not, then the distance matrix is also not).
//' @return Dataframe in square matrix form, rownames and colnames are labels. If the long table only contains \eqn{C_n^2} rows and `symmetric` is set to FALSE, then only lower triangle positions in the result is filled.
//' @examples
//' data <- as.data.frame(list(chars1=c("a","a","b"),chars2=c("b","c","c"),dist=c(1,2,3)))
//' mat <- long2squareform(data)
//[[Rcpp::export]]
DataFrame long2squareform(const DataFrame &data, bool symmetric = true)
{
    StringVector chars1_col = data[0];
    StringVector chars2_col = data[1];
    NumericVector dist_col = data[2];

    // str_vec chars = as<str_vect>(chars1_col);
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
        char2idx[unique_chars[i]] = i;
        names.push_back(unique_chars[i]);
        NumericVector col(unique_chars.size());
        result.push_back(col);
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

    result.attr("row.names") = names;
    result.attr("names") = names;
    return result;
}

List get_string_alignment(const str_vec &chars1, const str_vec &chars2, const cost_umap &cost, const String &delim = "")
{
    List report;
    std::string char1, char2, opt;
    int this_x, this_y, prev_x, prev_y;
    std::vector<std::vector<double_points>> dist = edit_dist_core_dp_record(chars1, chars2, cost);
    std::vector<std::vector<point>> paths = get_matched_paths(dist);
    double result = dist[chars2.size()][chars1.size()].first;
    report["distance"] = result;
    List path_dfs;

    for (auto &path : paths)
    {
        StringVector chars1_col, chars2_col, operation_col;
        NumericVector cost_col, cumcost_col;
        for (std::size_t j = 1; j < path.size(); j++)
        {
            this_x = path[j].first;
            this_y = path[j].second;
            prev_x = path[j - 1].first;
            prev_y = path[j - 1].second;
            char1 = this_y == prev_y ? EMPTY : chars1[this_y - 1];
            char2 = this_x == prev_x ? EMPTY : chars2[this_x - 1];
            opt = char1 == EMPTY ? "insert" : (char2 == EMPTY ? "delete" : (char1 == char2 ? "same" : "substitute"));
            chars1_col.push_back(char1);
            chars2_col.push_back(char2);
            operation_col.push_back(opt);
            cumcost_col.push_back(dist[this_x][this_y].first);
            cost_col.push_back(get_cost(char1, char2, cost));
        }
        DataFrame path_df = DataFrame::create(Named("Chars1") = chars1_col, Named("Chars2") = chars2_col, Named("Operation") = operation_col, Named("Cost") = cost_col, Named("Cumcost") = cumcost_col);
        path_dfs.push_back(path_df);
    }
    report["alignments"] = path_dfs;
    return report;
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
//' @return A list contains `distance` attribution storing the distance result. If `return_alignments` is TRUE, then a `alignments` attribution is present which is a list of dataframes with each storing a possible best alignment scenario.
//' @examples
//' cost.mat <- data.frame()
//' dist <- edit_dist_string("leaf","leaves")$distance
//' dist <- edit_dist_string("ph_l_i_z","p_l_i_s",cost_mat=cost.mat,delim="_")$distance
//' alignments <- edit_dist_string("ph_l_i_z","p_l_i_s",delim="_",return_alignments=TRUE)$alignments
//[[Rcpp::export]]
List edit_dist_string(const String &str1, const String &str2, Nullable<DataFrame> cost_mat = R_NilValue, const String &delim = "", bool return_alignments = false)
{
    cost_umap cost;
    if (cost_mat.isNotNull())
    {
        DataFrame cost_mat_ = DataFrame(cost_mat);
        cost = mat2umap(cost_mat_);
    }
    str_vec chars1 = split(str1, delim);
    str_vec chars2 = split(str2, delim);
    if (return_alignments)
    {
        return get_string_alignment(chars1, chars2, cost, delim);
    }
    else
    {
        double result = edit_dist_core_dp(chars1, chars2, cost);
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
//' @param symmetric Whether to the result matrix is symmetric. This depends on whether the `cost_mat` is symmetric.
//' @param parallel Whether to parallelize the computation.
//' @param n_threads The number of threads is used to parallelize the computation. Only meaningful if `parallel` is TRUE.
//' @return A dataframe in long table form if `squareform` is FALSE, otherwise in squareform. If `symmetric` is TRUE, the long table form has \eqn{C_n^2} rows otherwise \eqn{n^2} rows.
//' @examples
//' df <- as.data.frame(rbind(a=c("a_bc_d","d_bc_a"),b=c("b_bc_d","d_bc_a")))
//' cost.mat <- data.frame()
//' result <- edit_dist_df(df, cost_mat=cost.mat, delim="_")
//' result <- edit_dist_df(df, cost_mat=cost.mat, delim="_", squareform=TRUE)
//' result <- edit_dist_df(df, cost_mat=cost.mat, delim="_", parallel=TRUE, n_threads=4)
//[[Rcpp::export]]
DataFrame edit_dist_df(const DataFrame &data, Nullable<DataFrame> cost_mat = R_NilValue, const String &delim = "", bool squareform = false, bool symmetric = true, bool parallel = false, int n_threads = 2) // const DataFrame &cost_mat
{
    cost_umap cost;
    if (cost_mat.isNotNull())
    {
        DataFrame cost_mat_ = DataFrame(cost_mat);
        cost = mat2umap(cost_mat_);
    }
    int nrows = data.nrow();
    CharacterVector rownames = data.attr("row.names");
    StringVector lab1Col, lab2Col;

    std::vector<std::pair<int, int>> row_pairs;
    std::vector<std::vector<str_vec>> rows_vector = split_df(data, delim);
    if (symmetric)
    {
        for (int rowi = 0; rowi < nrows - 1; rowi++)
        {
            for (int rowj = rowi + 1; rowj < nrows; rowj++)
            {
                String label1 = rownames[rowi];
                String label2 = rownames[rowj];
                row_pairs.push_back(std::pair<int, int>{rowi, rowj});
                lab1Col.push_back(label1);
                lab2Col.push_back(label2);
            }
        }
    }
    else
    {
        for (int rowi = 0; rowi < nrows; rowi++)
        {
            for (int rowj = 0; rowj < nrows; rowj++)
            {
                String label1 = rownames[rowi];
                String label2 = rownames[rowj];
                row_pairs.push_back(std::pair<int, int>{rowi, rowj});
                lab1Col.push_back(label1);
                lab2Col.push_back(label2);
            }
        }
    }

    std::vector<double> dists(row_pairs.size());
    RcppThread::ProgressBar bar(row_pairs.size(), 1);
    if (parallel)
    {
        RcppThread::parallelFor(
            0,
            row_pairs.size(),
            [&](std::int32_t idx)
            {
                std::pair<int, int> &row_pair = row_pairs[idx];
                bar++;
                int rowi = row_pair.first;
                int rowj = row_pair.second;

                dists[idx] = edit_dist_row(rows_vector[rowi], rows_vector[rowj], cost);
            },
            n_threads);
    }
    else
    {
        for (std::size_t idx = 0; idx < row_pairs.size(); idx++)
        {
            std::pair<int, int> &row_pair = row_pairs[idx];
            bar++;
            int rowi = row_pair.first;
            int rowj = row_pair.second;

            dists[idx] = edit_dist_row(rows_vector[rowi], rows_vector[rowj], cost);
        }
    }
    DataFrame result = DataFrame::create(Named("lab1") = lab1Col, Named("lab2") = lab2Col, Named("dist") = NumericVector::import(dists.begin(), dists.end()));
    if (squareform)
    {
        result = long2squareform(result, symmetric);
    }
    return result;
}

str_vec get_all_unique_chars(const DataFrame &data, const String &delim = "")
{
    std::vector<std::vector<str_vec>> rows_vector = split_df(data, delim);
    str_vec all_chars;
    for (auto &row : rows_vector)
    {
        for (auto &chars : row)
        {
            all_chars.insert(all_chars.end(), chars.begin(), chars.end());
        }
    }

    std::set<std::string> unique_chars(all_chars.begin(), all_chars.end());
    str_vec all_unique_chars(unique_chars.begin(), unique_chars.end());
    return all_unique_chars;
}

//' Check whether there's missing characters in the cost matrix.
//'
//' Check whether there's missing characters in the cost matrix and return the missing characters.
//'
//' @param data DataFrame to be computed.
//' @param cost_mat Cost matrix to be checked.
//' @param delim The delimiter separating atomic symbols.
//' @return A string vector containing the missing characters, empty indicating there's no missing characters.
//' @examples
//' df <- as.data.frame(rbind(a=c("a_bc_d","d_bc_a"),b=c("b_bc_d","d_bc_a")))
//' cost.mat <- data.frame()
//' chars.not.found <- check_cost_defined(df, cost.mat, "_")
//[[Rcpp::export]]
StringVector check_cost_defined(const DataFrame &data, const DataFrame &cost_mat, const String &delim = "")
{

    str_vec all_chars = get_all_unique_chars(data, delim);
    all_chars.push_back(EMPTY);
    cost_umap umap = mat2umap(cost_mat);
    StringVector char_pair_not_found;
    for (auto &char1 : all_chars)
    {
        for (auto &char2 : all_chars)
        {
            if (char1 != char2)
            {
                std::string key = char1 + DELIM + char2;
                if (umap.find(key) == umap.end())
                {
                    char_pair_not_found.push_back(key);
                }
            }
        }
    }

    return char_pair_not_found;
}

//' Generate a default cost matrix
//'
//' generate a default cost matrix contains all possible characters in the raw data with all diagonal values set to 0 and others set to 1. This avoids you constructing the matrix from scratch.
//'
//' @param data DataFrame to be computed.
//' @param delim The delimiter separating atomic symbols.
//' @return Cost matrix contains all possible characters in the raw data with all diagonal values set to 0 and others set to 1.
//' @examples
//' df <- as.data.frame(rbind(a=c("a_bc_d","d_bc_a"),b=c("b_bc_d","d_bc_a")))
//' default.cost <- generate_default_cost_matrix(df, "_")
//[[Rcpp::export]]
DataFrame generate_default_cost_matrix(const DataFrame &data, const String &delim = "")
{
    DataFrame mat = DataFrame::create();
    str_vec all_chars = get_all_unique_chars(data, delim);
    all_chars.push_back(EMPTY);

    for (auto &char1 : all_chars)
    {
        NumericVector col;
        for (auto &char2 : all_chars)
        {
            col.push_back(char1 == char2 ? 0 : 1);
        }
        mat.push_back(col);
    }
    mat.attr("row.names") = wrap(all_chars);
    mat.attr("names") = wrap(all_chars);

    return mat;
}
