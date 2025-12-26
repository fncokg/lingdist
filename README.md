# lingdist: Fast Linguistic Distance and Alignment Computation

`lingdist` is an R package mainly for efficient computations of linguistic distances. It implements generalized edit distance and PMI (Pointwise Mutual Information) distance frameworks. The package originated from requirements in linguistic research but has the flexibility to be used in any suitable situations.

# Functionality

In general:
- **Fast Computation**: All algorithms as well as the parallel computation features and related tool functions are implemented in native `C++` leveraging the `Rcpp` package.
- **Parallel Computation for Large Datasets**: When pairwise distance matrix is required for large datasets, the parallel computation feature can be enabled as easy as setting the `parallel` parameter to `TRUE`.

Supported Algorithms:
- **Generalized Edit Distance**: Supports custom cost assignments, treating strings as single symbols (e.g., IPA transcriptions), and returning detailed alignment results.
- **PMI Distance**: Computes distances based on Pointwise Mutual Information, automatically learning costs from the data using an EM algorithm.
- **Weighted Jaccard Distance**: Computes distances for categorical data with hierarchical structure and multiple forms (e.g., "A_2_a#A_3").

# Installation

To install the stable version from CRAN, run:

```r
install.packages("lingdist")
```

For the development version from GitHub, run (NOTE: in this case, you need to compile the C++ source code, so make sure you have the necessary build tools installed):

```r
devtools::install_github("fncokg/lingdist")
```


# Usage

## Edit distance and alignment details for two strings

With `string_edit_dist` (cost_mat is optional and can be an empty dataframe):

```r
# cost.mat is optional
res <- string_edit_dist("pʰ_l_i_z̥", "p_l_i_s", cost_mat = NULL, delim = "_")
```

Inspect `res$distance` for distance result, and `res$alignments` for alignment details (if `return_alignments = TRUE`).

## Average pairwise edit distances for all rows in a dataframe

When you have a dataframe looks like this:

```r
df <- as.data.frame(rbind(
  albanian = c("pʰ_l_i_z̥", "k_o_l", "s_t_ɛ_l_aː"),
  azerbaijani = c("pʰ_l̥_i_z̥", "k_ɑ_lˠ", "s_t̪_ɛ_l_ə"),
  bengali = c("p_l_i_s", "k_o_l", "ə_s_t_e_l_ʌ")
))
```

And you want all pairwise distances between rows (i.e. albanian vs azerbaijani, bengali vs azerbaijani etc.), call:

```r
# cost.mat is optional
result <- pw_edit_dist(df, cost_mat = NULL, delim = "_", parallel = TRUE, n_threads = 4)
```

In this call, parallel computation is enabled and 4 threads are used in the computation.

## Pairwise PMI distances

To compute PMI distances which learn costs from the data:

```r
result <- pw_pmi_dist(df, delim = "_", parallel = TRUE, n_threads = 4)
```

## Pairwise Weighted Jaccard Distances

To compute Weighted Jaccard Distances for hierarchical categorical data:

```r
# Data with hierarchical structure (e.g. Category_Subcategory) and multiple forms (separated by #)
df_wjd <- as.data.frame(rbind(
  lang1 = c("A_1#B_2", "C_3"),
  lang2 = c("A_1", "C_3_x")
))

# Compute WJD
result <- pw_wjd(df_wjd, form_delim = "#", cate_delim = "_", parallel = TRUE)
```

## Utilities

- `generate_default_cost_matrix(data, delim)`: Generate a default cost matrix from data.
- `long2squareform(data, symmetric, default_diag)`: Convert long table distance result to square matrix.
- `check_cost_mat_symbols(cost_mat, data, delim)`: Check which symbols in `data` are missing from a cost matrix.

For more details and other APIs, refer to the documentation.

# How to cite

Forthcoming.