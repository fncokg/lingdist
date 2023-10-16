# lingdist: Fast Linguistic Distance and Alignment Computation

`lingdist` is an R package mainly for efficient computations of linguistic distances. In its current version, a generalized edit distance framework is implemented. More linguitic distance algorithms will be implemented in the later version. The package originated from requirements in linguistic researches but has the flexibility to be used in any suitable situations.

# Functionality

In general:
- **Fast Generalized Edit Distance Computation**: All alogrithms as well as the parallel computation features and related tool functions are implemented in native `C++` leveraging the `Rcpp` package.
- **Parallel Computation for Large Datasets**: When pairwise distance matrix is required for large datasets, the parallel computation feature can be enabled as easy as setting the `parallel` parameter to `TRUE`.

For the generalized edit distance algorithm:
- **Custom Cost Assignments**: You can specify all cost values for any character pairs. This is useful when you want specify different cost values for a-i substitution and a-k substitution.
- **String as a Character**: The package is able to treate a string as a single symbol when computing edit distance. This is useful when you want the algorithm to consider the IPA transcription \[kʰʷ\] as a single symbol, rather than 3 separate characters.
- **Alignment Results**: In addition to distance result, the package also provides you a detailed algiment result report including which two symbols are algined by the algorithm and how many alignment scenarios are possible.

# Installation

Just

```
install.packages("lingdist")
```

and then

```
library(lingdist)
```

You can also download the newest dev codes in Github and compile yourself (only for experienced users).

# Usage

## Edit distance and alignment details for two strings

With `edit_dist_string` (cost_mat is optional and can be an empty dataframe):

```
res <- edit_dist_string("pʰ_l_i_z̥", "p_l_i_s", cost_mat = cost.mat, delim = "_")
```

Inspect `res$distance` for distance result, and `res$alignments` for alignment details, which may look like this:

```
[[1]]
  Chars1 Chars2  Operation Cost Cumcost
1     pʰ      p substitute    1       1
2      l      l       same    0       1
3      i      i       same    0       1
4      z̥      s substitute    1       2
```

## Average pairwise edit distances for all rows in a dataframe

When you have a dataframe looks like this:

```
                    0      1           2
albanian    pʰ_l_i_z̥  k_o_l  s_t_ɛ_l_aː
azerbaijani pʰ_l̥_i_z̥ k_ɑ_lˠ   s_t̪_ɛ_l_ə
bengali      p_l_i_s  k_o_l ə_s_t_e_l_ʌ
```

And you want all pairwise distances between rows (i.e. albanian vs azerbaijani, bengali vs azerbaijani etc.), call:

```
result <- edit_dist_df(df, cost_mat = cost.mat, delim = "_", parallel = TRUE, n_threads = 4)
```

In this call, parallel computation is enabled and 4 threads are used in the computation. For more details and other APIs, refer to the [documentation](https://cran.r-project.org/web/packages/lingdist/lingdist.pdf).

# How to cite

Forthcoming.