# lingdist: Fast Linguistic Distance and Alignment Computation

`lingdist` is a fast generalized edit distance and string alignment computation mainly for linguistic aims. As a generalization to the classic edit distance algorithms, the package allows users to define custom cost for every symbol's insertion, deletion, and substitution. The package also allows character combinations in any length to be seen as a single symbol which is very useful for International Phonetic Alphabet (IPA) transcriptions with diacritics. In addition to edit distance result, users can get detailed alignment information such as all possible alignment scenarios between two strings which is useful for testing, illustration or any further usage. Either the distance matrix or its long table form can be obtained and tools to do such conversions are provided. All functions in the package are implemented in `C++` and the distance matrix computation is parallelized leveraging the `RcppThread` package.

# Functionality
- Fast Generalized Edit Distance Computation: All alogrithms as well as the parallel computation features are implemented in native C++.
- Parallel Computation for Large Datasets: When pairwise distance matrix is required for large datasets, the parallel computation feature can be enabled as easy as setting the `parallel` parameter to `TRUE`.
- Custom Cost Assignments: You can specify all cost values for any character pairs. This is useful when you want specify different cost values for a-i substitution and a-k substitution.
- String as a Character: The package is able to treate a string as a single symbol when computing edit distance. This is useful when you want the algorithm to consider the IPA transcription \[kʰʷ\] as a single symbol, rather than 3 separate characters.
- Alignment Results: In addition to distance result, the package also provides you a detailed algiment result report including which two symbols are algined by the algorithm and how many alignment scenarios are possible.

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

In this call, parallel computation is enabled and 4 threads are used in the computation. For more details and other APIs, refer to the [documentation](https://cran.r-project.org/web/packages/RcppThread/RcppThread.pdf).

# How to cite

Forthcoming.