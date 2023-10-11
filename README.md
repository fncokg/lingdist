# lingdist: A Fast Generalized Edit Distance Based Linguistic Distance and Alignment Computation R Package

`lingdist` is a fast generalized edit distance and string alignment computation mainly for linguistic aims. As a generalization to the classic edit distance algorithms, the package allows user to define custom cost for every symbol's insertion, deletion, and substitution. The package also allows character combinations in any length to be seen as a single symbol which is very useful for IPA transcriptions with diacritics. In addition to edit distance result, users can get detailed alignment information such as all possible alignment scenarios between two strings which is useful for testing ,illustration or any further usage such as PMI-based distance in linguistics. Either distance matrix or its long table form can be obtained and tools to do such conversions are provided. All functions in the package are implemented in C++ and the distance matrix computation is parallelized leveraging the `RcppThread` package.