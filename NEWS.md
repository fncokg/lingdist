# lingdist 2.3.0 2026-01-06

Fixed several bugs and improved performance.

Bug fixes:
- In `alignment.cpp`, fixed a bug that caused incorrect alignment results. All of the aligned sequences were accumulated into the last sequence.
- Fixed all possible issues on double types comparison.
- Fixed memory leak caused by RcppThread::ProgressBar not being destructed properly.
- In some corner cases, when nwords==0, all distance functions would return negative distances. Now they return `NA`.
- In `dist_wjd.cpp`, fixed a bug that caused incorrect WJD distance. `wjd_form` would never return 0.0 in the multi-category case.
- In `dist_wjd.cpp`, added checks for vector lengths of `multi_form_weights` and `multi_form_cats`.
- In `dist_wjd.cpp`, added checks for `multi_form_weights` sum being greater than 0.
- In `cost_table.cpp`, added checks for symbols existing in both row and column names of the cost matrix.

Performance improvements:
- In `alignment.cpp`, rewrited the DFS algorithm. Avoid insert to the front of a vector, which is slow.
- In `dist_pmi.cpp`, rewrited the alignment counting logic with a new structure `AlignmentCounter`. Now we use a dense vector to store counts instead of a hash map, and a `vector<string>` for string-to-index searching. This is faster when the number of unique symbols is small, which is usually the case.

# lingdist 2.2.1 2026-01-04

Add `quite` parameters to all distance functions to control printing messages and progress bars.

# lingdist 2.2.0 2025-12-26

Several improvements have been made:
- `long2squareform` fills diagonal with 0 now; this can be specified with `default_diag` parameter.
- `CostTable` is reconstructed to support both a normal and a fast implementation. In the fast case, `get_cost` will not lookup the cost table but directly return values to fasten classic edit distance computation.
- `pw_edit_dist` now:
  - checks missing symbols in cost matrix by default; this can be disabled with `check_missing_cost = FALSE`.
  - builds fast cost table when cost matrix is not provided; the fast table costs can be specified with `default_sub_cost` and `default_ins_del_cost` parameters.
- report more information in `pw_pmi_dist`.

# lingdist 2.1.0 2025-12-26

Added WJD distance supporting hierarchical categorical data with multiple forms.

# lingdist 2.0.0 2025-12-22

- Reconstructed the whole C++ codes for edit distance calculation.
- Added PMI distance supporting.

# lingdist 1.0.0

* Initial CRAN submission.
