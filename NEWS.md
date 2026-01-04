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
