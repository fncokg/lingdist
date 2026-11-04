#pragma once

#include "helpers.hpp"
#include "indexed_costs.hpp"

namespace lingdist
{
    // Import only costs in the supplied ID space, but validate the full matrix.
    // No intermediate string-indexed cost object or second dense copy is needed.
    IndexedCosts read_costs(const SymbolTable &symbols, const Rcpp::DataFrame &matrix);
    Rcpp::DataFrame write_costs(const IndexedCosts &costs);
    StrVec missing_cost_symbols(const SymbolTable &symbols, const Rcpp::DataFrame &matrix);
}
