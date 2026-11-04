#include "cost_io.hpp"

#include <cmath>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>

namespace
{
    // Get the row and column names of a square cost matrix,
    // validating that they are (1) present, (2) unique, and (3) identical.
    struct MatrixNames
    {
        lingdist::StrVec rows;
        lingdist::StrVec columns;
    };

    MatrixNames matrix_names(const Rcpp::DataFrame &matrix)
    {
        if (matrix.nrow() != matrix.ncol())
            throw std::invalid_argument("Cost matrix must be square.");
        const Rcpp::CharacterVector rows = matrix.attr("row.names");
        const Rcpp::CharacterVector columns = matrix.names();
        if (rows.size() != matrix.nrow() || columns.size() != matrix.ncol())
            throw std::invalid_argument("Cost matrix must have row and column names.");
        for (const auto &names : {rows, columns})
            for (const auto name : names)
                if (Rcpp::CharacterVector::is_na(name))
                    throw std::invalid_argument("Cost matrix names must not be missing.");
        MatrixNames result{Rcpp::as<lingdist::StrVec>(rows), Rcpp::as<lingdist::StrVec>(columns)};
        const std::set<std::string> row_set(result.rows.begin(), result.rows.end());
        const std::set<std::string> col_set(result.columns.begin(), result.columns.end());
        if (row_set.size() != result.rows.size() || col_set.size() != result.columns.size())
            throw std::invalid_argument("Cost matrix names must be unique.");
        if (row_set != col_set)
            throw std::invalid_argument("Cost matrix row and column symbols must match.");
        return result;
    }
}

namespace lingdist
{
    IndexedCosts read_costs(const SymbolTable &symbols, const Rcpp::DataFrame &matrix)
    {
        const auto names = matrix_names(matrix);
        const auto size = symbols.size();
        if (size > std::numeric_limits<std::size_t>::max() / size)
            throw std::length_error("Indexed cost dimensions are too large.");
        // Fill undefined operations once at the boundary: identical symbols 0,
        // all other pairs 1. The kernel never does string lookup or fallback.
        std::vector<double> values(size * size, 1.0);
        for (std::size_t i = 0; i < size; ++i)
            values[i * size + i] = 0.0;
        // A general-purpose matrix may include symbols absent from this input.
        // Project onto its ID space rather than growing the dense cost array.
        std::vector<std::optional<TokenId>> rows, columns;
        for (const auto &name : names.rows)
            rows.push_back(symbols.find(name));
        for (const auto &name : names.columns)
            columns.push_back(symbols.find(name));
        for (std::size_t c = 0; c < columns.size(); ++c)
        {
            const Rcpp::NumericVector column = Rcpp::as<Rcpp::NumericVector>(matrix[c]);
            for (std::size_t r = 0; r < rows.size(); ++r)
            {
                // Even unused rows/columns must satisfy the cost contract.
                if (!std::isfinite(column[r]))
                    throw std::invalid_argument("Edit costs must be finite.");
                if (rows[r] && columns[c])
                    values[static_cast<std::size_t>(*columns[c]) * size + *rows[r]] = column[r];
            }
        }
        return IndexedCosts(symbols, std::move(values));
    }

    Rcpp::DataFrame write_costs(const IndexedCosts &costs)
    {
        Rcpp::List columns(costs.nsyms());
        for (std::size_t target = 0; target < costs.nsyms(); ++target)
        {
            Rcpp::NumericVector column(costs.nsyms());
            for (std::size_t source = 0; source < costs.nsyms(); ++source)
                column[source] = costs.get(source, target);
            columns[target] = column;
        }
        columns.attr("names") = Rcpp::wrap(costs.symbol_table().names());
        columns.attr("row.names") = Rcpp::wrap(costs.symbol_table().names());
        columns.attr("class") = "data.frame";
        return Rcpp::DataFrame(columns);
    }

    StrVec missing_cost_symbols(const SymbolTable &symbols, const Rcpp::DataFrame &matrix)
    {
        const auto declared = matrix_names(matrix).rows;
        const std::set<std::string> present(declared.begin(), declared.end());
        StrVec missing;
        for (const auto &name : symbols.names())
            if (!present.count(name))
                missing.push_back(name);
        return missing;
    }
}
